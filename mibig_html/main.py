# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from datetime import datetime
import json
import logging
import os
import shutil
import tempfile
from typing import Any, cast, Dict, List, Optional, Union

import antismash
from antismash.common import logs
from antismash.common.record_processing import (
    AntismashInputError,
    records_contain_shotgun_scaffolds,
    generate_unique_id,
    sanitise_sequence,
    strip_record,
    _strict_parse,
)
from antismash.common.secmet.locations import (
    FeatureLocation,
    locations_overlap,
)
from antismash.detection import (
    cluster_hmmer,
    genefunctions,
    hmm_detection,
    nrps_pks_domains,
)
from antismash.main import (
    add_antismash_comments,
    AntismashModule,
    canonical_base_filename,
    ConfigType,
    DetectionResults,
    ModuleResults,
    run_module,
    serialiser,
    SeqIO,
    svg,
)
from mibig.converters.read.top import Everything

from mibig_html import annotations, html
from mibig_html.common.secmet import Record


__version__ = "3.0alpha"


def get_all_modules() -> List[AntismashModule]:
    """ Return a list of all modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    all_modules = []
    for modules in [
        antismash.main.get_detection_modules(),
        antismash.main.get_analysis_modules(),
        antismash.main.get_support_modules(),
    ]:
        all_modules.extend(modules)
    # add the pseudo-detection module for mibig annotations
    all_modules.append(annotations)
    # don't use "get_output_modules" as only one html module is relevant
    all_modules.append(html)
    return all_modules


def run_mibig_detection(record: Record, options: ConfigType,
                        module_results: Dict[str, Union[ModuleResults, Dict[str, Any]]]) -> Dict[str, float]:
    """ Detect different secondary metabolite clusters, PFAMs, and domains.

        Arguments:
            record: the Record to run detection over
            options: antiSMASH config
            module_results: a dictionary mapping a module's name to results from
                            a previous run on this module, as a ModuleResults subclass
                            or in JSON form

        Returns:
            the time taken by each detection module as a dictionary
    """
    timings: Dict[str, float] = {}

    logging.info("Loading MIBiG annotations")
    run_module(record, cast(AntismashModule, annotations), options, module_results, timings)
    results = module_results.get(annotations.__name__)
    if not results:
        raise ValueError("failed to load MiBIG annotations")
    assert isinstance(results, DetectionResults)
    for protocluster in results.get_predicted_protoclusters():
        record.add_protocluster(protocluster)
    for region in results.get_predicted_subregions():
        record.add_subregion(region)

    # annotate biosynthetic domains
    run_module(record, cast(AntismashModule, hmm_detection), options, module_results, timings)

    record.create_candidate_clusters()
    record.create_regions()

    if not record.get_regions():
        logging.info("No regions detected, skipping record")
        record.skip = "No regions detected"
        return timings

    logging.info("%d region(s) detected in record", len(record.get_regions()))

    # finally, run any detection limited to genes in regions
    for module in [nrps_pks_domains, cluster_hmmer, genefunctions]:
        run_module(record, cast(AntismashModule, module), options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            assert isinstance(results, ModuleResults)
            logging.debug("Adding detection results from %s to record", module.__name__)
            results.add_to_record(record)

    return timings


def write_outputs(results: serialiser.AntismashResults, options: ConfigType) -> None:
    """ Write output files (webpage, genbank files, etc) to the output directory

        Arguments:
            results: a serialiser.AntismashResults instance
            options: an antismash Config instance

        Returns:
            None
    """
    # don't use results for which the module no longer exists to regenerate/calculate
    module_results_per_record = []
    for record_results in results.results:
        record_result = {}
        for module_name, result in record_results.items():
            if isinstance(result, ModuleResults):
                record_result[module_name] = result
        module_results_per_record.append(record_result)

    logging.debug("Creating results page")
    html.write(results.records, module_results_per_record, options, get_all_modules())

    logging.debug("Creating results SVGs")
    svg.write(options, module_results_per_record)

    # convert records to biopython
    assert len(results.records) == 1
    record = results.records[0]
    bio_record = record.to_biopython()

    # add antismash meta-annotation to records
    add_antismash_comments(list(zip(results.records, [bio_record])), options)

    # write records to an aggregate output
    assert len(record.get_regions()) == 1
    region = record.get_regions()[0]
    base_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    filename = f"{base_filename}.gbk"
    logging.debug("Writing final genbank file to '%s'", filename)
    region.write_to_genbank(filename=filename, directory=options.output_dir, record=bio_record)

    zipfile = base_filename + ".zip"
    if os.path.exists(zipfile):
        os.remove(zipfile)
    if not options.skip_zip_file:
        logging.debug("Zipping output to '%s'", zipfile)
        with tempfile.NamedTemporaryFile(prefix="as_zip_tmp", suffix=".zip") as temp:
            shutil.make_archive(temp.name.replace(".zip", ""), "zip", root_dir=options.output_dir)
            shutil.copy(temp.name, zipfile)
            os.chmod(zipfile, 0o644)
        assert os.path.exists(zipfile)

    logging.debug("Saving mibig annotation file")
    annotation_filename = "annotations.json"
    shutil.copy(options.mibig_json, os.path.join(options.output_dir, annotation_filename))


def _get_mibig_acc(options: ConfigType) -> str:
    """Get the MIBiG accession from the input file."""
    return os.path.splitext(os.path.basename(options.mibig_json))[0]


def mibig_rename_records(records: List[Record], options: ConfigType) -> None:
    """ Rename records according to MIBiG identifiers. """
    mibig_acc = _get_mibig_acc(options)

    for record in records:
        record.id = "{}".format(mibig_acc)
        record.name = mibig_acc
        record.annotations['accessions'].insert(0, mibig_acc)


def pre_process_sequences(sequences: List[Record], options: ConfigType) -> List[Record]:
    """ Custom preprocess to remove restrictions from normal antiSMASH runs

        Arguments:
            sequences: the secmet.Record instances to process
            options: an antismash Config instance

        Returns:
            A list of altered secmet.Record
    """
    logging.debug("Preprocessing %d sequences", len(sequences))
    assert len(sequences) == 1

    # catch WGS master or supercontig entries
    if records_contain_shotgun_scaffolds(sequences):
        raise AntismashInputError("incomplete whole genome shotgun records are not supported")

    for i, seq in enumerate(sequences):
        seq.record_index = i + 1  # 1-indexed

    checking_required = not (options.reuse_results or options.skip_sanitisation)

    # keep sequences as clean as possible and make sure they're valid
    if checking_required:
        logging.debug("Sanitising record ids and sequences")
        # Ensure all records have unique names
        all_record_ids = {seq.id for seq in sequences}
        if len(all_record_ids) < len(sequences):
            all_record_ids = set()
            for record in sequences:
                if record.id in all_record_ids:
                    record.original_id = record.id
                    record.id = generate_unique_id(record.id, all_record_ids)[0]
                all_record_ids.add(record.id)
            assert len(all_record_ids) == len(sequences), "%d != %d" % (len(all_record_ids), len(sequences))
        sequences = [sanitise_sequence(sequences[0])]

    for record in sequences:
        if record.skip or not record.seq:
            logging.warning("Record %s has no sequence, skipping.", record.id)
        if not record.id:
            raise AntismashInputError("record has no name")

    if all(sequence.skip for sequence in sequences):
        raise AntismashInputError("all records skipped")

    mibig_rename_records(sequences, options)
    return sequences


def parse_input_sequence(filename: str, taxon: str = "bacteria",
                         start: int = 0, end: int = 0) -> List[Record]:
    """ Parse input records contained in a file

        Arguments:
            filename: the path of the file to read
            taxon: the taxon of the input, e.g. 'bacteria', 'fungi'
            start: a start location for trimming the sequence, or -1 to use all
            end: an end location for trimming the sequence, or -1 to use all

        Returns:
            A list of secmet.Record instances, one for each record in the file
    """
    logging.info('Parsing input sequence %r', filename)

    records = _strict_parse(filename)
    assert len(records) == 1
    record = records[0]

    if not Record.is_nucleotide_sequence(record.seq):
        raise AntismashInputError("protein records are not supported: %s" % record.id)

    # before conversion to secmet records, remove any irrelevant CDS features if possible
    if start > 0 or end != 0:
        if end == 0:
            end = len(record.seq)
        location = FeatureLocation(start, end)
        logging.critical(f"removing all CDS features outside area: {location}")
        features = []
        for feature in record.features:
            if feature.type != "CDS":  # keep all non-CDS features
                features.append(feature)
            elif locations_overlap(location, feature.location):
                features.append(feature)
        record.features = features

    # remove any previous or obselete antiSMASH annotations to minimise incompatabilities
    strip_record(record)

    logging.debug("Converting records from biopython to secmet")
    secmet_record = Record.from_biopython(record, taxon)
    # if parsable by secmet, it has a better context on what to strip, so run
    # the secmet stripping to ensure there's no surprises
    secmet_record.strip_antismash_annotations()

    return [secmet_record]


def _run_mibig(sequence_file: Optional[str], options: ConfigType) -> int:
    """ The real run_mibig, assumes logging is set up around it """
    options.all_enabled_modules = list(filter(lambda x: x.is_enabled(options), get_all_modules()))

    antismash.main.check_prerequisites(options.all_enabled_modules, options)

    # ensure the provided options are valid
    if not antismash.main.verify_options(options, options.all_enabled_modules):
        return 1

    assert annotations.is_enabled(options)

    start_time = datetime.now()

    with open(options.mibig_json) as handle:
        data = Everything(json.load(handle))
    loci = data.cluster.loci
    start = loci.start - 1 if loci.start else 0
    end = loci.end or 0

    if data.cluster.status == "retired":

        html.write_retired(data, options)

        running_time = datetime.now() - start_time
        logging.debug("MIBiG HTML generation finished at %s; runtime: %s",
                      datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(running_time))

        logging.info("MIBiG status: SUCCESS")
        return 0

    if sequence_file:
        records = parse_input_sequence(sequence_file, options.taxon, start, end)
        results = serialiser.AntismashResults(sequence_file.rsplit(os.sep, 1)[-1],
                                              records, [{} for i in records],
                                              __version__, taxon=options.taxon)
    else:
        orig_record = antismash.common.serialiser.Record
        antismash.common.serialiser.Record = Record
        results = antismash.main.read_data(sequence_file, options)
        antismash.common.serialiser.Record = orig_record

    # reset module timings
    results.timings_by_record.clear()

    antismash.main.prepare_output_directory(options.output_dir, sequence_file or options.reuse_results)

    results.records = pre_process_sequences(results.records, options)

    mibig_rename_records(results.records, options)
    analysis_modules = [cast(AntismashModule, antismash.modules.clusterblast)]

    for record, module_results in zip(results.records, results.results):
        # skip if we're not interested in it
        if record.skip:
            continue
        logging.info("Analysing record: %s", record.id)
        timings = run_mibig_detection(record, options, module_results)
        assert record.get_regions()
        analysis_timings = antismash.main.analyse_record(record, options, analysis_modules, module_results)
        timings.update(analysis_timings)
        results.timings_by_record[record.id] = timings

    # Write results
    json_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    json_filename += ".json"
    logging.debug("Writing json results to '%s'", json_filename)
    results.write_to_file(json_filename)

    # now that the json is out of the way, annotate the record
    # otherwise we could double annotate some areas
    antismash.main.annotate_records(results)

    # create relevant output files
    write_outputs(results, options)

    running_time = datetime.now() - start_time

    logging.debug("MIBiG HTML generation finished at %s; runtime: %s",
                  datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(running_time))

    logging.info("MIBiG status: SUCCESS")
    return 0


def run_mibig(sequence_file: Optional[str], options: ConfigType) -> int:
    """ Reads in data, runs detection and analysis modules over any records
        found, then outputs the results to file.

        Arguments:
            sequence_file: the sequence file to read in records from, can be
                            None if reusing results
            options: command line options

        Returns:
            0 if requested operations completed succesfully, otherwise 1
            Exceptions may also be raised
    """

    with logs.changed_logging(logfile=options.logfile, verbose=options.verbose,
                              debug=options.debug):
        orig_analysis_modules = antismash.main.get_analysis_modules
        antismash.main.get_analysis_modules = lambda: [antismash.modules.clusterblast]
        result = _run_mibig(sequence_file, options)
        antismash.main.get_analysis_modules = orig_analysis_modules
    return result
