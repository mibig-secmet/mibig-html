# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
import shutil
import tempfile
from typing import Any, cast, Dict, List, Optional, Union

import antismash
from antismash.common import logs
from antismash.common.record_processing import (
    pre_process_sequences as _as_pre_process_sequences,
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
    bio_records = [record.to_biopython() for record in results.records]

    # add antismash meta-annotation to records
    add_antismash_comments(list(zip(results.records, bio_records)), options)

    logging.debug("Writing cluster-specific genbank files")
    for record, bio_record in zip(results.records, bio_records):
        for region in record.get_regions():
            region.write_to_genbank(directory=options.output_dir, record=bio_record)

    # write records to an aggregate output
    base_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    combined_filename = base_filename + ".gbk"
    logging.debug("Writing final genbank file to '%s'", combined_filename)
    SeqIO.write(bio_records, combined_filename, "genbank")

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
    annotation_filename = "{}.json".format(_get_mibig_acc(options))
    shutil.copy(options.mibig_json, os.path.join(options.output_dir, annotation_filename))


def _get_mibig_acc(options: ConfigType) -> str:
    """Get the MIBiG accession from the input file."""
    return os.path.splitext(os.path.basename(options.mibig_json))[0]


def mibig_rename_records(records: List[Record], options: ConfigType) -> None:
    """ Rename records according to MIBiG identifiers. """
    mibig_acc = _get_mibig_acc(options)

    for record in records:
        record.id = "{}.1".format(mibig_acc)
        record.name = mibig_acc
        record.annotations['accessions'].insert(0, mibig_acc)


def pre_processor_wrapper(sequences: List[Record], options: ConfigType, genefinding: AntismashModule) -> List[Record]:
    sequences = _as_pre_process_sequences(sequences, options, genefinding)
    mibig_rename_records(sequences, options)
    return sequences


antismash.main.write_outputs = write_outputs
antismash.main.record_processing.pre_process_sequences = pre_processor_wrapper
antismash.main.record_processing.Record = Record
antismash.main.run_detection = run_mibig_detection


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
        result = antismash.run_antismash(sequence_file, options)
        antismash.main.get_analysis_modules = orig_analysis_modules
    return result
