# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Responsible for creating the single web page results """

import string
import os
import re
from typing import Any, Dict, List, Set, Tuple, Union, cast

from mibig.converters.read.cluster import GeneAnnotation

from antismash.common import html_renderer, path
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.json import JSONOrf
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.custom_typing import AntismashModule
from antismash.detection import hmm_detection
from antismash.config import ConfigType
from antismash.outputs.html import js
from antismash.outputs.html.generator import (
    find_plugins_for_cluster,
    generate_searchgtr_htmls,
    write_regions_js,
    VISUALISERS,
)

from mibig.converters.read.top import Everything

from mibig_html import annotations
from mibig_html.common.layers import OptionsLayer


def convert_categories(categories: List[str]) -> List[str]:
    """ Converts a list of MIBiG biosynthetic classes to antiSMASH categories """
    translations = {
        "NRP": "NRPS",
        "Polyketide": "PKS",
    }
    known = {category.name for category in hmm_detection.categories.get_rule_categories()}

    new = []
    for category in categories:
        if category in known:
            new.append(category)
            continue
        if category.lower() in known:
            new.append(category.lower())
            continue
        replacement = translations.get(category)
        if not replacement:
            raise ValueError(f"cannot translate MIBiG class {category} into a valid antiSMASH category")
        new.append(replacement)
    return new


def build_json_data(records: List[Record], results: List[Dict[str, ModuleResults]],
                    options: ConfigType, all_modules: List[AntismashModule],
                    categories: Set[str]) -> Tuple[
                        List[Dict[str, Any]],
                        List[Dict[str, Union[str, List[JSONOrf]]]],
                        Dict[str, Dict[str, Dict[str, Any]]]
                    ]:
    """ Builds JSON versions of records and domains for use in drawing SVGs with
        javascript.

        Arguments:
            records: a list of Records to convert
            results: a dictionary mapping record id to a list of ModuleResults to convert
            options: antiSMASH options
            categories: a list of antiSMASH compatible category strings

        Returns:
            a tuple of
                a list of JSON-friendly dicts representing records
                a list of JSON-friendly dicts representing domains
    """
    js_records = js.convert_records(records, results, options)

    js_domains: List[Dict[str, Union[str, List[JSONOrf]]]] = []
    js_results = {}
    assert len(records) == 1
    assert len(records[0].get_regions()) == 1
    mibig_results = results[0][annotations.__name__]
    assert isinstance(mibig_results, annotations.mibig.MibigAnnotations)

    for i, record in enumerate(records):
        json_record = js_records[i]
        # replace antismash cds_detail with mibig's one
        try:
            cds_annotations = mibig_results.data.cluster.genes.annotations
        except AttributeError:
            cds_annotations = []
        update_cds_description(record, json_record, cds_annotations, mibig_results)

        json_record['seq_id'] = "".join(char for char in json_record['seq_id'] if char in string.printable)
        for region, json_region in zip(record.get_regions(), json_record['regions']):
            json_region["product_categories"] = sorted(categories)
            handlers = find_plugins_for_cluster(all_modules, json_region)
            region_results = {}
            for handler in handlers:
                # if there's no results for the module, don't let it try
                if handler.__name__ not in results[i]:
                    continue
                if "generate_js_domains" in dir(handler):
                    domains_by_region = handler.generate_js_domains(region, record)
                    if domains_by_region:
                        js_domains.append(domains_by_region)
                if hasattr(handler, "generate_javascript_data"):
                    data = handler.generate_javascript_data(record, region, results[i][handler.__name__])
                    region_results[handler.__name__] = data

            for aggregator in VISUALISERS:
                if not hasattr(aggregator, "generate_javascript_data"):
                    continue
                if aggregator.has_enough_results(record, region, results[i]):
                    data = aggregator.generate_javascript_data(record, region, results[i])
                    region_results[aggregator.__name__] = data

            if region_results:
                js_results[RegionLayer.build_anchor_id(region)] = region_results

    return js_records, js_domains, js_results


def generate_html_sections(record: RecordLayer, results: Dict[str, ModuleResults],
                           options: ConfigType, modules: List[AntismashModule],
                           categories: Set[str]) -> List[HTMLSections]:
    """ Generates a mapping of record->region->HTMLSections for each record, region and module

        Arguments:
            records: a list of RecordLayers to pass through to the modules
            results: a dictionary mapping record name to
                        a dictionary mapping each module name to its results object
            options: the current antiSMASH config
            modules: modules that may have results for the BGC
            categories: a set of antiSMASH-compatible product categories

        Returns:
            a list of HTMLSections, one for each relevant module
    """
    assert len(record.regions) == 1
    region = record.regions[0]
    sections = []
    # don't use the region/record layer handlers for generating sections,
    # as they don't use the converted categories for checking
    for module in modules:
        if module.__name__ not in results:
            continue
        if hasattr(module, "will_handle") and module.will_handle([], categories):
            sections.append(module.generate_html(region, results[module.__name__], record, options))
    return sections


def generate_webpage(record: Record, result: Dict[str, ModuleResults],
                     options: ConfigType, all_modules: List[AntismashModule]) -> None:
    """ Generates and writes the HTML itself """
    mibig_results = result[annotations.__name__]
    assert isinstance(mibig_results, annotations.MibigAnnotations)
    categories = set(convert_categories(mibig_results.data.cluster.biosynthetic_class))

    # bring mibig module to the front of the module list
    all_modules.pop(all_modules.index(annotations))
    all_modules.insert(0, cast(AntismashModule, annotations))

    generate_searchgtr_htmls([record], options)
    json_records, js_domains, js_results = build_json_data([record], [result], options, all_modules, categories)
    write_regions_js(json_records, options.output_dir, js_domains, js_results)

    template = FileTemplate(path.get_full_path(__file__, "templates", "overview.html"))

    options_layer = OptionsLayer(options, all_modules)
    record_layer = RecordLayer(record, None, options_layer)

    mibig_id = os.path.splitext(os.path.basename(options.mibig_json))[0]
    annotation_filename = "{}.json".format(mibig_id)

    sections = generate_html_sections(record_layer, result, options, all_modules, categories)

    svg_tooltip = ("Shows the layout of the region, marking coding sequences and areas of interest. "
                   "Clicking a gene will select it and show any relevant details. "
                   "Clicking an area feature (e.g. a candidate cluster) will select all coding "
                   "sequences within that area. Double clicking an area feature will zoom to that area. "
                   "Multiple genes and area features can be selected by clicking them while holding the Ctrl key."
                   )
    region = record_layer.regions[0]
    aux = template.render(records=[record_layer], options=options_layer,
                          version=options.version, extra_data=js_domains,
                          regions_written=1, sections={record.id: {1: sections}},
                          config=options, page_title=mibig_id,
                          svg_tooltip=svg_tooltip,
                          record=record_layer, region=region, cluster=mibig_results.data.cluster,
                          annotation_filename=annotation_filename, mibig_id=mibig_id)
    with open(os.path.join(options.output_dir, 'index.html'), 'w', encoding="utf_8") as result_file:
        result_file.write(aux)


def generate_retired_page(data: Everything, options: ConfigType) -> None:
    template = FileTemplate(path.get_full_path(__file__, "templates", "retired.html"))

    options_layer = OptionsLayer(options, [])
    mibig_id = data.cluster.mibig_accession

    aux = template.render(options=options_layer,
                          reasons=data.cluster.retirement_reasons, see_also=data.cluster.see_also,
                          page_title=mibig_id, mibig_id=mibig_id)
    with open(os.path.join(options.output_dir, 'index.html'), 'w', encoding="utf_8") as result_file:
        result_file.write(aux)


def update_cds_description(record: Record, js_record: Dict[str, Any],
                           cds_annotations: List[GeneAnnotation], results: annotations.MibigAnnotations) -> None:
    original_accession = results.data.cluster.loci.accession

    id_to_annotations = {}
    for annotation in cds_annotations:
        if annotation.id:
            assert annotation.id not in id_to_annotations, annotation.id
            id_to_annotations[annotation.id] = annotation
        if annotation.name:
            id_to_annotations[annotation.name] = annotation

    template = html_renderer.FileTemplate(path.get_full_path(__file__, "templates", "cds_detail.html"))

    for reg_idx, js_region in enumerate(js_record["regions"]):
        for cds_idx, js_cds in enumerate(js_region["orfs"]):
            locus = js_cds["locus_tag"]
            cds = record.get_cds_by_name(locus)
            annotation = None
            for name in [cds.locus_tag, cds.protein_id, cds.gene]:
                annotation = id_to_annotations.get(name)
                if annotation:
                    break
            rendered = template.render(annotation=annotation)
            js_cds["description"] = re.sub(r"(\(total: (.+) nt\)<br>)", r"\1{}".format(rendered), str(js_cds["description"]))
            if original_accession.startswith("MIBIG"):
                js_cds["description"] = re.sub(r"<a[^>]+>View genomic context</a>", "", js_cds["description"])
            else:
                js_cds["description"] = re.sub(record.id, original_accession, js_cds["description"])
