# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Responsible for creating the single web page results """

import string
import os
import re
from typing import Any, Dict, List, Tuple, Union, cast

from mibig.converters.read.cluster import GeneAnnotation

from antismash.common import html_renderer, path
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.json import JSONOrf
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.custom_typing import AntismashModule
from antismash.detection import hmm_detection, nrps_pks_domains
from antismash.config import ConfigType
from antismash.outputs.html import js
from antismash.outputs.html.generator import (
    find_plugins_for_cluster,
    generate_searchgtr_htmls,
    write_regions_js,
    VISUALISERS,
)

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
                    options: ConfigType, all_modules: List[AntismashModule]) -> Tuple[
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
        update_cds_description(record, json_record, cds_annotations)

        json_record['seq_id'] = "".join(char for char in json_record['seq_id'] if char in string.printable)
        for region, json_region in zip(record.get_regions(), json_record['regions']):
            json_region["product_categories"] = convert_categories(mibig_results.data.cluster.biosynthetic_class)
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

def generate_html_sections(records: List[RecordLayer], results: Dict[str, Dict[str, ModuleResults]],
                           options: ConfigType) -> Dict[str, Dict[int, List[HTMLSections]]]:
    """ Generates a mapping of record->region->HTMLSections for each record, region and module

        Arguments:
            records: a list of RecordLayers to pass through to the modules
            results: a dictionary mapping record name to
                        a dictionary mapping each module name to its results object
            options: the current antiSMASH config

        Returns:
            a dictionary mapping record id to
                a dictionary mapping region number to
                    a list of HTMLSections, one for each module
    """
    details = {}
    for record in records:
        record_details = {}
        record_result = results[record.id]
        for region in record.regions:
            # work around mibig module not creating protoclusters with the expected types
            assert len(region.subregions) == 1 and region.subregions[0].tool == "mibig"
            if nrps_pks_domains.domain_drawing.has_domain_details(region.region_feature):
                region.handlers.append(cast(AntismashModule, nrps_pks_domains))

            sections = []
            for handler in region.handlers:
                if handler is nrps_pks_domains or handler.will_handle(region.products, region.product_categories):
                    handler_results = record_result.get(handler.__name__)
                    if handler_results is None:
                        continue
                    sections.append(handler.generate_html(region, handler_results, record, options))
            record_details[region.get_region_number()] = sections
        details[record.id] = record_details
    return details


def generate_webpage(records: List[Record], results: List[Dict[str, ModuleResults]],
                     options: ConfigType, all_modules: List[AntismashModule]) -> None:
    """ Generates and writes the HTML itself """
    # bring mibig module to the front of the module list
    all_modules = [cast(AntismashModule, annotations)] + [module for module in all_modules if module is not annotations]

    generate_searchgtr_htmls(records, options)
    json_records, js_domains, js_results = build_json_data(records, results, options, all_modules)
    write_regions_js(json_records, options.output_dir, js_domains, js_results)

    with open(os.path.join(options.output_dir, 'index.html'), 'w') as result_file:
        template = FileTemplate(path.get_full_path(__file__, "templates", "overview.html"))

        options_layer = OptionsLayer(options, all_modules)
        record_layers_with_regions = []
        record_layers_without_regions = []
        results_by_record_id: Dict[str, Dict[str, ModuleResults]] = {}
        for record, record_results in zip(records, results):
            if record.get_regions():
                record_layers_with_regions.append(RecordLayer(record, None, options_layer))
            else:
                record_layers_without_regions.append(RecordLayer(record, None, options_layer))
            results_by_record_id[record.id] = record_results

        regions_written = sum(len(record.get_regions()) for record in records)
        job_id = os.path.basename(options.output_dir)

        mibig_id = os.path.splitext(os.path.basename(options.mibig_json))[0]
        annotation_filename = "{}.json".format(mibig_id)
        page_title = mibig_id

        html_sections = generate_html_sections(record_layers_with_regions, results_by_record_id, options)

        svg_tooltip = ("Shows the layout of the region, marking coding sequences and areas of interest. "
                       "Clicking a gene will select it and show any relevant details. "
                       "Clicking an area feature (e.g. a candidate cluster) will select all coding "
                       "sequences within that area. Double clicking an area feature will zoom to that area. "
                       "Multiple genes and area features can be selected by clicking them while holding the Ctrl key."
                       )
        record_layer = record_layers_with_regions[0]
        region = record_layer.regions[0]
        mibig_results = results[0][annotations.__name__]
        assert isinstance(mibig_results, annotations.MibigAnnotations)
        aux = template.render(records=record_layers_with_regions, options=options_layer,
                              version=options.version, extra_data=js_domains,
                              regions_written=regions_written, sections=html_sections,
                              config=options, job_id=job_id, page_title=page_title,
                              records_without_regions=record_layers_without_regions,
                              svg_tooltip=svg_tooltip,
                              record=record_layer, region=region, cluster=mibig_results.data.cluster,
                              annotation_filename=annotation_filename, mibig_id=mibig_id)
        result_file.write(aux)


def update_cds_description(record: Record, js_record: Dict[str, Any], annotations: List[GeneAnnotation]) -> None:
    id_to_annotations = {}
    for annotation in annotations:
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
