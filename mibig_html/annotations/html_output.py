# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the MIBiG sideloader """

import os
from typing import Any, List, Set

from eutils import Client

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.html_renderer import HTMLSections, Markup
from antismash.common.layers import RegionLayer, RecordLayer
from mibig.converters.read.cluster import Publication

from mibig_html.common.html_renderer import FileTemplate
from mibig_html.common.layers import OptionsLayer

from .mibig import MibigAnnotations, PubmedCache, DoiCache


def will_handle(_products: List[str], _categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return True


def render_template(template_name: str, **kwargs: Any) -> Markup:
    """ Returns the markup from rendering the given template with the given keyword arguments

        Arguments:
            template_name: the filename of the template within the 'templates' subdirectory
            kwargs: keyword arguments to pass into the template for rendering

        Returns:
            the resulting Markup object from rendering
    """
    template = FileTemplate(path.get_full_path(__file__, "templates", template_name))
    return template.render(**kwargs)


def generate_html(region_layer: RegionLayer, results: ModuleResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    assert isinstance(results, MibigAnnotations)
    data = results.data
    tax = results.taxonomy
    # "class" is a reserved keyword in python, can't use it directly
    tax_class = getattr(tax, "class")

    html = HTMLSections("mibig-general")
    taxonomy_text = f"{tax.superkingdom} > {tax.kingdom} > {tax.phylum} > {tax_class} > {tax.order} > {tax.family} > {tax.name}"
    publications_links = ReferenceCollection(data.cluster.publications, results.pubmed_cache, results.doi_cache)

    general = render_template("general.html", data=results.data, taxonomy_text=taxonomy_text,
                              publications_links=publications_links.get_links())
    html.add_detail_section("General", general)

    compounds = render_template("compounds.html", compounds=results.data.cluster.compounds)
    html.add_detail_section("Compounds", compounds, class_name="mibig-compounds")

    genes = []
    annots = data.cluster.genes.annotations if data.cluster.genes else []
    for cds_feature in region_layer.cds_children:
        gene = {
            "locus_tag": cds_feature.locus_tag,
            "protein_id": cds_feature.protein_id,
            "gene": cds_feature.gene,
            "start": cds_feature.location.start + 1,
            "end": cds_feature.location.end,
            "strand": "+" if cds_feature.location.strand > 0 else "-",
            "product": cds_feature.product,
            "aa_seq": cds_feature.translation,
            "nt_seq": cds_feature.extract(record_layer.seq)
        }
        gene["functions"] = []
        for function in cds_feature.gene_functions:
            function_text = str(function.function)
            if function.tool != "mibig":
                continue
            if function.tool == "rule-based-clusters":
                function_text += " ({})".format(function.description)
            elif function.tool == "smcogs":
                function_text += " ({})".format(function.description.split(" (")[0])
            gene["functions"].append(function_text)
        annot_idx = -1
        for i, annot in enumerate(annots):
            annot_names = {annot.name, annot.id}
            feature_names = {cds_feature.locus_tag, cds_feature.protein_id, cds_feature.gene}
            if feature_names.intersection(annot_names):
                annot_idx = i
                break
        if annot_idx >= 0:
            annot = annots.pop(annot_idx)
            for function in annot.functions:
                function_text = function.category
                if annot.tailoring:
                    function_text += " ({}) ".format(", ".join(annot.tailoring))
                else:
                    function_text += " "
                gene["functions"].append(function_text)
                gene["evidences"] = sorted(set(function.evidence))
            if annot.mutation_phenotype:
                function_text = "Mutation phenotype: {}".format(annot.mutation_phenotype)
                gene["functions"].append(function_text)
            if annot.product:
                gene["product"] = annot.product
        genes.append(gene)
    for annot in annots:
        gene = {
            "locus_tag": annot.id,
            "protein_id": "None",
            "gene": annot.name,
            "product": annot.product or ""
        }
        gene["functions"] = []
        for function in annot.functions:
            function_text = function.category
            if annot.tailoring:
                function_text += " ({}) ".format(", ".join(annot.tailoring))
            else:
                function_text += " "
            gene["functions"].append(function_text)
            gene["evidences"] = sorted(set(function.evidence))
        genes.append(gene)
    html.add_detail_section("Genes", render_template("genes.html", genes=genes, record=record_layer),
                            class_name="mibig-genes")

    if data.cluster.polyketide:
        html.add_detail_section("Polyketide", render_template("polyketide.html", pk=results.data.cluster.polyketide, record=record_layer),
                                class_name="mibig-polyketide")

    if data.cluster.nrp:
        html.add_detail_section("NRP", render_template("nrp.html", nrp=results.data.cluster.nrp, record=record_layer),
                                class_name="mibig-nrp")

    if data.cluster.ripp:
        html.add_detail_section("RiPP", render_template("ripp.html", ripp=results.data.cluster.ripp, record=record_layer),
                                class_name="mibig-ripp")

    if data.cluster.saccharide:
        html.add_detail_section("Saccharide", render_template("saccharide.html", sac=results.data.cluster.saccharide, record=record_layer),
                                class_name="mibig-saccharide")

    if data.cluster.terpene:
        html.add_detail_section("Terpene", render_template("terpene.html", trp=results.data.cluster.terpene, record=record_layer),
                                class_name="mibig-terpene")

    logs = sorted(results.data.changelog, key=lambda log: log.mibig_version)
    html.add_detail_section("History", render_template("logs.html", logs=logs),
                            class_name="mibig-logs")

    return html


class ReferenceLink:
    """Keep track of a single reference link."""

    __slots__ = (
        'category',
        'ref',
        'title',
        'info',
    )

    def __init__(self, category: str, reference: str, title: str, info: str = None) -> None:
        self.category = category
        self.ref = reference
        self.title = title
        self.info = info


class ReferenceCollection:
    """Keep track of all references in a MIBiG entry."""

    __slots__ = (
        'client',
        'doi_cache',
        'pubmed_cache',
        'references',
    )

    def __init__(self, publications: List[Publication], pubmed_cache: PubmedCache,
                 doi_cache: DoiCache) -> None:
        self.client: Client = None
        self.references = {}
        self.pubmed_cache = pubmed_cache
        self.doi_cache = doi_cache
        pmids = []
        dois = []

        for publication in publications:
            if publication.category == "pubmed":
                if publication.content == "0":
                    continue
                reference = "https://www.ncbi.nlm.nih.gov/pubmed/{}".format(publication.content)
                pmids.append(publication.content)
            elif publication.category == "patent":
                reference = "https://patents.google.com/patent/{}".format(publication.content)
            elif publication.category == "doi":
                dois.append(publication.content)
                reference = "https://dx.doi.org/{}".format(publication.content)
            elif publication.category == "url":
                reference = publication.content

            self.references[publication.content] = ReferenceLink(
                publication.category, reference, publication.content)

        self._resolve_pmids(pmids)
        self._resolve_dois(dois)

    def get_links(self) -> List[ReferenceLink]:
        return list(self.references.values())

    def _resolve_pmids(self, pmids: List[str]) -> None:
        if not pmids:
            return

        missing = self.pubmed_cache.get_missing(pmids)

        if missing:
            if self.client is None:
                self.client = Client(api_key=os.environ.get("NCBI_API_KEY", None))
            articles = self.client.efetch(db="pubmed", id=missing)
            for article in articles:
                self.pubmed_cache.add(article.title, article.authors,
                                      article.year, article.jrnl, article.pmid)

        for pmid in pmids:
            entry = self.pubmed_cache.get(pmid)
            self.references[pmid].title = entry.title
            self.references[pmid].info = entry.info

    def _resolve_dois(self, dois: List[str]) -> None:
        for doi in dois:
            entry = self.doi_cache.get(doi)
            self.references[doi].title = entry.title
            self.references[doi].info = entry.info
