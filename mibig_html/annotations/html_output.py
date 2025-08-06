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
from mibig.converters.shared.common import Citation

from mibig_html.common.html_renderer import FileTemplate, build_short_form_citation_links
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
    entry = results.entry

    html = HTMLSections("mibig-general")
    publications_links = ReferenceCollection(entry.references, results.pubmed_cache, results.doi_cache)

    general = render_template("general.html", entry=results.entry,
                              publications_links=publications_links.get_links())
    html.add_detail_section("General", general)

    compounds = render_template("compounds.html", compounds=results.entry.compounds)
    html.add_detail_section("Compounds", compounds, class_name="mibig-compounds")

    genes = []
    annots = entry.genes.annotations if entry.genes and entry.genes.annotations else []
    genes_have_comments = False
    for cds_feature in region_layer.cds_children:
        gene = {
            "real_name": cds_feature.get_name(),
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
        gene["functions"] = {}
        for function in cds_feature.gene_functions:
            function_text = str(function.function)
            if function.tool != "mibig":
                continue
            if function.tool == "rule-based-clusters":
                function_text += " ({})".format(function.description)
            elif function.tool == "smcogs":
                function_text += " ({})".format(function.description.split(" (")[0])
            gene["functions"][function_text] = []
        annot_idx = -1
        for i, annot in enumerate(annots):
            annot_names = {str(annot.id)}
            if annot.name:
                annot_names.add(str(annot.name))
            feature_names = {cds_feature.locus_tag, cds_feature.protein_id, cds_feature.gene}
            if feature_names.intersection(annot_names):
                annot_idx = i
                break
        if annot_idx >= 0:
            annot = annots.pop(annot_idx)
            for function in annot.functions if annot.functions else []:
                function_text = function.function
                references_by_method = {}
                for evidence in function.evidence:
                    references_by_method.setdefault(evidence.method, []).extend(function.references)
                links = {}
                for method, references in references_by_method.items():
                    links[method] = build_short_form_citation_links(references)
                gene["functions"][function_text] = links
            for tailoring in annot.tailoring_functions if annot.tailoring_functions else []:
                if not tailoring.db_reference:
                    continue
                ref = tailoring.db_reference.split(":")[-1]
                function_text = f"{tailoring.function}"
                links = {}
                if tailoring.references:
                    links[ref] = build_short_form_citation_links(tailoring.references)
                gene["functions"][function_text] = links
            if annot.mutation_phenotype:
                function_text = "Mutation phenotype: {}".format(annot.mutation_phenotype)
                gene["functions"][function_text] = []
            if annot.product:
                gene["product"] = annot.product
            if annot.comment:
                gene["comment"] = annot.comment
                genes_have_comments = True
        genes.append(gene)
    # everything should've been popped
    assert not annots, [str(a) for a in annots]

    html.add_detail_section("Genes", render_template("genes.html", genes=genes, genes_have_comments=genes_have_comments, record=record_layer),
                            class_name="mibig-genes")

    if entry.biosynthesis:
        domains_with_substrates = set()
        for module in entry.biosynthesis.modules:
            for domain in module.extra_info.get_domains():
                if domain.substrates:
                    domains_with_substrates.add(domain)
        rendered = render_template(
            "biosynthesis.html", biosynthesis=entry.biosynthesis,
            domains_with_substrates=domains_with_substrates, record=record_layer,
        )
        html.add_detail_section("Biosynthesis", rendered, class_name="mibig-biosynthesis")

    if entry.biosynthesis.operons:
        rendered = render_template(
            "operons.html", operons=entry.biosynthesis.operons, record=record_layer,
        )
        html.add_detail_section("Operons", rendered, class_name="mibig-operons")

    #if entry.cluster.polyketide:
    #    html.add_detail_section("Polyketide", render_template("polyketide.html", pk=results.entry.cluster.polyketide, record=record_layer),
    #                            class_name="mibig-polyketide")

    #if entry.cluster.nrp:
    #    html.add_detail_section("NRP", render_template("nrp.html", nrp=results.entry.cluster.nrp, record=record_layer),
    #                            class_name="mibig-nrp")

    #if entry.cluster.ripp:
    #    html.add_detail_section("RiPP", render_template("ripp.html", ripp=results.entry.cluster.ripp, record=record_layer),
    #                            class_name="mibig-ripp")

    #if entry.cluster.saccharide:
    #    html.add_detail_section("Saccharide", render_template("saccharide.html", sac=results.entry.cluster.saccharide, record=record_layer),
    #                            class_name="mibig-saccharide")

    #if entry.cluster.terpene:
    #    html.add_detail_section("Terpene", render_template("terpene.html", trp=results.entry.cluster.terpene, record=record_layer),
    #                            class_name="mibig-terpene")

    html.add_detail_section("History", render_template("logs.html", releases=entry.changelog.releases),
                            class_name="mibig-logs")

    return html


class ReferenceLink:
    """Keep track of a single reference link."""

    __slots__ = (
        'title',
        'citation',
        'info',
    )

    def __init__(self, citation: Citation, *, title: str = "", info: str = "",
                 ) -> None:
        self.citation = citation
        self.title = title
        self.info = info

    @property
    def category(self) -> str:
        return self.citation.database

    @property
    def ref(self) -> str:
        return self.citation.to_url()


class ReferenceCollection:
    """Keep track of all references in a MIBiG entry."""

    __slots__ = (
        'client',
        'doi_cache',
        'pubmed_cache',
        'references',
    )

    def __init__(self, publications: list[Citation], pubmed_cache: PubmedCache,
                 doi_cache: DoiCache) -> None:
        self.client: Client = None
        self.references = {}
        self.pubmed_cache = pubmed_cache
        self.doi_cache = doi_cache
        pmids = []
        dois = []

        for publication in publications:
            if publication.database == "pubmed":
                if publication.value == "0":
                    continue
                pmids.append(publication.value)
            elif publication.database == "doi":
                dois.append(publication.value)

            self.references[publication.value] = ReferenceLink(publication)

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
