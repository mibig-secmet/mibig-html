# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Overrides some antiSMASH HTML rendering functions """

from typing import Any, List

from antismash.common.html_renderer import (
    collapser_start,
    collapser_end,
    FileTemplate as _FileTemplate,
    Markup,
    switch,
    cds_selector_span,
)
from antismash.common.secmet import Record
from mibig.converters.shared.common import Citation, Evidence

_TOOLTIP_COUNTER = 0


def help_tooltip(text: str, name: str, inline: bool = False) -> Markup:
    """ Constructs a help icon with tooltip, each will have a unique ID generated
        based on the given name.

        Arguments:
            text: the content of the tooltip
            name: a prefix for id generation

        Returns:
            A Markup instance with the constructed HTML
    """
    global _TOOLTIP_COUNTER  # pylint: disable=global-statement
    _TOOLTIP_COUNTER += 1
    unique_id = "%s-help-%d" % (name, _TOOLTIP_COUNTER)
    return Markup(('<div class="help-container{2}">'
                   ' <div class="help-icon" data-id="{0}"></div>'
                   ' <span class="help-tooltip" id="{0}">{1}</span>'
                   '</div>').format(unique_id, text, "-inline" if inline else ""))


def clickable_gene(name: str, record: Record, force_current: bool = False, real_name: str = None,
                   allow_missing: bool = False,
                   ) -> Markup:
    """ A template for HTML that will highlight the relevant CDS in the overview,
        using the "gene" name if possible.

        Arguments:
            name: any identifier for the CDS
            record: the containing Record
            force_current: if True, uses the given name instead of the gene name in display
            real_name: the unique name for the CDS if known, otherwise it will be looked up
            allow_missing: allows missing gene names to default back to the text instead
                           of raising an error

        Returns:
            A Markup instance with the constructed HTML
    """

    if not real_name:
        try:
            real_name = record.get_real_cds_name(name)
        except ValueError as err:
            assert "unknown CDS" in str(err)
            return Markup(f"<span>{name}</span>")
    assert real_name is not None
    if force_current:
        gene_name = name
    else:
        cds = record.get_cds_by_name(real_name)
        gene_name = cds.gene or real_name
    return cds_selector_span(real_name, display_name=gene_name)


def clickable_gene_list(names: List[str], record: Record,
                        force_current: bool = False, separator: str = " ", allow_missing: bool = False,
                        ) -> Markup:
    """ A template for HTML that will highlight relevant CDS features in the overview,
        using the "gene" name if possible.

        Arguments:
            names: any identifier for the CDS features
            record: the containing Record
            force_current: if True, uses the given name instead of the gene name in display
            separator: the separator to use in the output
            allow_missing: allows missing gene names to default back to the text instead
                           of raising an error

        Returns:
            A Markup instance with the constructed HTML
    """
    return Markup(separator.join(clickable_gene(name, record, force_current=force_current) for name in names))


def build_evidence_list(evidence: list[Evidence], use_superscript: bool = False) -> Markup:
    """ Builds markup for a list of evidence methods, including any supporting references.

        Arguments:
            evidence: the evidence to build markup from

        Returns:
            a Markup instance with the constructed HTML
    """
    sections = []
    for section in evidence:
        sections.append(f"{section.method}{build_short_form_citation_links(section.references, use_superscript=use_superscript)}")
    return Markup(", ".join(sections))


def build_short_form_citation_links(citations: list[Citation], use_superscript: bool = False) -> Markup:
    """ Builds markup for a list of citations in short form as superscript.

        Arguments:
            citations: the citations to use in the list

        Returns:
            a Markup instance with the constructed HTML
    """
    assert all(citation.short_id for citation in citations)
    citations = sorted(set(citations))
    chunks = [f' <a href="{citation.to_url()}">[{citation.short_id}]</a>' for citation in citations]
    return Markup(f'{"<sup>" if use_superscript else ""}{", ".join(chunks)}{"</sup>" if use_superscript else ""}')


class FileTemplate(_FileTemplate):  # pylint: disable=too-few-public-methods
    """ A template renderer for file templates """
    def render(self, **kwargs: Any) -> Markup:
        """ Renders the template HTML, providing any given arguments to the renderer """
        if not self.template:
            raise ValueError("attempting to render without a template")

        defaults = {
            "build_evidence_list": build_evidence_list,
            "build_short_form_citation_links": build_short_form_citation_links,
            "collapser_start": collapser_start,
            "collapser_end": collapser_end,
            "help_tooltip": help_tooltip,
            "switch": switch,
            "clickable_gene": clickable_gene,
            "clickable_gene_list": clickable_gene_list,
        }
        defaults.update(kwargs)
        return Markup(self.template.render(**defaults))
