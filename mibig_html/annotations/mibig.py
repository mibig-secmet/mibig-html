# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" MiBIG specific sideloading """

import json
from typing import Any, Dict, List, Optional
import logging

from antismash.common.errors import AntismashInputError
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, SubRegion
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    location_contains_other,
)

from mibig.converters.shared.mibig import MibigEntry

from mibig_html.common.secmet import Record

from .references import DoiCache, PubmedCache


class MibigAnnotations(DetectionResults):
    def __init__(self, record_id: str, area: SubRegion, entry: MibigEntry, cache_file: str,
                 pubmed_cache_file: str, doi_cache_file: str) -> None:
        super().__init__(record_id)
        self.entry = entry  # holds the original annotation json data
        # save calculated loci (relative to record), not annotated ones
        self.record_id = record_id
        self.area = area

        self.taxonomy = entry.taxonomy

        self.pubmed_cache = PubmedCache(pubmed_cache_file)
        self.doi_cache = DoiCache(doi_cache_file)

    def get_predicted_subregions(self) -> List[SubRegion]:
        return [self.area]

    def to_json(self) -> Dict[str, Any]:
        # save only information critical for deciding reusability
        locus = self.entry.loci[0]
        annotations = []
        for annot in (self.entry.genes.annotations if self.entry.genes and self.entry.genes.annotations else []):
            annotations.append(annot.to_json())
        to_add = []
        for extra_gene in (self.entry.genes.to_add if self.entry.genes and self.entry.genes.to_add else []):
            to_add.append(extra_gene.to_json())
        to_delete = []
        for gene_to_remove in (self.entry.genes.to_delete if self.entry.genes and self.entry.genes.to_delete else []):
            to_delete.append(gene_to_remove.to_json())
        return {
            "record_id": self.record_id,
            "genbank_accession": locus.accession,
            "coords": (locus.location.begin or -1, locus.location.end or -1),
            "gene_annotations": annotations,
            "extra_genes": to_add,
        }

    @staticmethod
    def from_json(prev: Dict[str, Any], record: Record, annotations_file: str,
                  cache_file: str, pubmed_cache_file: str, doi_cache_file: str) -> Optional["MibigAnnotations"]:
        with open(annotations_file) as handle:
            raw = json.load(handle)
            entry = MibigEntry.from_json(raw)

        # compare old vs new annotation, decide if we can 'reuse'
        can_reuse = True
        locus = entry.loci[0]
        gene_annotations = entry.genes.annotations if entry.genes and entry.genes.annotations else []
        to_add = entry.genes.to_add if entry.genes and entry.genes.to_add else []
        to_delete = entry.genes.to_delete if entry.genes and entry.genes.to_delete else []
        if locus.accession != prev["genbank_accession"]:
            logging.debug("Previous result's genbank_accession is not the same as the new one")
            can_reuse = False
        elif record.id != prev["record_id"]:
            logging.debug("Previous result's record_id is not the same as the new one")
            can_reuse = False
        elif (locus.location.begin or -1) != prev["coords"][0] or (locus.location.end or -1) != prev["coords"][1]:
            logging.debug("Previous result's start/end coordinate is not the same as the new one")
            can_reuse = False
        elif len(gene_annotations) != len(prev["gene_annotations"]):
            logging.debug("Gene annotations have changed")
            can_reuse = False
        elif len(to_add) != len(prev["extra_genes"]):
            logging.debug("Additional genes have changed")
            can_reuse = False
        elif len(to_delete) != len(prev["to_delete"]):
            logging.debug("Genes to delete have changed")
            can_reuse = False

        # if we can't reuse, stop running antismash, because CDS annotations won't be correct
        if can_reuse:
            product = ", ".join([bc.class_name.value for bc in entry.biosynthesis.classes])
            loci_region = FeatureLocation(
                locus.location.begin - 1 if locus.location.begin else 0,
                locus.location.end or len(record.seq)
            )
            area = SubRegion(loci_region, tool="mibig", label=product)
            return MibigAnnotations(record.id, area, entry, cache_file, pubmed_cache_file, doi_cache_file)
        else:
            logging.error("Can't reuse MIBiG annotation.")
            raise AntismashInputError(
                "Genbank record or gene annotations are updated, can't reuse result")


def mibig_loader(annotations_file: str, cache_file: str, pubmed_cache_file: str,
                 doi_cache_file: str, record: Record) -> MibigAnnotations:
    """This method will be called only when not reusing data"""
    with open(annotations_file) as handle:
        raw = json.load(handle)
        entry = MibigEntry.from_json(raw)

    product = ", ".join([bc.class_name.value for bc in entry.biosynthesis.classes])
    locus = entry.loci[0]
    loci_region = FeatureLocation(
        locus.location.begin - 1 if locus.location.begin else 0,
        locus.location.end or len(record.seq)
    )
    area = SubRegion(loci_region, tool="mibig", label=product)

    # re-annotate CDSes in the Record
    if entry.genes:
        # extra genes
        for gene in entry.genes.to_add if entry.genes.to_add else []:
            if gene.id and gene.location:
                exons = [FeatureLocation(exon.begin - 1, exon.end, strand=gene.location.strand)
                         for exon in gene.location.exons]
                location = CompoundLocation(exons) if len(exons) > 1 else exons[0]
                if not location_contains_other(area.location, location):
                    raise ValueError(f"additional gene {gene.id} lies outside cluster")
                translation = gene.translation or record.get_aa_translation_from_location(location)
                cds_feature = CDSFeature(
                    location=location, locus_tag=gene.id, translation=translation)
                record.add_cds_feature(cds_feature, auto_deduplicate=False)
                record.add_alteration(f"{cds_feature.get_name()} was added")
        # re-annotation
        for cds_feature in record.get_cds_features_within_location(area.location):
            locus_tag = cds_feature.locus_tag
            protein_id = cds_feature.protein_id
            name = cds_feature.gene
            for annot in (entry.genes.annotations if entry.genes and entry.genes.annotations else []):
                if locus_tag and annot.id == locus_tag:
                    pass
                elif protein_id and annot.id == protein_id:
                    pass
                elif name and annot.name == name:
                    pass
                else:
                    continue
                if annot.product:
                    cds_feature.product = annot.product

    existing = set()
    for cds in record.get_cds_features_within_location(area.location, with_overlapping=False):
        for name in [cds.locus_tag, cds.protein_id, cds.gene]:
            if name is not None:
                existing.add(name)

    referenced = set()
    if entry.genes and entry.genes.annotations:
        for gene in entry.genes.annotations:
            gene_id = str(gene.id)
            referenced.add(gene_id)
            if gene_id in existing:
                existing.add(str(gene.name))
    referenced.update([str(r) for r in entry.biosynthesis.genes_referenced])

    missing = referenced.difference(existing)
    if missing:
        raise ValueError(
            f"{entry.accession} refers to missing genes: {', '.join(sorted(missing))}")

    ambiguous = referenced.intersection(record.get_renames())
    if ambiguous:
        raise ValueError(
            f"{entry.accession} uses reference name(s) referring to multiple genes : {', '.join(sorted(ambiguous))}")

    return MibigAnnotations(record.id, area, entry, cache_file, pubmed_cache_file, doi_cache_file)
