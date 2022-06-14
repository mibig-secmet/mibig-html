# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" MiBIG specific sideloading """

import json
from os import path
from typing import Any, Dict, List, Optional
import logging

from Bio import Entrez
from xml.dom import minidom
import time

from antismash.common.errors import AntismashInputError
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, SubRegion, Record
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    location_contains_other,
)

from urllib.error import HTTPError
from mibig.converters.read.top import Everything


class MibigAnnotations(DetectionResults):
    def __init__(self, record_id: str, area: SubRegion, data: Everything, cache_file: str) -> None:
        super().__init__(record_id)
        self.data = data  # holds the original annotation json data
        # save calculated loci (relative to record), not annotated ones
        self.record_id = record_id
        self.area = area
        # fetch/update cached information (for taxonomy, etc.)
        cached = load_cached_information(data, cache_file)
        # save extra information from cache
        self.taxonomy = cached["taxonomy"][data.cluster.ncbi_tax_id]

    def get_predicted_subregions(self) -> List[SubRegion]:
        return [self.area]

    def to_json(self) -> Dict[str, Any]:
        # save only information critical for deciding reusability
        loci = self.data.cluster.loci
        annotations = []
        for annot in (self.data.cluster.genes.annotations if self.data.cluster.genes else []):
            annotations.append(annot.to_json())
        extra_genes = []
        for extra_gene in (self.data.cluster.genes.extra_genes if self.data.cluster.genes else []):
            extra_genes.append(extra_gene.to_json())
        return {
            "record_id": self.record_id,
            "genbank_accession": loci.accession,
            "coords": (loci.start or -1, loci.end or -1),
            "gene_annotations": annotations,
            "extra_genes": extra_genes,
        }

    @staticmethod
    def from_json(prev: Dict[str, Any], record: Record, annotations_file: str,
                  cache_file: str) -> Optional["MibigAnnotations"]:
        with open(annotations_file) as handle:
            raw = json.load(handle)
            data = Everything(raw)

        # compare old vs new annotation, decide if we can 'reuse'
        can_reuse = True
        loci = data.cluster.loci
        gene_annotations = data.cluster.genes.annotations if data.cluster.genes else []
        extra_genes = data.cluster.genes.extra_genes if data.cluster.genes else []
        if loci.accession != prev["genbank_accession"]:
            logging.debug("Previous result's genbank_accession is not the same as the new one")
            can_reuse = False
        elif record.id != prev["record_id"]:
            logging.debug("Previous result's record_id is not the same as the new one")
            can_reuse = False
        elif (loci.start or -1) != prev["coords"][0] or (loci.end or -1) != prev["coords"][1]:
            logging.debug("Previous result's start/end coordinate is not the same as the new one")
            can_reuse = False
        elif len(gene_annotations) != len(prev["gene_annotations"]):
            logging.debug("Gene annotations have changed")
            can_reuse = False
        elif len(extra_genes) != len(prev["extra_genes"]):
            logging.debug("Additional genes have changed")
            can_reuse = False

        # if we can't reuse, stop running antismash, because CDS annotations won't be correct
        if can_reuse:
            product = ", ".join(data.cluster.biosynthetic_class)
            loci_region = FeatureLocation(
                loci.start - 1 if loci.start else 0,
                loci.end or len(record.seq)
            )
            area = SubRegion(loci_region, tool="mibig", label=product)
            return MibigAnnotations(record.id, area, data, cache_file)
        else:
            logging.error("Can't reuse MIBiG annotation.")
            raise AntismashInputError("Genbank record or gene annotations are updated, can't reuse result")


def mibig_loader(annotations_file: str, cache_file: str, record: Record) -> MibigAnnotations:
    """This method will be called only when not reusing data"""
    with open(annotations_file) as handle:
        raw = json.load(handle)
        data = Everything(raw)

    product = ", ".join(data.cluster.biosynthetic_class)
    loci = data.cluster.loci
    loci_region = FeatureLocation(
        loci.start - 1 if loci.start else 0,
        loci.end or len(record.seq)
    )
    area = SubRegion(loci_region, tool="mibig", label=product)

    # re-annotate CDSes in the Record
    if data.cluster.genes:
        # extra genes
        for gene in data.cluster.genes.extra_genes:
            if gene.id and gene.location:
                exons = [FeatureLocation(exon.start - 1, exon.end, strand=gene.location.strand) for exon in gene.location.exons]
                location = CompoundLocation(exons) if len(exons) > 1 else exons[0]
                if not location_contains_other(area.location, location):
                    raise ValueError(f"additional gene {gene.id} lies outside cluster")
                translation = gene.translation or record.get_aa_translation_from_location(location)
                cds_feature = CDSFeature(location=location, locus_tag=gene.id, translation=translation)
                record.add_cds_feature(cds_feature, auto_deduplicate=False)
                record.add_alteration(f"{cds_feature.get_name()} was added")
        # re-annotation
        for cds_feature in record.get_cds_features_within_location(area.location):
            locus_tag = cds_feature.locus_tag
            protein_id = cds_feature.protein_id
            name = cds_feature.gene
            for annot in (data.cluster.genes.annotations if data.cluster.genes else []):
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
    for cds in record.get_cds_features():
        for name in [cds.locus_tag, cds.protein_id, cds.gene]:
            if name is not None:
                existing.add(name)

    referenced = set()
    if data.cluster.genes and data.cluster.genes.annotations:
        for gene in data.cluster.genes.annotations:
            referenced.add(gene.id)
    if data.cluster.nrp:
        for nrps_gene in data.cluster.nrp.nrps_genes:
            referenced.add(nrps_gene.gene_id)
        for thio in data.cluster.nrp.thioesterases:
            referenced.add(thio.gene)
    if data.cluster.polyketide:
        for synthase in data.cluster.polyketide.synthases:
            referenced.update(set(synthase.genes))
            for module in synthase.modules:
                referenced.update(set(module.genes))
    if data.cluster.saccharide:
        for transferase in data.cluster.saccharide.glycosyltransferases:
            referenced.add(transferase.gene_id)
    for compound in data.cluster.compounds:
        for moiety in compound.chem_moieties:
            if moiety.subcluster:
                referenced.update(set(moiety.subcluster))

    missing = referenced.difference(existing)
    if missing:
        raise ValueError(f"{data.cluster.mibig_accession} refers to missing genes: {', '.join(sorted(missing))}")

    return MibigAnnotations(record.id, area, data, cache_file)


def load_cached_information(annotations: Everything, cache_json_path: str, update: bool = True) -> Dict[str, Any]:
    """"""
    if len(cache_json_path) > 0 and (path.exists(cache_json_path)):
        with open(cache_json_path) as handle:
            cached = json.load(handle)
    else:
        cached = {}

    ncbi_email = "mibig@secondarymetabolites.org"

    assert isinstance(cached, dict)

    # fetch taxonomy information
    if "taxonomy" not in cached:
        cached["taxonomy"] = {}
    ncbi_tax_id = annotations.cluster.ncbi_tax_id
    if ncbi_tax_id not in cached["taxonomy"]:
        cached["taxonomy"][ncbi_tax_id] = get_ncbi_taxonomy(ncbi_tax_id, ncbi_email)

    # fetch BibTex for publications
    # ....

    if update:
        # update cache file
        save_cached_information(cached, cache_json_path)

    return cached


def save_cached_information(cached: Dict[str, Any], cache_json_path: str) -> None:
    """ Saves taxonomy information for later use as JSON

        Arguments:
            cached: the information to save, in a JSON-friendly format
            cache_json_path: the file path to save the cached information

        Returns:
            None
    """
    with open(cache_json_path, "w") as handle:
        handle.write(json.dumps(cached, indent=4, separators=(',', ': '), sort_keys=True))


def get_ncbi_taxonomy(tax_id: str, email: str) -> List[Dict[str, Any]]:
    """fetch taxonomy information from ncbi_tax_id"""
    taxonomy = []
    Entrez.email = email
    num_try = 1
    while num_try < 6:
        try:
            logging.debug("Fetching taxonomy information from NCBI for tax_id:{}...".format(tax_id))
            dom = minidom.parse(Entrez.efetch(db="taxonomy", id=tax_id))
            for dom_taxon in dom.getElementsByTagName('Taxon'):
                taxid = dom_taxon.getElementsByTagName("TaxId")[0].firstChild.nodeValue
                name = dom_taxon.getElementsByTagName("ScientificName")[0].firstChild.nodeValue
                rank = dom_taxon.getElementsByTagName("Rank")[0].firstChild.nodeValue
                taxonomy.append({"name": name, "taxid": taxid, "rank": rank})
            break
        except HTTPError:
            logging.error("Failed to query NCBI taxonomy database, retrying...")
            pass
        num_try += 1
        time.sleep(5)
    if len(taxonomy) > 1:  # shuffle species to the end of the list
        taxonomy.append(taxonomy.pop(0))

    return taxonomy
