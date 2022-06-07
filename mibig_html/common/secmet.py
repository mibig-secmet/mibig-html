# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from typing import Any, Dict, List, Tuple, Type, TypeVar, Union

from antismash.common.secmet.record import (
    CDSFeature,
    location_bridges_origin,
    Record as ASRecord,
    Seq,
    SeqRecord,
)

T = TypeVar("T", bound="Record")


def _get_biopython_cds_name(feature: SeqFeature) -> str:
    name = ""
    for qual in ["locus_tag", "gene", "protein_id"]:
        name = feature.qualifiers.get(qual, [""])[0]
        if name:
            break
    return name


class Record(ASRecord):
    def __init__(self, seq: Union[Seq, str] = "", transl_table: int = 1, **kwargs: Any) -> None:
        super().__init__(seq, transl_table=transl_table, **kwargs)

        self._altered_from_input: List[str] = []
        self._deduplicated_cds_names: Dict[str, List[str]] = {}

    def __getattr__(self, attr: str) -> Any:
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name", "annotations", "dbxrefs"]:
            return getattr(self._record, attr)
        if attr in self.__slots__:  # changed here, as ASRecord referred by name to itself
            return self.__getattribute__(attr)
        raise AttributeError("Record has no attribute '%s'" % attr)

    def add_alteration(self, description: str) -> None:
        """ Adds a description of a fundamental change to input values to the record """
        assert description
        self._altered_from_input.append(description)

    def get_alterations(self) -> Tuple[str, ...]:
        """ Returns any alterations made from the original input,
            not including any made by loading with biopython
        """
        return tuple(self._altered_from_input)

    def add_cds_feature(self, cds: CDSFeature) -> None:
        try:
            super().add_cds_feature(cds)
            return
        except ValueError as err:
            if "same name for mapping" not in str(err):
                raise
        # the remaining code here is only reached if adding the CDS raised an error,
        # but it was a duplicate CDS error
        original_name = cds.get_name()
        if original_name not in self._deduplicated_cds_names:
            self._deduplicated_cds_names[original_name] = []
        new_name = "%s_rename%d" % (original_name, len(self._deduplicated_cds_names[original_name]) + 1)
        self._deduplicated_cds_names[original_name].append(new_name)
        if cds.locus_tag:
            cds.locus_tag = new_name
        elif cds.gene:
            cds.gene = new_name
        elif cds.protein_id:
            cds.protein_id = new_name
        assert new_name not in self._cds_by_name
        self.add_cds_feature(cds)
        self._altered_from_input.append("CDS with name %s renamed to %s to avoid duplicates" % (
                                        original_name, new_name))

    @classmethod
    def from_biopython(cls: Type[T], seq_record: SeqRecord, taxon: str) -> T:
        # handle some entry reference records being mislabeled as RNA (e.g. BGCs 488, 720, 1166)
        molecule_type = seq_record.annotations.get("molecule_type", "DNA")
        if not molecule_type.upper().endswith("DNA"):
            seq_record.annotations["molecule_type"] = "DNA"

        can_be_circular = taxon == "bacteria"
        names = set()
        for feature in seq_record.features:
            if can_be_circular and location_bridges_origin(feature.location, allow_reversing=False):
                names.add(_get_biopython_cds_name(feature))
        record = super().from_biopython(seq_record, taxon)
        for name in sorted(names):
            record.add_alteration(f"{name} crossed the origin and was split into two features")
        return record

    @classmethod
    def from_antismash_record(cls, existing: ASRecord) -> "Record":
        new = cls(existing.seq)
        for key in existing.__slots__:
            try:
                setattr(new, key, getattr(existing, key))
            except AttributeError:
                continue
        return new
