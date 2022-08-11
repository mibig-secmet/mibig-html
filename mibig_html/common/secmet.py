# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple, Type, TypeVar, Union

from antismash.common.secmet.record import (
    CDSFeature,
    location_bridges_origin,
    Record as ASRecord,
    SecmetInvalidInputError,
    Seq,
    SeqFeature,
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
        self._deduplicated_cds_names: Dict[str, List[str]] = defaultdict(list)
        self._alternative_names: Dict[str, Set[str]] = defaultdict(set)

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

    def get_renames(self) -> Set[str]:
        """ Returns a set of original gene names that were altered """
        return set(self._deduplicated_cds_names)

    def add_cds_feature(self, cds_feature: CDSFeature, auto_deduplicate: bool = True) -> None:
        def add_alternative_names() -> None:
            real_name = cds_feature.get_name()
            for alternative in [cds_feature.locus_tag, cds_feature.gene, cds_feature.protein_id]:
                if alternative:
                    self._alternative_names[alternative].add(real_name)
            assert real_name in self._alternative_names

        if not auto_deduplicate:
            add_alternative_names()
            super().add_cds_feature(cds_feature)
            return

        original_name = cds_feature.get_name()
        duplicate_name = original_name in self._cds_by_name
        duplicate_location = str(cds_feature.location) in self._cds_by_location

        if not duplicate_name and not duplicate_location:
            add_alternative_names()
            super().add_cds_feature(cds_feature)
            return

        if duplicate_location and duplicate_name:
            self.add_alteration(f"removed an exact duplicate of CDS feature {cds_feature.get_name()}")
            return

        if duplicate_location:
            existing = self._cds_by_location[str(cds_feature.location)]
            self.add_alteration(f"removed {cds_feature.get_name()} as a duplicate of {existing.get_name()}")
            return

        # which leaves only names matching, so a rename is appropriate
        count = len(self._deduplicated_cds_names[original_name]) + 1
        new_name = f"{original_name}_rename{count}"

        # update the original name field with the new name
        if cds_feature.locus_tag == original_name:
            cds_feature.locus_tag = new_name
        elif cds_feature.gene == original_name:
            cds_feature.gene = new_name
        elif cds_feature.protein_id == original_name:
            cds_feature.protein_id = new_name

        add_alternative_names()

        # then add the modified feature and log the alteration
        super().add_cds_feature(cds_feature)
        self._deduplicated_cds_names[original_name].append(new_name)
        self.add_alteration(f"renamed CDS with name {original_name} at {cds_feature.location} to {new_name} to avoid duplicates")

    def get_real_cds_name(self, name: str) -> str:
        """ Gets the unique name for a CDS, if possible, from one of the names
            it was annotated with. If collisions would occur, an error is raised.
        """
        results = self._alternative_names.get(name, set())
        if not results:
            try:
                if self.get_cds_by_name(name):
                    return name
            except KeyError as err:
                raise ValueError(f"unknown CDS/gene name: {name}") from err
        if len(results) > 1:
            raise ValueError(f"multiple features map to alternative CDS/gene name {name}")
        return list(results)[0]

    def add_biopython_feature(self, feature: SeqFeature) -> None:
        try:
            super().add_biopython_feature(feature)
        except SecmetInvalidInputError as err:
            if "translation longer than location allows" not in str(err):
                raise
            # regenerate the translation
            feature.qualifiers.pop("translation")
            super().add_biopython_feature(feature)
            name = _get_biopython_cds_name(feature)
            self.add_alteration(f"the translation of {name} was too long and was regenerated")

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
        assert isinstance(record, cls)
        return record
