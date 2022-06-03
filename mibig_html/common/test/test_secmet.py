# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.errors import SecmetInvalidInputError
from Bio.SeqRecord import SeqRecord

from mibig_html.common.secmet import (
    ASRecord,
    Record,
)

class TestRecord(unittest.TestCase):
    def test_molecule_type_cleanup(self):
        bio = SeqRecord(seq="ATGC")
        for mol in ["RNA", "mRNA"]:
            bio.annotations["molecule_type"] = mol
            with self.assertRaisesRegex(SecmetInvalidInputError, f"{mol} records are not supported"):
                ASRecord.from_biopython(bio, taxon="bacteria")
            record = Record.from_biopython(bio, taxon="bacteria")
            assert record.annotations["molecule_type"] == "DNA"
