# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from mibig_html.html.generator import (
    convert_categories,
)

class TestCategories(unittest.TestCase):
    def test_case(self):
        assert convert_categories(["terpene"]) == ["terpene"]
        assert convert_categories(["Terpene"]) == ["terpene"]

    def test_alternatives(self):
        assert convert_categories(["NRP"]) == ["NRPS"]
        assert convert_categories(["Polyketide"]) == ["PKS"]

    def test_combo(self):
        categories = ["NRP", "Saccharide", "terpene"]
        result = convert_categories(categories)
        assert result == ["NRPS", "saccharide", "terpene"]

    def test_unknown(self):
        for bad in ["1", "invalid category name"]:
            with self.assertRaisesRegex(ValueError, "cannot translate MIBiG class"):
                convert_categories([bad])
