# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from antismash.common.layers import (
    Markup,
    OptionsLayer as _OptionsLayer,
)


class OptionsLayer(_OptionsLayer):
    def base_url(self) -> str:
        return "https://mibig.secondarymetabolites.org/"

    def get_name(self) -> str:
        """ Returns the ID of a record, with any extra notation included """
        name = self.seq_record.id
        if self.seq_record.has_multiple_sources():
            sources = self.seq_record.get_sources()
            name += f" (combined with {len(sources) - 1} other{'s' if len(sources) > 2 else ''})"
        return name

    def get_from_record(self) -> Markup:
        """ Returns the text to be displayed in the HTML overview page """

        return Markup(f"<strong>{self.get_name()}</strong>")

