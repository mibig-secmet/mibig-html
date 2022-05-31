# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Overrides some antiSMASH HTML rendering functions """

from typing import Any

from antismash.common.html_renderer import (
    collapser_start,
    collapser_end,
    FileTemplate as _FileTemplate,
    Markup,
    switch,
)

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


class FileTemplate(_FileTemplate):  # pylint: disable=too-few-public-methods
    """ A template renderer for file templates """
    def render(self, **kwargs: Any) -> Markup:
        """ Renders the template HTML, providing any given arguments to the renderer """
        if not self.template:
            raise ValueError("attempting to render without a template")

        defaults = {
            "collapser_start": collapser_start,
            "collapser_end": collapser_end,
            "help_tooltip": help_tooltip,
            "switch": switch,
        }
        defaults.update(kwargs)
        return Markup(self.template.render(**defaults))
