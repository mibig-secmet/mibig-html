# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Reference storage, look up, and caching """

import json
import logging
import os
from typing import Any, Dict, List, Set, Type


class ReferenceEntry:
    def __init__(self, title: str, authors: List[str], year: str, journal: str, identifier: str) -> None:
        # make title styling consistent
        if not title.endswith("."):
            title += "."
        self.title = title
        self.authors = authors
        self.year = year
        self.journal = journal
        self.identifier = identifier

    @property
    def info(self) -> str:
        extras = " et al." if len(self.authors) > 1 else ""
        return f"{self.authors[0]}{extras}, {self.journal} ({self.year})"

    def to_json(self) -> Dict[str, Any]:
        return {
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
            "journal": self.journal,
            "identifier": self.identifier,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "ReferenceEntry":
        return cls(data["title"], data["authors"], data["year"], data["journal"], data["identifier"])


class PubmedEntry(ReferenceEntry):
    @property
    def info(self) -> str:
        return f"{super().info} PMID:{self.identifier}"

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "ReferenceEntry":
        return cls(data["title"], data["authors"], data["year"], data["journal"], data.get("pmid", data.get("identifier")))

    def to_json(self) -> Dict[str, Any]:
        result = super().to_json()
        result["pmid"] = result.pop("identifier")
        return result


class ReferenceCache:
    def __init__(self, cache_file: str, entry_class: Type[ReferenceEntry], label: str) -> None:
        self.mappings: Dict[str, ReferenceEntry] = {}
        self.cache_file = cache_file
        self.entry_class = entry_class
        self.label = label
        if cache_file and os.path.exists(cache_file):
            with open(cache_file, 'r') as handle:
                entries = json.load(handle)
            for entry_id, entry_values in entries.items():
                self.mappings[entry_id] = entry_class.from_json(entry_values)
        self._new_updates: Set[str] = set()

    def add(self, title: str, authors: List[str], year: str, journal: str, identifier: str) -> None:
        self.mappings[identifier] = self.entry_class(title, authors, year, journal, identifier)
        self._new_updates.add(identifier)

    def get(self, identifier: str) -> ReferenceEntry:
        return self.mappings[identifier]

    def get_missing(self, want: List[str]) -> List[str]:
        have = self.mappings.keys()
        return sorted(list(want - have))

    def save(self) -> None:
        logging.debug("Updating %s cache file with %d new entries: %s",
                      self.label, len(self._new_updates), self.cache_file)
        with open(self.cache_file, "w", encoding="utf_8") as handle:
            json.dump({key: val.to_json() for key, val in self.mappings.items()}, handle)

    def __del__(self) -> None:
        # if the object was never fully instantiated, don't do anything with it
        if not hasattr(self, "_new_updates"):  # ideally the last thing set in __init__
            return
        # then ensure there's something to save and somewhere to save to
        if not self._new_updates or not self.cache_file:
            return
        # then ensure that python hasn't already cleaned up some core functionality
        try:
            _ = (open, logging, json)
        except NameError:
            return
        # finally, save the updates
        self.save()


class PubmedCache(ReferenceCache):
    def __init__(self, cache_file: str) -> None:
        super().__init__(cache_file, PubmedEntry, "PubMed")
