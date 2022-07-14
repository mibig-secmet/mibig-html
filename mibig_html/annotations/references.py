# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Reference storage, look up, and caching """

import json
import logging
import os
import re
import requests
from typing import Any, Dict, List, Optional, Set, Type


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


class DoiEntry(ReferenceEntry):
    def __init__(self, title: str, authors: List[str], year: str, journal: str, identifier: str) -> None:
        super().__init__(title, authors, year, journal, identifier)

    @property
    def info(self) -> str:
        return f"{super().info} DOI:{self.identifier}"


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


class DoiCache(ReferenceCache):
    def __init__(self, cache_file: str) -> None:
        super().__init__(cache_file, DoiEntry, "DOI")

    def add_entry(self, entry: DoiEntry) -> None:
        self.mappings[entry.identifier] = entry
        self._new_updates.add(entry.identifier)

    def get(self, identifier: str) -> ReferenceEntry:
        if identifier not in self.mappings:
            self.resolve(identifier)
        return self.mappings[identifier]

    def resolve(self, identifier: str) -> ReferenceEntry:
        accepts = {
            "json": "application/vnd.citationstyles.csl+json; charset=utf-8",
            "bibtex": "application/x-bibtex; charset=utf-8",
        }

        def simple_request(doi: str, accept: str) -> str:
            url = f"https://dx.doi.org/{doi}"
            req = requests.get(url, headers={"Accept": accept})
            if req.status_code == 404:
                raise ValueError(f"Invalid DOI: {url}")
            elif req.status_code == 204:  # valid DOI, but no json friendly version
                return ""
            elif req.status_code != 200:
                raise ValueError(f"Error fetching reference data {url}: {req.text}")
            assert req.status_code == 200
            return req.text

        def request_json(doi: str) -> Optional[DoiEntry]:
            def simplify_author(author: Dict[str, Any]) -> str:
                if "given" in author:
                    return f"{author['family']}, {author['given'][0]}"
                else:
                    return author["family"]

            raw = simple_request(doi, accepts["json"])
            if not raw:
                return None

            values = json.loads(raw)
            authors = [simplify_author(author) for author in values["author"]]
            year = values["published"]["date-parts"][0]
            if isinstance(year, list):
                year = year[0]
            assert isinstance(year, int)
            journal_titles = values.get("container-title", "")
            if isinstance(journal_titles, list):
                if journal_titles:
                    assert False, journal_titles
                journal_titles = ""
            journal = journal_titles.replace("&amp;", "&")
            # if the journal is missing, check if it's a preprint
            if not journal:
                if values.get("subtype") == "preprint":
                    journal = "preprint"
            # and if still not set, it might be a bioRxiv preprint
            if not journal:
                institutions = values.get("institution", [])
                for inst in institutions:
                    if inst.get("name") == "bioRxiv":
                        journal = "preprint"
                        break
            assert journal

            # strip HTML with odd spacing from the title
            title = values["title"]
            title = re.sub("\n +<[^>]*>", " ", title)  # opening tags
            title = re.sub("<[^>]*>\n +", " ", title)  # closing tags
            return DoiEntry(title, authors, str(year), journal, doi)

        def request_bibtex(doi: str) -> Optional[DoiEntry]:
            raw = simple_request(doi, accepts["bibtex"])
            if not raw:
                return None
            try:
                lines = (line for line in raw.splitlines() if "=" in line)
                pairs = {}
                for line in lines:
                    k, v = line.split(" = ", 1)
                    k = k.strip()
                    pairs[k] = v.strip().strip("{},")

                authors = pairs["author"].split(" and ")
                # is it the brief form reversed?
                if len(authors[0].split()[1].strip(".")) == 1:
                    authors = [author.strip(".") for author in authors]
                else:
                    reformatted_authors = []
                    for author in authors:
                        parts = author.split()
                        first = parts[0][0]
                        last = parts[1]
                        reformatted_authors.append(f"{last} {first}")

                return DoiEntry(pairs["title"], authors, pairs["year"], pairs["journal"], doi)
            except Exception as err:
                print(raw)
                raise err

        entry = request_json(identifier)
        if entry is None:
            entry = request_bibtex(identifier)
        if entry is None:
            raise ValueError(f"exhausted DOI metadata options for: {identifier}")
        self.add_entry(entry)
        return entry
