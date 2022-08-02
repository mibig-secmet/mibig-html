#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import glob
import json
import os
import sys
from typing import List

from mibig_html.annotations.references import DoiCache, DoiEntry

SPECIAL = {
    "10.12211/2096-8280.2021-024": {  # times out for anything but HTML
        "title": "Genome mining for novel natural products in Sorangium cellulosum So0157-2 by heterologous expression",
        "authors": ["Zhou, H" "Shen, Q", "Chen, H", "Wang, Z", "Li, Y", "Zhang, Y", "Bian, X"],
        "year": "2021",
        "journal": "Synthetic Biology Journal",
        "identifier": "10.12211/2096-8280.2021-024",
    },
}


def fetch_all(cache_file: str, files: List[str]) -> None:
    doi_cache = DoiCache(cache_file)
    for filename in files:
        with open(filename) as handle:
            data = json.load(handle)
        if "publications" not in data["cluster"]:
            continue
        for ref in data["cluster"]["publications"]:
            if not ref.startswith("doi:"):
                continue
            doi = ref.split(":", 1)[1]
            if doi in SPECIAL:
                doi_cache.add_entry(DoiEntry.from_json(SPECIAL[doi]))
            else:
                try:
                    doi_cache.get(doi)
                    doi_cache.save()
                except ValueError as err:
                    print("failed to import DOIs from", filename, err)
                    raise
    doi_cache.save()


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} input_dir cache_file")
        sys.exit(1)
    cache = "doi_cache.json"
    if len(sys.argv) == 3:
        cache = sys.argv[2]
    fetch_all(cache, glob.glob(f"{sys.argv[1]}/*.json"))
