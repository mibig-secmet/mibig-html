#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import glob
import json
import os
import sys

from eutils import Client

from mibig_html.annotations.references import PubmedCache, PubmedEntry

def extract_pmids(files: list[str]) -> list[str]:
    """Extract all unique pmids from mibig jspn files"""
    pmids: set[str] = set()
    for filename in files:
        with open(filename) as handle:
            data = json.load(handle)
        if "publications" not in data["cluster"]:
            continue
        for ref in data["cluster"]["publications"]:
            if not ref.startswith("pubmed:"):
                continue
            pmid = ref.split(":", 1)[1]
            if pmid == "0":
                continue
            pmids.add(pmid)
    return sorted(list(pmids))


def fetch_all(cache_file: str, pmids: list[str]) -> None:
    pubmed_cache = PubmedCache(cache_file)
    client = Client(api_key=os.environ.get("NCBI_API_KEY", None))

    missing = pubmed_cache.get_missing(pmids)
    while missing:
        articles = client.efetch(db="pubmed", id=missing)
        for article in articles:
            pubmed_cache.add(article.title, article.authors,
                                article.year, article.jrnl, article.pmid)
        missing = pubmed_cache.get_missing(pmids)
        
    pubmed_cache.save()


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} input_dir cache_file")
        sys.exit(1)
    cache = "pubmed_cache.json"
    if len(sys.argv) == 3:
        cache = sys.argv[2]
    pmids = extract_pmids(glob.glob(f"{sys.argv[1]}/*.json"))
    fetch_all(cache, pmids)
