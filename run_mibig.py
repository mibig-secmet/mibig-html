#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to generate mibig html output from mibig json file."""

import argparse
from datetime import datetime
import json
import os
from shutil import rmtree
import sys
from tempfile import NamedTemporaryFile
import traceback
from typing import List

import antismash
from antismash.common.subprocessing import execute
from mibig.converters.read.top import Everything
from mibig_taxa import TaxonCache  # pylint:disable=no-name-in-module

import mibig_html
from mibig_html import annotations
from mibig_html.common.secmet import Record


def run(commands: List[str]) -> bool:
    return execute(commands, stdout=sys.stdout, stderr=sys.stderr).successful()


def write_log(text: str, file_path: str) -> None:
    with open(file_path, "a") as o:
        o.write("[{}] {}\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), text))


def _main(json_path: str, gbk_folder: str, cache_path: str, output_folder: str,
          log_file_path: str, mibig_version: str, mibig_only: bool, pubmed_cache: str,
          doi_cache: str) -> int:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(json_path, "r") as json_file:
        raw = json.load(json_file)
    data = Everything(raw)
    mibig_acc = data.cluster.mibig_accession
    gbk_acc = data.cluster.loci.accession
    gbk_path = os.path.join(gbk_folder, "{}.gbk".format(gbk_acc))
    output_path = os.path.join(output_folder, mibig_acc)
    reusable_json_path = os.path.join(output_path, "{}.json".format(mibig_acc))

    # load/populate cache in advance, because it is needed to fetch taxonomy information
    cache = TaxonCache(cache_path)
    try:
        tax_id = int(data.cluster.ncbi_tax_id)
        taxon = cache.get_antismash_taxon(tax_id)
    except ValueError as err:
        try:
            entry = cache.get(tax_id, True)
            write_log(
                f"Outdated taxon {mibig_acc}: {tax_id} is now {entry.tax_id} ({err})", log_file_path)
        except ValueError as err:
            write_log(f"Unrecongnisable taxon {mibig_acc}: {err}", log_file_path)
        return 1

    # TODO: Properly support running on plants for MIBiG
    orig_taxon = taxon
    if taxon == "plants":
        taxon = "fungi"

    args = [
        "-v",
        "--cb-known",
        "--minlength", "1",
        "--taxon", taxon,
        "--genefinding-tool", "none",
        "--allow-long-headers",
        "--mibig-json", json_path,
        "--mibig-cache-json", cache_path,
        "--mibig-pubmed-json", pubmed_cache,
        "--mibig-doi-json", doi_cache,
        "--output-dir", output_path,
        "--output-basename", f"{mibig_acc}",
    ]
    all_modules = mibig_html.get_all_modules()
    assert mibig_html.annotations in all_modules
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    print("Generating MIBiG output for {}".format(mibig_acc))
    could_reuse = False
    operation = "generated"
    if os.path.exists(reusable_json_path):
        options = antismash.config.build_config(
            args + ["--reuse-results", reusable_json_path], parser=parser, modules=all_modules)
        try:
            mibig_html.run_mibig("", options)
            could_reuse = True
            operation = "reused"
        except Exception as err:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    if not could_reuse:
        options = antismash.config.build_config(args, parser=parser, modules=all_modules)
        if os.path.exists(output_path):
            # remove output path, proceed with caution!
            rmtree(output_path)
        try:
            mibig_html.run_mibig(gbk_path, options)
        except Exception as err:
            write_log(f"Failed to generate MIBiG page for {mibig_acc}: {err}", log_file_path)
            if os.path.exists(output_path):
                rmtree(output_path)
            raise
            return 1
    write_log(f"Successfully {operation} MIBiG page for {mibig_acc}", log_file_path)

    if mibig_only or orig_taxon == "plants":
        return 0

    print("Generating antiSMASH output for {}".format(mibig_acc))
    with open(os.path.join(output_path, "{}.json".format(mibig_acc)), "r") as result_json_txt:
        result_json = json.load(result_json_txt)
        assert len(result_json["records"]
                   ) == 1 and annotations.__name__ in result_json["records"][0]["modules"]
    prefix = mibig_acc
    region_gbk_path = os.path.join(output_path, f"{prefix}.gbk")
    output_path = os.path.join(output_path, "generated")
    reusable_as5_json_path = os.path.join(output_path, f"{prefix}.json")
    region_length = len(Record.from_genbank(region_gbk_path, taxon=taxon)[0].seq)

    if taxon in ["bacteria", "fungi"]:
        sideload_data = {
            "tool": {
                "name": "MIBiG",
                "version": mibig_version,
            },
            "records": [
                {
                    "name": mibig_acc,
                    "subregions": [
                        {
                            "start": 0,
                            "end": region_length,
                            "label": mibig_acc,
                        }
                    ]
                }
            ],
        }
        command = ["antismash"]
        args = [
            "-v",
            "--minlength", "1",
            "--html-title", "{}".format(mibig_acc),
            "--cb-known",
            "--cc-mibig",
            "--taxon", taxon,
            "--output-dir", output_path,
        ]
        reuse_as5_success = False
        run_success = False
        if os.path.exists(reusable_as5_json_path):
            reuse_as5_success = run(command + args + ["--reuse-results", reusable_as5_json_path])
            if not reuse_as5_success:
                write_log("Failed to reuse antiSMASH page for {}".format(mibig_acc), log_file_path)
                rmtree(output_path)
        if not reuse_as5_success:
            with NamedTemporaryFile() as temp:
                with open(temp.name, "w") as handle:
                    json.dump(sideload_data, handle)
                args.extend(["--sideload", temp.name])
                run_success = run(command + args + [region_gbk_path])
        if reuse_as5_success or run_success:
            write_log("Successfully generated antiSMASH page for {}".format(mibig_acc), log_file_path)
            # finally, ensure the freshly generated genbank is loadable
            try:
                antismash.common.secmet.Record.from_genbank(
                    os.path.join(output_path, f"{prefix}.gbk"), taxon=taxon)
            except antismash.common.secmet.errors.SecmetInvalidInputError as err:
                write_log(f"antiSMASH genbank for {mibig_acc} cannot be re-parsed", log_file_path)
                return 1
        else:
            if os.path.exists(output_path):
                rmtree(output_path)
            write_log("Failed to generate antiSMASH page for {}".format(mibig_acc), log_file_path)
            return 1
    elif taxon == "plants":
        write_log("Failed to generate antiSMASH page for plant BGC {}".format(mibig_acc), log_file_path)
        return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("json", type=str, help="The JSON file for the entry to run")
    parser.add_argument("genbanks", type=str, help="The directory containing genbank files")
    parser.add_argument("output", type=str, help="The directory to save results in")
    parser.add_argument("logfile", type=str, help="The path of the log file to use")
    parser.add_argument("mibig_version", type=str,
                        help="The version of mibig to display in results")
    parser.add_argument("-c", "--cache", type=str, default="tax_cache.json",
                        help="The cache file containing the mibig-taxa cache")
    parser.add_argument("-p", "--pubmed-cache", type=str, default="pubmed_cache.json",
                        help="The cache file containing the pubmed cache")
    parser.add_argument("-d", "--doi-cache", type=str, default="doi_cache.json",
                        help="The DOI cache file")
    parser.add_argument("-m", "--mibig-only", action="store_true",
                        help="Only run the MIBiG generation, skip full antiSMASH run")
    args = parser.parse_args()
    if _main(args.json, args.genbanks, args.cache, args.output, args.logfile,
             args.mibig_version, args.mibig_only, args.pubmed_cache, args.doi_cache):
        print("Errors were encountered, see log file for details")
        sys.exit(1)
    sys.exit(0)
