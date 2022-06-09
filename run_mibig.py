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
from typing import List

import antismash
from antismash.common.subprocessing import execute
from mibig.converters.read.top import Everything

import mibig_html
from mibig_html import annotations
from mibig_html.annotations.mibig import load_cached_information
from mibig_html.common.secmet import Record


def run(commands: List[str]) -> bool:
    return execute(commands, stdout=sys.stdout, stderr=sys.stderr).successful()


def write_log(text: str, file_path: str) -> None:
    with open(file_path, "a") as o:
        o.write("[{}] {}\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), text))


def _main(json_path: str, gbk_folder: str, cache_folder: str, output_folder: str,
          log_file_path: str, mibig_version: str) -> int:
    for path in [cache_folder, output_folder]:
        if not os.path.exists(path):
            os.makedirs(path)

    with open(json_path, "r") as json_file:
        raw = json.load(json_file)
    data = Everything(raw)
    mibig_acc = data.cluster.mibig_accession
    gbk_acc = data.cluster.loci.accession
    gbk_path = os.path.join(gbk_folder, "{}.gbk".format(gbk_acc))
    cache_json_path = os.path.join(cache_folder, "{}.cache.json".format(mibig_acc))
    output_path = os.path.join(output_folder, mibig_acc)
    reusable_json_path = os.path.join(output_path, "{}.1.json".format(mibig_acc))

    # load/populate cache in advance, because it is needed to fetch taxonomy information
    cached = load_cached_information(data, cache_json_path, True)
    taxonomy = [tax_obj["name"] for tax_obj in cached["taxonomy"][data.cluster.ncbi_tax_id]]

    if "Bacteria" in taxonomy:
        taxon = "bacteria"
    elif "Fungi" in taxonomy:
        taxon = "fungi"
    elif "Viridiplantae" in taxonomy:
        taxon = "plants"
    else:
        write_log("Unrecognizable taxons {} ({})".format(mibig_acc, ":".join(taxonomy)), log_file_path)
        return 1

    args = [
        "-v",
        "--cb-known",
        "--minlength", "1",
        "--taxon", taxon,
        "--genefinding-tool", "none",
        "--allow-long-headers",
        "--mibig-json", json_path,
        "--mibig-cache-json", cache_json_path,
        "--output-dir", output_path,
        "--output-basename", f"{mibig_acc}.1",
    ]
    all_modules = mibig_html.get_all_modules()
    assert mibig_html.annotations in all_modules
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    print("Generating MIBiG output for {}".format(mibig_acc))
    success = False
    if os.path.exists(reusable_json_path):
        options = antismash.config.build_config(args + ["--reuse-results", reusable_json_path], parser=parser, modules=all_modules)
        try:
            mibig_html.run_mibig("", options)
            success = True
            write_log("Successfully reused JSON file {}".format(reusable_json_path), log_file_path)
        except Exception as err:
            write_log(f"Failed to reuse JSON file {reusable_json_path}: {err}", log_file_path)
            success = False

    if not success:
        options = antismash.config.build_config(args, parser=parser, modules=all_modules)
        if os.path.exists(output_path):
            # remove output path, proceed with caution!
            rmtree(output_path)
            write_log("Removed {}".format(output_path), log_file_path)
        try:
            mibig_html.run_mibig(gbk_path, options)
        except Exception as err:
            write_log(f"Failed to generate MIBiG page for {mibig_acc}: {err}", log_file_path)
            if os.path.exists(output_path):
                rmtree(output_path)
            raise
            return 1
    write_log("Successfully generated MIBiG page for {}".format(mibig_acc), log_file_path)

    print("Generating antiSMASH output for {}".format(mibig_acc))
    with open(os.path.join(output_path, "{}.1.json".format(mibig_acc)), "r") as result_json_txt:
        result_json = json.load(result_json_txt)
        assert len(result_json["records"]) == 1 and annotations.__name__ in result_json["records"][0]["modules"]
    prefix = f"{mibig_acc}.1"
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
                    "name": f"{mibig_acc}.1",
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
    parser.add_argument("cache", type=str, help="The cache directory for the MIBiG module")
    parser.add_argument("output", type=str, help="The directory to save results in")
    parser.add_argument("logfile", type=str, help="The path of the log file to use")
    parser.add_argument("mibig_version", type=str, help="The version of mibig to display in results")

    args = parser.parse_args()
    sys.exit(_main(args.json, args.genbanks, args.cache, args.output, args.logfile, args.mibig_version))
