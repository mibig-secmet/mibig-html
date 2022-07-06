# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""HTML output format module

"""

import glob
import logging
import os
import shutil
from typing import Dict, List, Optional

import scss

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.custom_typing import AntismashModule
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from mibig_html.html.generator import generate_webpage

NAME = "mibig_html"
SHORT_DESCRIPTION = "HTML output for mibig-mode"


def get_arguments() -> ModuleArgs:
    """ Builds the arguments for the HMTL output module """
    return ModuleArgs("Output options", "html")


def prepare_data(_logging_only: bool = False) -> List[str]:
    """ Rebuild any dynamically buildable data """
    with path.changed_directory(path.get_full_path(__file__, "css")):
        built_file = os.path.abspath("mibig.css")

        if path.is_outdated(built_file, glob.glob("*.scss")):
            logging.info("CSS files out of date, rebuilding")

            result = scss.Compiler(output_style="expanded").compile("mibig.scss")
            assert result
            with open("mibig.css", "w") as out:
                out.write(result)
    return []


def check_prereqs(_options: ConfigType) -> List[str]:
    """ Check prerequisites """
    return prepare_data()


def check_options(_options: ConfigType) -> List[str]:
    """ Check options, but none to check """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Is the HMTL module enabled (currently always enabled) """
    return True


def write(records: List[Record], results: List[Dict[str, ModuleResults]],
          options: ConfigType, all_modules: List[AntismashModule]) -> None:
    """ Writes all results to a webpage, where applicable. Writes to options.output_dir

        Arguments:
            records: the list of Records for which results exist
            results: a list of dictionaries containing all module results for records
            options: antismash config object

        Returns:
            None
    """
    output_dir = options.output_dir

    copy_template_dir('css', output_dir, pattern="mibig.css")
    copy_template_dir('js', output_dir)
    copy_template_dir('images', output_dir)

    assert len(records) == 1
    assert len(results) == 1
    generate_webpage(records[0], results[0], options, all_modules)


def copy_template_dir(template: str, output_dir: str, pattern: Optional[str] = None) -> None:
    """ Copy files from a template directory to the output directory, removes
        any existing directory first. If pattern is supplied, only files within
        the template directory that match the template will be copied.

        Arguments:
            template: the source directory
            output_dir: the target directory
            pattern: a pattern to restrict to, if given

        Returns:
            None
    """
    target_dir = os.path.join(output_dir, template)
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    if pattern:
        os.makedirs(target_dir)
        for filename in glob.glob(path.get_full_path(__file__, template, pattern)):
            if os.path.isdir(filename):
                shutil.copytree(filename, target_dir)
            else:
                shutil.copy2(filename, target_dir)
    else:
        shutil.copytree(path.get_full_path(__file__, template), target_dir)
