#   hutspot - a DNAseq variant calling pipeline
#   Copyright (C) 2017-2019, Sander Bollen, Leiden University Medical Center
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
fastqc_stats.py
~~~~~~~~~~~~~~~

:copyright: (c) 2017-2019 Sander Bollen
:copyright: (c) 2017-2019 Leiden University Medical Center
:license: AGPL-3.0
"""

import argparse
import json
from pathlib import Path
from typing import List, Union
from copy import deepcopy
from warnings import warn
import zipfile


DO_NOT_TWO_COLUMN = ["Sequence_Length_Distribution",
                     "Per_sequence_quality_scores"]


def try_numeric(val: str) -> Union[str, int, float]:
    try:
        value = int(val)
    except ValueError:
        try:
            value = float(val)
        except ValueError:
            value = val
    return value


class FastqcModule(object):
    """Helper class for Fastqc Modules"""

    def __init__(self, name):
        self.name = name
        self.header = None
        self.rows = []

    def add_header(self, headerstring: str):
        self.header = headerstring.strip().split("\t")

    def add_row(self, rowstring: str):
        self.rows.append(rowstring.strip().split("\t"))

    def row_as_dict(self, row) -> dict:
        return {k.replace(" ", "_"): try_numeric(row[i])
                for i, k in enumerate(self.header)}

    def to_dict_list(self) -> Union[List, List[dict]]:
        # The semantics herein _ONLY_ works for fastqc version 0.11.5
        # If the version in the pipeline changes, we MUST revisit this section
        if self.header is None:
            return []
        if (all([len(x) == 2 for x in self.rows])
            and len(self.header) == 2
            and self.name not in DO_NOT_TWO_COLUMN):
            # two-column data is returned as a single list
            # second column is assumed to contain all the data
            return [try_numeric(x[1]) for x in self.rows]
        return [self.row_as_dict(x) for x in self.rows]


def extract_data_txt(zip_path: Path, encoding='utf-8') -> str:
    """
    Extract text of data.txt from fastqc zip file

    Returns emtpy string when file is completely empty
    """
    if zip_path.stat().st_size == 0:
        warn(f"File {str(zip_path)} is empty")
        return ""
    with zipfile.ZipFile(zip_path) as zp:
        iflist = zp.infolist()
        data_info = next(x for x in iflist
                         if x.filename.endswith("fastqc_data.txt"))
        with zp.open(data_info.filename) as data_handle:
            return data_handle.read().decode(encoding)


def make_data_modules(data: str) -> List[FastqcModule]:
    """Make fastqc modules from data"""
    modules = []
    in_module = False
    cur_module = None
    for line in data.split("\n"):
        if not in_module and line.startswith(">>"):
            name = line.split(">>")[1].split("\t")[0].replace(" ", "_")
            cur_module = FastqcModule(name)
            in_module = True
        elif in_module and line.startswith(">>END_MODULE"):
            modules.append(deepcopy(cur_module))
            cur_module = None
            in_module = False
        elif in_module and line.startswith("#"):
            cur_module.add_header(line.split("#")[1])
        elif in_module:
            cur_module.add_row(line)
    return modules


def data_to_dict(data: str, exclusions: List[str]) -> dict:
    """Create dictionary from data"""
    modules = filter(lambda x: x.name not in exclusions,
                     make_data_modules(data))
    return {m.name: m.to_dict_list() for m in modules}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--preqc-r1",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for forward strand of "
                             "merged fastq file before pre-processing")
    parser.add_argument("--preqc-r2",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for reverse strand of "
                             "merged fastq file before pre-processing")
    parser.add_argument("--postqc-r1",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for forward strand of "
                             "fastq file after pre-processing")
    parser.add_argument("--postqc-r2",
                        type=Path,
                        required=True,
                        help="Fastqc zip for file reverse strand of "
                             "fastq file after pre-processing")
    parser.add_argument("--exclude-modules",
                        action="append",
                        default=["Per_tile_sequence_quality",
                                 "Basic_Statistics"],
                        help="Fastqc modules to exclude")
    args = parser.parse_args()

    excl = args.exclude_modules

    data_prr1 = extract_data_txt(args.preqc_r1)
    data_prr2 = extract_data_txt(args.preqc_r2)
    data_por1 = extract_data_txt(args.postqc_r1)
    data_por2 = extract_data_txt(args.postqc_r2)

    d = {
        "preqc_r1": data_to_dict(data_prr1, excl),
        "preqc_r2": data_to_dict(data_prr2, excl),
        "postqc_r1": data_to_dict(data_por1, excl),
        "postqc_r2": data_to_dict(data_por2, excl)
    }

    print(json.dumps(d))
