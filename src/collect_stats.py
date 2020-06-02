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
collect_stats.py
~~~~~~~~~~~~~~~~

:copyright: (c) 2017-2019 Sander Bollen
:copyright: (c) 2017-2019 Leiden University Medical Center
:license: AGPL-3.0
"""
import argparse
import json

from os.path import basename


def parse_json_file(path):
    with open(path) as handle:
        d = json.load(handle)
    return d


def parse_num_file(path):
    with open(path) as handle:
        line = handle.readline()
    return int(line.strip())


def determine_gender(covstat, fthresh):
    """Determine gender from a covstat json """
    cv = covstat['stats']['coverage']
    all = cv['_all']['median']

    if 'chrX' in cv:
        x = cv['chrX']['median']
    elif 'X' in cv:
        x = cv['X']['median']
    else:
        return "NA"

    if all != 0:
        rat = x/all
    else:
        return "NA"

    if rat >= fthresh:
        return "female"
    return "male"


def main(args):
    mpnum = parse_num_file(args.mapped_num)
    mpbnum = parse_num_file(args.mapped_basenum)
    unum = parse_num_file(args.unique_num)
    ubnum = parse_num_file(args.usable_basenum)
    cutadapt = parse_json_file(args.cutadapt)
    

    d = {
        "sample_name": args.sample_name,
        "preqc_reads": cutadapt["preqc_reads"],
        "preqc_bases": cutadapt["preqc_bases"],
        "postqc_reads": cutadapt["postqc_reads"],
        "postqc_bases": cutadapt["postqc_bases"],
        "n_mapped_reads": mpnum,
        "n_mapped_bases": mpbnum,
        "n_usable_reads": unum,
        "n_usable_bases": ubnum
    }

    # "." is used to pass an 'empty' file from snakemake, since all snakemake
    # inputs must be files or folders which exist
    if args.covstats != ".":
        # Read the json file
        covstats = parse_json_file(args.covstats)
        # Format the coverage data and determine the gender
        cov_data = {
            "name": basename(args.covstats),
            "gender": determine_gender(covstats, args.female_threshold),
            "covstats": covstats
        }
        # Add the coverage data
        d["covstats"] = cov_data


    print(json.dumps(d))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample-name", required=True,
                        help="Sample name")
    parser.add_argument("--mapped-num", required=True,
                        help="Mapped num file")
    parser.add_argument("--mapped-basenum", required=True,
                        help="Mapped basenum file")
    parser.add_argument("--unique-num", required=True,
                        help="Unique num file")
    parser.add_argument("--usable-basenum", required=True,
                        help="Usable basenum")
    parser.add_argument("--female-threshold", default=0.6, 
                        help="Female threshold of X/all cov")
    parser.add_argument("--cutadapt", required=True,
                        help="Cutadapt summary output")
    parser.add_argument("covstats",
                        help="Coverage statistics")

    args = parser.parse_args()
    main(args)
