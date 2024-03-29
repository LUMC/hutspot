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
    cutadapt = parse_json_file(args.cutadapt)

    d = {
        "sample_name": args.sample_name,
        "preqc_reads": cutadapt["preqc_reads"],
        "preqc_bases": cutadapt["preqc_bases"],
        "postqc_reads": cutadapt["postqc_reads"],
        "postqc_bases": cutadapt["postqc_bases"]
    }

    # If a covstats file was specified
    if args.covstats:
        # Read the json file
        covstats = parse_json_file(args.covstats)

        # Determine the gender from the coverage data
        gender = determine_gender(covstats, args.female_threshold)
        d['gender'] = gender

        # Add the per chromosome coverage to the stats file
        chromosome_coverage = covstats['stats']['coverage']
        d['coverage'] = chromosome_coverage

    print(json.dumps(d))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample-name", required=True,
                        help="Sample name")
    parser.add_argument("--female-threshold", default=0.6, type=float,
                        help="Female threshold of X/all cov")
    parser.add_argument("--cutadapt", required=True,
                        help="Cutadapt summary output")
    parser.add_argument("--covstats", nargs="?",
                        help="Coverage statistics")

    args = parser.parse_args()
    main(args)
