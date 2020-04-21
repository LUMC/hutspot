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


def parse_json(path):
    with open(path) as handle:
        return json.load(handle)


def main(vcfstats, collectstats):
    data = dict()
    data["sample_stats"] = list()

    for vcf, stats in zip(vcfstats, collectstats):
        v = parse_json(vcf)
        cs = parse_json(stats)
        cs['vcfstats'] = v
        data["sample_stats"].append(cs)
    print(json.dumps(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcfstats',
                        nargs='+',
                        required=True,
                        help='Path to the vcfstats json for each sample')
    parser.add_argument('--collectstats',
                        nargs='+',
                        required=True,
                        help='Path to the collected stats for each sample')
    args = parser.parse_args()
    assert len(args.vcfstats) == len(args.collectstats)
    main(args.vcfstats, args.collectstats)
