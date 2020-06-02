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


def add_picard_insertSize(data, filename):
    """ Add the picard insertSize for each sample to data """
    insert = parse_json(filename)

    for sample in insert.values():
        name = sample.pop('SAMPLE_NAME')
        for d in data['sample_stats']:
            if d['sample_name'] == name:
                d['picard_insertSize'] = sample
                break
        else:
            raise RuntimeError(f"Unknown sample {name}")


def add_picard_HsMetrics(data, filename):
    """ Add the picard HsMetrics for each sample to data """
    HsMetrics = parse_json(filename)

    for sample in HsMetrics:
        for d in data['sample_stats']:
            if d['sample_name'] == sample:
                d['picard_HsMetrics'] = HsMetrics[sample]
                break
        else:
            raise RuntimeError(f"Unknown sample {sample}")


def main(args):
    data = dict()
    data["sample_stats"] = list()

    for stats in args.collectstats:
        cs = parse_json(stats)
        data["sample_stats"].append(cs)

    if args.picard_insertSize:
        add_picard_insertSize(data, args.picard_insertSize)

    if args.picard_HsMetrics:
        add_picard_HsMetrics(data, args.picard_HsMetrics)
    print(json.dumps(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--collectstats',
                        nargs='+',
                        required=True,
                        help='Path to the collected stats for each sample')
    parser.add_argument('--picard-insertSize',
                        required=False,
                        help=('Path to multiQC json summary for picard '
                        'insertSize'))
    parser.add_argument('--picard-HsMetrics',
                        required=False,
                        nargs='?',
                        help=('Path to multiQC json summary for picard '
                        'HsMetrics'))
    args = parser.parse_args()
    main(args)
