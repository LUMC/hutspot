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
stats_to_tsv.py
~~~~~~~~~~~~~~~

:copyright: (c) 2017-2019 Sander Bollen
:copyright: (c) 2017-2019 Leiden University Medical Center
:license: AGPL-3.0
"""
import argparse
import json

from collections import OrderedDict
from pathlib import Path


def get_covstats(cov_d):
    s_d = cov_d['_all']
    return {
        'median_coverage': s_d['median'],
        'mean_coverage': s_d['mean'],
        'modal_coverage': s_d['mode'],
        'horizontal_coverage': s_d['horizontal'],
        'coverage_frac_min_100x': s_d['frac_min_100x'],
        'coverage_frac_min_10x': s_d['frac_min_10x'],
        'coverage_frac_min_20x': s_d['frac_min_20x'],
        'coverage_frac_min_30x': s_d['frac_min_30x'],
        'coverage_frac_min_40x': s_d['frac_min_40x'],
        'coverage_frac_min_50x': s_d['frac_min_50x'],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to input JSON", type=Path)

    args = parser.parse_args()

    with args.input.open() as handle:
        orig_dict = json.load(handle)

    sdicts = []

    for sample in orig_dict['sample_stats']:
        sname = sample['sample_name']
        sample_dict = OrderedDict()
        sample_dict.update({
            "sample_name": sname,
            "preqc_reads" : sample['preqc_reads'],
            "preqc_bases" : sample['preqc_bases'],
            "postqc_reads": sample['postqc_reads'],
            "postqc_bases": sample['postqc_bases'],
            "mapped_reads": int(sample['picard_AlignmentSummaryMetrics']['PF_HQ_ALIGNED_READS']),
            "mapped_bases": int(sample['picard_AlignmentSummaryMetrics']['PF_HQ_ALIGNED_BASES'])
        })
        if 'coverage' in sample:
            sample_dict.update(get_covstats(sample['coverage']))
        sdicts.append(sample_dict)

    lens = [len(list(x.keys())) for x in sdicts]
    longest_dict = sdicts[lens.index(max(lens))]

    header = "\t".join(longest_dict.keys())
    print(header)
    for sd in sdicts:
        vals = []
        for k in longest_dict.keys():
            if k in sd:
                vals.append(sd[k])
            else:
                vals.append("NA")
        print("\t".join(map(str, vals)))
