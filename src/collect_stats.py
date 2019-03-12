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
import click
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


@click.command()
@click.option("--sample-name",
              type=click.STRING,
              required=True,
              help="Sample name")
@click.option("--pre-qc-fastq",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="pre-qc json from fastq-count")
@click.option("--post-qc-fastq",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Post-qc json from fastq-count")
@click.option("--mapped-num",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Mapped num file")
@click.option("--mapped-basenum",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Mapped basenum file")
@click.option("--unique-num",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Unique num file")
@click.option("--usable-basenum",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Usable basenum")
@click.option("--female-threshold",
              type=click.FLOAT,
              default=0.6,
              help="Female threshold of X/all cov")
@click.option("--fastqc-stats",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              help="Path to fastqc stats json")
@click.argument("covstats",
                type=click.Path(dir_okay=False, exists=True, readable=True),
                nargs=-1)
def main(sample_name, pre_qc_fastq, post_qc_fastq, mapped_num, mapped_basenum,
         unique_num, usable_basenum, female_threshold, fastqc_stats, covstats):

    preqcd = parse_json_file(pre_qc_fastq)
    posqcd = parse_json_file(post_qc_fastq)
    fastqc = parse_json_file(fastqc_stats)

    mpnum = parse_num_file(mapped_num)
    mpbnum = parse_num_file(mapped_basenum)
    unum = parse_num_file(unique_num)
    ubnum = parse_num_file(usable_basenum)

    covl = []
    for c in covstats:
        cd = parse_json_file(c)
        cdd = {
            "name": basename(c),
            "gender": determine_gender(cd, female_threshold),
            "covstats": cd
        }
        covl.append(cdd)

    d = {
        "sample_name": sample_name,
        "pre_qc_fastq_count": preqcd,
        "post_qc_fastq_count": posqcd,
        "n_mapped_reads": mpnum,
        "n_mapped_bases": mpbnum,
        "n_usable_reads": unum,
        "n_usable_bases": ubnum,
        "fastqc": fastqc,
        "covstats": covl
    }

    print(json.dumps(d))


if __name__ == "__main__":
    main()
