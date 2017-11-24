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
    all = cv['_all']

    if 'chrX' in cv:
        x = cv['chrX']
    elif 'X' in cv:
        x = cv['X']
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
              type=click.INT,
              default=0.6,
              help="Female threshold of X/all cov")
@click.option("--covstats",
              type=click.Path(dir_okay=False, exists=True, readable=True),
              required=True,
              nargs=-1,
              help="Covstat json files")
def main(pre_qc_fastq, post_qc_fastq, mapped_num, mapped_basenum, unique_num,
         usable_basenum, female_threshold, covstats):

    preqcd = parse_json_file(pre_qc_fastq)
    posqcd = parse_json_file(post_qc_fastq)

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
        "pre_qc_fastq_count": preqcd,
        "post_qc_fastq_count": posqcd,
        "n_mapped_reads": mpnum,
        "n_mapped_bases": mpbnum,
        "n_usable_reads": unum,
        "n_usable_bases": ubnum,
        "covstats": covl
    }

    json.dumps(covl)