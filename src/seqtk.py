"""
Little script from running seqtk with conda

Conda directives can't be used with a run directive,
so must be combined with script directive in stead.

This script assumes the following:
 - a `snakemake` object exists,
 - this object has the following attributes:
    - input: a list of two items:
        1. output of fastq-count as path to json file
        2. a fastq file to be sub-sampled
    - output: a list of one item containing path to output file
    - params: a list of one item containing the max number of bases
 - a `shell` function exists

This will _not_ work outside of a snakemake context.
"""
import json


def subsample(json_path, fastq_path, opath, max_bases):
    with open(json_path) as handle:
        bases = json.load(handle)['bases']
    if max_bases == "":
        frac = 100
    else:
        frac = max_bases / float(bases)

    if frac > 1:
        snakemake.shell("ln -s {0} {1}".format(fastq_path, opath))
    else:
        snakemake.shell("seqtk sample -s100 {0} {1} | gzip -c > {2}".format(fastq_path,
                                                                  frac,
                                                                  opath))


subsample(snakemake.input[0], snakemake.input[1],
          snakemake.output[0], snakemake.params[0])



