# Hutspot

This is a multi sample DNA variant calling pipeline based on Snakemake, bwa and
the GATK HaplotypeCaller.

## Features
* Any number of samples is supported
* Whole-genome calling, regardless of wet-lab library preparation.
* Follows modern best practices
    * Each sample is individually called as as a GVCF.
    * A VCF is then produced by genotyping the individual GVCFs separately
      for each sample.
* Data parallelization for calling and genotyping steps.
    * Using the `scatter_size` setting in the configuration file, the reference
      genome is split into chunks, and each chunk can be processed
      independenly. The default value of 1 billon will scatter the human
      reference genoom into 6 chunks.
* Reasonably fast.
    * 96 exomes in < 24 hours.
* No unnecessary jobs
* Calculate coverage metrics if a `bedfile` is specified.
* Fully containerized rules through singularity and biocontainers. Legacy
conda environments are no long available.

# Installation

To run this pipeline you will need the following at minimum:

* python 3.6
* snakemake 5.2.0 or newer

This repository contains a [conda](https://conda.io/docs/)
environment file that you can use to install all minimum dependencies in a
conda environment:

```bash
conda env create -f environment.yml
```

Alternatively, you can set up a python virtualenv and run

```bash
pip install -r requirements.txt
```

## Singularity

We highly recommend the user of the containerized rules through
[singularity](https://www.sylabs.io/singularity/).

This option does require you to install singularity on your system. As this
usually requires administrative privileges, singularity is not contained
within our provided conda environment file.

If you want to use singularity, make sure you install version 3 or higher.

### Debian
If you happen to use Debian buster, singularity 3.0.3 comes straight out
of the box with a simple:

```bash
sudo apt install singularity-container
```

### Docker

You can run singularity within a docker container. Please note that
the container **MUST** run in privileged mode for this to work.

We have provided our own container that includes singularity and snakemake
[here](https://hub.docker.com/r/lumc/singularity-snakemake).

### Manual install

If you don't use Debian buster and cannot run a privileged docker container,
you - unfortunately :-( - will have to install singularity manually.
Please see the installation instructions
[here](https://github.com/sylabs/singularity/blob/master/INSTALL.md) on how
to do that.


## Operating system

Hutspot was tested on Ubuntu 16.04 only.
It should reasonably work on most modern Linux distributions.

# Requirements

For every sample you wish to analyze, we require one or more paired end
readgroups in fastq format. They must be compressed with either `gzip` or
`bgzip`.

The configuration must be passed to the pipeline through a configuration file.
This is a json file listing the samples and their associated readgroups
as well as the other settings to be used.
An example config json can be found [here](config/example.json), and a
json schema describing the configuration file can be found [here](config/schema.json).
This json schema can also be used to validate your configuration file.

## Reference files

The following reference files **must** be provided in the configuration:

1. `reference`: A reference genome, in fasta format. Must be indexed with
   `samtools faidx`.
2. `dbsnp`: A dbSNP VCF file
3. `known_sites`: One ore more VCF files with known sites for base
    recalibration

The following reference files **may** be provided:

1. `targetsfile`: Bed file of the targets of the capture kit. Used to calculate coverage.
2. `baitsfile`: Bed file of the baits of the capture kit. Used to calculate picard HsMetric.
3. `refflat`: A refFlat file to calculate coverage over transcripts.
4. `scatter_size`: Size of the chunks to split the variant calling into.
5. `female_threshold`: Fraction of reads between X and the autosomes to call as
    female.


# How to run

After installing and activating the main conda environment, as described above,
the pipeline can be started with:

```bash
snakemake -s Snakefile \
--use-singularity \
--configfile tests/data/config/sample_config.json
```

This would start all jobs locally. Obviously this is not what one would
regularly do for a normal pipeline run. How to submit jobs on a cluster is
described later. Let's first move on to the necessary configuration values.

## Configuration values
The required and optional outputs are specified in the json schema located in
`config/schema.json`. Before running, the content of the `--configfile` is
validated against this schema.

The following configuration values are **required**:

| configuration | description |
| ------------- | ----------- |
| `reference` | Absolute path to fasta file |
| `samples` | One or more samples, with associated fastq files |
| `dbsnp` | Path to dbSNP VCF file|
| `known_sites` | Path to one or more VCF files with known sites. Can be the same as the `dbsnp` file|


The following configuration options are **optional**:

| configuration | description |
| ------------- | ----------- |
| `targetsfile` | Bed file of the targets of the capture kit. Used to calculate coverage |
| `baitsfile` | Bed file of the baits of the capture kit. Used to calculate picard HsMetrics |
| `female_threshold` | Float between 0 and 1 that signifies the threshold of the ratio between coverage on X/overall coverage that 'calls' a sample as female. Default = 0.6 |
| `scatter_size` | The size of chunks to divide the reference into for parallel execution. Default = 1000000000 |
| `coverage_threshold` | One or more threshold coverage values. For each value, a sample specific bed file will be created that contains the regions where the coverage is above the threshold |


## Cluster configuration

To run on a cluster, snakemake needs to be called with some extra arguments.
Additionally, it needs a cluster yaml file describing resources per job.

If you run on a cluster with drmaa support,an environment variable named
`DRMAA_LIBRARY_PATH` must be in the executing shell environment. This variable
points to the `.so` file of the DRMAA library.

An sge-cluster.yml is bundled with this pipeline in the cluster directory.
It is optimized for SGE clusters, where the default vmem limit is 4G.
If you run SLURM, or any other cluster system, you will have to write your own
cluster yaml file. Please see the [snakemake documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration)
for details on how to do so. Given the provided sge-cluster.yml, activating the
cluster mode can be done as follows:

```bash
snakemake -s Snakefile \
--cluster-config cluster/sge-cluster.yml
--drmaa ' -pe <PE_NAME> {cluster.threads} -q all.q -l h_vmem={cluster.vmem} -cwd -V -N hutspot' \
```

## Limitations
Sample names should be unique, and not overlap (such as `sample1` and
`sample10`). This is due to the way output files are parsed by multiQC,
when sample names overlap, the json output for picard DuplicationMetrics cannot
be parsed unambiguously.

## Binding additional directories under singularity

In singularity mode, snakemake binds the location of itself in the container.
The current working directory is also visible directly in the container.

In many cases, this is not enough, and will result in `FileNotFoundError`s.
E.g., suppose you run your pipeline in `/runs`, but your fastq files live in
`/fastq` and your reference genome lives in `/genomes`. We would have to bind
`/fastq` and `/genomes` in the container.

This can be accomplished with `--singularity-args`, which accepts a simple
string of arguments passed to singularity. E.g. in the above example,
we could do:

```bash
snakemake -S Snakefile \
--use-singularity  \
--singularity-args ' --bind /fastq:/fastq --bind /genomes:/genomes '
```

## Summing up

To sum up, a full pipeline run under a cluster would be called as:

```bash
snakemake -s Snakefile \
--use-singularity \
--singularity-args ' --bind /some_path:/some_path ' \
--cluster-config cluster/sge-cluster.yml \
--drmaa ' -pe <PE_NAME> {cluster.threads} -q all.q -l h_vmem={cluster.vmem} -cwd -V -N hutspot' \
--rerun-incomplete \
--jobs 200 \
-w 120 \
--max-jobs-per-second 30 \
--restart-times 2 \
--configfile config.json
```

# Graph

Below you can see the rule graph of the pipeline. The main variant calling flow
is highlighted in red. This only shows dependencies
between rules, and not between jobs. The actual job graph is considerably
more complex, as nearly all rules are duplicated by sample and some
(the scatter jobs) additionally by chunk.

As a rough estimate of the total number of jobs in pipeline you can use
the following formula:

<a href="https://www.codecogs.com/eqnedit.php?latex=n_{jobs}&space;=&space;4&plus;(21*n_{samples})&plus;(n_{samples}*n_{beds})&plus;(n_{samples}*n_{chunks})&plus;n_{chunks}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n_{jobs}&space;=&space;4&plus;(21*n_{samples})&plus;(n_{samples}*n_{beds})&plus;(n_{samples}*n_{chunks})&plus;n_{chunks}" title="n_{jobs} = 4+(21*n_{samples})+(n_{samples}*n_{beds})+(n_{samples}*n_{chunks})+n_{chunks}" /></a>

<!---
Note: math doesn't work on github. The following _does_ work in gitlab
```math
jobs = 4+(21*n_{samples})+(1*n_{samples}*n_{beds})+(1*n_{samples}*n_{chunks})+(1*n_{chunks})
```
--->

This gives about 12,000 jobs for a 96-sample run with 2 bed files and 100 chunks.

NOTE: the graph will only render if your markdown viewer supports `plantuml`.
Having trouble viewing the graph? See [this](img/rulegraph.svg) static SVG instead.

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.30 0.6 0.85", style="rounded"];
	1[label = "multiqc", color = "0.60 0.6 0.85", style="rounded"];
	2[label = "merge_stats", color = "0.17 0.6 0.85", style="rounded"];
	3[label = "bai", color = "0.09 0.6 0.85", style="rounded"];
	4[label = "genotype_gather\nsample: micro", color = "0.06 0.6 0.85", style="rounded"];
	5[label = "gvcf_gather\nsample: micro", color = "0.32 0.6 0.85", style="rounded"];
	6[label = "fastqc_raw\nsample: micro", color = "0.00 0.6 0.85", style="rounded"];
	7[label = "fastqc_merged", color = "0.11 0.6 0.85", style="rounded"];
	8[label = "fastqc_postqc", color = "0.02 0.6 0.85", style="rounded"];
	9[label = "stats_tsv", color = "0.45 0.6 0.85", style="rounded"];
	10[label = "collectstats", color = "0.24 0.6 0.85", style="rounded"];
	11[label = "vcfstats\nsampel: micro", color = "0.52 0.6 0.85", style="rounded"];
	12[label = "markdup", color = "0.47 0.6 0.85", style="rounded"];
	13[label = "scatterregions", color = "0.56 0.6 0.85", style="rounded"];
	14[label = "merge_r1\nsample: micro", color = "0.65 0.6 0.85", style="rounded"];
	15[label = "merge_r2\nsample: micro", color = "0.26 0.6 0.85", style="rounded"];
	16[label = "cutadapt", color = "0.22 0.6 0.85", style="rounded"];
	17[label = "fqcount_preqc", color = "0.37 0.6 0.85", style="rounded"];
	18[label = "fqcount_postqc", color = "0.58 0.6 0.85", style="rounded"];
	19[label = "mapped_reads_bases", color = "0.43 0.6 0.85", style="rounded"];
	20[label = "unique_reads_bases", color = "0.34 0.6 0.85", style="rounded"];
	21[label = "fastqc_stats", color = "0.13 0.6 0.85", style="rounded"];
	22[label = "covstats", color = "0.39 0.6 0.85", style="rounded"];
	23[label = "align", color = "0.49 0.6 0.85", style="rounded"];
	24[label = "create_markdup_tmp", color = "0.41 0.6 0.85", style="rounded,dashed"];
	25[label = "sickle", color = "0.19 0.6 0.85", style="rounded"];
	26[label = "genome", color = "0.62 0.6 0.85", style="rounded"];
	1 -> 0
	2 -> 0
	3 -> 0
	4 -> 0
	5 -> 0
	6 -> 0
	7 -> 0
	8 -> 0
	9 -> 1
	10 -> 2
	11 -> 2
	12 -> 3
	13 -> 4
	13 -> 5
	14 -> 7
	15 -> 7
	16 -> 8
	2 -> 9
	17 -> 10
	18 -> 10
	19 -> 10
	20 -> 10
	21 -> 10
	22 -> 10
	4 -> 11
	23 -> 12
	24 -> 12
	25 -> 16
	14 -> 17
	15 -> 17
	16 -> 18
	23 -> 19
	12 -> 20
	7 -> 21
	8 -> 21
	12 -> 22
	26 -> 22
	16 -> 23
	24 -> 23
	14 -> 25
	15 -> 25
}
```

LICENSE
=======

AGPL-3.0
