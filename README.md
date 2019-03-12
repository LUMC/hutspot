# Hutspot

This is a multisample DNA variant calling pipeline based on Snakemake, bwa and the
GATK HaplotypeCaller.  

## Features 
* Any number of samples is supported
* Whole-genome calling, regardless of wet-lab library preparation. 
* Follows modern best practices
    * Each sample is individually called as as a GVCF. 
    * A multisample VCF is then produced by genotyping the collection of GVCFs.
* Data parallelization for calling and genotyping steps.
    * Using ~100 chunks, we call an entire exome in ~15 minutes!
* Reasonably fast.
    * 96 exomes in < 24 hours.
* No unnecessary jobs
* Coverage metrics for any number of bed files.
* Separate conda environments for **every** step. No more dependency hell!
Every job can potentially use different versions of the same package.
* Optionally sub-sample inputs when number of bases exceeds a user-defined
threshold.

# Installation

We recommend the use of [conda](https://conda.io/docs/) for installing all
dependencies. All rules have a separate conda environment, which guarantees
every tool can use its own dependencies.

To install the base environment containing snakemake itself, activate conda
and run the following in your terminal:

`conda env create -f environment.yml`

Subsequently running the pipeline with `--use-conda` will make sure 
the correct conda environments get created. This requires a working
internet connection. If you do not want conda environment to be created for
each pipeline run, use the `--conda-prefix` argument. See the
[snakemake documentation](http://snakemake.readthedocs.io/en/stable/executable.html)
for more information. 

## GATK

For license reasons, conda cannot fully install the GATK. The JAR 
must be registered by running `gatk-register` after the environment is
created, which conflicts with the automated environment creation.
 
For this reason, hutspot **requires** you to manually specify the path to
the GATK executable JAR via `--config GATK=/path/to/gatk.jar`.

## Fastq-count

Several steps in the pipeline collect fastq metrics via [fastq-count](https://github.com/sndrtj/fastq-count).
This is a small tool implemented in Rust for speed reasons. As this tool
is not yet in conda, it must be compiled on the user's system before 
running the pipeline. When compiled, the path to the executable can be
supplied via `--config FASTQ_COUNT=/path/to/fastq-count`.

A drop-in replacement implemented in python exists in this repository.
Not specifing the `FASTQ_COUNT` config value will use the python replacement.
Do note that this python replacement is an order of magnitude slower.

## Operating system

Hutspot was tested on Ubuntu 16.04 only.
It should reasonably work on most modern Linux distributions. 
   

# Requirements

For every sample you wish to analyze, we require one or more paired end
readgroups in fastq format. They must be compressed with either `gzip` or
`bgzip`.

Samples must be passed to the pipeline through a config file. This is a
simple json file listing the samples and their associated readgroups/libraries.
An example config json can be found [here](config/example.json), and a
json schema describing the configuration file can be found [here](config/schema.json). 
This json schema can also be used to validate your configuration file.

## Reference files

The following reference files **must** be provided:

1. A reference genome, in fasta format. Must be indexed with `samtools faidx`.
2. A dbSNP VCF file
3. A VCF file from 1000Genomes
4. A VCF file from the HapMap project.

The following reference files **may** be provided:

1. Any number of BED files to calculate coverage on.


# How to run

After installing and activating the main conda environment, as described above,
the pipeline can be started with:

```bash
snakemake -s Snakefile \
--use-conda \
-T \
--config <CONFIGURATION VALUES>
```

This would start all jobs locally. Obviously this is not what one would
regularly do for a normal pipeline run. How to submit jobs on a cluster is
described later. Let's first move on to the necessary configuration values.

## Configuration values

The following configuration values are **required**:

| configuration | description |
| ------------- | ----------- |
| `OUT_DIR` | Absolute path to output directory |
| `REFERENCE` | Absolute path to fasta file |
| `SAMPLE_CONFIG` | Path to config file as described above |
| `GATK` | Path to GATK jar. **Must** be version 3.7  |
| `DBSNP` | Path to dbSNP VCF |
| `ONETHOUSAND` | Path to 1000Genomes VCF |
| `HAPMAP` | Path to HapMap VCF |

The following configuration options are **optional**:

| configuration | description |
| ------------- | ----------- |
| `BED` | Comma-separate list of paths to BED files of interest |
| `FEMALE_THRESHOLD` | Float between 0 and 1 that signifies the threshold of the ratio between coverage on X/overall coverage that 'calls' a sample as female. Default = 0.6 |
| `FASTQ_COUNT` | Path to `fastq-count` executable |
| `MAX_BASES` | Maximum allowed number of bases per sample before subsampling. Default = None (no subsampling) |


## Cluster configuration

To run on a cluster, snakemake needs to be called with some extra arguments.
Additionally, it needs a cluster yaml file describing resources per job.

In all cases, an environment variable named `DRMAA_LIBRARY_PATH` must be
in the executing shell environment. This variable points to the `.so` file
of the DRMAA library.

A cluster.yml is bundled with this pipeline. It is optimized for SGE clusters,
where the default vmem limit is 4G. If you run SLURM, or any other cluster
system, you will have to write your own cluster yaml file. Please see the
[snakemake documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration)
for details on how to do so. Given the provided cluster.yml, activating the
cluster mode can be done as follows:

```bash
snakemake -s Snakefile \
--cluster-config cluster.yml
--drmaa ' -pe <PE_NAME> {cluster.threads} -q all.q -l h_vmem={cluster.vmem} -cwd -V -N hutspot' \
```

## Summing up

To sum up, a full pipeline run under a cluster would be called as:

```bash
snakemake -s Snakefile \
--use-conda \
--cluster-config cluster.yml \
--drmaa ' -pe <PE_NAME> {cluster.threads} -q all.q -l h_vmem={cluster.vmem} -cwd -V -N hutspot' \
--rerun-incomplete \
--jobs 200 \
-w 120 \
--max-jobs-per-second 30 \
--restart-times 2 \
-T \
--config SAMPLE_CONFIG=samples.json \
OUTPUT_DIR=/path/to/odir \
REFERENCE=/path/to/genome.fasta \
GATK=/path/to/GenomeAnalysisTK.jar \
DBSNP=/path/to/dbsnp.vcf.gz \
ONETHOUSAND=/path/to/onekg.vcf \
HAPMAP=/path/to/hapmap.vcf \
FASTQ_COUNT=/path/to/fastq-count \
BED=/path/to/interesting_region.bed
```

# Graph

Below you can see the rulegraph of the pipeline. The main variant calling flow
is highlighted in red. This only shows dependencies
between rules, and not between jobs. The actual job graph is considerably
more complex, as nearly all rules are duplicated by sample and some
(the scatter jobs) additionally by chunk. 

As a rough estimate of the total number of jobs in pipeline you can use
the following formula:

```math
jobs = 4+(21*n_{samples})+(1*n_{samples}*n_{beds})+(1*n_{samples}*n_{chunks})+(1*n_{chunks})
``` 

This gives about 12,000 jobs for a 96-sample run with 2 bed files and 100 chunks.

NOTE: the graph will only render if your markdown viewer supports `plantuml`.
Having trouble viewing the graph? See [this](img/rulegraph.svg) static SVG in stead.

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    rankdir=LR;
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.32 0.6 0.85", style="rounded"];
	1[label = "fastqc_postqc", color = "0.48 0.6 0.85", style="rounded"];
	2[label = "fastqc_raw", color = "0.05 0.6 0.85", style="rounded"];
	3[label = "fastqc_merged", color = "0.57 0.6 0.85", style="rounded"];
	4[label = "merge_stats", color = "0.21 0.6 0.85", style="rounded"];
	5[label = "genotype_gather", color = "0.37 0.6 0.85", style="rounded"];
	6[label = "cutadapt", color = "0.30 0.6 0.85", style="rounded"];
	7[label = "merge_r2", color = "0.60 0.6 0.85", style="rounded"];
	8[label = "merge_r1", color = "0.23 0.6 0.85", style="rounded"];
	9[label = "collectstats", color = "0.02 0.6 0.85", style="rounded"];
	10[label = "vcfstats", color = "0.07 0.6 0.85", style="rounded"];
	11[label = "genotype_scatter", color = "0.16 0.6 0.85", style="rounded"];
	12[label = "sickle", color = "0.46 0.6 0.85", style="rounded"];
	13[label = "mapped_num", color = "0.64 0.6 0.85", style="rounded"];
	14[label = "fqcount_postqc", color = "0.00 0.6 0.85", style="rounded"];
	15[label = "unique_num", color = "0.11 0.6 0.85", style="rounded"];
	16[label = "covstats", color = "0.41 0.6 0.85", style="rounded"];
	17[label = "fqcount_preqc", color = "0.62 0.6 0.85", style="rounded"];
	18[label = "usable_basenum", color = "0.14 0.6 0.85", style="rounded"];
	19[label = "mapped_basenum", color = "0.09 0.6 0.85", style="rounded"];
	20[label = "gvcf_gather", color = "0.55 0.6 0.85", style="rounded"];
	21[label = "seqtk_r2", color = "0.51 0.6 0.85", style="rounded"];
	22[label = "seqtk_r1", color = "0.28 0.6 0.85", style="rounded"];
	23[label = "align", color = "0.25 0.6 0.85", style="rounded"];
	24[label = "markdup", color = "0.18 0.6 0.85", style="rounded"];
	25[label = "genome", color = "0.53 0.6 0.85", style="rounded"];
	26[label = "gvcf_scatter", color = "0.34 0.6 0.85", style="rounded"];
	27[label = "printreads", color = "0.39 0.6 0.85", style="rounded"];
	28[label = "baserecal", color = "0.44 0.6 0.85", style="rounded"];
	4 -> 0
	2 -> 0
	1 -> 0
	5 -> 0 [color="red"]
	3 -> 0
	6 -> 1
	8 -> 3
	7 -> 3
	10 -> 4
	9 -> 4
	11 -> 5 [color="red"]
	12 -> 6 [color="red"]
	14 -> 9
	19 -> 9
	18 -> 9
	16 -> 9
	15 -> 9
	17 -> 9
	13 -> 9
	5 -> 10
	20 -> 11 [color="red"]
	22 -> 12 [color="red"]
	21 -> 12 [color="red"]
	23 -> 13
	6 -> 14
	24 -> 15
	24 -> 16
	25 -> 16
	8 -> 17 [color="red"]
	7 -> 17 [color="red"]
	24 -> 18
	23 -> 19
	26 -> 20 [color="red"]
	17 -> 21 [color="red"]
	7 -> 21 [color="red"]
	17 -> 22 [color="red"]
	8 -> 22 [color="red"]
	6 -> 23 [color="red"]
	23 -> 24 [color="red"]
	27 -> 26 [color="red"]
	24 -> 27 [color="red"]
	28 -> 27 [color="red"]
	24 -> 28 [color="red"]
}  
```

LICENSE
=======

AGPL-3.0