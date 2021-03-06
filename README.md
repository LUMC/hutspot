[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3251553.svg)](https://doi.org/10.5281/zenodo.3251553)


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
* Fully containerized rules through singularity and biocontainers. Legacy 
conda environments are available as well.  
* Optionally sub-sample inputs when number of bases exceeds a user-defined
threshold.

# Installation

To run this pipeline you will need the following at minimum:

* python 3.6
* snakemake 5.2.0 or newer
* pyfaidx 

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

This option does, however,
require you to install singularity on your system. As this usually requires 
administrative privileges, singularity is not contained within our provided
conda environment file.

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


## GATK

For license reasons, conda and singularity cannot fully install the GATK. The JAR 
must be registered by running `gatk-register` after the environment is
created, which conflicts with the automated environment/container creation.
 
For this reason, hutspot **requires** you to manually specify the path to
the GATK executable JAR via `--config GATK=/path/to/gatk.jar`.

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
--use-singularity \
--config <CONFIGURATION VALUES>
```

This would start all jobs locally. Obviously this is not what one would
regularly do for a normal pipeline run. How to submit jobs on a cluster is
described later. Let's first move on to the necessary configuration values.

## Configuration values

The following configuration values are **required**:

| configuration | description |
| ------------- | ----------- |
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
--config SAMPLE_CONFIG=samples.json \
REFERENCE=/path/to/genome.fasta \
GATK=/path/to/GenomeAnalysisTK.jar \
DBSNP=/path/to/dbsnp.vcf.gz \
ONETHOUSAND=/path/to/onekg.vcf \
HAPMAP=/path/to/hapmap.vcf \
FASTQ_COUNT=/path/to/fastq-count \
BED=/path/to/interesting_region.bed
```

## Using conda instead of singularity

Legacy conda environments are also available for each and every rule. 
Simply use `--use-conda` instead of `--use-singularity` to enable conda
environments.

As dependency conflicts can and do arise with conda, it is recommended to 
combine this flag with `--conda-prefix`, such that you only have to 
build the environments once.

The conda environments use the same versions of tools as the singularity
containers, bar one:

* `fastqc` uses version 0.11.5 on conda, but 0.11.7 on singularity.    

# Graph

Below you can see the rulegraph of the pipeline. The main variant calling flow
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
    rankdir=LR;
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.62 0.6 0.85", style="rounded"];
	1[label = "genotype_gather", color = "0.31 0.6 0.85", style="rounded"];
	2[label = "multiqc", color = "0.14 0.6 0.85", style="rounded"];
	3[label = "bai", color = "0.41 0.6 0.85", style="rounded"];
	4[label = "split_vcf", color = "0.53 0.6 0.85", style="rounded"];
	5[label = "fastqc_raw", color = "0.63 0.6 0.85", style="rounded"];
	6[label = "fastqc_merged", color = "0.24 0.6 0.85", style="rounded"];
	7[label = "fastqc_postqc", color = "0.26 0.6 0.85", style="rounded"];
	8[label = "vtools_coverage", color = "0.58 0.6 0.85", style="rounded"];
	9[label = "merge_stats", color = "0.36 0.6 0.85", style="rounded"];
	10[label = "genotype_scatter", color = "0.09 0.6 0.85", style="rounded"];
	11[label = "genotype_chunkfile", color = "0.29 0.6 0.85", style="rounded"];
	12[label = "stats_tsv", color = "0.51 0.6 0.85", style="rounded"];
	13[label = "markdup", color = "0.55 0.6 0.85", style="rounded"];
	14[label = "genotype_gather_tbi", color = "0.19 0.6 0.85", style="rounded"];
	15[label = "merge_r1", color = "0.60 0.6 0.85", style="rounded"];
	16[label = "merge_r2", color = "0.10 0.6 0.85", style="rounded"];
	17[label = "cutadapt", color = "0.17 0.6 0.85", style="rounded"];
	18[label = "gvcf_gather", color = "0.32 0.6 0.85", style="rounded"];
	19[label = "gvcf_gather_tbi", color = "0.27 0.6 0.85", style="rounded"];
	20[label = "collectstats", color = "0.03 0.6 0.85", style="rounded"];
	21[label = "vcfstats", color = "0.00 0.6 0.85", style="rounded"];
	22[label = "align", color = "0.05 0.6 0.85", style="rounded"];
	23[label = "create_markdup_tmp", color = "0.44 0.6 0.85", style="rounded"];
	24[label = "sickle", color = "0.39 0.6 0.85", style="rounded"];
	25[label = "gvcf_scatter", color = "0.02 0.6 0.85", style="rounded"];
	26[label = "gvcf_chunkfile", color = "0.56 0.6 0.85", style="rounded"];
	27[label = "fqcount_preqc", color = "0.38 0.6 0.85", style="rounded"];
	28[label = "fqcount_postqc", color = "0.12 0.6 0.85", style="rounded"];
	29[label = "mapped_num", color = "0.50 0.6 0.85", style="rounded"];
	30[label = "mapped_basenum", color = "0.43 0.6 0.85", style="rounded"];
	31[label = "unique_num", color = "0.65 0.6 0.85", style="rounded"];
	32[label = "usable_basenum", color = "0.22 0.6 0.85", style="rounded"];
	33[label = "fastqc_stats", color = "0.46 0.6 0.85", style="rounded"];
	34[label = "covstats", color = "0.07 0.6 0.85", style="rounded"];
	35[label = "seqtk_r1", color = "0.34 0.6 0.85", style="rounded"];
	36[label = "seqtk_r2", color = "0.21 0.6 0.85", style="rounded"];
	37[label = "baserecal", color = "0.48 0.6 0.85", style="rounded"];
	38[label = "genome", color = "0.15 0.6 0.85", style="rounded"];
	9 -> 0
	4 -> 0 [color = "red"]
	3 -> 0
	6 -> 0
	7 -> 0
	1 -> 0
	8 -> 0
	2 -> 0
	5 -> 0
	11 -> 1 [color = "red"]
	10 -> 1 [color = "red"]
	12 -> 2
	13 -> 3
	1 -> 4 [color = "red"]
	14 -> 4 [color = "red"]
	16 -> 6
	15 -> 6
	17 -> 7
	19 -> 8
	18 -> 8
	20 -> 9
	21 -> 9
	19 -> 10 [color = "red"]
	18 -> 10 [color = "red"]
	9 -> 12
	23 -> 13 [color = "red"]
	22 -> 13 [color = "red"]
	1 -> 14  [color = "red"]
	24 -> 17 [color = "red"]
	25 -> 18 [color = "red"]
	26 -> 18 [color = "red"]
	18 -> 19 [color = "red"]
	28 -> 20
	27 -> 20
	32 -> 20
	30 -> 20
	33 -> 20
	34 -> 20
	29 -> 20
	31 -> 20
	1 -> 21
	14 -> 21
	17 -> 22 [color = "red"]
	36 -> 24 [color = "red"]
	35 -> 24 [color = "red"]
	37 -> 25 [color = "red"]
	13 -> 25 [color = "red"]
	16 -> 27 [color = "red"]
	15 -> 27 [color = "red"]
	17 -> 28
	22 -> 29
	22 -> 30
	13 -> 31
	13 -> 32
	7 -> 33
	6 -> 33
	38 -> 34
	13 -> 34
	27 -> 35 [color = "red"]
	15 -> 35 [color = "red"]
	27 -> 36 [color = "red"]
	16 -> 36 [color = "red"]
	13 -> 37 [color = "red"]
}
```

LICENSE
=======

AGPL-3.0
