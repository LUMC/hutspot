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
Main Snakefile for the pipeline.

:copyright: (c) 2017-2019 Sander Bollen
:copyright: (c) 2017-2019 Leiden University Medical Center
:license: AGPL-3.0
"""

import json
from functools import partial
from os.path import join, basename
from pathlib import Path

from pyfaidx import Fasta

REFERENCE = config.get("REFERENCE")
if REFERENCE is None:
    raise ValueError("You must set --config REFERENCE=<path>")
if not Path(REFERENCE).exists():
    raise FileNotFoundError("Reference file {0} "
                            "does not exist.".format(REFERENCE))

GATK = config.get("GATK")
if GATK is None:
    raise ValueError("You must set --config GATK=<path>")
if not Path(GATK).exists():
    raise FileNotFoundError("{0} does not exist.".format(GATK))

# these are all optional
BED = config.get("BED", "")  # comma-separated list of BED files
REFFLAT = config.get("REFFLAT", "")  # comma-separated list of refFlat files
FEMALE_THRESHOLD = config.get("FEMALE_THRESHOLD", 0.6)
FASTQ_COUNT = config.get("FASTQ_COUNT")
MAX_BASES = config.get("MAX_BASES", "")

# Make sure all files with known sites exist
known_sites = config.get("KNOWN_SITES").split(',')
for filename in known_sites:
    if not Path(filename).exists():
        raise FileNotFoundError("{0} does not exist".format(DBSNP))

# Generate the input string for basrecalibration
known_sites_argument = ''
for argument, filename in zip(repeat('-knownSites'), known_sites):
    known_sites_argument +=' {argument} {filename}'.format(argument, filename)

def fsrc_dir(*args):
    """Wrapper around snakemake's srcdir to work like os.path.join"""
    if len(args) == 1:
        return srcdir(args[0])
    return srcdir(join(*args))

covpy = fsrc_dir("src", "covstats.py")
colpy = fsrc_dir("src", "collect_stats.py")
mpy = fsrc_dir("src", "merge_stats.py")
seq = fsrc_dir("src", "seqtk.sh")
fqpy = fsrc_dir("src", "fastqc_stats.py")
tsvpy = fsrc_dir("src", "stats_to_tsv.py")
fqsc = fsrc_dir("src", "safe_fastqc.sh")

# sample config parsing
SCONFIG = config.get("SAMPLE_CONFIG")
if SCONFIG is None:
    raise ValueError("You must set --config SAMPLE_CONFIG=<path>")
if not Path(SCONFIG).exists():
    raise FileNotFoundError("{0} does not exist".format(SCONFIG))

with open(config.get("SAMPLE_CONFIG")) as handle:
    SAMPLE_CONFIG = json.load(handle)
SAMPLES = SAMPLE_CONFIG['samples'].keys()

if BED != "":
    BEDS = BED.split(",")
else:
    BEDS = []

if REFFLAT != "":
    REFFLATS = REFFLAT.split(",")
else:
    REFFLATS = []

BASE_BEDS = [basename(x) for x in BEDS]
BASE_REFFLATS = [basename(x) for x in REFFLATS]


def split_genome(ref, approx_n_chunks=100):
    """
    Split genome in chunks.

    Chunks are strings in the format: `<ctg>:<start>-<end>`
    These follow the region string format as used by htslib,
    which uses _1_-based indexing.
    See: http://www.htslib.org/doc/tabix.html
    """
    fa = Fasta(ref)
    tot_size = sum([len(x) for x in fa.records.values()])
    chunk_size = tot_size//approx_n_chunks
    chunks = []
    for chrom_name, chrom_value in fa.records.items():
        pos = 1
        while pos <= len(chrom_value):
            end = pos+chunk_size
            if end <= len(chrom_value):
                chunk = "{0}:{1}-{2}".format(chrom_name, pos, end)
            else:
                chunk = "{0}:{1}-{2}".format(chrom_name, pos, len(chrom_value))
            chunks.append(chunk)
            pos = end
    return chunks

CHUNKS = split_genome(REFERENCE)



def get_r(strand, wildcards):
    """Get fastq files on a single strand for a sample"""
    s = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    rs = []
    for l in sorted(s['libraries'].keys()):
        rs.append(s['libraries'][l][strand])
    return rs

get_r1 = partial(get_r, "R1")
get_r2 = partial(get_r, "R2")


def get_bedpath(wildcards):
    """Get absolute path of a bed file"""
    return next(x for x in BEDS if basename(x) == wildcards.bed)


def get_refflatpath(wildcards):
    """Get absolute path of a refflat file"""
    return next(x for x in REFFLATS if basename(x) == wildcards.ref)


def sample_gender(wildcards):
    """Get sample gender"""
    sam = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    return sam.get("gender", "null")


def metrics(do_metrics=True):
    if not do_metrics:
        return ""

    fqcr = expand("{sample}/pre_process/raw_fastqc/.done.txt",
                  sample=SAMPLES)
    fqcm = expand("{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip",
                  sample=SAMPLES)
    fqcp = expand("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip",
                  sample=SAMPLES)
    if len(REFFLATS) >= 1:
        coverage_stats = expand("{sample}/coverage/{ref}.coverages.tsv",
                                sample=SAMPLES, ref=BASE_REFFLATS)
    else:
        coverage_stats = []
    stats = "stats.json"
    return  fqcr + fqcm + fqcp + coverage_stats + [stats]


localrules: gvcf_chunkfile, genotype_chunkfile


rule all:
    input:
        combined="multisample/genotyped.vcf.gz",
        report="multiqc_report/multiqc_report.html",
        bais=expand("{sample}/bams/{sample}.markdup.bam.bai",sample=SAMPLES),
        vcfs=expand("{sample}/vcf/{sample}_single.vcf.gz",sample=SAMPLES),
        stats=metrics()


rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output: directory("tmp")
    singularity: "docker://debian:buster-slim"
    shell: "mkdir -p {output}"

rule genome:
    """Create genome file as used by bedtools"""
    input: REFERENCE
    output: "current.genome"
    singularity: "docker://debian:buster-slim"
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule merge_r1:
    """Merge all forward fastq files into one"""
    input: get_r1
    output: temp("{sample}/pre_process/{sample}.merged_R1.fastq.gz")
    singularity: "docker://debian:buster-slim"
    shell: "cat {input} > {output}"

rule merge_r2:
    """Merge all reverse fastq files into one"""
    input: get_r2
    output: temp("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    singularity: "docker://debian:buster-slim"
    shell: "cat {input} > {output}"


rule seqtk_r1:
    """Either subsample or link forward fastq file"""
    input:
        stats="{sample}/pre_process/{sample}.preqc_count.json",
        fastq="{sample}/pre_process/{sample}.merged_R1.fastq.gz",
        seqtk=seq
    params:
        max_bases=str(MAX_BASES)
    output:
        fastq=temp("{sample}/pre_process/{sample}.sampled_R1.fastq.gz")
    conda: "envs/seqtk.yml"
    singularity: "docker://quay.io/biocontainers/mulled-v2-13686261ac0aa5682c680670ff8cda7b09637943:d143450dec169186731bb4df6f045a3c9ee08eb6-0"
    shell: "bash {input.seqtk} {input.stats} {input.fastq} {output.fastq} "
           "{params.max_bases}"


rule seqtk_r2:
    """Either subsample or link reverse fastq file"""
    input:
        stats = "{sample}/pre_process/{sample}.preqc_count.json",
        fastq = "{sample}/pre_process/{sample}.merged_R2.fastq.gz",
        seqtk=seq
    params:
        max_bases =str(MAX_BASES)
    output:
        fastq = temp("{sample}/pre_process/{sample}.sampled_R2.fastq.gz")
    conda: "envs/seqtk.yml"
    singularity: "docker://quay.io/biocontainers/mulled-v2-13686261ac0aa5682c680670ff8cda7b09637943:d143450dec169186731bb4df6f045a3c9ee08eb6-0"
    shell: "bash {input.seqtk} {input.stats} {input.fastq} {output.fastq} "
           "{params.max_bases}"


# contains original merged fastq files as input to prevent them from being prematurely deleted
rule sickle:
    """Trim fastq files"""
    input:
        r1 = "{sample}/pre_process/{sample}.sampled_R1.fastq.gz",
        r2 = "{sample}/pre_process/{sample}.sampled_R2.fastq.gz",
        rr1 = "{sample}/pre_process/{sample}.merged_R1.fastq.gz",
        rr2 = "{sample}/pre_process/{sample}.merged_R2.fastq.gz"
    output:
        r1 = temp("{sample}/pre_process/{sample}.trimmed_R1.fastq"),
        r2 = temp("{sample}/pre_process/{sample}.trimmed_R2.fastq"),
        s = "{sample}/pre_process/{sample}.trimmed_singles.fastq"
    singularity: "docker://quay.io/biocontainers/sickle-trim:1.33--ha92aebf_4"
    conda: "envs/sickle.yml"
    shell: "sickle pe -f {input.r1} -r {input.r2} -t sanger -o {output.r1} "
           "-p {output.r2} -s {output.s}"

rule cutadapt:
    """Clip fastq files"""
    input:
        r1 = "{sample}/pre_process/{sample}.trimmed_R1.fastq",
        r2 = "{sample}/pre_process/{sample}.trimmed_R2.fastq"
    output:
        r1 = temp("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2 = temp("{sample}/pre_process/{sample}.cutadapt_R2.fastq")
    singularity: "docker://quay.io/biocontainers/cutadapt:1.14--py36_0"
    conda: "envs/cutadapt.yml"
    shell: "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 1 -o {output.r1} "
           "{input.r1} -p {output.r2} {input.r2}"

rule align:
    """Align fastq files"""
    input:
        r1 = "{sample}/pre_process/{sample}.cutadapt_R1.fastq",
        r2 = "{sample}/pre_process/{sample}.cutadapt_R2.fastq",
        ref = REFERENCE,
        temp = ancient("tmp")
    params:
        rg = "@RG\\tID:{sample}_lib1\\tSM:{sample}\\tPL:ILLUMINA"
    output: temp("{sample}/bams/{sample}.sorted.bam")
    singularity: "docker://quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0"
    conda: "envs/bwa.yml"
    shell: "bwa mem -t 8 -R '{params.rg}' {input.ref} {input.r1} {input.r2} "
           "| picard -Xmx4G SortSam CREATE_INDEX=TRUE TMP_DIR={input.temp} "
           "INPUT=/dev/stdin OUTPUT={output} SORT_ORDER=coordinate"

rule markdup:
    """Mark duplicates in BAM file"""
    input:
        bam = "{sample}/bams/{sample}.sorted.bam",
        tmp = ancient("tmp")
    output:
        bam = "{sample}/bams/{sample}.markdup.bam",
        bai = "{sample}/bams/{sample}.markdup.bai",
        metrics = "{sample}/bams/{sample}.markdup.metrics"
    singularity: "docker://quay.io/biocontainers/picard:2.14--py36_0"
    conda: "envs/picard.yml"
    shell: "picard -Xmx4G MarkDuplicates CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "INPUT={input.bam} OUTPUT={output.bam} "
           "METRICS_FILE={output.metrics} "
           "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

rule bai:
    """Copy bai files as some genome browsers can only .bam.bai files"""
    input:
        bai = "{sample}/bams/{sample}.markdup.bai"
    output:
        bai = "{sample}/bams/{sample}.markdup.bam.bai"
    singularity: "docker://debian:buster-slim"
    shell: "cp {input.bai} {output.bai}"

rule baserecal:
    """Base recalibrated BAM files"""
    input:
        bam = "{sample}/bams/{sample}.markdup.bam",
        gatk = GATK,
        ref = REFERENCE,
        dbsnp = DBSNP,
        one1kg = ONETHOUSAND,
        hapmap = HAPMAP
    output:
        grp = "{sample}/bams/{sample}.baserecal.grp"
    params:
        knownsites = known_sites_argument
    singularity: "docker://quay.io/biocontainers/gatk:3.7--py36_1"
    conda: "envs/gatk.yml"
    shell: "java -XX:ParallelGCThreads=1 -jar {input.gatk} -T "
           "BaseRecalibrator -I {input.bam} -o {output.grp} -nct 8 "
           "-R {input.ref} -cov ReadGroupCovariate -cov QualityScoreCovariate "
           "-cov CycleCovariate -cov ContextCovariate {params.knownsites}"

rule gvcf_scatter:
    """Run HaplotypeCaller in GVCF mode by chunk"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        bqsr="{sample}/bams/{sample}.baserecal.grp",
        dbsnp=DBSNP,
        ref=REFERENCE,
        gatk=GATK
    params:
        chunk="{chunk}"
    output:
        gvcf=temp("{sample}/vcf/{sample}.{chunk}.part.vcf.gz"),
        gvcf_tbi=temp("{sample}/vcf/{sample}.{chunk}.part.vcf.gz.tbi")
    singularity: "docker://quay.io/biocontainers/gatk:3.7--py36_1"
    conda: "envs/gatk.yml"
    shell: "java -jar -Xmx4G -XX:ParallelGCThreads=1 {input.gatk} "
           "-T HaplotypeCaller -ERC GVCF -I "
           "{input.bam} -R {input.ref} -D {input.dbsnp} "
           "-L '{params.chunk}' -o '{output.gvcf}' "
           "-variant_index_type LINEAR -variant_index_parameter 128000 "
           "-BQSR {input.bqsr}"


rule gvcf_chunkfile:
    """
    Create simple text file with paths to chunks for GVCF.
    
    This uses a run directive in stead of a shell directive because
    the amount of chunks may be so large the shell would error out with
    an "argument list too long" error. 
    See https://unix.stackexchange.com/a/120842 for more info
    
    This also means this rule lives outside of singularity/conda and is
    executed in snakemake's own environment. 
    """
    params:
        chunkfiles = expand("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz",
                            chunk=CHUNKS)
    output:
        file = "{sample}/vcf/chunkfile.txt"
    run:
        with open(output.file, "w") as ohandle:
            for filename in params.chunkfiles:
                corrected = filename.format(sample=wildcards.sample)
                ohandle.write(corrected + "\n")

rule gvcf_gather:
    """Gather all GVCF scatters"""
    input:
        gvcfs = expand("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz",
                       chunk=CHUNKS),
        chunkfile = "{sample}/vcf/chunkfile.txt"
    output:
        gvcf = "{sample}/vcf/{sample}.g.vcf.gz"
    conda: "envs/bcftools.yml"
    singularity: "docker://quay.io/biocontainers/bcftools:1.9--ha228f0b_4"
    shell: "bcftools concat -f {input.chunkfile} -n > {output.gvcf}"


rule gvcf_gather_tbi:
    """Index GVCF"""
    input:
        gvcf = "{sample}/vcf/{sample}.g.vcf.gz"
    output:
        tbi = "{sample}/vcf/{sample}.g.vcf.gz.tbi"
    conda: "envs/tabix.yml"
    singularity: "docker://quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    shell: "tabix -pvcf {input.gvcf}"


rule genotype_scatter:
    """Run GATK's GenotypeGVCFs by chunk"""
    input:
        gvcfs = expand("{sample}/vcf/{sample}.g.vcf.gz", sample=SAMPLES),
        tbis = expand("{sample}/vcf/{sample}.g.vcf.gz.tbi",
                      sample=SAMPLES),
        ref=REFERENCE,
        gatk=GATK
    params:
        li=" -V ".join(expand("{sample}/vcf/{sample}.g.vcf.gz",
                              sample=SAMPLES)),
        chunk="{chunk}"
    output:
        vcf=temp("multisample/genotype.{chunk}.part.vcf.gz"),
        vcf_tbi=temp("multisample/genotype.{chunk}.part.vcf.gz.tbi")
    singularity: "docker://quay.io/biocontainers/gatk:3.7--py36_1"
    conda: "envs/gatk.yml"
    shell: "java -jar -Xmx15G -XX:ParallelGCThreads=1 {input.gatk} -T "
           "GenotypeGVCFs -R {input.ref} "
           "-V {params.li} -L '{params.chunk}' -o '{output.vcf}'"


rule genotype_chunkfile:
    """
    Create simple text file with paths to chunks for genotyping
    
    This uses a run directive in stead of a shell directive because
    the amount of chunks may be so large the shell would error out with
    an "argument list too long" error. 
    See https://unix.stackexchange.com/a/120842 for more info
    
    This also means this rule lives outside of singularity/conda and is
    executed in snakemake's own environment. 
    """
    params:
        vcfs = expand("multisample/genotype.{chunk}.part.vcf.gz",
                      chunk=CHUNKS)
    output:
        file = "multisample/chunkfile.txt"
    run:
        with open(output.file, "w") as ohandle:
            for filename in params.vcfs:
                ohandle.write(filename + "\n")


rule genotype_gather:
    """Gather all genotyping VCFs"""
    input:
        vcfs = expand("multisample/genotype.{chunk}.part.vcf.gz",
                      chunk=CHUNKS),
        chunkfile = "multisample/chunkfile.txt"
    output:
        vcf = "multisample/genotyped.vcf.gz"
    conda: "envs/bcftools.yml"
    singularity: "docker://quay.io/biocontainers/bcftools:1.9--ha228f0b_4"
    shell: "bcftools concat -f {input.chunkfile} -n > {output.vcf}"


rule genotype_gather_tbi:
    """Index genotyped vcf file"""
    input:
        vcf = "multisample/genotyped.vcf.gz"
    output:
        tbi = "multisample/genotyped.vcf.gz.tbi"
    conda: "envs/tabix.yml"
    singularity: "docker://quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    shell: "tabix -pvcf {input.vcf}"


rule split_vcf:
    """Split multisample VCF in single samples"""
    input:
        vcf="multisample/genotyped.vcf.gz",
        tbi = "multisample/genotyped.vcf.gz.tbi",
        gatk=GATK,
        ref=REFERENCE
    params:
        s="{sample}"
    output:
        splitted="{sample}/vcf/{sample}_single.vcf.gz"
    singularity: "docker://quay.io/biocontainers/gatk:3.7--py36_1"
    conda: "envs/gatk.yml"
    shell: "java -Xmx15G -XX:ParallelGCThreads=1 -jar {input.gatk} "
           "-T SelectVariants -sn {params.s} -env -R {input.ref} -V "
           "{input.vcf} -o {output.splitted}"


## bam metrics

rule mapped_num:
    """Calculate number of mapped reads"""
    input:
        bam="{sample}/bams/{sample}.sorted.bam"
    output:
        num="{sample}/bams/{sample}.mapped.num"
    singularity: "docker://quay.io/biocontainers/samtools:1.6--he673b24_3"
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | wc -l > {output.num}"


rule mapped_basenum:
    """Calculate number of mapped bases"""
    input:
        bam="{sample}/bams/{sample}.sorted.bam"
    output:
        num="{sample}/bams/{sample}.mapped.basenum"
    singularity: "docker://quay.io/biocontainers/samtools:1.6--he673b24_3"
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | cut -f10 | wc -c > {output.num}"


rule unique_num:
    """Calculate number of unique reads"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam"
    output:
        num="{sample}/bams/{sample}.unique.num"
    singularity: "docker://quay.io/biocontainers/samtools:1.6--he673b24_3"
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | wc -l > {output.num}"


rule usable_basenum:
    """Calculate number of bases on unique reads"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam"
    output:
        num="{sample}/bams/{sample}.usable.basenum"
    singularity: "docker://quay.io/biocontainers/samtools:1.6--he673b24_3"
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | cut -f10 | wc -c > "
           "{output.num}"


## fastqc

rule fastqc_raw:
    """
    Run fastqc on raw fastq files
    NOTE: singularity version uses 0.11.7 in stead of 0.11.5 due to 
    perl missing in the container of 0.11.5
    """
    input:
        r1=get_r1,
        r2=get_r2
    params:
        odir="{sample}/pre_process/raw_fastqc"
    output:
        aux="{sample}/pre_process/raw_fastqc/.done.txt"
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    conda: "envs/fastqc.yml"
    shell: "fastqc --nogroup -o {params.odir} {input.r1} {input.r2} "
           "&& echo 'done' > {output.aux}"


rule fastqc_merged:
    """
    Run fastqc on merged fastq files    
    NOTE: singularity version uses 0.11.7 in stead of 0.11.5 due to 
    perl missing in the container of 0.11.5
    """
    input:
        r1="{sample}/pre_process/{sample}.merged_R1.fastq.gz",
        r2="{sample}/pre_process/{sample}.merged_R2.fastq.gz",
        fq=fqsc
    params:
        odir="{sample}/pre_process/merged_fastqc"
    output:
        r1="{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip",
        r2="{sample}/pre_process/merged_fastqc/{sample}.merged_R2_fastqc.zip"
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    conda: "envs/fastqc.yml"
    shell: "bash {input.fq} {input.r1} {input.r2} "
           "{output.r1} {output.r2} {params.odir}"


rule fastqc_postqc:
    """
    Run fastqc on fastq files post pre-processing
    NOTE: singularity version uses 0.11.7 in stead of 0.11.5 due to 
    perl missing in the container of 0.11.5     
    """
    input:
        r1="{sample}/pre_process/{sample}.cutadapt_R1.fastq",
        r2="{sample}/pre_process/{sample}.cutadapt_R2.fastq",
        fq=fqsc
    params:
        odir="{sample}/pre_process/postqc_fastqc"
    output:
        r1="{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip",
        r2="{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R2_fastqc.zip"
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    conda: "envs/fastqc.yml"
    shell: "bash {input.fq} {input.r1} {input.r2} "
           "{output.r1} {output.r2} {params.odir}"


## fastq-count

rule fqcount_preqc:
    """Calculate number of reads and bases before pre-processing"""
    input:
        r1="{sample}/pre_process/{sample}.merged_R1.fastq.gz",
        r2="{sample}/pre_process/{sample}.merged_R2.fastq.gz"
    output:
        "{sample}/pre_process/{sample}.preqc_count.json"
    singularity: "docker://quay.io/biocontainers/fastq-count:0.1.0--h14c3975_0"
    conda: "envs/fastq-count.yml"
    shell: "fastq-count {input.r1} {input.r2} > {output}"


rule fqcount_postqc:
    """Calculate number of reads and bases after pre-processing"""
    input:
        r1="{sample}/pre_process/{sample}.cutadapt_R1.fastq",
        r2="{sample}/pre_process/{sample}.cutadapt_R2.fastq"
    output:
        "{sample}/pre_process/{sample}.postqc_count.json"
    singularity: "docker://quay.io/biocontainers/fastq-count:0.1.0--h14c3975_0"
    conda: "envs/fastq-count.yml"
    shell: "fastq-count {input.r1} {input.r2} > {output}"


# fastqc stats
rule fastqc_stats:
    """Collect fastq stats for a sample in json format"""
    input:
        preqc_r1="{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip",
        preqc_r2="{sample}/pre_process/merged_fastqc/{sample}.merged_R2_fastqc.zip",
        postqc_r1="{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip",
        postqc_r2="{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R2_fastqc.zip",
        sc=fqpy
    singularity: "docker://python:3.6-slim"
    conda: "envs/collectstats.yml"
    output:
        "{sample}/pre_process/fastq_stats.json"
    shell: "python {input.sc} --preqc-r1 {input.preqc_r1} "
           "--preqc-r2 {input.preqc_r2} "
           "--postqc-r1 {input.postqc_r1} "
           "--postqc-r2 {input.postqc_r2} > {output}"

## coverages

rule covstats:
    """Calculate coverage statistics on bam file"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        genome="current.genome",
        covpy=covpy,
        bed=get_bedpath
    params:
        subt="Sample {sample}"
    output:
        covj="{sample}/coverage/{bed}.covstats.json",
        covp="{sample}/coverage/{bed}.covstats.png"
    singularity: "docker://quay.io/biocontainers/mulled-v2-3251e6c49d800268f0bc575f28045ab4e69475a6:4ce073b219b6dabb79d154762a9b67728c357edb-0"
    conda: "envs/covstat.yml"
    shell: "bedtools coverage -sorted -g {input.genome} -a {input.bed} "
           "-b {input.bam} -d  | python {input.covpy} - --plot {output.covp} "
           "--title 'Targets coverage' --subtitle '{params.subt}' "
           "> {output.covj}"


rule vtools_coverage:
    """Calculate coverage statistics per transcript"""
    input:
        gvcf="{sample}/vcf/{sample}.g.vcf.gz",
        tbi = "{sample}/vcf/{sample}.g.vcf.gz.tbi",
        ref=get_refflatpath
    output:
        tsv="{sample}/coverage/{ref}.coverages.tsv"
    singularity: "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
    conda: "envs/vcfstats.yml"
    shell: "vtools-gcoverage -I {input.gvcf} -R {input.ref} > {output.tsv}"


## vcfstats

rule vcfstats:
    """Calculate vcf statistics"""
    input:
        vcf="multisample/genotyped.vcf.gz",
        tbi = "multisample/genotyped.vcf.gz.tbi"
    output:
        stats="multisample/vcfstats.json"
    singularity: "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
    conda: "envs/vcfstats.yml"
    shell: "vtools-stats -i {input.vcf} > {output.stats}"


## collection
if len(BASE_BEDS) >= 1:
    rule collectstats:
        """Collect all stats for a particular sample with beds"""
        input:
            preqc="{sample}/pre_process/{sample}.preqc_count.json",
            postq="{sample}/pre_process/{sample}.postqc_count.json",
            mnum="{sample}/bams/{sample}.mapped.num",
            mbnum="{sample}/bams/{sample}.mapped.basenum",
            unum="{sample}/bams/{sample}.unique.num",
            ubnum="{sample}/bams/{sample}.usable.basenum",
            fastqc="{sample}/pre_process/fastq_stats.json",
            cov=expand("{{sample}}/coverage/{bed}.covstats.json",
                       bed=BASE_BEDS),
            colpy=colpy
        params:
            sample_name="{sample}",
            fthresh=FEMALE_THRESHOLD
        output:
            "{sample}/{sample}.stats.json"
        singularity: "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
        conda: "envs/collectstats.yml"
        shell: "python {input.colpy} --sample-name {params.sample_name} "
               "--pre-qc-fastq {input.preqc} --post-qc-fastq {input.postq} "
               "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} "
               "--unique-num {input.unum} --usable-basenum {input.ubnum} "
               "--female-threshold {params.fthresh} "
               "--fastqc-stats {input.fastqc} {input.cov} > {output}"
else:
    rule collectstats:
        """Collect all stats for a particular sample without beds"""
        input:
            preqc = "{sample}/pre_process/{sample}.preqc_count.json",
            postq = "{sample}/pre_process/{sample}.postqc_count.json",
            mnum = "{sample}/bams/{sample}.mapped.num",
            mbnum = "{sample}/bams/{sample}.mapped.basenum",
            unum = "{sample}/bams/{sample}.unique.num",
            ubnum = "{sample}/bams/{sample}.usable.basenum",
            fastqc="{sample}/pre_process/fastq_stats.json",
            colpy = colpy
        params:
            sample_name = "{sample}",
            fthresh = FEMALE_THRESHOLD
        output:
            "{sample}/{sample}.stats.json"
        singularity: "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
        conda: "envs/collectstats.yml"
        shell: "python {input.colpy} --sample-name {params.sample_name} "
               "--pre-qc-fastq {input.preqc} --post-qc-fastq {input.postq} "
               "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} "
               "--unique-num {input.unum} --usable-basenum {input.ubnum} "
               "--female-threshold {params.fthresh} "
               "--fastqc-stats {input.fastqc}  > {output}"

rule merge_stats:
    """Merge all stats of all samples"""
    input:
        cols=expand("{sample}/{sample}.stats.json", sample=SAMPLES),
        vstat="multisample/vcfstats.json",
        mpy=mpy
    output:
        stats="stats.json"
    singularity: "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
    conda: "envs/collectstats.yml"
    shell: "python {input.mpy} --vcfstats {input.vstat} {input.cols} "
           "> {output.stats}"


rule stats_tsv:
    """Convert stats.json to tsv"""
    input:
        stats="stats.json",
        sc=tsvpy
    output:
        stats="stats.tsv"
    singularity: "docker://python:3.6-slim"
    conda: "envs/collectstats.yml"
    shell: "python {input.sc} -i {input.stats} > {output.stats}"


rule multiqc:
    """
    Create multiQC report
    Depends on stats.tsv to forcefully run at end of pipeline
    """
    input:
        stats="stats.tsv"
    params:
        odir=".",
        rdir="multiqc_report"
    output:
        report="multiqc_report/multiqc_report.html"
    singularity: "docker://quay.io/biocontainers/multiqc:1.5--py36_0"
    conda: "envs/multiqc.yml"
    shell: "multiqc -f -o {params.rdir} {params.odir} || touch {output.report}"
