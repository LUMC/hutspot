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
import itertools

from pyfaidx import Fasta

REFERENCE = config.get("REFERENCE")
if REFERENCE is None:
    raise ValueError("You must set --config REFERENCE=<path>")
if not Path(REFERENCE).exists():
    raise FileNotFoundError("Reference file {0} "
                            "does not exist.".format(REFERENCE))
DBSNP = config.get("DBSNP")
if DBSNP is None:
    raise ValueError("You must set --config DBSNP=<path>")
if not Path(DBSNP).exists():
    raise FileNotFoundError("{0} does not exist".format(DBSNP))

SCONFIG = config.get("SAMPLE_CONFIG")
if SCONFIG is None:
    raise ValueError("You must set --config SAMPLE_CONFIG=<path>")
if not Path(SCONFIG).exists():
    raise FileNotFoundError("{0} does not exist".format(SCONFIG))

KNOWN_SITES = config.get("KNOWN_SITES")
if KNOWN_SITES is None:
    raise ValueError("You must set --config KNOWN_SITES=<path>")

# The user can specify multiple known sites
KNOWN_SITES = KNOWN_SITES.split(',')
for filename in KNOWN_SITES:
    if not Path(filename).exists():
        raise FileNotFoundError("{0} does not exist".format(filename))


# these are all optional
BED = config.get("BED", "")  # comma-separated list of BED files
REFFLAT = config.get("REFFLAT", "")  # comma-separated list of refFlat files
FEMALE_THRESHOLD = config.get("FEMALE_THRESHOLD", 0.6)
MAX_BASES = config.get("MAX_BASES", "")

# Generate the input string for basrecalibration
BSQR_known_sites = ''
for argument, filename in zip(itertools.repeat('-knownSites'), KNOWN_SITES):
    BSQR_known_sites +=' {} {}'.format(argument, filename)

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
pywc = fsrc_dir("src", "pywc.py")

# sample config parsing
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

containers = {
    "bcftools": "docker://quay.io/biocontainers/bcftools:1.9--ha228f0b_4",
    "bedtools-2.26-python-2.7": "docker://quay.io/biocontainers/mulled-v2-3251e6c49d800268f0bc575f28045ab4e69475a6:4ce073b219b6dabb79d154762a9b67728c357edb-0",
    "biopet-scatterregions": "docker://quay.io/biocontainers/biopet-scatterregions:0.2--0",
    "bwa-0.7.17-picard-2.18.7": "docker://quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0",
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:1.14--py36_0",
    "debian": "docker://debian:buster-slim",
    "fastq-count": "docker://quay.io/biocontainers/fastq-count:0.1.0--h14c3975_0",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.11.7--4",
    "gatk": "docker://broadinstitute/gatk3:3.7-0",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.5--py36_0",
    "picard-2.14": "docker://quay.io/biocontainers/picard:2.14--py36_0",
    "python3": "docker://python:3.6-slim",
    "samtools-1.7-python-3.6": "docker://quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:1abf1824431ec057c7d41be6f0c40e24843acde4-0",
    "seqtk-1.2-jq-1.6": "docker://quay.io/biocontainers/mulled-v2-13686261ac0aa5682c680670ff8cda7b09637943:d143450dec169186731bb4df6f045a3c9ee08eb6-0",
    "sickle": "docker://quay.io/biocontainers/sickle-trim:1.33--ha92aebf_4",
    "tabix": "docker://quay.io/biocontainers/tabix:0.2.6--ha92aebf_0",
    "vtools": "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
}

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


rule all:
    input:
        report="multiqc_report/multiqc_report.html",
        bais=expand("{sample}/bams/{sample}.markdup.bam.bai",sample=SAMPLES),
        vcfs=expand("{sample}/vcf/{sample}.vcf.gz",sample=SAMPLES),
        stats=metrics()


rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output: directory("tmp")
    singularity: containers["debian"]
    shell: "mkdir -p {output}"

rule genome:
    """Create genome file as used by bedtools"""
    input: REFERENCE
    output: "current.genome"
    singularity: containers["debian"]
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule merge_r1:
    """Merge all forward fastq files into one"""
    input: get_r1
    output: temp("{sample}/pre_process/{sample}.merged_R1.fastq.gz")
    singularity: containers["debian"]
    shell: "cat {input} > {output}"

rule merge_r2:
    """Merge all reverse fastq files into one"""
    input: get_r2
    output: temp("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    singularity: containers["debian"]
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
    singularity: containers["seqtk-1.2-jq-1.6"]
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
    singularity: containers["seqtk-1.2-jq-1.6"]
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
    singularity: containers["sickle"]
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
    singularity: containers["cutadapt"]
    shell: "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 1 -o {output.r1} "
           "{input.r1} -p {output.r2} {input.r2}"

rule align:
    """Align fastq files"""
    input:
        r1 = "{sample}/pre_process/{sample}.cutadapt_R1.fastq",
        r2 = "{sample}/pre_process/{sample}.cutadapt_R2.fastq",
        ref = REFERENCE,
        tmp = ancient("tmp")
    params:
        rg = "@RG\\tID:{sample}_lib1\\tSM:{sample}\\tPL:ILLUMINA"
    output: temp("{sample}/bams/{sample}.sorted.bam")
    singularity: containers["bwa-0.7.17-picard-2.18.7"]
    shell: "bwa mem -t 8 -R '{params.rg}' {input.ref} {input.r1} {input.r2} "
           "| picard -Xmx4G -Djava.io.tmpdir={input.tmp} SortSam "
           "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
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
    singularity: containers["picard-2.14"]
    shell: "picard -Xmx4G -Djava.io.tmpdir={input.tmp} MarkDuplicates "
           "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "INPUT={input.bam} OUTPUT={output.bam} "
           "METRICS_FILE={output.metrics} "
           "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

rule bai:
    """Copy bai files as some genome browsers can only .bam.bai files"""
    input:
        bai = "{sample}/bams/{sample}.markdup.bai"
    output:
        bai = "{sample}/bams/{sample}.markdup.bam.bai"
    singularity: containers["debian"]
    shell: "cp {input.bai} {output.bai}"

rule baserecal:
    """Base recalibrated BAM files"""
    input:
        bam = "{sample}/bams/{sample}.markdup.bam",
        ref = REFERENCE,
    output:
        grp = "{sample}/bams/{sample}.baserecal.grp"
    params:
        known_sites = BSQR_known_sites
    singularity: containers["gatk"]
    shell: "java -XX:ParallelGCThreads=1 -jar /usr/GenomeAnalysisTK.jar -T "
           "BaseRecalibrator -I {input.bam} -o {output.grp} -nct 8 "
           "-R {input.ref} -cov ReadGroupCovariate -cov QualityScoreCovariate "
           "-cov CycleCovariate -cov ContextCovariate {params.known_sites}"

checkpoint scatterregions:
    """Scatter the reference genome"""
    input:
        ref = REFERENCE,
    output:
        directory("scatter")
    singularity: containers["biopet-scatterregions"]
    shell: "mkdir -p {output} && "
           "biopet-scatterregions "
           "--referenceFasta {input.ref} --scatterSize 1000000000 "
           "--outputDir scatter"

rule gvcf_scatter:
    """Run HaplotypeCaller in GVCF mode by chunk"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        bqsr="{sample}/bams/{sample}.baserecal.grp",
        dbsnp=DBSNP,
        ref=REFERENCE,
        region="scatter/scatter-{chunk}.bed"
    output:
        gvcf=temp("{sample}/vcf/{sample}.{chunk}.g.vcf.gz"),
        gvcf_tbi=temp("{sample}/vcf/{sample}.{chunk}.g.vcf.gz.tbi")
    wildcard_constraints:
        chunk="[0-9]+"
    singularity: containers["gatk"]
    shell: "java -jar -Xmx4G -XX:ParallelGCThreads=1 /usr/GenomeAnalysisTK.jar "
           "-T HaplotypeCaller -ERC GVCF -I "
           "{input.bam} -R {input.ref} -D {input.dbsnp} "
           "-L '{input.region}' -o '{output.gvcf}' "
           "-variant_index_type LINEAR -variant_index_parameter 128000 "
           "-BQSR {input.bqsr}"


rule genotype_scatter:
    """Run GATK's GenotypeGVCFs by chunk"""
    input:
        gvcf = "{sample}/vcf/{sample}.{chunk}.g.vcf.gz",
        tbi = "{sample}/vcf/{sample}.{chunk}.g.vcf.gz.tbi",
        ref=REFERENCE
    output:
        vcf="{sample}/vcf/{sample}.{chunk}.vcf.gz",
        vcf_tbi="{sample}/vcf/{sample}.{chunk}.vcf.gz.tbi"
    wildcard_constraints:
        chunk="[0-9]+"
    singularity: containers["gatk"]
    shell: "java -jar -Xmx15G -XX:ParallelGCThreads=1 /usr/GenomeAnalysisTK.jar -T "
           "GenotypeGVCFs -R {input.ref} "
           "-V {input.gvcf} -o '{output.vcf}'"


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.vcf.gz",
           i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)


rule genotype_gather:
    """Gather all genotyping VCFs"""
    input:
        vcfs = aggregate_input
    output:
        vcf = "{sample}/vcf/{sample}.vcf.gz"
    singularity: containers["bcftools"]
    shell: "bcftools concat {input.vcfs} --allow-overlaps --output {output.vcf} "
           "--output-type z"


rule genotype_gather_tbi:
    """Index genotyped vcf file"""
    input:
        vcf = "{sample}/vcf/{sample}.vcf.gz"
    output:
        tbi = "{sample}/vcf/{sample}.vcf.gz.tbi"
    singularity: containers["tabix"]
    shell: "tabix -pvcf {input.vcf}"


## bam metrics
rule mapped_reads_bases:
    """Calculate number of mapped reads"""
    input:
        bam="{sample}/bams/{sample}.sorted.bam",
        pywc=pywc
    output:
        reads="{sample}/bams/{sample}.mapped.num",
        bases="{sample}/bams/{sample}.mapped.basenum"
    singularity: containers["samtools-1.7-python-3.6"]
    shell: "samtools view -F 4 {input.bam} | cut -f 10 | python {input.pywc} "
           "--reads {output.reads} --bases {output.bases}"


rule unique_reads_bases:
    """Calculate number of unique reads"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        pywc=pywc
    output:
        reads="{sample}/bams/{sample}.unique.num",
        bases="{sample}/bams/{sample}.usable.basenum"
    singularity: containers["samtools-1.7-python-3.6"]
    shell: "samtools view -F 4 -F 1024 {input.bam} | cut -f 10 | python {input.pywc} "
           "--reads {output.reads} --bases {output.bases}"


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
    singularity: containers["fastqc"]
    shell: "fastqc --threads 4 --nogroup -o {params.odir} {input.r1} {input.r2} "
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
    singularity: containers["fastqc"]
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
    singularity: containers["fastqc"]
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
    singularity: containers["fastq-count"]
    shell: "fastq-count {input.r1} {input.r2} > {output}"


rule fqcount_postqc:
    """Calculate number of reads and bases after pre-processing"""
    input:
        r1="{sample}/pre_process/{sample}.cutadapt_R1.fastq",
        r2="{sample}/pre_process/{sample}.cutadapt_R2.fastq"
    output:
        "{sample}/pre_process/{sample}.postqc_count.json"
    singularity: containers["fastq-count"]
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
    singularity: containers["python3"]
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
    singularity: containers["bedtools-2.26-python-2.7"]
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
    singularity: containers["vtools"]
    shell: "vtools-gcoverage -I {input.gvcf} -R {input.ref} > {output.tsv}"


## vcfstats

rule vcfstats:
    """Calculate vcf statistics"""
    input:
        vcf="{sample}/vcf/{sample}.vcf.gz",
        tbi = "{sample}/vcf/{sample}.vcf.gz.tbi"
    output:
        stats="{sampel}/vcf/{sample}.vcfstats.json"
    singularity: containers["vtools"]
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
        singularity: containers["vtools"]
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
        singularity: containers["vtools"]
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
        vstat=expand("{sample}/vcf/{sample}.vcfstats.json", sample=SAMPLES),
        mpy=mpy
    output:
        stats="stats.json"
    singularity: containers["vtools"]
    shell: "python {input.mpy} --vcfstats {input.vstat} {input.cols} "
           "> {output.stats}"


rule stats_tsv:
    """Convert stats.json to tsv"""
    input:
        stats="stats.json",
        sc=tsvpy
    output:
        stats="stats.tsv"
    singularity: containers["python3"]
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
    singularity: containers["multiqc"]
    shell: "multiqc -f -o {params.rdir} {params.odir} || touch {output.report}"
