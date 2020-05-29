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
import jsonschema
from functools import partial
from os.path import join, basename
from pathlib import Path
import itertools

# Read the json schema
with open(srcdir('config/schema.json'), 'rt') as fin:
    schema = json.load(fin)

# Validate the config against the schema
try:
    jsonschema.validate(config, schema)
except jsonschema.ValidationError as e:
    raise jsonschema.ValidationError(f'Invalid CONFIG_JSON: {e}')


# Set default values
def set_default(key, value):
    """ Set default config values """
    if key not in config:
        config[key] = value

# Set the default config
set_default('scatter_size', 1000000000)
set_default('female_threshold', 0.6)

# Set the script paths
set_default("covstats", srcdir("src/covstats.py"))
set_default("collect_stats", srcdir("src/collect_stats.py"))
set_default("merge_stats", srcdir("src/merge_stats.py"))
set_default("stats_to_tsv", srcdir("src/stats_to_tsv.py"))
set_default("py_wordcount", srcdir("src/pywc.py"))
set_default("cutadapt_summary", srcdir("src/cutadapt_summary.py"))

containers = {
    "bcftools": "docker://quay.io/biocontainers/bcftools:1.9--ha228f0b_4",
    "bedtools-2.26-python-2.7": "docker://quay.io/biocontainers/mulled-v2-3251e6c49d800268f0bc575f28045ab4e69475a6:4ce073b219b6dabb79d154762a9b67728c357edb-0",
    "biopet-scatterregions": "docker://quay.io/biocontainers/biopet-scatterregions:0.2--0",
    "bwa-0.7.17-picard-2.18.7": "docker://quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0",
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:2.9--py37h516909a_0",
    "debian": "docker://debian:buster-slim",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.11.7--4",
    "gatk": "docker://broadinstitute/gatk3:3.7-0",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.8--py_2",
    "picard": "docker://quay.io/biocontainers/picard:2.22.8--0",
    "python3": "docker://python:3.6-slim",
    "samtools-1.7-python-3.6": "docker://quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:1abf1824431ec057c7d41be6f0c40e24843acde4-0",
    "vtools": "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
}

def get_forward(wildcards):
    """ Get the forward fastq file from the config """
    return (
        config["samples"][wildcards.sample]["read_groups"]
                [wildcards.read_group]["R1"]
    )

def get_reverse(wildcards):
    """ Get the reverse fastq file from the config """
    return (
        config["samples"][wildcards.sample]["read_groups"]
            [wildcards.read_group]["R2"]
    )

def get_readgroup(wildcards):
    return config["samples"][wildcards.sample]["read_groups"]

def get_readgroup_per_sample():
    for sample in config["samples"]:
        for rg in config["samples"][sample]["read_groups"]:
            yield rg, sample


def coverage_stats(wildcards):
    files = expand("{sample}/coverage/refFlat_coverage.tsv", sample=config["samples"])
    return files if "refflat" in config else []

rule all:
    input:
        multiqc = "multiqc_report/multiqc_report.html",
        stats = "stats.json",
        stats_tsv = "stats.tsv",
        bam = expand("{sample}/bams/{sample}.bam", sample=config["samples"]),
        vcfs = expand("{sample}/vcf/{sample}.vcf.gz", sample=config["samples"]),
        vcf_tbi = expand("{sample}/vcf/{sample}.vcf.gz.tbi", sample=config["samples"]),
        gvcfs = expand("{sample}/vcf/{sample}.g.vcf.gz", sample=config["samples"]),
        gvcf_tbi = expand("{sample}/vcf/{sample}.g.vcf.gz.tbi", sample=config["samples"]),
        fastqc_raw = (f"{sample}/pre_process/raw-{sample}-{read_group}/" for read_group, sample in get_readgroup_per_sample()),
        fastqc_trimmed = (f"{sample}/pre_process/trimmed-{sample}-{read_group}/" for read_group, sample in get_readgroup_per_sample()),
        cutadapt = (f"{sample}/pre_process/{sample}-{read_group}.txt" for read_group, sample in get_readgroup_per_sample()),
        coverage_stats = coverage_stats,


rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output: directory("tmp")
    singularity: containers["debian"]
    shell: "mkdir -p {output}"

rule genome:
    """Create genome file as used by bedtools"""
    input: config["reference"]
    output: "current.genome"
    singularity: containers["debian"]
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule cutadapt:
    """Clip fastq files"""
    input:
        r1 = get_forward,
        r2 = get_reverse
    output:
        r1 = "{sample}/pre_process/{sample}-{read_group}_R1.fastq.gz",
        r2 = "{sample}/pre_process/{sample}-{read_group}_R2.fastq.gz",
    log:
        "{sample}/pre_process/{sample}-{read_group}.txt"
    singularity: containers["cutadapt"]
    shell: "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG "
           "--minimum-length 1 --quality-cutoff=20,20 "
           "--output {output.r1} --paired-output {output.r2} -Z "
           "{input.r1} {input.r2} "
           "--report=minimal > {log}"

rule align:
    """Align fastq files"""
    input:
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2,
        ref = config["reference"],
        tmp = ancient("tmp")
    params:
        rg = "@RG\\tID:{read_group}\\tSM:{sample}\\tPL:ILLUMINA"
    output: "{sample}/bams/{sample}-{read_group}.sorted.bam"
    singularity: containers["bwa-0.7.17-picard-2.18.7"]
    shell: "bwa mem -t 8 -R '{params.rg}' {input.ref} {input.r1} {input.r2} "
           "| picard -Xmx4G -Djava.io.tmpdir={input.tmp} SortSam "
           "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "INPUT=/dev/stdin OUTPUT={output} SORT_ORDER=coordinate"

def markdup_bam_input(wildcards):
    """ Generate the INPUT for each bam file """
    return ["INPUT={sample}/bams/{sample}-{read_group}.sorted.bam".format(sample=wildcards.sample,
            read_group=rg) for rg in get_readgroup(wildcards)]

rule markdup:
    """Mark duplicates in BAM file"""
    input:
        bam = lambda wildcards:
        ("{sample}/bams/{sample}-{read_group}.sorted.bam".format(sample=wildcards.sample,
        read_group=rg) for rg in get_readgroup(wildcards)),
        tmp = ancient("tmp")
    output:
        bam = "{sample}/bams/{sample}.bam",
        bai = "{sample}/bams/{sample}.bai",
        metrics = "{sample}/bams/{sample}.metrics"
    params:
        bams=markdup_bam_input
    singularity: containers["picard"]
    shell: "picard -Xmx4G -Djava.io.tmpdir={input.tmp} MarkDuplicates "
           "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "{params.bams} OUTPUT={output.bam} "
           "METRICS_FILE={output.metrics} "
           "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

rule bai:
    """Copy bai files as some genome browsers can only .bam.bai files"""
    input:
        bai = rules.markdup.output.bai
    output:
        bai = "{sample}/bams/{sample}.bam.bai"
    singularity: containers["debian"]
    shell: "cp {input.bai} {output.bai}"


def bqsr_bam_input(wildcards):
    """ Generate the bam input string for each read group for BQSR"""
    return " ".join(["-I {sample}/bams/{sample}-{read_group}.sorted.bam".format(sample=wildcards.sample,
            read_group=rg) for rg in get_readgroup(wildcards)])

rule baserecal:
    """Base recalibrated BAM files"""
    input:
        bam = lambda wildcards:
        ("{sample}/bams/{sample}-{read_group}.sorted.bam".format(sample=wildcards.sample,
        read_group=rg) for rg in get_readgroup(wildcards)),
        ref = config["reference"],
        vcfs = config["known_sites"]
    output: "{sample}/bams/{sample}.baserecal.grp"
    params:
        known_sites = " ".join(expand("-knownSites {vcf}", vcf=config["known_sites"])),
        bams = bqsr_bam_input
    singularity: containers["gatk"]
    shell: "java -XX:ParallelGCThreads=1 -jar /usr/GenomeAnalysisTK.jar -T "
           "BaseRecalibrator {params.bams} -o {output} -nct 8 "
           "-R {input.ref} -cov ReadGroupCovariate -cov QualityScoreCovariate "
           "-cov CycleCovariate -cov ContextCovariate {params.known_sites}"

checkpoint scatterregions:
    """Scatter the reference genome"""
    input:
        ref = config["reference"],
    params:
        size = config['scatter_size']
    output:
        directory("scatter")
    singularity: containers["biopet-scatterregions"]
    shell: "mkdir -p {output} && "
           "biopet-scatterregions -Xmx24G "
           "--referenceFasta {input.ref} --scatterSize {params.size} "
           "--outputDir scatter"


rule gvcf_scatter:
    """Run HaplotypeCaller in GVCF mode by chunk"""
    input:
        bam = rules.markdup.output.bam,
        bqsr = rules.baserecal.output,
        dbsnp = config["dbsnp"],
        ref = config["reference"],
        region = "scatter/scatter-{chunk}.bed"
    output:
        gvcf = temp("{sample}/vcf/{sample}.{chunk}.g.vcf.gz"),
        gvcf_tbi = temp("{sample}/vcf/{sample}.{chunk}.g.vcf.gz.tbi")
    wildcard_constraints:
        chunk = "[0-9]+"
    singularity: containers["gatk"]
    shell: "java -jar -Xmx4G -XX:ParallelGCThreads=1 /usr/GenomeAnalysisTK.jar "
           "-T HaplotypeCaller -ERC GVCF -I "
           "{input.bam} -R {input.ref} -D {input.dbsnp} "
           "-L '{input.region}' -o '{output.gvcf}' "
           "-variant_index_type LINEAR -variant_index_parameter 128000 "
           "-BQSR {input.bqsr} "
           "--GVCFGQBands 20 --GVCFGQBands 40 --GVCFGQBands 60 "
           "--GVCFGQBands 80 --GVCFGQBands 100 "

def aggregate_gvcf(wildcards):
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.g.vcf.gz",
           i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

def aggregate_gvcf_tbi(wildcards):
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.g.vcf.gz.tbi",
           i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

rule gvcf_gather:
    """ Join the gvcf files together """
    input:
        gvcfs = aggregate_gvcf,
        tbis = aggregate_gvcf_tbi,
    output:
        gvcf = "{sample}/vcf/{sample}.g.vcf.gz",
        gvcf_tbi = "{sample}/vcf/{sample}.g.vcf.gz.tbi"
    singularity: containers["bcftools"]
    shell: "bcftools concat {input.gvcfs} --allow-overlaps --output {output.gvcf} "
           "--output-type z && "
           "bcftools index --tbi --output-file {output.gvcf_tbi} {output.gvcf}"

rule genotype_scatter:
    """Run GATK's GenotypeGVCFs by chunk"""
    input:
        gvcf = rules.gvcf_scatter.output.gvcf,
        tbi = rules.gvcf_scatter.output.gvcf_tbi,
        ref= config["reference"]
    output:
        vcf = temp("{sample}/vcf/{sample}.{chunk}.vcf.gz"),
        vcf_tbi = temp("{sample}/vcf/{sample}.{chunk}.vcf.gz.tbi")
    wildcard_constraints:
        chunk="[0-9]+"
    singularity: containers["gatk"]
    shell: "java -jar -Xmx15G -XX:ParallelGCThreads=1 /usr/GenomeAnalysisTK.jar -T "
           "GenotypeGVCFs -R {input.ref} "
           "-V {input.gvcf} -o '{output.vcf}'"


def aggregate_vcf(wildcards):
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.vcf.gz",
           i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)


def aggregate_vcf_tbi(wildcards):
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.vcf.gz.tbi",
           i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)


rule genotype_gather:
    """Gather all genotyping VCFs"""
    input:
        vcfs = aggregate_vcf,
        vcfs_tbi = aggregate_vcf_tbi
    output:
        vcf = "{sample}/vcf/{sample}.vcf.gz",
        vcf_tbi = "{sample}/vcf/{sample}.vcf.gz.tbi"
    singularity: containers["bcftools"]
    shell: "bcftools concat {input.vcfs} --allow-overlaps --output {output.vcf} "
           "--output-type z && "
           "bcftools index --tbi --output-file {output.vcf_tbi} {output.vcf}"


## bam metrics
rule mapped_reads_bases:
    """Calculate number of mapped reads"""
    input:
        bam = rules.markdup.output.bam,
        pywc = config["py_wordcount"]
    output:
        reads="{sample}/bams/{sample}.mapped.num",
        bases="{sample}/bams/{sample}.mapped.basenum"
    singularity: containers["samtools-1.7-python-3.6"]
    shell: "samtools view -F 4 {input.bam} | cut -f 10 | python {input.pywc} "
           "--reads {output.reads} --bases {output.bases}"


rule unique_reads_bases:
    """Calculate number of unique reads"""
    input:
        bam = rules.markdup.output.bam,
        pywc=config["py_wordcount"]
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
    """
    input:
        r1=get_forward,
        r2=get_reverse
    output:
        directory("{sample}/pre_process/raw-{sample}-{read_group}/")
    singularity: containers["fastqc"]
    shell: "fastqc --threads 4 --nogroup -o {output} {input.r1} {input.r2} "


rule fastqc_postqc:
    """
    Run fastqc on fastq files post pre-processing
    """
    input:
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2
    output:
        directory("{sample}/pre_process/trimmed-{sample}-{read_group}/")
    singularity: containers["fastqc"]
    shell: "fastqc --threads 4 --nogroup -o {output} {input.r1} {input.r2} "


## coverage
rule covstats:
    """Calculate coverage statistics on bam file"""
    input:
        bam = rules.markdup.output.bam,
        genome="current.genome",
        covpy=config["covstats"],
        bed=config.get("bedfile","")
    params:
        subt="Sample {sample}"
    output:
        covj="{sample}/coverage/covstats.json",
        covp="{sample}/coverage/covstats.png"
    singularity: containers["bedtools-2.26-python-2.7"]
    shell: "bedtools coverage -sorted -g {input.genome} -a {input.bed} "
           "-b {input.bam} -d  | python {input.covpy} - --plot {output.covp} "
           "--title 'Targets coverage' --subtitle '{params.subt}' "
           "> {output.covj}"


rule vtools_coverage:
    """Calculate coverage statistics per transcript"""
    input:
        gvcf = rules.gvcf_gather.output.gvcf,
        tbi = rules.gvcf_gather.output.gvcf_tbi,
        ref = config.get('refflat', "")
    output:
        tsv="{sample}/coverage/refFlat_coverage.tsv"
    singularity: containers["vtools"]
    shell: "vtools-gcoverage -I {input.gvcf} -R {input.ref} > {output.tsv}"


rule collect_cutadapt_summary:
    """ Colect cutadapt summary from each readgroup per sample """
    input:
        cutadapt = lambda wildcards:
        ("{sample}/pre_process/{sample}-{read_group}.txt".format(sample=wildcards.sample,
        read_group=read_group) for read_group in get_readgroup(wildcards)),
        cutadapt_summary= config["cutadapt_summary"]
    output: "{sample}/cutadapt.json"
    singularity: containers["python3"]
    shell: "python {input.cutadapt_summary} --sample {wildcards.sample} "
           "--cutadapt-summary {input.cutadapt} > {output}"


## collection
if "bedfile" in config:
    rule collectstats:
        """Collect all stats for a particular sample with beds"""
        input:
            mnum = rules.mapped_reads_bases.output.reads,
            mbnum = rules.mapped_reads_bases.output.bases,
            unum = rules.unique_reads_bases.output.reads,
            ubnum = rules.unique_reads_bases.output.bases,
            cov = rules.covstats.output.covj,
            cutadapt = rules.collect_cutadapt_summary.output,
            colpy=config["collect_stats"]
        params:
            sample_name="{sample}",
            fthresh=config["female_threshold"]
        output: "{sample}/{sample}.stats.json"
        singularity: containers["vtools"]
        shell: "python {input.colpy} --sample-name {params.sample_name} "
               "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} "
               "--unique-num {input.unum} --usable-basenum {input.ubnum} "
               "--female-threshold {params.fthresh} "
               "--cutadapt {input.cutadapt} "
               "{input.cov} > {output}"
else:
    rule collectstats:
        """Collect all stats for a particular sample without beds"""
        input:
            mnum = rules.mapped_reads_bases.output.reads,
            mbnum = rules.mapped_reads_bases.output.bases,
            unum = rules.unique_reads_bases.output.reads,
            ubnum = rules.unique_reads_bases.output.bases,
            cutadapt = rules.collect_cutadapt_summary.output,
            colpy=config["collect_stats"]
        params:
            sample_name = "{sample}",
            fthresh = config["female_threshold"]
        output: "{sample}/{sample}.stats.json"
        singularity: containers["vtools"]
        shell: "python {input.colpy} --sample-name {params.sample_name} "
               "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} "
               "--unique-num {input.unum} --usable-basenum {input.ubnum} "
               "--female-threshold {params.fthresh} "
               "--cutadapt {input.cutadapt} "
               "> {output}"

rule merge_stats:
    """Merge all stats of all samples"""
    input:
        cols=expand("{sample}/{sample}.stats.json", sample=config['samples']),
        mpy=config["merge_stats"]
    output: "stats.json"
    singularity: containers["vtools"]
    shell: "python {input.mpy} --collectstats {input.cols} "
           "> {output}"


rule stats_tsv:
    """Convert stats.json to tsv"""
    input:
        stats = rules.merge_stats.output,
        sc=config["stats_to_tsv"]
    output:
        stats="stats.tsv"
    singularity: containers["python3"]
    shell: "python {input.sc} -i {input.stats} > {output.stats}"


rule multiple_metrics:
    """Run picard CollectMultipleMetrics"""
    input:
        bam = rules.markdup.output.bam,
        ref = config["reference"],
    params:
        prefix="{sample}/bams/{sample}",
    output:
        alignment="{sample}/bams/{sample}.alignment_summary_metrics",
        insert="{sample}/bams/{sample}.insert_size_metrics"
    singularity: containers["picard"]
    shell: "picard CollectMultipleMetrics "
           "I={input.bam} O={params.prefix} "
           "R={input.ref} "
           "PROGRAM=CollectAlignmentSummaryMetrics "
           "PROGRAM=CollectInsertSizeMetrics "

rule multiqc:
    """
    Create multiQC report
    Depends on stats.tsv to forcefully run at end of pipeline
    """
    input:
        stats="stats.tsv",
        bam=expand("{sample}/bams/{sample}.bam", sample=config["samples"]),
        metric=expand("{sample}/bams/{sample}.metrics", sample=config["samples"]),
        alignment_metrics=expand("{sample}/bams/{sample}.alignment_summary_metrics", sample=config["samples"]),
        insert_metrics=expand("{sample}/bams/{sample}.insert_size_metrics", sample=config["samples"]),
    params:
        odir=".",
        rdir="multiqc_report"
    output: "multiqc_report/multiqc_report.html"
    singularity: containers["multiqc"]
    shell: "multiqc --data-format json -f -o {params.rdir} {params.odir} || touch {output}"
