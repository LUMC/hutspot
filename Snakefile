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

# Read the config json
with open(config['CONFIG_JSON'], 'rt') as fin:
    settings = json.load(fin)

# Read the json schema
with open('config/schema.json', 'rt') as fin:
    schema = json.load(fin)

# Validate the settings against the schema
try:
    jsonschema.validate(settings, schema)
except jsonschema.ValidationError as e:
    raise jsonschema.ValidationError(f'Invalid CONFIG_JSON: {e}')


# Set default values
def set_default(key, value):
    """ Set default config values """
    if key not in settings:
        settings[key] = value

# Set the default settings
set_default('scatter_size', 1000000000)
set_default('female_threshold', 0.6)

# Set the script paths
set_default("covstats", "src/covstats.py")
set_default("collect_stats", "src/collect_stats.py")
set_default("merge_stats", "src/merge_stats.py")
set_default("fastq_stats", "src/fastqc_stats.py")
set_default("stats_to_tsv", "src/stats_to_tsv.py")
set_default("py_wordcount", "src/pywc.py")
set_default("collect_metrics", "src/collect_metrics.py")
set_default("merge_metrics", "src/merge_metrics.py")

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
    "picard-2.14": "docker://quay.io/biocontainers/picard:2.14--py36_0",
    "python3": "docker://python:3.6-slim",
    "samtools-1.7-python-3.6": "docker://quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:1abf1824431ec057c7d41be6f0c40e24843acde4-0",
    "sickle": "docker://quay.io/biocontainers/sickle-trim:1.33--ha92aebf_4",
    "tabix": "docker://quay.io/biocontainers/tabix:0.2.6--ha92aebf_0",
    "vtools": "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
}

def get_forward(wildcards):
    """ Get the forward fastq file from the config """
    return (
        settings["samples"][wildcards.sample]["libraries"]
                [wildcards.read_group]["R1"]
    )

def get_reverse(wildcards):
    """ Get the reverse fastq file from the config """
    return (
        settings["samples"][wildcards.sample]["libraries"]
            [wildcards.read_group]["R2"]
    )

def get_readgroup(wildcards):
    return settings["samples"][wildcards.sample]["libraries"]

def get_readgroup_per_sample():
    for sample in settings["samples"]:
        for rg in settings["samples"][sample]["libraries"]:
            yield rg, sample


def coverage_stats(wildcards):
    files = expand("{sample}/coverage/refFlat_coverage.tsv", sample=settings['samples'])
    return files if "refflat" in settings else []

rule all:
    input:
        multiqc="multiqc_report/multiqc_report.html",
        stats = "stats.json",
        stats_tsv = "stats.tsv",
        bais=expand("{sample}/bams/{sample}.markdup.bam.bai", sample=settings['samples']),
        vcfs=expand("{sample}/vcf/{sample}.vcf.gz", sample=settings['samples']),
        vcf_tbi=expand("{sample}/vcf/{sample}.vcf.gz.tbi", sample=settings['samples']),
        gvcfs=expand("{sample}/vcf/{sample}.g.vcf.gz", sample=settings['samples']),
        gvcf_tbi=expand("{sample}/vcf/{sample}.g.vcf.gz.tbi", sample=settings['samples']),
        fastqc_raw = (f"{sample}/pre_process/raw-{sample}-{read_group}/" for read_group, sample in get_readgroup_per_sample()),
        fastqc_trimmed = (f"{sample}/pre_process/trimmed-{sample}-{read_group}/" for read_group, sample in get_readgroup_per_sample()),
        cutadapt = (f"{sample}/pre_process/{sample}-{read_group}.txt" for read_group, sample in get_readgroup_per_sample()),
        sample_metrics = expand("{sample}/metrics.json", sample=settings['samples']),
        metrics_json = "metrics.json",
        metrics_tsv = "metrics.tsv",
        coverage_stats = coverage_stats,
        #covstat_png=expand("{sample}/coverage/covstats.png", sample=settings['samples'])


rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output: directory("tmp")
    singularity: containers["debian"]
    shell: "mkdir -p {output}"

rule genome:
    """Create genome file as used by bedtools"""
    input: settings["reference"]
    output: "current.genome"
    singularity: containers["debian"]
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule cutadapt:
    """Clip fastq files"""
    input:
        r1=get_forward,
        r2=get_reverse
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
        r1 = "{sample}/pre_process/{sample}-{read_group}_R1.fastq.gz",
        r2 = "{sample}/pre_process/{sample}-{read_group}_R2.fastq.gz",
        ref = settings["reference"],
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
        bam = "{sample}/bams/{sample}.markdup.bam",
        bai = "{sample}/bams/{sample}.markdup.bai",
        metrics = "{sample}/bams/{sample}.markdup.metrics"
    params:
        bams=markdup_bam_input
    singularity: containers["picard-2.14"]
    shell: "picard -Xmx4G -Djava.io.tmpdir={input.tmp} MarkDuplicates "
           "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "{params.bams} OUTPUT={output.bam} "
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
        ref = settings["reference"],
        vcfs = settings["known_sites"]
    output:
        grp = "{sample}/bams/{sample}.baserecal.grp"
    params: " ".join(expand("-knownSites {vcf}", vcf=settings["known_sites"]))
    singularity: containers["gatk"]
    shell: "java -XX:ParallelGCThreads=1 -jar /usr/GenomeAnalysisTK.jar -T "
           "BaseRecalibrator -I {input.bam} -o {output.grp} -nct 8 "
           "-R {input.ref} -cov ReadGroupCovariate -cov QualityScoreCovariate "
           "-cov CycleCovariate -cov ContextCovariate {params}"

checkpoint scatterregions:
    """Scatter the reference genome"""
    input:
        ref = settings["reference"],
    params:
        size = settings['scatter_size']
    output:
        directory("scatter")
    singularity: containers["biopet-scatterregions"]
    shell: "mkdir -p {output} && "
           "biopet-scatterregions "
           "--referenceFasta {input.ref} --scatterSize {params.size} "
           "--outputDir scatter"

rule gvcf_scatter:
    """Run HaplotypeCaller in GVCF mode by chunk"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        bqsr="{sample}/bams/{sample}.baserecal.grp",
        dbsnp=settings["dbsnp"],
        ref=settings["reference"],
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
        gvcf = "{sample}/vcf/{sample}.{chunk}.g.vcf.gz",
        tbi = "{sample}/vcf/{sample}.{chunk}.g.vcf.gz.tbi",
        ref= settings["reference"]
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
        bam="{sample}/bams/{sample}.markdup.bam",
        pywc=settings["py_wordcount"]
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
        pywc=settings["py_wordcount"]
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
        r1="{sample}/pre_process/{sample}-{read_group}_R1.fastq.gz",
        r2="{sample}/pre_process/{sample}-{read_group}_R2.fastq.gz",
    output:
        directory("{sample}/pre_process/trimmed-{sample}-{read_group}/")
    singularity: containers["fastqc"]
    shell: "fastqc --threads 4 --nogroup -o {output} {input.r1} {input.r2} "


## coverages

rule covstats:
    """Calculate coverage statistics on bam file"""
    input:
        bam="{sample}/bams/{sample}.markdup.bam",
        genome="current.genome",
        covpy=settings["covstats"],
        bed=settings.get("bedfile","")
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
        gvcf="{sample}/vcf/{sample}.g.vcf.gz",
        tbi = "{sample}/vcf/{sample}.g.vcf.gz.tbi",
        ref = settings.get('refflat', "")
    output:
        tsv="{sample}/coverage/refFlat_coverage.tsv"
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


# collection
if "bedfile" in settings:
    rule collectstats:
        """Collect all stats for a particular sample with beds"""
        input:
            mnum="{sample}/bams/{sample}.mapped.num",
            mbnum="{sample}/bams/{sample}.mapped.basenum",
            unum="{sample}/bams/{sample}.unique.num",
            ubnum="{sample}/bams/{sample}.usable.basenum",
            cov="{sample}/coverage/covstats.json",
            cutadapt = "{sample}/metrics.json",
            colpy=settings["collect_stats"]
        params:
            sample_name="{sample}",
            fthresh=settings["female_threshold"]
        output:
            "{sample}/{sample}.stats.json"
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
            mnum = "{sample}/bams/{sample}.mapped.num",
            mbnum = "{sample}/bams/{sample}.mapped.basenum",
            unum = "{sample}/bams/{sample}.unique.num",
            ubnum = "{sample}/bams/{sample}.usable.basenum",
            cutadapt = "{sample}/metrics.json",
            colpy = settings["collect_stats"]
        params:
            sample_name = "{sample}",
            fthresh = settings["female_threshold"]
        output:
            "{sample}/{sample}.stats.json"
        singularity: containers["vtools"]
        shell: "python {input.colpy} --sample-name {params.sample_name} "
               "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} "
               "--unique-num {input.unum} --usable-basenum {input.ubnum} "
               "--female-threshold {params.fthresh} "
               "--cutadapt {input.cutadapt} "
               "> {output}"

rule collect_metrics:
    """ Colect metrics for each sample """
    input:
        cutadapt = lambda wildcards:
        ("{sample}/pre_process/{sample}-{read_group}.txt".format(sample=wildcards.sample,
        read_group=read_group) for read_group in get_readgroup(wildcards)),
        collect_metrics = settings["collect_metrics"]
    output: "{sample}/metrics.json"
    singularity: containers["python3"]
    shell: "python {input.collect_metrics} --sample {wildcards.sample} "
           "--cutadapt-summary {input.cutadapt} > {output}"

rule merge_metrics:
   """ Merge per sample metrics files """
   input:
        metrics = expand("{sample}/metrics.json", sample=settings["samples"]),
        merge_metrics = settings["merge_metrics"]
   output:
        json = "metrics.json",
        tsv = "metrics.tsv"
   singularity: containers["python3"]
   shell: "python {input.merge_metrics} --metrics {input.metrics} "
          "--json {output.json} --tsv {output.tsv}"

rule merge_stats:
    """Merge all stats of all samples"""
    input:
        cols=expand("{sample}/{sample}.stats.json", sample=settings['samples']),
        vstat=expand("{sample}/vcf/{sample}.vcfstats.json", sample=settings['samples']),
        mpy=settings["merge_stats"]
    output:
        stats="stats.json"
    singularity: containers["vtools"]
    shell: "python {input.mpy} --vcfstats {input.vstat} --collectstats {input.cols} "
           "> {output.stats}"


rule stats_tsv:
    """Convert stats.json to tsv"""
    input:
        stats="stats.json",
        sc=settings["stats_to_tsv"]
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
        stats=expand("{sample}/metrics.json", sample=settings["samples"])
    params:
        odir=".",
        rdir="multiqc_report"
    output:
        "multiqc_report/multiqc_report.html"
    singularity: containers["multiqc"]
    shell: "multiqc -f -o {params.rdir} {params.odir} --ignore metrics.tsv || touch {output}"
