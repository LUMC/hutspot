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
set_default("safe_fastqc", "src/safe_fastqc.sh")
set_default("py_wordcount", "src/pywc.py")

# Generate the input string for basrecalibration
BSQR_known_sites = ''
for argument, filename in zip(itertools.repeat('-knownSites'), settings["known_sites"]):
    BSQR_known_sites +=' {} {}'.format(argument, filename)

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
    "sickle": "docker://quay.io/biocontainers/sickle-trim:1.33--ha92aebf_4",
    "tabix": "docker://quay.io/biocontainers/tabix:0.2.6--ha92aebf_0",
    "vtools": "docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0"
}

def get_r(strand, wildcards):
    """Get fastq files on a single strand for a sample"""
    s = settings['samples'].get(wildcards.sample)
    rs = []
    for l in sorted(s['libraries'].keys()):
        rs.append(s['libraries'][l][strand])
    return rs

get_r1 = partial(get_r, "R1")
get_r2 = partial(get_r, "R2")


def get_bedpath(wildcards):
    """Get absolute path of a bed file"""
    return next(x for x in BEDS if basename(x) == wildcards.bed)


def sample_gender(wildcards):
    """Get sample gender"""
    sam = settings['samples'].get(wildcards.sample)
    return sam.get("gender", "null")


def metrics(do_metrics=True):
    if not do_metrics:
        return ""

    fqcr = expand("{sample}/pre_process/raw_fastqc/.done.txt",
                  sample=settings['samples']),
    fqcm = expand("{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip",
                  sample=settings['samples']),
    fqcp = expand("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip",
                  sample=settings['samples']),
    if "refflat" in settings:
        coverage_stats = tuple(expand("{sample}/refFlat_coverage.tsv", sample=settings['samples']))
    else:
        coverage_stats = tuple()
    stats = "stats.json",
    print(coverage_stats)
    return  fqcr + fqcm + fqcp + coverage_stats + stats


rule all:
    input:
        multiqc="multiqc_report/multiqc_report.html",
        bais=expand("{sample}/bams/{sample}.markdup.bam.bai", sample=settings['samples']),
        vcfs=expand("{sample}/vcf/{sample}.vcf.gz", sample=settings['samples']),
        vcf_tbi=expand("{sample}/vcf/{sample}.vcf.gz.tbi", sample=settings['samples']),
        gvcfs=expand("{sample}/vcf/{sample}.g.vcf.gz", sample=settings['samples']),
        gvcf_tbi=expand("{sample}/vcf/{sample}.g.vcf.gz.tbi", sample=settings['samples']),
        stats=metrics()


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

# contains original merged fastq files as input to prevent them from being prematurely deleted
rule sickle:
    """Trim fastq files"""
    input:
        r1 = "{sample}/pre_process/{sample}.merged_R1.fastq.gz",
        r2 = "{sample}/pre_process/{sample}.merged_R2.fastq.gz"
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
        ref = settings["reference"],
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
        ref = settings["reference"],
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
           "-BQSR {input.bqsr}"

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
        bam="{sample}/bams/{sample}.sorted.bam",
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
        fq=settings["safe_fastqc"]
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
        fq=settings["safe_fastqc"]
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
        sc=settings["fastq_stats"]
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
        covpy=settings["covstats"],
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
        ref = settings.get('refflat', "")
    output:
        tsv="{sample}/refFlat_coverage.tsv"
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
if "bedfile" in settings:
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
            cov=expand("{{sample}}/coverage/{bed}.covstats.json", bed=settings["bedfile"]),
            colpy=settings["collect_stats"]
        params:
            sample_name="{sample}",
            fthresh=settings["female_threshold"]
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
            colpy = settings["collect_stats"]
        params:
            sample_name = "{sample}",
            fthresh = settings["female_threshold"]
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
        cols=expand("{sample}/{sample}.stats.json", sample=settings['samples']),
        vstat=expand("{sample}/vcf/{sample}.vcfstats.json", sample=settings['samples']),
        mpy=settings["merge_stats"]
    output:
        stats="stats.json"
    singularity: containers["vtools"]
    shell: "python {input.mpy} --vcfstats {input.vstat} {input.cols} "
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
        stats="stats.tsv"
    params:
        odir=".",
        rdir="multiqc_report"
    output:
        "multiqc_report/multiqc_report.html"
    singularity: containers["multiqc"]
    shell: "multiqc -f -o {params.rdir} {params.odir} || touch {output}"
