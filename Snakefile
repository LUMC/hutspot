import json
from os.path import join
from os import mkdir

from pyfaidx import Fasta

OUT_DIR = config.get("OUTPUT_DIR")
REFERENCE = config.get("REFERENCE")
JAVA = config.get("JAVA")
GATK = config.get("GATK")
DBSNP = config.get("DBSNP")
ONETHOUSAND = config.get("ONETHOUSAND")
HAPMAP = config.get("HAPMAP")
QUEUE = config.get("QUEUE")
BED = config.get("BED")
REFFLAT = config.get("REFFLAT")
FEMALE_THRESHOLD = config.get("FEMALE_THRESHOLD", 0.6)

_this_dir = workflow.current_basedir


env_dir = join(_this_dir, "envs")
main_env = join(_this_dir, "environment.yml")

settings_template = join(join(_this_dir, "templates"), "pipeline_settings.md.j2")

with open(config.get("SAMPLE_CONFIG")) as handle:
    SAMPLE_CONFIG = json.load(handle)
SAMPLES = SAMPLE_CONFIG['samples'].keys()


def split_genome(ref, approx_n_chunks=100):
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


def out_path(path):
    return join(OUT_DIR, path)


try:
    mkdir(out_path("tmp"))
except OSError:
    pass


def get_r1(wildcards):
    s = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    r1 = [x['R1'] for _, x in s['libraries'].items()]
    return r1


def get_r2(wildcards):
    s = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    r2 = [x['R2'] for _, x in s['libraries'].items()]
    return r2


def sample_gender(wildcards):
    sam = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    return sam.get("gender", "null")


rule all:
    input:
        combined=out_path("multisample/genotyped.vcf.gz")

rule genome:
    input: REFERENCE
    output: out_path("current.genome")
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule merge_r1:
    input: get_r1
    output: temp(out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"))
    shell: "cat {input} > {output}"

rule merge_r2:
    input: get_r2
    output: temp(out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz"))
    shell: "cat {input} > {output}"

rule sickle:
    input:
        r1 = out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        r2 = out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    output:
        r1 = temp(out_path("{sample}/pre_process/{sample}.trimmed_R1.fastq")),
        r2 = temp(out_path("{sample}/pre_process/{sample}.trimmed_R2.fastq")),
        s = out_path("{sample}/pre_process/{sample}.trimmed_singles.fastq"),
    conda: "envs/sickle.yml"
    shell: "sickle pe -f {input.r1} -r {input.r2} -t sanger -o {output.r1} " \
           "-p {output.r2} -s {output.s}"

rule cutadapt:
    input:
        r1 = out_path("{sample}/pre_process/{sample}.trimmed_R1.fastq"),
        r2 = out_path("{sample}/pre_process/{sample}.trimmed_R2.fastq")
    output:
        r1 = temp(out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq")),
        r2 = temp(out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq"))
    conda: "envs/cutadapt.yml"
    shell: "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 1 -o {output.r1} " \
           "{input.r1} -p {output.r2} {input.r2}"

rule align:
    input:
        r1 = out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2 = out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq"),
        ref = REFERENCE
    params:
        rg = "@RG\\tID:{sample}_lib1\\tSM:{sample}\\tPL:ILLUMINA"
    output: temp(out_path("{sample}/bams/{sample}.sorted.bam"))
    conda: "envs/bwa.yml"
    shell: "bwa mem -t 8 -R '{params.rg}' {input.ref} {input.r1} {input.r2} " \
           "| picard SortSam CREATE_INDEX=TRUE TMP_DIR=null " \
           "INPUT=/dev/stdin OUTPUT={output} SORT_ORDER=coordinate"

rule markdup:
    input:
        bam = out_path("{sample}/bams/{sample}.sorted.bam"),
    params:
        tmp = out_path("tmp")
    output:
        bam = temp(out_path("{sample}/bams/{sample}.markdup.bam")),
        metrics = out_path("{sample}/bams/{sample}.markdup.metrics")
    conda: "envs/picard.yml"
    shell: "picard MarkDuplicates CREATE_INDEX=TRUE TMP_DIR={params.tmp} " \
           "INPUT={input.bam} OUTPUT={output.bam} " \
           "METRICS_FILE={output.metrics} " \
           "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

rule baserecal:
    input:
        bam = out_path("{sample}/bams/{sample}.markdup.bam"),
        java = JAVA,
        gatk = GATK,
        ref = REFERENCE,
        dbsnp = DBSNP,
        one1kg = ONETHOUSAND,
        hapmap = HAPMAP
    output:
        grp = out_path("{sample}/bams/{sample}.baserecal.grp")
    conda: "envs/gatk.yml"
    shell: "{input.java} -jar {input.gatk} -T BaseRecalibrator " \
           "-I {input.bam} -o {output.grp} -nct 8 -R {input.ref} " \
           "-cov ReadGroupCovariate -cov QualityScoreCovariate " \
           "-cov CycleCovariate -cov ContextCovariate -knownSites " \
           "{input.dbsnp} -knownSites {input.one1kg} " \
           "-knownSites {input.hapmap}"

rule printreads:
    input:
        grp=out_path("{sample}/bams/{sample}.baserecal.grp"),
        bam=out_path("{sample}/bams/{sample}.markdup.bam"),
        java=JAVA,
        gatk=GATK,
        ref=REFERENCE
    output:
        bam=out_path("{sample}/bams/{sample}.baserecal.bam"),
        bai=out_path("{sample}/bams/{sample}.baserecal.bai")
    conda: "envs/gatk.yml"
    shell: "{input.java} -jar {input.gatk} -T PrintReads -I {input.bam} "\
           "-o {output.bam} -R {input.ref} -BQSR {input.grp}"


rule gvcf_scatter:
    input:
        bam=out_path("{sample}/bams/{sample}.baserecal.bam"),
        dbsnp=DBSNP,
        ref=REFERENCE,
        gatk=GATK
    params:
        chunk="{chunk}"
    output:
        gvcf=out_path("{sample}/vcf/{sample}.{chunk}.part.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -jar {input.gatk} -T HaplotypeCaller -ERC GVCF -I "\
           "{input.bam} -R {input.ref} -D {input.dbsnp} "\
           "-L {params.chunk} -o {output.gvcf}"


rule gvcf_gather:
    input:
        gvcfs=expand(out_path("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz"),
                     chunk=CHUNKS),
        ref=REFERENCE,
        gatk=GATK
    params:
        gvcfs=" -V ".join(expand(out_path("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz"),
                                 chunk=CHUNKS))
    output:
        gvcf=out_path("{sample}/vcf/{sample}.g.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "\
           "-R {input.ref} -V {params.gvcfs} -output {output.gvcf} "\
           "-assumeSorted"


rule genotype_scatter:
    input:
        gvcfs = expand(out_path("{sample}/vcf/{sample}.g.vcf.gz"),
                       sample=SAMPLES),
        ref=REFERENCE,
        gatk=GATK
    params:
        li=" -V ".join(expand(out_path("{sample}/vcf/{sample}.g.vcf.gz"),
                              sample=SAMPLES)),
        chunk="{chunk}"
    output:
        vcf=out_path("multisample/genotype.{chunk}.part.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -jar {input.gatk} -T GenotypeGVCFs -R {input.ref} "\
           "-V {params.li} -L {params.chunk} -o {output.vcf}"


rule genotype_gather:
    input:
        vcfs=expand(out_path("multisample/genotype.{chunk}.part.vcf.gz"),
                    chunk=CHUNKS),
        ref=REFERENCE,
        gatk=GATK
    params:
        vcfs=" -V ".join(expand(out_path("multisample/genotype.{chunk}.part.vcf.gz"),
                                chunk=CHUNKS))
    output:
        combined=out_path("multisample/genotyped.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -cp {input.gatk} org.broadinstitute.gatk.tool.CatVariants "\
           "-R {input.ref} -V {params.vcfs} -output {output.combined} "\
           "-assumeSorted"


