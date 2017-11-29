import json
from os.path import join, basename
from os import mkdir

from pyfaidx import Fasta

OUT_DIR = config.get("OUTPUT_DIR")
REFERENCE = config.get("REFERENCE")
JAVA = config.get("JAVA")
GATK = config.get("GATK")
DBSNP = config.get("DBSNP")
ONETHOUSAND = config.get("ONETHOUSAND")
HAPMAP = config.get("HAPMAP")
BED = config.get("BED", "")  # comma-separated list of BED files
REFFLAT = config.get("REFFLAT", "")  # comma-separated list of refFlat files
FEMALE_THRESHOLD = config.get("FEMALE_THRESHOLD", 0.6)
FASTQ_COUNT = config.get("FASTQ_COUNT")
MAX_BASES = config.get("MAX_BASES", "")

_this_dir = workflow.current_basedir


env_dir = join(_this_dir, "envs")
main_env = join(_this_dir, "environment.yml")

settings_template = join(join(_this_dir, "templates"), "pipeline_settings.md.j2")
covpy = join(join(_this_dir, "src"), "covstats.py")
colpy = join(join(_this_dir, "src"), "collect_stats.py")
vs_py = join(join(_this_dir, "src"), "vcfstats.py")
mpy = join(join(_this_dir, "src"), "merge_stats.py")

with open(config.get("SAMPLE_CONFIG")) as handle:
    SAMPLE_CONFIG = json.load(handle)
SAMPLES = SAMPLE_CONFIG['samples'].keys()

BEDS = BED.split(",")
REFFLATS = REFFLAT.split(",")

BASE_BEDS = [basename(x) for x in BEDS]
BASE_REFFLATS = [basename(x) for x in BEDS]


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
    r1 = []
    for l in sorted(s['libraries'].keys()):
        r1.append(s['libraries'][l]['R1'])
    return r1


def get_r2(wildcards):
    s = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    r2 = []
    for l in sorted(s['libraries'].keys()):
        r2.append(s['libraries'][l]['R2'])
    return r2


def get_bedpath(wildcards):
    return [x for x in BEDS if basename(x) == wildcards.bed][0]


def get_refflatpath(wildcards):
    return [x for x in REFFLATS if basename(x) == wildcards.refflat][0]


def sample_gender(wildcards):
    sam = SAMPLE_CONFIG['samples'].get(wildcards.sample)
    return sam.get("gender", "null")


def metrics(do_metrics=True):
    if not do_metrics:
        return ""

    fqcr = expand(out_path("{sample}/pre_process/raw_fastqc/.done.txt"),
                  sample=SAMPLES)
    fqcm = expand(out_path("{sample}/pre_process/merged_fastqc/.done.txt"),
                  sample=SAMPLES)
    fqcp = expand(out_path("{sample}/pre_process/postqc_fastqc/.done.txt"),
                  sample=SAMPLES)
    stats = out_path("stats.json")
    return  fqcr + fqcm + fqcp + [stats]


rule all:
    input:
        combined=out_path("multisample/genotyped.vcf.gz"),
        stats=metrics()

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


rule seqtk_r1:
    input:
        stats=out_path("{sample}/pre_process/{sample}.preqc_count.json"),
        fastq=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz")
    params:
        max_bases=MAX_BASES
    output:
        fastq=temp(out_path("{sample}/pre_process/{sample}.sampled_R1.fastq.gz"))
    conda: "envs/seqtk.yml"
    script: "src/seqtk.py"


rule seqtk_r2:
    input:
        stats = out_path("{sample}/pre_process/{sample}.preqc_count.json"),
        fastq = out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    params:
        max_bases = MAX_BASES
    output:
        fastq = temp(out_path("{sample}/pre_process/{sample}.sampled_R2.fastq.gz"))
    conda: "envs/seqtk.yml"
    script: "src/seqtk.py"


# contains original merged fastq files as input to prevent them from being prematurely deleted
rule sickle:
    input:
        r1 = out_path("{sample}/pre_process/{sample}.sampled_R1.fastq.gz"),
        r2 = out_path("{sample}/pre_process/{sample}.sampled_R2.fastq.gz"),
        rr1 = out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        rr2 = out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
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
    shell: "java -jar -Xmx4G {input.gatk} -T HaplotypeCaller -ERC GVCF -I "\
           "{input.bam} -R {input.ref} -D {input.dbsnp} "\
           "-L {params.chunk} -o {output.gvcf} "\
           "-variant_index_type LINEAR -variant_index_parameter 128000"


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
    shell: "java -Xmx4G -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "\
           "-R {input.ref} -V {params.gvcfs} -out {output.gvcf} "\
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
    shell: "java -jar -Xmx4G {input.gatk} -T GenotypeGVCFs -R {input.ref} "\
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
    shell: "java -Xmx4G -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "\
           "-R {input.ref} -V {params.vcfs} -out {output.combined} "\
           "-assumeSorted"


## bam metrics

rule mapped_num:
    input:
        bam=out_path("{sample}/bams/{sample}.sorted.bam")
    output:
        num=out_path("{sample}/bams/{sample}.mapped.num")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | wc -l > {output.num}"


rule mapped_basenum:
    input:
        bam=out_path("{sample}/bams/{sample}.sorted.bam")
    output:
        num=out_path("{sample}/bams/{sample}.mapped.basenum")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | wc -c > {output.num}"


rule unique_num:
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam")
    output:
        num=out_path("{sample}/bams/{sample}.unique.num")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | wc -l > {output.num}"


rule usable_basenum:
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam")
    output:
        num=out_path("{sample}/bams/{sample}.usable.basenum")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | cut -f10 | wc -c > {output.num}"


## fastqc

rule fastqc_raw:
    input:
        r1=get_r1,
        r2=get_r2
    params:
        odir=out_path("{sample}/pre_process/raw_fastqc")
    output:
        aux=out_path("{sample}/pre_process/raw_fastqc/.done.txt")
    conda: "envs/fastqc.yml"
    shell: "fastqc -o {params.odir} {input.r1} {input.r2} && echo 'done' > {output.aux}"


rule fastqc_merged:
    input:
        r1=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        r2=out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    params:
        odir=out_path("{sample}/pre_process/merged_fastqc")
    output:
        aux=out_path("{sample}/pre_process/merged_fastqc/.done.txt")
    conda: "envs/fastqc.yml"
    shell: "fastqc -o {params.odir} {input.r1} {input.r2} && echo 'done' > {output.aux}"


rule fastqc_postqc:
    input:
        r1=out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2=out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq")
    params:
        odir=out_path("{sample}/pre_process/postqc_fastqc")
    output:
        aux=out_path("{sample}/pre_process/postqc_fastqc/.done.txt")
    conda: "envs/fastqc.yml"
    shell: "fastqc -o {params.odir} {input.r1} {input.r2} && echo 'done' > {output.aux}"


## fastq-count

rule fqcount_preqc:
    input:
        r1=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        r2=out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    params:
        fastqcount=FASTQ_COUNT
    output:
        out_path("{sample}/pre_process/{sample}.preqc_count.json")
    shell: "{params.fastqcount} {input.r1} {input.r2} > {output}"


rule fqcount_postqc:
    input:
        r1=out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2=out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq")
    params:
        fastqcount=FASTQ_COUNT
    output:
        out_path("{sample}/pre_process/{sample}.postqc_count.json")
    shell: "{params.fastqcount} {input.r1} {input.r2} > {output}"


## coverages

rule covstats:
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam"),
        genome=out_path("current.genome"),
        covpy=covpy,
        bed=get_bedpath
    params:
        subt="Sample {sample}"
    output:
        covj=out_path("{sample}/coverage/{bed}.covstats.json"),
        covp=out_path("{sample}/coverage/{bed}.covstats.png")
    conda: "envs/covstat.yml"
    shell: "bedtools coverage -sorted -g {input.genome} -a {input.bed} -b {input.bam} " \
           "-d  | python {input.covpy} - --plot {output.covp} " \
           "--title 'Targets coverage' --subtitle '{params.subt}' > {output.covj}"


## vcfstats

rule vcfstats:
    input:
        vcf=out_path("multisample/genotyped.vcf.gz"),
        vs_py=vs_py
    output:
        stats=out_path("multisample/vcfstats.json")
    conda: "envs/vcfstats.yml"
    shell: "python {input.vs_py} -i {input.vcf} > {output.stats}"


## collection
rule collectstats:
    input:
        preqc=out_path("{sample}/pre_process/{sample}.preqc_count.json"),
        postq=out_path("{sample}/pre_process/{sample}.postqc_count.json"),
        mnum=out_path("{sample}/bams/{sample}.mapped.num"),
        mbnum=out_path("{sample}/bams/{sample}.mapped.basenum"),
        unum=out_path("{sample}/bams/{sample}.unique.num"),
        ubnum=out_path("{sample}/bams/{sample}.usable.basenum"),
        cov=expand(out_path("{{sample}}/coverage/{bed}.covstats.json"), bed=BASE_BEDS),
        colpy=colpy
    params:
        sample_name="{sample}",
        fthresh=FEMALE_THRESHOLD
    output:
        out_path("{sample}/{sample}.stats.json")
    conda: "envs/collectstats.yml"
    shell: "python {input.colpy} --sample-name {params.sample_name} " \
           "--pre-qc-fastq {input.preqc} --post-qc-fastq {input.postq} " \
           "--mapped-num {input.mnum} --mapped-basenum {input.mbnum} " \
           "--unique-num {input.unum} --usable-basenum {input.ubnum} " \
           "--female-threshold {params.fthresh} {input.cov} > {output}"

rule merge_stats:
    input:
        cols=expand(out_path("{sample}/{sample}.stats.json"), sample=SAMPLES),
        vstat=out_path("multisample/vcfstats.json"),
        mpy=mpy
    output:
        stats=out_path("stats.json")
    conda: "envs/collectstats.yml"
    shell: "python {input.mpy} --vcfstats {input.vstat} {input.cols} > {output.stats}"
