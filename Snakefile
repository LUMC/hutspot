import json
from functools import partial
from os.path import join, basename

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

if FASTQ_COUNT is None:
    fqc = "python {0}".format(fsrc_dir("src", "fastq-count.py"))
else:
    fqc = FASTQ_COUNT

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


def out_path(path):
    return join(OUT_DIR, path)

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

    fqcr = expand(out_path("{sample}/pre_process/raw_fastqc/.done.txt"),
                  sample=SAMPLES)
    fqcm = expand(out_path("{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip"),
                  sample=SAMPLES)
    fqcp = expand(out_path("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip"),
                  sample=SAMPLES)
    if len(REFFLATS) >= 1:
        coverage_stats = expand(out_path("{sample}/coverage/{ref}.coverages.tsv"),
                                sample=SAMPLES, ref=BASE_REFFLATS)
    else:
        coverage_stats = []
    stats = out_path("stats.json")
    return  fqcr + fqcm + fqcp + coverage_stats + [stats]


rule all:
    input:
        combined=out_path("multisample/genotyped.vcf.gz"),
        bais=expand(out_path("{sample}/bams/{sample}.markdup.bam.bai"),
                    sample=SAMPLES),
        stats=metrics()


rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output: ancient(out_path("tmp"))
    shell: "mkdir -p {output}"

rule genome:
    """Create genome file as used by bedtools"""
    input: REFERENCE
    output: out_path("current.genome")
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output}"

rule merge_r1:
    """Merge all forward fastq files into one"""
    input: get_r1
    output: temp(out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"))
    shell: "cat {input} > {output}"

rule merge_r2:
    """Merge all reverse fastq files into one"""
    input: get_r2
    output: temp(out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz"))
    shell: "cat {input} > {output}"


rule seqtk_r1:
    """Either subsample or link forward fastq file"""
    input:
        stats=out_path("{sample}/pre_process/{sample}.preqc_count.json"),
        fastq=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        seqtk=seq
    params:
        max_bases=str(MAX_BASES)
    output:
        fastq=temp(out_path("{sample}/pre_process/{sample}.sampled_R1.fastq.gz"))
    conda: "envs/seqtk.yml"
    shell: "bash {input.seqtk} {input.stats} {input.fastq} {output.fastq} {params.max_bases}"


rule seqtk_r2:
    """Either subsample or link reverse fastq file"""
    input:
        stats = out_path("{sample}/pre_process/{sample}.preqc_count.json"),
        fastq = out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz"),
        seqtk=seq
    params:
        max_bases =str(MAX_BASES)
    output:
        fastq = temp(out_path("{sample}/pre_process/{sample}.sampled_R2.fastq.gz"))
    conda: "envs/seqtk.yml"
    shell: "bash {input.seqtk} {input.stats} {input.fastq} {output.fastq} {params.max_bases}"


# contains original merged fastq files as input to prevent them from being prematurely deleted
rule sickle:
    """Trim fastq files"""
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
    shell: "sickle pe -f {input.r1} -r {input.r2} -t sanger -o {output.r1} "
           "-p {output.r2} -s {output.s}"

rule cutadapt:
    """Clip fastq files"""
    input:
        r1 = out_path("{sample}/pre_process/{sample}.trimmed_R1.fastq"),
        r2 = out_path("{sample}/pre_process/{sample}.trimmed_R2.fastq")
    output:
        r1 = temp(out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq")),
        r2 = temp(out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq"))
    conda: "envs/cutadapt.yml"
    shell: "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 1 -o {output.r1} "
           "{input.r1} -p {output.r2} {input.r2}"

rule align:
    """Align fastq files"""
    input:
        r1 = out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2 = out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq"),
        ref = REFERENCE
    params:
        rg = "@RG\\tID:{sample}_lib1\\tSM:{sample}\\tPL:ILLUMINA"
    output: temp(out_path("{sample}/bams/{sample}.sorted.bam"))
    conda: "envs/bwa.yml"
    shell: "bwa mem -t 8 -R '{params.rg}' {input.ref} {input.r1} {input.r2} "
           "| picard SortSam CREATE_INDEX=TRUE TMP_DIR=null "
           "INPUT=/dev/stdin OUTPUT={output} SORT_ORDER=coordinate"

rule markdup:
    """Mark duplicates in BAM file"""
    input:
        bam = out_path("{sample}/bams/{sample}.sorted.bam"),
        tmp = ancient(out_path("tmp"))
    output:
        bam = out_path("{sample}/bams/{sample}.markdup.bam"),
        bai = out_path("{sample}/bams/{sample}.markdup.bai"),
        metrics = out_path("{sample}/bams/{sample}.markdup.metrics")
    conda: "envs/picard.yml"
    shell: "picard MarkDuplicates CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
           "INPUT={input.bam} OUTPUT={output.bam} "
           "METRICS_FILE={output.metrics} "
           "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

rule bai:
    """Copy bai files as some genome browsers can only .bam.bai files"""
    input:
        bai = out_path("{sample}/bams/{sample}.markdup.bai")
    output:
        bai = out_path("{sample}/bams/{sample}.markdup.bam.bai")
    shell: "cp {input.bai} {output.bai}"

rule baserecal:
    """Base recalibrated BAM files"""
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
    shell: "{input.java} -jar {input.gatk} -T BaseRecalibrator "
           "-I {input.bam} -o {output.grp} -nct 8 -R {input.ref} "
           "-cov ReadGroupCovariate -cov QualityScoreCovariate "
           "-cov CycleCovariate -cov ContextCovariate -knownSites "
           "{input.dbsnp} -knownSites {input.one1kg} "
           "-knownSites {input.hapmap}"

rule gvcf_scatter:
    """Run HaplotypeCaller in GVCF mode by chunk"""
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam"),
        bqsr=out_path("{sample}/bams/{sample}.baserecal.grp"),
        dbsnp=DBSNP,
        ref=REFERENCE,
        gatk=GATK
    params:
        chunk="{chunk}"
    output:
        gvcf=out_path("{sample}/vcf/{sample}.{chunk}.part.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -jar -Xmx4G {input.gatk} -T HaplotypeCaller -ERC GVCF -I "
           "{input.bam} -R {input.ref} -D {input.dbsnp} "
           "-L '{params.chunk}' -o '{output.gvcf}' "
           "-variant_index_type LINEAR -variant_index_parameter 128000 "
           "-BQSR {input.bqsr}"


rule gvcf_gather:
    """Gather all gvcf scatters"""
    input:
        gvcfs=expand(out_path("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz"),
                     chunk=CHUNKS),
        ref=REFERENCE,
        gatk=GATK
    params:
        gvcfs="' -V '".join(expand(out_path("{{sample}}/vcf/{{sample}}.{chunk}.part.vcf.gz"),
                                 chunk=CHUNKS))
    output:
        gvcf=out_path("{sample}/vcf/{sample}.g.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -Xmx4G -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "
           "-R {input.ref} -V '{params.gvcfs}' -out {output.gvcf} "
           "-assumeSorted"


rule genotype_scatter:
    """Run GATK's GenotypeGVCFs by chunk"""
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
    shell: "java -jar -Xmx4G {input.gatk} -T GenotypeGVCFs -R {input.ref} "
           "-V {params.li} -L '{params.chunk}' -o '{output.vcf}'"


rule genotype_gather:
    """Gather all genotyping scatters"""
    input:
        vcfs=expand(out_path("multisample/genotype.{chunk}.part.vcf.gz"),
                    chunk=CHUNKS),
        ref=REFERENCE,
        gatk=GATK
    params:
        vcfs="' -V '".join(expand(out_path("multisample/genotype.{chunk}.part.vcf.gz"),
                                chunk=CHUNKS))
    output:
        combined=out_path("multisample/genotyped.vcf.gz")
    conda: "envs/gatk.yml"
    shell: "java -Xmx4G -cp {input.gatk} org.broadinstitute.gatk.tools.CatVariants "
           "-R {input.ref} -V '{params.vcfs}' -out {output.combined} "
           "-assumeSorted"


## bam metrics

rule mapped_num:
    """Calculate number of mapped reads"""
    input:
        bam=out_path("{sample}/bams/{sample}.sorted.bam")
    output:
        num=out_path("{sample}/bams/{sample}.mapped.num")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | wc -l > {output.num}"


rule mapped_basenum:
    """Calculate number of mapped bases"""
    input:
        bam=out_path("{sample}/bams/{sample}.sorted.bam")
    output:
        num=out_path("{sample}/bams/{sample}.mapped.basenum")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 {input.bam} | cut -f10 | wc -c > {output.num}"


rule unique_num:
    """Calculate number of unique reads"""
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam")
    output:
        num=out_path("{sample}/bams/{sample}.unique.num")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | wc -l > {output.num}"


rule usable_basenum:
    """Calculate number of bases on unique reads"""
    input:
        bam=out_path("{sample}/bams/{sample}.markdup.bam")
    output:
        num=out_path("{sample}/bams/{sample}.usable.basenum")
    conda: "envs/samtools.yml"
    shell: "samtools view -F 4 -F 1024 {input.bam} | cut -f10 | wc -c > {output.num}"


## fastqc

rule fastqc_raw:
    """Run fastqc on raw fastq files"""
    input:
        r1=get_r1,
        r2=get_r2
    params:
        odir=out_path("{sample}/pre_process/raw_fastqc")
    output:
        aux=out_path("{sample}/pre_process/raw_fastqc/.done.txt")
    conda: "envs/fastqc.yml"
    shell: "fastqc --nogroup -o {params.odir} {input.r1} {input.r2} && echo 'done' > {output.aux}"


rule fastqc_merged:
    """Run fastqc on merged fastq files"""
    input:
        r1=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        r2=out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    params:
        odir=out_path("{sample}/pre_process/merged_fastqc")
    output:
        r1=out_path("{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip"),
        r2=out_path("{sample}/pre_process/merged_fastqc/{sample}.merged_R2_fastqc.zip")
    conda: "envs/fastqc.yml"
    shell: "fastqc --nogroup -o {params.odir} {input.r1} {input.r2}"


rule fastqc_postqc:
    """Run fastqc on fastq files post pre-processing"""
    input:
        r1=out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2=out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq")
    params:
        odir=out_path("{sample}/pre_process/postqc_fastqc")
    output:
        r1=out_path("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip"),
        r2=out_path("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R2_fastqc.zip")
    conda: "envs/fastqc.yml"
    shell: "fastqc --nogroup -o {params.odir} {input.r1} {input.r2}"


## fastq-count

rule fqcount_preqc:
    """Calculate number of reads and bases before pre-processing"""
    input:
        r1=out_path("{sample}/pre_process/{sample}.merged_R1.fastq.gz"),
        r2=out_path("{sample}/pre_process/{sample}.merged_R2.fastq.gz")
    params:
        fastqcount=fqc
    output:
        out_path("{sample}/pre_process/{sample}.preqc_count.json")
    shell: "{params.fastqcount} {input.r1} {input.r2} > {output}"


rule fqcount_postqc:
    """Calculate number of reads and bases after pre-processing"""
    input:
        r1=out_path("{sample}/pre_process/{sample}.cutadapt_R1.fastq"),
        r2=out_path("{sample}/pre_process/{sample}.cutadapt_R2.fastq")
    params:
        fastqcount=fqc
    output:
        out_path("{sample}/pre_process/{sample}.postqc_count.json")
    shell: "{params.fastqcount} {input.r1} {input.r2} > {output}"


# fastqc stats
rule fastqc_stats:
    """Collect fastq stats for a sample in json format"""
    input:
        preqc_r1=out_path("{sample}/pre_process/merged_fastqc/{sample}.merged_R1_fastqc.zip"),
        preqc_r2=out_path("{sample}/pre_process/merged_fastqc/{sample}.merged_R2_fastqc.zip"),
        postqc_r1=out_path("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R1_fastqc.zip"),
        postqc_r2=out_path("{sample}/pre_process/postqc_fastqc/{sample}.cutadapt_R2_fastqc.zip"),
        sc=fqpy
    conda: "envs/collectstats.yml"
    output:
        out_path("{sample}/pre_process/fastq_stats.json")
    shell: "python {input.sc} --preqc-r1 {input.preqc_r1} "
           "--preqc-r2 {input.preqc_r2} "
           "--postqc-r1 {input.postqc_r1} "
           "--postqc-r2 {input.postqc_r2} > {output}"

## coverages

rule covstats:
    """Calculate coverage statistics on bam file"""
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
    shell: "bedtools coverage -sorted -g {input.genome} -a {input.bed} -b {input.bam} "
           "-d  | python {input.covpy} - --plot {output.covp} "
           "--title 'Targets coverage' --subtitle '{params.subt}' > {output.covj}"


rule vtools_coverage:
    """Calculate coverage statistics per transcript"""
    input:
        gvcf=out_path("{sample}/vcf/{sample}.g.vcf.gz"),
        ref=get_refflatpath
    output:
        tsv=out_path("{sample}/coverage/{ref}.coverages.tsv")
    conda: "envs/vcfstats.yml"
    shell: "vtools-gcoverage -I {input.gvcf} -R {input.ref} > {output.tsv}"


## vcfstats

rule vcfstats:
    """Calculate vcf statistics"""
    input:
        vcf=out_path("multisample/genotyped.vcf.gz")
    output:
        stats=out_path("multisample/vcfstats.json")
    conda: "envs/vcfstats.yml"
    shell: "vtools-stats -i {input.vcf} > {output.stats}"


## collection
if len(BASE_BEDS) >= 1:
    rule collectstats:
        """Collect all stats for a particular sample with beds"""
        input:
            preqc=out_path("{sample}/pre_process/{sample}.preqc_count.json"),
            postq=out_path("{sample}/pre_process/{sample}.postqc_count.json"),
            mnum=out_path("{sample}/bams/{sample}.mapped.num"),
            mbnum=out_path("{sample}/bams/{sample}.mapped.basenum"),
            unum=out_path("{sample}/bams/{sample}.unique.num"),
            ubnum=out_path("{sample}/bams/{sample}.usable.basenum"),
            fastqc=out_path("{sample}/pre_process/fastq_stats.json"),
            cov=expand(out_path("{{sample}}/coverage/{bed}.covstats.json"),
                       bed=BASE_BEDS),
            colpy=colpy
        params:
            sample_name="{sample}",
            fthresh=FEMALE_THRESHOLD
        output:
            out_path("{sample}/{sample}.stats.json")
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
            preqc = out_path("{sample}/pre_process/{sample}.preqc_count.json"),
            postq = out_path("{sample}/pre_process/{sample}.postqc_count.json"),
            mnum = out_path("{sample}/bams/{sample}.mapped.num"),
            mbnum = out_path("{sample}/bams/{sample}.mapped.basenum"),
            unum = out_path("{sample}/bams/{sample}.unique.num"),
            ubnum = out_path("{sample}/bams/{sample}.usable.basenum"),
            fastqc=out_path("{sample}/pre_process/fastq_stats.json"),
            colpy = colpy
        params:
            sample_name = "{sample}",
            fthresh = FEMALE_THRESHOLD
        output:
            out_path("{sample}/{sample}.stats.json")
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
        cols=expand(out_path("{sample}/{sample}.stats.json"), sample=SAMPLES),
        vstat=out_path("multisample/vcfstats.json"),
        mpy=mpy
    output:
        stats=out_path("stats.json")
    conda: "envs/collectstats.yml"
    shell: "python {input.mpy} --vcfstats {input.vstat} {input.cols} > {output.stats}"


rule multiqc:
    """
    Create multiQC report
    Depends on stats.json to forcefully run at end of pipeline
    """
    input:
        stats=out_path("stats.json")
    params:
        odir=out_path(".")
    output:
        report=out_path("multiqc_report/multiqc_report.html")
    shell: "multiqc -o {output.report} {params.odir}"
