#   hutspot - a DNAseq variant calling pipeline
#   Copyright (C) 2017-2019, Leiden University Medical Center
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

include: "common.smk"

process_config()

localrules:
    collectstats,
    create_markdup_tmp,
    cutadapt_summary,
    merge_stats,
    multiqc,
    stats_tsv

rule all:
    input:
        multiqc = "multiqc_report/multiqc_report.html",
        stats = "stats.json",
        stats_tsv = "stats.tsv",
        bam = expand("{s}/bams/{s}.bam", s=config["samples"]),
        vcfs = expand("{s}/vcf/{s}.vcf.gz", s=config["samples"]),
        vcf_tbi = expand("{s}/vcf/{s}.vcf.gz.tbi", s=config["samples"]),
        gvcfs = expand("{s}/vcf/{s}.g.vcf.gz", s=config["samples"]),
        gvcf_tbi = expand("{s}/vcf/{s}.g.vcf.gz.tbi", s=config["samples"]),
        coverage_stats = coverage_stats,
        coverage_files = coverage_files

rule create_markdup_tmp:
    """Create tmp directory for mark duplicates"""
    output:
        directory("tmp")
    log:
        "log/create_markdup_tmp.log"
    container:
        containers["debian"]
    shell:
        "mkdir -p {output} 2> {log}"

rule genome:
    """Create genome file as used by bedtools"""
    input:
        config["reference"]
    output:
        "current.genome"
    log:
        "log/genome.log"
    container:
        containers["debian"]
    shell:
        "awk -v OFS='\t' {{'print $1,$2'}} {input}.fai > {output} 2> {log}"

rule cutadapt:
    """Clip fastq files"""
    input:
        r1 = lambda wc: (config['samples'][wc.sample]['read_groups'][wc.read_group]['R1']),
        r2 = lambda wc: (config['samples'][wc.sample]['read_groups'][wc.read_group]['R2'])
    output:
        r1 = "{sample}/pre_process/{sample}-{read_group}_R1.fastq.gz",
        r2 = "{sample}/pre_process/{sample}-{read_group}_R2.fastq.gz"
    log:
        "{sample}/pre_process/{sample}-{read_group}.txt"
    container:
        containers["cutadapt"]
    threads:
        8
    shell:
        "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG "
        "--minimum-length 1 --quality-cutoff=20,20 "
        "--compression-level=1 -j 8 "
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
    output:
        "{sample}/bams/{sample}-{read_group}.sorted.bam"
    params:
        compression_level = 1,
        rg = "@RG\\tID:{sample}-library-{read_group}\\tSM:{sample}\\tLB:library\\tPL:ILLUMINA"

    log:
        bwa = "log/{sample}/align.{read_group}.bwa.log",
        samtools = "log/{sample}/align.{read_group}.samtools.log"
    container:
        containers["bwa-0.7.17-samtools-1.10"]
    threads:
        8
    shell:
        "set -eo pipefail;"
        "bwa mem -t {threads} -R '{params.rg}' {input.ref} "
        "{input.r1} {input.r2} 2> {log.bwa} | "
        "samtools sort "
        "-l {params.compression_level} "
        "- -o {output} 2> {log.samtools};"
        "samtools index {output}"

rule markdup:
    """Mark duplicates in BAM file"""
    input:
        bam = sample_bamfiles,
        tmp = ancient("tmp")
    output:
        bam = "{sample}/bams/{sample}.bam",
        bai = "{sample}/bams/{sample}.bai",
        metrics = "{sample}/bams/{sample}.metrics"
    params:
        bams = lambda wc: expand("INPUT={bam}", bam=sample_bamfiles(wc))
    log:
        "log/{sample}/markdup.log"
    container:
        containers["picard"]
    shell:
        "picard -Xmx4G -Djava.io.tmpdir={input.tmp} MarkDuplicates "
        "CREATE_INDEX=TRUE TMP_DIR={input.tmp} "
        "{params.bams} OUTPUT={output.bam} "
        "METRICS_FILE={output.metrics} "
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 2> {log}"

rule baserecal:
    """Base recalibrated BAM files"""
    input:
        bam = sample_bamfiles,
        ref = config["reference"],
        vcfs = config["known_sites"]
    output:
        "{sample}/bams/{sample}.baserecal.grp"
    params:
        known_sites = expand("-knownSites {vcf}", vcf=config["known_sites"]),
        region = f"-L {config['restrict_BQSR']}" if "restrict_BQSR" in config else "",
        gatk_jar = config["gatk_jar"],
        bams = lambda wc: expand("-I {bam}", bam=sample_bamfiles(wc))
    log:
        "log/{sample}/baserecal.log"
    container:
        containers["gatk"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 6000,
    threads:
        8
    shell:
        "java -XX:ParallelGCThreads=1 -jar {params.gatk_jar} -T "
        "BaseRecalibrator {params.bams} -o {output} -nct {threads} "
        "-R {input.ref} -cov ReadGroupCovariate -cov QualityScoreCovariate "
        "-cov CycleCovariate -cov ContextCovariate {params.known_sites} "
        "{params.region} 2> {log}"

checkpoint scatterregions:
    """Scatter the reference genome"""
    input:
        ref = config["reference"],
    output:
        directory("scatter")
    params:
        size = config["scatter_size"]
    log:
        "log/scatterregions.log"
    container:
        containers["biopet-scatterregions"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000
    shell:
        "mkdir -p {output} && "
        "biopet-scatterregions -Xmx24G "
        "--referenceFasta {input.ref} --scatterSize {params.size} "
        "--outputDir scatter 2> {log}"

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
    params:
        gatk_jar = config["gatk_jar"],
    log:
        "log/{sample}/gvcf_scatter.{chunk}.log"
    container:
        containers["gatk"]
    shell:
        "java -jar -Xmx4G -XX:ParallelGCThreads=1 {params.gatk_jar} -T "
        "HaplotypeCaller -ERC GVCF -I "
        "{input.bam} -R {input.ref} -D {input.dbsnp} "
        "-L '{input.region}' -o '{output.gvcf}' "
        "-variant_index_type LINEAR -variant_index_parameter 128000 "
        "-BQSR {input.bqsr} "
        "--GVCFGQBands 20 --GVCFGQBands 40 --GVCFGQBands 60 "
        "--GVCFGQBands 80 --GVCFGQBands 100 2> {log}"

rule gvcf_gather:
    """Join the gvcf files together"""
    input:
        gvcfs = gather_gvcf,
        tbis = gather_gvcf_tbi
    output:
        gvcf = "{sample}/vcf/{sample}.g.vcf.gz",
        gvcf_tbi = "{sample}/vcf/{sample}.g.vcf.gz.tbi"
    log:
        concat = "log/{sample}/gvcf_gather_concat.log",
        index = "log/{sample}/gvcf_gather_index.log"
    container:
        containers["bcftools"]
    shell:
        "bcftools concat {input.gvcfs} --allow-overlaps "
        "--output {output.gvcf} --output-type z 2> {log.concat} && "
        "bcftools index --tbi --output-file {output.gvcf_tbi} "
        "{output.gvcf} 2> {log.index}"

rule genotype_scatter:
    """Run GATK's GenotypeGVCFs by chunk"""
    input:
        gvcf = rules.gvcf_scatter.output.gvcf,
        tbi = rules.gvcf_scatter.output.gvcf_tbi,
        ref = config["reference"]
    output:
        vcf = temp("{sample}/vcf/{sample}.{chunk}.vcf.gz"),
        vcf_tbi = temp("{sample}/vcf/{sample}.{chunk}.vcf.gz.tbi")
    params:
        gatk_jar = config["gatk_jar"]
    log:
        "log/{sample}/genotype_scatter.{chunk}.log"
    container:
        containers["gatk"]
    resources:
        mem_mb = 6000
    shell:
        "java -jar -Xmx15G -XX:ParallelGCThreads=1 {params.gatk_jar} -T "
        "GenotypeGVCFs -R {input.ref} "
        "-V {input.gvcf} -o {output.vcf} 2> {log}"

rule genotype_gather:
    """Gather all genotyping VCFs"""
    input:
        vcfs = gather_vcf,
        vcfs_tbi = gather_vcf_tbi
    output:
        vcf = "{sample}/vcf/{sample}.vcf.gz",
        vcf_tbi = "{sample}/vcf/{sample}.vcf.gz.tbi"
    log:
        "log/{sample}/genotype_gather.log"
    container:
        containers["bcftools"]
    shell:
        "bcftools concat {input.vcfs} --allow-overlaps "
        "--output {output.vcf} --output-type z 2> {log} && "
        "bcftools index --tbi --output-file {output.vcf_tbi} {output.vcf}"

rule fastqc:
    """Run fastqc on fastq files post pre-processing"""
    input:
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2
    output:
        done = "{sample}/pre_process/trimmed-{sample}-{read_group}/.done"
    log:
        "log/{sample}/fastqc.{read_group}.log"
    container:
        containers["fastqc"]
    threads:
        4
    shell:
        "fastqc --threads {threads} --nogroup -o $(dirname {output.done}) "
        "{input.r1} {input.r2} 2> {log} && "
        "touch {output.done}"

rule covstats:
    """Calculate coverage statistics on bam file"""
    input:
        bam = rules.markdup.output.bam,
        genome = "current.genome",
        covstats = srcdir("src/covstats.py"),
        bed = config.get("targetsfile", "")
    output:
        covj = "{sample}/coverage/covstats.json",
        covp = "{sample}/coverage/covstats.png"
    params:
        subt = "Sample {sample}"
    log:
        bedtools = "log/{sample}/covstats_bedtools.log",
        covstats = "log/{sample}/covstats_covstats.log"
    container:
        containers["bedtools-2.26-python-2.7"]
    resources:
        mem_mb = 6000
    threads:
        2
    shell:
        "bedtools coverage -sorted -g {input.genome} -a {input.bed} "
        "-b {input.bam} -d  2> {log.bedtools} | python {input.covstats} - "
        "--plot {output.covp} --title 'Targets coverage' "
        "--subtitle '{params.subt}' > {output.covj} 2> {log.covstats}"

rule vtools_coverage:
    """Calculate coverage statistics per transcript"""
    input:
        gvcf = rules.gvcf_gather.output.gvcf,
        tbi = rules.gvcf_gather.output.gvcf_tbi,
        ref = config.get("refflat", "")
    output:
        "{sample}/coverage/refFlat_coverage.tsv"
    log:
        "log/{sample}/vtools_coverage.log"
    container:
        containers["vtools"]
    shell:
        "vtools-gcoverage -I {input.gvcf} "
        "-R {input.ref} > {output} 2> {log}"

rule cutadapt_summary:
    """Colect cutadapt summary from each readgroup per sample """
    input:
        cutadapt = sample_cutadapt_files,
        cutadapt_summary = srcdir("src/cutadapt_summary.py")
    output:
        "{sample}/cutadapt.json"
    log:
        "log/{sample}/cutadapt_summary.log"
    container:
        containers["python3"]
    shell:
        "python {input.cutadapt_summary} --sample {wildcards.sample} "
        "--cutadapt-summary {input.cutadapt} > {output}"

rule collectstats:
    """Collect all stats for a particular sample"""
    input:
        cov = rules.covstats.output.covj if "targetsfile" in config else [],
        cutadapt = rules.cutadapt_summary.output,
        collect_stats = srcdir("src/collect_stats.py")
    output:
        "{sample}/{sample}.stats.json"
    params:
        fthresh = config["female_threshold"]
    log:
        "log/{sample}/collectstats.log"
    container:
        containers["python3"]
    shell:
        "python {input.collect_stats} --sample-name {wildcards.sample} "
        "--female-threshold {params.fthresh} "
        "--cutadapt {input.cutadapt} "
        "--covstats {input.cov} > {output} 2> {log}"

rule multiple_metrics:
    """Run picard CollectMultipleMetrics"""
    input:
        bam = rules.markdup.output.bam,
        ref = config["reference"]
    output:
        alignment = "{sample}/bams/{sample}.alignment_summary_metrics",
        insertMetrics = "{sample}/bams/{sample}.insert_size_metrics"
    params:
        prefix = lambda wildcards, output: output.alignment[:-26]
    log:
        "log/{sample}/multiple_metrics.log"
    container:
        containers["picard"]
    shell:
        "picard -Xmx8G -XX:CompressedClassSpaceSize=256m CollectMultipleMetrics "
        "I={input.bam} O={params.prefix} "
        "R={input.ref} "
        "PROGRAM=CollectAlignmentSummaryMetrics "
        "PROGRAM=CollectInsertSizeMetrics 2> {log}"

rule bed_to_interval:
    """Run picard BedToIntervalList

    This is needed to convert the bed files for the capture kit to a format
    picard can read
    """
    input:
        targets = config.get("targetsfile", ""),
        baits = config.get("baitsfile", ""),
        ref = config["reference"]
    output:
        target_interval = "target.interval",
        baits_interval = "baits.interval"
    log:
        target = "log/bed_to_interval_target.log",
        baits = "log/bed_to_interval_target.log"
    container:
        containers["picard"]
    shell:
        "picard -Xmx4G BedToIntervalList "
        "I={input.targets} O={output.target_interval} "
        "SD={input.ref} 2> {log.target} && "
        "picard BedToIntervalList "
        "I={input.baits} O={output.baits_interval} "
        "SD={input.ref} 2> {log.baits}"

rule hs_metrics:
    """Run picard CollectHsMetrics"""
    input:
        bam = rules.markdup.output.bam,
        ref = config["reference"],
        targets = rules.bed_to_interval.output.target_interval,
        baits = rules.bed_to_interval.output.baits_interval
    output:
        "{sample}/bams/{sample}.hs_metrics.txt"
    log:
        "log/{sample}/hs_metrics.log"
    container:
        containers["picard"]
    shell:
        "picard -Xmx8G -XX:CompressedClassSpaceSize=256m CollectHsMetrics "
        "I={input.bam} O={output} "
        "R={input.ref} "
        "BAIT_INTERVALS={input.baits} "
        "TARGET_INTERVALS={input.targets}"

rule multiqc:
    """
    Create multiQC report
    Depends on stats.tsv to forcefully run at end of pipeline
    """
    input:
        bam = expand("{s}/bams/{s}.bam", s=config["samples"]),
        metric = expand("{s}/bams/{s}.metrics", s=config["samples"]),
        alignment_metrics = expand("{s}/bams/{s}.alignment_summary_metrics", s=config["samples"]),
        insert_metrics = expand("{s}/bams/{s}.insert_size_metrics", s=config["samples"]),
        fastqc = all_trimmed_fastqc,
        hs_metric = expand("{s}/bams/{s}.hs_metrics.txt", s=config["samples"]) if "baitsfile" in config else []
    output:
        html = "multiqc_report/multiqc_report.html",
        insertSize = "multiqc_report/multiqc_data/multiqc_picard_insertSize.json",
        AlignmentMetrics = "multiqc_report/multiqc_data/multiqc_picard_AlignmentSummaryMetrics.json",
        DuplicationMetrics = "multiqc_report/multiqc_data/multiqc_picard_dups.json",
        HsMetrics = "multiqc_report/multiqc_data/multiqc_picard_HsMetrics.json" if "baitsfile" in config else []
    log:
        "log/multiqc.log"
    container:
        containers["multiqc"]
    shell:
        "multiqc --data-format json --force "
        "--outdir multiqc_report . 2> {log} "
        "|| touch {output}"

rule merge_stats:
    """Merge all stats of all samples"""
    input:
        cols = expand("{sample}/{sample}.stats.json", sample=config["samples"]),
        merge_stats = srcdir("src/merge_stats.py"),
        insertSize = rules.multiqc.output.insertSize,
        AlignmentMetrics = rules.multiqc.output.AlignmentMetrics,
        DuplicationMetrics = rules.multiqc.output.DuplicationMetrics,
        HsMetrics = rules.multiqc.output.HsMetrics
    output:
        "stats.json"
    log:
        "log/stats.tsv.log"
    container:
        containers["vtools"]
    shell:
        "python {input.merge_stats} --collectstats {input.cols} "
        "--picard-insertSize {input.insertSize} "
        "--picard-AlignmentMetrics {input.AlignmentMetrics} "
        "--picard-DuplicationMetrics {input.DuplicationMetrics} "
        "--picard-HsMetrics {input.HsMetrics} > {output} 2> {log}"

rule stats_tsv:
    """Convert stats.json to tsv"""
    input:
        stats = rules.merge_stats.output,
        stats_to_tsv = srcdir("src/stats_to_tsv.py")
    output:
        "stats.tsv"
    log:
        "log/stats.tsv.log"
    container:
        containers["python3"]
    shell:
        "python {input.stats_to_tsv} -i {input.stats} > {output} 2> {log}"

rule gvcf2coverage:
    """ Determine coverage from gvcf files """
    input:
        rules.gvcf_gather.output.gvcf
    output:
        "{sample}/vcf/{sample}_{threshold}.bed"
    log:
        "log/{sample}/gvcf2coverage.{threshold}.log"
    container:
        containers["gvcf2coverage"]
    shell:
        "gvcf2coverage -t {wildcards.threshold} < {input} 2> {log} | cut -f 1,2,3 > {output}"
