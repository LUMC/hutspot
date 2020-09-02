import itertools
import json
import jsonschema
import os

containers = {
   'bcftools': 'docker://quay.io/biocontainers/bcftools:1.9--ha228f0b_4',
   'bedtools-2.26-python-2.7': 'docker://quay.io/biocontainers/mulled-v2-3251e6c49d800268f0bc575f28045ab4e69475a6:4ce073b219b6dabb79d154762a9b67728c357edb-0',
   'biopet-scatterregions': 'docker://quay.io/biocontainers/biopet-scatterregions:0.2--0',
   'bwa-0.7.17-samtools-1.10': 'docker://quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0',
   'cutadapt': 'docker://quay.io/biocontainers/cutadapt:2.9--py37h516909a_0',
   'debian': 'docker://debian:buster-slim',
   'fastqc': 'docker://quay.io/biocontainers/fastqc:0.11.7--4',
   'gatk': 'docker://broadinstitute/gatk3:3.7-0',
   'gvcf2coverage': 'docker://lumc/gvcf2coverage:0.1-dirty-2',
   'multiqc': 'docker://quay.io/biocontainers/multiqc:1.8--py_2',
   'picard': 'docker://quay.io/biocontainers/picard:2.22.8--0',
   'python3': 'docker://python:3.6-slim',
   'vtools': 'docker://quay.io/biocontainers/vtools:1.0.0--py37h3010b51_0'
}

def process_config():
    """ Process the config file and set the default values """

    def set_default(key, value):
        """Set default config values"""
        if key not in config:
            config[key] = value

    # Read the json schema
    with open(srcdir('config/schema.json'), 'rt') as fin:
        schema = json.load(fin)

    # Validate the config against the schema
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError as e:
        raise jsonschema.ValidationError(f'Invalid --configfile: {e.message}')

    # If you specify a baitsfile, you also have to specify a targets file for
    # picard
    if 'baitsfile' in config and 'targetsfile' not in config:
        msg = 'Invalid --configfile: "baitsfile" specified without "targetsfile"'
        raise jsonschema.ValidationError(msg)

    # If you specify a target file but no baitsfile, we use the targets as
    # baits. This is needed because picard HsMetrics needs both a baitfile and
    # targets file as input
    if 'targetsfile' in config and 'baitsfile' not in config:
        set_default('baitsfile', config['targetsfile'])

    # A sample name cannot be a substring of another sample, since that breaks picard
    # metrics parsing by multiqc
    msg = 'Invalid --configfile: sample names should not overlap ("{s1}" is contained in "{s2}")'
    for s1, s2 in itertools.permutations(config['samples'], 2):
        if s1 in s2:
            raise jsonschema.ValidationError(msg.format(s1=s1, s2=s2))

    # Set the default config values
    set_default('scatter_size', 1000000000)
    set_default('female_threshold', 0.6)

    # Hide the absolute path so the snakemake linter doesn't cry about it
    set_default('gatk_jar', os.path.join(os.path.sep,'usr','GenomeAnalysisTK.jar'))

def coverage_stats(wildcards):
    files = expand("{sample}/coverage/refFlat_coverage.tsv",
                   sample=config["samples"])
    return files if "refflat" in config else []

def coverage_files(wildcards):
    """ Return a list of all coverage files

    The coverage is calculated for each sample, for each specified threshold
    """

    # We only calculate the coverage when this is specified in the
    # configuration
    if 'coverage_threshold' not in config:
        return list()

    # Fetch the values we need from the configuration
    samples = config['samples']
    thresholds = config['coverage_threshold']

    files = list()
    for sample, threshold in itertools.product(samples, thresholds):
        files.append(f'{sample}/vcf/{sample}_{threshold}.bed')
    return files

def sample_bamfiles(wildcards):
    """ Determine the bam files for a sample (one for each readgroup)
    """
    files = list()
    sample = config['samples'][wildcards.sample]
    sample_name = wildcards.sample
    for read_group in sample['read_groups']:
        files.append(f'{sample_name}/bams/{sample_name}-{read_group}.sorted.bam')
    return files

def gather_gvcf(wildcards):
    """ Gather the gvcf files based on the scatterregions checkpoint

    This is depends on the 'scatter_size' parameter and the reference genome
    used
    """
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.g.vcf.gz",
       i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

def gather_gvcf_tbi(wildcards):
    """ Gather the gvcf index files based on the scatterregions checkpoint
    This is depends on the 'scatter_size' parameter and the reference genome
    used
    """
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.g.vcf.gz.tbi",
       i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

def gather_vcf(wildcards):
    """ Gather the vcf files based on the scatterregions checkpoint
    This is depends on the 'scatter_size' parameter and the reference genome
    used
    """
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.vcf.gz",
       i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

def gather_vcf_tbi(wildcards):
    """ Gather the vcf index files based on the scatterregions checkpoint
    This is depends on the 'scatter_size' parameter and the reference genome
    used
    """
    checkpoint_output = checkpoints.scatterregions.get(**wildcards).output[0]
    return expand("{{sample}}/vcf/{{sample}}.{i}.vcf.gz.tbi",
       i=glob_wildcards(os.path.join(checkpoint_output, 'scatter-{i}.bed')).i)

def sample_cutadapt_files(wildcards):
    """ Determine the cutadapt log files files for a sample (one for each
    readgroup).
    """
    files = list()
    sample = config['samples'][wildcards.sample]
    sample_name = wildcards.sample
    for read_group in sample['read_groups']:
        files.append(f'{sample_name}/pre_process/{sample_name}-{read_group}.txt')
    return files

def all_trimmed_fastqc(wildcards):
    """ Determine the trimmed fastq files for each sample """
    fastq_files = list()
    for sample in config['samples']:
        for read_group in config['samples'][sample]['read_groups']:
            fastq_files.append(f"{sample}/pre_process/trimmed-{sample}-{read_group}/.done")
    return fastq_files
