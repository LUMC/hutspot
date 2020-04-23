The `micro_R1.fq.gz` and `micro_R2.fq.gz` fastq files were generated using 
data downloaded from ENA with project number
 [ERR2438055](https://www.ebi.ac.uk/ena/data/view/ERR2438055).

The full ENA dataset was aligned to hg19, after which about 150k reads 
aligning to the mitochondrial genome were selected and converted back to
fastq.

The `ref.fa` fasta file contains human mitochondrial genome (NC_012920).

The `database.vcf.gz` VCF file contains a single artificial variant. 

The `micro_rg1_R1.fq.gz` and `micro_rg1_R2.fq.gz` were created by taking the
first 15440 lines from `micro_R1.fq.gz` and `micro_R2.fq.gz, respectively.

The `micro_rg2_R1.fq.gz` and `micro_rg2_R2.fq.gz` were created by taking the
last 15440 lines from `micro_R1.fq.gz` and `micro_R2.fq.gz, respectively.
