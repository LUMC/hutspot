Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v2.0.1
---------------------------
+ `multisample_vcf` now acts on the scatters, instead of on the merged g.vcf
files.
+ The multisample output is located in `multisample/multisample.vcf.gz`.
+ Intermediate .bam, .bai and fastq files are automatically removed when no
longer needed.
+ Switch to using chunked-scatter

v2.0.0
---------------------------
+ Add an environment.yml file for conda.
+ Greatly simplified the snakemake workflow.
+ All statistics are now calculated using existing tools.
+ Add option `multisample_vcf` to enable joint variantcalling.