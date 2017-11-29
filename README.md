# hutspot


# graph

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.32 0.6 0.85", style="rounded"];
	1[label = "fastqc_postqc", color = "0.48 0.6 0.85", style="rounded"];
	2[label = "fastqc_raw", color = "0.05 0.6 0.85", style="rounded"];
	3[label = "fastqc_merged", color = "0.57 0.6 0.85", style="rounded"];
	4[label = "merge_stats", color = "0.21 0.6 0.85", style="rounded"];
	5[label = "genotype_gather", color = "0.37 0.6 0.85", style="rounded"];
	6[label = "cutadapt", color = "0.30 0.6 0.85", style="rounded"];
	7[label = "merge_r2", color = "0.60 0.6 0.85", style="rounded"];
	8[label = "merge_r1", color = "0.23 0.6 0.85", style="rounded"];
	9[label = "collectstats", color = "0.02 0.6 0.85", style="rounded"];
	10[label = "vcfstats", color = "0.07 0.6 0.85", style="rounded"];
	11[label = "genotype_scatter", color = "0.16 0.6 0.85", style="rounded"];
	12[label = "sickle", color = "0.46 0.6 0.85", style="rounded"];
	13[label = "mapped_num", color = "0.64 0.6 0.85", style="rounded"];
	14[label = "fqcount_postqc", color = "0.00 0.6 0.85", style="rounded"];
	15[label = "unique_num", color = "0.11 0.6 0.85", style="rounded"];
	16[label = "covstats", color = "0.41 0.6 0.85", style="rounded"];
	17[label = "fqcount_preqc", color = "0.62 0.6 0.85", style="rounded"];
	18[label = "usable_basenum", color = "0.14 0.6 0.85", style="rounded"];
	19[label = "mapped_basenum", color = "0.09 0.6 0.85", style="rounded"];
	20[label = "gvcf_gather", color = "0.55 0.6 0.85", style="rounded"];
	21[label = "seqtk_r2", color = "0.51 0.6 0.85", style="rounded"];
	22[label = "seqtk_r1", color = "0.28 0.6 0.85", style="rounded"];
	23[label = "align", color = "0.25 0.6 0.85", style="rounded"];
	24[label = "markdup", color = "0.18 0.6 0.85", style="rounded"];
	25[label = "genome", color = "0.53 0.6 0.85", style="rounded"];
	26[label = "gvcf_scatter", color = "0.34 0.6 0.85", style="rounded"];
	27[label = "printreads", color = "0.39 0.6 0.85", style="rounded"];
	28[label = "baserecal", color = "0.44 0.6 0.85", style="rounded"];
	4 -> 0
	2 -> 0
	1 -> 0
	5 -> 0 [color="red"]
	3 -> 0
	6 -> 1
	8 -> 3
	7 -> 3
	10 -> 4
	9 -> 4
	11 -> 5 [color="red"]
	12 -> 6 [color="red"]
	14 -> 9
	19 -> 9
	18 -> 9
	16 -> 9
	15 -> 9
	17 -> 9
	13 -> 9
	5 -> 10
	20 -> 11 [color="red"]
	22 -> 12 [color="red"]
	21 -> 12 [color="red"]
	23 -> 13
	6 -> 14
	24 -> 15
	24 -> 16
	25 -> 16
	8 -> 17 [color="red"]
	7 -> 17 [color="red"]
	24 -> 18
	23 -> 19
	26 -> 20 [color="red"]
	17 -> 21 [color="red"]
	7 -> 21 [color="red"]
	17 -> 22 [color="red"]
	8 -> 22 [color="red"]
	6 -> 23 [color="red"]
	23 -> 24 [color="red"]
	27 -> 26 [color="red"]
	24 -> 27 [color="red"]
	28 -> 27 [color="red"]
	24 -> 28 [color="red"]
}  
```