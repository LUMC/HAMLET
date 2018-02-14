# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples.


## Running

    snakemake -p -T --use-conda -s Snakefile --configfile config.yml \
        --rerun-incomplete --restart-times 3 \
        --cluster-config config-cluster.yml --jobname 'hamlet.{jobid}' \
        --jobs 100 -w 120 --max-jobs-per-second 3 \
        --drmaa ' -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V' \
        --drmaa-log-dir .drmaa-logs


## Steps

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    rankdir=LR;
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "combine_plots", color = "0.06 0.6 0.85", style="rounded"];
	1[label = "subset_sf", color = "0.01 0.6 0.85", style="rounded"];
	2[label = "intersect_fusions", color = "0.03 0.6 0.85", style="rounded"];
	3[label = "plot_cp", color = "0.04 0.6 0.85", style="rounded"];
	4[label = "count_bases_gene", color = "0.59 0.6 0.85", style="rounded"];
	5[label = "table_vars_hi", color = "0.10 0.6 0.85", style="rounded"];
	6[label = "count_fragments", color = "0.09 0.6 0.85", style="rounded"];
	7[label = "plot_vars", color = "0.41 0.6 0.85", style="rounded"];
	8[label = "sample_stats", color = "0.12 0.6 0.85", style="rounded"];
	9[label = "detect_flt3", color = "0.13 0.6 0.85", style="rounded"];
	10[label = "plot_sf", color = "0.58 0.6 0.85", style="rounded"];
	11[label = "align_flt3", color = "0.43 0.6 0.85", style="rounded"];
	12[label = "rg_stats", color = "0.39 0.6 0.85", style="rounded"];
	13[label = "plot_itd_flt3", color = "0.44 0.6 0.85", style="rounded"];
	14[label = "plot_itd_kmt2a", color = "0.62 0.6 0.85", style="rounded"];
	15[label = "fusioncatcher", color = "0.15 0.6 0.85", style="rounded"];
	16[label = "plot_fc", color = "0.46 0.6 0.85", style="rounded"];
	17[label = "fastqc_processed", color = "0.47 0.6 0.85", style="rounded"];
	18[label = "align_vars", color = "0.16 0.6 0.85", style="rounded"];
	19[label = "rna_stats", color = "0.18 0.6 0.85", style="rounded"];
	20[label = "clip_trim_sync", color = "0.19 0.6 0.85", style="rounded"];
	21[label = "call_vars", color = "0.49 0.6 0.85", style="rounded"];
	22[label = "align_kmt2a", color = "0.21 0.6 0.85", style="rounded"];
	23[label = "star_fusion_cp", color = "0.22 0.6 0.85", style="rounded"];
	24[label = "table_vars_all", color = "0.50 0.6 0.85", style="rounded"];
	25[label = "merge_fastqs_raw", color = "0.52 0.6 0.85", style="rounded"];
	26[label = "aln_stats", color = "0.31 0.6 0.85", style="rounded"];
	27[label = "fusioncatcher_cp", color = "0.53 0.6 0.85", style="rounded"];
	28[label = "count_bases_exon", color = "0.25 0.6 0.85", style="rounded"];
	29[label = "star_fusion", color = "0.27 0.6 0.85", style="rounded"];
	30[label = "count_raw_bases", color = "0.24 0.6 0.85", style="rounded"];
	31[label = "insert_stats", color = "0.30 0.6 0.85", style="rounded"];
	32[label = "plot_isect", color = "0.55 0.6 0.85", style="rounded"];
	33[label = "reorder_aln_header", color = "0.00 0.6 0.85", style="rounded"];
	34[label = "annotate_vars", color = "0.56 0.6 0.85", style="rounded"];
	35[label = "detect_kmt2a", color = "0.34 0.6 0.85", style="rounded"];
	36[label = "fastqc_raw", color = "0.07 0.6 0.85", style="rounded"];
	37[label = "all", color = "0.36 0.6 0.85", style="rounded"];
	38[label = "merge_fastqs", color = "0.37 0.6 0.85", style="rounded"];
	39[label = "extract_vars", color = "0.61 0.6 0.85", style="rounded"];
	40[label = "calc_exon_ratios", color = "0.64 0.6 0.85", style="rounded"];
	41[label = "idsort_aln", color = "0.65 0.6 0.85", style="rounded"];
	3 -> 0
	2 -> 1
	23 -> 1
	27 -> 2
	23 -> 2
	32 -> 3
	16 -> 3
	10 -> 3
	30 -> 4
	39 -> 5
	41 -> 6
	39 -> 7
	12 -> 8
	11 -> 9
	23 -> 10
	38 -> 11
	17 -> 12
	36 -> 12
	9 -> 13
	35 -> 14
	25 -> 15
	27 -> 16
	20 -> 17
	38 -> 18
	33 -> 19
	36 -> 20
	33 -> 21
	38 -> 22
	29 -> 23
	39 -> 24
	33 -> 26
	15 -> 27
	30 -> 28
	38 -> 29
	33 -> 30
	33 -> 31
	1 -> 32
	18 -> 33
	21 -> 34
	22 -> 35
	0 -> 37
	1 -> 37
	2 -> 37
	3 -> 37
	33 -> 37
	6 -> 37
	5 -> 37
	8 -> 37
	9 -> 37
	19 -> 37
	22 -> 37
	23 -> 37
	28 -> 37
	31 -> 37
	26 -> 37
	35 -> 37
	38 -> 37
	7 -> 37
	11 -> 37
	13 -> 37
	24 -> 37
	27 -> 37
	34 -> 37
	14 -> 37
	4 -> 37
	40 -> 37
	20 -> 38
	34 -> 39
	28 -> 40
	33 -> 41
}
```
