# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples. It is based on the
[Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/). Four distinct analysi modules
comprise Hamlet:

  1. `snv-indels`, for small variant detection
  2. `fusion`, for fusion gene detection
  3. `itd`, for tandem duplication detection
  4. `expression`, for expression analysis

There is also a `qc-seq` module which does quality control for input sequence files. Everything is tied together by
a main [`Snakefile`](https://snakemake.readthedocs.io/en/stable/executable.html#useful-command-line-arguments) using
the [`include`](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes) statement.

Configurations and input file definitions are supplied in a single YAML configuration file. Here you can specify
arbitrary number of samples, each containing an arbitrary number of paired-end files. The merging of these files and
separation of output per-sample are taken care of automatically by Hamlet. Additionally, if you mayalso supply your own
cluster-specific YAML configuration to determine how much resource such as memory or CPU cores should be allocated to
the jobs.

Hamlet tries to make use of isolated conda environment per Snakemake rule whenever possible. The base execution
environment is also defined by an `environment.yml` file.

In addition to the raw output files, Hamlet also generates a PDF report containing an overview of the essential results
and a zipped file containing this report and the essential result files.

The following diagram shows an overview of the Snakemake rules that are executed by Hamlet:

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.38 0.6 0.85", style="rounded"];
	1[label = "merge_fastqs_r2", color = "0.01 0.6 0.85", style="rounded"];
	2[label = "merge_fastqs_r1", color = "0.11 0.6 0.85", style="rounded"];
	3[label = "create_summary", color = "0.28 0.6 0.85", style="rounded"];
	4[label = "generate_report", color = "0.14 0.6 0.85", style="rounded"];
	5[label = "package_results", color = "0.57 0.6 0.85", style="rounded"];
	6[label = "reorder_aln_header", color = "0.36 0.6 0.85", style="rounded"];
	7[label = "annotate_vars", color = "0.60 0.6 0.85", style="rounded"];
	8[label = "table_vars_all", color = "0.29 0.6 0.85", style="rounded"];
	9[label = "table_vars_hi", color = "0.45 0.6 0.85", style="rounded"];
	10[label = "plot_vars", color = "0.55 0.6 0.85", style="rounded"];
	11[label = "star_fusion_cp", color = "0.46 0.6 0.85", style="rounded"];
	12[label = "plot_cp\next: star-fusion", color = "0.16 0.6 0.85", style="rounded"];
	13[label = "combine_plots", color = "0.33 0.6 0.85", style="rounded"];
	14[label = "fusioncatcher_cp", color = "0.06 0.6 0.85", style="rounded"];
	15[label = "plot_cp\next: fusioncatcher", color = "0.16 0.6 0.85", style="rounded"];
	16[label = "intersect_fusions", color = "0.30 0.6 0.85", style="rounded"];
	17[label = "plot_cp\next: sf-isect", color = "0.16 0.6 0.85", style="rounded"];
	18[label = "subset_sf", color = "0.63 0.6 0.85", style="rounded"];
	19[label = "count_fragments", color = "0.31 0.6 0.85", style="rounded"];
	20[label = "count_bases_gene", color = "0.07 0.6 0.85", style="rounded"];
	21[label = "count_bases_exon", color = "0.23 0.6 0.85", style="rounded"];
	22[label = "calc_exon_ratios", color = "0.47 0.6 0.85", style="rounded"];
	23[label = "sample_stats", color = "0.18 0.6 0.85", style="rounded"];
	24[label = "aln_stats", color = "0.17 0.6 0.85", style="rounded"];
	25[label = "rna_stats", color = "0.21 0.6 0.85", style="rounded"];
	26[label = "insert_stats", color = "0.39 0.6 0.85", style="rounded"];
	27[label = "exon_cov", color = "0.37 0.6 0.85", style="rounded"];
	28[label = "align_flt3", color = "0.48 0.6 0.85", style="rounded"];
	29[label = "detect_flt3", color = "0.56 0.6 0.85", style="rounded"];
	30[label = "plot_itd_flt3", color = "0.03 0.6 0.85", style="rounded"];
	31[label = "align_kmt2a", color = "0.54 0.6 0.85", style="rounded"];
	32[label = "detect_kmt2a", color = "0.44 0.6 0.85", style="rounded"];
	33[label = "plot_itd_kmt2a", color = "0.34 0.6 0.85", style="rounded"];
	34[label = "cutadapt\nread_group: rg_1\nsample: Sample1", color = "0.62 0.6 0.85", style="rounded"];
	35[label = "sort_bamfile", color = "0.22 0.6 0.85", style="rounded"];
	36[label = "call_vars", color = "0.52 0.6 0.85", style="rounded"];
	37[label = "extract_vars", color = "0.25 0.6 0.85", style="rounded"];
	38[label = "star_fusion", color = "0.10 0.6 0.85", style="rounded"];
	39[label = "plot_sf", color = "0.64 0.6 0.85", style="rounded"];
	40[label = "fusioncatcher", color = "0.43 0.6 0.85", style="rounded"];
	41[label = "plot_fc", color = "0.02 0.6 0.85", style="rounded"];
	42[label = "plot_isect", color = "0.13 0.6 0.85", style="rounded"];
	43[label = "idsort_aln", color = "0.05 0.6 0.85", style="rounded"];
	44[label = "count_raw_bases", color = "0.61 0.6 0.85", style="rounded"];
	45[label = "rg_stats", color = "0.08 0.6 0.85", style="rounded"];
	46[label = "exon_cov_ref", color = "0.24 0.6 0.85", style="rounded"];
	47[label = "genome_txt", color = "0.59 0.6 0.85", style="rounded"];
	48[label = "align_vars", color = "0.41 0.6 0.85", style="rounded"];
	49[label = "merge_fastqs_raw_r1\nsample: Sample1", color = "0.66 0.6 0.85", style="rounded"];
	50[label = "merge_fastqs_raw_r2\nsample: Sample1", color = "0.00 0.6 0.85", style="rounded"];
	51[label = "fastqc_raw\npair: R1\nread_group: rg_1\nsample: Sample1", color = "0.49 0.6 0.85", style="rounded"];
	52[label = "fastqc_raw\npair: R2\nread_group: rg_1\nsample: Sample1", color = "0.49 0.6 0.85", style="rounded"];
	53[label = "fastqc_processed\npair: R1", color = "0.09 0.6 0.85", style="rounded"];
	54[label = "fastqc_processed\npair: R2", color = "0.09 0.6 0.85", style="rounded"];
	1 -> 0
	2 -> 0
	3 -> 0
	4 -> 0
	5 -> 0
	6 -> 0
	7 -> 0
	8 -> 0
	9 -> 0
	10 -> 0
	11 -> 0
	12 -> 0
	13 -> 0
	14 -> 0
	15 -> 0
	16 -> 0
	17 -> 0
	18 -> 0
	19 -> 0
	20 -> 0
	21 -> 0
	22 -> 0
	23 -> 0
	24 -> 0
	25 -> 0
	26 -> 0
	27 -> 0
	28 -> 0
	29 -> 0
	30 -> 0
	31 -> 0
	32 -> 0
	33 -> 0
	34 -> 1
	34 -> 2
	23 -> 3
	24 -> 3
	25 -> 3
	26 -> 3
	7 -> 3
	27 -> 3
	10 -> 3
	9 -> 3
	13 -> 3
	30 -> 3
	33 -> 3
	29 -> 3
	32 -> 3
	22 -> 3
	3 -> 4
	3 -> 5
	8 -> 5
	9 -> 5
	10 -> 5
	13 -> 5
	19 -> 5
	20 -> 5
	21 -> 5
	22 -> 5
	29 -> 5
	30 -> 5
	32 -> 5
	33 -> 5
	4 -> 5
	35 -> 6
	36 -> 7
	37 -> 8
	37 -> 9
	37 -> 10
	38 -> 11
	39 -> 12
	12 -> 13
	15 -> 13
	17 -> 13
	40 -> 14
	41 -> 15
	11 -> 16
	14 -> 16
	42 -> 17
	11 -> 18
	16 -> 18
	43 -> 19
	44 -> 20
	44 -> 21
	21 -> 22
	45 -> 23
	6 -> 24
	6 -> 25
	6 -> 26
	6 -> 27
	46 -> 27
	47 -> 27
	2 -> 28
	1 -> 28
	28 -> 29
	29 -> 30
	2 -> 31
	1 -> 31
	31 -> 32
	32 -> 33
	48 -> 35
	6 -> 36
	7 -> 37
	2 -> 38
	1 -> 38
	11 -> 39
	49 -> 40
	50 -> 40
	14 -> 41
	18 -> 42
	6 -> 43
	6 -> 44
	51 -> 45
	52 -> 45
	53 -> 45
	54 -> 45
	2 -> 48
	1 -> 48
	34 -> 53
	34 -> 54
}
```

# Installation

The dependencies required for running the pipeline are listed in the provided `environment.yml` file. To use it, first
set up a Conda virtual environment and then update it:

```bash
# Set up and activate your conda environment.
# Install the dependencies
$ conda env update -f environment.yml
```

When running the pipeline, Snakemake will then install all the required tools via Conda, whenever possible.

Unfortunately, not all the required tools are available in Conda (or if they are present, they may not have been
compiled using certain optimizations). There are three tools that fall under this category, and you may need to do
additional steps to have them installed:

1. GSNAP, which is used by the alignment step before variant calling. On some certain architectures, the version of
   GSNAP that Hamlet uses can be compiled to utilize SIMD instructions. As of inclusion of this tool in the pipeline,
   the optimization is not available in Conda unfortunately. Here you may supply your own compiled GSNAP using the
   `gsnap_exe` settings. If you do, Hamlet will use that GSNAP executable. If you do not, Hamlet will use the
   GSNAP available from Conda and may require more time to run.

2. FusionCatcher, which is used by the fusion detection module. The results from fusioncatcher will be intersected
   with those of STAR-Fusion, and both the independend and the intersected results will be reported.

3. VEP, which is used to annotate the variant calling results. Hamlet was tested using an older version of VEP (77),
   which could not be made to work with Conda at the time of development. You will need to supply the path to the main
   executable via the `vep_exe` settings.

4. bgzip, which needs to be supplied via the `bgzip_exe` settings. This may be replaced with a proper Conda installation
   in the future.


# Usage

## Input files

Hamlet requires gzipped, paired-end mRNA-seq files. Any number of samples can be processed in a single execution, and
each of them may have differing number of pair sets. These files may have arbitrary names, but they *must* be supplied
correctly in the configuration files.

## Execution

After installation of all the required tools, you will need to fill in the required settings and configurations in
several YAML files.

For the runtime settings, use the provided `config-base.yml` file as template and fill in the required values as
instructed. Fill also your sample names and paths to their input files as instructed in the same YAML file. Let's call
this file `config.yml`.

If running in a cluster, you may also want to define the resource configurations in another YAML file. Read more about
this type of configuration on the official [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration). For this
file, let's call it `config-cluster.yml`

You can then run the pipeline by invoking Snakemake, for example:

```bash
    $ snakemake -s Snakefile \
        --configfile config.yml --cluster-config config-cluster.yml \
        --rerun-incomplete
        # ... other flags
```

Here is another example that uses more flags (see the Snakemake help for explanation of these flags):

```bash
    $ snakemake -p -T -s Snakefile \
        --configfile config.yml --cluster-config config-cluster.yml \
        --rerun-incomplete
        --jobname 'hamlet.{jobid}' \
        --jobs 100 -w 120 --max-jobs-per-second 3 \
        --drmaa ' -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V' \
        --drmaa-log-dir .drmaa-logs
```

## Output files

Assuming the output directory is set to `/path/to/output`, Hamlet will create `/path/to/output/{sample_name}` for each
sample present in the config file. Inside the directory, all the essential results are packaged in a zip file called
`hamlet_results.{sample_name}.zip`. This includes a PDF report called `hamlet_report.{sample_name}.pdf` which contains
the overview of the essential results.


## Notes

1. You can run Hamlet from anywhere, but preferrably this is done outside of the repository. This way, the temporary
   Snakemake files are written elsewhere and does not pollute the repository.

2. You can direct Hamlet to create the output directory anywhere. This is a configuration value that is supplied in the
   config file via `output_dir`.
