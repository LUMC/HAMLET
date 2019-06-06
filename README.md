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
    rankdir=LR;
    node[fontname=sans, fontsize=10, label="\N", penwidth=2, shape=box, style=rounded];
    edge[color=grey, penwidth=2];
    0[color="0.13 0.6 0.85", label=all];
    1[color="0.38 0.6 0.85", label=fusioncatcher_cp];
    10[color="0.05 0.6 0.85", label=intersect_fusions];
    34[color="0.29 0.6 0.85", label=plot_fc];
    2[color="0.21 0.6 0.85", label=annotate_vars];
    19[color="0.44 0.6 0.85", label=create_summary];
    38[color="0.26 0.6 0.85", label=extract_vars];
    3[color="0.09 0.6 0.85", label=count_bases_exon];
    13[color="0.60 0.6 0.85", label=package_results];
    25[color="0.59 0.6 0.85", label=calc_exon_ratios];
    4[color="0.47 0.6 0.85", label=detect_flt3];
    28[color="0.63 0.6 0.85", label=plot_itd_flt3];
    5[color="0.16 0.6 0.85", label=plot_cp];
    17[color="0.54 0.6 0.85", label=combine_plots];
    6[color="0.52 0.6 0.85", label=star_fusion_cp];
    16[color="0.56 0.6 0.85", label=subset_sf];
    39[color="0.27 0.6 0.85", label=plot_sf];
    7[color="0.50 0.6 0.85", label=generate_report];
    8[color="0.08 0.6 0.85", label=align_flt3];
    9[color="0.24 0.6 0.85", label=insert_stats];
    1 -> 0
    1 -> 10
    1 -> 34
    2 -> 0
    2 -> 19
    2 -> 38
    3 -> 0
    3 -> 13
    3 -> 25
    4 -> 0
    4 -> 13
    4 -> 28
    5 -> 0
    5 -> 17
    6 -> 0
    6 -> 10
    6 -> 16
    6 -> 39
    7 -> 0
    7 -> 13
    8 -> 0
    8 -> 4
    9 -> 0
    9 -> 19
    10 -> 0
    10 -> 16
    11[color="0.34 0.6 0.85", label=align_kmt2a];
    11 -> 0
    12[color="0.22 0.6 0.85", label=detect_kmt2a];
    11 -> 12
    12 -> 0
    12 -> 13
    23[color="0.12 0.6 0.85", label=plot_itd_kmt2a];
    12 -> 23
    13 -> 0
    14[color="0.35 0.6 0.85", label=rna_stats];
    14 -> 0
    14 -> 19
    15[color="0.04 0.6 0.85", label=merge_fastqs];
    15 -> 0
    15 -> 8
    15 -> 11
    35[color="0.55 0.6 0.85", label=star_fusion];
    15 -> 35
    44[color="0.18 0.6 0.85", label=align_vars];
    15 -> 44
    16 -> 0
    41[color="0.17 0.6 0.85", label=plot_isect];
    16 -> 41
    17 -> 0
    17 -> 13
    17 -> 19
    18[color="0.61 0.6 0.85", label=aln_stats];
    18 -> 0
    18 -> 19
    19 -> 0
    19 -> 7
    19 -> 13
    20[color="0.14 0.6 0.85", label=sample_stats];
    20 -> 0
    20 -> 19
    21[color="0.65 0.6 0.85", label=table_vars_all];
    21 -> 0
    21 -> 13
    22[color="0.00 0.6 0.85", label=table_vars_hi];
    22 -> 0
    22 -> 13
    23 -> 0
    23 -> 13
    23 -> 19
    24[color="0.58 0.6 0.85", label=plot_vars];
    24 -> 0
    24 -> 13
    24 -> 19
    25 -> 0
    25 -> 13
    25 -> 19
    26[color="0.43 0.6 0.85", label=count_bases_gene];
    26 -> 0
    26 -> 13
    27[color="0.07 0.6 0.85", label=count_fragments];
    27 -> 0
    27 -> 13
    28 -> 0
    28 -> 13
    28 -> 19
    29[color="0.30 0.6 0.85", label=exon_cov];
    30[color="0.64 0.6 0.85", label=reorder_aln_header];
    29 -> 0
    29 -> 19
    30 -> 0
    30 -> 9
    30 -> 14
    30 -> 18
    30 -> 29
    30 -> 32
    30 -> 33
    30 -> 40
    31 -> 1
    32 -> 2
    33 -> 3
    33 -> 26
    34 -> 5
    35 -> 6
    36 -> 15
    36 -> 47
    37 -> 20
    38 -> 21
    38 -> 22
    38 -> 24
    39 -> 5
    40 -> 27
    41 -> 5
    42 -> 29
    43 -> 29
    44 -> 30
    45 -> 31
    46 -> 36
    46 -> 37
    47 -> 37
    32[color="0.01 0.6 0.85", label=call_vars];
    33[color="0.48 0.6 0.85", label=count_raw_bases];
    40[color="0.42 0.6 0.85", label=idsort_aln];
    31[color="0.20 0.6 0.85", label=fusioncatcher];
    36[color="0.33 0.6 0.85", label=clip_trim_sync];
    47[color="0.39 0.6 0.85", label=fastqc_processed];
    37[color="0.03 0.6 0.85", label=rg_stats];
    42[color="0.41 0.6 0.85", label=genome_txt];
    43[color="0.46 0.6 0.85", label=exon_cov_ref];
    45[color="0.37 0.6 0.85", label=merge_fastqs_raw];
    46[color="0.25 0.6 0.85", label=fastqc_raw];
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

2. FusionCatcher, which is used by the fusion detection module. Unfortunately, [this tool is not and will probably never
   be available in Conda](https://github.com/ndaniel/fusioncatcher/issues/20). If you would like to use FusionCatcher
   in the fusion detection steps, we recommend installing version 0.99.5a and then supply the path to the main
   executable using the `fusioncatcher_exe` settings. You may also keep this configuration value empty, though that will
   result in the fusion detection module only using STAR-Fusion and skipping any possible intersection of the results.

3. rose-dt, which is used to detect the FLT3 ITD and KMT2A PTD. This tool may be available in Conda soon, but for now
   it must be compiled and installed directly [from source](https://git.lumc.nl/hem/rose-dt). You can then supply
   the path to the executable, the plot script, and the plot script environment using the `rose_dt_exe` and `plot_itd`
   settings.

4. VEP, which is used to annotate the variant calling results. Hamlet was tested using an older version of VEP (77),
   which could not be made to work with Conda at the time of development. You will need to supply the path to the main
   executable via the `vep_exe` settings.

5. bgzip, which needs to be supplied via the `bgzip_exe` settings. This may be replaced with a proper Conda installation
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
