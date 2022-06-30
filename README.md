[![Continuous Integration](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml/badge.svg)](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml)

# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples. Please use the
[public github repository](https://github.com/LUMC/HAMLET) to open issues or pull requests.


Four distinct analysis modules comprise Hamlet, which can be run independently:

  1. `snv-indels`, for small variant detection
  2. `fusion`, for fusion gene detection
  3. `itd`, for tandem duplication detection
  4. `expression`, for expression analysis

There is also a `qc-seq` module which does quality control for input sequence files. Everything is tied together by
a main `Snakefile` using
the [`include`](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes) statement.

Configurations and input file definitions are supplied in a single YAML configuration file. Here you can specify
arbitrary number of samples, each containing an arbitrary number of paired-end files. The merging of these files and
separation of output per-sample are taken care of automatically by Hamlet. Additionally, if you may also supply your own
cluster-specific YAML configuration to determine how much resource such as memory or CPU cores should be allocated to
the jobs.

Hamlet tries to make use of isolated Conda environment per Snakemake rule whenever possible. The base execution
environment is also defined by an `environment.yml` file.

In addition to the raw output files, Hamlet also generates a PDF report containing an overview of the essential results
and a zipped file containing this report and the essential result files.

# Which version should you use?
The version of HAMLET on the 'devel' branch is the latest version, and is probably the version most people are
interested in. The tagged versions on the 'master' branch have been validated more extensively and are used by the LUMC
for clinical purposes. Although we can give **no guarantees** that the output of HAMLET is correct or complete, this is
probably the version you should use (**at your own risk**) if you plan to analyse clinical samples.

# Installation

The dependencies required for running the pipeline are listed in the provided `environment.yml` file. To use it, first
make sure that you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.
Then, set up a Conda virtual environment and then update it:

```bash
# Set up and activate your conda environment.
# Install the dependencies
conda env create -f environment.yml

# Activate the conda environment
conda activate HAMLET
```

Additionally, `singularity` version 3 or greater should be installed on the system.

## Data files
HAMLET requires around 100GB of reference files to run. You can download the data files [here](https://barmsijs.lumc.nl/HAMLET/deps.tar.gz),
the md5sum for the archive is `5ca99cde00955cae44cb63ef3f7d3fd4`.
Please contact the author or open an issue if the link is not working.

## Testing
The following commands can be used to test different aspects of HAMLET. First, activate the HAMLET Conda environment
which was created in the previous step.
```bash
conda activate HAMLET
```

To test if all dependencies of HAMLET have been installed, use
```bash
pytest --tag sanity
```

If any of the tests fail, append `--keep-workflow-wd-on-fail` or `--kwdof` to
the pytest command and inspect the `log.err` and `log.out` files in the run
folder.

To test if HAMLET can parse the Snakemake files and find the appropriate output files, use
```bash
pytest --tag dry-run
```

To test if HAMLET can run the quality control part of the pipeline, using example data, use
```bash
pytest --tag functional
```

To test the full behaviour of HAMLET, you can use
```bash
pytest --tag functional
```

**Important: pytest copies the current directory to /tmp to run the tests.  Therefore, do not place large reference
or sample files inside the HAMLET root folder when running tests, or these will be copied over dozens of times.**

If you want to manually test HAMLET without using pytest-workflow, you can run the following command. Please make sure
you have updated the paths in `test/data/config/test-hamlet-chrM.json` to point to the copy of the HAMLET
reference files.

```bash
snakemake -rp --snakefile Snakefile --configfile test/data/config/test-chrM.json --use-singularity
```

# Usage
## Input files

Hamlet requires gzipped, paired-end mRNA-seq files. Any number of samples can be processed in a single execution, and
each of them may have differing number of pair sets. These files may have arbitrary names, but they *must* be supplied
correctly in the configuration files.

## Execution

After installation of all the required tools, you will need to fill in the required settings and configurations in
several YAML files.

For the runtime settings, use the provided `test/data/config/test-chrM.json` file as template and
update the values as required. Let's call this file `config.json`.

If running in a cluster, you may also want to define the resource configurations in another YAML file. Read more about
this type of configuration on the official [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration). For this
file, let's call it `config-cluster.yml`

### Example command
```bash
$ snakemake -s Snakefile \
    --configfile config.json \
    --cluster-config config-cluster.yml \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args ' --containall --bind /exports/:/exports/' \
    # ... other flags
```

### Explanation for the various flags
| flag | description | required |
| ---- | ----------- | -------- |
| --configfile config.json | The configuration file for the pipeline | Yes |
| --cluster-config | A cluster configuration file, only relevant when you are running HAMLET on a cluster | No |
| --rerun-incomplete | Re-run jobs if the output appears incomplete | No |
| --use-singularity | Use Singularity images to fetch all required dependencies. | Yes |
| --singularity-args | Arguments to pass to singularity. Use --bind to specify which folders on your system should be accessible inside the container. This should at least be the folders where your samples and reference files are located | Yes |

See `test/test_hamlet.yml` for a working example of the flags required to run HAMLET with Singularity.

## Output files

Assuming the output directory is set to `/path/to/output`, Hamlet will create `/path/to/output/{sample_name}` for each
sample present in the config file. Inside the directory, all the essential results are packaged in a zip file called
`hamlet_results.{sample_name}.zip`. This includes a PDF report called `hamlet_report.{sample_name}.pdf` which contains
the overview of the essential results.

## Notes

1. You can run Hamlet from anywhere, but preferably this is done outside of the repository. This way, the temporary
Snakemake files are written elsewhere and does not pollute the repository.

2. You can direct Hamlet to create the output directory anywhere. This is a configuration value that is supplied in the
config file via `output_dir`.

# Citation
If you use HAMLET in your research, please cite the [HAMLET publication](https://www.nature.com/articles/s41375-020-0762-8).

# Common issues
## Using singularity
If you forget the `--use-singularity` flag for Snakemake, you will find that many rules break due to the required tools
not being available on your system.

## Snakemake errors about reserved keywords
If you install Snakemake manually instead of using Conda and the provided `environment.yml` file, you might get errors
about reserved keyword that are used in the Snakefiles. Please use the Snakemake version specified in the
`environment.yml` file.
