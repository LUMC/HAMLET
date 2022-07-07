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
[modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).

HAMLET is build to use Singularity to run every Snakemake rule inside its own container. The base execution
environment is for HAMLET defined by an `environment.yml` file.

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
HAMLET requires around 100GB of reference files to run. You can download the data files [here](https://barmsijs.lumc.nl/HAMLET/deps-1.0.0.tar.gz),
the md5sum for the archive is `5541718e8bc17bcd00ec90ff23ebcfae`.
Please contact the author or open an issue if the link is not working.

## Testing
The following commands can be used to test different aspects of HAMLET. If any
of the tests fail, you can inspect the `log.err` and `log.out` files in the run
folder.

Activate the HAMLET conda environment you installed above.
```bash
conda activate HAMLET
```

To test if all dependencies of HAMLET have been installed, use
```bash
pytest --kwd --tag sanity
```


To test if HAMLET can parse the example configurations and find the appropriate output files, use
```bash
pytest --kwd --tag dry-run
```

To test the full behaviour of HAMLET, you can use
```bash
pytest --kwd --tag functional
```

If you want to manually test HAMLET without using pytest-workflow, you can run the following command.

```bash
snakemake --cores 1 --configfile test/data/config/chrM.json --config pepfile=test/pep/chrM.csv --use-singularity
```

# Usage
## Input files
HAMLET requires two separate input files. Firstly, a `json` file that contains
the settings and reference files for the pipeline, you can see an example
[here](test/data/config/chrM.json).

Secondly, HAMLET requires a [Portable Encapsulated
Project](http://pep.databio.org/en/2.1.0/) configuration that specifies the
samples and their associated gzipped, paired-end mRNA-seq files. For simple use
cases, this can be a csv file with one line per read-pair, as can be seen
[here](test/pep/chrM-trio-subsamples.csv).

Any number of samples can be processed in a single execution, and each sample
may have any number of read pairs, and HAMLET will handle those properly.

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
