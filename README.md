[![Continuous Integration](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml/badge.svg)](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pep.databio.org)
![GitHub release](https://img.shields.io/github/v/release/LUMC/HAMLET)
![Commits since latest release](https://img.shields.io/github/commits-since/LUMC/HAMLET/latest)

# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples. Please use the
[public github repository](https://github.com/LUMC/HAMLET) to open issues or pull requests.


Four distinct analysis modules comprise Hamlet, which can be run independently and have their own documentation:

  1. [qc-seq](docs/qc-seq.md), for adapter trimming and quality control
  2. [snv-indels](docs/snv-indels.md), for small variant detection
  3. [fusion](docs/fusion.md), for fusion gene detection
  4. [itd](docs/itd.md), for tandem duplication detection

Everything is tied together by a main `Snakefile` using
[modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).

HAMLET is build to use Singularity to run every Snakemake rule inside its own container. The base execution
environment for HAMLET defined by an `environment.yml` file.

In addition to the raw output files, Hamlet also generates a PDF report containing an overview of the essential results
and a JSON file containing the underlying data that are shown in the report.

# Installation
The dependencies required for running the pipeline are listed in the provided `environment.yml` file. To use it, first
make sure that you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.
Then, set up a Conda virtual environment and activate it:

```bash
# Set up and activate your conda environment.
# Install the dependencies
conda env create -f environment.yml

# Activate the conda environment
conda activate HAMLET
```

Additionally, [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) version 3 or greater should be installed on the system.

## Data files
Automatically generate the required reference files for the HAMLET pipeline
in the `HAMLET-data` folder with
```bash
snakemake \
    --snakefile utilities/deps/Snakefile \
    --use-singularity \
    --singularity-args '--cleanenv --bind /tmp' \
    --directory HAMLET-data
```

Next, you can automatically generate a configuration file with the following helper script
```bash
python3 utilities/create-config.py HAMLET-data
```

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

# Usage
## Input files
HAMLET requires two separate input files. Firstly, a `json` file that contains
the settings and reference files for the pipeline, see above.

Secondly, HAMLET requires a [Portable Encapsulated
Project](http://pep.databio.org/en/2.1.0/) configuration that specifies the
samples and their associated gzipped, paired-end mRNA-seq files. For simple use
cases, this can be a csv file with one line per read-pair, as can be seen
[here](test/pep/chrM-trio-subsamples.csv).

Any number of samples can be processed in a single execution, and each sample
may have any number of read pairs, and HAMLET will handle those properly.

## Execution
If running in a cluster, you may also want to define the resource configurations in another YAML file. Read more about
this type of configuration on the official [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration). For this
file, let's call it `config-cluster.yml`

### Example command
```bash
$ snakemake -s Snakefile \
    --configfile config.json \
    --config pepfile=sample_sheet.csv \
    --cluster-config config-cluster.yml \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args ' --containall' \
    # ... other flags
```

### Explanation for the various flags
| flag | description | required |
| ---- | ----------- | -------- |
| --configfile config.json | The configuration file for the pipeline | Yes |
| --cluster-config | A cluster configuration file, only relevant when you are running HAMLET on a cluster | No |
| --config pepfile=sample_sheet.csv | A PEP configuration file that contains all samples, can be CSV | Yes |
| --rerun-incomplete | Re-run jobs if the output appears incomplete | No |
| --use-singularity | Use Singularity images to fetch all required dependencies. | Yes |
| --singularity-args | Arguments to pass to singularity. Use --bind to specify which folders on your system should be accessible inside the container. This should at least be the folders where your samples and reference files are located | Yes |

## Output files

Assuming the output directory is set to `/path/to/output`, Hamlet will create
`/path/to/output/{sample_name}` for each sample present in the config file.
Inside the directory, there will be a PDF report called
`hamlet_report.{sample_name}.pdf` which contains the overview of the essential
results. The same data is also present in the JSON file called `{sample_name}.summary.json`.

## Grouping results from multiple samples
If you analysed multiple samples using HAMLET, you can generate an overview of
multiple samples using the `utilities/hamlet_table.py` script, rather than
reading many PDF files. This script uses the `{sample_name}.summary.json` files
which are generated as part of the default HAMLET output.

```bash
$ python3 utilities/hamlet_table.py --help
usage: hamlet_table.py [-h] [--itd-gene ITD_GENE] {variant,fusion,itd} json_files [json_files ...]

positional arguments:
  {variant,fusion,itd}  Table to output
  json_files

options:
  -h, --help            show this help message and exit
  --itd-gene ITD_GENE
```

## Notes

1. You can run Hamlet from anywhere, but preferably this is done outside of the repository. This way, the temporary
Snakemake files are written elsewhere and does not pollute the repository.

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
