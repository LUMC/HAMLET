[![Continuous Integration](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml/badge.svg)](https://github.com/LUMC/HAMLET/actions/workflows/ci.yml)

# Hamlet

Hamlet is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples. Please use the
[public github repository](https://github.com/LUMC/HAMLET) to open issues or pull requests.


Four distinct analysis modules comprise Hamlet, which can be run independently:

  1. `snv-indels`, for small variant detection
  2. `fusion`, for fusion gene detection
  3. `itd`, for tandem duplication detection

There is also a `qc-seq` module which does quality control for input sequence files. Everything is tied together by
a main `Snakefile` using
[modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).

HAMLET is build to use Singularity to run every Snakemake rule inside its own container. The base execution
environment for HAMLET defined by an `environment.yml` file.

In addition to the raw output files, Hamlet also generates a PDF report containing an overview of the essential results
and a zipped file containing this report and the essential result files.

# Which version should you use?
The version of HAMLET on the 'devel' branch is the latest version, and is probably the version most people are
interested in.

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
python3 python3 utilities/create-config.py HAMLET-data
```

For the old **master** branch of HAMLET, you can download the data files
[here](https://barmsijs.lumc.nl/HAMLET/deps-1.0.0.tar.gz), the md5sum for the
archive is `5541718e8bc17bcd00ec90ff23ebcfae`. Please contact the author or
open an issue if the link is not working.

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
the settings and reference files for the pipeline. You can use the example
below, just update the path to the `deps-1.0.0` folder.
```json
{
  "kmt2a_fasta": "deps-1.0.0/data/bwa-kmt2a-index/kmt2a-213.fa",
  "genome_star_fusion_lib": "deps-1.0.0/data/GRCh38_gencode_v31_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir",
  "flt3_start": 1787,
  "gtf": "deps-1.0.0/data/ucsc_gencode.gtf",
  "vcf_gonl": "deps-1.0.0/data/gonl.grch38.sorted.filtered.snps_indels.r5.vcf.gz",
  "kmt2a_start": 406,
  "annotation_refflat": "deps-1.0.0/data/ucsc_gencode.refFlat",
  "kmt2a_end": 4769,
  "kmt2a_name": "KMT2A-213",
  "bed_variant_hotspots": "deps-1.0.0/data/hotspots_genome.bed",
  "exon_names": [
    "MECOM:169146722-169147734"
  ],
  "exon_min_ratio": 0.1,
  "flt3_end": 2024,
  "genome_fasta": "deps-1.0.0/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
  "fusioncatcher_data": "deps-1.0.0/fusioncatcher/human_v98",
  "flt3_name": "FLT3-001",
  "cache_vep": "deps-1.0.0/vep/ensembl-vep/cache_dir/",
  "ref_id_mapping": "deps-1.0.0/data/id_mappings.tsv",
  "flt3_fasta": "deps-1.0.0/data/bwa-flt3-index/flt3-001.fa",
  "rrna_refflat": "deps-1.0.0/data/ucsc_rrna.refFlat",
  "genome_gmap_index": "deps-1.0.0/data/gmap-genome-index/GCA_000001405.15_GRCh38_no_alt_analysis_set",
  "run_name": "104372-045",
  "vcf_1kg": "deps-1.0.0/data/1KG.phase3.GRCh38.with_chr.vcf.gz",
  "transcripts_bed": "deps-1.0.0/data/transcripts.bed"
}
```

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
    --singularity-args ' --containall --bind /exports/:/exports/' \
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
