*****
Usage
*****


Input files
===========
HAMLET requires two separate input files. Firstly, a ``json`` file that contains
the settings and reference files for the pipeline, which can be generated with
the ``utilities/create-config.py`` script.

Secondly, HAMLET requires a `Portable Encapsulated
Project <http://pep.databio.org/en/2.1.0/>`_ configuration that specifies the
samples and their associated gzipped, paired-end mRNA-seq files. For simple use
cases, this can be a csv file with one line per read-pair, as can be seen below.

.. csv-table:: Example sample specification for HAMLET
   :delim: ,
   :file: ../../test/pep/chrM-trio-subsamples.csv

Any number of samples can be processed in a single execution, and each sample
may have any number of read pairs, and HAMLET will handle those properly.

**Note that spaces in the file paths are supported, but not in sample names**

Execution
=========

If running in a cluster, you may also want to define the resource
configurations in another YAML file. Read more about this type of configuration
on the official `Snakemake documentation
<https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration>`_.
For this file, let's call it ``config-cluster.yml``

Run locally
-----------
To run HAMLET on your local PC, supply the configuration file ``config.json`` and the sample sheet, as well as the flags which indicate to Snakemake that singularity should be used.

.. code:: bash

  $ snakemake -s Snakefile \
      --configfile config.json \
      --config pepfile=sample_sheet.csv \
      --cluster-config config-cluster.yml \
      --rerun-incomplete \
      --use-singularity \
      --singularity-args ' --containall' \
      # ... other flags

.. list-table:: Configuration options
  :widths: 25 80 20
  :header-rows: 1

  * - Option
    - Description
    - Required
  * - --configfile
    - The configuration file for HAMLET
    - yes
  * - --config pepfile
    - A PEP configuration file that contains all samples (can be CSV)
    - yes
  * - --use-singularity
    - Use Singularity or Apptainer to fetch all required dependencies
    - yes
  * - --singularity-args
    - Arguments to pass to singularity. Use --bind to specify which folders on
      your system should be accessible inside the container. This should at
      least be the folders where your samples and reference files are located
    - yes
  * - --cluster-config
    - A cluster configuration file
    - no
  * - --rerun-incomplete
    - Re-run jobs if the output appears incomplete
    - no

Run on Slurm cluster
--------------------


Output files
============
HAMLET will create a separate folder for every sample in the current directory.
Files which are shared across samples will be created once in the current
folder. You can run HAMLET from anywhere, but preferably this is done outside
of the HAMLET folder. This way, the temporary Snakemake files are written
elsewhere and does not pollute the repository.

Inside each sample directory, there will be a PDF report called
``hamlet_report.{sample_name}.pdf`` which contains the overview of the essential
results. The same data is also present in the JSON file called
``{sample_name}.summary.json``.

HAMLET will also run `MultiQC <https://docs.seqera.io/multiqc>`_ and generate a
single html output file which contains quality control metrics for every
sample. This can be used to assess the quality of each individual sample and
find outliers in your sample set.

Grouping results from multiple samples
--------------------------------------

If you analysed multiple samples using HAMLET, you can generate an overview of
multiple samples using the ``utilities/hamlet_table.py`` script, rather than
relying on individual PDF files. This script uses the
``{sample_name}.summary.json`` files which are generated as part of the default
HAMLET output. Simply specify the results you are interested in (``variant``,
``fusion`` or ``itd``) as shown below.

.. code:: bash

  usage: hamlet_table.py [-h] [--itd-gene ITD_GENE]
                         {variant,fusion,itd} json_files [json_files ...]

  positional arguments:
    {variant,fusion,itd}  Table to output
    json_files

  options:
    -h, --help            show this help message and exit
    --itd-gene ITD_GENE
