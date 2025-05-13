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
cases, this can be a ``CSV`` file with one line per read-pair, as can be seen below.

.. csv-table:: Example sample specification for HAMLET
   :delim: ,
   :file: ../../test/pep/chrM-trio-subsamples.csv

Any number of samples can be processed in a single execution, and each sample
may have any number of read pairs, and HAMLET will handle those properly.

**Note that spaces in the file paths are supported, but not in sample names**

Execution
=========

To run the HAMLET pipeline, you need to supply the input files, as well as a
`Snakemake profile
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_,
which configures Snakemake to run the HAMLET pipeline. The example profile,
located in ``HAMLET/cfg/config.v8+.yaml`` is shown below.

Snakemake profile
-----------------
.. literalinclude:: ../../cfg/config.v8+.yaml
   :language: yaml

Please consult the `Snakemake documentation
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#>`_ for an
explanation of all settings.

Make sure to modify the Singularity settings to your specific situation. In
particular, the `--bind` directive determines which parts of the file system
will be visible to HAMLET. In the example, only ``/home`` and ``/tmp`` will be
visible. Make sure that the locations of HAMLET itself, the HAMLET-data as well
as the samples are included here, or HAMLET will not be able to find the
required files.

Since HAMLET includes many tools, the singularity cache will grow to multiple
gigabytes. If you have limited space in your home folder, modify
``singularity-prefix`` to a location with more available space.

The resource requirements will depend on the characteristics of your samples.
The example configuration is based on poly-A captured RNAseq, with up to 200
million reads per sample.

Running HAMLET
--------------
Since all settings can be set in the Snakemake profile, the actual command to
run HAMLET is quite simple.

.. code:: bash

  $ snakemake -s HAMLET/Snakefile \
      --profile HAMLET/cfg \
      --configfile config.json \
      --config pepfile=sample_sheet.csv

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
