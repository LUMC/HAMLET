************
Installation
************

Download HAMLET
===============
You can download and unpack the latest HAMLET release from the `HAMLET release page <https://github.com/LUMC/HAMLET/releases>`_, or clone the github repository directly:

.. code:: bash

   git clone git@github.com:LUMC/HAMLET.git

   # Alternatively, via https
   https://github.com/LUMC/HAMLET.git


Setup the environment
=====================
The dependencies required for running the pipeline are listed in the provided ``environment.yml`` file. To use it, first make sure that you have `Conda <https://docs.conda.io/en/latest/miniconda.html>`_ installed on your system. Then, navigate to the HAMLET folder and set up the Conda virtual environment.

.. code:: bash

   # Navigate to the HAMLET folder, e.g.
   cd ~/Downloads/HAMLET

   # Install the conda environment
   conda env create -f environment.yml

   # Activate the conda environment
   conda activate HAMLET


Additionally, `Singularity <https://docs.sylabs.io/guides/3.0/user-guide/installation.html>`_ version 3 or greater, or `Apptainer <https://apptainer.org/docs/admin/main/installation.html>`_ should be installed on your system.


Download the reference files
============================
HAMLET uses a large number of reference files, in total around 35GB in size. These can be automatically downloaded using the following helper pipeline. Indexing the reference genome for use with STAR will take 60GB RAM and around two hours, so preferably run this step on a workstation or cluster and go on to the next step. The reference files will be placed in the HAMLET-data folder.

.. code:: bash

   # Generate the reference files
   snakemake \
    --snakefile utilities/deps/Snakefile \
    --profile cfg \
    --directory HAMLET-data

  # Create a HAMLET configuration file
  python3 utilities/create-config.py HAMLET-data


Run the test suite
==================
HAMLET comes with a test suite for almost every part of the pipeline, which can be run without the reference files. While the download of the reference files is in progress, you can use the tests to ensure that HAMLET can run on your system. Make sure you are still in the HAMLET conda environment.

.. code:: bash

   # Test if all the HAMLET dependencies have been installed properly
   pytest --kwd --tag sanity

   # Test if HAMLET can parse the example configurations
   pytest --kwd --tag dry-run

   # Run the full HAMLET pipeline on the most complete example data
   #
   # Feel free to investigate the output files produced by the test:
   #
   # - SRR8615409/hamlet_report.SRR8615409.pdf, the PDF report
   # - multiqc_hamlet.html, MultiQC report for all samples
   # - SRR8615409.summary.json, detailed JSON file
   pytest --kwd --tag test-hamlet-chrM

   # Run the HAMLET pipeline for all tests
   pytest --kwd --tag functional

