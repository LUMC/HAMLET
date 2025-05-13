itd module
==========

The ``itd`` module is responsible for finding Internal Tandem Duplications in select genes, specifically *FLT3* and *KMT2A*.

Tools
-----
First, this module uses `bwa <https://github.com/lh3/bwa>`_ to align the trimmed reads to a custom reference, which contains the transcript sequence of *FLT3* and *KMT2A*. Next, a custom tool, `rose-dt <https://git.lumc.nl/hem/rose-dt>`_, is used to detect and visualise Internal Tandem Duplications, using evidence from soft-clipped reads.

Input
-----
The input for this module is a single pair of FastQ files per sample, specified in a PEP configuration file, as is shown below.

.. csv-table:: Example input for the itd module
  :delim: ,
  :file: ../../test/pep/itd.csv

Output
------
The output of this module are a JSON file with an overview of the most important results, as well as a number of other output files:

* For both *FLT3 and *KMT2A*, a .csv file with the detected tandem duplications.
* For both *FLT3* and *KMT2A*, a figure to visualise the detected tandem duplications.

Configuration
-------------
The configuration for this module is tailored to the provided reference files, be very careful if you want to modify any of these settings.
You can automatically generate a configuration for the fusion module using the ``utilities/create-config.py`` script.

.. literalinclude:: ../../test/data/config/itd.json
   :language: json

Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options

  * - Option
    - Description
    - Required
  * - fasta
    - The fasta file containing the trasncript sequence for FLT3 and KMT2A
    - yes
  * - flt3_name
    - The name of the FLT3 sequence in the fasta file
    - yes
  * - flt3_start
    - The start of the FLT3 region to investigate
    - yes
  * - flt3_end
    - The end of the FLT3 region to investigate
    - yes
  * - kmt2a_name
    - The name of the KMT2A sequence
    - yes
  * - kmt2a_start
    - The start of the KMT2A region to investigate
    - yes
  * - kmt2a_end
    - The end of the KMT2A region to investigate
    - yes
