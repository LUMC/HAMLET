snv-indels module
=================

The `snv-indels` module is responsible for aligning the reads to the reference, and calling variants. The bam and count files produced by this module are used in the fusion and gene expression modules.

Tools
-----
This module uses `STAR <https://github.com/alexdobin/STAR>`_ to align the reads to the reference using twopass mode. `VarDict <https://github.com/AstraZeneca-NGS/VarDictJava>`_ is used to call variants, which are annotated using `VEP <https://www.ensembl.org/info/docs/tools/vep/index.html>`_.
Variants are filtered based on the criteria defined in `filter_criteria`, and annotated based on `annotation_criteria`.

The variants annotated by VEP are then filtered based on a number of different criteria:

1. Variants that are present on the `blacklist` are excluded.
2. Only variants that are present on one of the specified transcripts in `ref_id_mapping` are included.
3. Only variants that match one of the consequences defined in `vep_include_consequence` are included.
4. Variant that have a population frequency of more than 1% in the `gnomADe` population are excluded.

Picard is used to generate various alignment statistics.

Input
-----
The input for this module is a single pair of FastQ files per sample, specified in a PEP configuration file, as is shown below.

.. csv-table:: Example input for the snv-indels module
  :delim: ,
  :file: ../../test/pep/targetted.csv

Output
------
The output of this module are a JSON file with an overview of the most important results, as well as a number of other output files:

* A .bam and .bai per sample, which contain the aligned reads.
* A VEP output file (`vep_high`), which contains the final set of filtered variants.
* A VEP output file (`vep_target`), which contains the variants on the transcripts of interest. These variants have not been filtered on `vep_include_consequence` terms.

Configuration
-------------
You can automatically generate a configuration for the fusion module using the `utilities/create-config.py` script.

Example
^^^^^^^
.. literalinclude:: ../../test/data/config/snv-indels.json
   :language: json

Note that the `vep-cache` entry is missing for this example file, which means
that the online API of VEP will be used. For the best performance, please
specify a `vep-cache` folder as well.

Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options
  :widths: 30 80 15
  :header-rows: 1

  * - Option
    - Description
    - Required
  * - forward_adapter
    - The forward adapter sequence
    - yes
  * - reverse_adapter
    - The reverse adapter sequence
    - yes
  * - genome_fasta
    - Reference genome, in FASTA format
    - yes
  * - genome_fai
    - .fai index file for the reference fasta
    - yes
  * - genome_dict
    - .dict index file for the reference fasta
    - yes
  * - star_index
    - STAR index database
    - yes
  * - ref_id_mapping
    - File of transcripts of interest
    - yes
  * - filter_criteria
    - Criteria file to filter variants
    - yes
  * - annotation_criteria
    - Criteria file to annotate variants
    - yes
  * - rrna_refflat
    - File of rRNA transcripts
    - yes
  * - gtf
    - GTF file with transcripts, used by STAR
    - yes
  * - annotation_refflat
    - File used to determine exon coverage
    - yes
  * - blacklist
    - File of blacklisted variants
    - yes
  * - vep-cache
    - Folder containing the VEP cache
    - no
  * - vep_include_consequence
    - List of VEP consequences to report
    - yes
  * - variant_allele_frequency
    - Minimum variant allele frequency in the sample to call a variant

      (default=0.05)
    - no 
