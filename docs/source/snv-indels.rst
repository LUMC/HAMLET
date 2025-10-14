snv-indels module
=================

The ``snv-indels`` module is responsible for aligning the reads to the reference, and calling variants. The bam and count files produced by this module are used in the fusion and gene expression modules.

Tools
-----
This module uses `STAR <https://github.com/alexdobin/STAR>`_ to align the reads
to the reference using two-pass mode. `VarDict
<https://github.com/AstraZeneca-NGS/VarDictJava>`_ is used to call variants,
which are annotated using ``VEP``. Variants are filtered based on the criteria
defined in ``filter_criteria``, and annotated based on ``annotation_criteria``.

The variants annotated by VEP are then filtered based on a number of different criteria:


1. Variant that have a population frequency of more than 1% in the ``gnomADe`` population are excluded.
2. Variants which are specified in ``known_variants`` will be included.
3. Variants that match at least one criteria in ``filter_criteria`` or ``annotation_criteria`` are included.

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
* The filtered VEP output file (``filter_vep``), which contains the final set of filtered and annotated variants.
* The ``counts`` file produced by STAR, which contains the coverage per gene.

Configuration
-------------
You can automatically generate a configuration for the fusion module using the ``utilities/create-config.py`` script.

Example
^^^^^^^
.. literalinclude:: ../../test/data/config/snv-indels.json
   :language: json

Note that the ``vep-cache`` entry is missing for this example file, which means
that VEP will be run with only the ``fasta`` and ``gtf`` files as input. For the best performance, please
specify a ``vep-cache`` folder as well.

Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options
  :widths: 30 80 15
  :header-rows: 1

  * - Option
    - Description
    - Required
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
  * - vep-cache
    - Folder containing the VEP cache
    - no
  * - variant_allele_frequency
    - Minimum variant allele frequency to call a variant
    - no (default=0.05)
  * - min_variant_depth
    - Minimum read depth to call a variant
    - no (default=2)
  * - known_variants
    - File containing known variants and their annotation
    - no

Filter and annotation criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
HAMLET include the ability to specify separate filter criteria for every
transcript, based on the position and the VEP consequence of the variant. The
criteria are used both the filter which variants will be part of the output
(``filter_criteria``), and also annotate the identified variants
(``annotation_criteria``).

The used columns are ``transcript_id``, ``consequence``, ``start``, ``end``
and ``frame``. For annotating variants, the ``annotation`` column is used.
Every column except for ``transcript_id`` can be empty.

.. csv-table:: Example ``filter_criteria`` file, from the HAMLET tests
  :delim: U+0009
  :file: ../../test/data/config/filter_criteria.tsv

.. csv-table:: Example ``annotation_criteria`` file, from the HAMLET tests
  :delim: U+0009
  :file: ../../test/data/config/annotation_criteria.tsv

Known variant annotations
^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to the annotation criteria desribed above, it is also possible to
supply HAMLET with annotations for specific variants via the ``known_variants``
file. Annotations from this file have a higher priority than the annotations
specified in ``annotation_criteria``.

The used columns are ``variant`` and ``annotation``. These columns cannot be 
empty.

.. csv-table:: Example ``known_variants`` file, from the HAMLET tests
  :delim: U+0009
  :file: ../../test/data/config/known_variants.tsv
