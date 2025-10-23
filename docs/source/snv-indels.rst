snv-indels module
=================

The ``snv-indels`` module is responsible for aligning the reads to the reference, and calling variants. The bam and count files produced by this module are used in the fusion and gene expression modules.

Tools
-----
This module uses `STAR <https://github.com/alexdobin/STAR>`_ to align the reads
to the reference using two-pass mode. `VarDict
<https://github.com/AstraZeneca-NGS/VarDictJava>`_ is used to call variants,
which are annotated using ``VEP``. Variants are filtered based on the criteria
defined in ``inclusion_criteria``, and annotated based on ``annotation_criteria``.

The variants annotated by VEP are then filtered and annotated in the following order:

1. Variants that have a population frequency of more than 1% in the ``gnomADe``
   population are removed.
2. Only variant which match the ``inclusion_criteria`` will be included.
3. If a variant is present in ``known_variants``, the annotation from that file
   will be added to the variant.
4. If a variant is not in ``known_variants``, it will be checked against the
   criteria in ``annotation_criteria``. The annotation from the first matching
   definition will be added to the variant.

Picard is used to generate various alignment statistics.

Variant annotations
-------------------
By default, HAMLET comes with variant filters and annotations which are tuned towards diagnosing AML. When using the default variant filters and annotations, HAMLET uses the following definitions:

.. list-table:: Configuration options
  :widths: 30 70
  :header-rows: 1

  * - Annotation
    - Definition
  * - Known pathogenic
    - This variant is known to be associated with AML
  * - Pathogenic
    - All evidence point to this variant being pathogenic for AML
  * - Likely pathogenic
    - This variant should be considered pathogenic, unless there is evidence to
      the contrary (*e.g.* it is a known benign variant)
  * - Possible pathogenic
    - This variant should **not** be considered pathogenic, unless there is
      additional evidence (*e.g.* it is a known pathogenic variant)
  * - Discard
    - This variant should not be considere pathogenic
  * - Artifact
    - This variant is most likely an artifact produced by the pipeline, *i.e.*
      the variant is not truly present in the sample

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
* Various quality control metrics produced by MultiQC.

Configuration
-------------
You can automatically generate a configuration for the fusion module using the ``utilities/create-config.py`` script.

Example
^^^^^^^

.. code:: bash

   $ python3 utilities/create-config.py --module snv-indels HAMLET-data

    {
     "annotation_criteria": "HAMLET-data/annotation_criteria.tsv",
     "annotation_refflat": "HAMLET-data/ucsc_gencode.refFlat",
     "inclusion_criteria": "HAMLET-data/filter_criteria.tsv",
     "genome_dict": "HAMLET-data/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
     "genome_fai": "HAMLET-data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
     "genome_fasta": "HAMLET-data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
     "gtf": "HAMLET-data/Homo_sapiens.GRCh38.115.chr.gtf",
     "known_variants": "HAMLET-data/known_variants.tsv",
     "min_variant_depth": 2,
     "rrna_refflat": "HAMLET-data/ucsc_rrna.refFlat",
     "star_index": "HAMLET-data/star-index",
     "variant_allele_frequency": 0.05,
     "vep_cache": "HAMLET-data"
    }


Note that although the ``vep-cache`` entry is optional, for the best
performance, please specify a ``vep-cache`` folder as well.

Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options
  :widths: 30 70 25
  :header-rows: 1

  * - Option
    - Description
    - Required
  * - annotation_criteria
    - Criteria file to annotate variants
    - yes
  * - annotation_refflat
    - File used to determine exon coverage
    - yes
  * - inclusion_criteria
    - Criteria file to filter variants
    - yes
  * - genome_dict
    - .dict index file for the reference fasta
    - yes
  * - genome_fai
    - .fai index file for the reference fasta
    - yes
  * - genome_fasta
    - Reference genome, in FASTA format
    - yes
  * - gtf
    - GTF file with transcripts, used by STAR
    - yes
  * - known_variants
    - File containing known variants and their annotation
    - no
  * - min_variant_depth
    - Minimum read depth to call a variant
    - no (default=2)
  * - rrna_refflat
    - File of rRNA transcripts
    - yes
  * - star_index
    - STAR index database
    - yes
  * - variant_allele_frequency
    - Minimum variant allele frequency to call a variant
    - no (default=0.05)
  * - vep-cache
    - Folder containing the VEP cache
    - no

Filter and annotation criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
HAMLET include the ability to specify separate filter criteria for every
transcript, based on the position and the VEP consequence of the variant. The
criteria are used both to filter which variants will be part of the output
(``inclusion_criteria``), and also annotate the identified variants
(``annotation_criteria``).

The used columns are ``transcript_id``, ``consequence``, ``start``, ``end``
and ``frame``. For annotating variants, the ``annotation`` column is used.
Every column except for ``transcript_id`` can be empty.

.. csv-table:: Example ``inclusion_criteria`` file, from the HAMLET tests
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
