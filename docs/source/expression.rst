expression module
=================

The `expression` module is responsible for determining gene expression levels
from STAR bam and count files. For the highest accuracy the strandedness of the
library preparation can be specified. By default, the module assumes the data
is unstranded.

The expression of a configurable set of housekeeping genes is used to normalize
the expression of the genes of interest.

Tools
-----
This module relies on the STAR count files in combination with a set of
housekeeping genes to normalize gene expression levels.

Input
-----
The minimal input for this module is one BAM file and one STAR count table
specified in a PEP configuration file, as is shown below.

.. csv-table:: Minimal input for the expression module
   :delim: ,
   :file: ../../test/pep/expression.csv

For more accurate results, it is possible to specify the strandedness of your
RNA library prep (`unstranded`, `forward` and `reverse`). If strandedness is
not specified, all samples will be treated as unstranded.

.. csv-table:: Sample configuration with strandedness
   :delim: ,
   :file: ../../test/pep/expression_strandedness.csv

Output
------

* The genes specified under the `report` section of the configuration will be
  included in the HAMLET pdf repor.
* All genes from the `bed` file and the `genes_of_interest` will be included in
  the MultiQC report. Stranded (forward and reverse) and unstranded samples
  will be listed separately, since their values cannot be compared directly.

Configuration
-------------
The following options are available for the `expression` module


Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options

  * - Option
    - Description
    - Required
  * - gtf
    - A GTF file, to look up the ENSG for the housekeeping genes
    - yes
  * - housekeeping
    - A list of genes to use for normalizing the expression
    - yes
  * - bed
    - A BED file with genomic regions (genes) to quantify
    - no
  * - genes_of_interest
    - A list of gene names to quantify (must be present in the gtf)
    - no
  * - report
    - Genes to include in the PDF report, can include names from the bed file
    - no

Example
^^^^^^^
.. code-block:: json

  {
    "gtf": "test/data/reference/hamlet-ref.gtf",
    "housekeeping": [
      "MT-CO2"
    ],
    "bed": "path/to/bed/file.bed",
    "genes_of_interest": [
      "MT-ND3",
      "MT-ND2"
    ],
    "report": [
      "MT-ND3"
    ],
  }
