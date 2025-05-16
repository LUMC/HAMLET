expression module
=================

The ``expression`` module is responsible for determining gene expression levels
from STAR bam and count files. For the highest accuracy the strandedness of the
library preparation can be specified. By default, the module assumes the data
is unstranded.

The expression of a configurable set of housekeeping genes is used to normalize
the expression of the genes of interest.

Tools
-----
This module relies on the STAR count files in combination with a set of
housekeeping genes to normalize gene expression levels.

`seAMLess <https://github.com/eonurk/seAMLess>`_ is used to determine the cell
type composition of the sample, as well as the resistance to Venetoclax. Please
see `Karakaslar et al. <https://www.nature.com/articles/s41698-024-00596-9>`_
for details.

Input
-----
The minimal input for this module is one BAM file and one STAR count table
specified in a PEP configuration file, as is shown below.

.. csv-table:: Minimal input for the expression module
   :delim: ,
   :widths: 25 53 77
   :header-rows: 1
   :file: ../../test/pep/expression.csv

For more accurate results, it is possible to specify the strandedness of your
RNA library prep (``unstranded``, ``forward`` and ``reverse``). If strandedness is
not specified, all samples will be treated as unstranded.

.. csv-table:: Sample configuration with strandedness
   :delim: ,
   :widths: 25 24 52 50
   :header-rows: 1
   :file: ../../test/pep/expression_strandedness.csv

Output
------

* The genes specified under the ``report`` section of the configuration will be
  included in the HAMLET pdf report.
* All genes from the ``bed`` file and the ``genes_of_interest`` will be included in
  the MultiQC report. Stranded (forward and reverse) and unstranded samples
  will be listed separately, since their values cannot be compared directly.
* The cell type composition for each sample wil be included in the HAMLET pdf
  report. In addition, and overview of the composition of all samples is
  included in the MultiQC report.

Note that the overexpression and cell type composition results can be exported
to csv directly from the MultiQC html file.

Configuration
-------------
The following options are available for the ``expression`` module


Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options
  :widths: 25 80 20
  :header-rows: 1

  * - Option
    - Description
    - Required
  * - seamless_ref
    - Path to seAMLess reference expression csv file (genes x samples)
    - yes
  * - seamless_meta
    - Path to seAMLess metadata csv file (samples x features)
    - yes
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
    "seamless_ref": "test/data/seAMLess/seAMLess_expression_reference.csv",
    "seamless_meta": "test/data/seAMLess/seAMLess_expression_metadata.csv",
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
