expression
==========

The `expression` module is responsible for determining gene expression levels
from STAR bam and count files. Although the strandedness of the library
preparation is important when determining, the module itself is strand
agnostic. Instead, we take inspiration from STAR and produce output files for
unstranded, forward stranded and reverse stranded libraries, and leave it to
the user to select the relevant output for their samples.

Tools
-----
This module relies on the STAR count files in combination with a set of
housekeeping genes to normalize gene expression levels.

Input
-----
The minimal input for this module is one BAM file and one STAR count table
specified in a PEP configuration file, as is shown below.

.. csv-table:: Example input for the expression module
   :delim: ,
   :file: ../../test/pep/expression.csv

For more accurate results, it is possible to specify the strandedness of your
RNA library prep (`unstranded`, `forward` and `reverse`). If strandedness is
not specified, all samples will be treated as unstranded.

.. csv-table:: Example input for the expression module
   :delim: ,
   :file: ../../test/pep/expression_forward.csv

Output
------

* Three files with the normalized gene expression levels, one for each strandedness.
* A single MultiQC report which contains the same data.

Configuration
-------------
The following options are available for the `expression` module


Example
^^^^^^^
.. literalinclude:: ../../test/data/config/expression.json
   :language: json

Configuration options
^^^^^^^^^^^^^^^^^^^^^
.. list-table:: Configuration options

  * - Option
    - Description
    - Required
  * - housekeeping
    - A list of genes to use for normalizing the expression
    - yes
  * - gtf
    - A GTF file, to look up the ENSG for the housekeeping genes
    - yes
  * - bed
    - A BED file with genomic regions (genes) to quantify
    - yes
