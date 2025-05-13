==========
Background
==========

HAMLET is a pipeline for analysis of human acute myeloid leukemia RNA-seq samples. Please use the `public github repository <https://github.com/LUMC/HAMLET>`_ to open issues or pull requests. It consists of a number of distinct modules, which can be run independently and have their own documentation. Everything is tied together by the main ``Snakefile`` using `modules <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules>`_.

HAMLET is build to use Singularity or Apptainer to run every Snakemake rule inside its own container. The base execution environment for HAMLET is defined by an ``environment.yml`` file.

In addition to the raw output files for each sample, HAMLET also generates a PDF report containing an overview of the essential results. Additional information is available from the summary json file.

Citation
========
If you use HAMLET in your research, please cite the `HAMLET publication <https://www.nature.com/articles/s41375-020-0762-8>`_.
