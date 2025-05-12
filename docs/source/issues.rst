*************
Common issues
*************

Using singularity
=================
If you forget the ``--use-singularity`` flag for Snakemake, you will find that
many rules break due to the required tools not being available on your system.

Snakemake errors about reserved keywords
========================================
If you install Snakemake manually instead of using Conda and the provided
``environment.yml`` file, you might get errors about reserved keyword that are
used in the Snakefiles. Please use the Snakemake version specified in the
``environment.yml`` file.
