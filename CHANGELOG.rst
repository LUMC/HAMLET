#########
Changelog
#########

.. Newest changes should be on top.

..  This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

**********
v2.3.1-dev
**********
* **Breaking change**: Update the json output format
* **Breaking change**: Update Snakemake to version 8
* Fix a bug with the Java runtime environment for Picard
* Fix a bug with caching of report assets introduced in snakemake 8
* Fix a bug with structural variants missing effect prediction
* Add an option to specify `variant_allele_frequency` for snv-indels
* Speed up VarDict by running with 8 threads
* Speed up VarDict by sorting the call_regions.bed file
* Run VEP with either vep_cache or just the gtf/fasta files as database
* Update picard from 2.27.4 to 3.3.0
* Update cutadapt from 4.6 to 5.0
* Update sequali from 0.9.1 to 0.12.0
* Update star from 2.7.10b to 2.7.11b
* Update MultiQC to 1.27.1

**********
v2.2.1
**********

Breaking changes
================
* The `bed_variant_call_regions` option has been removed, variants are now
  called for all genes present in the `gtf` file.
* Add graphviz/`dot` as a dependency (developer only).
* Please create a new HAMLET configuration file with `create-config.py` script.
* To use the latest hotspot regions and artifact blacklist, please recreate the
  HAMLET reference data.

Novel module
============
* Add novel module, **expression**, which analyzes gene expression.
    * Add optional input `strandedness` to the sample configuration.
    * Add json output file for the expression module.

Bugfixes
========
* Fix a rare bug where different modules use the same MultiQC file list.
* Fix a bug with filtering VEP records that contain multiple population.
  frequency records for a single variant.

Updates
=======
* Add ability to generate configurations for each module using the
  `utilities/create-config.py` script.
* Update the hotspot regions reference file.
* Update the blacklist of known artifacts.
* Remove various superfluous plots from the MultiQC report.

**********
v2.1.3
**********
* Add `pysam` to the conda environment (developer only)
* Add exon number to variant table
* Add chromosomes to fusion table

**********
v2.1.2
**********
* Set the maximum population frequency to 1%. This was accidentally set to 5% in v2.1.1

**********
v2.1.1
**********

Bugfixes
========
* Fix a bug where VEP removed rare variants

**********
v2.1.0
**********

Breaking changes
================
* Remove the JSON output for the qc-seq module (this has been replaced by a
  MultiQC report)
* Add sample name to STAR counts table

Changes
=======
* Automatically remove _STAR temporary folders
* Modified PDF formatting
    * Change cover image
    * Add bookmarks under chapter variant
    * Sort the genes of interest alphabetically
    * Remove the "Sequencing Results" section from the report (this has been replaced by a MultiQC report)
* Replace FastQC with Sequali

Bugfixes
========
* Fix a bug where the trimmed FastQ files are not removed when no longer needed

Updates
=======
* Update Cutadapt to 4.6
* Update MultiQC to 1.22
* Update snakefmt to 0.10.0 (developer only)
* Update black to 24.3.0 (developer only)


**********
v2.0.5
**********
* Increase space for the HGVS description in "Results Overview" table

**********
v2.0.4
**********
* Automatically check the release tag is set correctly

**********
v2.0.3
**********

Bugfixes
========
* Fix a bug where long HGVS descriptions make the "Results Overview" table overflow the page

**********
v2.0.2
**********
* Include the sample name in the final BAM file

**********
v2.0.1
**********

Bugfixes
========
* Update version number in HAMLET report

**********
v2.0.0
**********

Bugfixes
========
* Fix a bug with inconsistent config setting 'blacklist' in snv*indels
* Fix a bug where unmapped reads are not included in STAR output file

Tool changes
============
* Replace StarFusion and FusionCatcher with Arriba
* Replace VarScan variant caller with VarDict
* Replace GSNAP aligner with STAR

Tool updates
============
* Update VEP to 108.2
* Update Picard to 2.27.4
* Update FastQC to 0.11.9
* Update Cutadapt to 4.1

Speed improvements
==================
* Use multiple threads for Cutadapt, and reduce the compression of output files

Changes
=======
* Remove run name from the report
* Deprecate option `fusion*partners`, in favour of `report_genes`, which points
  to a list of fusion genes to report
* Show allele frequency as a percentage in the pfd report
* Add additional genes of interest
    - SRSF2
    - SF3B1
    - U2AF1
    - BCOR
    - STAG2
    - ZRSR2
    - EZH2
* Filter fusion results based on fusion partners
* Add fusion plots from Arriba
* Add default blacklist with common false*positive variants
* Add support for variant blacklist in VEP hgvsc format
* Add script to generate a configuration file
* Add pipeline to generate reference files
* Add per*module configuration options
* Add support for PEP sample configuration
* Add support for Snakemake 7.8.5
* Remove variants plots
* Use MANE select transcript for all genes
