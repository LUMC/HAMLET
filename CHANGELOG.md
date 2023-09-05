# Changelog

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

## v2.0.0-dev

### Bugfixes
- Fix a bug with inconsistent config setting 'blacklist' in snv-indels 

### Tool changes
- Replace StarFusion and FusionCatcher with Arriba
- Replace VarScan variant caller with VarDict
- Replace GSNAP aligner with STAR

### Tool updates
- Update VEP to 108.2
- Update Picard to 2.27.4
- Update FastQC to 0.11.9
- Update Cutadapt to 4.1

### Speed improvements
- Use multiple threads for Cutadapt, and reduce the compression of output files

### Changes
- Add additional genes of interest
    - SRSF2
    - SF3B1
    - U2AF1
    - BCOR
    - STAG2
    - ZRSR2
    - EZH2
- Filter fusion results based on fusion partners
- Add fusion plots from Arriba
- Add default blacklist with common false-positive variants
- Add support for variant blacklist in VEP hgvsc format
- Add script to generate a configuration file
- Add pipeline to generate reference files
- Add per-module configuration options
- Add support for PEP sample configuration
- Add support for Snakemake 7.8.5
- Remove variants plots
