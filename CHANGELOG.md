# Changelog

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

## v2.0.0-dev

### Tool changes
- Replace GSNAP aligner with STAR

### Tool updates
- Update VEP to 108.2
- Update Picard to 2.27.4
- Update FastQC to 0.11.9
- Update Cutadapt to 4.1

### Speed improvements
- Use multiple threads for Cutadapt, and reduce the compression of output files

### Changes
- Add support for PEP sample configuration
- Add support for Snakemake 7.8.5
- Remove variants plots
