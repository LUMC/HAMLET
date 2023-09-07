# qc-seq

The `qc-seq` module is responsible for removing adapter sequences and low
quality reads, and generating read-level statistics. It also merges the FastQ
files per sample, so they can be used by the other modules. Every set of FastQ
files can be analysed in parallel.

## Tools
This module uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove adapter sequences and low quality bases.
[FastQC](https://cutadapt.readthedocs.io/en/stable/) is used to generate detailed quality statistics.

## Input
The input for this module is one or more pairs of FastQ files per sample, specified in a PEP configuration file, as is shown [here](../test/pep/chrM-trio-subsamples.csv).

## Output
The output of this module are one set of merged FastQ files per sample, as well as a JSON file with statistics.

## configuration
The only configurable option for this module is adapter sequences for
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove.

| Option            | description                  | required |
| ---------------   | ---------------------------- | -------- |
| `forward_adapter` | The forward adapter sequence | yes      |
| `reverse_adapter` | The reverse adapter sequence | yes      |
