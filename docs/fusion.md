# fusion

The `fusion` module uses [Arriba](https://github.com/suhrig/arriba) to call fusion events.

## Tools
This module uses the bam file from [STAR](https://github.com/alexdobin/STAR) to
call fusion events.

The fusion events are filtered based on the `blacklist` from Arriba itself. Also, only fusions where at least one of the involved genes is in `report_genes` will be included in the final output.

For each fusion event that remains after filtering, we also generate a figure using the `draw_fusions.R` script provided by Arriba.

## Input
The input for this module is a single bam file, generated by STAR per sample, specified in a PEP configuration file, as is shown [here](../test/pep/chrM-bam.csv).

## Output
The output of this module are a JSON file with an overview of the most important results, as well as a number of other output files:
- The final Arriba output file, after filtering.
- One figure per fusion event

## configuration
| Option                      | description                             | required |
| --------------------------- | --------------------------------------- | -------- |
| `genome_fasta`              | Reference genome, in FASTA format       | yes      |
| `gtf`                       | GTF file with transcript information    | yes      |
| `blacklist`                 | File of blacklisted variants            | yes      |      
| `known_fusions`             | A file of known fusion events           | yes      |
| `report_genes`              | A file of genes to report fusions for   | yes      |
| `cytobands`                 | A file with cytoband information        | yes      |
| `protein_domains`           | A file with protein domains             | yes      |
