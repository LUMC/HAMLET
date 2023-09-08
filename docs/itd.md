# itd

The `itd` module is responsible for finding Internal Tandem Duplications in select genes, specifically *FLT3* and *KMT2A*.

## Tools
First, this module uses [bwa]() to align the trimmed reads to a custom reference, which contains the transcript sequence of *FLT3* and *KMT2A*. Next, a custom tool, [rose-dt](https://git.lumc.nl/hem/rose-dt),
is used to detect and visualise Internal Tandem Duplications, using evindence from soft-clipped reads.

## Input
The input for this module is a single pair of FastQ files per sample, specified in a PEP configuration file, as is shown [here](../test/pep/itd.csv).

## Output
The output of this module are a JSON file with an overview of the most important results, as well as a number of other output files:
- For both *FLT3 and *KMT2A*, a .csv file with the detected tandem duplications.
- For both *FLT3* and *KMT2A*, a figure to visualise the detected tandem duplications.

## configuration
The configuration for this module is tailored to the provided reference files, be very careful if you want to modify any of these settings.

| Option                      | description                                   | required |
| --------------------------- | --------------------------------------------- | -------- |
| `fasta`                     | The fasta file, which contains FLT3 and KMT2A | yes      |
| `flt3_name`                 | The name of the FLT3 sequence                 | yes      |
| `flt3_start`                | The start of the FLT3 region to investigate   | yes      |
| `flt3_end`                  | The end of the FLT3 region to investigate     | yes      |
| `kmt2a_name`                | The name of the KMT2A sequence                | yes      |
| `kmt2a_start`               | The start of the KMT2A region to investigate  | yes      |
| `kmt2a_end`                 | The end of the KMT2A region to investigate    | yes      |
