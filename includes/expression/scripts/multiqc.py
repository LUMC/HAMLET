#!/usr/bin/env python3

import argparse
from collections import defaultdict


def read_expression(fname, strand):
    columns = {"unstranded": 1, "forward": 2, "reverse": 3}
    column = columns[strand]
    with open(fname) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")
            yield spline[0], spline[column]


def main(samples, countfiles, strandedness):
    # Group the samples by strand
    grouped = defaultdict(list)
    for sample, fname, strand in zip(samples, countfiles, strandedness):
        grouped[strand].append((sample, fname))

    for strand in ["unstranded", "forward", "reverse"]:
        samples = [x[0] for x in grouped[strand]]
        fnames = [x[1] for x in grouped[strand]]

        write_multiqc(samples, fnames, strand)


def write_multiqc(samples, countfiles, strand):
    outfile = f"merged_expression_{strand}_mqc.tsv"

    with open(outfile, "wt") as fout:
        # If there are no samples, we don't even write the header
        if not samples:
            return
        # Did we already write the gene header?
        gene_header = None

        # First, we print the multiqc header
        multiqc_header = f"""# plot_type: "table"
# id: mqc_expression_{strand}
# section_name: "{strand.capitalize()}"
# description: "Normalized gene expression assuming <b>{strand} library preparation</b>." """

        print(multiqc_header, file=fout)

        for sample, fname in zip(samples, countfiles):
            # Store the normalized expression for each gene per strand
            expression = dict()
            for values in read_expression(fname, strand):
                gene, level = values
                expression[gene] = level

            # IF this is the first time, write the header
            if not gene_header:
                gene_header = list(expression.keys())
                print("Sample", *gene_header, sep="\t", file=fout)

            print(
                sample, *(expression[gene] for gene in gene_header), sep="\t", file=fout
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--samples", required=True, nargs="+", help="Sample names")
    parser.add_argument(
        "--counts",
        required=True,
        nargs="+",
        help="Files with normalized expression levels",
    )
    parser.add_argument(
        "--strandedness",
        required=True,
        nargs="+",
        choices=["unstranded", "forward", "reverse"],
    )

    args = parser.parse_args()

    if len(args.samples) != len(args.counts) or len(args.samples) != len(
        args.strandedness
    ):
        msg = (
            "--counts, --samples and --strandedness must have the same number of items"
        )
        raise RuntimeError(msg)

    main(args.samples, args.counts, args.strandedness)
