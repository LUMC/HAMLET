#!/usr/bin/env python3

import argparse


def read_expression(fname, strand):
    columns = {"unstranded": 1, "forward": 2, "reverse": 3}
    column = columns[strand]
    with open(fname) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")
            yield spline[0], spline[column]


def main(samples, countfiles, strand):
    # Did we already write the gene header?
    gene_header = None

    # First, we print the multiqc header
    multiqc_header = f"""# plot_type: "table"
# id: mqc_expression_{strand}
# section_name: "{strand.capitalize()}"
# description: "Normalized gene expression assuming <b>{strand} library preparation</b>." """

    print(multiqc_header)

    for sample, fname in zip(samples, countfiles):
        # Store the normalized expression for each gene per strand
        expression = dict()
        for values in read_expression(fname, args.strand):
            gene, level = values
            expression[gene] = level

        # IF this is the first time, write the header
        if not gene_header:
            gene_header = list(expression.keys())
            print("Sample", *gene_header, sep="\t")

        print(sample, *(expression[gene] for gene in gene_header), sep="\t")


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
        "--strand", required=True, choices=["unstranded", "forward", "reverse"]
    )

    args = parser.parse_args()

    if len(args.samples) != len(args.counts):
        msg = "--counts and --samples must have the same number of items"
        raise RuntimeError(msg)

    main(args.samples, args.counts, args.strand)
