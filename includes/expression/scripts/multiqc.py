#!/usr/bin/env python3

import argparse


def read_expression(fname, strand):
    """Read the expression data from the correct column for the specified strand"""
    columns = {"unstranded": 1, "forward": 2, "reverse": 3}
    column = columns[strand]

    expression = dict()
    with open(fname) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")
            gene = spline[0]
            expression[gene] = spline[column]
    return expression


def read_data(samples, countfiles, strandedness):
    """Read the data for every sample"""
    data = dict()
    for sample, fname, strand in zip(samples, countfiles, strandedness):
        data[sample] = read_expression(fname, strand)
    return data


def main(samples, countfiles, strandedness):
    # Read all the data
    data = read_data(samples, countfiles, strandedness)

    # Group the samples by strandedness
    unstranded = dict()
    stranded = dict()
    for sample, strand in zip(samples, strandedness):
        if strand == "unstranded":
            unstranded[sample] = data[sample]
        elif strand =="forward" or strand == "reverse":
            stranded[sample] = data[sample]
        else:
            raise RuntimeError

    write_multiqc(unstranded, "unstranded")
    write_multiqc(stranded, "stranded")


def write_multiqc(data, strand):
    outfile = f"merged_expression_{strand}_mqc.tsv"

    with open(outfile, "wt") as fout:
        # If there are no data for this strand, do nothing
        if not data:
            return
        # Did we already write the gene header?
        gene_header = None

        # First, we print the multiqc header
        multiqc_header = f"""# plot_type: "table"
# id: mqc_expression_{strand}
# section_name: "{strand.capitalize()}"
# description: "Normalized gene expression assuming <b>{strand} library preparation</b>." """

        print(multiqc_header, file=fout)

        for sample, genes in data.items():
            # IF this is the first time, write the header
            if not gene_header:
                gene_header = list(genes.keys())
                print("Sample", *gene_header, sep="\t", file=fout)

            print(
                sample, *(genes[gene] for gene in gene_header), sep="\t", file=fout
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
