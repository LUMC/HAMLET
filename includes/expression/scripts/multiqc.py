#!/usr/bin/env python3

import argparse
import json


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

def write_cell_types(sample_json):
    outfile="merged_expression_cell_types_mqc.tsv"
    url="https://github.com/eonurk/seAMLess>seAMLess"

    with open(outfile, "wt") as fout:
        multiqc_header = f"""# plot_type: "table"
# id: mqc_expression_cell_types
# section_name: Cell type composition
# description: Estimated cell type composition based on <a href={url}>seAMLess</a>."""

        print(multiqc_header, file=fout)
        header = None

        for fname in sample_json:
            with open(fname) as fin:
                js = json.load(fin)

            # First time
            if header is None:
                header = list(js["expression"]["cell-types"].keys())
                header.insert(0, "Sample")
                print(*header, sep='\t', file=fout)

            # Get the cell type data
            data = js["expression"]["cell-types"]
            # Add the sample
            sample = js["expression"]["metadata"]["sample_name"]
            data["Sample"] = sample

            # Print
            print(*(data[field] for field in header), sep="\t", file=fout)

def main(samples, countfiles, strandedness, sample_json):
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
    write_cell_types(sample_json)


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
    parser.add_argument(
        "--sample-json",
        required=True,
        nargs="+",
        help="Per sample expression json files"
    )

    args = parser.parse_args()

    if len(args.samples) != len(args.counts) or len(args.samples) != len(
        args.strandedness
    ):
        msg = (
            "--counts, --samples and --strandedness must have the same number of items"
        )
        raise RuntimeError(msg)

    main(args.samples, args.counts, args.strandedness, args.sample_json)
