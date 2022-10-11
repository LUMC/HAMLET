#!/usr/bin/env python

import json
import sys

import click
import pandas as pd


GENE_SEP = ":"
COORD_SEP = "-"


def process_table(relative_gene, exons, df, min_threshold):
    res = {}
    _, sample = df.columns.values
    df.columns = ("exon", "count")
    gene_df = df[df["exon"].str.startswith(relative_gene + GENE_SEP)]
    gene_total_bases = gene_df["count"].sum()

    exons_df = df[df.exon.apply(lambda v: v in exons)]
    exons_df.insert(0, "sample_name", sample)
    exons_df = exons_df.assign(ratio=exons_df["count"] / float(gene_total_bases))
    exons_df = exons_df.assign(above_threshold=exons_df["ratio"] >= min_threshold)
    exons_df = exons_df.assign(divisor_gene=relative_gene)
    exons_df = exons_df.assign(divisor_exp=gene_total_bases)
    #exons_df["above_threshold"] = exons_df["above_threshold"].map({False: "no", True: "yes"})
    exons_df.reset_index(inplace=True)
    exons_df.drop("index", axis=1, inplace=True)

    return exons_df


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("table", type=str)
@click.argument("relative_gene", type=str)
@click.argument("exons", type=str, nargs=-1)
@click.option("-r", "--min-ratio", type=float,
              default=0.1)
@click.option("--sample-id", type=str,
              help="Sample ID.")
def main(table, relative_gene, exons, min_ratio, sample_id):
    if table == "-":
        fh = sys.stdin
    else:
        fh = open(table, "r")

    df = pd.read_table(fh)
    edf = process_table(relative_gene, exons, df, min_ratio)
    result = edf.to_json(index=False, orient="table")

    # Get the data
    ratios = json.loads(result)["data"]

    # Clean up the data
    for r in ratios:
        # Remove sample name
        r.pop("sample_name")
        # If there is no ratio, set it to zero
        if r["ratio"] is None:
            r["ratio"] = 0

    # Write to stdout
    json.dump({"expression": ratios}, sys.stdout, indent=2)

if __name__ == "__main__":
    main()
