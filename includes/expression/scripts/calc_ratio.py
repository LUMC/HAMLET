#!/usr/bin/env python

import argparse
import json
import sys

import pandas as pd


GENE_SEP = ":"


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
    exons_df.reset_index(inplace=True)
    exons_df.drop("index", axis=1, inplace=True)

    return exons_df


def main(table, relative_gene, exons, min_ratio):
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-ratio", type=float, default=0.1)
    parser.add_argument("table", type=str)
    parser.add_argument("relative_gene", type=str)
    parser.add_argument("exons", type=str, nargs='+')

    args = parser.parse_args()
    main(args.table, args.relative_gene, args.exons, args.min_ratio)
