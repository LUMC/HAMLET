#!/usr/bin/env python3

import argparse
import statistics
from gtf import gene_id_name


def get_coverage(fname):
    cov = dict()
    with open(fname) as fin:
        for line in fin:
            name, *coverage = line.strip("\n").split("\t")
            cov[name] = [int(x) for x in coverage]
    return cov


def get_normalizer_values(fname, genes):
    # Store the counts for the housekeeping genes, for each strand
    expression_counts = list()

    with open(fname) as fin:
        # Skip the headers
        for _ in range(4):
            next(fin)

        for line in fin:
            spline = line.strip("\n").split("\t")
            ensg = spline[0]
            expression = [int(spline[column]) for column in [1, 2, 3]]
            if ensg not in genes:
                continue
            expression_counts.append(expression)

    unstranded = statistics.median(row[0] for row in expression_counts)
    forward = statistics.median(row[1] for row in expression_counts)
    reverse = statistics.median(row[2] for row in expression_counts)

    return unstranded, forward, reverse


def get_names_ensg(fname):
    with open(fname) as fin:
        ensg_to_name = gene_id_name(fin)
    return {v: k for k, v in ensg_to_name.items()}


def main(coverage_file, counts_file, housekeeping_genes, gtf_file):
    # Get the coverage from the bam file
    coverage = get_coverage(coverage_file)

    # Determine the ENSG for the housekeeping genes
    names_to_ensg = get_names_ensg(gtf_file)

    # This will crash if the specified gene is not in the GTF
    housekeeping_ensg = [names_to_ensg[x] for x in housekeeping_genes]

    # Get the mediant expression of the housekeeping genes
    median_housekeeping = get_normalizer_values(counts_file, housekeeping_ensg)

    for gene, cov in coverage.items():
        norm_values = [count / norm for count, norm in zip(cov, median_housekeeping)]
        print(gene, *norm_values, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--coverage", required=True, help="Coverage file for genes of interest"
    )
    parser.add_argument("--counts", required=True, help="STAR counts file")
    parser.add_argument(
        "--housekeeping-genes",
        required=True,
        nargs="+",
        help="Genes to use to normalize the expression",
    )
    parser.add_argument(
        "--gtf", required=True, help="Needed to lookup the ENSG for the spcified genes"
    )

    args = parser.parse_args()

    main(args.coverage, args.counts, args.housekeeping_genes, args.gtf)
