#!/usr/bin/env python3

import argparse

def main(coverage_file, counts_file, strand, housekeeping_genes):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--coverage", required=True, help="Coverage file for genes of interest"
    )
    parser.add_argument("--counts", required=True, help="STAR counts file")
    parser.add_argument("--strand", required=True, help="Strandedness of the sample")
    parser.add_argument("--housekeeping-genes", required=True, help="Genes to use to normalize the expression")

    args = parser.parse_args()

    main(args.coverage, args.counts, args.strand, args.housekeeping_genes)
