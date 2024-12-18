#!/usr/bin/env python

import json
import argparse
import os

from multiqc import read_expression

def fusion_results(fname):
    # Initialise the results dictionary
    with open(fname) as fin:
        return json.load(fin)

def main(args):
    """ Create json output of expression results """
    
    # Create results dictionary
    results = {
            "metadata": {
                "sample_name": args.sample
            },
    }
    
    # Read the normalized expression for the specified strand
    raw_coverage = read_expression(args.coverage, args.strandedness)
    norm_coverage = read_expression(args.norm_coverage, args.strandedness)

    # Extract the genes of interest
    genes = dict()
    for gene in args.genes:
        # If all housekeeping genes have 0 expression, the normalized expression is None
        if norm_coverage[gene] == "None":
            norm = None
        else:
            norm = float(norm_coverage[gene])

        genes[gene] = {
            "raw": int(raw_coverage[gene]),
            "normalized": norm
        }

    results["gene-expression"] = genes
    print(json.dumps({"expression": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--coverage', help='Raw coverage file')
    parser.add_argument('--norm-coverage', help='Normalized coverage file')
    parser.add_argument('--sample', help = 'Sample name')
    parser.add_argument('--strandedness', help='Strandedness of the sample')
    parser.add_argument('--genes', nargs='*', default=list(), help='genes to include')

    args = parser.parse_args()
    main(args)

