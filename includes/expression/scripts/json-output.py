#!/usr/bin/env python

import json
import argparse

from multiqc import read_expression

def parse_deconvolution(fname):
    with open(fname) as fin:
        header = next(fin).strip("\n").replace('"','').split(',')[1:]
        data = next(fin).strip("\n").split(",")
        data = [float(x) for x in data[1:]]

        assert len(data) == len(header)
        d = {k:v for k, v in zip(header, data)}
    return d

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
    results["cell-types"] = parse_deconvolution(args.deconvolution)
    print(json.dumps({"expression": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--coverage', help='Raw coverage file')
    parser.add_argument('--norm-coverage', help='Normalized coverage file')
    parser.add_argument('--sample', help = 'Sample name')
    parser.add_argument('--strandedness', help='Strandedness of the sample')
    parser.add_argument('--genes', nargs='*', default=list(), help='genes to include')
    parser.add_argument('--deconvolution', help="seAMLess deconvolution results")

    args = parser.parse_args()
    main(args)

