#!/usr/bin/env python

import json
import argparse
from pathlib import Path
from typing import Union

from multiqc import read_expression


def parse_deconvolution(fname: str) -> dict[str, float]:
    with open(fname) as fin:
        header = next(fin).strip("\n").replace('"', "").split(",")[1:]
        spline = next(fin).strip("\n").split(",")
        data = [float(x) for x in spline[1:]]

        assert len(data) == len(header)
        d = {k: v for k, v in zip(header, data)}
    return d


AMLmapRType = dict[str, Union[str, float, bool]]


def parse_amlmapr(fname: str) -> AMLmapRType:
    with open(fname) as fin:
        header = next(fin).strip("\n").replace('"', "").split(",")
        data = next(fin).strip("\n").replace('"', "").split(",")

    fix_type: list[Union[str, float, bool]] = []

    # Convert to float
    value: Union[float, str]
    for x in data:
        try:
            value = float(x)
        except ValueError:
            value = x
        fix_type.append(value)
    # Convert to boolean
    fix_type[-2] = fix_type[-2] == "TRUE"

    assert len(header) == len(fix_type)
    return {k: v for k, v in zip(header, fix_type)}


def main(args: argparse.Namespace) -> None:
    """Create json output of expression results"""

    # Create results dictionary
    results = {
        "metadata": {"sample_name": args.sample},
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

        genes[gene] = {"raw": int(raw_coverage[gene]), "normalized": norm}

    results["gene-expression"] = genes
    results["cell-types"] = dict()
    results["cell-types"]["data"] = parse_deconvolution(args.deconvolution)
    results["cell-types"]["plot"] = str(Path(args.cell_types).resolve())
    results["subtype"] = parse_amlmapr(args.subtype)
    print(json.dumps({"expression": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", help="Raw coverage file")
    parser.add_argument("--norm-coverage", help="Normalized coverage file")
    parser.add_argument("--sample", help="Sample name")
    parser.add_argument("--strandedness", help="Strandedness of the sample")
    parser.add_argument("--genes", nargs="*", default=list(), help="genes to include")
    parser.add_argument("--deconvolution", help="seAMLess deconvolution results")
    parser.add_argument("--cell-types", help="seAMLess cell type bar chart")
    parser.add_argument("--subtype", help="AMLmapR subtype prediction results")

    args = parser.parse_args()
    main(args)
