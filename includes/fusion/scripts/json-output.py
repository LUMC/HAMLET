#!/usr/bin/env python

import json
import argparse

from crimson import star_fusion
from crimson import fusioncatcher

def transform_star_fusion(event):
    """ Transform a star_fusion event to a simplified form """
    return {
        "jr_count": event["nJunctionReads"],
        "name": event["fusionName"],
        "sf_count": event["nSpanningFrags"],
        "type": event["spliceType"]
    }

def transform_fusioncatcher(event):
    """ Transform a fusioncatcher event to a simplified form """
    return {
        "jr_count": None,
        "name": f"{event['5end']['geneSymbol']}--{event['3end']['geneSymbol']}",
        "sf_count": event["nSpanningPairs"],
        "type": None
    }

def fusion_results(args):
    # Initialise the results dictionary
    results = {
            "intersected": args.intersected,
            "plots": dict(),
            "tables": dict()
    }

    # Add the results for each tool
    for tool, plot, table in zip(args.tools, args.plots, args.tables):
        results["plots"][tool] = plot

        if tool == "fusioncatcher":
            events = fusioncatcher.parse(table)
            results["tables"][tool] = {
                "path": table,
                "top20": [transform_fusioncatcher(event) for event in events[:20]]
            }
        else:
            events = star_fusion.parse(table)
            results["tables"][tool] = {
                "path": table,
                "top20": [transform_star_fusion(event) for event in events[:20]]
            }

    return results


def main(args):
    """ Create json output of fusion results """
    results = fusion_results(args)
    print(json.dumps({"fusion": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--intersected', default=False, action='store_true')
    parser.add_argument('--tools', required=True, nargs='+', help='Tools')
    parser.add_argument('--plots', required=True, nargs='+', help='PNG plots, in the same order as --tools')
    parser.add_argument('--tables', required=True, nargs='+', help='Result tables, in the same order as --tools')

    args = parser.parse_args()
    assert len(args.tools) == len(args.plots)
    assert len(args.tools) == len(args.tables)
    main(args)

