#!/usr/bin/env python

import json
import argparse
import os


def fusion_results(args):
    # Initialise the results dictionary
    with open(args.arriba) as fin:
        return json.load(fin)

    # Add the results for each tool
    for tool, plot, table in zip(args.tools, args.plots, args.tables):
        results["plots"][tool] = os.path.abspath(plot)

        if tool == "fusioncatcher":
            events = fusioncatcher.parse(table)
            results["tables"][tool] = {
                "path": os.path.abspath(table),
                "top20": [transform_fusioncatcher(event) for event in events[:20]]
            }
        else:
            events = star_fusion.parse(table)
            results["tables"][tool] = {
                "path": os.path.abspath(table),
                "top20": [transform_star_fusion(event) for event in events[:20]]
            }

    return results


def main(args):
    """ Create json output of fusion results """
    results = fusion_results(args)
    print(json.dumps({"fusion": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arriba', help='Arriba output converted to JSON')

    args = parser.parse_args()
    main(args)

