#!/usr/bin/env python

import json
import argparse

def parse_top20(tp):
    with open(tp, "r") as src:
        sft = []
        for lineno, line in enumerate(src):
            if lineno == 0:
                continue
            name, jr_count, sf_count, fusion_type, *_ = line.split("\t")
            sft.append({"name": name, "jr_count": int(jr_count),
                        "sf_count": int(sf_count), "type": fusion_type,})
            if lineno > 20:
                break
    return sft

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
        results["tables"][tool] = {
            "path": table,
            "top20": parse_top20(table)
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

