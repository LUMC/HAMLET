#!/usr/bin/env python

import json
import argparse
import os


def fusion_results(fname):
    # Initialise the results dictionary
    with open(fname) as fin:
        return json.load(fin)

def main(args):
    """ Create json output of fusion results """
    results = fusion_results(args.arriba)
    for i, fusion in enumerate(results, 1):
        # Here, we magically know the name of the png file.
        plot = os.path.join(args.plots, f"fusion-{i}.png")
        if not os.path.exists(plot):
            raise RuntimeError(f"Missing fusion figure: {plot}")

        fusion["plot"] = os.path.abspath(plot)

    print(json.dumps({"fusion": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arriba', help='Arriba output converted to JSON')
    parser.add_argument('plots', help='Arriba fusion plots')

    args = parser.parse_args()
    main(args)

