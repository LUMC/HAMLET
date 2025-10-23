#!/usr/bin/env python

import json
import argparse
import os
from typing import Any, cast


ArribaResults = list[dict[str, Any]]


def fusion_results(fname: str) -> ArribaResults:
    # Initialise the results dictionary
    with open(fname) as fin:
        return cast(ArribaResults, json.load(fin))


def main(args: argparse.Namespace) -> None:
    """Create json output of fusion results"""
    results = fusion_results(args.arriba)

    # The padding in the png filename is dependent on the total number of
    # fusions
    nr_fusions = len(results)
    padding = len(str(nr_fusions))

    for count, fusion in enumerate(results, 1):
        # Add padding to the fusion count
        i = str(count).zfill(padding)
        # Here, we magically know the name of the png file.
        plot = os.path.join(args.plots, f"fusion-{i}.png")
        if not os.path.exists(plot):
            raise RuntimeError(f"Missing fusion figure: {plot}")

        fusion["plot"] = os.path.abspath(plot)

    data = {"fusion": {"events": results, "metadata": {"sample_name": args.sample}}}
    print(json.dumps(data, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("arriba", help="Arriba output converted to JSON")
    parser.add_argument("plots", help="Arriba fusion plots")
    parser.add_argument("--sample")

    args = parser.parse_args()
    main(args)
