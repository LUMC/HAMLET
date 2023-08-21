#!/usr/bin/env python

import json
import argparse
import os


def fusion_results(args):
    # Initialise the results dictionary
    with open(args.arriba) as fin:
        return json.load(fin)

def main(args):
    """ Create json output of fusion results """
    results = fusion_results(args)
    print(json.dumps({"fusion": results}, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('arriba', help='Arriba output converted to JSON')

    args = parser.parse_args()
    main(args)

