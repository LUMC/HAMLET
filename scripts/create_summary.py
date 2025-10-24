#!/usr/bin/env python

import argparse
import json

# Important: This scripts runs in a container that has Python 2.7


def parse_idm(fh):
    idms = []
    for lineno, line in enumerate(fh):
        if lineno == 0:
            continue
        gid, gsym, raw_tids = line.strip().split("\t")
        tids = raw_tids.split(",")
        idms.append({"gene_id": gid, "gene_symbol": gsym, "transcript_ids": tids})
    return idms


def main(args):
    """Helper script for combining multiple stats files into one JSON."""
    with open(args.id_mappings_path) as fin:
        idm = parse_idm(fin)
    combined = {
        "metadata": {
            "pipeline_version": args.pipeline_version,
            "sample_name": args.sample_name,
            "genes_of_interest": idm,
        },
        "modules": dict(),
    }

    # Read and add module json files
    for m in args.module:
        with open(m) as fin:
            module = json.load(fin)
        for key, data in module.items():
            combined["modules"][key] = data

    print(json.dumps(combined, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("id_mappings_path")
    parser.add_argument(
        "--sample-name", help="Name of the sample from which the stats were generated."
    )
    parser.add_argument("--pipeline-version", help="Version string of the pipeline.")
    parser.add_argument(
        "--module", action="append", help="JSON outputs from various modules"
    )

    args = parser.parse_args()
    main(args)
