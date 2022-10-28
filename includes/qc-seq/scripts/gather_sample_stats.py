#!/usr/bin/env python

import argparse
import json

def main(rg_stats):
    rgs = []
    for rg in rg_stats:
        with open(rg, "r") as src:
            rgs.append(json.load(src))

    stats = {
        "all_read_groups": {
            "raw": {
                "num_reads_r1": sum([x["raw"]["num_reads_r1"] for x in rgs]),
                "num_reads_r2": sum([x["raw"]["num_reads_r2"] for x in rgs])
            },
            "proc": {
                "num_reads_r1": sum([x["proc"]["num_reads_r1"] for x in rgs]),
                "num_reads_r2": sum([x["proc"]["num_reads_r2"] for x in rgs])
            },
        },
        "per_read_group": rgs,
    }
    print(json.dumps({ "qc_seq": stats}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rg_stats", nargs="+")
    args = parser.parse_args()
    main(args.rg_stats)
