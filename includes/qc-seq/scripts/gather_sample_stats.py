#!/usr/bin/env python

import json
import sys

import click


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("rg_stats", nargs=-1,
                type=click.Path(exists=True, dir_okay=False))
@click.option("--name", type=str)
def main(rg_stats, name):
    rgs = []
    for rg in rg_stats:
        with open(rg, "r") as src:
            rgs.append(json.load(src))

    stats = {
        "name": name,
        "raw": {
            "R1": {
                "num_seq": sum([x["raw"]["R1"]["num_seq"] for x in rgs]),
            },
            "R2": {
                "num_seq": sum([x["raw"]["R2"]["num_seq"] for x in rgs]),
            },
        },
        "proc": {
            "R1": {
                "num_seq": sum([x["proc"]["R1"]["num_seq"] for x in rgs]),
            },
            "R2": {
                "num_seq": sum([x["proc"]["R2"]["num_seq"] for x in rgs]),
            },
        },
        "read_groups": rgs,
    }
    json.dump(stats, sys.stdout, separators=(",", ":"))


if __name__ == "__main__":
    main()
