#!/usr/bin/env python

import json
import sys

import click


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("rg_stats", nargs=-1,
                type=click.Path(exists=True, dir_okay=False))
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
    json.dump({ "qc_seq": stats}, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
