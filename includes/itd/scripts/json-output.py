#!/usr/bin/env python

import json
from pathlib import Path

import click


def add_itd_table(csv_fname):
    rv = []
    with open(csv_fname, "r") as src:
        header_cols = next(src).strip().split("\t")
        for line in (l.strip() for l in src):
            d = dict(zip(header_cols, line.split("\t")))

            # Integer fields
            int_fields = [
                "fuzziness", "rose_end_anchor_pos", "rose_end_count",
                "rose_end_pos", "rose_start_anchor_pos", "rose_start_count",
                "rose_start_pos"
            ]

            # Convert values to int
            for field in int_fields:
                d[field] = int(d[field])

            # Convert to list of ints
            d["td_ends"] = [int(x) for x in d["td_ends"].split(",")]
            d["td_starts"] = [int(x) for x in d["td_starts"].split(",")]

            rv.append(d)
    return rv


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("flt3_csv",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("flt3_plot",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("kmt2a_csv",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("kmt2a_plot",
                type=click.Path(exists=True, dir_okay=False))
def main(flt3_csv, flt3_plot, kmt2a_csv, kmt2a_plot):
    """Helper script for combining multiple stats files into one JSON."""
    combined = dict()
    combined["itd"] = {
        "flt3": {"path": str(Path(flt3_plot).resolve()),
                 "table": add_itd_table(flt3_csv)},
        "kmt2a": {"path": str(Path(kmt2a_plot).resolve()),
                  "table": add_itd_table(kmt2a_csv)},
    }
    print(json.dumps(combined, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()

