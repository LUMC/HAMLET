#!/usr/bin/env python

import argparse
import json
from pathlib import Path
from typing import Any


ITDresult = dict[str, Any]


def add_itd_table(csv_fname: str) -> list[ITDresult]:
    rv = []
    with open(csv_fname, "r") as src:
        header_cols = next(src).strip().split("\t")
        for line in (l.strip() for l in src):
            d: ITDresult = dict(zip(header_cols, line.split("\t")))

            # Integer fields
            int_fields = [
                "fuzziness",
                "rose_end_anchor_pos",
                "rose_end_count",
                "rose_end_pos",
                "rose_start_anchor_pos",
                "rose_start_count",
                "rose_start_pos",
            ]

            # Convert values to int
            for field in int_fields:
                d[field] = int(d[field])

            # Convert to list of ints
            d["td_ends"] = [int(x) for x in d["td_ends"].split(",")]
            d["td_starts"] = [int(x) for x in d["td_starts"].split(",")]

            rv.append(d)
    return rv


def main(
    flt3_csv: str, flt3_plot: str, kmt2a_csv: str, kmt2a_plot: str, sample_name: str
) -> None:
    """Helper script for combining multiple stats files into one JSON."""
    combined = dict()
    combined["itd"] = {
        "flt3": {
            "path": str(Path(flt3_plot).resolve()),
            "table": add_itd_table(flt3_csv),
        },
        "kmt2a": {
            "path": str(Path(kmt2a_plot).resolve()),
            "table": add_itd_table(kmt2a_csv),
        },
        "metadata": {"sample_name": sample_name},
    }
    print(json.dumps(combined, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("flt3_csv")
    parser.add_argument("flt3_plot")
    parser.add_argument("kmt2a_csv")
    parser.add_argument("kmt2a_plot")
    parser.add_argument("--sample")

    args = parser.parse_args()
    main(args.flt3_csv, args.flt3_plot, args.kmt2a_csv, args.kmt2a_plot, args.sample)
