#!/usr/bin/env python

# Script for converting raw FLT3 JSON data to CSV.

import csv
import json
import sys

import click
import pandas as pd


START = 1786 # 0-based
END = 2024


def within_range(pos, start, end):
    if start is None and end is None:
        return True
    return start <= pos < end


def load_json(fname, sample_name, start, end):
    with open(fname, "r") as src:
        data = json.load(src)

    pileups_by_pos = {x["pos"]: {"count": x["count"]} for x in data["pileups"]}
    sc_by_pos = {x["pos"]: {k: x[k] for k in ("altPosCount", "count")}
                 for x in data["scs"] if within_range(x["pos"], start, end)}
    inserts_by_pos = {}
    for x in data["inserts"]:
        pos = x["pos"]
        if not within_range(pos, start, end):
            continue
        if pos not in inserts_by_pos:
            inserts_by_pos[pos] = []
        inserts_by_pos[pos].append({"count": x["count"], "seq": x["seq"]})

    dps = []
    for pos in range(start, end):
        pileup_count = pileups_by_pos.get(pos, {"count": 0})["count"]
        sc = sc_by_pos.get(pos, {"altPosCount": [], "count": 0})
        sc_count = sc["count"]
        try:
            sc_pct = sc_count * 100.0 / pileup_count
        except ZeroDivisionError:
            sc_pct = float("nan")

        links = [x for x in sc["altPosCount"] if within_range(x["pos"], start, end)]
        links_count = sum(x["count"] for x in links)
        try:
            lsc_count = sum(x["count"]
                            for x in sc["altPosCount"] if within_range(x["pos"], start, end))
            lsc_pct = lsc_count * 100.0 / pileup_count
        except ZeroDivisionError:
            lsc_pct = float("nan")

        inserts = inserts_by_pos.get(pos) or []
        inserts_count = sum(x["count"] for x in inserts)

        dps.append({
            "sample": sample_name,
            "pos": pos,
            "pileup_count": pileup_count,
            "sc_count": sc_count,
            "sc_pct": sc_pct,
            "lsc_pct": lsc_pct,
            "lsc_count": lsc_count,
            "insert_count": inserts_count,
            "links_count": links_count,

            "links": links,
            "inserts": inserts,
        })

    df = pd.DataFrame(dps).set_index(["pos", "sample"], drop=False,
                                     verify_integrity=True)
    return df


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("input",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("sample_id", type=str)
@click.option("--start", type=int, default=START)
@click.option("--end", type=int, default=END)
def main(input, sample_id, start, end):

    df = load_json(input, sample_id, start=start, end=end)
    df_cols = ["sample", "pos", "pileup_count", "sc_count", "lsc_count",
               "sc_pct"]
    header = ["sample", "position", "coverage", "n_soft_clips",
              "n_alt_soft_clip_aln", "pct_soft_clips"]
    # Convert coordinate to 1-based, fully-closed.
    df["pos"] = df["pos"].apply(lambda p: p + 1)
    df[df_cols].to_csv(sys.stdout, index=False, header=header)


if __name__ == "__main__":
    main()
