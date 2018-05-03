#!/usr/bin/env python

import json
import sys
import statistics as stats
from collections import namedtuple
from typing import Dict, Iterable, Optional, TextIO, Tuple

import click


class Row(namedtuple("Row",
                     ["chrom", "start", "end",
                      "feature", "exon_num", "exon_pos", "cov"])):

    @classmethod
    def from_raw_line(cls, line: str) -> "Row":
        cols = line.strip().split("\t")
        return cls(cols[0], int(cols[1]), int(cols[2]), cols[3], cols[4],
                   int(cols[5]) - 1, int(cols[6]))

    @property
    def key(self) -> Tuple[str, int]:
        return (self.feature, self.exon_num)

    @property
    def genome_pos(self) -> Tuple[str, int]:
        return (self.chrom, self.start + self.exon_pos)


def aggr_covs_entry(entry: dict, cov_limits: Iterable[int]) -> dict:
    covs = entry.pop("covs")
    metrics = {k: None for k in ("min", "max", "avg", "median", "stdev")}
    metrics["count"] = 0
    metrics["frac_cov_at_least"] = {f"{k}x": 0 for k in cov_limits}
    metrics["len"] = entry["end"] - entry["start"]

    if covs:
        covs.sort()
        count = len(covs)
        avg = stats.mean(covs)
        median = stats.median(covs)
        try:
            stdev = stats.stdev(covs, xbar=avg)
        except stats.StatisticsError:
            stdev = None

        limit_counts = {k: 0 for k in cov_limits}
        for cov in covs:
            for limit in cov_limits:
                if cov >= limit:
                    limit_counts[limit] += 1

        metrics = {
            "count": count,
            "min": covs[0],
            "max": covs[-1],
            "avg": avg,
            "median": median,
            "stdev": stdev,
            "frac_cov_at_least": {f"{k}x": v / count
                                  for k, v in limit_counts.items()}
        }

    entry["metrics"] = metrics
    return entry


def parse_idm(idm_fh: TextIO) -> Dict[str, str]:
    idms = {}
    for lineno, line in enumerate(idm_fh):
        if lineno == 0:
            continue
        gid, gsym, raw_tids = line.strip().split("\t")
        tids = raw_tids.split(",")
        for tid in raw_tids.split(","):
            assert tid not in idms, tid
            idms[tid] = gsym
    return idms


def group_per_exon(input_fh: TextIO, idm_fh: Optional[TextIO]=None,
                   cov_limits: Iterable[int]=(8, 10, 20, 30, 40, 50)):
    idms = {} if idm_fh is None else parse_idm(idm_fh)
    grouped = {}

    def include_row(row):
        if idms:
            return row.feature in idms
        return True

    for idx, line in enumerate(input_fh, start=1):
        row = Row.from_raw_line(line)

        if include_row(row):
            if row.key not in grouped:
                grouped[row.key] = {
                    "chrom": row.chrom, "start": row.start, "end": row.end,
                    "gx": idms[row.feature], "trx": row.feature,
                    "exon_num": row.exon_num,
                    "covs": []}
            grouped[row.key]["covs"].append(row.cov)

        if idx % 1_000_000 == 0 or idx == 1:
            show = f"{idx // 1000000}M lines" if idx != 1 else f"{idx} line"
            print(f"processed {show} ...", file=sys.stderr)

    print(f"processed {idx:,} lines in total", file=sys.stderr)
    print("aggregating coverage values ...", file=sys.stderr)

    return {k: aggr_covs_entry(v, cov_limits)
            for k, v in grouped.items()}


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("input_tsv", type=click.File("r"))
@click.argument("output", type=click.File("w"), default="-")
@click.option("-m", "--id-mapping", type=click.File("r"),
              help="Path to ID mapping file.")
@click.option("-cm", "--cov-limit", type=int, multiple=True,
              default=[8, 10, 20, 30, 40, 50],
              help="Values at which fraction coverage will be calculated.")
def main(input_tsv, output, id_mapping, cov_limit):
    """Calculates exon-level coverage metrics.

    The input to this script is a TSV file with the following columns:
        1. chromosome name
        2. start coordinate (0-based, half open)
        3. end coordinate (0-based, half open)
        4. transcript name
        5. exon number (1-based)
        6. position in exon (1-based)
        7. depth of coverage at the position

    If using bedtools coverage (as of v2.27.1), such a file can be obtained
    using the following command (assuming a position-sorted BED and BAM):

        bedtools coverage -d -sorted -a <bed> -b <bam> | cut -f1,2,3,4,5,8,7

    This script works with streaming input. Consequently, it assumes the input
    is grouped by chromosome and then sorted by start and stop coordinates,
    in ascending order.

    """
    grouped = group_per_exon(input_tsv, id_mapping, cov_limit)

    def serialize_key(row_key):
        return f"{row_key[0]}|{row_key[1]}"

    json.dump({serialize_key(k): v for k, v in grouped.items()},
              output, indent=None, separators=(",", ":"))


if __name__ == "__main__":
    main()

