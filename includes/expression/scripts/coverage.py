#!/usr/bin/env python3

import argparse
from typing import Any
import pysam

from collections import defaultdict
from dataclasses import dataclass


@dataclass
class Bed:
    """Class for BED records"""

    chrom: str
    chromStart: int
    chromEnd: int
    name: str
    score: int
    strand: str


def parse_bed(fname: str):
    """Parse BED file and yield records"""
    with open(fname) as fin:
        for line in fin:
            spline: list[Any] = line.strip("\n").split("\t")

            # To int
            for field in [1, 2, 4]:
                spline[field] = int(spline[field])

            B = Bed(*spline)

            yield B


def get_reads(Bed, bamfile):
    """Get reads that are not supplementary or secondary or not proper pairs"""
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for read in samfile.fetch(
            contig=Bed.chrom, start=Bed.chromStart, end=Bed.chromEnd
        ):
            if read.is_supplementary or read.is_secondary or not read.is_proper_pair:
                continue
            yield read


def read_names(Bed, bamfile):
    """Extract the read names in Bed region, per strandedness"""

    # Store the read name to ensure we don't count forward/reverse double
    by_strand = {"unstranded": set(), "forward": set(), "reverse": set()}

    for read in get_reads(Bed, bamfile):
        orientation = forward_orientation(read)
        read_name = read.query_name
        if orientation == Bed.strand:
            by_strand["forward"].add(read_name)
        else:
            by_strand["reverse"].add(read_name)

        # For unstranded, we count both orientations
        by_strand["unstranded"].add(read_name)

    return by_strand


def forward_orientation(read: pysam.AlignedSegment) -> str:
    """Determine the orientation of the forward read of the pair"""
    if read.is_read1 and read.is_forward:
        return "+"
    elif read.is_read2 and read.is_forward:
        return "-"
    elif read.is_read1 and read.is_reverse:
        return "-"
    elif read.is_read2 and read.is_reverse:
        return "+"

    msg = f"Unexpected read encountered: {read.query_name}"
    raise RuntimeError(msg)


def main(bed: str, bamfile: str):
    # If there are multiple regions in the BED file with the same name, we add
    # the read names together so we don't count double
    by_name = defaultdict(
        lambda: {"unstranded": set(), "forward": set(), "reverse": set()}
    )

    for Bed in parse_bed(bed):
        names = read_names(Bed, bamfile)
        for strand in ["unstranded", "forward", "reverse"]:
            by_name[Bed.name][strand].update(names[strand])

    for name, reads in by_name.items():
        u = reads["unstranded"]
        f = reads["forward"]
        r = reads["reverse"]

        print(name, len(u), len(f), len(r), sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bed", required=True, help="BED file with regions to determine coverage for"
    )
    parser.add_argument("--bam", required=True, help="Bam file")

    args = parser.parse_args()

    main(args.bed, args.bam)
