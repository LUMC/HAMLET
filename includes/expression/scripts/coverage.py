#!/usr/bin/env python3

import argparse
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
            spline = line.strip("\n").split("\t")

            # To int
            for field in [1, 2, 4]:
                spline[field] = int(spline[field])

            B = Bed(*spline)

            yield B


def count_reads(Bed, bamfile):
    """Count the number of reads in region"""

    # Store the read name to ensure we don't count forward/reverse double
    names = set()

    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for read in samfile.fetch(contig=Bed.chrom, start=Bed.chromStart, end=Bed.chromEnd):
            names.add(read.qname)

    return len(names)


def main(bed: str, bamfile: str):
    # If there are multiple regions in the BED file with the same name, we add
    # the counts together
    counts = defaultdict(int)

    for Bed in parse_bed(bed):
        count = count_reads(Bed, bamfile)
        counts[Bed.name] += count

    for name, count in counts.items():
        print(name, count, sep=",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bed", required=True, help="BED file with regions to determine coverage for"
    )
    parser.add_argument("--bam", required=True, help="Bam file")

    args = parser.parse_args()

    main(args.bed, args.bam)
