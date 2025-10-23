#!/usr/bin/env python3

import sys


class Bed:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.size = self.end - self.start

    # def __repr__(self):
    #     chrom = self.chrom
    #     start = self.start
    #     end = self.end
    #     size = self.size
    #     return f"Bed({chrom=}, {start=}, {end=}, {size=})"

    def __str__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}"

    def __gt__(self, other):
        return self.size > other.size


def main():
    # Read all bed records
    records = list()
    for line in sys.stdin:
        row = line.strip("\n").split("\t")
        record = Bed(*row[:3])
        records.append(record)

    for record in sorted(records, reverse=True):
        print(record)


if __name__ == "__main__":
    main()
