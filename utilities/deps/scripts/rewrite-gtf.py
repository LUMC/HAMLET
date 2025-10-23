#!/usr/bin/env python3

import sys


def main():
    infile = sys.argv[1]

    # How to rename the regular chromosomes
    mapping = {str(x): f"chr{x}" for x in range(1, 23)}
    mapping["MT"] = "chrM"
    mapping["X"] = "chrX"
    mapping["Y"] = "chrY"

    with open(infile) as fin:
        for line in fin:
            if line.startswith("#"):
                print(line.strip())
                continue

            spline = line.strip().split("\t")
            chrom = spline[0]
            if chrom in mapping:
                print(mapping[chrom], *spline[1:], sep="\t")
            else:
                print(line, end="")


if __name__ == "__main__":
    main()
