#!/usr/bin/env python3

import sys

def unplaced(chrom):
    """ Guess if a chromosome is unplaced, from the name """
    return chrom.startswith("GL") or chrom.startswith("KI")


def main():
    refflat = sys.argv[1]

    # Put the gene name in the first column
    order = [11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    with open(refflat) as fin:
        for line in fin:
            spline = line.strip('\n').split('\t')
            if not spline[11]:
                print(f"Skipping {spline[0]}, no gene name", file=sys.stderr)
                continue
            chromosome = spline[1]
            if unplaced(chromosome):
                print(f"Skipping {spline[0]}, unplaced chromosome", file=sys.stderr)
                continue
            print(*[spline[f] for f in order ], sep='\t')

if __name__ == '__main__':
    main()
