#!/usr/bin/env python

import sys


if __name__ == "__main__":

    fuma_fname = sys.argv[1]
    sf_fname = sys.argv[2]

    isect = set([])
    with open(fuma_fname, "r") as src:
        src.readline()
        for line in (raw_line.strip() for raw_line in src):
            cols = filter(lambda x: len(x) > 0, line.split("\t"))
            if len(cols) < 5:
                continue
            lefts = cols[0].split(":")
            rights = cols[1].split(":")
            for l in lefts:
                for r in rights:
                    isect.add(tuple(sorted([l, r])))

    with open(sf_fname, "r") as src:
        print(src.readline().strip())
        for line in (raw_line.strip() for raw_line in src):
            cols = line.split("\t")
            fusion = tuple(sorted([cols[4].split("^")[0], cols[6].split("^")[0]]))
            if fusion in isect:
                print(line)
