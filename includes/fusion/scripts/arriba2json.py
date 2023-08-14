#!/usr/bin/env python3

import sys
import json

def to_list(d, field, sep=','):
    """Split values in d[field] into a list"""
    if not d[field]:
        return
    d[field] = d[field].split(sep)


def parse_arriba(fin):
    header = next(fin)[1:-1].split('\t')
    for line in fin:
        d = {k: v for k, v in zip(header, line.strip('\n').split('\t'))}

        # Convert '.' to None
        for key, value in d.items():
            d[key] = value if value != "." else None

        # Convert lists
        to_list(d, "read_identifiers")
        to_list(d, "filters")

        # Convert to int
        for field in ["split_reads1", "split_reads2", "coverage1", "coverage2", "discordant_mates"]:
            d[field] = int(d[field])

        # Update filters, if there are any
        if d["filters"]:
            filters = dict()
            for f in d["filters"]:
                name, removed_reads = f.split("(")
                filters[name] = int(removed_reads[:-1])
            d["filters"] = filters
        yield d


if __name__ == "__main__":
    with open(sys.argv[1]) as fin:
        # We want to keep the fusions in the same order as the input file
        fusions = list()
        for record in parse_arriba(fin):
            fusions.append(record)
        print(json.dumps(fusions, indent=True))
