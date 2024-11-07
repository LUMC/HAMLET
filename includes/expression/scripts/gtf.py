#!/usr/bin/env python3

import sys
import json

def read_attributes(fin):
    for line in fin:
        if line.startswith("#"):
            continue
        spline = line.strip().split(maxsplit=8)
        kvals = spline[-1]

        d = dict()
        for pair in kvals.split(";"):
            if not pair:
                continue
            pair = pair.strip(" ")
            k, v = pair.split(" ", maxsplit=1)
            d[k] = v.replace('"','')
        yield d


def gene_id_name(fin) -> dict[str, str]:
    """Return a two-way mapping between gene names and ID's"""
    mapping = dict()
    for record in read_attributes(fin):
        if "gene_id" in record and "gene_name" in record:
            # From geneID (ENSG) to gene name
            mapping[record["gene_id"]] = record["gene_name"]
    return mapping


if __name__ == "__main__":
    with open(sys.argv[1]) as fin:
        mapping = gene_id_name(fin)
    print(json.dumps(mapping, indent=True))
