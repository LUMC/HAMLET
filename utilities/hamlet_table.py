#!/usr/bin/env python3

import argparse
import json


class HAMLET_V1:
    def __init__(self, data):
        self.json = data
        self.variant_fields = [
            "CHROM",
            "POS",
            "HGVSc",
            "HGVSp",
            "PVAL",
            "SYMBOL",
            "genotype",
            "Existing_variation",
            "FREQ",
            "is_in_hotspot",
        ]

    @property
    def variants(self):
        for variants in self.json["results"]["var"]["overview"].values():
            for var in variants:
                d = {f: var[f] for f in self.variant_fields}
                d["POS"] = int(d["POS"])
                d["Existing_variation"] = var["Existing_variation"].split(",")
                d["FREQ"] = float(var["FREQ"][:-1]) / 100
                d["is_in_hotspot"] = var["is_in_hotspot"] == "yes"
                yield d

    @property
    def sample(self):
        return self.json["metadata"]["sample_name"]


def main(args):
    if args.version == "v1.0":
        HAMLET = HAMLET_V1
    else:
        raise NotImplementedError

    if args.table == "variant":
        print_variant_table(HAMLET, args.json_files)

def print_variant_table(HAMLET, json_files):
    """ Print variant tables """
    # Did we print the header already
    header_printed = False

    # Process every json file
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
        H = HAMLET(data)

        if not header_printed:
            print("sample", *H.variant_fields, sep="\t")
            header_printed = True

        for variant in H.variants:
            print(
                H.sample, *[variant[field] for field in H.variant_fields], sep="\t"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="v1.0", help="HAMLET version")
    parser.add_argument(
        "table", choices=["variant"], default="variant", help="Table to output"
    )
    parser.add_argument("json_files", nargs="+")

    args = parser.parse_args()

    main(args)
