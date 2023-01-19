#!/usr/bin/env python3

import argparse
import json


class HAMLET_V1:
    def __init__(self, data):
        self.json = data
        self.variant_fields = [
            "Gene",
            "CHROM",
            "POS",
            "HGVSc",
            "HGVSp",
            "REF",
            "genotype",
            "PVAL",
            "Existing_variation",
            "FREQ",
            "is_in_hotspot",
        ]
        self.fusion_fields = ["jr_count", "name", "sf_count", "type"]

        self.overexpression_fields = [
            "exon",
            "divisor_gene",
            "count",
            "divisor_exp",
            "ratio",
            "above_threshold",
        ]

    @property
    def variants(self):
        for gene, variants in self.json["results"]["var"]["overview"].items():
            for var in variants:
                var["Gene"] = gene
                d = {f: var[f] for f in self.variant_fields}
                d["POS"] = int(d["POS"])
                d["Existing_variation"] = var["Existing_variation"].split(",")
                d["FREQ"] = float(var["FREQ"][:-1]) / 100
                d["is_in_hotspot"] = var["is_in_hotspot"] == "yes"
                yield d

    @property
    def sample(self):
        return self.json["metadata"]["sample_name"]

    @property
    def fusions(self):
        for fusion in self.json["results"]["fusion"]["tables"]["intersection"]["top20"]:
            fusion["jr_count"] = int(fusion["jr_count"])
            fusion["sf_count"] = int(fusion["sf_count"])
            yield fusion

    @property
    def overexpression(self):
        for expression in self.json["results"]["expr"]:
            d = {f: expression[f] for f in self.overexpression_fields}
            d["count"] = int(d["count"])
            d["ratio"] = float(d["ratio"])
            d["above_threshold"] = d["above_threshold"] == "yes"
            d["divisor_exp"] = int(d["divisor_exp"])
            yield d


def main(args):
    if args.version == "v1":
        HAMLET = HAMLET_V1
    else:
        raise NotImplementedError

    if args.table == "variant":
        print_variant_table(HAMLET, args.json_files)
    elif args.table == "fusion":
        print_fusion_table(HAMLET, args.json_files)
    elif args.table == "overexpression":
        print_overexpression_table(HAMLET, args.json_files)


def print_variant_table(HAMLET, json_files):
    """ Print variant table """
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
            print(H.sample, *[variant[field] for field in H.variant_fields], sep="\t")


def print_fusion_table(HAMLET, json_files):
    """ Print fusion table """
    # Did we print the header already
    header_printed = False

    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
        H = HAMLET(data)

        if not header_printed:
            print("sample", *H.fusion_fields, sep="\t")
            header_printed = True

        for fusion in H.fusions:
            print(H.sample, *[fusion[field] for field in H.fusion_fields], sep="\t")


def print_overexpression_table(HAMLET, json_files):
    """ Print overexpression table """
    # Did we print the header already
    header_printed = False

    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
        H = HAMLET(data)

        if not header_printed:
            print("sample", *H.overexpression_fields, sep="\t")
            header_printed = True

        for expr in H.overexpression:
            print(
                H.sample, *[expr[field] for field in H.overexpression_fields], sep="\t"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="v1", help="HAMLET version")
    parser.add_argument(
        "table",
        choices=["variant", "fusion", "overexpression"],
        default="variant",
        help="Table to output",
    )
    parser.add_argument("json_files", nargs="+")

    args = parser.parse_args()

    main(args)
