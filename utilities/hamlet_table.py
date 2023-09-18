#!/usr/bin/env python3

import argparse
import json


def main(args):
    if args.table == "variant":
        print_variant_table(args.json_files)
    elif args.table == "fusion":
        print_fusion_table( args.json_files)
    elif args.table == "overexpression":
        print_overexpression_table(args.json_files)
    elif args.table == "itd":
        print_itd_table(args.json_files, args.itd_gene)


def print_variant_table(json_files):
    """Print variant table"""
    def parse_genes_v1(_, genes):
        """ Parse HAMLET variant results in v1 format"""
        for gene in genes:
            for variant in genes[gene]:
                yield variant

    def parse_genes_v2(sample, genes):
        """Parse HAMLET variant results in V2 format"""
        for gene, variants in genes.items():
            for variant in variants:
                for transcript_consequence in variant["transcript_consequences"]:
                    values = parse_variant(sample, gene, variant)
                    # Copy over hgvs
                    for field in ["hgvsc", "hgvsp", "impact"]:
                        values[field] = transcript_consequence[field]
                    # Copy over consequence terms
                    values["consequence_terms"] = ",".join(transcript_consequence["consequence_terms"])
                    # Throw out the 'input' field, since that is just the VCF
                    # information again
                    values.pop("input")
                    yield values

    def parse_colocated_variants(coloc):
        return ','.join((var["id"] for var in coloc))

    def parse_variant(sample, gene, variant):
        """ Parse a single variant"""
        var = dict()
        var["sample"] = sample
        var["gene"] = gene

        # Add headers for the VCF fields
        vcf_header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT FORMAT_DATA".split()
        var.update({field: value for field, value in zip(vcf_header, variant["input"].split("\t"))})

        # Copy over all simple values
        for field in variant:
            if isinstance(variant[field], (str, float, int)):
                var[field] = variant[field]

        # Copy over all colocated variants
        var["colocated_variants"] = parse_colocated_variants(variant.get("colocated_variants", list()))

        return var
    # The headers in the csv file are based on the json output
    header = None

    # Process every json file
    for fname in json_files:
        with open(fname) as fin:
            data = json.load(fin)

        # See if the data is from HAMLET 1.0 or 2.0
        if "modules" in data: # HAMLET 2.0
            genes = data["modules"]["snv_indels"]["genes"]
            parse = parse_genes_v2
        elif "results" in data: # HAMLET 1.0
            genes = data["results"]["var"]["overview"]
            parse = parse_genes_v1

        # Extract the sample name
        sample = data["metadata"]["sample_name"]

        # Print the headers if this is the first time
        # Next print every variant line
        for var_line in parse(sample, genes):
            if header is None:
                header = list(var_line.keys())
                print(*header, sep="\t")
            else:
                current_keys = list(var_line.keys())
                if current_keys != header:
                    msg = f"{current_keys}\n{dynamic_keys}\n do not match!"
                    raise RuntimeError(msg)
            print(*(var_line[field] for field in header), sep="\t")


def print_fusion_table(HAMLET, json_files):
    """Print fusion table"""
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
    """Print overexpression table"""
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


def print_itd_table(HAMLET, json_files, itd_gene):
    # Did we print the header already
    header_printed = False

    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
        H = HAMLET(data)

        if not header_printed:
            print("sample", *H.itd_fields, sep="\t")
            header_printed = True

        for expr in H.itd(itd_gene):
            print(H.sample, *[expr[field] for field in H.itd_fields], sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "table",
        choices=["variant", "fusion", "overexpression", "itd"],
        default="variant",
        help="Table to output",
    )
    parser.add_argument("json_files", nargs="+")
    parser.add_argument("--itd-gene", required=False)

    args = parser.parse_args()

    if args.table == "itd" and not args.itd_gene:
        raise parser.error("Please specify an itd gene with --itd-gene")

    main(args)
