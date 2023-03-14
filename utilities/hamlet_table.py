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
            "ALT",
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

        self.itd_fields = [
            "td_starts",
            "td_ends",
            "rose_start_count",
            "rose_end_count",
            "rose_start_pos",
            "rose_start_anchor_pos",
            "rose_end_pos",
            "rose_end_anchor_pos",
            "boundary_type",
            "fuzziness",
        ]

    @staticmethod
    def get_alt(ref, genotype):
        """Extract alt from genotype string"""
        gt = genotype.split("/")
        # Remove the reference call from the genotype
        alt = set(gt)-{ref}
        # Make sure that there is only a single ALT allele
        if len(alt) > 1:
            raise NotImplementedError(genotype)
        return alt.pop()


    @property
    def variants(self):
        for gene, variants in self.json["results"]["var"]["overview"].items():
            for var in variants:
                var["Gene"] = gene
                var["ALT"] = self.get_alt(var["REF"], var["genotype"])
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

    def itd(self, gene):
        for event in self.json["results"]["itd"][gene]["table"]:
            yield event


class HAMLET_V2(HAMLET_V1):
    def __init__(self, data):
        super().__init__(data)
        # Translate VEP json values to old HAMLET format
        self.vep_lookup = {
            "Gene": "Gene",
            "CHROM": "seq_region_name",
            "POS": "start",
            "HGVSc": "HGVSc",
            "HGVSp": "HGVSp",
            "REF": "REF",
            "ALT": "ALT",
            "genotype": "allele_string",
            "PVAL": "PVAL",
            "Existing_variation": "Existing_variation",
            "FREQ": "FREQ",
            "is_in_hotspot": "is_in_hotspot",
        }

    @property
    def variants(self):
        for gene, variants in self.json["modules"]["snv_indels"]["genes"].items():
            for v in variants:
                # Split variants so they always contain (at most) a single
                # transcript consequence, so we can print it as a table
                for var in self.split_by_consequence(v):
                    # Include reference base in variant
                    self.rewrite_indel(var)
                    # Get the transcript_consequence object
                    cons = var["transcript_consequences"][0]
                    # Format var to match the existing HAMLET output format
                    var["Gene"] = gene
                    var["POS"] = self.vcf_pos(var)
                    var["HGVSc"] = cons["hgvsc"]
                    var["HGVSp"] = cons["hgvsp"]
                    var["REF"] = cons["used_ref"]
                    var["ALT"] = cons["variant_allele"]
                    var["PVAL"] = var["FORMAT"]["PVAL"]
                    existing = [
                        x["id"] for x in var.get("colocated_variants", []) if "id" in x
                    ]
                    var["Existing_variation"] = existing
                    var["FREQ"] = var["FORMAT"]["FREQ"]
                    var["is_in_hotspot"] = None
                    d = {f: var[self.vep_lookup[f]] for f in self.variant_fields}
                    yield d

    @staticmethod
    def split_by_consequence(variant):
        """Yield a variant with a single transcript_consequence, for each
        transcript_consequence"""
        if "transcript_consequences" not in variant:
            yield variant
        for i in range(len(variant["transcript_consequences"])):
            newvar = variant.copy()
            newvar["transcript_consequences"] = [variant["transcript_consequences"][i]]
            yield newvar

    @staticmethod
    def rewrite_indel(var):
        """Rewrite an insertion/deletion to include the reference base

        VEP does not include the reference base in the allele_description from
        the JSON output, but we need this to be compatible with older version
        of the HAMLET output.

        Insertion:
        ---------
        VEP:
        allele_string: -/T
        start: 105234996
        end:   105234996

        VCF:
        REF: T
        ALT: TT

        Deletion:
        ---------
        VEP:
        allele_string: T/-
        start: 32431640
        end:   32431639

        VCF:
        REF: AT
        ALT: A
        POS: 32431639
        """
        supported_classes = {"SNV", "insertion", "deletion"}
        if var["variant_class"] not in supported_classes:
            raise NotImplementedError

        # We do not have to re-write SNVs
        if var["variant_class"] == "SNV":
            return

        # Get the reference_sequence from the VCF field
        ref_seq = var["input"].split("\t")[3]

        # Get the ref/alt call from VEP out of the allele_string
        ref, alt = var["allele_string"].split("/")

        if var["variant_class"] == "insertion":
            # Prepend the reference base to the allele description
            var["allele_string"] = f"{ref_seq}/{ref_seq}{alt}"
        if var["variant_class"] == "deletion":
            # Include the reference base in the allele description
            var["allele_string"] = f"{ref_seq}/{ref_seq[0]}"

        # Shift the start/end positions
        var["start"] -= 1
        var["end"] -= 1

    @staticmethod
    def vcf_pos(var):
        c = var["variant_class"]
        if c == "SNV":
            return var["start"]
        elif c == "insertion" or c == "deletion" or c == "sequence variation":
            return var["start"] - 1

def main(args):
    if args.version == "v1":
        HAMLET = HAMLET_V1
    elif args.version == "v2":
        HAMLET = HAMLET_V2
    else:
        raise NotImplementedError

    if args.table == "variant":
        print_variant_table(HAMLET, args.json_files)
    elif args.table == "fusion":
        print_fusion_table(HAMLET, args.json_files)
    elif args.table == "overexpression":
        print_overexpression_table(HAMLET, args.json_files)
    elif args.table == "itd":
        print_itd_table(HAMLET, args.json_files, args.itd_gene)


def print_variant_table(HAMLET, json_files):
    """Print variant table"""
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
    parser.add_argument("--version", default="v2", help="HAMLET version")
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
