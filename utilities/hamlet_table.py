#!/usr/bin/env python3

import argparse
import json


def main(args):
    if args.table == "variant":
        print_variant_table(args.json_files)
    elif args.table == "fusion":
        print_fusion_table(args.json_files)
    elif args.table == "itd":
        print_itd_table(args.json_files, args.itd_gene)
    elif args.table == "expression":
        print_expression_table(args.json_files)


def print_expression_table(json_files):
    """Print gene expression table"""
    genes = None
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
            expression = data["modules"]["expression"]["gene-expression"]
            sample = sample_name(data)
            if genes is None:
                genes = list(expression.keys())

                # Print the header
                header = ["sample"]
                #  First add all normalized columns
                for gene in genes:
                    header.append(f"{gene}-normalized")
                # Then add all raw columncs
                for gene in genes:
                    header.append(f"{gene}-raw")

                print(*header, sep="\t")

            # Print the data
            row = [sample]
            # First we append the normalized data in order
            for gene in genes:
                row.append(expression[gene]["normalized"])
            for gene in genes:
                row.append(expression[gene]["raw"])

            print(*row, sep="\t")

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
                        values[field] = transcript_consequence.get(field, '')
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
                header = get_fields(var_line)
                print(*header, sep="\t")
            else:
                current_keys = get_fields(var_line)
                if current_keys != header:
                    msg = f"\n{current_keys}\n{header}\n do not match!"
                    raise RuntimeError(msg)
            print(*(var_line[field] for field in header), sep="\t")

def get_fields(var_line):
    fields = list(var_line.keys())
    # minimised is not set for every variant, so we remove it
    try:
        fields.remove('minimised')
    except ValueError:
        pass
    return fields

def print_fusion_table(json_files):
    """Print fusion table"""
    header = ["sample", "gene1", "gene2", "strand1(gene/fusion)", "strand2(gene/fusion)", "breakpoint1", "breakpoint2", "site1", "site2", "type", "split_reads1", "split_reads2", "discordant_mates", "coverage1", "coverage2", "confidence", "reading_frame", "tags", "retained_protein_domains", "closest_genomic_breakpoint1", "closest_genomic_breakpoint2", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "direction1", "direction2", "filters", "fusion_transcript", "peptide_sequence"]

    print(*header, sep="\t")
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)

        sample = sample_name(data)
        if "modules" in data: # HAMLET 2.0
            fusions = data["modules"]["fusion"]["events"]
        elif "results" in data: # HAMLET 1.0
            fusions = data["results"]["fusion"]["tables"]["intersection"]["top20"]
        else:
            raise ValueError

        for fusion in fusions:
            fusion["sample"] = sample
            print(*(fusion[key] for key in header), sep="\t")

def sample_name(data):
    return data["metadata"]["sample_name"]

def print_itd_table(json_files, itd_gene):
    def join_list(positions):
        """Join a list of positions

        Also handle the case where we don't get a list at all
        """
        if isinstance(positions, list):
            return ",".join((str(x) for x in positions))
        else:
            return positions

    # Did we print the header already
    header = None

    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)

        if "modules" in data:  # HAMLET 2.0
            itd_table = data["modules"]["itd"][itd_gene]["table"]
        elif "results" in data: # HAMLET 1.0
            itd_table = data["results"]["itd"][itd_gene]["table"]

        for event in itd_table:
            event["td_ends"] = join_list(event["td_ends"])
            event["td_starts"] = join_list(event["td_starts"])

            # Sneakily put the sample name first
            new_event = {"sample": sample_name(data)}
            new_event.update(event)
            event = new_event

            if header is None:
                header = list(event.keys())
                print(*header, sep="\t")
            else:
                new_header = list(event.keys())
                if new_header != header:
                    raise RuntimError()
            print(*(event[field] for field in header), sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "table",
        choices=["variant", "fusion", "itd", "expression"],
        default="variant",
        help="Table to output",
    )
    parser.add_argument("json_files", nargs="+")
    parser.add_argument("--itd-gene", required=False, choices=["flt3", "kmt2a"])

    args = parser.parse_args()

    if args.table == "itd" and not args.itd_gene:
        raise parser.error("Please specify an itd gene with --itd-gene")

    main(args)
