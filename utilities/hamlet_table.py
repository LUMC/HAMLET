#!/usr/bin/env python3

import argparse
from collections.abc import Sequence
import json
import os
import functools
from typing import Any, Dict, cast


def main(args: argparse.Namespace) -> None:
    functions: Dict[str, Any] = {
        "variant": print_variant_table,
        "fusion": print_fusion_table,
        "itd": print_itd_table,
        "expression": print_expression_table,
        "celltype": print_celltype_table,
        "aml_subtype": print_aml_subtype_table,
    }
    if args.table == "itd":
        functions[args.table](args.json_files, args.itd_gene)
    elif args.table in functions:
        functions[args.table](args.json_files)
    elif args.table == "all":
        os.makedirs(args.output, exist_ok=True)

        modules = "variant fusion expression celltype aml_subtype".split()
        for m in modules:
            fname = f"{args.output}/{m}.tsv"
            with open(fname, "wt") as fout:
                write = functools.partial(print, file=fout)
                functions[m](args.json_files, write)

        for gene in ["flt3", "kmt2a"]:
            fname = f"{args.output}/{gene}_itd.tsv"
            with open(fname, "wt") as fout:
                write = functools.partial(print, file=fout)
                functions["itd"](args.json_files, gene, write)

    else:
        raise NotImplementedError(args.table)


def print_aml_subtype_table(json_files: Sequence[str], write: Any = print) -> None:
    header = None
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
            subtype = data["modules"]["expression"]["subtype"]
            subtype["sample"] = subtype.pop("sample_id")

            if header is None:
                header = ["sample", "prediction", "pass_cutoff"]
                clusters = [field for field in subtype if field not in header]

                header += clusters
                write(*header, sep="\t")

            write(*(subtype[field] for field in header), sep="\t")


def print_celltype_table(json_files: Sequence[str], write: Any = print) -> None:
    header = None
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)
            celltypes = data["modules"]["expression"]["cell-types"]
            sample = sample_name(data)

            if header is None:
                header = list(celltypes["data"].keys())
                write("sample", *header, sep="\t")

            write(
                sample, *(celltypes["data"][celltype] for celltype in header), sep="\t"
            )


def print_expression_table(json_files: Sequence[str], write: Any = print) -> None:
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

                write(*header, sep="\t")

            # Print the data
            row = [sample]
            # First we append the normalized data in order
            for gene in genes:
                row.append(expression[gene]["normalized"])
            for gene in genes:
                row.append(expression[gene]["raw"])

            write(*row, sep="\t")


def print_variant_table(json_files: Sequence[str], write: Any = print) -> None:
    """Print variant table"""

    # Print the header
    header = "sample gene_name hgvsc hgvsp hgvsg database vaf exon annotation ref_depth alt_depth total_depth".split()
    write(*header, sep="\t")

    for fname in json_files:
        with open(fname) as fin:
            js = json.load(fin)

        if "modules" in js:
            genes = js["modules"]["snv_indels"]["genes"]
            name = sample_name(js)
        elif "snv_indels" in js:
            genes = js["snv_indels"]["genes"]
            name = sample_name(js["snv_indels"])
        else:
            raise RuntimeError("Unknown json format")

        for gene, variants in genes.items():
            for variant in variants:
                # Dict to store all the column values we will print
                ref, alt = variant["FORMAT"]["AD"].split(",")

                db_ids = ",".join(
                    (var["id"] for var in variant.get("colocated_variants", []))
                )
                to_print = {
                    "sample": name,
                    "gene_name": gene,
                    "total_depth": variant["FORMAT"]["DP"],
                    "vaf": variant["FORMAT"]["AF"],
                    "ref_depth": ref,
                    "alt_depth": alt,
                    "database": db_ids,
                }

                for transcript in variant["transcript_consequences"]:
                    to_print["hgvsc"] = transcript.get("hgvsc", "")
                    to_print["hgvsp"] = transcript.get("hgvsp", "")
                    to_print["hgvsg"] = transcript.get("hgvsg", "")

                    to_print["exon"] = transcript.get("exon", "")
                    to_print["annotation"] = transcript.get("annotation", "")
                write(*(to_print[field] for field in header), sep="\t")


def print_fusion_table(json_files: Sequence[str], write: Any = print) -> None:
    """Print fusion table"""
    header = [
        "sample",
        "gene1",
        "gene2",
        "strand1(gene/fusion)",
        "strand2(gene/fusion)",
        "breakpoint1",
        "breakpoint2",
        "site1",
        "site2",
        "type",
        "split_reads1",
        "split_reads2",
        "discordant_mates",
        "coverage1",
        "coverage2",
        "confidence",
        "reading_frame",
        "tags",
        "retained_protein_domains",
        "closest_genomic_breakpoint1",
        "closest_genomic_breakpoint2",
        "gene_id1",
        "gene_id2",
        "transcript_id1",
        "transcript_id2",
        "direction1",
        "direction2",
        "filters",
        "fusion_transcript",
        "peptide_sequence",
    ]

    write(*header, sep="\t")
    for js in json_files:
        with open(js) as fin:
            data = json.load(fin)

        if "modules" in data:  # HAMLET 2.0
            fusions = data["modules"]["fusion"]["events"]
            sample = sample_name(data)
        elif "results" in data:  # HAMLET 1.0
            fusions = data["results"]["fusion"]["tables"]["intersection"]["top20"]
            sample = sample_name(data)
        elif "fusion" in data:  # Output of the fusion module
            fusions = data["fusion"]["events"]
            sample = sample_name(data["fusion"])
        else:
            raise ValueError

        for fusion in fusions:
            fusion["sample"] = sample
            write(*(fusion[key] for key in header), sep="\t")


def sample_name(data: Dict[str, Any]) -> str:
    return cast(str, data["metadata"]["sample_name"])


def print_itd_table(
    json_files: Sequence[str], itd_gene: str, write: Any = print
) -> None:
    def join_list(positions: Sequence[Any]) -> Any:
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
        elif "results" in data:  # HAMLET 1.0
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
                write(*header, sep="\t")
            else:
                new_header = list(event.keys())
                if new_header != header:
                    raise RuntimeError()
            write(*(event[field] for field in header), sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "table",
        choices=[
            "variant",
            "fusion",
            "itd",
            "expression",
            "celltype",
            "aml_subtype",
            "all",
        ],
        default="variant",
        help="Table to output",
    )
    parser.add_argument("json_files", nargs="+")
    parser.add_argument("--itd-gene", required=False, choices=["flt3", "kmt2a"])
    parser.add_argument("--output", help="Output folder when using 'all'")

    args = parser.parse_args()

    if args.table == "itd" and not args.itd_gene:
        raise parser.error("Please specify an itd gene with --itd-gene")

    if args.table == "all" and not args.output:
        raise parser.error("Please specify an --output folder")

    main(args)
