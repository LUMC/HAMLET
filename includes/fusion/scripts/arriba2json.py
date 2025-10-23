#!/usr/bin/env python3

from typing import Any, Dict, Set
import argparse
import json

arriba_header = [
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
    "read_identifiers",
]


def arriba_to_json(header, line):
    """Convert an arriba line to json"""
    d = {k: v for k, v in zip(header, line.split("\t"))}

    # Convert '.' to None
    for key, value in d.items():
        d[key] = value if value != "." else None

    # Convert lists
    d["read_identifiers"] = d["read_identifiers"].split(",")

    # Convert to int, if possible
    for field in d.keys():
        try:
            d[field] = int(d[field])
        except (ValueError, TypeError):
            pass

    return d


def json_to_arriba(header, data):
    """Convert arriba data to a list of strings"""
    # Join read identifiers into single string
    data["read_identifiers"] = ",".join(data["read_identifiers"])

    # Convert 'None' values to a dot
    for field, value in data.items():
        if value is None:
            data[field] = "."

    return "\t".join((str(data[field]) for field in header))


def parse_arriba(fin) -> Dict[str, Any]:
    header = next(fin)[1:-1].split("\t")

    if header != arriba_header:
        raise RuntimeError("Unexpected file header")

    for line in fin:
        yield arriba_to_json(header, line.strip("\n"))


def read_genes(fname: str) -> Set[str]:
    """Read the genes we want to report"""
    report_genes = set()
    with open(fname) as fin:
        for line in fin:
            report_genes.add(line.strip("\n"))
    return report_genes


def main(fusion_file: str, report_genes_file: str) -> None:
    with open(fusion_file) as fin:
        # We want to keep the fusions in the same order as the input file
        fusions = list()

        # If no filter list was specified
        if report_genes_file is None:
            fusions = list(parse_arriba(fin))

        else:
            # Get the genes we want to report
            report_genes = read_genes(report_genes_file)

            for record in parse_arriba(fin):
                if record["gene1"] in report_genes or record["gene2"] in report_genes:
                    fusions.append(record)
        print(json.dumps(fusions, indent=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fusions")
    parser.add_argument("--report-genes", required=False)

    args = parser.parse_args()

    main(args.fusions, args.report_genes)
