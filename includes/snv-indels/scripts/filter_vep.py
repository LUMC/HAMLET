#!/usr/bin/env python3

import argparse
import json

from functools import partial


# From most to least severe, taken from the ensembl website
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
severity = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]


def gene_of_interest(cons, genes):
    """Is a VEP consequence applicable to a gene of interest"""
    return cons["gene_id"] in genes


def transcript_of_interest(cons, transcripts):
    """Is a VEP consequence applicable to a transcript of interest"""
    return cons["transcript_id"] in transcripts


def consequence_of_interest(cons, genes, transcripts, impact=None):
    """Is a VEP consequence of interest for both gene and transcript"""
    return all(
        [
            gene_of_interest(cons, genes),
            transcript_of_interest(cons, transcripts),
            # Only filter on impact if one was specified
            cons["impact"] == impact if impact else True,
        ]
    )


def consequences_of_interest(cons, genes, transcripts, impact=None):
    """Filter consequences to only those of interest"""
    return [
        con for con in cons if consequence_of_interest(con, genes, transcripts, impact)
    ]


def update_most_severe(vep):
    """The most severe consequence for all genes and transcript is stored at
    the top level of the VEP object. After filtering the transcript
    consequences, we need to update this field
    """
    # Gather all consequences
    cons = set()
    for consequence in vep["transcript_consequences"]:
        cons.update(consequence["consequence_terms"])

    # Set the first matching consequence (they are sorted based on severity)
    for term in severity:
        if term in cons:
            vep["most_severe_consequence"] = term
            break


def vep_of_interest(vep, genes, transcripts, impact=None):
    """Return a new VEP object which only contains consequences of interest

    1. Restrict the transcript_consequences to only include genes/transcripts
    of interest, and the impact, if specified.
    2. Rewrite the most_severe_consequence based on the remaining transcripts.
    """
    # Copy the VEP object
    new_vep = vep.copy()

    # Extract the transcript consequences
    cons = vep.get("transcript_consequences", list())

    # Extract only the transcripts of interest
    cons_int = consequences_of_interest(cons, genes, transcripts, impact)

    # Replace the transcripts
    new_vep["transcript_consequences"] = cons_int

    # Set the most severe consequence
    update_most_severe(new_vep)

    return new_vep


def read_goi_file(fname):
    goi = set()
    toi = set()
    with open(fname) as fin:
        header = next(fin).strip("\n").split("\t")
        for line in fin:
            d = {k: v for k, v in zip(header, line.strip("\n").split("\t"))}
            # There is only a single goi
            goi.add(d["GOI_ID"])
            # There can be multiple toi's
            toi.update(d["TOI_IDS"].split(","))
    return goi, toi


def parse_vep_json(vep_file):
    """Parse the VEP 'json' output file, each line contains a JSON entry"""
    with open(vep_file) as fin:
        for line in fin:
            yield json.loads(line)


def main(vep_file, goi_file, impact):
    # Get genes and transcripts of interest
    goi, toi = read_goi_file(goi_file)

    for variant in parse_vep_json(vep_file):
        vep = vep_of_interest(variant, goi, toi, impact)
        # If there is no consequence of interest
        if not vep["transcript_consequences"]:
            continue
        print(json.dumps(vep, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract genes (and transcript) of interest from VEP output"
    )

    parser.add_argument("vep", help="VEP json output file")
    parser.add_argument("goi", help="Genes of interest")
    parser.add_argument(
        "--impact", choices=["HIGH", "MODERATE", "MODIFIER"], default=None
    )

    args = parser.parse_args()

    main(args.vep, args.goi, args.impact)
