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


def consequence_of_interest(cons, genes, transcripts):
    """Is a VEP consequence of interest for both gene and transcript"""
    return gene_of_interest(cons, genes) and transcript_of_interest(cons, transcripts)


def consequences_of_interest(cons, genes, transcripts):
    """Filter consequences to only those of interest"""
    return [con for con in cons if consequence_of_interest(con, genes, transcripts)]


def vep_of_interest(vep, genes, transcripts):
    """Return a new VEP object which only contains consequences of interest

    1. Restrict the transcript_consequences to only include genes/transcripts
    of interest.
    2. Rewrite the most_severe_consequence based on the remaining transcripts.
    """
    # Copy the VEP object
    new_vep = vep.copy()

    # Extract the transcript consequences
    cons = vep["transcript_consequences"]

    # Extract only the transcripts of interest
    cons_int = consequences_of_interest(cons, genes, transcripts)

    # Replace the transcripts
    new_vep["transcript_consequences"] = cons_int

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


def main(vep_file, goi_file):
    # Get genes and transripts of interest
    goi, toi = read_goi_file(goi_file)

    for variant in parse_vep_json(vep_file):
        # Store the consequences of interest
        cons = list()
        for consequence in variant.get("transcript_consequences", list()):
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract genes (and transcript) of interest from VEP output"
    )

    parser.add_argument("vep", help="VEP json output file")
    parser.add_argument("goi", help="Genes of interest")

    args = parser.parse_args()

    main(args.vep, args.goi)