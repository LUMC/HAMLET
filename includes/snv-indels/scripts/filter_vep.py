#!/usr/bin/env python3

import argparse
import json

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


class VEP(dict):
    """Class to work with VEP objects"""
    def filter_transcript_id(self, transcripts):
        """Filter transcript consequences by transcript_id"""
        tc = self.get("transcript_consequences", list())
        tc = [x for x in tc if x["transcript_id"] in transcripts]
        self["transcript_consequences"] = tc
        self.update_most_severe()


    def filter_consequence_term(self, consequences):
        """Filter transcript consequences by consequence_term"""
        tc = self["transcript_consequences"]
        tc = [x for x in tc if not set(x["consequence_terms"]).isdisjoint(consequences)]
        self["transcript_consequences"] = tc
        self.update_most_severe()


    def update_most_severe(self):
        """The most severe consequence for all genes and transcript is stored at
        the top level of the VEP object. After filtering the transcript
        consequences, we need to update this field
        """
        # Gather all consequences
        cons = set()
        for consequence in self["transcript_consequences"]:
            cons.update(consequence["consequence_terms"])

        # Set the first matching consequence (they are sorted based on severity)
        for term in severity:
            if term in cons:
                self["most_severe_consequence"] = term
                break


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
            yield VEP(json.loads(line))


def get_hotspot(fname):
    hotspots = set()
    with open(fname) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            else:
                hotspots.add(line.strip('\n'))
    return hotspots

def main(vep_file, goi_file, consequences, hotspot_file):
    # Get genes and transcripts of interest
    goi, toi = read_goi_file(goi_file)

    # Get the hotspot mutations
    if hotspot_file:
        hotspot = get_hotspot(hotspot_file)
    else:
        hotspot = set()

    for vep in parse_vep_json(vep_file):
        # Filter on transcript of interest
        vep.filter_transcript_id(toi)
        # Filter on consequences of interest
        vep.filter_consequence_term(consequences)
        # Add is_in_hotspot field
        vep["is_in_hotspot"] = vep["input"] in hotspot

        # If there is no consequence of interest left
        if not vep["transcript_consequences"]:
            continue
        print(json.dumps(vep, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract genes (and transcript) of interest from VEP output"
    )


    parser.add_argument("vep", help="VEP json output file")
    parser.add_argument("goi", help="Genes of interest")
    parser.add_argument( "--consequences", nargs='*',
            type=str, default=list())
    parser.add_argument("--hotspot", help="VCF file with hotspot variants")

    args = parser.parse_args()

    main(args.vep, args.goi, args.consequences, args.hotspot)
