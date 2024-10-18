#!/usr/bin/env python3

import argparse
import gzip
import json

from typing import Any, Dict, Set, Tuple, Iterator

# Type for the frequencies entry from VEP
FrequenciesType = Dict[str, Dict[str, float]]
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


class VEP(dict[str, Any]):
    """Class to work with VEP objects"""

    def filter_transcript_id(self, transcripts: Set[str]) -> None:
        """Filter transcript consequences by transcript_id"""
        tc = self.get("transcript_consequences", list())
        tc = [x for x in tc if x["transcript_id"] in transcripts]
        self["transcript_consequences"] = tc
        self.update_most_severe()

    def filter_consequence_term(self, consequences: Set[str]) -> None:
        """Filter transcript consequences by consequence_term"""
        # If there are not consequences to filter on, we do nothing
        if not consequences:
            return
        tc = self["transcript_consequences"]
        tc = [x for x in tc if not set(x["consequence_terms"]).isdisjoint(consequences)]
        self["transcript_consequences"] = tc
        self.update_most_severe()

    def filter_hgvsc_blacklist(self, blacklist: Set[str]) -> None:
        """Filter transcript consequences by hgvsc field"""
        if not blacklist:
            return
        tc = self["transcript_consequences"]
        tc = [x for x in tc if not x.get("hgvsc") in blacklist]
        self["transcript_consequences"] = tc
        self.update_most_severe()

    def update_most_severe(self) -> None:
        """The most severe consequence for all genes and transcript is stored
        at the top level of the VEP object. After filtering the transcript
        consequences, we need to update this field
        """
        # Gather all consequences
        cons: Set[str] = set()
        for consequence in self["transcript_consequences"]:
            cons.update(consequence["consequence_terms"])

        # Set the most severe consequence (the list of consequences is sorted
        # by severity)
        for term in severity:
            if term in cons:
                self["most_severe_consequence"] = term
                break

    def _extract_frequencies(self) -> FrequenciesType:
        """
        Extract the population allele frequencies from a VEP record

        The structure of the VEP record is quite complicated. I've looked at
        191k VEP entries from HAMLET, and found the following:
        - Each record has zero or more "colocated_variants" entries
        - Each colocated_variant entry can have a "frequencies" section
        - Only one of the colocated_variant entries can have a "frequencies"
            section
        """
        frequencies: FrequenciesType = dict()
        if "colocated_variants" not in self:
            return frequencies

        nr_freq_entries = 0
        for var in self["colocated_variants"]:
            if "frequencies" in var:
                nr_freq_entries += 1
                frequencies = var["frequencies"]

        if nr_freq_entries > 1:
            msg = "Multiple colocated variants with 'frequencies' entry encountered"
            raise RuntimeError(msg)

        if len(frequencies) > 1:
            msg = "'frequencies' entry from VEP should only contain a single key"
            raise RuntimeError(msg)

        return frequencies

    @classmethod
    def _extract_population(
        self, population: str, frequencies: FrequenciesType
    ) -> float:
        """
        Extract the frequency from the specified population from the VEP record

        TODO check that there is only a single key in frequencies
        TODO make private function
        TODO write top_level function taht uses extract_frequencies and extract
             population on self
        """
        # If frequencies is empty, return 0
        if not frequencies:
            return 0

        if len(frequencies) > 1:
            msg = "'frequencies' entry from VEP should only contain a single key"
            raise RuntimeError(msg)

        # Get the single key from frequencies
        key = next(iter(frequencies))
        return frequencies[key].get(population, 0)

    def extract_population_frequency(self, population: str) -> float:
        """
        Extract the specified population frequency from the VEP record
        """
        frequencies = self._extract_frequencies()
        return self._extract_population(population, frequencies)


def read_goi_file(fname: str) -> Tuple[Set[str], Set[str]]:
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


def parse_vep_json(vep_file: str) -> Iterator[VEP]:
    """Parse the VEP 'json' output file, each line contains a JSON entry"""
    with gzip.open(vep_file, "rt") as fin:
        for line in fin:
            yield VEP(json.loads(line))


def get_hotspot(fname: str) -> Set[str]:
    hotspots = set()
    with open(fname) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            else:
                hotspots.add(line.strip("\n"))
    return hotspots


def get_blacklist(fname: str) -> Set[str]:
    blacklist = set()
    with open(fname) as fin:
        for line in fin:
            blacklist.add(line.strip("\n"))
    return blacklist


def main(
    vep_file: str,
    goi_file: str,
    consequences: Set[str],
    hotspot_file: str,
    blacklist_file: str,
) -> None:
    # Get genes and transcripts of interest
    goi, toi = read_goi_file(goi_file)

    # Get the hotspot mutations
    hotspot = get_hotspot(hotspot_file) if hotspot_file else set()

    # Get the blacklisted variant
    blacklist = get_blacklist(blacklist_file) if blacklist_file else set()

    for vep in parse_vep_json(vep_file):
        # Filter on transcript of interest
        vep.filter_transcript_id(toi)
        # Filter on consequences of interest
        vep.filter_consequence_term(consequences)
        # Add is_in_hotspot field
        vep["is_in_hotspot"] = vep["input"] in hotspot
        # Filter on hgvsc blacklist
        vep.filter_hgvsc_blacklist(blacklist)

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
    parser.add_argument("--consequences", nargs="*", type=str, default=list())
    parser.add_argument("--hotspot", help="VCF file with hotspot variants")
    parser.add_argument(
        "--blacklist",
        help=(
            "File with blacklisted variants, one per line. Format "
            "should match the 'hgvsc' field of VEP"
        ),
    )

    args = parser.parse_args()

    main(args.vep, args.goi, args.consequences, args.hotspot, args.blacklist)
