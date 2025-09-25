from typing import Any, Iterator, Dict, Tuple
import functools
from collections import namedtuple
from mutalyzer_hgvs_parser import to_model  # type: ignore

# Type for the frequencies entry from VEP
FrequenciesType = Dict[str, Dict[str, float]]

# Tuple to store locations of HGVS points
Location = namedtuple("Location", ["downstream", "position", "offset"])

# Tuple to store a region
Region = namedtuple("Region", ["start", "end"])

class VEP(dict[str, Any]):
    """Class to work with VEP objects"""
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

    def filter_criteria(self, criteria: list["Criterion"]) -> None:
        filtered = list()
        for tc in self.get("transcript_consequences", list()):
            hgvsc = tc.get("hgvsc")
            if hgvsc is None:
                continue

            variant = Variant(hgvsc, tc["consequence_terms"])
            for crit in criteria:
                if crit.match(variant):
                    filtered.append(tc)
                    break
        self["transcript_consequences"] = filtered
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
        for term in self.severity:
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
        - There can be multiple colocated_variants entries with a "frequencies"
            section, as long as they are all identical
        """
        frequencies: FrequenciesType = dict()
        if "colocated_variants" not in self:
            return frequencies

        for var in self["colocated_variants"]:
            if "frequencies" in var:
                if not frequencies:
                    frequencies = var["frequencies"]
                elif var["frequencies"] != frequencies:
                    location = self.location
                    msg = f"Multiple colocated variants with 'frequencies' entry encountered on {location}"
                    raise RuntimeError(msg)

        if len(frequencies) > 1:
            location = self.location
            msg = f"'frequencies' entry with multiple keys encountered on {location}"
            raise RuntimeError(msg)

        return frequencies

    @classmethod
    def _extract_population(
        cls, population: str, frequencies: FrequenciesType
    ) -> float:
        """
        Extract the frequency from the specified population from the VEP record
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

    def population_frequency(self, population: str) -> float:
        """
        Return the specified population frequency from the VEP record

        Returns 0 if the data is missing
        """
        frequencies = self._extract_frequencies()
        return self._extract_population(population, frequencies)

    def above_population_threshold(self, population: str, threshold: float) -> bool:
        """
        Determine if the VEP population frequency is above the specified threshold
        """
        return self.population_frequency(population) > threshold

    @property
    def location(self) -> str:
        """
        Return a representation the location of the VEP record on the genome
        """
        input = self.get("input")
        if input is None:
            return "unknown location"

        chrom, pos = input.split("\t")[:2]
        return f"{chrom}:{pos}"


class Variant:
    def __init__(self, hgvs: str, consequences: list[str]):
        self.hgvs = hgvs
        self.consequences = consequences

    @classmethod
    def from_VEP(cls, vep: VEP) -> Iterator["Variant"]:
        """Yield all Variants from a VEP record"""
        for transcript in vep.get("transcript_consequences", list()):
            consequences = transcript["consequence_terms"]
            for hgvs in ["hgvsg", "hgvsc", "hgvsp"]:
                if hgvs in transcript:
                    yield Variant(transcript[hgvs], consequences)

    def __str__(self) -> str:
        return f"Variant({self.hgvs}, {self.consequences})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError
        return str(self) == str(other)


class Criterion:
    def __init__(
        self,
        identifier: str,
        coordinate: str = "c",
        consequence: str | None = None,
        start: str | None = None,
        end: str | None = None,
    ) -> None:
        self.identifier = identifier
        self.coordinate = coordinate
        self.consequence = consequence

        # Define the type of start and end
        self.start: Location | None
        self.end: Location | None

        if start is not None:
            # if end is not specified, set end equal to start
            if end is None:
                end = start
            # Make a fake hgvs description to determine the start and end Location
            hgvs = f"{identifier}:c.{start}_{end}"
            self.start, self.end = get_position(hgvs)
        else:
            self.start = None
            self.end = None

    def __str__(self) -> str:
        return (
            f"Criterion(identifier={self.identifier}, "
            f"coordinate={self.coordinate}, "
            f"consequence={self.consequence}, "
            f"start={self.start}, "
            f"end={self.end})"
        )

    def __repr__(self) -> str:
        return str(self)

    def match(self, variant: Variant) -> bool:
        return (self.match_id(variant) and
                self.match_coordinate(variant) and
                self.match_consequence(variant) and
                self.match_region(variant)
                )

    def match_id(self, variant: Variant) -> bool:
        """Determine if the identifier matches the variant

        Raises a runtime error if only the version doesn't match, since this
        could indicate the use of incompatible reference seqeunces
        """

        def split_version(variant_id: str) -> Tuple[str, str]:
            """Split the version from the identifier

            Set version to 0 if there is no version
            """
            try:
                id_, version = variant_id.split(".")
            # If there is no version (e.g. chr5), set version to 0
            except ValueError:
                id_ = variant_id
                version = "0"

            return id_, version

        # Get the variant id
        variant_id = variant.hgvs.split(":")[0]
        criterion_id = self.identifier

        id1, v1 = split_version(variant_id)
        id2, v2 = split_version(criterion_id)

        if id1 != id2:
            return False

        if v1 == v2:
            return True
        else:
            msg = f"Version mismatch between {variant} and {self}"
            raise RuntimeError(msg)

    def match_coordinate(self, variant: Variant) -> bool:
        variant_coordinate = variant.hgvs.split(":")[1].split(".")[0]
        return variant_coordinate == self.coordinate

    def match_consequence(self, variant: Variant) -> bool:
        if self.consequence is None:
            return True
        return self.consequence in variant.consequences

    def match_region(self, variant: Variant) -> bool:
        if self.start is None:
            return True
        var_region = get_position(variant.hgvs)
        crit_region = Region(self.start, self.end)

        return region_overlap(var_region, crit_region)

@functools.lru_cache(maxsize=1000)
def get_position(hgvs: str) -> Region:
    # Get the variant part of the hgvs description
    coordinate, var = hgvs.split(":")[1].split(".")

    if coordinate == "p":
        model = to_model(var, "p_variant")
    else:
        model = to_model(var, "variant")

    location = model["location"]
    if location["type"] == "point":
        start = point_to_tuple(location)
        end = start
    elif location["type"] == "range":
        start = point_to_tuple(location["start"])
        end = point_to_tuple(location["end"])

        # make sure "end" is downstream of "start"
        if start > end:
            s, e = start, end
            start = e
            end = s
    else:
        raise RuntimeError

    return Region(start, end)


def point_to_tuple(point: dict[str, Any]) -> Location:
    """Convert a 'point' dictionary from parse_hgvs into a Location"""

    # Default values
    downstream = 0
    position = 0
    offset = 0

    position = point["position"]
    if "outside_cds" in point:
        if point["outside_cds"] == "upstream":
            position *= -1
        elif point["outside_cds"] == "downstream":
            downstream = 1
        else:
            raise RuntimeError
    if "offset" in point:
        offset = point["offset"]["value"]

    return Location(downstream, position, offset)


def region_overlap(region1: Region, region2: Region) -> bool:
    return any(
        [
            region1[0] >= region2[0] and region1[0] <= region2[1],
            region1[1] >= region2[0] and region1[1] <= region2[1],
            region1[0] < region2[0] and region1[1] > region2[1],
        ]
    )
