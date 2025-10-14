from typing import Any, Iterator, Dict, Sequence, Tuple, Set
import functools
from collections import namedtuple
from mutalyzer_hgvs_parser import to_model  # type: ignore
from itertools import zip_longest
from collections import OrderedDict

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

    def filter_criteria(self, criteria: Sequence["Criterion"]) -> None:
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

    def filter_annotate_transcripts(self, known_variants: dict[str, str], criteria: dict["Criterion", str], ) -> None:
        """Filter and annotate the transcripts"""

        filtered = list()
        for transcript in self.get("transcript_consequences", list()):
            # Make sure hgvsc is defined
            hgvsc = transcript.get("hgvsc")

            # If there is no HGVSC, we do nothing
            if hgvsc is None:
                continue

            # If hgvsc is a known variant, we update the annotation
            if hgvsc in known_variants:
                transcript["annotation"] = known_variants[hgvsc]
                filtered.append(transcript)
                continue

            # Otherwise, we check the criteria
            variant = Variant(hgvsc, transcript["consequence_terms"])
            for crit, annotation in criteria.items():
                if crit.match(variant):
                    transcript["annotation"] = annotation
                    filtered.append(transcript)
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

    @staticmethod
    def _location_size(location: dict[str, Any]) -> int:
        """Determine the size of a 'location' dict from mutalyzer"""
        if location["type"] == "point":
            size = 1
        elif location["type"] == "range":
            start = int(location["start"]["position"])
            end = int(location["end"]["position"])
            size = end - start + 1
        else:
            raise NotImplementedError

        return size

    def _outside_cds(self, location: dict[str, Any]) -> bool:
        """Determine if a location object from Mutalyzer is outside the CDS"""
        if "outside_cds" in location:
            return True
        if "offset" in location:
            return True
        if "start" in location and "outside_cds" in location["start"]:
            return True
        if "end" in location and "outside_cds" in location["end"]:
            return True
        if "start" in location and "offset" in location["start"]:
            return True
        if "end" in location and "offset" in location["end"]:
            return True

        return False

    def size(self) -> int:
        """Return the size of the variant, inserted-deleted"""
        model = to_model(self.hgvs)
        variants = model["variants"]

        # No variant
        if not model["variants"]:
            return 0

        if len(variants) > 1:
            raise NotImplementedError

        variant = variants[0]

        # SNP
        if "deleted" in variant and "inserted" in variant:
            d = variant["deleted"][0]["sequence"]
            i = variant["inserted"][0]["sequence"]
            return len(i) - len(d)

        # A deletion
        elif variant["type"] == "deletion":
            return -self._location_size(variant["location"])

        # An insertion
        elif variant["type"] == "insertion":
            return len(variant["inserted"][0]["sequence"])

        # A delins
        elif variant["type"] == "deletion_insertion":
            del_size = self._location_size(variant["location"])
            ins_size = len(variant["inserted"][0]["sequence"])
            return ins_size - del_size
        elif variant["type"] == "inversion":
            return 0
        elif variant["type"] == "duplication":
            return self._location_size(variant["location"])
        else:
            raise NotImplementedError

    def frame(self) -> int:
        model = to_model(self.hgvs)
        # Determine if the coordinate system supports detecting the frame
        coordinate = model["coordinate_system"]
        if coordinate != "c":
            raise ValueError("Determining the frame is only supported for c. variants")

        variants = model["variants"]
        # No variant
        if not model["variants"]:
            return 0

        if len(variants) > 1:
            raise NotImplementedError

        variant = variants[0]
        # Next, determine if the variant supports detecting the frame
        if self._outside_cds(variant["location"]):
            raise ValueError(
                "Determing the frame is only supported in the coding region"
            )

        return self.size() % 3


class Criterion:
    def __init__(
        self,
        identifier: str,
        coordinate: str = "c",
        consequence: str | None = None,
        start: str | None = None,
        end: str | None = None,
        frame: int | None = None,
    ) -> None:
        self.identifier = identifier
        self.coordinate = coordinate
        self.consequence = consequence

        # Define the type of start and end
        self.start: Location | None
        self.end: Location | None

        if end is not None and start is None:
            raise ValueError("Please specify a start and end")

        if start is not None:
            # if end is not specified, set end equal to start
            if end is None:
                end = start

            # Make a fake hgvs description to determine the start and end Location
            hgvs = f"{identifier}:c.{start}_{end}"
            self.start, self.end = get_position(hgvs)
            if (
                self.start is not None
                and self.end is not None
                and self.start > self.end
            ):
                raise ValueError(f"{start=} cannot be after {end=}")

        else:
            self.start = None
            self.end = None

        self.frame = frame

    def __str__(self) -> str:
        return (
            f"Criterion(identifier={self.identifier}, "
            f"coordinate={self.coordinate}, "
            f"consequence={self.consequence}, "
            f"start={self.start}, "
            f"end={self.end}, "
            f"frame={self.frame})"
        )

    def __repr__(self) -> str:
        return str(self)

    def _contains_identifier(self, other: "Criterion") -> bool:
        if self.identifier is None:
            return True

        # Check that the identifier versions match
        id1, v1 = self.split_version(self.identifier)
        id2, v2 = self.split_version(other.identifier)

        # Identifiers do not match
        if id1 != id2:
            return False

        # Version numbers do not match
        if v1 != v2:
            msg = f"Version mismatch between {self} and {other}"
            raise ValueError(msg)

        return self.identifier == other.identifier

    def _contains_consequence(self, other: "Criterion") -> bool:
        if self.consequence is None:
            return True
        else:
            return self.consequence == other.consequence

    def _contains_region(self, other: "Criterion") -> bool:
        r1 = Region(self.start, self.end)
        r2 = Region(other.start, other.end)
        if r1 == (None, None):
            return True
        else:
            return region_contains(r1, r2)

    def _contains_frame(self, other: "Criterion") -> bool:
        if self.frame is None:
            return True
        else:
            return self.frame == other.frame

    def contains(self, other: "Criterion") -> bool:
        """Determine if other falls within self

        In this context, this means that any Variant that matches other will
        also match self.
        """
        if not isinstance(other, Criterion):
            raise NotImplementedError


        return (
            self._contains_identifier(other)
            and self.coordinate == other.coordinate
            and self._contains_consequence(other)
            and self._contains_region(other)
            and self._contains_frame(other)
        )

    def match(self, variant: Variant) -> bool:
        return (
            self.match_id(variant)
            and self.match_coordinate(variant)
            and self.match_consequence(variant)
            and self.match_region(variant)
            and self.match_frame(variant)
        )

    def split_version(self, identifier: str) -> Tuple[str, str]:
        """Split the version from the identifier

        Set version to 0 if there is no version
        """
        try:
            id_, version = identifier.split(".")
        # If there is no version (e.g. chr5), set version to 0
        except ValueError:
            id_ = identifier
            version = "0"

        return id_, version

    def match_id(self, variant: Variant) -> bool:
        """Determine if the identifier matches the variant

        Raises a runtime error if only the version doesn't match, since this
        could indicate the use of incompatible reference seqeunces
        """

        # Get the variant id
        variant_id = variant.hgvs.split(":")[0]
        criterion_id = self.identifier

        id1, v1 = self.split_version(variant_id)
        id2, v2 = self.split_version(criterion_id)

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

    def match_frame(self, variant: Variant) -> bool:
        if self.frame is None:
            return True
        else:
            return variant.frame() == self.frame


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

def region_contains(region1: Region, region2: Region) -> bool:
    """Determine of region1 contains region2

        Note that both start and end of both regions can be None
    """
    # Determine if the start positions are set for region1 and region2
    if region1.start is None and region2.start is not None:
        return False
    elif region1.start is not None and region2.start is None:
        return False

    # Determine if the end positions are set for self and other
    if region1.end is None and region2.end is not None:
        return False
    elif region1.end is not None and region2.end is None:
        return False

    # If both regions are simply None
    if region1 == (None, None) and region2 == (None, None):
        return True

    # Here, we know both regions do not contain any None
    return bool(region1.start <= region2.start and region1.end >= region2.end)

def read_criteria_file(criteria_file: str) -> OrderedDict[Criterion, str]:
    """Read the criteria and annotations from the criteria file

    If the annotations column does not exist, all criteria will get an empty
    string as annotation
    """

    annotations: OrderedDict[Criterion, str] = OrderedDict()

    header = None
    with open(criteria_file) as fin:
        for line in fin:
            if line.startswith("#"):
                continue

            spline = line.strip("\n").split("\t")
            if header is None:
                header = spline
                continue

            # Read into dict, convert '' to None
            d = {k: v if v else None for k, v in zip_longest(header, spline)}

            # Check that at least the transcript id is set
            transcript_id = d.get("transcript_id")
            assert transcript_id is not None

            # Get the frame
            frame = d.get("frame")
            if frame is not None:
                frame = int(frame)

            c = Criterion(
                identifier=transcript_id,
                coordinate="c",
                consequence=d["consequence"],
                start=d["start"],
                end=d["end"],
                frame=frame,
            )
            annotation = str(d.get("annotation", ""))

            annotations[c] = annotation

    return annotations

def read_known_variants(fname: str) -> dict[str, str]:
    known_variants = dict()
    header = None
    with open(fname) as fin:
        for line in fin:
            if line.startswith("#"):
                continue

            spline = line.strip("\n").split("\t")
            if header is None:
                header = spline
                continue

            # Read into dict, convert '' to None
            d = {k: v if v else None for k, v in zip_longest(header, spline)}

            # Check that the variant column is filled
            variant = d.get("variant")
            assert variant is not None

            # Check that the annotation column is filled
            annotation = d.get("annotation")
            assert annotation is not None

            known_variants[variant] = annotation
    return known_variants
