#!/usr/bin/env python3

# Authors: Anne van der Grinten, Redmar van den Berg

from typing import Sequence
from typing import Any, Dict, List, Tuple
import json
import pytest

from utils import (
    VEP,
    FrequenciesType,
    Criterion,
    Variant,
    Location,
    Region,
    region_overlap,
    get_position,
)


class TestFilterAnnotate:
    """
    Tests for the filtering/annotation of transcripts for variants, based on Criteria and known variants


    The test setup is as follows:
    - We have a single VEP record, which contains 3 transcript consequences
    - We have 3 filter/annotation Criteria, which match each of the 3 transcript
      consequences
    - We have 3 known variants, which match each of the 3 transcript consequences

    The tests combine the Criteria with the known_variants to ensure that the
    annotations and filtering works correctly, with the expected priority.
    """

    @pytest.fixture
    def vep(self) -> VEP:
        minimal_vep_record = {
            "transcript_consequences": [
                {
                    "hgvsc": "ENST1.1:c.100del",
                    "consequence_terms": ["frameshift"],
                },
                {
                    "hgvsc": "ENST2.1:c.100A>T",
                    "consequence_terms": ["missense", "splice site"],
                },
                {
                    "hgvsc": "ENST3.1:c.100+100A>T",
                    "consequence_terms": [],
                },
            ]
        }
        return VEP(minimal_vep_record)

    def criteria(self):
        return [
            Criterion("ENST1.1", consequence="frameshift"),
            Criterion("ENST2.1", start="0", end="200", consequence="splice site"),
            Criterion("ENST3.1", start="100+20", end="101-20"),
        ]

    @pytest.fixture
    def inclusion_criteria(self) -> Sequence[Criterion]:
        return self.criteria()

    @pytest.fixture
    def annotation_criteria(self) -> Sequence[Criterion]:
        """The annotation criteria are the same as the inclusion criteria"""
        return self.criteria()

    @pytest.fixture
    def known_variants(self) -> Sequence[str]:
        return ["ENST1.1:c.100del", "ENST2.1:c.100A>T", "ENST3.1:c.100+100A>T"]

    def test_filter_annotate_variants_empty(self, vep: VEP) -> None:
        vep.filter_annotate_transcripts(dict(), list(), dict())
        assert not vep["transcript_consequences"]

    @pytest.mark.parametrize(
        "inclusion_index, var_index, annotation_index, expected_annotations",
        [
            # No inclusion means all transcript consequences are filtered out
            ("", "012", "012", []),
            # If all three inclusion and annotation criteria are set
            ("012", "", "012", ["crit", "crit", "crit"]),
            # If two inclusion and three annotation criteria are set
            ("01", "", "012", ["crit", "crit"]),
            # If one inclusion and three annotation criteria are set
            ("0", "", "012", ["crit"]),
            # If all three inclusion and known variants are set. Variants have
            # priorty over annotations
            ("012", "012", "012", ["var", "var", "var"]),
            # Idem, but with two inclusion criteria set
            ("01", "012", "012", ["var", "var"]),
            # Idem, but with one inclusion criteria set
            ("0", "012", "012", ["var"]),
            # Idem, but without inclusion criteria set
            ("", "012", "012", []),
            # The second variant is not in known variants
            ("012", "02", "012", ["var", "crit", "var"]),
            ("012", "02", "1", ["var", "crit", "var"]),
            ("012", "02", "", ["var", "", "var"]),
            # Everything included, but nothing annotated
            # Note that the annotation exist, but is empty
            ("012", "", "", ["", "", ""]),
        ],
    )
    def test_filter_annotate_variants(
        self,
        vep: VEP,
        inclusion_criteria: Sequence[Criterion],
        known_variants: Sequence[str],
        annotation_criteria: Sequence[Criterion],
        inclusion_index: str,
        var_index: str,
        annotation_index: str,
        expected_annotations: list[str],
    ) -> None:
        # Create the list for the inclusion criteria
        inclusion = list()
        for index in inclusion_index:
            inclusion.append(inclusion_criteria[int(index)])

        # Create the dict for the known variants
        variants = dict()
        for index in var_index:
            variants[known_variants[int(index)]] = "var"

        # Create the dict for the annotation criteria
        annotation = dict()
        for index in annotation_index:
            annotation[annotation_criteria[int(index)]] = "crit"

        vep.filter_annotate_transcripts(inclusion, variants, annotation)

        annotations = [
            transcript["annotation"] for transcript in vep["transcript_consequences"]
        ]

        assert annotations == expected_annotations


@pytest.fixture
def vep() -> VEP:
    with open("test/data/output/v2/vep.json") as fin:
        js = json.load(fin)
    return VEP(js[0])


CRITERIA = [
    # Criteria, readout of gene names
    ([Criterion("ENTS0123.1")], ["gene1"]),
    ([Criterion("ENTS0123.1"), Criterion("ENTS0124.1")], ["gene1", "gene2"]),
    ([Criterion("ENST999.9")], []),
    ([Criterion("ENTS0123.1", consequence="transcript_ablation")], ["gene1"]),
    ([Criterion("ENTS0125.1", consequence="stop_gained")], ["gene3"]),
    ([Criterion("ENTS0125.1", consequence="inframe_insertion")], ["gene3"]),
]


@pytest.mark.parametrize(["criteria", "readout"], CRITERIA)
def test_filter_criteria(
    vep: VEP, criteria: list[Criterion], readout: list[str]
) -> None:
    """Test filtering by consequence_term"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_annotate_transcripts(criteria, dict(), dict())
    genes = [tc["gene_id"] for tc in vep["transcript_consequences"]]
    assert genes == readout


def test_filter_criteria_no_match(vep: VEP) -> None:
    """Test that we get an empty list if no criteria matches"""
    vep.filter_annotate_transcripts([Criterion("no_such_transcript")], dict(), dict())
    assert not vep["transcript_consequences"]


def test_filter_no_criteria(vep: VEP) -> None:
    """Check that we filter nothing if there are no criteria"""
    vep.filter_annotate_transcripts([], dict(), dict())
    assert not vep["transcript_consequences"]


def test_vep_of_interest_one_transcript(vep: VEP) -> None:
    """Test rewriting the VEP object when there is a single transcript of
    interest"""
    c = Criterion("ENTS0125.1")

    # Restrict to transcript of interest
    vep.filter_annotate_transcripts([c], dict(), dict())

    # Get the consequences after rewriting the VEP obejct
    cons = vep["transcript_consequences"]

    # Test that only gene3 is left
    assert len(cons) == 1
    assert cons[0]["gene_id"] == "gene3"
    assert cons[0]["transcript_id"] == "transcript3"
    assert cons[0]["consequence_terms"] == ["inframe_insertion", "stop_gained"]

    # Test that we set the most severe consequence for transcript3
    assert vep["most_severe_consequence"] == "stop_gained"


def test_update_most_severe(vep: VEP) -> None:
    vep.update_most_severe()
    assert vep["most_severe_consequence"] == "transcript_ablation"

    # Remove gene1
    vep["transcript_consequences"].pop(0)
    vep.update_most_severe()
    assert vep["most_severe_consequence"] == "splice_acceptor_variant"

    # Remove gene2
    vep["transcript_consequences"].pop(0)
    vep.update_most_severe()
    assert vep["most_severe_consequence"] == "stop_gained"


def test_no_transcript_consequence() -> None:
    """Test there are no crashes on empty VEP objects"""
    # Create an empty VEP object
    empty = VEP(dict())
    c = Criterion("ENTS0123.1")
    empty.filter_annotate_transcripts([c], dict(), dict())
    empty.update_most_severe()


# fmt: off
COLOCATED_VARIANTS :List[Tuple[Any,Any]]= [
    # No "colocated_variants" entry
    (dict(), dict()),
    # No "frequencies" entry in colocated variants
    ({"colocated_variants": []}, dict()),
    # A single frequencies entry
    (
        {
            "colocated_variants": [
                {
                    "frequencies": {"T": dict()}
                }
            ]
        },
        {"T": dict()}
    ),
    # A single frequencies entry in the second colocated variant
    (
        {
            "colocated_variants": [
                {},
                {
                    "frequencies": {"T": dict()}
                }
            ]
        },
        {"T": dict()}
    ),
    # A single frequencies entry in the first colocated variant
    (
        {
            "colocated_variants": [
                {
                    "frequencies": {"T": dict()}
                },
                {}
            ]
        },
        {"T": dict()}
    ),
    # Multiple identical frequencies entries
    (
        {
            "colocated_variants": [
                {
                    "frequencies": {"T": {"af": 0.5}}
                },
                {
                    "frequencies": {"T": {"af": 0.5}}
                }
            ]
        },
        {"T": {"af": 0.5}}
    ),
]
# fmt: on


@pytest.mark.parametrize(["vep", "frequencies"], COLOCATED_VARIANTS)
def test_extract_frequencies(vep: VEP, frequencies: FrequenciesType) -> None:
    V = VEP(vep)
    assert V._extract_frequencies() == frequencies


def test_error_extract_multiple_frequencies() -> None:
    """
    GIVEN a VEP record with multiple differing "frequencies" entries
    WHEN we extract the frequency values
    THEN we get an error
    """
    # fmt: off
    data: Dict[str, Any] = {
            "colocated_variants": [
                {
                    "frequencies": {"af": 0.9}
                },
                {
                    "frequencies": {"af": 0.8}
                }
            ]
    }
    # fmt: on
    V = VEP(data)
    with pytest.raises(RuntimeError):
        V._extract_frequencies()


def test_error_extract_frequencies_multiple_entries() -> None:
    """
    GIVEN a VEP record with a single "frequencies" entry, which contains multiple keys
    WHEN we extract the frequency values
    THEN we get an error
    """
    # fmt: off
    data: Dict[str, Any] = {"colocated_variants": [
            {"frequencies":
                {
                    "T": dict(),
                    "A": dict()
                }
            }
        ]
    }
    # fmt: on
    V = VEP(data)
    with pytest.raises(RuntimeError):
        V._extract_frequencies()


POPULATION = [
    # If frequencies is empty, return 0
    ({}, 0),
    # If the specified population is missing, return 0
    ({"T": {}}, 0),
    # Otherwise, return the requested population frequency
    ({"T": {"gnomAD": 0.5}}, 0.5),
]


@pytest.mark.parametrize(["frequencies", "expected"], POPULATION)
def test_extract_population(frequencies: FrequenciesType, expected: float) -> None:
    """
    GIVEN a VEP frequencies record
    WHEN we extract the specified population
    THEN we should get the corresponding population frequency
    """
    V = VEP()
    assert V._extract_population("gnomAD", frequencies) == expected


def test_extract_population_frequency() -> None:
    """
    GIVEN a VEP record
    WHEN we extract the specified population frequency
    THEN we should get the corresponding value

    This is an end to end test which relies on
        VEP_.extract_frequencies
        VEP._extract_population
    """
    # fmt: off
    data = {
        "colocated_variants": [
            # An empty colocated variant
            dict(),
            # Another colocated variant with some random annotations
            {"some": "nonsense"},
            # The colocated variant which contains the frequencies
            {
                "frequencies": {
                    "T": {
                        "gnomAD": 0.42
                    }
                 }
            }
        ]
    }
    # fmt: on

    V = VEP(data)
    assert V.population_frequency("gnomAD") == 0.42


THRESHOLD = [
    ("gnomAD", 0.41999, True),
    ("gnomAD", 42, False),
    ("gnomAD", 43, False),
    # If the population is missing, assume the frequency is 0, which means
    # nothing is above the threshold
    ("Missing", 1, False),
]


@pytest.mark.parametrize(["population", "threshold", "expected"], THRESHOLD)
def test_above_population_threshold(
    population: str, threshold: float, expected: bool
) -> None:
    # fmt: off
    data = {
        "colocated_variants": [
            # An empty colocated variant
            dict(),
            # Another colocated variant with some random annotations
            {"some": "nonsense"},
            # The colocated variant which contains the frequencies
            {
                "frequencies": {
                    "T": {
                        "gnomAD": 0.42
                    }
                 }
            }
        ]
    }
    # fmt: on
    V = VEP(data)

    assert V.above_population_threshold(population, threshold) is expected


@pytest.fixture
def variant() -> Variant:
    return Variant("ENST123.5:c.10A>T", consequences=["missense", "inframe"])


@pytest.fixture
def minimal_vep() -> VEP:
    # fmt: off
    minimal_vep_record = {
        "transcript_consequences": [
            {
                "consequence_terms": ["inframe_insertion"],
                "hgvsg": "chr13:g.28034118_28034147dup",
                "hgvsc": "ENST00000241453.1:c.1772_1801dup",
                "hgvsp": "ENSP00000241453.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            },
            {
                "consequence_terms": ["inframe_insertion", "NMD_transcript_variant"],
                "hgvsg": "chr13:g.28034118_28034147dup",
                "hgvsc": "ENST00000380987.1:c.1772_1801dup",
                "hgvsp": "ENSP00000370374.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            },
        ]
    }
    # fmt: on

    return VEP(minimal_vep_record)


def test_create_criterion() -> None:
    C = Criterion("ENST123.5")
    assert C.identifier == "ENST123.5"


def test_create_criterion_Locations() -> None:
    """
    WHEN we initialize a criterion
    THEN start and end must be Location tuples
    """
    C = Criterion("ENST123.5", start="10", end="*20-5")

    assert C.start == Location(0, 10, 0)
    assert C.end == Location(1, 20, -5)


def test_create_variant(variant: Variant) -> None:
    assert variant.hgvs == "ENST123.5:c.10A>T"
    assert variant.consequences == ["missense", "inframe"]


def test_variant_match(variant: Variant) -> None:
    """not
    GIVEN a variant and a criterion with different transcripts
    WHEN we check for a match
    THEN we should not find a match
    """
    C = Criterion("ENST120.5")

    assert not C.match(variant)


def test_transcript_match(variant: Variant) -> None:
    """
    GIVEN two identical transcripts
    WHEN we check for a match
    THEN we should find a match
    """
    C = Criterion("ENST123.5")
    assert C.match_id(variant)


def test_transcript_no_match(variant: Variant) -> None:
    """
    GIVEN two different transcripts
    WHEN we check for a match
    THEN we should not find a match
    """
    C = Criterion("ENST120.5")
    assert not C.match_id(variant)


def test_transcript_version_no_match(variant: Variant) -> None:
    """
    GIVEN a version mismatch between two transcripts
    WHEN we check for a match
    THEN we should get a runtime error
    """
    C = Criterion("ENST123.4")
    with pytest.raises(RuntimeError):
        C.match_id(variant)


def test_transcript_id_no_version(variant: Variant) -> None:
    """
    GIVEN a Transcript and Criterion without version number
    WHEN we compare them
    THEN they should match
    """
    V = Variant("Chr12:g.10A>T", consequences=[])
    C = Criterion("Chr12", coordinate="g")

    assert C.match_id(V)


def test_transcript_version_no_match_error(variant: Variant) -> None:
    """
    GIVEN a version mismatch between criterion and variant
    WHEN we check for a match
    THEN we should get a runtime error
    """
    C = Criterion("ENST123.4")
    with pytest.raises(RuntimeError):
        C.match(variant)


CONS = [(None, True), ("missense", True), ("stop lost", False)]


@pytest.mark.parametrize("consequence, expected", CONS)
def test_match_consequence(
    consequence: str | None, expected: bool, variant: Variant
) -> None:
    """
    GIVEN a variant and a criterion with identical transcript and a consequence
    WHEN we check for a match
    THEN we should get the answer that is expected
    """
    C = Criterion("ENST123.5", "c", consequence)
    assert C.match(variant) == expected


# fmt: off
POSITION = [
    ("-20del",(
        Location(0,-20,0),
        Location(0,-20,0)
    )),
    ("-10+5del",(
        Location(0,-10,5),
        Location(0,-10,5)
    )),
    ("-9-5del",(
        Location(0,-9,-5),
        Location(0,-9,-5)
    )),
    ("1del",(
        Location(0,1,0),
        Location(0,1,0)
    )),
    ("10+5del",(
        Location(0,10,+5),
        Location(0,10,+5)
    )),
    ("11-5del",(
        Location(0,11,-5),
        Location(0,11,-5)
    )),
    ("*1del",(
        Location(1,1,0),
        Location(1,1,0)
    )),
    ("*5+5del",(
        Location(1,5,5),
        Location(1,5,5)
    )),
    ("*6-5del",(
        Location(1,6,-5),
        Location(1,6,-5)
    )),
    ("10_15del",(
        Location(0,10,0),
        Location(0,15,0)
    )),
]
# fmt: on


@pytest.mark.parametrize("mutation, expected_position", POSITION)
def test_variant_position_calculation(
    mutation: str, expected_position: tuple[Location, Location]
) -> None:
    """
    GIVEN a mutation
    WHEN we extract the position from hgvs
    THEN we should get a 3-integer tuple for start and end
    """
    hgvs = "ENST123.5:c." + mutation
    V1 = Variant(hgvs, consequences=["missense"])
    assert get_position(V1.hgvs) == expected_position


def test_variant_position_no_change_protein() -> None:
    """
    GIVEN a protein variant with no amino acid change
    WHEN we extract the position from the HGVS
    THEN we should not get an error
    """
    hgvs = "ENSP123.5:p.Leu615="
    assert get_position(hgvs) == (Location(0, 615, 0), Location(0, 615, 0))


# fmt: off
POSITION_CRITERION = [
    ("10","15",(
        Location(0,10,0),
        Location(0,15,0)
    )),
    ("15",None,(
        Location(0,15,0),
        Location(0,15,0)
    ))
]
# fmt: on


@pytest.mark.parametrize("start,end,expected_position_crit", POSITION_CRITERION)
def test_criterion_position_calculation(
    start: str, end: str, expected_position_crit: tuple[Location, Location]
) -> None:
    """
    GIVEN a criterion
    WHEN we convert start and end into positions
    THEN we should get a 3-integer tuple for start and end
    """
    C1 = Criterion("ENST123.5", start=start, end=end)
    assert (C1.start, C1.end) == expected_position_crit


# fmt: off
REGION = [
    # region1 lies before region2
    ((1,2),(5,9),False),
    # r1 ends on start r2
    ((3,5),(5,9),True),
    # r1 ends inside r2
    ((3,6),(5,9),True),
    # r1 ends on end r2
    ((3,9),(5,9),True),
    # r2 is inside r1
    ((3,10),(5,9),True),
    # r1 starts on start r2, ends inside r2
    ((5,6),(5,9),True),
    # r1 starts on start r2, ends on end r2
    ((5,9),(5,9),True),
    # r1 starts on start r2, ends outside r2
    ((5,10),(5,9),True),
    # r1 starts and ends inside r2
    ((6,7),(5,9),True),
    # r1 starts in r2, ends on end r2
    ((6,9),(5,9),True),
    # r1 starts in r2, ends outside r2
    ((6,10),(5,9),True),
    # r1 starts at end r2, ends outside r2
    ((9,10),(5,9),True),
    # r1 lies after r2
    ((10,11),(5,9),False)
]
# fmt: on


@pytest.mark.parametrize("region1, region2, overlap", REGION)
def test_region_overlap(region1: Region, region2: Region, overlap: bool) -> None:
    """
    GIVEN two regions
    WHEN we call the overlap
    THEN we should get whether the regions overlap

    Important: assumes inclusive positions, so (1,5) and (5, 6) overlap
    """

    assert region_overlap(region1, region2) == overlap


@pytest.mark.parametrize("region1, region2, overlap", REGION)
def test_region_overlap_reverse(
    region1: Region, region2: Region, overlap: bool
) -> None:
    """
    GIVEN two regions
    WHEN we call the overlap
    THEN we should get whether the regions overlap

    Important: assumes inclusive positions, so (1,5) and (5, 6) overlap
    """

    assert region_overlap(region2, region1) == overlap


def test_match_region(variant: Variant) -> None:
    """
    GIVEN a criterion and a variant in the specified region
    WHEN we check the match
    THEN we should see that the variant matches the criterion
    """
    c = Criterion("ENST123.5", consequence="inframe", start="10", end="20")
    assert c.match(variant) == True


def test_match_coordinate(variant: Variant) -> None:
    """
    GIVEN a Variant and Criterion
    WHEN we test if the coordinates match
    THEN they should only match when they are identical
    """
    c = Criterion("ENST123.5", coordinate="c")
    g = Criterion("ENST123.5", coordinate="g")

    assert c.match_coordinate(variant)
    assert not g.match_coordinate(variant)


def test_criterion_coordinate_mismatch(variant: Variant) -> None:
    """
    GIVEN a Variant and Criterion where only the coordinates don't match
    WHEN we match them
    THEN the Criterion should not match
    """
    c = Criterion("ENST123.5", coordinate="r", consequence="missense")
    assert not c.match(variant)


def test_Variant_from_VEP(minimal_vep: VEP) -> None:
    """
    GIVEN a VEP record
    WHEN we generate Variants from the VEP record
    THEN we should get a Variant for every hgvs description in every transcript_consequences
    """
    expected = [
        Variant("chr13:g.28034118_28034147dup", ["inframe_insertion"]),
        Variant("ENST00000241453.1:c.1772_1801dup", ["inframe_insertion"]),
        Variant(
            "ENSP00000241453.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            ["inframe_insertion"],
        ),
        Variant(
            "chr13:g.28034118_28034147dup",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
        Variant(
            "ENST00000380987.1:c.1772_1801dup",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
        Variant(
            "ENSP00000370374.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
    ]
    assert list(Variant.from_VEP(minimal_vep)) == expected
