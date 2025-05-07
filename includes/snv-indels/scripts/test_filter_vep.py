#!/usr/bin/env python3

from filter_vep import VEP, FrequenciesType
from filter_variants import Criterion

import json
import pytest

from typing import Any, Dict, Set, List, Tuple


@pytest.fixture
def vep() -> VEP:
    with open("test/data/output/v2/vep.json") as fin:
        js = json.load(fin)
    return VEP(js[0])


# transcripts, transcript_consequences lenght
FILTER_TRANSCRIPTS = [
        ({"ENTS0123.1"}, 1),
        ({"ENTS0123.1", "ENTS0124.1"}, 2),
        ({"ENST999.9"}, 0),
]

# consequence_terms, first gene_id
FILTER_CONSEQUENCE = [
        ({"transcript_ablation"}, "gene1"),
        ({"stop_gained"}, "gene3"),
        ({"transcript_ablation", "inframe_insertion"}, "gene1"),
]

@pytest.mark.parametrize(["transcripts", "length"], FILTER_TRANSCRIPTS)
def test_filter_transcript_id(vep: VEP, transcripts: Set[str],
                              length: int) -> None:
    """Test filtering by transcript_id"""
    criteria = [Criterion(id) for id in transcripts]
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_criterion(criteria)
    assert len(vep["transcript_consequences"]) == length


@pytest.mark.parametrize(["consequences", "gene"], FILTER_CONSEQUENCE)
def test_filter_consequence_term(vep: VEP, consequences: Set[str],
                                 gene: str) -> None:
    """Test filtering by consequence_term"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_consequence_term(consequences)
    assert vep["transcript_consequences"][0]["gene_id"] == gene


def test_filter_consequence_emtpy(vep: VEP) -> None:
    """Test that we get an empty list if no consequence matches"""
    vep.filter_consequence_term({"no_such_consequence"})
    assert not vep["transcript_consequences"]


def test_filter_consequence_none(vep: VEP) -> None:
    """Check that we filter nothing if the consequences are empty"""
    before = vep["transcript_consequences"]
    vep.filter_consequence_term(set())
    after = vep["transcript_consequences"]
    assert before == after


def test_vep_of_interest_one_transcript(vep: VEP) -> None:
    """Test rewriting the VEP object when there is a single transcript of
    interest"""
    transcripts = {"transcript3"}

    # Restrict to transcript of interest
    vep.filter_transcript_id(transcripts)

    # Get the consequences after rewriting the VEP obejct
    cons = vep["transcript_consequences"]

    # Test that only gene3 is left
    assert len(cons) == 1
    assert cons[0]["gene_id"] == "gene3"
    assert cons[0]["transcript_id"] == "transcript3"
    assert cons[0]["consequence_terms"] == ["inframe_insertion","stop_gained"]

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
    transcripts = {"transcript1"}
    empty.filter_transcript_id(transcripts)
    empty.update_most_severe()


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
    V = VEP(data)
    with pytest.raises(RuntimeError):
        V._extract_frequencies()

def test_error_extract_frequencies_multiple_entries() -> None:
    """
    GIVEN a VEP record with a single "frequencies" entry, which contains multiple keys
    WHEN we extract the frequency values
    THEN we get an error
    """
    data: Dict[str, Any] = {"colocated_variants": [
            {"frequencies":
                {
                    "T": dict(),
                    "A": dict()
                }
            }
        ]
    }
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
def test_above_population_threshold(population: str, threshold:float, expected: bool) -> None:
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
    V =VEP(data)

    assert V.above_population_threshold(population, threshold) is expected
