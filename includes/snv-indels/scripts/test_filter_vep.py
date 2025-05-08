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


CRITERIA= [
        # Criteria, readout of gene names
        ([Criterion("ENTS0123.1")], ["gene1"]),
        ([Criterion("ENTS0123.1"), Criterion("ENTS0124.1")], ["gene1", "gene2"]),
        ([Criterion("ENST999.9")], []),
        ([Criterion("ENTS0123.1", consequence="transcript_ablation")], ["gene1"]),
        ([Criterion("ENTS0125.1", consequence="stop_gained")], ["gene3"]),
        ([Criterion("ENTS0125.1", consequence="inframe_insertion")], ["gene3"]),
]

@pytest.mark.parametrize(["criteria", "readout"], CRITERIA)
def test_filter_criteria(vep: VEP, criteria: list[Criterion],
                                 readout: list[str]) -> None:
    """Test filtering by consequence_term"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_criteria(criteria)
    genes = [tc["gene_id"] for tc in vep["transcript_consequences"]]
    assert genes == readout


def test_filter_criteria_no_match(vep: VEP) -> None:
    """Test that we get an empty list if no criteria matches"""
    vep.filter_criteria([Criterion("no_such_transcript")])
    assert not vep["transcript_consequences"]


def test_filter_no_criteria(vep: VEP) -> None:
    """Check that we filter nothing if there are no criteria"""
    vep.filter_criteria([])
    assert not vep["transcript_consequences"]


def test_vep_of_interest_one_transcript(vep: VEP) -> None:
    """Test rewriting the VEP object when there is a single transcript of
    interest"""
    c = Criterion("ENTS0125.1")

    # Restrict to transcript of interest
    vep.filter_criteria([c])

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
    c = Criterion("ENTS0123.1")
    empty.filter_criteria([c])
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
