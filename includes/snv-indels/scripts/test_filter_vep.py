#!/usr/bin/env python3

from filter_vep import VEP

import json
import pytest

from typing import Dict, Set


@pytest.fixture
def vep() -> VEP:
    with open("test/data/output/v2/vep.json") as fin:
        js = json.load(fin)
    return VEP(js[0])


# transcripts, transcript_consequences lenght
FILTER_TRANSCRIPTS = [
        ({"transcript1"}, 1),
        ({"transcript1", "transcript2"}, 2),
        ({"transcript9"}, 0),
]

# consequence_terms, first gene_id
FILTER_CONSEQUENCE = [
        ({"transcript_ablation"}, "gene1"),
        ({"stop_gained"}, "gene3"),
        ({"transcript_ablation", "inframe_insertion"}, "gene1"),
]

# hgvsc variant, transcript_consequence length
FILTER_BLACKLIST = [
    (set(), 3),
    ({"ENTS0123.1:c.40A>C"}, 2),
    ({"ENTS0123.1:c.40A>C", "ENTS0124.1:c.40A>C", "ENTS0125.1:c.40A>C"}, 1),
]


@pytest.mark.parametrize(["transcripts", "length"], FILTER_TRANSCRIPTS)
def test_filter_transcript_id(vep: VEP, transcripts: Set[str],
                              length: int) -> None:
    """Test filtering by transcript_id"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_transcript_id(transcripts)
    assert len(vep["transcript_consequences"]) == length


@pytest.mark.parametrize(["consequences", "gene"], FILTER_CONSEQUENCE)
def test_filter_consequence_term(vep: VEP, consequences: Set[str],
                                 gene: str) -> None:
    """Test filtering by consequence_term"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_consequence_term(consequences)
    vep["transcript_consequences"][0] == gene


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
    cons[0]["gene_id"] == "gene3"
    cons[0]["transcript_id"] == "transcript"
    cons[0]["consequence_terms"] == ["stop_gained"]

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


@pytest.mark.parametrize(["blacklist", "length"], FILTER_BLACKLIST)
def test_filter_transcript_blacklist(vep: VEP, blacklist: Set[str],
                                     length: int) -> None:
    """Test filtering on hgvsc field"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_hgvsc_blacklist(blacklist)
    assert len(vep["transcript_consequences"]) == length


def test_most_severe_transcript_blacklist(vep: VEP) -> None:
    """Test if most severe consequence is updated"""
    vep.filter_hgvsc_blacklist({"ENTS0123.1:c.40A>C"})
    assert vep["most_severe_consequence"] == "splice_acceptor_variant"


COLOCATED_VARIANTS = [
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
]

@pytest.mark.parametrize(["vep", "frequencies"], COLOCATED_VARIANTS)
def test_extract_frequencies(vep: VEP, frequencies: Dict[str, float]) -> None:
    V = VEP(vep)
    assert V.extract_frequencies() == frequencies


def test_error_extract_multiple_frequencies() -> None:
    """
    GIVEN a VEP record with multiple "frequencies" entries
    WHEN we extract the frequency values
    THEN we get an error
    """
    data = {
            "colocated_variants": [
                {
                    "frequencies": dict()
                },
                {
                    "frequencies": dict()
                }
            ]
    }
    V = VEP(data)
    with pytest.raises(RuntimeError):
        V.extract_frequencies()

def test_error_extract_frequencies_multiple_entries() -> None:
    """
    GIVEN a VEP record with a single "frequencies" entry, which contains multiple keys
    WHEN we extract the frequency values
    THEN we get an error
    """
    data = {"colocated_variants": [
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
        V.extract_frequencies()
