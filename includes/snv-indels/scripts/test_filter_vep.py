#!/usr/bin/env python3

from filter_vep import VEP

import json
import pytest


@pytest.fixture
def vep():
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

@pytest.mark.parametrize(["transcripts", "length"], FILTER_TRANSCRIPTS)
def test_filter_transcript_id(vep, transcripts, length):
    """Test filtering by transcript_id"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_transcript_id(transcripts)
    assert len(vep["transcript_consequences"]) == length


@pytest.mark.parametrize(["consequences", "gene"], FILTER_CONSEQUENCE)
def test_filter_consequence_term(vep, consequences, gene):
    """Test filtering by consequence_term"""
    assert len(vep["transcript_consequences"]) == 3
    vep.filter_consequence_term(consequences)
    vep["transcript_consequences"][0] == gene


def test_filter_consequence_emtpy(vep):
    """Test that we get an empty list if no consequence matches"""
    vep.filter_consequence_term({"no_such_consequence"})
    assert not vep["transcript_consequences"]


def test_vep_of_interest_one_transcript(vep):
    """Test rewriting the VEP object when there is a single transcript of
    interest"""
    transcripts = ["transcript3"]

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


def test_update_most_severe(vep):
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
