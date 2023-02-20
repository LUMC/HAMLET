#!/usr/bin/env python3

from vep_goi import (
    gene_of_interest,
    transcript_of_interest,
    consequence_of_interest,
    consequences_of_interest,
    vep_of_interest,
    update_most_severe,
)

import pytest


def make_consequence(gene, transcript, consequence_terms=["stop_lost"]):
    return {
            "gene_id": gene,
            "transcript_id": transcript,
            "consequence_terms": consequence_terms
    }


@pytest.fixture
def gene1():
    return make_consequence("gene1", "transcript1")


@pytest.fixture
def gene2():
    return make_consequence("gene2", "transcript2")


@pytest.fixture
def vep():
    return {
        "most_severe_consequence": None,
        "transcript_consequences": [
            make_consequence("gene1", "transcript1", ["transcript_ablation"]),
            make_consequence("gene2", "transcript2", ["splice_acceptor_variant"]),
            make_consequence("gene3", "transcript3", ["inframe_insertion", "stop_gained"]),
        ]
    }


def test_is_gene_of_interest(gene1):
    """Test recognition of gene of interest"""
    assert gene_of_interest(gene1, {"gene1"})
    assert gene_of_interest(gene1, {"gene2", "gene3", "gene1"})


def test_is_not_gene_of_interest(gene1):
    """Test recognition of gene not of interest"""
    assert not gene_of_interest(gene1, {"gene9"})


def test_is_transcript_of_interest(gene1):
    """Test recognition of transcript of interest"""
    assert transcript_of_interest(gene1, {"transcript1"})
    assert transcript_of_interest(gene1, {"transcript2", "transcript1"})


def test_is_not_transcript_of_interest(gene1):
    """Test recognition of transcript not of interest"""
    assert not transcript_of_interest(gene1, {"transcript12"})


def test_consequence_of_interest(gene1):
    """Test recognition of consequence of interest"""
    assert consequence_of_interest(gene1, ["gene1"], ["transcript1"])


def test_consequences_of_interest(gene1, gene2):
    """Test filtering only consequences of interest"""
    consequences = [gene1, gene2]
    genes = {"gene1"}
    transcripts = {"transcript1"}
    filtered = consequences_of_interest(consequences, genes, transcripts)
    assert len(filtered) == 1
    con = filtered[0]
    assert con["gene_id"] == "gene1"
    assert con["transcript_id"] == "transcript1"


def test_multiple_consequences_of_interest(gene1, gene2):
    """Test returning multiple consequences of interest"""
    consequences = [gene1, gene2]
    genes = {"gene1", "gene2"}
    transcripts = {"transcript1", "transcript2"}
    filtered = consequences_of_interest(consequences, genes, transcripts)
    assert consequences == filtered


def test_vep_of_interest_one_transcript(vep):
    """Test rewriting the VEP object when there is a single transcript of
    interest"""
    genes = ["gene3"]
    transcripts = ["transcript3"]

    # Restrict
    new_vep = vep_of_interest(vep, genes, transcripts)

    # Get the consequences after rewriting the VEP obejct
    cons = new_vep["transcript_consequences"]

    # Test that only gene3 is left
    assert len(cons) == 1
    cons[0]["gene_id"] == "gene3"
    cons[0]["transcript_id"] == "transcript"
    cons[0]["consequence_terms"] == ["stop_gained"]

    # Test that we set the most severe consequence for transcript3
    assert new_vep["most_severe_consequence"] == "stop_gained"


def test_update_most_severe(vep):
    update_most_severe(vep)
    assert vep["most_severe_consequence"] == "transcript_ablation"

    # Remove gene1
    vep["transcript_consequences"].pop(0)
    update_most_severe(vep)
    assert vep["most_severe_consequence"] == "splice_acceptor_variant"

    # Remove gene2
    vep["transcript_consequences"].pop(0)
    update_most_severe(vep)
    assert vep["most_severe_consequence"] == "stop_gained"
