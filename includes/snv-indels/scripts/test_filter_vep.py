#!/usr/bin/env python3

from filter_vep import (
    gene_of_interest,
    transcript_of_interest,
    consequence_of_interest,
    consequences_of_interest,
    vep_of_interest,
    update_most_severe,
)

import pytest


def make_consequence(
    gene, transcript, consequence_terms=["stop_lost"], impact="MODIFIER"
):
    return {
        "gene_id": gene,
        "transcript_id": transcript,
        "consequence_terms": consequence_terms,
        "impact": impact,
    }


@pytest.fixture
def gene1():
    return make_consequence("gene1", "transcript1", impact="HIGH")


@pytest.fixture
def gene2():
    return make_consequence("gene2", "transcript2", impact="HIGH")


@pytest.fixture
def vep():
    return {
        "most_severe_consequence": None,
        "transcript_consequences": [
            make_consequence("gene1", "transcript1", ["transcript_ablation"]),
            make_consequence(
                "gene2", "transcript2", ["splice_acceptor_variant"], "HIGH"
            ),
            make_consequence(
                "gene3", "transcript3", ["inframe_insertion", "stop_gained"]
            ),
        ],
    }

# Compare vs gene1
GENES_OF_INTEREST = [
        ({"gene1"}, True),
        ({"gene2", "gene3", "gene1"}, True),
        ({"gene9"}, False),
]


# Compare vs gene1
TRANSCRIPTS_OF_INTEREST = [
        ({"transcript1"}, True),
        ({"transcript2", "transcript1"}, True),
        ({"transcript9"}, False),
]

@pytest.mark.parametrize(["goi", "boolean"], GENES_OF_INTEREST)
def test_gene_of_interest(gene1, goi, boolean):
    """Test if gene1 is a gene of interest

    Relative to genes of interest in goi
    """
    assert gene_of_interest(gene1, goi) == boolean


@pytest.mark.parametrize(["toi", "boolean"], TRANSCRIPTS_OF_INTEREST)
def test_transcript_of_interest(gene1, toi, boolean):
    """Test if gene1 has transcript of interest

    Relative to transcripts of interest toi
    """
    assert transcript_of_interest(gene1, toi) == boolean


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


def test_consequence_of_interest_impact():
    genes = [f"gene{i}" for i in [1, 2, 3]]
    trans = [f"transcript{i}" for i in [1, 2, 3]]

    c1 = make_consequence("gene1", "transcript1")
    c2 = make_consequence("gene2", "transcript2", impact="HIGH")
    c3 = make_consequence("gene3", "transcript3", impact="MODERATE")

    assert consequence_of_interest(c1, genes, trans, "MODIFIER")
    assert not consequence_of_interest(c1, genes, trans, "HIGH")

    assert consequence_of_interest(c2, genes, trans, "HIGH")
    assert not consequence_of_interest(c2, genes, trans, "MODIFIER")

    # Test without an impact specified
    assert consequence_of_interest(c3, genes, trans, "MODERATE")
    assert consequence_of_interest(c3, genes, trans, None)


def test_multiple_consequences_of_interest(gene1, gene2):
    """Test returning multiple consequences of interest"""
    consequences = [gene1, gene2]
    genes = {"gene1", "gene2"}
    transcripts = {"transcript1", "transcript2"}
    filtered = consequences_of_interest(consequences, genes, transcripts)
    assert consequences == filtered


def test_multiple_consequences_of_interest_HIGH(vep):
    genes = "gene1 gene2 gene3".split()
    transcript = "transcript1 transcript2 transcript3".split()

    cons = consequences_of_interest(
        vep["transcript_consequences"], genes, transcript, "HIGH"
    )
    assert len(cons) == 1

    con = cons[0]
    assert con["gene_id"] == "gene2"


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
