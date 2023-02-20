#!/usr/bin/env python3

from vep_goi import (
    gene_of_interest,
    transcript_of_interest,
    consequence_of_interest,
    consequences_of_interest,
)

import pytest


def make_consequence(gene, transcript):
    return {"gene_id": gene, "transcript_id": transcript}


@pytest.fixture
def gene1():
    return make_consequence("gene1", "transcript1")


@pytest.fixture
def gene2():
    return make_consequence("gene2", "transcript2")


@pytest.fixture
def vep():
    return {
        "transcript_consequences": [
            make_consequence("gene1", "transcript1"),
            make_consequence("gene2", "transcript2"),
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
