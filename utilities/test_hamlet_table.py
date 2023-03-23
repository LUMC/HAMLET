#!/usr/bin/env python3

import pytest
import json

from hamlet_table import HAMLET_V1
from hamlet_table import HAMLET_V2


@pytest.fixture
def v1():
    v1_output = "test/data/output/v1/SRR8615409.summary.json"
    with open(v1_output) as fin:
        return HAMLET_V1(json.load(fin))


def test_v1_variants(v1):
    variants = list(v1.variants)
    assert len(variants) == 1
    var = variants[0]
    assert var["CHROM"] == "chr1"
    assert var["POS"] == 114716123
    assert var["FREQ"] == pytest.approx(0.46, 0.1)
    assert var["is_in_hotspot"]
    assert isinstance(var["Existing_variation"], list)


def test_v1_fusion(v1):
    fusions = list(v1.fusions)
    fusion = fusions[0]

    assert fusion["jr_count"] == 33
    assert fusion["name"] == "SFPQ--ZFP36L1"
    assert fusion["sf_count"] == 37


def test_v1_overexpression(v1):
    overexpression = list(v1.overexpression)
    exp = overexpression[0]

    assert exp["divisor_gene"] == "HMBS"
    assert exp["count"] == 0
    assert exp["ratio"] == 0.0
    assert exp["above_threshold"] == False
    assert exp["divisor_exp"] == 208030


def test_v1_itd_flt3(v1):
    events = list(v1.itd("flt3"))
    assert len(events) == 2


def test_split_by_consequence_none():
    """Test splitting VEP variant object without consequence"""
    var = {"some_value": 1}

    # Test that the variant is unchanged
    assert next(HAMLET_V2.split_by_consequence(var)) == var


def test_split_by_consequence_one():
    """Test splitting VEP variant object with a single consequence"""
    var = {"some_value": 1, "transcript_consequences": [1]}


    # Test that the variant is unchanged
    assert next(HAMLET_V2.split_by_consequence(var)) == var


def test_split_by_consequence_two():
    """TEST splitting VEP variant object with multiple consequences"""
    var = {"some_value": 1, "transcript_consequences": [1, 2]}
    for i, split in enumerate(HAMLET_V2.split_by_consequence(var), 1):
        assert split["transcript_consequences"] == [i]


def test_extract_alt_heterozygous():
    ref = "A"
    genotype = "A/T"
    assert HAMLET_V1.get_alt(ref, genotype) == "T"


def test_extract_alt_hom_alt():
    ref = "A"
    genotype = "T/T"
    assert HAMLET_V1.get_alt(ref, genotype) == "T"


def test_two_alts():
    """Not supported"""
    ref = "A"
    genotype = "T/C"
    with pytest.raises(NotImplementedError):
        HAMLET_V1.get_alt(ref, genotype)


def test_vcf_pos_SNV():
    """For SNV's, the vcf POS is the same as the VEP start and end"""
    var = {"variant_class": "SNV", "start": 10}
    assert HAMLET_V2.vcf_pos(var) == 10


def test_vcf_pos_insertion():
    """For insertions, the vcf POS start - 1"""
    var = {"variant_class": "insertion", "start": 11}
    assert HAMLET_V2.vcf_pos(var) == 10


def test_vcf_pos_deletion():
    """For deletions, the vcf POS start - 1"""
    var = {"variant_class": "deletion", "start": 11}
    assert HAMLET_V2.vcf_pos(var) == 10


def test_vcf_pos_sequence_variation():
    """For sequence variation, the vcf POS start - 1"""
    var = {"variant_class": "sequence variation", "start": 11}
    assert HAMLET_V2.vcf_pos(var) == 10
