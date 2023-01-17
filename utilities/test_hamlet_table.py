#!/usr/bin/env python3

import pytest
import json

from hamlet_table import HAMLET_V1

@pytest.fixture
def v1():
    v1_output ='test/data/output/v1/SRR8615409.summary.json'
    with open(v1_output) as fin:
        return HAMLET_V1(json.load(fin))


def test_v1_variants(v1):
    variants = list(v1.variants)
    assert len(variants) == 1
    var = variants[0]
    assert var['CHROM'] == 'chr1'
    assert var['POS'] == 114716123
    assert var['FREQ'] == pytest.approx(0.46, 0.1)
    assert var['is_in_hotspot']
    assert isinstance(var['Existing_variation'], list)


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
