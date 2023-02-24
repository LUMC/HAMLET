#!/usr/bin/env python3

import importlib
import pytest
from collections import defaultdict

json_output = importlib.import_module("json-output")
update_variant_overview = json_output.update_variant_overview
get_gene_id = json_output.get_gene_id
get_gene_symbol = json_output.get_gene_symbol
idf_to_gene_symbol = json_output.idf_to_gene_symbol
get_format = json_output.get_format
get_info = json_output.get_info


@pytest.fixture
def id_mapping():
    """Dummy gene to transcript ID mapping"""
    gene1 = { "gene_id": "gene1",
              "gene_symbol": "GENE1",
              "transcript_ids": ["transcript1"]
            }
    gene2 = { "gene_id": "gene2",
              "gene_symbol": "GENE2",
              "transcript_ids": ["transcript1", "transcript2"]
            }
    return [gene1, gene2]


@pytest.fixture
def vcf_line():
    """The VCF input is included in the VEP json object"""
    return "chrM\t13650\t.\tC\tT\t.\tPASS\tADP=17;WT=0;HET=0;HOM=1;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t1/1:81:17:17:0:15:100%:6.4467E-9:0:40:0:0:10:5"

@pytest.fixture
def mapping():
    return {"gene1": "GENE1", "gene2": "GENE2"}


@pytest.fixture
def vep_single():
    """A single vep entry, with a single transcript"""
    transcript = { "gene_id": "gene1",
            "transcript_id": "transcript1"}
    vep = { "transcript_consequences": [transcript] }

    return vep


def test_get_gene_id():
    assert get_gene_id({"gene_id": "gene1"}) == "gene1"


def test_idf_to_gene_symbol(id_mapping):
    gene_symbol = idf_to_gene_symbol(id_mapping)
    assert gene_symbol == {"gene1": "GENE1", "gene2": "GENE2"}


def test_get_gene_symbol(vep_single, mapping):
    assert get_gene_symbol(vep_single, mapping) == "GENE1"


def test_update_variant_overview(mapping, vep_single):
    overview = defaultdict(list)
    update_variant_overview(mapping, vep_single, overview)
    assert "GENE1" in overview


def test_get_format(vcf_line):
    FORMAT = get_format(vcf_line)
    assert FORMAT["GT"] == "1/1"
    assert FORMAT["ADR"] == "5"


def test_get_info(vcf_line):
    INFO = get_info(vcf_line)
    assert INFO["ADP"] == "17"
    assert INFO["HOM"] == "1"
