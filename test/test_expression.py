import json
import pathlib
import pytest

@pytest.mark.workflow('test-expression-chrM')
def test_output_targetted_seq_against_schema(workflow_dir):
    """ Test content of exon ratio json file """
    sample = "SRR8615409"
    exon_ratio = pathlib.Path(workflow_dir, f"{sample}/expression/{sample}.exon_ratios.json")
    with open(exon_ratio) as fin:
        js = json.load(fin)["expression"][0]
    assert js["above_threshold"]
    assert js["count"] == 219
