import json
import pytest
import pathlib


@pytest.mark.workflow('test-snv-indels-chrM')
def test_is_in_hotspot2(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/snv-indels/{sample}.vep.high.txt")

    # Expected results, there should be 4 variants in the vep.high.txt file
    positions = [8390, 8566, 8701, 8860]
    is_in_hotspot = [True, False, True, True]

    with open(output_file) as fin:
        for pos, hotspot, line in zip(positions, is_in_hotspot, fin):
            d = json.loads(line)
            assert d["start"] == pos
            assert d["is_in_hotspot"] == hotspot
