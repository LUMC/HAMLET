import json
import pytest
import pathlib

HOTSPOTS = [
        (8390, True),
        (8566, False),
        (8701, True),
        (8860, True),
]

@pytest.mark.workflow('test-snv-indels-chrM')
@pytest.mark.parametrize(["pos", "hotspot"], HOTSPOTS)
def test_is_in_hotspot2(workflow_dir, pos, hotspot):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/snv-indels/{sample}.vep.high.txt")

    # Read the file, and pick out the line we want based on pos
    with open(output_file) as fin:
        for line in fin:
            d = json.loads(line)
            if d["start"] == pos:
                assert d["is_in_hotspot"] == hotspot
                break
        else:
            raise RuntimeError(f"Position {pos} not found in {output_file}")
