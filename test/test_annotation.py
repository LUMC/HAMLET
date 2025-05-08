import gzip
import json
import pytest
import pathlib

HOTSPOTS = [
        (8390),
        (8860),
]

@pytest.mark.workflow('test-snv-indels-chrM')
@pytest.mark.parametrize("pos", HOTSPOTS)
def test_annotation(workflow_dir, pos):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/snv-indels/{sample}.vep.annotated.txt.gz")

    # Read the file, and pick out the line we want based on pos
    with gzip.open(output_file, "rt") as fin:
        for line in fin:
            d = json.loads(line)
            if d["start"] == pos:
                ts = d["transcript_consequences"][0]
                assert ts["annotation"] == "Hotspot"
                break
        else:
            raise RuntimeError(f"Position {pos} not found in {output_file}")
