import gzip
import json
import pytest
import pathlib

HOTSPOTS = [
    (8390, "Hotspot"),
    (8860, "known variant"),
    (8566, "artifact"),
    (8701, "artifact"),
]


@pytest.mark.workflow("test-snv-indels-chrM")
@pytest.mark.parametrize("pos, annotation", HOTSPOTS)
def test_annotation(workflow_dir: str, pos: int, annotation: str) -> None:
    sample = "SRR8615409"
    output_file = pathlib.Path(
        workflow_dir, f"{sample}/snv-indels/{sample}.vep.annotated.txt.gz"
    )

    # Read the file, and pick out the line we want based on pos
    with gzip.open(output_file, "rt") as fin:
        for line in fin:
            d = json.loads(line)
            if d["start"] == pos:
                ts = d["transcript_consequences"][0]
                assert ts["annotation"] == annotation
                break
        else:
            raise RuntimeError(f"Position {pos} not found in {output_file}")
