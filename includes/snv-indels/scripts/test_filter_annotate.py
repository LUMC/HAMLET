from typing import Sequence
import pytest

from utils import (
    Criterion,
    VEP,
)

"""
Tests for the filtering/annotation of transcripts for variants, based on Criteria and known variants


The test setup is as follows:
- We have a single VEP record, which contains 3 transcript consequences
- We have 3 filter/annotation Criteria, which match each of the 3 transcript
  consequences
- We have 3 known variants, which match each of the 3 transcript consequences

The tests combine the Criteria with the known_variants to ensure that the
annotations and filtering works correctly, with the expected priority.
"""
@pytest.fixture
def vep() -> VEP:
    minimal_vep_record = {
        "transcript_consequences": [
            {
                "hgvsc": "ENST1.1:c.100del",
                "consequence_terms": ["frameshift"],
            },
            {
                "hgvsc": "ENST2.1:c.100A>T",
                "consequence_terms": ["missense", "splice site"],
            },
            {
                "hgvsc": "ENST3.1:c.100+100A>T",
                "consequence_terms": [],
            },
        ]
    }
    return VEP(minimal_vep_record)


@pytest.fixture
def criteria() -> Sequence[Criterion]:
    return [
        Criterion("ENST1.1", consequence = "frameshift"),
        Criterion("ENST2.1", start="0", end="200", consequence="splice site"),
        Criterion("ENST3.1", start="100+20", end="101-20")
    ]

@pytest.fixture
def known_variants() -> Sequence[str]:
    return [
        "ENST1.1:c.100del",
        "ENST2.1:c.100A>T",
        "ENST3.1:c.100+100A>T"
    ]


def test_filter_annotate_variants_empty(vep: VEP) -> None:
    vep.filter_annotate_transcripts(dict(), dict())
    assert not vep["transcript_consequences"]

@pytest.mark.parametrize(
    "v, c, expected_annotations",
    [
        # No annotations means all transcript consequences are filtered out
        ("", "", []),
        # All three Criteria
        ("", "012", ["crit", "crit", "crit"]),
        # All three known variants
        ("012", "", ["var", "var", "var"]),
        # Both criteria and known variants, variants take precedence
        ("012", "012", ["var", "var", "var"]),
        # The second variant is not in known variants
        ("02", "012", ["var", "crit", "var"]),
        ("02", "1", ["var", "crit", "var"]),
        ("02", "", ["var", "var"]),
    ]
)
def test_filter_annotate_variants(vep: VEP, criteria: Sequence[Criterion], known_variants: Sequence[str], c: str, v: str, expected_annotations: list[str]) -> None:
    # Create the dict for the criteria
    criteria_dict = dict()
    for index in c:
        criteria_dict[criteria[int(index)]] = "crit"

    # Create the dict for the known variants
    known_variants_dict = dict()
    for index in v:
        known_variants_dict[known_variants[int(index)]] = "var"


    vep.filter_annotate_transcripts(known_variants_dict, criteria_dict)

    annotations = [transcript["annotation"] for transcript in vep["transcript_consequences"]]

    assert annotations == expected_annotations
