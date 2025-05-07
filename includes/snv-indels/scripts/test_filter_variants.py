import pytest

# Authors: Anne van der Grinten, Redmar van den Berg

from filter_variants import (
    Criterion,
    Variant,
    VEP,
    Location,
    Region,
    region_overlap,
    get_position,
)


@pytest.fixture
def variant() -> Variant:
    return Variant("ENST123.5:c.10A>T", consequences=["missense", "inframe"])


@pytest.fixture
def vep() -> VEP:
    minimal_vep_record = {
        "transcript_consequences": [
            {
                "consequence_terms": ["inframe_insertion"],
                "hgvsg": "chr13:g.28034118_28034147dup",
                "hgvsc": "ENST00000241453.1:c.1772_1801dup",
                "hgvsp": "ENSP00000241453.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            },
            {
                "consequence_terms": ["inframe_insertion", "NMD_transcript_variant"],
                "hgvsg": "chr13:g.28034118_28034147dup",
                "hgvsc": "ENST00000380987.1:c.1772_1801dup",
                "hgvsp": "ENSP00000370374.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            },
        ]
    }

    return minimal_vep_record


def test_create_criterion() -> None:
    C = Criterion("ENST123.5")
    assert C.identifier == "ENST123.5"


def test_create_criterion_Locations() -> None:
    """
    WHEN we initialize a criterion
    THEN start and end must be Location tuples
    """
    C = Criterion("ENST123.5", start="10", end="*20-5")

    assert C.start == Location(0, 10, 0)
    assert C.end == Location(1, 20, -5)


def test_create_variant(variant: Variant) -> None:
    assert variant.hgvs == "ENST123.5:c.10A>T"
    assert variant.consequences == ["missense", "inframe"]


def test_variant_match(variant: Variant) -> None:
    """not
    GIVEN a variant and a criterion with different transcripts
    WHEN we check for a match
    THEN we should not find a match
    """
    C = Criterion("ENST120.5")

    assert not C.match(variant)


def test_transcript_match(variant: Variant) -> None:
    """
    GIVEN two identical transcripts
    WHEN we check for a match
    THEN we should find a match
    """
    C = Criterion("ENST123.5")
    assert C.match_id(variant)


def test_transcript_no_match(variant: Variant) -> None:
    """
    GIVEN two different transcripts
    WHEN we check for a match
    THEN we should not find a match
    """
    C = Criterion("ENST120.5")
    assert not C.match_id(variant)


def test_transcript_version_no_match(variant: Variant) -> None:
    """
    GIVEN a version mismatch between two transcripts
    WHEN we check for a match
    THEN we should get a runtime error
    """
    C = Criterion("ENST123.4")
    with pytest.raises(RuntimeError):
        C.match_id(variant)


def test_transcript_id_no_version(variant: Variant) -> None:
    """
    GIVEN a Transcript and Criterion without version number
    WHEN we compare them
    THEN they should match
    """
    V = Variant("Chr12:g.10A>T", consequences=[])
    C = Criterion("Chr12", coordinate="g")

    assert C.match_id(V)


def test_transcript_version_no_match_error(variant: Variant) -> None:
    """
    GIVEN a version mismatch between criterion and variant
    WHEN we check for a match
    THEN we should get a runtime error
    """
    C = Criterion("ENST123.4")
    with pytest.raises(RuntimeError):
        C.match(variant)


CONS = [(None, True), ("missense", True), ("stop lost", False)]


@pytest.mark.parametrize("consequence, expected", CONS)
def test_match_consequence(
    consequence: str | None, expected: bool, variant: Variant
) -> None:
    """
    GIVEN a variant and a criterion with identical transcript and a consequence
    WHEN we check for a match
    THEN we should get the answer that is expected
    """
    C = Criterion("ENST123.5", "c", consequence)
    assert C.match(variant) == expected


# fmt: off
POSITION = [
    ("-20del",(
        Location(0,-20,0),
        Location(0,-20,0)
    )),
    ("-10+5del",(
        Location(0,-10,5),
        Location(0,-10,5)
    )),
    ("-9-5del",(
        Location(0,-9,-5),
        Location(0,-9,-5)
    )),
    ("1del",(
        Location(0,1,0),
        Location(0,1,0)
    )),
    ("10+5del",(
        Location(0,10,+5),
        Location(0,10,+5)
    )),
    ("11-5del",(
        Location(0,11,-5),
        Location(0,11,-5)
    )),
    ("*1del",(
        Location(1,1,0),
        Location(1,1,0)
    )),
    ("*5+5del",(
        Location(1,5,5),
        Location(1,5,5)
    )),
    ("*6-5del",(
        Location(1,6,-5),
        Location(1,6,-5)
    )),
    ("10_15del",(
        Location(0,10,0),
        Location(0,15,0)
    )),
    #in case start and end are specified antisense in criterion
    ("15_10del",(
        Location(0,10,0),
        Location(0,15,0)
    )),
]
# fmt: on


@pytest.mark.parametrize("mutation, expected_position", POSITION)
def test_variant_position_calculation(
    mutation: str, expected_position: tuple[Location, Location]
) -> None:
    """
    GIVEN a mutation
    WHEN we extract the position from hgvs
    THEN we should get a 3-integer tuple for start and end
    """
    hgvs = "ENST123.5:c." + mutation
    V1 = Variant(hgvs, consequences=["missense"])
    assert get_position(V1.hgvs) == expected_position


def test_variant_position_no_change_protein() -> None:
    """
    GIVEN a protein variant with no amino acid change
    WHEN we extract the position from the HGVS
    THEN we should not get an error
    """
    hgvs = "ENSP123.5:p.Leu615="
    assert get_position(hgvs) == (Location(0, 615, 0), Location(0, 615, 0))


# fmt: off
POSITION_CRITERION = [
    ("10","15",(
        Location(0,10,0),
        Location(0,15,0)
    )),
    ("15","10",(
        Location(0,10,0),
        Location(0,15,0)
    )),
    ("15",None,(
        Location(0,15,0),
        Location(0,15,0)
    ))
]
# fmt: on


@pytest.mark.parametrize("start,end,expected_position_crit", POSITION_CRITERION)
def test_criterion_position_calculation(
    start: str, end: str, expected_position_crit: tuple[Location, Location]
) -> None:
    """
    GIVEN a criterion
    WHEN we convert start and end into positions
    THEN we should get a 3-integer tuple for start and end
    """
    C1 = Criterion("ENST123.5", start=start, end=end)
    assert (C1.start, C1.end) == expected_position_crit


# fmt: off
REGION = [
    # region1 lies before region2
    ((1,2),(5,9),False),
    # r1 ends on start r2
    ((3,5),(5,9),True),
    # r1 ends inside r2
    ((3,6),(5,9),True),
    # r1 ends on end r2
    ((3,9),(5,9),True),
    # r2 is inside r1
    ((3,10),(5,9),True),
    # r1 starts on start r2, ends inside r2
    ((5,6),(5,9),True),
    # r1 starts on start r2, ends on end r2
    ((5,9),(5,9),True),
    # r1 starts on start r2, ends outside r2
    ((5,10),(5,9),True),
    # r1 starts and ends inside r2
    ((6,7),(5,9),True),
    # r1 starts in r2, ends on end r2
    ((6,9),(5,9),True),
    # r1 starts in r2, ends outside r2
    ((6,10),(5,9),True),
    # r1 starts at end r2, ends outside r2
    ((9,10),(5,9),True),
    # r1 lies after r2
    ((10,11),(5,9),False)
]
# fmt: on


@pytest.mark.parametrize("region1, region2, overlap", REGION)
def test_region_overlap(region1: Region, region2: Region, overlap: bool) -> None:
    """
    GIVEN two regions
    WHEN we call the overlap
    THEN we should get whether the regions overlap

    Important: assumes inclusive positions, so (1,5) and (5, 6) overlap
    """

    assert region_overlap(region1, region2) == overlap


@pytest.mark.parametrize("region1, region2, overlap", REGION)
def test_region_overlap_reverse(
    region1: Region, region2: Region, overlap: bool
) -> None:
    """
    GIVEN two regions
    WHEN we call the overlap
    THEN we should get whether the regions overlap

    Important: assumes inclusive positions, so (1,5) and (5, 6) overlap
    """

    assert region_overlap(region2, region1) == overlap


def test_match_region(variant: Variant) -> None:
    """
    GIVEN a criterion and a variant in the specified region
    WHEN we check the match
    THEN we should see that the variant matches the criterion
    """
    c = Criterion("ENST123.5", consequence="inframe", start="10", end="20")
    assert c.match(variant) == True


def test_match_coordinate(variant: Variant) -> None:
    """
    GIVEN a Variant and Criterion
    WHEN we test if the coordinates match
    THEN they should only match when they are identical
    """
    c = Criterion("ENST123.5", coordinate="c")
    g = Criterion("ENST123.5", coordinate="g")

    assert c.match_coordinate(variant)
    assert not g.match_coordinate(variant)


def test_criterion_coordinate_mismatch(variant: Variant) -> None:
    """
    GIVEN a Variant and Criterion where only the coordinates don't match
    WHEN we match them
    THEN the Criterion should not match
    """
    c = Criterion("ENST123.5", coordinate="r", consequence="missense")
    assert not c.match(variant)


def test_Variant_from_VEP(vep: VEP) -> None:
    """
    GIVEN a VEP record
    WHEN we generate Variants from the VEP record
    THEN we should get a Variant for every hgvs description in every transcript_consequences
    """
    expected = [
        Variant("chr13:g.28034118_28034147dup", ["inframe_insertion"]),
        Variant("ENST00000241453.1:c.1772_1801dup", ["inframe_insertion"]),
        Variant(
            "ENSP00000241453.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            ["inframe_insertion"],
        ),
        Variant(
            "chr13:g.28034118_28034147dup",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
        Variant(
            "ENST00000380987.1:c.1772_1801dup",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
        Variant(
            "ENSP00000370374.1:p.Asp600_Leu601insHisValAspPheArgGluTyrGluTyrAsp",
            ["inframe_insertion", "NMD_transcript_variant"],
        ),
    ]
    assert list(Variant.from_VEP(vep)) == expected
