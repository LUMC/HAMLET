import pysam
import pytest
import dataclasses


@dataclasses.dataclass
class FakeRead():
    name: str
    flag: str 
    ref_name: str = "chr1"
    ref_pos: str = "10"
    next_ref_name: str = "="
    next_ref_pos: str = "20"
    seq: str = "ATCG"
    qual: str = "@@@@"
    map_quality: str = "255"
    cigar: str = "4M"
    length: str = "4"
    tags = []


@pytest.fixture
def reads():

    reads = [
            FakeRead("read1", "99")
    ]

    header = {
            "SQ": [ {"SN": "chr1", "LN": 100}]
    }

    H = pysam.AlignmentHeader.from_dict(header)

    return [pysam.AlignedSegment.from_dict(dataclasses.asdict(x), H) for x in reads]


def test_true(reads):
    assert True

    print()
    for read in reads:
        print(read)
