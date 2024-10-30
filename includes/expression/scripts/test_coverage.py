import pysam
import pytest
import dataclasses

from typing import Any

from coverage import forward_orientation


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


def to_pysam_read(fake_read: FakeRead, header: pysam.AlignmentHeader) -> pysam.AlignedSegment:
    dataclass_dict = dataclasses.asdict(fake_read)
    return pysam.AlignedSegment.from_dict(dataclass_dict, header)


@pytest.fixture
def reads():

    reads = [
            # paired, proper pair, mate reverse, first in pair
            FakeRead("read1", "99"),
            # paired, proper pair, mate reverse, second in pair
            FakeRead("read2", "163"),
            # paired, proper pair, read reverse, first in pair
            FakeRead("read3", "83"),
            # paired, proper pair, read reverse, second in pair
            FakeRead("read4", "147"),
    ]

    header = {
            "SQ": [{"SN": "chr1", "LN": 100}]
    }

    H = pysam.AlignmentHeader.from_dict(header)

    return [to_pysam_read(read, H) for read in reads]


def test_true(reads):
    assert True

    # Orientation of the first read of the pair, for the reads in the fixture
    expected = ["+", "-", "-", "+"]

    for e, r in zip(expected, reads):
        assert forward_orientation(r) == e

    print()
    for read in reads:
        print(read)
