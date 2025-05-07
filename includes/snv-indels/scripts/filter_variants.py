#!/usr/bin/python3

# Authors: Anne van der Grinten, Redmar van den Berg

from mutalyzer_hgvs_parser import to_model  # type: ignore
from collections import namedtuple
from collections.abc import Iterator
import functools

from typing import Any, Tuple

# Tuple to store locations of HGVS points
Location = namedtuple("Location", ["downstream", "position", "offset"])

# Tuple to store a region
Region = namedtuple("Region", ["start", "end"])

VEP = dict[str, Any]


class Variant:
    def __init__(self, hgvs: str, consequences: list[str]):
        self.hgvs = hgvs
        self.consequences = consequences

    @classmethod
    def from_VEP(cls, vep: VEP) -> Iterator["Variant"]:
        """Yield all Variants from a VEP record"""
        for transcript in vep.get("transcript_consequences", list()):
            consequences = transcript["consequence_terms"]
            for hgvs in ["hgvsg", "hgvsc", "hgvsp"]:
                if hgvs in transcript:
                    yield Variant(transcript[hgvs], consequences)

    def __str__(self) -> str:
        return f"Variant({self.hgvs}, {self.consequences})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError
        return str(self) == str(other)


class Criterion:
    def __init__(
        self,
        identifier: str,
        coordinate: str = "c",
        consequence: str | None = None,
        start: str | None = None,
        end: str | None = None,
    ) -> None:
        self.identifier = identifier
        self.coordinate = coordinate
        self.consequence = consequence

        # Define the type of start and end
        self.start: Location | None
        self.end: Location | None

        if start is not None:
            # if end is not specified, set end equal to start
            if end is None:
                end = start
            # Make a fake hgvs description to determine the start and end Location
            hgvs = f"{identifier}:c.{start}_{end}"
            self.start, self.end = get_position(hgvs)
        else:
            self.start = None
            self.end = None

    def __str__(self) -> str:
        return (
            f"Criterion(identifier={self.identifier}, "
            f"coordinate={self.coordinate}, "
            f"consequence={self.consequence}, "
            f"start={self.start}, "
            f"end={self.end})"
        )

    def __repr__(self) -> str:
        return str(self)

    def match(self, variant: Variant) -> bool:
        return all(
            [
                self.match_id(variant),
                self.match_coordinate(variant),
                self.match_consequence(variant),
                self.match_region(variant),
            ]
        )

    def match_id(self, variant: Variant) -> bool:
        """Determine if the identifier matches the variant

        Raises a runtime error if only the version doesn't match, since this
        could indicate the use of incompatible reference seqeunces
        """

        def split_version(variant_id: str) -> Tuple[str, str]:
            """Split the version from the identifier

            Set version to 0 if there is no version
            """
            try:
                id_, version = variant_id.split(".")
            # If there is no version (e.g. chr5), set version to 0
            except ValueError:
                id_ = variant_id
                version = "0"

            return id_, version

        # Get the variant id
        variant_id = variant.hgvs.split(":")[0]
        criterion_id = self.identifier

        id1, v1 = split_version(variant_id)
        id2, v2 = split_version(criterion_id)

        if id1 != id2:
            return False

        if v1 == v2:
            return True
        else:
            msg = f"Version mismatch between {variant} and {self}"
            raise RuntimeError(msg)

    def match_coordinate(self, variant: Variant) -> bool:
        variant_coordinate = variant.hgvs.split(":")[1].split(".")[0]
        return variant_coordinate == self.coordinate

    def match_consequence(self, variant: Variant) -> bool:
        if self.consequence is None:
            return True
        return self.consequence in variant.consequences

    def match_region(self, variant: Variant) -> bool:
        if self.start is None:
            return True
        var_region = get_position(variant.hgvs)
        crit_region = Region(self.start, self.end)

        return region_overlap(var_region, crit_region)


@functools.lru_cache(maxsize=1000)
def get_position(hgvs: str) -> Region:
    # Get the variant part of the hgvs description
    coordinate, var = hgvs.split(":")[1].split(".")

    if coordinate == "p":
        model = to_model(var, "p_variant")
    else:
        model = to_model(var, "variant")

    location = model["location"]
    if location["type"] == "point":
        start = point_to_tuple(location)
        end = start
    elif location["type"] == "range":
        start = point_to_tuple(location["start"])
        end = point_to_tuple(location["end"])

        # make sure "end" is downstream of "start"
        if start > end:
            s, e = start, end
            start = e
            end = s
    else:
        raise RuntimeError

    return Region(start, end)


def point_to_tuple(point: dict[str, Any]) -> Location:
    """Convert a 'point' dictionary from parse_hgvs into a Location"""

    # Default values
    downstream = 0
    position = 0
    offset = 0

    position = point["position"]
    if "outside_cds" in point:
        if point["outside_cds"] == "upstream":
            position *= -1
        elif point["outside_cds"] == "downstream":
            downstream = 1
        else:
            raise RuntimeError
    if "offset" in point:
        offset = point["offset"]["value"]

    return Location(downstream, position, offset)


def region_overlap(region1: Region, region2: Region) -> bool:
    return any(
        [
            region1[0] >= region2[0] and region1[0] <= region2[1],
            region1[1] >= region2[0] and region1[1] <= region2[1],
            region1[0] < region2[0] and region1[1] > region2[1],
        ]
    )
