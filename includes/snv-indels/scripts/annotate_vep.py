#!/usr/bin/env python3

import argparse
import json
from utils import Variant, VEP, read_criteria_file
from itertools import zip_longest
from collections import OrderedDict

from typing import Iterator


def parse_vep_json(vep_file: str) -> Iterator[VEP]:
    """Parse the VEP 'json' output file, each line contains a JSON entry"""
    with open(vep_file) as fin:
        for line in fin:
            yield VEP(json.loads(line))

def read_known_variants(fname: str) -> dict[str, str]:
    raise NotImplementedError

def main(vep_file: str, annotations_file: str, known_variants_file: str | None):
    criteria = read_criteria_file(annotations_file)
    known_variants = read_known_variants(known_variants_file) if known_variants_file else dict()

    for vep in parse_vep_json(vep_file):
        for transcript in vep["transcript_consequences"]:
            # Make sure hgvsc is defined
            hgvsc = transcript.get("hgvsc")
            # If there is no HGVSC, we do nothing
            if hgvsc is None:
                continue
            # If hgvsc is a known variant, we update the annotation
            if hgvsc in known_variants:
                transcript["annotation"] = known_variants[hgvsc]
                continue
            # Otherwise, we check the criteria
            variant = Variant(hgvsc, transcript["consequence_terms"])
            for crit, annotation in criteria.items():
                if crit.match(variant):
                    transcript["annotation"] = annotation
                    break
        print(json.dumps(vep, sort_keys=True))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract genes (and transcript) of interest from VEP output"
    )

    parser.add_argument("--vep", help="VEP json output file", required=True)
    parser.add_argument("--annotations", help="Criteria file with annotations", required=True)
    parser.add_argument("--known-variants", help="File with known pathogenic variants", nargs='?')

    args = parser.parse_args()

    main(
        args.vep,
        args.annotations,
        args.known_variants
    )
