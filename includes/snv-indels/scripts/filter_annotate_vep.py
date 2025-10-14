#!/usr/bin/env python3

# Authors: Anne van der Grinten, Redmar van den Berg

import argparse
import json

from collections.abc import Iterator


from utils import VEP, read_criteria_file, read_known_variants


def parse_vep_json(vep_file: str) -> Iterator[VEP]:
    """Parse the VEP 'json' output file, each line contains a JSON entry"""
    with open(vep_file) as fin:
        for line in fin:
            yield VEP(json.loads(line))


def main(
    vep_file: str,
    criteria_file: str,
    annotation_file: str,
    known_variants_file: str,
    population: str,
    frequency: float,
) -> None:
    # Get genes and transcripts of interest
    annotations = read_criteria_file(annotation_file)
    filter_criteria = read_criteria_file(criteria_file)

    # Add the filter criteria to the end of the annotations, so they have a
    # lower priority
    annotations.update(filter_criteria)

    # TODO, read known variants file
    known_variants = read_known_variants(known_variants_file) if known_variants_file else dict()

    for vep in parse_vep_json(vep_file):
        # Skip variants that are above the specified population frequency
        if vep.above_population_threshold(population, frequency):
            continue

        # Filter and annotate transcripts based on known variants and the annotations
        # The known variants have the higher priority
        vep.filter_annotate_transcripts(known_variants, annotations)

        # If there is no consequence of interest left
        if not vep["transcript_consequences"]:
            continue

        print(json.dumps(vep, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract genes (and transcript) of interest from VEP output"
    )

    parser.add_argument("--vep", help="VEP json output file")
    parser.add_argument("--filter-criteria", help="File with filter criteria")
    parser.add_argument("--annotation-criteria", help="File with annotation criteria")
    parser.add_argument("--known-variants", help="File with known variants")
    parser.add_argument(
        "--population", help="Population to use for variant frequency", default="gnomAD"
    )
    parser.add_argument(
        "--frequency",
        help="Variants with a population frequency above this threshold will be filtered out",
        default=0.05,
        type=float,
    )

    args = parser.parse_args()

    main(
        args.vep,
        args.filter_criteria,
        args.annotation_criteria,
        args.known_variants,
        args.population,
        args.frequency,
    )
