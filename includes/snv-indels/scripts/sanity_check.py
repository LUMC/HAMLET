#!/usr/bin/env python3

import argparse

from utils import VEP, read_criteria_file, Variant
from annotate_vep import read_known_variants

def main(filter_criteria_file, annotation_criteria_file, known_variants_file):
    filter_criteria = read_criteria_file(filter_criteria_file)
    annotation_criteria = read_criteria_file(annotation_criteria_file)

    known_variants = [Variant(hgvs, consequences=[]) for hgvs in read_known_variants(known_variants_file).keys()]

    annotation_errors = list()
    # Test that each annotation criteria falls in at least one filter criteria
    for annotation in annotation_criteria:
        for filter in filter_criteria:
            if filter.contains(annotation):
                print("CONTAINED")
                print("annotation:", annotation)
                print("filter    :", filter)
                print()
                break
        else:
            annotation_errors.append(annotation)

    if annotation_errors:
        print("The following annotation criteria can never be met:")
        print("\n".join([str(criteria) for criteria in annotation_errors]))

    variant_errors = list()
    for variant in known_variants:
        for filter in filter_criteria:
            if filter.match(variant):
                print("MATCH")
                print(f"Variant {variant} matches {filter}")
                break
        else:
            variant_errors.append(variant)

    if variant_errors:
        print("The following variants can never be found:")
        print("\n".join([str(variant) for variant in variant_errors]))

    if annotation_errors or variant_errors:
        exit(-1)
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "Check the internal consistency of the filter and annotation criteria"
    )

    parser.add_argument("--filter-criteria", help="File with filter criteria", required=True)
    parser.add_argument("--annotation-criteria", help="File with annotation criteria", required=True)
    parser.add_argument("--known-variants", help="File with known variants")
    args = parser.parse_args()

    main(
        args.filter_criteria,
        args.annotation_criteria,
        args.known_variants
    )
