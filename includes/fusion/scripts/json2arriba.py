#!/usr/bin/env python3

from arriba2json import arriba_header, json_to_arriba
import argparse
import json


def main(fusion_file: str) -> None:
    with open(fusion_file) as fin:
        fusions = json.load(fin)

    # Prin the header
    print("#", "\t".join(arriba_header), sep="")

    for fusion in fusions:
        print(json_to_arriba(arriba_header, fusion))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fusions", help="Arriba fusions in json format")

    args = parser.parse_args()

    main(args.fusions)
