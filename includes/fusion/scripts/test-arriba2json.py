#!/usr/bin/env python3

from arriba2json import arriba_to_json, json_to_arriba

header = [
    "gene1",
    "gene2",
    "read_identifiers",
    "split_reads_1",
    "closest_genomic_breakpoint1",
]

# A fusion event in json format
json_format = {
    "gene1": "BCR",
    "gene2": "ABL",
    "read_identifiers": ["read1", "read2"],
    "split_reads_1": 15,
    "closest_genomic_breakpoint1": None,
}

# A fusion event in Arriba tsv format
arriba_format = "BCR\tABL\tread1,read2\t15\t."


def test_arriba_to_json() -> None:
    assert arriba_to_json(header, arriba_format) == json_format


def test_json_to_arriba() -> None:
    assert json_to_arriba(header, json_format) == arriba_format
