#!/usr/bin/env python3

from dataclasses import dataclass


@dataclass
class Mapping:
    gene_id: str
    gene_name: str
    transcript_ids: set[str]


def read_attributes(fin):
    for line in fin:
        if line.startswith("#"):
            continue
        spline = line.strip().split(maxsplit=8)
        kvals = spline[-1]

        d = dict()
        for pair in kvals.split(";"):
            if not pair:
                continue
            pair = pair.strip(" ")
            k, v = pair.split(" ", maxsplit=1)
            d[k] = v.replace('"', "")
        yield d


def create_mapping(gtf_file: str, transcripts: set[str]) -> dict[str, Mapping]:
    """Create the mapping for each transcript in transcripts"""
    results = dict()
    with open(gtf_file) as fin:
        for record in read_attributes(fin):
            transcript_id = record.get("transcript_id")
            if transcript_id not in transcripts or transcript_id is None:
                continue
            gene_id = record["gene_id"]
            gene_name = record["gene_name"]

            if gene_id in results:
                results[gene_id].transcript_ids.add(transcript_id)
            else:
                results[gene_id] = Mapping(gene_id, gene_name, {transcript_id})
    return results


def read_transcripts(filter_file: str) -> set[str]:
    """Read the transcripts of interest from the filter criteria file"""
    transcripts = set()
    header = None
    with open(filter_file) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")

            if header is None:
                header = spline
                continue

            d = {k: v for k, v in zip(header, spline)}
            transcripts.add(d["transcript_id"])
    return transcripts


def main(gtf_file: str, filter_file: str) -> None:

    transcripts = read_transcripts(filter_file)
    results = create_mapping(gtf_file, transcripts)

    # Make sure we didn't miss any transcripts
    found_transcripts = set()
    for mapping in results.values():
        found_transcripts.update(mapping.transcript_ids)

    if missing := transcripts - found_transcripts:
        raise RuntimeError(f"Unable to find transcripts: {', '.join(missing)}")

    # Print header
    print("GOI_ID", "GOI_SYMBOL", "TOI_IDS", sep="\t")
    for r in results.values():
        print(r.gene_id, r.gene_name, ",".join(r.transcript_ids), sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("gtf", help="gtf file")
    parser.add_argument(
        "--filter-file", help="File with filter criteria", required=True
    )
    # parser.add_argument(
    #     "--transcripts", nargs="+", help="Transcripts of interest", required=True
    # )

    args = parser.parse_args()
    main(args.gtf, args.filter_file)
