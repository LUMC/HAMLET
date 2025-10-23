#!/usr/bin/env python3

from dataclasses import dataclass
import sys


@dataclass
class Mapping:
    gene_id: str
    gene_name: str
    transcript_ids: set[str]

    def __str__(self):
        return f"{self.gene_id}\t{self.gene_name}\t{','.join(self.transcript_ids)}"


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
            if (transcript_id := record.get("transcript_id")) is None:
                continue

            # Transscript can have a version number, or not
            if transcript_id not in transcripts:
                continue
            gene_id = record["gene_id"]
            gene_name = record["gene_name"]

            if gene_id in results:
                results[gene_id].transcript_ids.add(transcript_id)
            else:
                results[gene_id] = Mapping(gene_id, gene_name, {transcript_id})
    return results


def read_transcripts(filter_file: str) -> set[str]:
    """Read the transcripts of interest from the filter criteria file

    Ignore the transcript version
    """
    transcripts = set()
    header = None
    with open(filter_file) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")

            if header is None:
                header = spline
                continue

            d = {k: v for k, v in zip(header, spline)}
            # Strip the version from the transcript id
            transcript = d["transcript_id"].split(".")[0]
            transcripts.add(transcript)
    return transcripts

def known_variants(fname: str) -> set[str]:
    """Read the transcripts from the known variants file"""
    transcripts = set()
    header = None
    with open(fname) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")
            if header is None:
                header = spline
                continue
            d = {k: v for k, v in zip(header, spline)}
            # Get the transcript identifier from the HGVS
            transcript_id = d["variant"].split(":c.")[0]
            # Remove the version number
            transcript = transcript_id.split(".")[0]
            transcripts.add(transcript)
    return transcripts

def main(gtf_file: str, filter_file: str, annotation_file: str | None, known_variants_file: str |None) -> None:

    transcripts = read_transcripts(filter_file)

    if annotation_file:
        transcripts.update(read_transcripts(annotation_file))

    if known_variants_file:
        transcripts.update(known_variants(known_variants_file))

    results = create_mapping(gtf_file, transcripts)

    # Make sure we didn't miss any transcripts
    found_transcripts = set()
    for mapping in results.values():
        found_transcripts.update(mapping.transcript_ids)

    if missing := transcripts - found_transcripts:
        print("Missing transcripts:", file=sys.stderr)
        print(*missing, sep="\n", file=sys.stderr)
        raise RuntimeError(f"Unable to find {len(missing)} transcripts")

    # Print header
    print("GOI_ID", "GOI_SYMBOL", "TOI_IDS", sep="\t")
    for mapping in results.values():
        print(mapping)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("gtf", help="gtf file")
    parser.add_argument(
        "--filter-file", help="File with filter criteria", required=True
    )
    parser.add_argument("--annotation-file", help="File with annotation criteria")
    parser.add_argument("--known-variants", help="File with known variants")

    args = parser.parse_args()
    main(args.gtf, args.filter_file, args.annotation_file, args.known_variants)
