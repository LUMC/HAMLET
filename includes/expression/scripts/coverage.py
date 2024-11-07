#!/usr/bin/env python3

import argparse
from typing import Any
import pysam
from statistics import median

from collections import defaultdict
from dataclasses import dataclass

from typing import Iterator, Optional
from gtf import gene_id_name


@dataclass
class Bed:
    """Class for BED records"""

    chrom: str
    chromStart: int
    chromEnd: int
    name: str
    score: int
    strand: str


@dataclass
class Coverage:
    """Class to store coverage information by strand"""

    unstranded: float
    forward: float
    reverse: float


def parse_bed(fname: str):
    """Parse BED file and yield records"""
    with open(fname) as fin:
        for line in fin:
            spline: list[Any] = line.strip("\n").split("\t")

            # To int
            for field in [1, 2, 4]:
                spline[field] = int(spline[field])

            B = Bed(*spline)

            yield B


def get_reads(record: Bed, bamfile: str) -> Iterator[pysam.AlignedSegment]:
    """Get reads that overlap Bed record

    To mimic the results from the STAR counts table, ignore reads that are
     - supplementary
     - secondary
     - not in proper pairs
    """
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for read in samfile.fetch(
            contig=record.chrom, start=record.chromStart, end=record.chromEnd
        ):
            if read.is_supplementary or read.is_secondary or not read.is_proper_pair:
                continue
            yield read


def read_names(Bed: Bed, bamfile):
    """Extract the read names in Bed region, per strandedness"""

    # Store the read name to ensure we don't count forward/reverse double
    by_strand: dict[str, set[str]] = {
        "unstranded": set(),
        "forward": set(),
        "reverse": set(),
    }

    for read in get_reads(Bed, bamfile):
        orientation = orientation_first(read)
        read_name = read.query_name
        if orientation == Bed.strand:
            by_strand["forward"].add(read_name)
        else:
            by_strand["reverse"].add(read_name)

        # For unstranded, we count both orientations
        by_strand["unstranded"].add(read_name)

    return by_strand


def get_bed_coverage(bedfile: str, bamfile: str) -> dict[str, Coverage]:
    """Get coverage for the regions specified in the bedfile"""
    # If there are multiple regions in the BED file with the same name, we add
    # the read names together so we don't count double
    by_name: dict[str, dict[str, set]] = defaultdict(
        lambda: {"unstranded": set(), "forward": set(), "reverse": set()}
    )

    # Final coverage counts
    coverage = dict()

    for record in parse_bed(bedfile):
        names = read_names(record, bamfile)
        for strand in ["unstranded", "forward", "reverse"]:
            by_name[record.name][strand].update(names[strand])

    for name, reads in by_name.items():
        u = reads["unstranded"]
        f = reads["forward"]
        r = reads["reverse"]

        coverage[name] = Coverage(len(u), len(f), len(r))

    return coverage


def orientation_first(read: pysam.AlignedSegment) -> str:
    """Determine the orientation of the first read of the pair"""
    if read.is_read1 and read.is_forward:
        return "+"
    elif read.is_read2 and read.is_forward:
        return "-"
    elif read.is_read1 and read.is_reverse:
        return "-"
    elif read.is_read2 and read.is_reverse:
        return "+"

    msg = f"Unexpected read encountered: {read.query_name}"
    raise RuntimeError(msg)


def read_STAR_counts(fname):
    """Read the STAR counts file into Coverage by strand"""
    with open(fname) as fin:
        # Skip the headers
        for _ in range(4):
            next(fin)

        for line in fin:
            spline = line.strip("\n").split("\t")
            ensg = spline[0]
            expression = [int(spline[column]) for column in [1, 2, 3]]
            yield ensg, Coverage(*expression)


def counts_by_name(countsfile: str, gtffile: str) -> dict[str, Coverage]:
    """Read the STAR counts by gene name"""
    counts = dict()
    # Read the gft files to get the ENSG to name mapping
    with open(gtffile) as fin:
        ensg_to_name = gene_id_name(fin)

    # Read the STAR counts file
    for ensg, coverage in read_STAR_counts(countsfile):
        # If the ENSG has no name in the GTF, keep it
        name = ensg_to_name.get(ensg, ensg)
        counts[name] = coverage

    return counts


def get_normalizer_values(
    STAR_counts: dict[str, Coverage], genes: list[str]
) -> Coverage:
    """Calculate the normalizer values for each strand"""

    unstranded = [STAR_counts[gene].unstranded for gene in genes]
    forward = [STAR_counts[gene].forward for gene in genes]
    reverse = [STAR_counts[gene].reverse for gene in genes]

    return Coverage(median(unstranded), median(forward), median(reverse))


def main(
    bamfile: str,
    countsfile: str,
    gtffile: str,
    housekeeping: list[str],
    bedfile: Optional[str],
    genes: list[str],
    raw: bool,
) -> None:
    if not bedfile and not genes:
        print("Nothing to do")
        exit(0)

    coverage: dict[str, Coverage] = dict()
    # Read the counts file, by gene name instead of ENSG
    STAR_counts = counts_by_name(countsfile, gtffile)

    # Calculate the normalizer for the housekeeping genes
    normalizer = get_normalizer_values(STAR_counts, housekeeping)

    if bedfile:
        coverage = get_bed_coverage(bedfile, bamfile)

    # Get the genes of interest from the STAR_counts
    for gene in genes:
        coverage[gene] = STAR_counts[gene]

    # Normalize the coverage
    if not raw:
        for gene in coverage:
            coverage[gene].unstranded /= normalizer.unstranded
            coverage[gene].forward /= normalizer.forward
            coverage[gene].reverse /= normalizer.reverse

    for gene, cov in coverage.items():
        print(gene, cov.unstranded, cov.forward, cov.reverse, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("--bam", required=True, help="BAM file")
    parser.add_argument("--counts", required=True, help="STAR counts file")
    parser.add_argument(
        "--housekeeping",
        required=True,
        help="List of housekeeping gene names",
        nargs="+",
    )
    parser.add_argument(
        "--gtf", required=True, help="GTF file, used to find the ENSG for each gene"
    )

    # Optional arguments
    parser.add_argument(
        "--bed",
        nargs="?",
        default="",
        const="",
        help="BED file with regions of interest",
    )
    parser.add_argument(
        "--genes", nargs="*", default=list(), help="List of genes of interest"
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        default=False,
        help="Do not normalize the expression counts",
    )

    args = parser.parse_args()

    main(
        args.bam,
        args.counts,
        args.gtf,
        args.housekeeping,
        args.bed,
        args.genes,
        args.raw,
    )
