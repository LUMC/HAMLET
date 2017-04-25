#!/usr/bin/env python

import concurrent.futures
import json
import os
import sys
from collections import Counter, namedtuple
from enum import Enum
from math import log
from os import path
from pathlib import Path

import click
import pysam
from align.calign import aligner
from Bio import SeqIO

__author__ = ["Wibowo Arindrarto", "Daniel Borras"]
__contact__ = "w.arindrarto@lumc.nl"


INSERT_SC = 1
SPLICE_SC = 3
CIGAR_SC = 4
REF_SC = (0, 2, 3, 6, 7, 8)


class Region(namedtuple("Region", ["contig", "start", "end", "name"])):

    """Helper class for representing a region."""

    __slots__ = ()

    @property
    def length(self):
        return self.end - self.start


class SCType(Enum):

    """Enumeration of possible soft clip location relative to a read."""

    start = 0
    end = 1


def get_bam_sample(bam):
    """Returns the BAM sample name based on the RG:SM tag, or its filename."""
    return bam.header.get("RG", [{}])[0].get(
        "SM",
        path.basename(bam.filename.decode()))


def calc_cigar_bit(cigar_op):
    """Given a cigar operation integer, return the cigar bit."""
    # taken from htslib bam_cigar_type function
    return 0x3c1a7 >> (cigar_op << 1) & 3


def consumes_query(cigar_op):
    """Given a cigar operation integer, returns whether it consumes the
    query sequence."""
    return calc_cigar_bit(cigar_op) & 1


def consumes_ref(cigar_op):
    """Given a cigar operation integer, returns whether it consumes the
    reference sequence."""
    return calc_cigar_bit(cigar_op) & 2


def overlaps(reg, target_reg):
    """Given a region and another target region, returns whether the region
    overlaps the target region."""
    return (target_reg.start <= reg.start < target_reg.end) \
        or (target_reg.start <= reg.end < target_reg.end)


def envelops(reg, target_reg):
    """Given a region and another target region, returns whether the region
    is enveloped within the target region."""
    return (target_reg.start <= reg.start < target_reg.end) \
        and (target_reg.start <= reg.end < target_reg.end)


def advance_pos(init_ref_pos, init_query_pos, cigar_ops):
    """Advance `ref` and `query` positions based on the cigar string."""
    for op, op_len in cigar_ops:
        if consumes_query(op):
            init_query_pos += op_len
        if consumes_ref(op):
            init_ref_pos += op_len
    return init_ref_pos, init_query_pos


def extract_coord(reg_str):
    """Given a SAM-compatible genome coordinate, extract the values."""
    reg_str = reg_str.replace(",", "")
    try:
        contig, reg_str = reg_str.rsplit(":", 1)
    except ValueError:
        # No start and end specified
        return reg_str, None, None
    try:
        start, end = reg_str.rsplit("-", 1)
    except ValueError:
        # Only start specified
        return contig, int(reg_str) - 1, None
    # Start and end specified.
    # Convert coord to zero-based, half open.
    start, end = int(start) - 1, int(end)
    if start < 0:
        raise click.BadParameter("Start position must be at least 1.")
    if start > end:
        raise click.BadParameter(
            "Invalid interval: {0} - {1}.".format(start + 1, end))
    return contig, start, end


def within_fragment(read, sc_type):
    """Returns whether a soft-clipped region occurs within a fragment
    or not."""
    if sc_type not in (SCType.start, SCType.end):
        raise ValueError("Invalid soft clip sc_typeation " + sc_type + ".")
    # no mate or mate unmapped always false
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    # read has mate and is the leftmost
    if read.tlen > 0:
        return sc_type == SCType.end
    if read.tlen < 0:
        return sc_type == SCType.start
    msg = "Unexpected fragment configuration: for read at {0}".format(read.pos)
    print(msg, file=sys.stderr)


def get_inserts(read, min_insertion_length):
    """Returns a list of insert sequences and where they occur (ref-wise)
    inside the read"""
    ref_pos = read.pos
    read_pos = 0
    cigar_ops = read.cigartuples
    seq = read.seq
    inserts = []
    for op, op_len in cigar_ops:
        if consumes_query(op):
            if op == INSERT_SC and op_len >= min_insertion_length:
                insert_seq = seq[read_pos:read_pos + op_len]
                inserts.append((ref_pos, insert_seq))
            read_pos += op_len
        if consumes_ref(op):
            ref_pos += op_len
    return inserts


def start_is_soft_clipped(read, min_length=0):
    """Returns whether the beginning of a read is soft-clipped or not."""
    cigar = read.cigartuples
    if not cigar:
        msg = "Unexpected read at {0}: no cigar string".format(read.pos)
        print(msg, file=sys.stderr)
        return
    return cigar[0][0] == CIGAR_SC and cigar[0][1] >= min_length


def end_is_soft_clipped(read, min_length=0):
    """Returns whether the end of a read is soft-clipped or not."""
    cigar = read.cigartuples
    if not cigar:
        msg = "Unexpected read at {0}: no cigar string".format(read.pos)
        print(msg, file=sys.stderr)
        return
    return cigar[-1][0] == CIGAR_SC and cigar[-1][1] >= min_length


def count_pileups(bam, contig, start, end):
    """Returns the number of piled-up bases between the given
    start and end positions in the given contig."""
    raw_counts = Counter({x.pos: x.n
                          for x in bam.pileup(contig, start, end)
                          if start <= x.pos <= end})
    return [{"pos": pos, "count": raw_counts[pos]}
            for pos in range(start, end)]


def get_alt_sc_coords(sc_seq, ref, start_ref, end_ref, sc_type):
    """Returns a set of coordinates where the given soft-clipped sequence
    may align."""
    sc_len = len(sc_seq)
    # Only consider soft clips whose length is at least 1 + log4 of
    # the potential candidate region length to avoid getting hits
    # by chance.
    if sc_len < round(log(end_ref - start_ref, 4) + 1):
        return set([])
    alns = [aln
            for aln in aligner(sc_seq, ref[start_ref:end_ref],
                               matrix="DNAFULL", method="glocal",
                               gap_open=-7, gap_extend=-1,
                               gap_double=-7, max_hits=None)
            # Filter for alignments with at most 2 hits and at most
            # 10% mismatches, rounded up.
            if aln.n_gaps1 <= 2 and aln.n_mismatches <= round(0.1 * sc_len)]
    # Also ensure we are using regular Python types from here on.
    res = {int(aln.end2) + start_ref for aln in alns} \
        if sc_type == SCType.start else \
        {int(aln.start2) + start_ref - 1 for aln in alns}
    if len(res) > 1:
        return set([])
    return res


def process_read(read, target_reg, ref, min_sc_length, min_insertion_length):
    """Counts the insertions and soft clips present in the given read."""
    cigar = read.cigartuples
    query_pos = read.pos
    inserts = get_inserts(read, min_insertion_length)
    scs = []

    if start_is_soft_clipped(read):
        sc_reg = Region(target_reg.contig, query_pos - cigar[0][1],
                        query_pos, read.qname)
        # cigar[[0][1] denotes length of sc region
        if envelops(sc_reg, target_reg) and sc_reg.length >= min_sc_length:
            alt_sc_coords = get_alt_sc_coords(
                read.seq[:sc_reg.length],
                ref, query_pos, target_reg.end, SCType.start)
            for asc in (alt_sc_coords or {None}):
                scs.append((query_pos - 1, asc))

    if end_is_soft_clipped(read):
        adv_ref_pos, adv_query_pos = advance_pos(query_pos, 0, cigar[:-1])
        sc_reg = Region(target_reg.contig, adv_ref_pos,
                        adv_ref_pos + cigar[-1][1], read.qname)
        if envelops(sc_reg, target_reg) and sc_reg.length >= min_sc_length:
            alt_sc_coords = get_alt_sc_coords(
                read.seq[adv_query_pos:],
                ref, target_reg.start, adv_ref_pos, SCType.end)
            for asc in (alt_sc_coords or {None}):
                scs.append((adv_ref_pos, asc))

    return inserts, scs


def integrate_per_read_result(raw_inserts, raw_scs, counts_i, counts_sc):
    # counts_i: [(pos, seq)] of inserts
    for ipos, iseq in counts_i:
        raw_inserts[(ipos, iseq.upper())] += 1
    # counts_sc: [(pos, altPos)] of scs
    for scpos, altscpos in counts_sc:
        if scpos not in raw_scs:
            raw_scs[scpos] = {"count": 0, "altPosCount": {}}
        raw_scs[scpos]["count"] += 1
        if altscpos is not None:
            if altscpos not in raw_scs[scpos]["altPosCount"]:
                raw_scs[scpos]["altPosCount"][altscpos] = 0
            raw_scs[scpos]["altPosCount"][altscpos] += 1


def process_region(aln, ref, contig, start, end, nt, min_sc_length=3,
                   min_insertion_length=3, sample_id=None, output_zeros=True):
    raw_inserts, raw_scs = Counter(), {}
    result = {
        "bamFile": str(Path(aln.filename.decode()).resolve()),
        "sampleName": sample_id,
        "region": {
            "start": start,
            "end": end,
            "contig": contig,
        },
        "pileups": count_pileups(aln, contig, start, end),
    }
    # Add 1 bp padding to capture reads whose soft clip occurs just after the
    # boundary ~ we are interested in these but pysam doesn't consider
    # them overlaps (for the right reasons).
    reads = aln.fetch(contig, max(start - 1, 0), end + 1)

    target_reg = Region(contig, start, end, None)
    with concurrent.futures.ThreadPoolExecutor(max_workers=nt) as executor:
        futures = \
            {executor.submit(process_read, read, target_reg, ref,
                             min_sc_length, min_insertion_length): read.qname
             for read in reads}

        for future in concurrent.futures.as_completed(futures):
            try:
                counts_i, counts_sc = future.result()
            except Exception as exc:
                read_name = futures[future]
                print("Error when processing read {0!r}: {1}."
                      "".format(read_name, exc), file=sys.stderr)
            integrate_per_read_result(raw_inserts, raw_scs, counts_i,
                                      counts_sc)

    result["inserts"] = [{"pos": pos, "count": count, "seq": seq}
                         for (pos, seq), count in raw_inserts.items()]
    result["scs"] = [{"pos": pos, "count": item["count"],
                      "altPosCount": [{"pos": k, "count": v}
                                      for k, v in item["altPosCount"].items()]}
                     for pos, item in raw_scs.items()]

    return result


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("fasta",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("bam",
                type=click.Path(exists=True, dir_okay=False))
@click.option("--region", type=str,
              help="Region of soft clip and insertion counting.")
@click.option("--min-sc-length", type=int, default=3,
              help="Minimum length of soft-clipped region to count.")
@click.option("--min-insertion-length", type=int, default=3,
              help="Minimum length of insertion to count.")
@click.option("--sample-id", type=str,
              help="Name of the sample to which the reads belong."
                   " If not given, the sample name will be parsed"
                   " from the BAM header. If this is not possible,"
                   " the BAM file basename will be used.")
@click.option("--nthreads", type=int, default=os.cpu_count(),
              help="Number of threads to use.")
def main(fasta, bam, region, min_sc_length, min_insertion_length, sample_id,
         nthreads):
    """
    Counts soft clip and insertion events in a region of an indexed BAM file.

    \b
    Input
    =====

    \b
        * Region of a contig in the BAM file from which the soft clip
          and insertion events are counted. The coordinate format is
          '<contig>:<start>-<end>', where the first base is numbered
          '1' and the end coordinate is included in the region.
        * FASTA file of the BAM reference sequence.
        * Paired-end BAM alignment, position sorted and indexed.


    \b
    Output
    ======

    \b
        * JSON with the following pseudoschema (all coordinates are
          zero-based, half open):

    \b
        {
            "bamFile": <path to BAM file>,
            "sampleName": <name of the sample>,
            "region": {
                "contig": <contig name>,
                "start": <start coordinate>,
                "end": <end coordinate>
            },
            "scs": [{
                "pos": <base position>,
                "count": <number of reads with soft clipping at the position>,
                "altPosCount": [{
                    "pos": <reference position where soft clip may be aligned>,
                    "count": <number of soft clip alignable to the position>
                }],
            }],
            "inserts": [{
                "pos": <reference position where insertion occurs>,
                "sequence": <sequence of the insertion>,
                "count": <how many times the insertion occurs>
            }],
            "pileups": [
                "pos": <base position>,
                "count": <number of piled-up bases at the position>
            ]
        }

    Entries in the `altPosCount` object each denotes the position in the
    reference sequence where the soft-clipped sequence may be aligned. The
    position always refers to the reference position where the soft clip base
    is closest to the read. In other words, if the soft clip occurs at the 5'
    end of a read, this position refers to the 3' end of the sequence, and
    vice versa.


    Copyright (c) 2016 Leiden University Medical Center

    All rights reserved.

    """
    aln = pysam.AlignmentFile(bam)

    if region is not None:
        contig, start, end = extract_coord(region)
    else:
        try:
            contig, = aln.references
            start, end = None, None
        except ValueError:
            raise click.BadParameter("Contig is not specified and there"
                                     " is not exactly one contig in the"
                                     " alignment file.")

    if not aln.has_index():
        raise click.BadParameter("Alignment file is not indexed.")
    if contig not in aln.references:
        raise click.BadParameter("Contig {0!r} is not in the alignment file."
                                 "".format(contig))

    fa_recs = [r for r in SeqIO.parse(fasta, "fasta")
               if r.id == contig]
    if not fa_recs:
        raise click.BadParameter("Reference FASTA does not contain the contig"
                                 " {0!r}.".format(contig))
    elif len(fa_recs) > 1:
        raise click.BadParameter("Reference FASTA contains multiple contigs"
                                 " with the name {0!r}.".format(contig))
    ref = str(fa_recs.pop().seq)

    counts = process_region(
        aln, ref, contig,
        start or 0,
        end or aln.lengths[aln.references.index(contig)],
        nthreads, min_sc_length, min_insertion_length,
        sample_id or get_bam_sample(aln),
        True)

    json.dump(counts, sys.stdout, sort_keys=True)


if __name__ == "__main__":
    main()
