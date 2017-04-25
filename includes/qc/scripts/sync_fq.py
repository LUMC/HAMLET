#!/usr/bin/env python

"""
Script for syncing paired-end FASTQs. Adapted from here[1].


[1] https://github.com/martijnvermaat/bio-playground/blob/master/sync-paired-end-reads/sync_paired_end_reads.py
"""

import argparse
import gzip
import json
import re
import sys


_RE_ID = re.compile(rb"[_/][12]$")


def sync_paired_reads(ori, in1, in2, synced1, synced2):

    def next_record(fh):
        return [fh.readline().strip() for i in range(4)]

    def get_header(record):
        return _RE_ID.split(record[0].split(b" ", 1)[0])[0]

    headers = (get_header([line.strip()])
               for idx, line in enumerate(ori) if not (idx % 4))

    n_filtered1 = n_filtered2 = n_kept = 0

    rec1, rec2 = next_record(in1), next_record(in2)

    for header in headers:

        header1, header2 = get_header(rec1), get_header(rec2)

        if header == header1 and header2 != header:
            rec1 = next_record(in1)
            n_filtered1 += 1

        if header == header2 and header1 != header:
            rec2 = next_record(in2)
            n_filtered2 += 1

        if header == header1 == header2:
            synced1.write(b"\n".join(rec1) + b"\n")
            synced2.write(b"\n".join(rec2) + b"\n")
            rec1, rec2 = next_record(in1), next_record(in2)
            n_kept += 1

    return n_filtered1, n_filtered2, n_kept


def make_fh(fname, mode="rb"):
    if fname.endswith(".gz") or fname.endswith(".gzip"):
        return gzip.open(fname, mode=mode)
    return open(fname, mode=mode)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--ori", type=str, required=True,
                        help="path to either R1 or R2 of the FASTQ file"
                             " containing a superset all R1 and R2 records")
    parser.add_argument("--i1", type=str, required=True,
                        help="path to input R1 FASTQ file")
    parser.add_argument("--i2", type=str, required=True,
                        help="path to input R2 FASTQ file")
    parser.add_argument("--o1", type=str, required=True,
                        help="path to output R1 FASTQ file")
    parser.add_argument("--o2", type=str, required=True,
                        help="path to output R2 FASTQ file")
    parser.add_argument("--stats", type=str, required=False,
                        help="path to output stats JSON file")

    args = parser.parse_args()

    (n1, n2, nk) = sync_paired_reads(
        make_fh(args.ori), make_fh(args.i1), make_fh(args.i2),
        make_fh(args.o1, mode="wb"), make_fh(args.o2, mode="wb"))

    statsd = {
        "num_read1_filtered": n1,
        "num_read2_filtered": n1,
        "num_kept": nk,
    }

    stats_fh = open(args.stats, "w") if args.stats is not None else sys.stderr
    json.dump(statsd, stats_fh)
