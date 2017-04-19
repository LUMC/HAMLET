#!/usr/bin/env python

"""Snakemake helper for clipping, trimming, and syncing paired-end FASTQ files."""

import argparse
import os
import shutil


def clip_sync_trim(input_r1, input_r2, input_fqc_r1, input_fqc_r2,
                   output_r1, output_r2, contams_fname):
    shutil.copy(input_r1, output_r1)
    shutil.copy(input_r2, output_r2)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--i1", type=str, required=True,
                        help="path to input R1 FASTQ file")
    parser.add_argument("--i2", type=str, required=True,
                        help="path to input R2 FASTQ file")
    parser.add_argument("--fqc1", type=str, required=True,
                        help="path to input FastQC JSON output from R1")
    parser.add_argument("--fqc2", type=str, required=True,
                        help="path to input FastQC JSON output from R2")
    parser.add_argument("--o1", type=str, required=True,
                        help="path to output R1 FASTQ file")
    parser.add_argument("--o2", type=str, required=True,
                        help="path to output R2 FASTQ file")
    parser.add_argument("--contaminants", type=str, required=True,
                        help="path to FastQC contaminants file")

    args = parser.parse_args()

    clip_sync_trim(args.i1, args.i2, args.fqc1, args.fqc2,
                   args.o1, args.o2, args.contaminants)
