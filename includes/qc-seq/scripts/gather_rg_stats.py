#!/usr/bin/env python

import json
import sys

import click
from crimson import fastqc


def parse_fastqc_stats(fastqc_dir):
    raw = fastqc.parse(fastqc_dir)
    stats = {
        "pct_gc": raw["Basic Statistics"]["contents"]["%GC"],
        "num_reads": raw["Basic Statistics"]["contents"]["Total Sequences"],
    }
    return stats


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("raw_fastqc_r1_dir",
                type=click.Path(exists=True, file_okay=False))
@click.argument("raw_fastqc_r2_dir",
                type=click.Path(exists=True, file_okay=False))
@click.argument("proc_fastqc_r1_dir",
                type=click.Path(exists=True, file_okay=False))
@click.argument("proc_fastqc_r2_dir",
                type=click.Path(exists=True, file_okay=False))
@click.option("--name", type=str)
def main(raw_fastqc_r1_dir, raw_fastqc_r2_dir,
         proc_fastqc_r1_dir, proc_fastqc_r2_dir, name):

    # Parse all fastq results
    raw_fastq_r1=parse_fastqc_stats(raw_fastqc_r1_dir)
    raw_fastq_r2=parse_fastqc_stats(raw_fastqc_r2_dir)

    proc_fastq_r1=parse_fastqc_stats(proc_fastqc_r1_dir)
    proc_fastq_r2=parse_fastqc_stats(proc_fastqc_r2_dir)

    stats = {
        "raw": {
            "num_reads_r1": raw_fastq_r1["num_reads"],
            "num_reads_r2": raw_fastq_r2["num_reads"],
            "pct_gc_r1": raw_fastq_r1["pct_gc"],
            "pct_gc_r2": raw_fastq_r2["pct_gc"]
        },
        "proc": {
            "num_reads_r1": proc_fastq_r1["num_reads"],
            "num_reads_r2": proc_fastq_r2["num_reads"],
            "pct_gc_r1": proc_fastq_r1["pct_gc"],
            "pct_gc_r2": proc_fastq_r2["pct_gc"]
        },
        "name": name
    }
    json.dump(stats, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
