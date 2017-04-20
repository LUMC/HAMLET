#!/usr/bin/env python

"""Snakemake helper for clipping, trimming, and syncing paired-end FASTQs."""

import argparse
import concurrent.futures
import os
import subprocess
import sys

from crimson import fastqc


# FastQC-sickle encoding name map
# from FastQC source
# (uk.ac.babraham.FastQC.Sequence.QualityEncoding.PhredEncoding.java)
# and sickle source (src/sickle.h)
FQ_SICKLE_ENC_MAP = {
    "Sanger / Illumina 1.9": ("sanger", 33),
    "Illumina <1.3": ("solexa", 59),
    "Illumina 1.3": ("illumina", 64),
    "Illumina 1.5": ("illumina", 64),
}


def parse_contam_file(contam_file, delimiter="\t"):
    """Returns a dictionary of contaminant sequences names and their
    sequences in the given contaminants file."""
    with open(contam_file, "r") as source:
        # read only lines not beginning with "#" and discard empty lines
        lines = filter(None, (line.strip() for line in source if not
                       line.startswith("#")))
        # parse contam seq lines into lists of [id, sequence]
        parsed = (filter(None, line.split(delimiter))
                  for line in lines)
        # and create a dictionary, key=sequence id and value=sequence
        contam_ref = {name: seq for name, seq in parsed}

    return contam_ref


def get_contams_present(fqcd, contamd):
    """Returns the full sequence of found FastQC contaminants."""
    contam_names = {x["Possible Source"]
                    for x in fqcd["Overrepresented sequences"]["contents"]}
    contams_present = {name: seq for name, seq in contamd.items()
                       if any([cid.startswith(name)
                               for cid in contam_names if cid != "No Hit"])}
    return contams_present


def construct_command(in_fname, out_fname, enc, enc_offset, adapters):
    fifos = []
    if adapters:
        cutadapt_stderr = out_fname + ".cutadapt"
        cutadapt_fifop = out_fname + ".cutadapt.fifo"
        os.mkfifo(cutadapt_fifop)
        fifos.append(cutadapt_fifop)

        cutadapt_toks = ["cutadapt"]
        for adapter in adapters:
            cutadapt_toks.extend(["-a", adapter])
        cutadapt_toks.extend(["-o", cutadapt_fifop])
        cutadapt_toks.append(in_fname)

        subprocess.Popen(cutadapt_toks, stderr=open(cutadapt_stderr, "w"))

        in_fname = cutadapt_fifop

    sickle_stdout = out_fname + ".sickle"
    sickle_fifop = out_fname + ".sickle.fifo"
    os.mkfifo(sickle_fifop)
    fifos.append(sickle_fifop)

    subprocess.Popen(
        ["sickle", "se", "-f", in_fname, "-o", sickle_fifop, "-t", "sanger"],
        stdout=open(sickle_stdout, "w"))

    gzip_proc = subprocess.Popen(
        ["gzip", "-c"], stdout=open(out_fname, "w"),
        stdin=open(sickle_fifop, "r"))

    return gzip_proc, fifos


def process_read(in_fname, in_fqc_dir, out_fname, contams_fname):
    fqcd = fastqc.parse(in_fqc_dir)
    contamd = parse_contam_file(contams_fname)

    raw_enc = fqcd["Basic Statistics"]["contents"]["Encoding"]
    enc, enc_offset = FQ_SICKLE_ENC_MAP[raw_enc]
    adapters = get_contams_present(fqcd, contamd)

    cmd, fifos = construct_command(in_fname, out_fname, enc, enc_offset,
                                   adapters)

    try:
        cmd.wait()
    except:
        raise
    finally:
        for fifo in fifos:
            os.unlink(fifo)


def clip_sync_trim(input_r1, input_r2, input_fqc_r1, input_fqc_r2,
                   output_r1, output_r2, contams_fname):

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:

        futures = {
            executor.submit(process_read, i, j, o, contams_fname): (i, j, o)
            for i, j, o in ((input_r1, input_fqc_r1, output_r1),
                            (input_r2, input_fqc_r2, output_r2))}

        for future in concurrent.futures.as_completed(futures):
            in_fname, _, out_fname = futures[future]
            try:
                future.result()
            except Exception as exc:
                print("processing of {!r} generated an exception: {}"
                      .format(in_fname, exc), file=sys.stderr)
                raise
            else:
                print("processed {!r} to {!r}".format(in_fname, out_fname),
                      file=sys.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--i1", type=str, required=True,
                        help="path to input R1 FASTQ file")
    parser.add_argument("--i2", type=str, required=True,
                        help="path to input R2 FASTQ file")
    parser.add_argument("--fqc1", type=str, required=True,
                        help="path to FastQC R1 output directory")
    parser.add_argument("--fqc2", type=str, required=True,
                        help="path to FastQC R2 output directory")
    parser.add_argument("--o1", type=str, required=True,
                        help="path to output R1 FASTQ file")
    parser.add_argument("--o2", type=str, required=True,
                        help="path to output R2 FASTQ file")
    parser.add_argument("--contaminants", type=str, required=True,
                        help="path to FastQC contaminants file")

    args = parser.parse_args()

    clip_sync_trim(args.i1, args.i2, args.fqc1, args.fqc2,
                   args.o1, args.o2, args.contaminants)
