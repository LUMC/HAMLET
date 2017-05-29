#!/usr/bin/env python

"""Snakemake helper for clipping, trimming, and syncing paired-end FASTQs."""

import argparse
import os
import subprocess
from os import path

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


FQ_EXTS = (
    ".fq.gzip", ".fastq.gzip",
    ".fq.gz", ".fastq.gz",
    ".fq", ".fastq",
)


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


def construct_command(in_fname, out_fname_unsynced, enc, enc_offset, adapters):
    os.mkfifo(out_fname_unsynced)
    fifos = []
    if adapters:
        cutadapt_stderr = out_fname_unsynced + ".cutadapt"
        cutadapt_fifop = out_fname_unsynced + ".cutadapt.fifo"
        os.mkfifo(cutadapt_fifop)
        fifos.append(cutadapt_fifop)

        cutadapt_toks = ["cutadapt", "-m", "20"]
        for adapter in adapters:
            cutadapt_toks.extend(["-a", adapter])
        cutadapt_toks.extend(["-o", cutadapt_fifop])
        cutadapt_toks.append(in_fname)

        subprocess.Popen(cutadapt_toks, stderr=open(cutadapt_stderr, "w"))

        in_fname = cutadapt_fifop

    sickle_stdout = out_fname_unsynced + ".sickle"

    subprocess.Popen(
        ["sickle", "se", "-f", in_fname, "-o", out_fname_unsynced,
         "-t", "sanger"],
        stdout=open(sickle_stdout, "w"))

    return fifos


def process_read(in_fname, in_fqc_dir, out_fname, contams_fname):
    fqcd = fastqc.parse(in_fqc_dir)
    contamd = parse_contam_file(contams_fname)

    raw_enc = fqcd["Basic Statistics"]["contents"]["Encoding"]
    enc, enc_offset = FQ_SICKLE_ENC_MAP[raw_enc]
    adapters = get_contams_present(fqcd, contamd)

    return construct_command(in_fname, out_fname, enc, enc_offset, adapters)


def splitext_fq(fname):
    fnamel = fname.lower()
    for pext in FQ_EXTS:
        if fnamel.endswith(pext):
            lpext = len(pext)
            return fname[:-lpext], fname[-lpext:]
    return fname


def mark_unsynced(fname):
    dname, bname_fq = path.dirname(fname), path.basename(fname)
    bname, _ = splitext_fq(bname_fq)
    return path.join(dname, bname + ".unsynced.fq")


def clip_sync_trim(input_r1, input_r2, input_fqc_r1, input_fqc_r2,
                   output_r1, output_r2, contams_fname, sync_scr_fname,
                   stats_fname):

    output_r1_unsynced = mark_unsynced(output_r1)
    output_r2_unsynced = mark_unsynced(output_r2)

    fifos = [output_r1_unsynced, output_r2_unsynced]
    for i, j, o in ((input_r1, input_fqc_r1, output_r1_unsynced),
                    (input_r2, input_fqc_r2, output_r2_unsynced)):
        fifo_handles = process_read(i, j, o, contams_fname)
        fifos.extend(fifo_handles)

    sync_toks = [
        "python", sync_scr_fname, "--ori", input_r1,
        "--i1", output_r1_unsynced, "--i2", output_r2_unsynced,
        "--o1", output_r1, "--o2", output_r2]
    if stats_fname is not None:
        sync_toks.extend(["--stats", stats_fname])

    sync_proc = subprocess.Popen(sync_toks)

    try:
        sync_proc.wait()
    except:
        raise
    finally:
        for temp_out in fifos:
            os.unlink(temp_out)


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
    parser.add_argument("--sync-scr", type=str, required=True,
                        help="path to sync script")
    parser.add_argument("--stats", type=str, required=False,
                        help="path to output stats JSON file")

    args = parser.parse_args()

    clip_sync_trim(args.i1, args.i2, args.fqc1, args.fqc2,
                   args.o1, args.o2, args.contaminants,
                   args.sync_scr, args.stats)
