#!/usr/bin/env python

import gzip
import locale
import re
import sys
from collections import namedtuple
from contextlib import contextmanager

import click


locale.setlocale(locale.LC_ALL, "en_US.utf8")

_ATTR_RE = re.compile(r' ?(.+) "(.+)"')


class GTFRecord(namedtuple("GTFRecord", ["raw"])):

    @property
    def columns(self):
        if not hasattr(self, "_columns"):
            self._columns = self.raw.split("\t")
        return self._columns

    @property
    def chrom(self):
        return self.columns[0]

    @property
    def feature(self):
        return self.columns[2]

    @property
    def start(self):
        return int(self.columns[3]) - 1

    @property
    def end(self):
        return int(self.columns[4])

    @property
    def strand(self):
        return self.columns[6]

    @property
    def attributes(self):
        if not hasattr(self, "_attributes"):
            self._attributes = [
                _ATTR_RE.match(attr).groups()
                for attr in filter(None, self.columns[-1].split(";"))]
        return self._attributes

    @property
    def gene_id(self):
        try:
            return [x[1] for x in self.attributes if x[0] == "gene_id"][0]
        except IndexError:
            return

    @property
    def gene_name(self):
        try:
            return [x[1] for x in self.attributes if x[0] == "gene_name"][0]
        except IndexError:
            return

    @property
    def chrom_is_canonical(self):
        try:
            int(self.chrom)
            return True
        except ValueError:
            return self.chrom.lower in {"x", "y"}

    def as_bed_line(self):
        bed_columns = (
            "chr" + self.chrom, self.start, self.end, self.gene_name or "-", 0,
            self.strand
        )
        return "\t".join([str(x) for x in bed_columns]) + "\n"

    def as_gtf_line(self):
        return "chr" + self.raw + "\n"


@contextmanager
def bopen(fname):
    lname = fname.lower()
    if lname.endswith("gz") or lname.endswith("gzip"):
        yield gzip.open(fname)
    else:
        yield open(fname, "rb")


def parse_id_mappings(fname, encoding):
    mappings = {}
    with bopen(fname) as src:
        for line in src:
            if line.startswith(b"#"):
                continue
            goi, goi_symbol, toi = line.decode(encoding).split("\t")
            tois = toi.split(",")
            if goi in mappings:
                raise ValueError("Gene ID {} can not be listed twice."
                                 .format(goi))
            mappings[goi] = {"symbol": goi_symbol, "tois": tois}
    return mappings


def gtf_line_is_nonrecord(line):
    return line.startswith("#") \
        or line.startswith("track ") \
        or line.startswith("browser ")


def parse_gtf(fname, encoding):
    with bopen(fname) as src:
        for line in src:
            line = line.decode(encoding)
            # skip comment lines and track lines
            if gtf_line_is_nonrecord(line):
                continue
            yield GTFRecord(line.strip())


def write_gtf_and_bed(goi, in_name_gtf, out_name_gtf, out_name_bed, encoding):
    with open(out_name_gtf, "w") as target_gtf, \
            open(out_name_bed, "w") as target_bed:

        for i, record in enumerate(parse_gtf(in_name_gtf, encoding), start=1):

            if i == 10 or i % 50000 == 0:
                fmted = locale.format("%d", i, grouping=True)
                print("Processed {} records ...".format(fmted),
                      file=sys.stderr)

            if record.gene_id not in goi:
                continue

            if record.feature.lower() == "gene":
                target_bed.write(record.as_bed_line())

            target_gtf.write(record.as_gtf_line())

        fmted = locale.format("%d", i, grouping=True)
        print("Processed a total of {} records ...".format(fmted),
              file=sys.stderr)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("id_mappings", type=click.Path(exists=True, dir_okay=False))
@click.argument("input_gtf", type=click.Path(exists=True, dir_okay=False))
@click.argument("output_gtf", type=click.Path(dir_okay=False))
@click.argument("output_bed", type=click.Path(dir_okay=False))
@click.option("--encoding", type=str, default="utf-8")
def main(id_mappings, input_gtf, output_gtf, output_bed, encoding):
    goi = parse_id_mappings(id_mappings, encoding)
    write_gtf_and_bed(goi, input_gtf, output_gtf, output_bed, encoding)


if __name__ == "__main__":
    main()
