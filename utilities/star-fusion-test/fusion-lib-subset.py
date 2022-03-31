#!/usr/bin/env python3

import argparse
import logging

import gzip

logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
)


def annot_to_dict(line):
    """ Parse the annotation field from gtf into a python dictionary """
    d = dict()
    for entry in line.split('; ')[:-1]:
        key, value = entry.split(' ')
        d[key] = value.replace('"', '')
    return d


def get_transcripts(filename):
    transcripts = set()
    with open(filename) as fin:
        for line in fin:
            annot = line.strip().split('\t')[8]
            d = annot_to_dict(annot)
            transcripts.add(d['gene_name'])
    return transcripts


def main(args):
    transcripts = get_transcripts(args.gtf)

    with gzip.open(args.input, 'rt') as fin:
        with gzip.open(args.output, 'wt') as fout:
            for line in fin:
                genes = line.strip().split('\t')[0]
                if '--' not in genes:
                    continue

                gene1, gene2 = genes.split('--')
                if gene1 in transcripts and gene2 in transcripts:
                    print(line, end='', file=fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--gtf', required=True, help='GTF file containing transcripts of interest')
    parser.add_argument('--input', required=True, help='Input fusion_lib from fusioncatcher')
    parser.add_argument('--output', required=True, help='Output fusion lib to write')
    arguments = parser.parse_args()
    main(arguments)
