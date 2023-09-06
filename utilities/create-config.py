#!/usr/bin/env python3

import argparse
import functools
import json
import os


fusion_partners = ['ABL', 'AFDN', 'AFF1', 'BCOR', 'BCR', 'CBFA2T3', 'CBFB',
        'CREBBP', 'DEK', 'ELL', 'ERG', 'ETV6', 'FIP1L1', 'FUS', 'GATA2',
        'GLIS2', 'IRF2BP2', 'KAT6A', 'KMT2A', 'MECOM', 'MLF1', 'MLLT1',
        'MLLT10', 'MLLT3', 'MNX1', 'MRTF1', 'MYC', 'MYH11', 'NPM1', 'NUP214',
        'NUP98', 'PICALM', 'PML', 'PRDM16', 'RARA', 'RBM15', 'RPN1',
        'RUNX1', 'RUNX1T1', 'STAT3B', 'STAT5B', 'TBL1XR1', 'TET1', 'ZBTB16',
        'KMD5A', 'NSD1']

def get_qc_config():
    return {"forward_adapter": "AGATCGGAAGAG", "reverse_adapter": "AGATCGGAAGAG"}


def get_reference(dirname):
    return os.path.join(dirname, "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")


def get_gtf(dirname):
    return os.path.join(dirname, "Homo_sapiens.GRCh38.104.chr.gtf")


def get_itd_config(dirname):
    return {
        "flt3_fasta": os.path.join(dirname, "flt3-001.fa"),
        "flt3_name": "FLT3-001",
        "flt3_start": 1787,
        "flt3_end": 2024,
        "kmt2a_fasta": os.path.join(dirname, "kmt2a-213.fa"),
        "kmt2a_name": "KMT2A-213",
        "kmt2a_start": 406,
        "kmt2a_end": 4769,
    }


def get_fusion_config(dirname):
    join = functools.partial(os.path.join, dirname)
    return {
        "genome_fasta": get_reference(dirname),
        "gtf": get_gtf(dirname),
        "blacklist": join("arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"),
        "cytobands": join("arriba/cytobands_hg38_GRCh38_v2.4.0.tsv"),
        "known_fusions": join("arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"),
        "protein_domains": join("arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3"),
        "fusion-partners": fusion_partners,
    }


def get_snv_indels_config(dirname):
    join = functools.partial(os.path.join, dirname)
    return {
        "genome_fasta": get_reference(dirname),
        "genome_fai": join("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"),
        "genome_dict": join("GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"),
        "star_index": join("star-index"),
        "ref_id_mapping": join("id_mappings.tsv"),
        "rrna_refflat": join("ucsc_rrna.refFlat"),
        "bed_variant_hotspots": join("hotspots_genome.bed"),
        "bed_variant_call_regions": join("call_regions.bed"),
        "gtf": get_gtf(dirname),
        "annotation_refflat": join("ucsc_gencode.refFlat"),
        "blacklist": join("blacklist.txt"),
        "vep_cache": dirname,
        "vep_include_consequence": [
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "inframe_insertion",
            "inframe_deletion",
            "protein_altering_variant",
            "missense_variant",
        ],
    }


def main(dirname):
    # Get the absolute path to the root reference folder
    dirname = os.path.abspath(dirname)

    config = dict()
    config["qc-seq"] = get_qc_config()
    config["itd"] = get_itd_config(dirname)
    config["fusion"] = get_fusion_config(dirname)
    config["snv-indels"] = get_snv_indels_config(dirname)

    print(json.dumps(config, indent=True, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("reference_dir", help="Path to the reference folder")

    args = parser.parse_args()
    main(args.reference_dir)
