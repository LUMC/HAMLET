#!/usr/bin/env python3

import argparse
import functools
import json
import os


def get_qc_config():
    return {"forward_adapter": "AGATCGGAAGAG", "reverse_adapter": "AGATCGGAAGAG"}


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
    return {
        "genome_star_fusion_lib": os.path.join(
            dirname,
            "GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
        ),
        "transcripts_bed": os.path.join(dirname, "transcripts.bed"),
        "fusioncatcher_data": os.path.join(dirname, "fusioncatcher/current"),
    }


def get_snv_indels_config(dirname):
    join = functools.partial(os.path.join, dirname)
    return {
        "genome_fasta": join("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"),
        "genome_fai": join("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"),
        "genome_dict": join("GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"),
        "star_index": join("star-index"),
        "ref_id_mapping": join("id_mappings.tsv"),
        "rrna_refflat": join("ucsc_rrna.refFlat"),
        "bed_variant_hotspots": join("hotspots_genome.bed"),
        "bed_variant_call_regions": join("call_regions.bed"),
        "gtf": join("Homo_sapiens.GRCh38.104.chr.gtf"),
        "annotation_refflat": join("ucsc_gencode.refFlat"),
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
