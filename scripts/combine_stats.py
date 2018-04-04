#!/usr/bin/env python

import json

import click
from crimson import picard, vep


def process_aln_stats(path):
    pd = picard.parse(path)
    raw = next(x for x in pd["metrics"]["contents"]
               if x["CATEGORY"] == "PAIR")

    aln_proper_pairs = (raw["READS_ALIGNED_IN_PAIRS"] -
                        raw["PF_READS_IMPROPER_PAIRS"])

    return {
        "num_aligned_bases": raw["PF_ALIGNED_BASES"],
        "num_aligned_reads": raw["PF_READS_ALIGNED"],
        "num_aligned_reads_proper_pairs": aln_proper_pairs,
        "num_total_reads": raw["TOTAL_READS"],
        "pct_adapter": raw["PCT_ADAPTER"],
        "pct_aligned_reads_from_total": (raw["PF_READS_ALIGNED"] * 100. /
                                         raw["TOTAL_READS"]),
        "pct_aligned_reads_proper_pairs": (aln_proper_pairs * 100. /
                                           raw["PF_READS_ALIGNED"]),
        "pct_chimeras": raw["PCT_CHIMERAS"],
        "rate_indel": raw["PF_INDEL_RATE"],
        "rate_mismatch": raw["PF_MISMATCH_RATE"],
        "strand_balance": raw["STRAND_BALANCE"],
    }


def process_rna_stats(path):
    pd = picard.parse(path)
    raw = pd["metrics"]["contents"]
    cov_sort_key = lambda x: x["normalized_position"]
    return {
        "median_3prime_bias": raw["MEDIAN_3PRIME_BIAS"],
        "median_5prime_bias": raw["MEDIAN_5PRIME_BIAS"],
        "median_5prime_to_3prime_bias": raw["MEDIAN_5PRIME_TO_3PRIME_BIAS"],
        "median_cv_coverage": raw["MEDIAN_CV_COVERAGE"],
        "num_coding_bases": raw["CODING_BASES"],
        "num_intergenic_bases": raw["INTERGENIC_BASES"],
        "num_intronic_bases": raw["INTRONIC_BASES"],
        "num_mrna_bases": raw["CODING_BASES"] + raw["UTR_BASES"],
        "num_ribosomal_bases": (raw["RIBOSOMAL_BASES"]
                                if raw["RIBOSOMAL_BASES"] != ""
                                else None),
        "num_total_bases": raw["PF_BASES"],
        "num_utr_bases": raw["UTR_BASES"],
        "pct_coding_bases": raw["PCT_CODING_BASES"],
        "pct_intergenic_bases": raw["PCT_INTERGENIC_BASES"],
        "pct_intronic_bases": raw["PCT_INTRONIC_BASES"],
        "pct_mrna_bases": raw["PCT_MRNA_BASES"],
        "pct_ribosomal_bases": (raw["RIBOSOMAL_BASES"] * 100. / raw["PF_ALIGNED_BASES"]
                                if raw["RIBOSOMAL_BASES"] != ""
                                else None),
        "pct_utr_bases": raw["PCT_UTR_BASES"],
        "normalized_cov": [item["All_Reads.normalized_coverage"]
                           for item in sorted(pd["histogram"]["contents"],
                                              key=cov_sort_key, reverse=False)]
    }


def process_seq_stats(path):
    with open(path) as src:
        raw = json.load(src)

    def f(rgi):
        return {
            "raw": {
                "num_reads_r1": rgi["raw"]["R1"]["num_seq"],
                "num_reads_r2": rgi["raw"]["R2"]["num_seq"],
                "pct_gc_r1": rgi["raw"]["R1"]["pct_gc"],
                "pct_gc_r2": rgi["raw"]["R2"]["pct_gc"],
            },
            "proc": {
                "num_reads_r1": rgi["proc"]["R1"]["num_seq"],
                "num_reads_r2": rgi["proc"]["R2"]["num_seq"],
                "pct_gc_r1": rgi["proc"]["R1"]["pct_gc"],
                "pct_gc_r2": rgi["proc"]["R2"]["pct_gc"],
            },
            "name": rgi["name"],
        }

    rv = {
        "all_read_groups": {
            "raw": {
                "num_reads_r1": raw["raw"]["R1"]["num_seq"],
                "num_reads_r2": raw["raw"]["R2"]["num_seq"],
            },
            "proc": {
                "num_reads_r1": raw["proc"]["R1"]["num_seq"],
                "num_reads_r2": raw["proc"]["R2"]["num_seq"],
            },
        },
        "per_read_group": [f(item) for item in raw["read_groups"]],

    }
    return rv


def process_insert_stats(path):
    pd = picard.parse(path)
    raw = pd["metrics"]["contents"]
    # Raw can be a list if there are more than one PAIR_ORIENTATION.
    if isinstance(raw, list):
        raw = next(r for r in raw if r["PAIR_ORIENTATION"] == "FR")
    return {
        "median_absolute_deviation": raw["MEDIAN_ABSOLUTE_DEVIATION"],
        "median_insert_size": raw["MEDIAN_INSERT_SIZE"],
        "min_insert_size": raw["MIN_INSERT_SIZE"],
        "max_insert_size": raw["MAX_INSERT_SIZE"],
    }


def process_var_stats(path):
    pd = vep.parse(path)
    return {
        "coding_consequences":
        {k: v for k, v in pd["Coding consequences"].items()},
        "num_deletions": pd["Variant classes"].get("deletion", 0),
        "num_insertions": pd["Variant classes"].get("insertion", 0),
        "num_snvs": pd["Variant classes"].get("SNV", 0),
        "per_chromosome":
        {k: v
         for k, v in pd["Variants by chromosome"].items()
         if k in {str(i) for i in range(1, 23)}.union({"X", "Y", "MT"})},
        "polyphen": {
            "num_benign_variants": pd["PolyPhen summary"].get("benign", 0),
            "num_possibly_damaging_variants":
            pd["PolyPhen summary"].get("possibly_damaging", 0),
            "num_probably_damaging_variants":
            pd["PolyPhen summary"].get("probably_damaging", 0),
            "num_unknown_variants": pd["PolyPhen summary"].get("unknown", 0),
        },
        "sift": {
            "num_deleterious_variants":
            pd["SIFT summary"].get("deleterious", 0),
            "num_tolerated_variants":
            pd["SIFT summary"].get("tolerated", 0),
        },
    }


def post_process(cs):
    cs["aln"]["num_total_bases"] = (cs["rna"]
                                          .pop("num_total_bases"))
    pct_bases_aln = (cs["aln"]["num_aligned_bases"] * 100. /
                     cs["aln"]["num_total_bases"])
    cs["aln"]["pct_aligned_bases_from_total"] = pct_bases_aln
    return cs


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("seq_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("aln_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("rna_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("insert_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("vep_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.option("-n", "--sample-name", type=str,
              help="Name of the sample from which the stats were generated.")
def main(seq_stats_path, aln_stats_path, rna_stats_path, insert_stats_path,
         vep_stats_path, sample_name):
    """Helper script for combining multiple stats files into one JSON."""
    combined = {
        "sample_name": sample_name,
        "stats": {
            "seq": process_seq_stats(seq_stats_path),
            "aln": process_aln_stats(aln_stats_path),
            "rna": process_rna_stats(rna_stats_path),
            "insert": process_insert_stats(insert_stats_path),
            "var": process_var_stats(vep_stats_path),
        },
    }
    combined["stats"] = post_process(combined["stats"])
    print(json.dumps(combined, separators=(",", ":"), sort_keys=True))


if __name__ == "__main__":
    main()

