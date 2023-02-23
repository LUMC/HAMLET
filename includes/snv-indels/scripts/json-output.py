#!/usr/bin/env python

import argparse
import json
import csv
from pathlib import Path
from collections import defaultdict

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

    # If there are no variants, insert an empty defaultdict so that any queried
    # value returns 0
    if "Variants by chromosome" not in pd:
        pd["Variants by chromosome"] = defaultdict(int)

    return {
        "coding_consequences": {k: v for k, v in pd["Coding consequences"].items()},
        "num_deletions": pd["Variant classes"].get("deletion", 0),
        "num_insertions": pd["Variant classes"].get("insertion", 0),
        "num_snvs": pd["Variant classes"].get("SNV", 0),
        "per_chromosome": {k: v
                             for k, v in pd["Variants by chromosome"].items()
                             if k in {str(i) for i in range(1, 23)}.union({"X", "Y", "MT"})},
        "polyphen": {
            "num_benign_variants": pd["PolyPhen summary"].get("benign", 0),
            "num_possibly_damaging_variants":
            pd["PolyPhen summary"].get("possibly damaging", 0),
            "num_probably_damaging_variants":
            pd["PolyPhen summary"].get("probably damaging", 0),
            "num_unknown_variants": pd["PolyPhen summary"].get("unknown", 0),
        },
        "sift": {
            "num_deleterious_variants":
            pd["SIFT summary"].get("deleterious", 0),
            "num_tolerated_variants":
            pd["SIFT summary"].get("tolerated", 0),
        },
    }


def process_exon_cov_stats(path, idm):
    with open(path, "r") as src:
        raw = json.load(src)

    tid_map = {item["gene_symbol"]: set(item["transcript_ids"])
               for item in idm}

    tempd = {}
    for val in raw.values():
        gene = val.pop("gx")
        if gene not in tempd:
            tempd[gene] = {}
        trx = val.pop("trx")
        if trx not in tid_map[gene]:
            continue
        if trx not in tempd[gene]:
            tempd[gene][trx] = {}
        exon = int(val["exon_num"])
        assert exon not in tempd[gene][trx], f"{gene}:{trx}:{exon}"
        tempd[gene][trx][exon] = val

    return {gid: {tid: sorted([exn for exn in tv.values()],
                               key=lambda x: int(x["exon_num"]))
                  for tid, tv in gv.items()}
            for gid, gv in tempd.items()}


def post_process(cs):
    cs["aln"]["num_total_bases"] = (cs["rna"]
                                          .pop("num_total_bases"))
    pct_bases_aln = (cs["aln"]["num_aligned_bases"] * 100. /
                     cs["aln"]["num_total_bases"])
    cs["aln"]["pct_aligned_bases_from_total"] = pct_bases_aln
    return cs


def parse_idm(fh):
    idms = []
    for lineno, line in enumerate(fh):
        if lineno == 0:
            continue
        gid, gsym, raw_tids = line.strip().split("\t")
        tids = raw_tids.split(",")
        idms.append({"gene_id": gid, "gene_symbol": gsym,
                     "transcript_ids": tids})
    return idms


def add_variant_overview(idm, fn_csv):
    idms = set([])
    for item in idm:
        gsym = item["gene_symbol"]
        assert gsym not in idms, item
        idms.add(gsym)

    rv = {}
    with open(fn_csv, "r") as src:
        reader = csv.DictReader(src)
        for row in reader:
            gs = row["gene_symbol"]
            if gs not in idms:
                continue
            if gs not in rv:
                rv[gs] = []

            # Fields to include in the variant overview
            report_fields = {
                "CHROM", "POS", "REF", "genotype", "HGVSc", "HGVSp",
                "Existing_variation", "FREQ", "is_in_hotspot", "PVAL", "RD",
                "AD", "DP"
            }
            # Extract relevant fields
            d = {k: v for k, v in row.items() if k in report_fields }

            # Replace NA values with an empty string
            d = {k: v if v != "NA" else "" for k, v in d.items()}

            # Update field types
            d["FREQ"] = float(d["FREQ"][:-1]) # Cut off %
            d["is_in_hotspot"] = d["is_in_hotspot"] == "yes"
            d["PVAL"] = float(d["PVAL"])

            # Convert to int
            for field in ["POS", "RD", "AD", "DP"]:
                d[field] = int(d[field])
            if not d["Existing_variation"]:
                d["Existing_variation"] = list()
            else:
                d["Existing_variation"] = d["Existing_variation"].split(",")

            rv[gs].append(d)

    return rv


def main(id_mappings_path, var_csv,
         aln_stats_path, rna_stats_path,
         insert_stats_path, exon_cov_stats_path, vep_stats_path):
    """Helper script for combining multiple stats files into one JSON."""
    with open(id_mappings_path) as fin:
        idm = parse_idm(fin)
    combined = {
        "snv_indels": {
            "genes": {},
            "stats": {
                "aln": process_aln_stats(aln_stats_path),
                "rna": process_rna_stats(rna_stats_path),
                "cov": process_exon_cov_stats(exon_cov_stats_path, idm),
                "ins": process_insert_stats(insert_stats_path),
                "var": process_var_stats(vep_stats_path)
            }
        },
    }

    combined["snv_indels"]["genes"] = add_variant_overview(
        idm, var_csv)
    combined["snv_indels"]["stats"] = post_process(combined["snv_indels"]["stats"])
    print(json.dumps(combined, sort_keys=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("id_mappings_path")
    parser.add_argument("var_csv")
    parser.add_argument("aln_stats_path")
    parser.add_argument("rna_stats_path")
    parser.add_argument("insert_stats_path")
    parser.add_argument("exon_cov_stats_path")
    parser.add_argument("vep_stats_path")

    args = parser.parse_args()
    main(
        args.id_mappings_path,
        args.var_csv,
        args.aln_stats_path,
        args.rna_stats_path,
        args.insert_stats_path,
        args.exon_cov_stats_path,args.vep_stats_path
    )

