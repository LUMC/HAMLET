#!/usr/bin/env python

import json
import csv
from pathlib import Path
from collections import defaultdict

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


def add_variant_plots(idm, var_plot_dir):
    idms = set([])
    for item in idm:
        gsym = item["gene_symbol"]
        assert gsym not in idms, item
        idms.add(gsym)

    plots = []
    vpd = Path(var_plot_dir)
    for png in vpd.glob("*.png"):
        stem = png.stem
        _, gene = png.stem.rsplit("_gene_", 1)
        if gene in idms:
            plots.append({"path": str(png.absolute()), "gene": gene})

    return plots

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
            rv[gs].append({k: v for k, v in row.items()})

    return rv


def add_fusion_results(fusion_results_dir):
    frd = Path(fusion_results_dir)
    sf_plot = next(frd.rglob("*star-fusion-circos/*.png"))
    sf_table = next(frd.glob("*.star-fusion"))
    intersected = any(True for _ in frd.glob("*.fusions-combined.svg"))
    rv = {
        "intersected": intersected,
        "plots": {
            "star-fusion": str(sf_plot.resolve())
        },
        "tables": {
            "star-fusion": {"path": str(sf_table.resolve())},
        }
    }

    def parse_top20(tp):
        with open(tp, "r") as src:
            sft = []
            for lineno, line in enumerate(src):
                if lineno == 0:
                    continue
                name, jr_count, sf_count, fusion_type, *_ = line.split("\t")
                sft.append({"name": name, "jr_count": jr_count,
                            "sf_count": sf_count, "type": fusion_type,})
                if lineno > 20:
                    break
        return sft

    rv["tables"]["star-fusion"]["top20"] = \
        parse_top20(rv["tables"]["star-fusion"]["path"])

    if intersected:
        fc_plot = next(frd.rglob("*fusioncatcher-circos/*.png"))
        fc_table = next(frd.glob("*.star-fusion"))
        rv["plots"]["fusioncatcher"] = str(fc_plot.resolve())
        rv["tables"]["fusioncatcher"] = {"path": str(fc_table.resolve())}
        rv["tables"]["fusioncatcher"]["top20"] = \
            parse_top20(rv["tables"]["fusioncatcher"]["path"])

        isect_plot = next(frd.rglob("*sf-isect-circos/*.png"))
        isect_table = next(frd.glob("*.sf-isect"))
        rv["plots"]["intersection"] = str(isect_plot.resolve())
        rv["tables"]["intersection"] = {"path": str(isect_table.resolve())}
        rv["tables"]["intersection"]["top20"] = \
            parse_top20(rv["tables"]["intersection"]["path"])

    return rv


def add_expr_results(exon_ratios_path):
    rv = []
    with open(exon_ratios_path, "r") as src:
        header_cols = next(src).strip().split("\t")
        for line in (l.strip() for l in src):
            d = dict(zip(header_cols, line.split("\t")))
            # Update types
            d["count"] = int(d["count"])
            d["divisor_exp"] = int(d["divisor_exp"])
            d["ratio" ] = float(d["ratio"])
            d["above_threshold"] = d["above_threshold"] == "yes"
            del d["sample_name"]
            rv.append(d)
    return rv


def add_itd_table(csv_fname):
    rv = []
    with open(csv_fname, "r") as src:
        header_cols = next(src).strip().split("\t")
        for line in (l.strip() for l in src):
            rv.append(dict(zip(header_cols, line.split("\t"))))
    return rv


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("id_mappings_path", type=click.File("r"))
@click.argument("var_plot_dir",
                type=click.Path(exists=True, file_okay=False))
@click.argument("var_csv",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("fusion_results_dir",
                type=click.Path(exists=True, file_okay=False))
@click.argument("flt3_csv",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("flt3_plot",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("kmt2a_csv",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("kmt2a_plot",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("exon_ratios_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("seq_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("aln_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("rna_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("insert_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("exon_cov_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("vep_stats_path",
                type=click.Path(exists=True, dir_okay=False))
@click.option("-r", "--run-name", type=str,
              help="Name of the run in which the stats were generated.")
@click.option("-n", "--sample-name", type=str,
              help="Name of the sample from which the stats were generated.")
@click.option("--pipeline-version", type=str,
              help="Version string of the pipeline.")
def main(id_mappings_path, var_plot_dir, var_csv, fusion_results_dir,
         flt3_csv, flt3_plot, kmt2a_csv, kmt2a_plot, exon_ratios_path,
         seq_stats_path, aln_stats_path, rna_stats_path,
         insert_stats_path, exon_cov_stats_path, vep_stats_path,
         run_name, sample_name, pipeline_version):
    """Helper script for combining multiple stats files into one JSON."""
    idm = parse_idm(id_mappings_path)
    combined = {
        "metadata": {
            "pipeline_version": pipeline_version,
            "run_name": run_name,
            "sample_name": sample_name,
            "genes_of_interest": idm,
        },
        "stats": {
            "seq": process_seq_stats(seq_stats_path),
            "aln": process_aln_stats(aln_stats_path),
            "rna": process_rna_stats(rna_stats_path),
            "cov": process_exon_cov_stats(exon_cov_stats_path, idm),
            "ins": process_insert_stats(insert_stats_path),
            "var": process_var_stats(vep_stats_path),
        },
        "results": {
            "var": {"plots": [], "overview": {}},
            "fusion": {},
        },
    }
    combined["results"]["var"]["plots"].extend(
        add_variant_plots(idm, var_plot_dir))
    combined["results"]["var"]["overview"] = add_variant_overview(
        idm, var_csv)
    combined["results"]["fusion"] = add_fusion_results(fusion_results_dir)
    combined["results"]["itd"] = {
        "flt3": {"path": str(Path(flt3_plot).resolve()),
                 "table": add_itd_table(flt3_csv)},
        "kmt2a": {"path": str(Path(kmt2a_plot).resolve()),
                  "table": add_itd_table(kmt2a_csv)},
    }
    combined["results"]["expr"] = add_expr_results(exon_ratios_path)
    combined["stats"] = post_process(combined["stats"])
    print(json.dumps(combined, separators=(",", ":"), sort_keys=True))


if __name__ == "__main__":
    main()

