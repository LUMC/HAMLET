#!/usr/bin/env python3

"""
Script for extracting a JSON payload for plotting from VarScan VCF files
annotated with VEP.

Only samples with calls are returned.

If a file containing Ensembl IDs (one per line) is supplied, only records
associated with the IDs are returned.

Requirements:
    * Python 3.x
    * PyVCF 0.6.7 <https://pyvcf.readthedocs.org>
    * Click 5.1 <http://click.pocoo.org/>

Input requirements:
    * VCF output from VarScan (may contain more than one samples).
    * May or may not be gzipped.
    * Must have been annotated using VEP with the `--allele_number` flag.

This tool was tested using VCF files generated using:
    * VarScan v2.3.7 (mpileup2cns).
    * VEP (Ensembl tools 77, GRCh 38 assembly).


Copyright (c) 2015 Leiden University Medical Center <http://lumc.nl>
All rights reserved.

"""

import json
import os
import re
import sys
from functools import partial
from os import path

import click
import vcf
from pybedtools import BedTool, Interval

__author__ = "Wibowo Arindrarto"
__contact__ = "w.arindrarto@lumc.nl"

__all__ = []

# Mapping of VEP consequence to its predicted impact
# Source: http://www.ensembl.org/info/genome/variation/predicted_data.html
VEP_IMPACTS = {
    "TFBS_amplification": "MODIFIER",
    "regulatory_region_amplification": "MODIFIER",
    "5_prime_UTR_variant": "MODIFIER",
    "regulatory_region_ablation": "MODERATE",
    "start_lost": "HIGH",
    "intron_variant": "MODIFIER",
    "inframe_insertion": "HIGH",
    "non_coding_transcript_exon_variant": "MODIFIER",
    "synonymous_variant": "LOW",
    "mature_miRNA_variant": "MODIFIER",
    "splice_donor_variant": "MODERATE",
    "3_prime_UTR_variant": "MODIFIER",
    "feature_truncation": "MODIFIER",
    "TF_binding_site_variant": "MODIFIER",
    "splice_acceptor_variant": "MODERATE",
    "transcript_amplification": "HIGH",
    "upstream_gene_variant": "MODIFIER",
    "stop_lost": "HIGH",
    "stop_retained_variant": "LOW",
    "inframe_deletion": "HIGH",
    "TFBS_ablation": "MODIFIER",
    "stop_gained": "HIGH",
    "regulatory_region_variant": "MODIFIER",
    "incomplete_terminal_codon_variant": "LOW",
    "intergenic_variant": "MODIFIER",
    "downstream_gene_variant": "MODIFIER",
    "splice_region_variant": "LOW",
    "transcript_ablation": "HIGH",
    "protein_altering_variant": "HIGH",
    "frameshift_variant": "HIGH",
    "feature_elongation": "MODIFIER",
    "NMD_transcript_variant": "MODIFIER",
    "coding_sequence_variant": "MODIFIER",
    "missense_variant": "HIGH",
    "non_coding_transcript_variant": "MODIFIER",
}
CSQ_NAME = "CSQ"


def af_below(threshold, key, varscan_d, vep_d):
    af_values = varscan_d.get(key)
    if af_values is None:
        # default to filtering in variants without AFs
        return True
    alleles = [x.upper() for x in varscan_d["ALT"]]
    allele_af_map = {al.upper(): af
                     for al, af in (v.split(":")
                     for v in af_values)}
    afs = [float(allele_af_map.get(x, 0)) for x in alleles]
    if any([af < threshold for af in afs]):
        return True
    return False


def all_af_subpop_below(threshold, varscan_d, vep_d):

    def is_subpop_af_key(key):
        """Returns whether the given VCF record key denotes a subpopulation
        AF value or not."""
        # 1KG_P3_AF denotes the combined 1KG Phase 3 subpopulation AF
        return (key.startswith("1KG_P3_") and key != "1KG_P3_AF") or \
            key == "GONL_AF"

    subpop_afs = {k: v for k, v in varscan_d.items() if is_subpop_af_key(k)}
    if not subpop_afs:
        return True
    alleles = [x.upper() for x in varscan_d["ALT"]]
    all_below_thresh = []
    for subpop, afs in subpop_afs.items():
        has_below_thresh = []
        for af in afs:
            k, af_values = af.split(":")
            k = k.upper()
            if k not in alleles:
                continue
            if float(af_values) < threshold:
                has_below_thresh.append(True)
            else:
                has_below_thresh.append(False)
        if has_below_thresh:
            all_below_thresh.append(any(has_below_thresh))
        else:
            all_below_thresh.append(True)

    return all(all_below_thresh)


af_1kg_below_1pct = partial(af_below, 0.01, "1KG_P3_AF")
af_1kg_below_5pct = partial(af_below, 0.05, "1KG_P3_AF")
all_af_subpop_below_1pct = partial(all_af_subpop_below, 0.01)
all_af_subpop_below_5pct = partial(all_af_subpop_below, 0.05)


def vcfrec2interval(record):
    """Given a VCF record object, return a pybedtools Interval object."""
    # NOTE: we need to do coordinate conversion manually
    return Interval(record.CHROM, record.POS - 1, record.POS)


def make_record_extractor(reader, csq_info_name=CSQ_NAME):
    """Creates a function for extracting the records of the given VCF."""

    # regex for checking whether a string can be converted to int
    is_int = re.compile(r'^([-+]?\d+)L?$')
    # regex for checking whether a string can be converted to float
    is_decimal = re.compile(r'^([-+]?\d*\.?\d+(?:[eE][-+]?[0-9]+)?)$')
    # samples in this VCF file
    samples = reader.samples
    split_attrs = set(["Consequence", "Existing_variation", "TREMBL"])
    # Assumes the VEP formatting is given as the last space-separated
    # string in the VCF header.
    vep_keys = reader.infos[csq_info_name].desc.split(" ")[-1].split("|")

    def convert_token(tok):
        """Given a string token, tries to convert it to int or float
        if possible. Empty strings are converted to None."""
        if isinstance(tok, str):
            if is_int.match(tok):
                return int(tok)
            if is_decimal.match(tok):
                return float(tok)
            if not tok:
                return
        return tok

    def get_s_alleles(sample_call, alleles):
        """Given a sample call and all the alleles present in a VCF record,
        return the alleles called in the sample and its number."""
        phase_char = sample_call.gt_phase_char()
        num_alleles = sample_call.data.GT.split(phase_char)
        allele_nums = [int(n) for n in num_alleles]
        return [alleles.__getitem__(n) for n in allele_nums], allele_nums

    def get_vep_impact(vep_cons):
        """Given a VEP consequence string, return its impact. If there are
        multiple consequences, returns the most severe impact."""
        if "&" in vep_cons:
            vep_conss = vep_cons.split("&")
            impacts = set([get_vep_impact(v) for v in vep_conss])
            if len(impacts) == 1:
                return impacts.pop()
            elif "HIGH" in impacts:
                return "HIGH"
            elif "MODERATE" in impacts:
                return "MODERATE"
            elif "LOW" in impacts:
                return "LOW"
            elif "MODIFIER" in impacts:
                return "MODIFIER"
            assert False
        return VEP_IMPACTS.get(vep_cons, "UNKNOWN")

    def split_if_exists(d, key, split_char="&"):
        if key in d:
            d[key] = d[key].split(split_char)
        return d

    def parse_raw_vep(raw_str, keep_ampersand):
        """Parses the given raw VEP string into a dictionary. If the number
        of fields and values do not match, None is returned."""
        values = [convert_token(x) for x in raw_str.split("|")]
        if len(vep_keys) == len(values):
            res = {k: v for k, v in zip(vep_keys, values) if v is not None}
            assert "impact" not in res
            res["impact"] = get_vep_impact(res["Consequence"])
            if not keep_ampersand:
                for attr in split_attrs:
                    res = split_if_exists(res, attr)
            return res
        raise click.ClickException("Unexpected VEP values in string '{0}'"
                                   "".format(raw_str))

    def extract_sample_data(record, sample, vep_data, hotspot_ivals):
        """Given a record, a sample name, and a parsed VEP annotation, return
        the sample data."""
        alleles = [record.REF] + [str(x) for x in record.ALT]
        call = record.genotype(sample)
        data = call.data

        if data.GT is None or not call.called:
            vep_data = []
            varscan_ok = {}
        else:
            s_alleles, s_allele_nums = get_s_alleles(call, alleles)
            vep_data = [v for v in vep_data
                        if int(v.get("ALLELE_NUM")) in s_allele_nums]
            varscan_data = [(k, convert_token(getattr(data, k)))
                            for k in data._fields]
            genotype = "{0}/{0}".format(s_alleles[0]) if len(s_alleles) == 1 \
                       else "/".join(s_alleles)
            varscan_data += [
                ("CHROM", record.CHROM),
                ("POS", record.POS),
                ("REF", record.REF),
                ("ALT", [str(a) for a in record.ALT]),
                ("alleles", s_alleles),
                ("genotype", genotype)]
            in_hotspot = None
            if hotspot_ivals is not None:
                ival = vcfrec2interval(record)
                in_hotspot = bool(hotspot_ivals.any_hits(ival))
            varscan_data.append(("is_in_hotspot", in_hotspot))
            # custom af keys
            af = {}
            for custom_info_key in ("GONL", "GONL_AF", "P3", "P3_AF",
                                    "P3_AFR_AF", "P3_AMR_AF", "P3_EAS_AF",
                                    "P3_EUR_AF", "P3_SAS_AF"):
                key_val = record.INFO.get(custom_info_key)
                if key_val is not None:
                    key_val = list(map(str, key_val))
                    custom_info_key = custom_info_key.replace("P3", "1KG_P3")
                    key_val = map(lambda x: ":".join(x),
                                  zip([str(x) for x in record.ALT], key_val))
                    af[custom_info_key] = list(key_val)
            for k, v in af.items():
                varscan_data.append((k, v))
            varscan_ok = {k: v for k, v in varscan_data}

            filters = []
            clean_bases_ratio = \
                (varscan_ok["RD"] + varscan_ok["AD"]) / varscan_ok["DP"]
            clean_filter_ok = clean_bases_ratio > 0.2
            if not af_1kg_below_5pct(varscan_ok, vep_data):
                filters.append("1KGAFAtLeast5Pct")
            if not all_af_subpop_below_5pct(varscan_ok, vep_data):
                filters.append("SubpopAFAtLeast5Pct")
            if not clean_filter_ok:
                filters.append("LowQualBases")

            assert "filters" not in varscan_ok
            varscan_ok["filters"] = filters

        return {
            "sample": sample,
            "vep": vep_data,
            "varscan": varscan_ok,
        }

    def extractor(gene_ids, keep_ampersand, filter_goi, hotspot_ivals,
                  record):
        """Function for extracting records into a dictionary."""

        toi_ids = set([toi_id for ids in gene_ids.values()
                       for toi_id in ids])
        onames_key = "Existing_variation"

        def annotate_toi(vep_data):
            if vep_data.get("Feature") in toi_ids:
                assert "is_toi" not in vep_data
                vep_data["is_toi"] = True
                assert "in_cosmic" not in vep_data
                if not keep_ampersand:
                    has_cosmic = any(["COSM" in name for name
                                     in vep_data.get(onames_key, [])])
                    vep_data["in_cosmic"] = has_cosmic
                else:
                    has_cosmic = "COSM" in vep_data.get(onames_key, "")
                    vep_data["in_cosmic"] = has_cosmic
            return vep_data

        csq_values = record.INFO[csq_info_name]
        vep_data = [parse_raw_vep(x, keep_ampersand) for x in csq_values]
        # select only variants affecting genes of interest
        if filter_goi:
            goi_data = [x for x in vep_data if x.get("Gene") in gene_ids]
        else:
            goi_data = [x for x in vep_data]
        toi_data = list(map(annotate_toi, goi_data))
        sample_data = [extract_sample_data(record, k, toi_data, hotspot_ivals)
                       for k in samples]
        return sample_data

    return extractor


def group(extracted_iter, sample_names, gene_ids):
    """Given the raw dictionary results, group into per-sample,
    per-gene dictionary."""
    samples = {s: {} for s in sample_names}
    dns = ("vep", "varscan")
    for lined in extracted_iter:
        for sampled in lined:
            sample = sampled["sample"]
            entry = {dn: sampled[dn] for dn in dns}
            # select only variants with VEP annotation
            if entry["vep"]:
                for gene_id in gene_ids:
                    if gene_id not in samples[sample]:
                        samples[sample][gene_id] = []
                    gene_varscan = entry["varscan"]
                    gene_vep = [x for x in entry["vep"]
                                if x["Gene"] == gene_id]
                    if gene_vep:
                        gene_entry = {
                            "vep": gene_vep,
                            "varscan": gene_varscan,
                            "sample": sample,
                            "gene": gene_id
                        }
                        samples[sample][gene_id].append(gene_entry)
    return samples


def parse_id_file(fh):
    """Parses the given ID file handle into a dictionary between ENSG IDs and
    ENST IDs."""
    # discard header line
    fh.readline()
    id_mapping = {gid: tid.split(",")
                  for gid, _, tid in (line.strip().split("\t") for line in fh)}
    return id_mapping


@click.command()
@click.argument("id_file",
                type=click.File())
@click.argument("input_vcf",
                type=click.Path(dir_okay=False))
@click.option("--keep-amp", default=False, is_flag=True,
              help="Whether to keep ampersand-separated VEP values as strings"
                   "or split them into a list.")
@click.option("--hotspots", default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help="Path to a BED file containing hotspots region. The "
                   "regions will be annotated in the output JSON file.")
@click.option("--sample-id", type=str,
              help="Set VCF sample name to the given value. If there are "
                   "more than one samples in the VCF, this flag is ignored.")
# TODO: add option for pretty output (default now is compact)
def main(id_file, input_vcf, keep_amp, hotspots, sample_id):
    if input_vcf == "-":
        reader = vcf.Reader(sys.stdin)
    elif not path.exists(input_vcf):
        raise click.BadParameter("Input file not found.")
    elif not os.access(input_vcf, os.R_OK):
        raise click.BadParameter("Input file can not be read.")
    else:
        reader = vcf.Reader(filename=input_vcf)

    if sample_id is not None:
        if len(reader.samples) <= 1:
            reader.samples = [sample_id]
            reader._sample_indexes = {sample_id: 0}
        else:
            msg = "Ignoring `--sample-id` flag since VCF contains " \
                "more than one sample."
            print(msg, file=sys.stderr)

    filter_goi = True
    gene_ids = parse_id_file(id_file)
    hotspot_intervals = BedTool(hotspots).as_intervalfile() if hotspots \
        is not None else None

    def gene_check(rec):
        for csq in rec.INFO["CSQ"]:
            gene = csq.split("|")[1]
            if gene in gene_ids:
                return True
        return False

    f = partial(make_record_extractor(reader), gene_ids, keep_amp, filter_goi,
                hotspot_intervals)
    raw_items = (f(rec) for rec in reader if gene_check(rec))
    samples = group(raw_items, reader.samples, gene_ids)
    meta = {
        "format_keys": list(reader.formats.keys()),
        "vep_keys": [
            x for x in reader.infos[CSQ_NAME].desc.split(" ")[-1].split("|")
        ]
    }
    json.dump({"_meta": meta, "samples": samples}, sys.stdout, sort_keys=True,
              separators=(",", ":"))


if __name__ == "__main__":

    main.__doc__ = __doc__
    main()
