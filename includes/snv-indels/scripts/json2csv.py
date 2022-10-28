#!/usr/bin/env python3

"""
Script for converting Hamlet VCF JSON payload to a CSV file.


Copyright (c) 2016 Leiden University Medical Center <http://lumc.nl>
All rights reserved.

"""

import argparse
import json
import os
import sys
from os import path

__author__ = "Wibowo Arindrarto"
__contact__ = "w.arindrarto@lumc.nl"

__all__ = []

# TODO: use the csv module


def make_csv(grouped, meta, fallback="NA", csq_info_name="CSQ",
             with_header=True, toi_filter=True, af_filter=False,
             impact_filter=False, dp_ratio_filter=False):
    vep_keys = ["is_toi", "impact", "in_cosmic"] + meta["vep_keys"]
    var_keys = ["CHROM", "POS", "REF", "alleles", "genotype"] + \
        ["is_in_hotspot", "filters", "GONL", "GONL_AF"] + \
        ["1KG_P3" + x for x in ["", "_AF", "_AFR_AF", "_AMR_AF",
                                "_EAS_AF", "_EUR_AF", "_SAS_AF"]] + \
        meta["format_keys"]
    header_cols = [
        '"{}"'.format(x) if isinstance(x, str) and "," in x else x
        for x in ["sample_id", "gene_symbol", "gene_id"] + var_keys + vep_keys]
    header_fmt = "{" + "},{".join(header_cols) + "}"

    def vep_ok(vep):
        cond = True
        if toi_filter:
            cond = cond and vep.get("is_toi", False)
        if impact_filter:
            cond = cond and vep.get("impact", "").upper() == "HIGH"
        return cond

    if with_header:
        yield ",".join(header_cols)

    for sample_id, sampled in grouped.items():
        for gene_id, gened in sampled.items():
            for variant in gened:
                varscan = variant["varscan"]
                if af_filter and \
                        "SubpopAFAtLeast5Pct" in varscan.get("filters", []):
                    continue
                if dp_ratio_filter and \
                        "LowQualBases" in varscan.get("filters", []):
                    continue
                for aff in (x for x in variant["vep"] if vep_ok(x)):
                    rowd = [("sample_id", sample_id),
                            ("gene_symbol", aff.get("SYMBOL")),
                            ("gene_id", gene_id)]
                    rowd += [(x, varscan.get(x)) for x in var_keys]
                    rowd += [(x, aff.get(x)) for x in vep_keys]
                    to_print = []
                    for k, v in rowd:
                        if isinstance(v, list):
                            to_print.append((k,
                                             ",".join(v) if v else fallback))
                        elif isinstance(v, bool):
                            to_print.append((k, "yes" if v else "no"))
                        elif v is None:
                            to_print.append((k, fallback))
                        else:
                            to_print.append((k, v))
                    row = header_fmt.format(**{
                        k: '"{}"'.format(v)
                           if isinstance(v, str) and "," in v else v
                        for k, v in to_print})
                    yield row


def main(input_json, hi):
    if input_json == "-":
        payload = json.load(sys.stdin)
    else:
        with open(input_json, "r") as src:
            payload = json.load(src)

    payload_meta = payload["_meta"]
    payload_samples = payload["samples"]

    if not payload_samples:
        raise RuntimeError("Can not process JSON file without samples.")
    elif len(payload_samples) > 1:
        msg = "Can not process JSON with more than one samples."
        raise RuntimeError(msg)

    sample_id, sampled = next(iter(payload_samples.items()))
    items = (x for x in sorted(sampled.items(), key=lambda x: x[0]) if x[1])
    for idx, ginfo in enumerate(items, start=1):
        gene_id, gened = ginfo
        faked = {sample_id: {gene_id: gened}}
        for row in make_csv(faked, payload_meta, toi_filter=True,
                            af_filter=hi, impact_filter=hi,
                            dp_ratio_filter=hi, with_header=(idx == 1)):
            print(row, file=sys.stdout)


if __name__ == "__main__":
    main.__doc__ = __doc__
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json")
    parser.add_argument("--hi", default=False, action='store_true')

    args = parser.parse_args()
    main(args.input_json, args.hi)
