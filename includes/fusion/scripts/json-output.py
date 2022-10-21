#!/usr/bin/env python

import json
import sys
import argparse

from pathlib import Path

def fusion_results(fusion_results_dir, intersected):
    frd = Path(fusion_results_dir)
    sf_plot = next(frd.rglob("*star-fusion-circos/*.png"))
    sf_table = next(frd.glob("*.star-fusion"))
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
                sft.append({"name": name, "jr_count": int(jr_count),
                            "sf_count": int(sf_count), "type": fusion_type,})
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


def main(args):
    """ Create json output of fusion results """
    results = fusion_results(args.fusion_results_dir, args.intersected)
    json.dump({"fusion": results}, sys.stdout, sort_keys=True, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--intersected', default=False, action='store_true')
    parser.add_argument('fusion_results_dir')

    args = parser.parse_args()
    main(args)

