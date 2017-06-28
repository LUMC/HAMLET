#!/usr/bin/env python

import json
import math

import click
import pandas as pd
import matplotlib; matplotlib.use("Agg")  # noqa
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.gridspec as gs
from matplotlib.path import Path as pth
from matplotlib.ticker import FuncFormatter


def to_percent(y, position):
    """Input function for ``FuncFormatter`` for percentage formatting."""
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)
    # The percent symbol needs escaping in latex
    if plt.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'


def calc_sc_ratio(sc_count, pileup_count):
    """Given the soft clip count and pileup count of a position,
    return the soft clip count ratio."""
    if sc_count is None:
        return math.nan
    return sc_count / (sc_count + pileup_count)


def calc_insert_ratio(insert_count, pileup_count):
    """Given the insert count and pileup count of a position,
    return the insert count ratio."""
    try:
        return insert_count / pileup_count
    except ZeroDivisionError:
        return math.nan


def make_sc_sample_df(countd, region, sc_bg=None, min_scl_count=2):
    """Converts the given count dictionary into a data frame of soft clip
    event counts in the given region.

    The returned data frame has the following columns:
        * sample: Name of the sample.
        * pos: Position of the event.
        * pileup_count: Number of reads mapped to the position.
        * sc_count: Number of soft clip events on the position.
        * sc_ratio: Ratio of soft clip events to the total
                    reads (including the soft clip events) on the position.
        * asc: List of position and counts where a a soft clip from the
               position can map.

    In addition, the returned data frame is indexed on the `position` column.

    :param dict countd: Dictionary of count data. See script documentation
                        for more information.
    :param range region: Range over which the count data will be extracted.
    :param int min_count: The minimum number of counts an event must have
                          for it to be considered. Events whose count is
                          lower than this will be ignored. Default: 2.
    :returns: A pandas DataFrame object.

    """
    sample = countd["sampleName"]
    raw_pileups = {item["pos"]: item["count"] for item in countd["pileups"]}
    raw_scs = {item["pos"]:
               {"count": item["count"], "altPosCount": item["altPosCount"]}
               for item in countd["scs"]}

    dps = []
    for p in region:
        pileup_count = raw_pileups.get(p, 0)
        sc_count = raw_scs.get(p, {}).get("count", None)
        if sc_bg is not None:
            try:
                bg = sc_bg.loc[p].threshold
            except KeyError:
                sc_ratio = math.nan
            else:
                sc_ratio = calc_sc_ratio(sc_count, pileup_count)
                if sc_ratio * 100.0 <= bg:
                    sc_ratio = math.nan
        else:
            sc_ratio = calc_sc_ratio(sc_count, pileup_count)
        dp = {
            "sample": sample,
            "pos": p,
            "pileup_count": raw_pileups.get(p, 0),
            "sc_count": raw_scs.get(p, {}).get("count"),
            "sc_ratio": sc_ratio,
            "asc": [item for item in raw_scs.get(p, {}).get("altPosCount", list())
                    if item["count"] >= min_scl_count],
        }
        dps.append(dp)

    df = pd.DataFrame(dps)
    df.set_index("pos", inplace=True, drop=False)
    return df


def make_insert_sample_df(countd, region, min_count=2):
    """Converts the given count dictionary into a data frame of insertion
    event counts in the given region.

    The returned data frame has the following columns:
        * sample: Name of the sample.
        * pos: Position of the event.
        * pileup_count: Number of reads mapped to the position.
        * insert_count: Number of insertion events on the position.
        * insert_ratio: Ratio of insertion events to the total reads
                        on the position.
        * insert_seq: Sequence of the insertion.

    :param dict countd: Dictionary of count data. See script documentation
                        for more information.
    :param range region: Range over which the count data will be extracted.
    :param int min_count: The minimum number of counts an event must have
                          for it to be considered. Events whose count is
                          lower than this will be ignored. Default: 2.
    :returns: A pandas DataFrame object.

    """
    sample = countd["sampleName"]
    raw_pileups = {x["pos"]: x["count"] for x in countd["pileups"]}
    inserts = [item for item in countd["inserts"]
               if region.start <= item["pos"] < region.stop and
               item["count"] >= min_count]
    insert_ratios = [calc_insert_ratio(item["count"], raw_pileups[item["pos"]])
                     for item in inserts]
    df = pd.DataFrame({
        "sample": [sample] * len(inserts),
        "pos": [item["pos"] for item in inserts],
        "insert_count": [item["count"] for item in inserts],
        "insert_ratio": insert_ratios,
        "insert_seq": [item["seq"] for item in inserts],
        "pileup_count": [raw_pileups[item["pos"]] for item in inserts],
    })
    df.sort_values(by=["insert_ratio"], inplace=True, ascending=False)
    return df


def adjust_pair_fuzz(pair, scd, fuzziness):
    """Adjusts the second element of a pair to point to a position within
    a given fuzziness region where there may be a reciprocating soft clip
    alignment.

    :param (int, int) pair: Pair of coordinate whose second element will
                            be adjusted.
    :param pd.DataFrame scd: Pandas data frame of soft clip event counts.
    :param int fuzziness: Length of 5' and 3' extension between which
                          candidate for pair will be searched.
    :returns: The coordinate whose second element may be adjusted.

    """
    # Two coordinates: pair[0] (soft clip), pair[1] (soft clip pair)
    # Strategy:
    #     * Get count_asc, which is the asc_count of pair[1]
    #     * In the fuzz region (pair[1] +- fuzziness), find pos whose count is
    #       closest to count_asc
    #     * That is the adjusted pair[1]
    count_asc, = [item["count"]
                  for item in scd.loc[pair[0]]["asc"]
                  if item["pos"] == pair[1]]
    fuzz_origin = range(pair[0]-fuzziness, pair[0]+fuzziness+1)
    fuzz_dest = range(pair[1]-fuzziness, pair[1]+fuzziness+1)
    fuzz_scs = scd.loc[[p for p in fuzz_dest]][["sc_count", "pos"]]
    cands = sorted([(p, count_asc - c)
                    for p, c in zip(fuzz_scs["pos"], fuzz_scs["sc_count"])
                    if not math.isnan(c) and
                    {item["pos"]
                     for item in scd.loc[p]["asc"]}.intersection(fuzz_origin)],
                   key=lambda x: x[1])
    if not cands:
        return
    return pair[0], cands[0][0]


def calc_scs_pairs(scd, region, fuzziness):
    """Returns all reciprocating soft clip event within the given region.

    :param pd.DataFrame scd: Pandas data frame of soft clip event counts.
    :param range region: Range over which the reciprocating soft clip event
                         will be searched.
    :param int fuzziness: Length of 5' and 3' extension between which
                          each reciprocating soft clip event will be
                          searched.
    :returns: A dictionary whose keys are the reciprocating soft clip
              event and whose values are the soft clip ratios of the
              keys.

    """
    ascs = [item for item in scd["asc"]]
    poss = [v for v in scd.index.values]
    all_alt_pairs = {(pa[0], asc["pos"])
                     for pa in zip(poss, ascs)
                     for asc in pa[1]
                     if region.start <= pa[0] < region.stop and
                     region.start <= asc["pos"] < region.stop}
    sc_ratios = {p: r
                 for p, r in zip((v for _, v in scd["pos"].iteritems()),
                                 (v for _, v in scd["sc_ratio"].iteritems()))
                 if isinstance(r, (int, float)) and not math.isnan(r)}

    # adjust arc pairs since the match can be fuzzy
    adj_arc_pairs = {}
    for ap in all_alt_pairs:
        np = adjust_pair_fuzz(ap, scd, fuzziness)
        if np is not None:
            try:
                adj_arc_pairs[np] = (sc_ratios[np[0]], sc_ratios[np[1]])
            except KeyError:
                pass

    # only return reciprocal ones
    larger_starts = [p for p in adj_arc_pairs.keys() if p[0] > p[1]]
    return {k: v for k, v in adj_arc_pairs.items()
            if (k[1], k[0]) in larger_starts}


def plot_sample_df(countd, region, sc_fuzziness, sc_bg,
                   min_insert_count, output_fname=None):
    """Plots the given count data over the given region.

    :param dict countd: Dictionary of soft clip and insertion event counts.
    :param range region: Range over which the plot will be made.
    :param int sc_fuzziness: Length of 5' and 3' extension between which
                             each reciprocating soft clip event will be
                             searched.
    # :param int min_sc_count: Minimum number of count for a soft clip event
    #                          to be considered.
    :param int min_insert_count: Mininum number of for an insertion event
                                 to be considered.
    :param str output_fname: Name of output file of the plot.

    """
    sample = countd["sampleName"]
    scd = make_sc_sample_df(countd, region, sc_bg=sc_bg)
    ind = make_insert_sample_df(countd, region, min_count=min_insert_count)
    scs_pairs = calc_scs_pairs(scd, region, sc_fuzziness)

    plt.rcParams["figure.figsize"] = (16 + 0.1 * len(ind["insert_seq"]), 4)
    plt.suptitle("Sample {0!r}\nFLT3 Soft Clip & Insertion Plot"
                 "".format(sample), y=1.01)
    grid = gs.GridSpec(1, 2, width_ratios=[6, 1])
    gx1 = plt.subplot(grid[0])
    gx2 = plt.subplot(grid[1], sharey=gx1)
    grid.update(wspace=0.025)
    axes = [gx1, gx2]

    has_scs = True
    if not [x for x in scd["sc_ratio"] if not math.isnan(x)]:
        axes[0].set_xlabel("No Soft Clips Found")
        axes[0].set_xticklabels([])
        axes[0].xaxis.grid(False)
        has_scs = False
    else:
        ax1 = scd.plot.scatter(x="pos", y="sc_ratio",
                               marker="o", color="#4C72B0",
                               s=25,
                               zorder=2, ax=axes[0])

        yspan = max(0.15, scd["sc_ratio"].max(), ind["insert_ratio"].max())
        xmin = min(scd["pos"])
        xmax = max(scd["pos"])
        xspan = xmax - xmin
        ax1.set_ylim(yspan * -0.1, yspan * 1.4)
        ax1.set_xlim(xmin - (xspan * 0.05), xmax + (xspan * 0.05))

        ax1.spines["bottom"].set_color("#444444")
        ax1.spines["left"].set_color("#444444")
        ax1.tick_params(direction="inout", length=6, top=False, right=False)
        ax1.yaxis.set_major_formatter(FuncFormatter(to_percent))
        ax1.xaxis.set_major_formatter(
            FuncFormatter(lambda x, p: format(int(x), ",")))

        ax1.set_xlabel("Soft Clip Positions in Transcript")
        ax1.set_ylabel("% of Total Coverage")

        for ((p0, p1), (v0, v1)) in scs_pairs.items():
            codes = [pth.MOVETO, pth.CURVE4, pth.CURVE4, pth.CURVE4]
            verts = [(p0, v0), (p0, v0+(yspan*0.3)),
                     (p1, v1+(yspan*0.3)), (p1, v1)]
            paths = pth(verts, codes)
            patch = ptc.PathPatch(paths, facecolor="none", edgecolor="#E24A33",
                                  lw=2, zorder=1)
            ax1.add_patch(patch)

    if not len(ind["insert_seq"]):
        axes[1].set_xlabel("No Insertions Found")
        axes[1].set_xticklabels([])
        if not has_scs:
            axes[1].set_yticklabels([])
        axes[1].xaxis.grid(False)
    else:
        ax2 = ind.plot.bar(x="insert_seq", y="insert_ratio",
                           width=0.9, color="#55A868",
                           edgecolor="none",
                           legend=False, ax=axes[1])

        ax2.spines["bottom"].set_color("#444444")
        ax2.set_xlabel("Insertions")

        rects = ax2.patches
        annots = ["{0}bp".format(len(x)) for x in ind["insert_seq"]]
        ax2.set_xticklabels(["pos:{0}".format(x) for x in ind["pos"]],
                            rotation="vertical")
        for rect, annot in zip(rects, annots):
            height = rect.get_height()
            ax2.text(rect.get_x() + rect.get_width() / 2, height, annot,
                     size=9, ha="center", va="bottom")

    if output_fname is not None:
        plt.savefig(output_fname, bbox_inches="tight")

    return grid


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("input",
                type=click.Path(exists=True, dir_okay=False))
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("output", type=str)
@click.option("--fuzziness", type=int, default=12,
              help="Length of 5' and 3' extension when considering soft"
                   " clip matching.")
@click.option("--sc-bg", type=click.Path(exists=True, dir_okay=False))
@click.option("--min-insert-count", type=int, default=2,
              help="Minimum count of insertion event to consider.")
@click.option("--padding", type=int, default=10,
              help="Length of 5' and 3' extension of the given region.")
def main(input, start, end, output, fuzziness, sc_bg,
         min_insert_count, padding):

    ext = output.lower().rsplit(".", 1)
    if not ext or not ext[-1] in ("png", "svg"):
        raise click.BadParameter("Output file must have a PNG or an SVG"
                                 " extension.")

    if padding < 0:
        raise click.BadParameter("Length of padding must be at least 0.")

    with open(input, "r") as src:
        countd = json.load(src)

    if (start - 1) < countd["region"]["start"]:
        raise click.BadParameter("Invalid start coordinate: {0!r}."
                                 "".format(start))
    if end > countd["region"]["end"]:
        raise click.BadParameter("Invalid end coordinate: {0!r}."
                                 "".format(end))
    region = range(max(0, start - padding),
                   min(countd["region"]["end"], end + padding + 1))

    plt.style.use("seaborn-colorblind")
    if sc_bg is not None:
        bg = pd.DataFrame.from_csv(sc_bg, parse_dates=False)
    else:
        bg = None
    plot_sample_df(countd, region, fuzziness, bg, min_insert_count,
                   output_fname=output)


if __name__ == "__main__":
    main()
