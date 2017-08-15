#!/usr/bin/env python

import click
import svgutils.transform as sg


def load_single_fig(fname, width="4300px", height="4300px"):
    """Loads the SVG of a single figure and adjusts its size."""
    fig = sg.fromfile(fname)
    fig.width = width
    fig.height = height

    return fig


def combine_figs(*labeled_figs, unit="px"):
    """Combines the given figures horizontally."""
    # Assumes all figures have the same height and use the same size units
    plots = []
    labels = []
    cur_offset = 0
    for fig, label in labeled_figs:
        plot = fig.getroot()
        plot.moveto(cur_offset, 0)
        plots.append(plot)

        width = int(fig.width.strip(unit))
        left_label_offset = 0.0125 * width
        top_label_offset = 0.0375 * width
        font_size = 0.03 * width

        text = sg.TextElement(left_label_offset + cur_offset,
                              top_label_offset, label, size=font_size)
        labels.append(text)

        cur_offset += width

    height = labeled_figs[0][0].height

    combined_fig = sg.SVGFigure(cur_offset, height)
    combined_fig.append(plots)
    combined_fig.append(labels)

    return combined_fig


def combine_files(*inputs, output):
    """Combines the given input figures into a single output figure."""
    fig = combine_figs(*[(load_single_fig(fn), label) for label, fn in inputs],
                       unit="px")
    fig.save(output)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("inputs", nargs=-1, type=str)
@click.argument("output", nargs=1,
                type=click.Path(dir_okay=False))
def main(inputs, output):
    labeled_inputs = [item.rsplit(":", 1) for item in inputs]
    combine_files(*labeled_inputs, output=output)


if __name__ == "__main__":
    main()
