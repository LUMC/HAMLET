import argparse
import json
import os
import shutil
from datetime import datetime as dt
from pathlib import Path
from tempfile import NamedTemporaryFile as NTF
from typing import Optional

from jinja2 import Environment, FileSystemLoader


class Report(object):

    """Analysis report of a single sample."""

    def __init__(self, summaryd: dict,
                 tpl_dir: str="templates",
                 imgs_dir: str="assets/img",
                 cover_tpl_fname: str="cover.html.j2",
                 contents_tpl_fname: str="contents.html.j2",
                 css_fname: str="assets/style.css",
                 toc_fname: str="assets/toc.xsl",
                 header_line: bool=True,
                 header_caption: Optional[str]=None,
                 footer_line: bool=True,
                 footer_lcaption: Optional[str]=None,
                 footer_rcaption: Optional[str]=None,
                 pdfkit_opts: Optional[dict]=None,
                 timestamp: Optional[dt]=None) -> None:
        sdm = summaryd["metadata"]
        self.summary = summaryd
        self.sample_name = sdm["sample_name"]
        self.run_name = sdm["run_name"]
        self.pipeline_version = sdm["pipeline_version"]
        self.cover_tpl_fname = cover_tpl_fname
        self.contents_tpl_fname = contents_tpl_fname
        self.css_fname = css_fname
        self.toc_fname = toc_fname
        self.imgs_dir = str(Path(imgs_dir).resolve())

        self.timestamp = timestamp or dt.now()

        hf_ctx = {
            "sample_name": self.sample_name,
            "run_name": self.run_name,
            "timestamp": self.timestamp,
        }
        pdfkit_opts = pdfkit_opts or {
            "title": header_caption,
            "encoding": "UTF-8",
            "quiet": None,

            "page-size": "A4",
            "page-offset": -1,
            "margin-top": "20",
            "margin-right": "16",
            "margin-bottom": "20",
            "margin-left": "16",

            "header-spacing": "5",
            "footer-spacing": "4",
        }
        if header_line:
            pdfkit_opts["header-line"] = ""
        if header_caption is not None:
            pdfkit_opts["header-center"] = header_caption.format(**hf_ctx)
        if footer_line:
            pdfkit_opts["footer-line"] = ""
        if footer_lcaption is not None:
            pdfkit_opts["footer-left"] = footer_lcaption.format(**hf_ctx)
        if footer_rcaption is not None:
            pdfkit_opts["footer-right"] = footer_rcaption.format(**hf_ctx)
        self.pdfkit_opts = pdfkit_opts

        def show_int(value):
            if value is None or value == "":
                return "?"
            return "{:,d}".format(int(value))

        def show_pct(value1, value2):
            if value2 == 0:
                return "undefined"
            elif any(v is None or v == "" for v in (value1, value2)):
                return "?"
            return "{:,.2f}%".format(value1 * 100. / value2)

        def show_float(value, spec=".3g"):
            if value is None or value == "":
                return "?"
            fmt = "{:," + spec + "}"
            return fmt.format(float(value))

        def as_pct(value):
            if value is None or value == "":
                return "?"
            return "{:,.2f}%".format(float(value) * 100.)

        def num_tids(idm):
            return sum([len(v["transcript_ids"]) for v in idm])

        env = Environment(loader=FileSystemLoader(tpl_dir))
        env.filters["show_int"] = show_int
        env.filters["show_pct"] = show_pct
        env.filters["show_float"] = show_float
        env.filters["as_pct"] = as_pct
        env.filters["num_tids"] = num_tids
        self.env = env
        self.cover_tpl = env.get_template(cover_tpl_fname)
        self.contents_tpl = env.get_template(contents_tpl_fname)

    def write(self, output_path: Path) -> None:
        """Writes the report to the given path."""
        tmp_prefix = str(Path.cwd()) + "/"
        cover_ctx = {
            "sample_name": self.sample_name,
            "run_name": self.run_name,
            "pipeline_version": self.pipeline_version,
            "timestamp": self.timestamp,
            "css_fname": self.css_fname,
            "imgs_dir": self.imgs_dir,
        }
        contents_ctx = self.summary
        contents_ctx["css_fname"] = self.css_fname
        contents_ctx["imgs_dir"] = self.imgs_dir

        with NTF(prefix=tmp_prefix, suffix=".html") as cov_fh:
            cov_txt = self.cover_tpl.render(**cover_ctx)
            cov_fh.write(cov_txt.encode("utf-8"))
            cov_fh.seek(0)

            con_txt = self.contents_tpl.render(**contents_ctx)
            with open(output_path, 'wt') as fout:
                print(con_txt, file=fout)


def localise_assets(sd, html_output, css_path):
    """ Localise the html assets to a folder next to the html output """

    # Create the report assets folder
    assets_folder = '.'.join(html_output.split('.')[:-1])
    os.makedirs(assets_folder, exist_ok=True)

    # Localise the fusion results, they all have the same name so we have to
    # include the containing folder in the file name
    for plot in sd['results']['fusion']['plots']:
        path = sd['results']['fusion']['plots'][plot]
        newpath = os.path.join(assets_folder, '_'.join(path.split('/')[-2:]))
        shutil.copy(path, newpath)
        local_html_path = '/'.join(newpath.split('/')[1:])
        sd['results']['fusion']['plots'][plot] = local_html_path

    # Localise itd
    for gene in sd['results']['itd']:
        path = sd['results']['itd'][gene]['path']
        new_path = os.path.join(assets_folder, os.path.basename(path))
        shutil.copy(path, new_path)
        local_html_path = '/'.join(new_path.split('/')[1:])
        sd['results']['itd'][gene]['path'] = local_html_path

    # Localise variant results
    for index, gene in enumerate(sd['results']['var']['plots']):
        path = gene['path']
        new_path = os.path.join(assets_folder, os.path.basename(path))
        shutil.copy(path, new_path)
        local_html_path = '/'.join(new_path.split('/')[1:])
        sd['results']['var']['plots'][index]['path'] = local_html_path

    # Copy the css file
    local_css = os.path.join(assets_folder, os.path.basename(css_path))
    shutil.copy(css_path, local_css)
    local_css_path = '/'.join(local_css.split('/')[1:])
    return local_css_path


def main(input_summary_path, css_path, templates_dir, imgs_dir, toc_path,
         html_output):
    """Script for generating PDF report of a sample analyzed with the Hamlet
    pipeline."""
    with open(input_summary_path) as src:
        sd = json.load(src)

    css_path = localise_assets(sd, html_output, css_path)

    sdm = sd["metadata"]
    sample_name = sdm["sample_name"]
    run_name = sdm["run_name"]
    header_caption = ("Hamlet Report -"
                      f" Sample {sample_name!r} of Run {run_name!r}")
    footer_lcaption = "Generated on {timestamp:%A, %d %B %Y at %H:%M}"
    footer_rcaption = "[page]/[toPage]"

    report = Report(sd,
                    tpl_dir=templates_dir,
                    imgs_dir=imgs_dir,
                    css_fname=css_path,
                    toc_fname=toc_path,
                    header_caption=header_caption,
                    footer_lcaption=footer_lcaption,
                    footer_rcaption=footer_rcaption)
    report.write(Path(html_output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-summary", required=True,
                        help="Input summary path")
    parser.add_argument("--css-path", default="assets/style.css",
                        help="Path to css file")
    parser.add_argument("--templates-dir", default="templates",
                        help="Path to templates directory")
    parser.add_argument("--imgs-dir", default="assets/img",
                        help="Path to images directory")
    parser.add_argument("--toc-path", default="assets/toc.xsl",
                        help="Path to table of content")
    parser.add_argument("--html-output", required=True,
                        help="Path to output HTML file")

    args = parser.parse_args()
    main(args.input_summary, args.css_path, args.templates_dir, args.imgs_dir,
         args.toc_path, args.html_output)
