import argparse
import json
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
            "page-offset":-1,
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

    def write(self, html, pdf) -> None:
        """Writes the report to the given path."""
        toc = {"xsl-style-sheet": self.toc_fname}
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

            if html:
                with open(html, "wt") as fout:
                    fout.write(cov_txt)
                    fout.write(con_txt)
            if pdf:
                import pdfkit
                pdfkit.from_string(con_txt, pdf,
                               options=self.pdfkit_opts, css=self.css_fname,
                               toc=toc, cover=cov_fh.name, cover_first=True)


def main(input_summary_path, css_path, templates_dir,
         imgs_dir, toc_path, html, pdf):
    """Script for generating PDF report of a sample analyzed with the Hamlet
    pipeline."""
    with open(input_summary_path) as src:
        sd = json.load(src)

    # Format p-values for readability
    for gene in sd["modules"]["snv_indels"]["genes"].values():
        for variant in gene:
            variant["FORMAT"]["PVAL"] = "{:0.1e}".format(float(variant["FORMAT"]["PVAL"]))

    sdm = sd["metadata"]
    sample_name = sdm["sample_name"]
    run_name = sdm["run_name"]
    pipeline_version = sdm["pipeline_version"]
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
    report.write(html, pdf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input_summary_path")
    parser.add_argument("--css-path", default="report/assets/style.css")
    parser.add_argument("--templates-dir", default="report/templates")
    parser.add_argument("--imgs-dir", default="report/assets/img")
    parser.add_argument("--toc-path", default="report/assets/toc.xsl")
    parser.add_argument("--html-output", dest="html")
    parser.add_argument("--pdf-output", dest="pdf")

    args = parser.parse_args()

    main(args.input_summary_path, args.css_path,
        args.templates_dir, args.imgs_dir, args.toc_path, args.html, args.pdf)
