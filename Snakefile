import os
from functools import partial
from os.path import dirname
from uuid import uuid4

import git
from rattle import Run


BASE_PIPELINE_VERSION = "0.1.0"
try:
    repo = git.Repo(path=srcdir(""), search_parent_directories=True)
except git.exc.InvalidGitRepositoryError:
    repo = None
    sha = "unknown"
    is_dirty = "?"
else:
    sha = repo.head.object.hexsha[:8]
    is_dirty = "*" if repo.is_dirty() else ""

PIPELINE_VERSION = f"{BASE_PIPELINE_VERSION}-{sha}{is_dirty}"


RUN = Run(config)

RUN_NAME = RUN.settings.get("run_name") or f"hamlet-{uuid4().hex[:8]}"


include: "includes/qc-seq/Snakefile"
include: "includes/snv-indels/Snakefile"
include: "includes/fusion/Snakefile"
include: "includes/expression/Snakefile"
include: "includes/itd/Snakefile"


def make_pattern(extension, dirname):
    """Helper function to create a wildcard-containing path for output files."""
    return f"{{sample}}/{dirname}/{{sample}}{extension}"


seqqc_output = partial(make_pattern, dirname="qc-seq")
var_output = partial(make_pattern, dirname="snv-indels")
fusion_output = partial(make_pattern, dirname="fusion")
expr_output = partial(make_pattern, dirname="expression")
itd_output = partial(make_pattern, dirname="itd")


OUTPUTS = dict(
    # Merged FASTQs, stats, and packaged results
    fqs="{sample}/{sample}-{pair}.fq.gz",
    summary="{sample}/{sample}.summary.json",
    report="{sample}/hamlet_report.{sample}.pdf",
    package="{sample}/hamlet_results.{sample}.zip",

    # Small variants
    smallvars_bam=var_output(".snv-indel.bam"),
    smallvars_vcf=var_output(".annotated.vcf.gz"),
    smallvars_csv_all=var_output(".variants_all.csv"),
    smallvars_csv_hi=var_output(".variants_hi.csv"),
    smallvars_plots="{sample}/snv-indels/variant_plots/.done",

    # Fusion
    star_fusion_txt=fusion_output(".star-fusion"),
    star_fusion_svg=fusion_output(".star-fusion.svg"),
    fusions_svg=fusion_output(".fusions-combined.svg"),

    # Expression
    count_fragments_per_gene=expr_output(".fragments_per_gene"),
    count_bases_per_gene=expr_output(".bases_per_gene"),
    count_bases_per_exon=expr_output(".bases_per_exon"),
    ratio_exons=expr_output(".exon_ratios"),

    # Stats
    seq_stats=seqqc_output(".seq_stats.json"),
    aln_stats=var_output(".aln_stats"),
    rna_stats=var_output(".rna_stats"),
    insert_stats=var_output(".insert_stats"),
    vep_stats=var_output(".vep_stats"),
    exon_cov_stats=var_output(".exon_cov_stats.json"),

    # ITD module
    flt3_bam=itd_output(".flt3.bam"),
    flt3_csv=itd_output(".flt3.csv"),
    flt3_bg_csv=itd_output(".flt3.bg.csv"),
    flt3_png=itd_output(".flt3.png"),
    kmt2a_bam=itd_output(".kmt2a.bam"),
    kmt2a_csv=itd_output(".kmt2a.csv"),
    kmt2a_bg_csv=itd_output(".kmt2a.bg.csv"),
    kmt2a_png=itd_output(".kmt2a.png"),
)

if "fusioncatcher_exe" in RUN.settings:
    OUTPUTS.update(
        dict(
            fusioncatcher_txt=fusion_output(".fusioncatcher"),
            fusioncatcher_svg=fusion_output(".fusioncatcher.svg"),
            fusions_txt=fusion_output(".fuma"),
            isect_svg=fusion_output(".sf-isect.svg"),
            isect_txt=fusion_output(".sf-isect"),
    ))


rule all:
    input:
        [expand(RUN.output(p), sample=RUN.samples, pair={"R1", "R2"})
         for p in OUTPUTS.values()]


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        seq_stats=RUN.output(OUTPUTS["seq_stats"]),
        aln_stats=RUN.output(OUTPUTS["aln_stats"]),
        rna_stats=RUN.output(OUTPUTS["rna_stats"]),
        insert_stats=RUN.output(OUTPUTS["insert_stats"]),
        vep_stats=RUN.output(OUTPUTS["vep_stats"]),
        exon_cov_stats=RUN.output(OUTPUTS["exon_cov_stats"]),
        idm=RUN.settings["ref_id_mapping"],
        var_plots=RUN.output(OUTPUTS["smallvars_plots"]),
        fusions_svg=RUN.output(OUTPUTS["fusions_svg"]),
        flt3_plot=RUN.output(OUTPUTS["flt3_png"]),
        kmt2a_plot=RUN.output(OUTPUTS["kmt2a_png"]),
        exon_ratios=RUN.output(OUTPUTS["ratio_exons"]),
        scr=srcdir("scripts/create_summary.py"),
    params:
        pipeline_ver=PIPELINE_VERSION,
        run_name=RUN_NAME,
    output:
        js=RUN.output(OUTPUTS["summary"])
    conda: srcdir("envs/create_summary.yml")
    shell:
        "python {input.scr}"
        " {input.idm}"
        " `dirname {input.var_plots}`"
        " `dirname {input.fusions_svg}`"
        " {input.flt3_plot} {input.kmt2a_plot} {input.exon_ratios}"
        " {input.seq_stats} {input.aln_stats} {input.rna_stats} {input.insert_stats}"
        " {input.exon_cov_stats} {input.vep_stats}"
        " --pipeline-version {params.pipeline_ver}"
        " --run-name {params.run_name}"
        " --sample-name {wildcards.sample}"
        " > {output.js}"


rule generate_report:
    """Generates a PDF report of the essential results."""
    input:
        summary=RUN.output(OUTPUTS["summary"]),
        css=srcdir("report/assets/style.css"),
        templates=srcdir("report/templates"),
        imgs=srcdir("report/assets/img"),
        toc=srcdir("report/assets/toc.xsl"),
        scr=srcdir("scripts/generate_report.py"),
    output:
        pdf=RUN.output(OUTPUTS["report"]),
    conda: srcdir("envs/create_report.yml")
    shell:
        "python {input.scr}"
        " --templates-dir {input.templates} --imgs-dir {input.imgs}"
        " --css-path {input.css} --toc-path {input.toc}"
        " {input.summary} {output.pdf}"


rule package_results:
    """Copies essential result files into one directory and zips it."""
    input:
        summary=RUN.output(OUTPUTS["summary"]),
        smallvars_csv_all=RUN.output(OUTPUTS["smallvars_csv_all"]),
        smallvars_csv_hi=RUN.output(OUTPUTS["smallvars_csv_hi"]),
        smallvars_plots_dir=RUN.output(dirname(OUTPUTS["smallvars_plots"])),
        fusions_svg=RUN.output(OUTPUTS["fusions_svg"]),
        count_fragments_per_gene=RUN.output(OUTPUTS["count_fragments_per_gene"]),
        count_bases_per_gene=RUN.output(OUTPUTS["count_bases_per_gene"]),
        count_bases_per_exon=RUN.output(OUTPUTS["count_bases_per_exon"]),
        ratio_exons=RUN.output(OUTPUTS["ratio_exons"]),
        flt_csv=RUN.output(OUTPUTS["flt3_csv"]),
        flt_bg_csv=RUN.output(OUTPUTS["flt3_bg_csv"]),
        flt_png=RUN.output(OUTPUTS["flt3_png"]),
        kmt_csv=RUN.output(OUTPUTS["kmt2a_csv"]),
        kmt_bg_csv=RUN.output(OUTPUTS["kmt2a_bg_csv"]),
        kmt_png=RUN.output(OUTPUTS["kmt2a_png"]),
        report=RUN.output(OUTPUTS["report"]),
    output:
        pkg=RUN.output(OUTPUTS["package"]),
    params:
        tmp="/tmp/hamlet-pkg.{sample}." + str(uuid4()) + "/hamlet_results.{sample}"
    conda: srcdir("envs/package_results.yml")
    shell:
        "(mkdir -p {params.tmp}"
        " && cp -r {input} {params.tmp}"
        " && cd $(dirname {params.tmp})"
        " && zip -9 -x *.done -r {output.pkg} *"
        " && cd -"
        " && rm -rf {params.tmp})"
        " || rm -rf {params.tmp}"
