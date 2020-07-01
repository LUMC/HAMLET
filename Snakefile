import os
import subprocess
from functools import partial
from os.path import dirname
from uuid import uuid4

include: "includes/qc-seq/Snakefile"
include: "includes/snv-indels/Snakefile"
include: "includes/fusion/Snakefile"
include: "includes/expression/Snakefile"
include: "includes/itd/Snakefile"

localrules: create_summary, generate_report, package_results

containers = {
    "fsnviz": "docker://quay.io/biocontainers/fsnviz:0.3.0--py_3",
    "hamlet-scripts": "docker://lumc/hamlet-scripts:0.3",
    "debian": "docker://debian:buster-slim",
    "zip": "docker://lumc/zip:3.0"
}

settings=config["settings"]

# Determine version of HAMLET from git
out = subprocess.Popen(['git', 'describe', '--tags'], stdout=subprocess.PIPE)
stdout, stderr = out.communicate()
PIPELINE_VERSION = stdout.strip().decode('utf-8')

RUN_NAME = settings.get("run_name") or f"hamlet-{uuid4().hex[:8]}"


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
    reportje="{sample}/hamlet_report.{sample}.pdf",
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

    # Fusioncatcher
    fusioncatcher_txt=fusion_output(".fusioncatcher"),
    fusioncatcher_svg=fusion_output(".fusioncatcher.svg"),
    fusions_txt=fusion_output(".fuma"),
    isect_svg=fusion_output(".sf-isect.svg"),
    isect_txt=fusion_output(".sf-isect"),

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
    vep_stats=var_output(".vep_stats.txt"),
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

rule all:
    input:
        [expand(p, sample=config["samples"], pair={"R1", "R2"})
         for p in OUTPUTS.values()]


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        seq_stats=OUTPUTS["seq_stats"],
        aln_stats=OUTPUTS["aln_stats"],
        rna_stats=OUTPUTS["rna_stats"],
        insert_stats=OUTPUTS["insert_stats"],
        vep_stats=OUTPUTS["vep_stats"],
        exon_cov_stats=OUTPUTS["exon_cov_stats"],
        idm=settings["ref_id_mapping"],
        var_plots=OUTPUTS["smallvars_plots"],
        var_csv=OUTPUTS["smallvars_csv_hi"],
        fusions_svg=OUTPUTS["fusions_svg"],
        flt3_plot=OUTPUTS["flt3_png"],
        kmt2a_plot=OUTPUTS["kmt2a_png"],
        flt3_csv=OUTPUTS["flt3_csv"],
        kmt2a_csv=OUTPUTS["kmt2a_csv"],
        exon_ratios=OUTPUTS["ratio_exons"],
        scr=srcdir("scripts/create_summary.py"),
    params:
        pipeline_ver=PIPELINE_VERSION,
        run_name=RUN_NAME,
    output:
        js=OUTPUTS["summary"]
    singularity: containers["fsnviz"]
    shell:
        "python {input.scr}"
        " {input.idm}"
        " `dirname {input.var_plots}` {input.var_csv}"
        " `dirname {input.fusions_svg}`"
        " {input.flt3_csv} {input.flt3_plot}"
        " {input.kmt2a_csv} {input.kmt2a_plot}"
        " {input.exon_ratios}"
        " {input.seq_stats} {input.aln_stats} {input.rna_stats} {input.insert_stats}"
        " {input.exon_cov_stats} {input.vep_stats}"
        " --pipeline-version {params.pipeline_ver}"
        " --run-name {params.run_name}"
        " --sample-name {wildcards.sample}"
        " > {output.js}"


rule generate_report:
    """Generates a PDF report of the essential results."""
    input:
        summary=OUTPUTS["summary"],
        css=srcdir("report/assets/style.css"),
        templates=srcdir("report/templates"),
        imgs=srcdir("report/assets/img"),
        toc=srcdir("report/assets/toc.xsl"),
        scr=srcdir("scripts/generate_report.py"),
    output:
        pdf=OUTPUTS["reportje"],
    singularity: containers["hamlet-scripts"]
    shell:
        "python3 {input.scr}"
        " --templates-dir {input.templates} --imgs-dir {input.imgs}"
        " --css-path {input.css} --toc-path {input.toc}"
        " {input.summary} {output.pdf}"


rule package_results:
    """Copies essential result files into one directory and zips it."""
    input:
        summary=OUTPUTS["summary"],
        smallvars_csv_all=OUTPUTS["smallvars_csv_all"],
        smallvars_csv_hi=OUTPUTS["smallvars_csv_hi"],
        smallvars_plots=OUTPUTS["smallvars_plots"],
        fusions_svg=OUTPUTS["fusions_svg"],
        count_fragments_per_gene=OUTPUTS["count_fragments_per_gene"],
        count_bases_per_gene=OUTPUTS["count_bases_per_gene"],
        count_bases_per_exon=OUTPUTS["count_bases_per_exon"],
        ratio_exons=OUTPUTS["ratio_exons"],
        flt_csv=OUTPUTS["flt3_csv"],
        flt_bg_csv=OUTPUTS["flt3_bg_csv"],
        flt_png=OUTPUTS["flt3_png"],
        kmt_csv=OUTPUTS["kmt2a_csv"],
        kmt_bg_csv=OUTPUTS["kmt2a_bg_csv"],
        kmt_png=OUTPUTS["kmt2a_png"],
        reportje=OUTPUTS["reportje"],
    output:
        pkg=OUTPUTS["package"],
    params:
        tmp="tmp/hamlet-pkg.{sample}." + str(uuid4()) + "/hamlet_results.{sample}"
    singularity: containers["zip"]
    shell:
        "(mkdir -p {params.tmp}"
        " && cp -r {input} {params.tmp}"
        " && cp -r $(dirname {input.smallvars_plots}) {params.tmp}"
        " && zip -9 -x *.done -r {output.pkg} {params.tmp}"
        " && rm -rf {params.tmp})"
        " || rm -rf {params.tmp}"
