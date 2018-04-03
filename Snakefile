from functools import partial

from rattle import Run


RUN = Run(config)


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
    # Merged FASTQs
    fqs="{sample}/{sample}-{pair}.fq.gz",

    # Small variants
    smallvars_bam=var_output(".snv-indel.bam"),
    smallvars_vcf=var_output(".annotated.vcf.gz"),
    smallvars_csv_all=var_output(".variants_all.csv"),
    smallvars_csv_hi=var_output(".variants_hi.csv"),
    smallvars_plots="{sample}/snv-indels/variant_plots/.done",

    # Fusion
    star_fusion_txt=fusion_output(".star-fusion"),
    fusioncatcher_txt=fusion_output(".fusioncatcher"),
    star_fusion_svg=fusion_output(".star-fusion.svg"),
    fusioncatcher_svg=fusion_output(".fusioncatcher.svg"),
    fusions_txt=fusion_output(".fuma"),
    isect_svg=fusion_output(".sf-isect.svg"),
    isect_txt=fusion_output(".sf-isect"),
    fusions_svg=fusion_output(".fusions-combined.svg"),

    # Expression
    count_fragments_per_gene=expr_output(".fragments_per_gene"),
    count_bases_per_gene=expr_output(".bases_per_gene"),
    count_bases_per_exon=expr_output(".bases_per_exon"),
    ratio_exons=expr_output(".exon_ratios"),

    # Stats
    fqs_stats=seqqc_output("-seq-stats.json"),
    rna_stats=var_output(".rna_stats"),
    aln_stats=var_output(".aln_stats"),
    insert_stats=var_output(".insert_stats"),

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
        [expand(RUN.output(p), sample=RUN.samples, pair={"R1", "R2"})
         for p in OUTPUTS.values()]
