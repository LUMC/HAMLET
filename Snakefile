from rattle import Run

RUN = Run(config)

include: "includes/qc-seq/Snakefile"
include: "includes/snv-indels/Snakefile"
include: "includes/fusion/Snakefile"
include: "includes/expression/Snakefile"
include: "includes/itd/Snakefile"

OUTPUTS = dict(
    # Merged FASTQs
    fqs="{sample}/{sample}-{pair}.fq.gz",

    # Small variants
    smallvars_bam="{sample}/snv-indels/{sample}.snv-indel.bam",
    smallvars_vcf="{sample}/snv-indels/{sample}.annotated.vcf.gz",
    smallvars_csv_all="{sample}/snv-indels/{sample}.variants_all.csv",
    smallvars_csv_hi="{sample}/snv-indels/{sample}.variants_hi.csv",
    smallvars_plots="{sample}/snv-indels/variant_plots/.done",

    # Fusion
    star_fusion_txt="{sample}/fusion/{sample}.star-fusion",
    fusioncatcher_txt="{sample}/fusion/{sample}.fusioncatcher",
    star_fusion_svg="{sample}/fusion/{sample}.star-fusion.svg",
    fusioncatcher_svg="{sample}/fusion/{sample}.fusioncatcher.svg",
    fusions_txt="{sample}/fusion/{sample}.fuma",
    isect_svg="{sample}/fusion/{sample}.sf-isect.svg",
    isect_txt="{sample}/fusion/{sample}.sf-isect",
    fusions_svg="{sample}/fusion/{sample}.fusions-combined.svg",

    # Expression
    count_fragments_per_gene="{sample}/expression/{sample}.fragments_per_gene",
    count_bases_per_gene="{sample}/expression/{sample}.bases_per_gene",
    count_bases_per_exon="{sample}/expression/{sample}.bases_per_exon",
    ratio_exons="{sample}/expression/{sample}.exon_ratios",

    # Stats
    fqs_stats="{sample}/qc-seq/{sample}-seq-stats.json",
    rna_stats="{sample}/snv-indels/{sample}.rna_stats",
    aln_stats="{sample}/snv-indels/{sample}.aln_stats",
    insert_stats="{sample}/snv-indels/{sample}.insert_stats",

    # ITD module
    flt3_bam="{sample}/itd/{sample}.flt3.bam",
    flt3_csv="{sample}/itd/{sample}.flt3.csv",
    flt3_bg_csv="{sample}/itd/{sample}.flt3.bg.csv",
    flt3_png="{sample}/itd/{sample}.flt3.png",
    kmt2a_bam="{sample}/itd/{sample}.kmt2a.bam",
    kmt2a_csv="{sample}/itd/{sample}.kmt2a.csv",
    kmt2a_bg_csv="{sample}/itd/{sample}.kmt2a.bg.csv",
    kmt2a_png="{sample}/itd/{sample}.kmt2a.png",
)


rule all:
    input:
        [expand(RUN.output(p), sample=RUN.samples, pair={"R1", "R2"})
         for p in OUTPUTS.values()]
