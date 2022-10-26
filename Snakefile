include: "common.smk"


localrules:
    create_summary,
    generate_report,
    package_results,


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
    smallvars_plots="{sample}/snv-indels/variant_plots",
    # Stats
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
        [expand(p, sample=samples, pair={"R1", "R2"}) for p in OUTPUTS.values()],


# Define HAMLET modules
module qc_seq:
    snakefile:
        "includes/qc-seq/Snakefile"
    config:
        config


use rule * from qc_seq as qc_seq_*


module itd:
    snakefile:
        "includes/itd/Snakefile"
    config:
        config


use rule * from itd as itd_*


# Connect the align_kmt2a rule to the output of qc-seq
use rule align_kmt2a from itd as itd_align_kmt2a with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        fasta=config["kmt2a_fasta"],


# Connect the align_flt3 rule to the output of qc-seq
use rule align_flt3 from itd as itd_align_flt3 with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        fasta=config["flt3_fasta"],


module align:
    snakefile:
        "includes/snv-indels/Snakefile"
    config:
        config


use rule * from align as align_*


# Connect the align rule to the output of qc-seq
use rule align_vars from align as align_align_vars with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        index=config.get("genome_gmap_index") or "gmap_index/reference",


module expression:
    snakefile:
        "includes/expression/Snakefile"
    config:
        config


use rule * from expression as expression_*


# Connect the idsort rule to the output of snv-indels
use rule idsort_aln from expression as expression_idsort_aln with:
    input:
        bam=align.module_output.bam,
    container:
        "docker://quay.io/biocontainers/picard:2.20.5--0"


# Connect the count_raw_bases rule to the output of snv-indels
use rule count_raw_bases from expression as expression_count_raw_bases with:
    input:
        bam=align.module_output.bam,
        bed=config["expression_bed"],
        count_script=config["base_count_script"],
    container:
        "docker://quay.io/biocontainers/mulled-v2-a9ddcbd438a66450297b5e0b61ac390ee9bfdb61:e60f3cfda0dfcf4a72f2091c6fa1ebe5a5400220-0"


module fusion:
    snakefile:
        "includes/fusion/Snakefile"
    config:
        config


use rule * from fusion as fusion_*


# Connect the star_fusion rule to the output of qc-seq
use rule star_fusion from fusion as fusion_star_fusion with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        lib=config["genome_star_fusion_lib"],
    container:
        "docker://quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1"


# Connect the fusioncather rule to the output of qc-seq. Fusioncatcher should
# use the raw (untrimmed) output files, according to the manual.
use rule fusioncatcher from fusion as fusion_fusioncatcher with:
    input:
        fq1=qc_seq.module_output.forward_raw,
        fq2=qc_seq.module_output.reverse_raw,
    container:
        "docker://quay.io/biocontainers/fusioncatcher:1.20--2"


# Fusioncatcher outputs
if config.get("fusioncatcher_data"):
    OUTPUTS["fusioncatcher_txt"] = fusion.module_output.optional.fusion_catcher
    OUTPUTS["fusions_txt"] = fusion.module_output.optional.intersect
    OUTPUTS["isect_txt"] = fusion.module_output.optional.subset_predictions


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        aln_stats=OUTPUTS["aln_stats"],
        rna_stats=OUTPUTS["rna_stats"],
        insert_stats=OUTPUTS["insert_stats"],
        vep_stats=OUTPUTS["vep_stats"],
        exon_cov_stats=OUTPUTS["exon_cov_stats"],
        idm=config["ref_id_mapping"],
        var_plots=OUTPUTS["smallvars_plots"],
        var_csv=OUTPUTS["smallvars_csv_hi"],
        flt3_plot=OUTPUTS["flt3_png"],
        kmt2a_plot=OUTPUTS["kmt2a_png"],
        flt3_csv=OUTPUTS["flt3_csv"],
        kmt2a_csv=OUTPUTS["kmt2a_csv"],
        exon_ratios=expression.module_output.json,
        qc_seq_json=qc_seq.module_output.json,
        fusion_json=fusion.module_output.json,
        scr=srcdir("scripts/create_summary.py"),
    params:
        pipeline_ver=PIPELINE_VERSION,
        run_name=RUN_NAME,
    output:
        js=OUTPUTS["summary"],
    log:
        "log/create_summary.{sample}.txt",
    container:
        containers["crimson"]
    shell:
        """
        python {input.scr} \
            {input.idm} \
            `dirname {input.var_plots}` \
            {input.var_csv} \
            {input.flt3_csv} {input.flt3_plot} \
            {input.kmt2a_csv} {input.kmt2a_plot} \
            {input.aln_stats} {input.rna_stats} {input.insert_stats} \
            {input.exon_cov_stats} {input.vep_stats} \
            --pipeline-version {params.pipeline_ver} \
            --run-name {params.run_name} \
            --sample-name {wildcards.sample} \
            --module {input.exon_ratios} \
            --module {input.fusion_json} \
            --module {input.qc_seq_json} > {output.js} 2>{log}
        """


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
    log:
        "log/generate_report.{sample}.txt",
    container:
        containers["hamlet-scripts"]
    shell:
        """
        python3 {input.scr} \
            --templates-dir {input.templates} \
            --imgs-dir {input.imgs} \
            --css-path {input.css} \
            --toc-path {input.toc} \
            {input.summary} \
            {output.pdf} 2> {log}
        """


rule package_results:
    """Copies essential result files into one directory and zips it."""
    input:
        summary=OUTPUTS["summary"],
        smallvars_csv_all=OUTPUTS["smallvars_csv_all"],
        smallvars_csv_hi=OUTPUTS["smallvars_csv_hi"],
        smallvars_plots=OUTPUTS["smallvars_plots"],
        count_fragments_per_gene=expression.module_output.fragments_per_gene,
        count_bases_per_gene=expression.module_output.bases_per_gene,
        count_bases_per_exon=expression.module_output.bases_per_exon,
        ratio_exons=expression.module_output.json,
        flt_csv=OUTPUTS["flt3_csv"],
        flt_bg_csv=OUTPUTS["flt3_bg_csv"],
        flt_png=OUTPUTS["flt3_png"],
        kmt_csv=OUTPUTS["kmt2a_csv"],
        kmt_bg_csv=OUTPUTS["kmt2a_bg_csv"],
        kmt_png=OUTPUTS["kmt2a_png"],
        reportje=OUTPUTS["reportje"],
    params:
        tmp=f"tmp/hamlet-pkg.{{sample}}.{uuid4()}/hamlet_results.{{sample}}",
    output:
        pkg=OUTPUTS["package"],
    log:
        "log/package_results.{sample}.txt",
    container:
        containers["zip"]
    shell:
        """
        mkdir -p {params.tmp} && \
        cp -r {input} {params.tmp} && \
        cp -r $(dirname {input.smallvars_plots}) {params.tmp} && \
        zip -9 -x *.done -r {output.pkg} {params.tmp} && \
        rm -rf {params.tmp} 2> {log} \
        \
        || rm -rf {params.tmp}
        """
