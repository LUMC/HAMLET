include: "common.smk"


localrules:
    align_exon_cov_ref,
    align_genome_txt,
    align_json_output,
    cleanup_multiqc,
    create_summary,
    generate_report,
    itd_detect_itd_flt3,
    itd_detect_itd_kmt2a,
    itd_json_output,
    itd_plot_itd_flt3,
    itd_plot_itd_kmt2a,


rule all:
    input:
        summary=expand("{sample}/{sample}.summary.json", sample=samples),
        report=expand("{sample}/hamlet_report.{sample}.pdf", sample=samples),
        multiqc="multiqc_hamlet.html",
        cleanup=".cleanup_multiqc",


# Add the PEP configuration to each submodule
config["qc-seq"]["pepfile"] = config["pepfile"]
config["snv-indels"]["pepfile"] = config["pepfile"]
config["itd"]["pepfile"] = config["pepfile"]
config["fusion"]["pepfile"] = config["pepfile"]
config["expression"]["pepfile"] = config["pepfile"]


# Define HAMLET modules
module qc_seq:
    snakefile:
        "includes/qc-seq/Snakefile"
    config:
        config["qc-seq"]


use rule * from qc_seq as qc_seq_*


# Make the final trimmed, merged FastQ files from the qc_seq module temporary
use rule cutadapt from qc_seq as qc_seq_cutadapt with:
    output:
        fq1=temp("{sample}/qc-seq/{sample}-R1.fq.gz"),
        fq2=temp("{sample}/qc-seq/{sample}-R2.fq.gz"),
        json="{sample}/qc-seq/{sample}.cutadapt.json",


# Mark qc-seq MultiQC files as temporary
use rule multiqc from qc_seq as qc_seq_multiqc with:
    output:
        html=temporary("multiqc_qc_seq.html"),
        folder=directory("multiqc_qc_seq_data"),
        parquet="multiqc_qc_seq_data/multiqc.parquet",
        filelist=temporary("multiqc_filelist_qc_seq.txt"),


module itd:
    snakefile:
        "includes/itd/Snakefile"
    config:
        config["itd"]


use rule * from itd as itd_*


# Connect the align_reads rule to the output of qc-seq
use rule align_reads from itd as itd_align_reads with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        fasta=config["itd"]["fasta"],


module align:
    snakefile:
        "includes/snv-indels/Snakefile"
    config:
        config["snv-indels"]


use rule * from align as align_*


# Connect the align rule to the output of qc-seq
use rule align_vars from align as align_align_vars with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        index=config["snv-indels"]["star_index"],
        gtf=config["snv-indels"]["gtf"],


# Mark align MultiQC files as temporary
use rule multiqc from align as align_multiqc with:
    output:
        html=temporary("multiqc_snv_indels.html"),
        folder=directory("multiqc_snv_indels_data"),
        parquet="multiqc_snv_indels_data/multiqc.parquet",
        filelist=temporary("multiqc_filelist_snv_indels.txt"),


module fusion:
    snakefile:
        "includes/fusion/Snakefile"
    config:
        config["fusion"]


use rule * from fusion as fusion_*


# Connect the output of snv-indels to the arriba rule
use rule arriba from fusion as fusion_arriba with:
    input:
        bam=align.module_output.bam,
        ref=config["fusion"]["genome_fasta"],
        gtf=config["fusion"]["gtf"],
        blacklist=config["fusion"]["blacklist"],
        known_fusions=config["fusion"]["known_fusions"],
        protein_domains=config["fusion"]["protein_domains"],
    container:
        fusion.containers["arriba"]


use rule plot_fusions from fusion as fusion_plot_fusions with:
    input:
        fusions="{sample}/fusion/arriba/fusions.tsv",
        bam=align.module_output.bam,
        bai=align.module_output.bai,
        gtf=config["fusion"]["gtf"],
        blacklist=config["fusion"]["blacklist"],
        cytobands=config["fusion"]["cytobands"],
        protein_domains=config["fusion"]["protein_domains"],
    container:
        fusion.containers["arriba"]


module expression:
    snakefile:
        "includes/expression/Snakefile"
    config:
        config["expression"]


use rule * from expression as expression_*


# Connect the output of snv-indels to expression
use rule normalized_coverage from expression as expression_normalized_coverage with:
    input:
        bam=align.module_output.bam,
        bai=align.module_output.bai,
        counts=align.module_output.counts,
        gtf=config["expression"]["gtf"],
        bed=config["expression"].get("bed", []),
        src=workflow.source_path("includes/expression/scripts/coverage.py"),


use rule transform_counts from expression as expression_transform_counts with:
    input:
        counts=align.module_output.counts,
        src=workflow.source_path("includes/expression/scripts/transform_counts.py"),


use rule multiqc from expression as expression_multiqc with:
    output:
        html=temporary("multiqc_expression.html"),
        folder=directory("multiqc_expression_data"),
        parquet="multiqc_expression_data/multiqc.parquet",
        filelist=temporary("multiqc_filelist_expression.txt"),


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        idm=align.module_output.id_mapping,
        fusion_json=fusion.module_output.json,
        snv_indels_json=align.module_output.json,
        itd_json=itd.module_output.json,
        expression_json=expression.module_output.json,
        scr=workflow.source_path("scripts/create_summary.py"),
    params:
        pipeline_ver=PIPELINE_VERSION,
    output:
        js="{sample}/{sample}.summary.json",
    log:
        "log/create_summary.{sample}.txt",
    container:
        containers["hamlet-scripts"]
    shell:
        """
        python {input.scr} \
            {input.idm} \
            --pipeline-version {params.pipeline_ver} \
            --sample-name {wildcards.sample} \
            --module {input.fusion_json} \
            --module {input.snv_indels_json} \
            --module {input.expression_json} \
            --module {input.itd_json} > {output.js} 2>{log}
        """


rule generate_report:
    """Generates a PDF report of the essential results."""
    input:
        summary=rules.create_summary.output.js,
        css=workflow.source_path("report/assets/style.css"),
        toc=workflow.source_path("report/assets/toc.xsl"),
        scr=workflow.source_path("scripts/generate_report.py"),
        # Ensure all report files are localised to the stupid Snakemake cache
        report_files=report_files,
    output:
        "{sample}/hamlet_report.{sample}.pdf",
    log:
        "log/generate_report.{sample}.txt",
    container:
        containers["hamlet-scripts"]
    shell:
        """
        # Find the cached location of the report folder
        css={input.css}
        report=${{css%assets/style.css}}
        templates="${{report}}templates"
        imgs="${{report}}assets/img"

        python3 {input.scr} \
            --templates-dir ${{templates}} \
            --imgs-dir ${{imgs}} \
            --css-path {input.css} \
            --toc-path {input.toc} \
            {input.summary} \
            --pdf-output {output} 2> {log}
        """


rule generate_html_report:
    """Generates a HTML report of the essential results, used for testing only"""
    input:
        summary=rules.create_summary.output.js,
        css=workflow.source_path("report/assets/style.css"),
        toc=workflow.source_path("report/assets/toc.xsl"),
        scr=workflow.source_path("scripts/generate_report.py"),
        # Ensure all report files are localised to the stupid Snakemake cache
        report_files=report_files,
    output:
        "{sample}/hamlet_report.{sample}.html",
    log:
        "log/generate_html_report.{sample}.txt",
    container:
        containers["hamlet-scripts"]
    shell:
        """
        # Find the cached location of the report folder
        css={input.css}
        report=${{css%assets/style.css}}
        templates="${{report}}templates"
        imgs="${{report}}assets/img"

        python3 {input.scr} \
            --templates-dir ${{templates}} \
            --imgs-dir ${{imgs}} \
            --css-path {input.css} \
            --toc-path {input.toc} \
            {input.summary} \
            --html-output {output} 2> {log}
        """


rule multiqc:
    input:
        qc_stats=qc_seq.module_output.multiqc_parquet,
        snv_indel_stats=align.module_output.multiqc_parquet,
        expression_stats=expression.module_output.multiqc_parquet,
        config=workflow.source_path("cfg/multiqc.yml"),
        background=config.get("background_samples", []),
    params:
        exclude="dedup",
    output:
        html="multiqc_hamlet.html",
    log:
        "log/multiqc.txt",
    container:
        containers["multiqc"]
    shell:
        """
        multiqc \
        --force \
        --config {input.config} \
        --exclude {params.exclude} \
        --filename {output.html} \
        {input.qc_stats} \
        {input.snv_indel_stats} \
        {input.expression_stats} \
        {input.background} \
        2> {log}
        """


rule cleanup_multiqc:
    """Rule to clean up leftover MultiQC output files and folders from the modules"""
    input:
        qc_stats=qc_seq.module_output.multiqc_parquet,
        snv_indel_stats=align.module_output.multiqc_parquet,
        expression_stats=expression.module_output.multiqc_parquet,
        run_after_multiqc_rule=rules.multiqc.output.html,
    output:
        target=".cleanup_multiqc",
    log:
        "log/cleanup_multiqc.txt",
    container:
        containers["multiqc"]
    shell:
        """
        rm -rfv $(dirname {input.qc_stats}) 2>&1 >> {log}
        rm -rfv $(dirname {input.snv_indel_stats}) 2>&1 >> {log}
        rm -rfv $(dirname {input.expression_stats}) 2>&1 >> {log}

        touch {output.target}
        """
