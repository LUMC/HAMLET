include: "common.smk"


localrules:
    align_exon_cov_ref,
    align_genome_txt,
    align_json_output,
    align_table_vars_all,
    align_table_vars_hi,
    create_summary,
    fusion_arriba_to_json,
    generate_report,
    itd_detect_itd_ftl3,
    itd_detect_itd_kmt2a,
    itd_json_output,
    itd_plot_itd_flt3,
    itd_plot_itd_kmt2a,
    package_results,


rule all:
    input:
        summary=expand("{sample}/{sample}.summary.json", sample=samples),
        report=expand("{sample}/hamlet_report.{sample}.pdf", sample=samples),


# Add the PEP configuration to each submodule
config["qc-seq"]["pepfile"] = config["pepfile"]
config["snv-indels"]["pepfile"] = config["pepfile"]
config["itd"]["pepfile"] = config["pepfile"]
config["fusion"]["pepfile"] = config["pepfile"]


# Define HAMLET modules
module qc_seq:
    snakefile:
        "includes/qc-seq/Snakefile"
    config:
        config["qc-seq"]


use rule * from qc_seq as qc_seq_*


# Make the trimmed, merged FastQ files temporary
use rule merge_fastqs_r1 from qc_seq as qc_seq_merge_fastqs_r1 with:
    output:
        merged=temp("{sample}/{sample}-R1.fq.gz"),


use rule merge_fastqs_r2 from qc_seq as qc_seq_merge_fastqs_r2 with:
    output:
        merged=temp("{sample}/{sample}-R2.fq.gz"),


module itd:
    snakefile:
        "includes/itd/Snakefile"
    config:
        config["itd"]


use rule * from itd as itd_*


# Connect the align_kmt2a rule to the output of qc-seq
use rule align_kmt2a from itd as itd_align_kmt2a with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        fasta=config["itd"]["kmt2a_fasta"],


# Connect the align_flt3 rule to the output of qc-seq
use rule align_flt3 from itd as itd_align_flt3 with:
    input:
        fq1=qc_seq.module_output.forward,
        fq2=qc_seq.module_output.reverse,
        fasta=config["itd"]["flt3_fasta"],


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


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        idm=config["snv-indels"]["ref_id_mapping"],
        qc_seq_json=qc_seq.module_output.json,
        fusion_json=fusion.module_output.json,
        snv_indels_json=align.module_output.json,
        itd_json=itd.module_output.json,
        scr=srcdir("scripts/create_summary.py"),
    params:
        pipeline_ver=PIPELINE_VERSION,
        run_name=RUN_NAME,
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
            --run-name {params.run_name} \
            --sample-name {wildcards.sample} \
            --module {input.fusion_json} \
            --module {input.snv_indels_json} \
            --module {input.itd_json} \
            --module {input.qc_seq_json} > {output.js} 2>{log}
        """


rule generate_report:
    """Generates a PDF report of the essential results."""
    input:
        summary=rules.create_summary.output.js,
        css=srcdir("report/assets/style.css"),
        templates=srcdir("report/templates"),
        imgs=srcdir("report/assets/img"),
        toc=srcdir("report/assets/toc.xsl"),
        scr=srcdir("scripts/generate_report.py"),
    output:
        pdf="{sample}/hamlet_report.{sample}.pdf",
        html="{sample}/hamlet_report.{sample}.html",
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
            --html-output {output.html} \
            --pdf-output {output.pdf} 2> {log}
        """
