include: "common.smk"


localrules:
    align_exon_cov_ref,
    align_genome_txt,
    align_json_output,
    align_table_vars_hi,
    create_summary,
    fusion_arriba_to_json,
    generate_report,
    itd_detect_itd_ftl3,
    itd_detect_itd_kmt2a,
    itd_json_output,
    itd_plot_itd_flt3,
    itd_plot_itd_kmt2a,


rule all:
    input:
        summary=expand("{sample}/{sample}.summary.json", sample=samples),
        report=expand("{sample}/hamlet_report.{sample}.pdf", sample=samples),
        multiqc="multiqc_hamlet.html",


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


# Make the final trimmed, merged FastQ files from the qc_seq module temporary
use rule cutadapt from qc_seq as qc_seq_cutadapt with:
    output:
        fq1=temp("{sample}/qc-seq/{sample}-R1.fq.gz"),
        fq2=temp("{sample}/qc-seq/{sample}-R2.fq.gz"),
        json="{sample}/qc-seq/{sample}.cutadapt.json",


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


rule create_summary:
    """Combines statistics and other info across modules to a single JSON file per sample."""
    input:
        idm=config["snv-indels"]["ref_id_mapping"],
        fusion_json=fusion.module_output.json,
        snv_indels_json=align.module_output.json,
        itd_json=itd.module_output.json,
        scr=srcdir("scripts/create_summary.py"),
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
            --module {input.itd_json} > {output.js} 2>{log}
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
        "{sample}/hamlet_report.{sample}.pdf",
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
            --pdf-output {output} 2> {log}
        """


rule generate_html_report:
    """Generates a HTML report of the essential results, used for testing only"""
    input:
        summary=rules.create_summary.output.js,
        css=srcdir("report/assets/style.css"),
        templates=srcdir("report/templates"),
        imgs=srcdir("report/assets/img"),
        toc=srcdir("report/assets/toc.xsl"),
        scr=srcdir("scripts/generate_report.py"),
    output:
        "{sample}/hamlet_report.{sample}.html",
    log:
        "log/generate_html_report.{sample}.txt",
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
            --html-output {output} 2> {log}
        """


rule multiqc:
    input:
        qc_stats=qc_seq.module_output.multiqc_files,
        snv_indel_stats=align.module_output.multiqc_files,
        config=srcdir("cfg/multiqc.yml"),
    params:
        filelist="multiqc_filelist.txt",
        depth=2,
    output:
        html="multiqc_hamlet.html",
    log:
        "log/multiqc.txt",
    container:
        containers["multiqc"]
    shell:
        """
        rm -f {params.filelist}

        for fname in {input.qc_stats} {input.snv_indel_stats}; do
            echo $fname >> {params.filelist}
        done

        multiqc \
        --force \
        --dirs \
        --dirs-depth {params.depth} \
        --fullnames \
        --fn_as_s_name \
        --file-list {params.filelist} \
        --config {input.config} \
        --filename {output.html} 2> {log}
        """
