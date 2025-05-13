from uuid import uuid4


pepfile: config["pepfile"]


samples = pep.sample_table["sample_name"]

for s in samples:
    if " " in s:
        raise RuntimeError(f'Spaces in samples are not supported ("{s}")')


containers = {
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.3",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0",
}


def report_files(wildcards):
    files = [
        "report/templates/contents_about.html.j2",
        "report/templates/contents_aln.html.j2",
        "report/templates/contents_basic.html.j2",
        "report/templates/contents_fusion.html.j2",
        "report/templates/contents_var.html.j2",
        "report/templates/cover.html.j2",
        "report/templates/contents.html.j2",
        "report/templates/contents_itd.html.j2",
        "report/templates/contents_rna.html.j2",
        "report/templates/contents_overview.html.j2",
        "report/assets/style.css",
        "report/assets/toc.xsl",
        "report/assets/img/hamlet-logo.jpg",
        "report/assets/img/lumc-logo.jpg",
        "report/assets/img/Myeloblast_logo.png",
    ]

    return [workflow.source_path(fname) for fname in files]


# The version of HAMLET
PIPELINE_VERSION = "v2.3.2-dev"
