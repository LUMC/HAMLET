pepfile: config["pepfile"]

samples = pep.sample_table["sample_name"]

containers = {
    "debian": "docker://debian:buster-slim",
    "fsnviz": "docker://quay.io/biocontainers/fsnviz:0.3.0--py_4",
    "fuma": "docker://quay.io/biocontainers/fuma:4.0.0--pyhb7b1952_0",
    "fusioncatcher": "docker://quay.io/biocontainers/fusioncatcher:1.20--2",
    "star-fusion": "docker://quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.2"
}

def set_default(key, value):
    """ Set default value for settings """
    if key not in config:
        config[key] = value

set_default("sf_subset_script", srcdir("scripts/subset_sf.py"))
set_default("plot_combined_script", srcdir("scripts/combine_svgs.py"))
set_default("sf_rewrite_script", srcdir("scripts/rewrite_star_fusion_header.py"))
set_default("fusioncatcher_data", False)

def get_fusioncatcher_outputs():
    # If no fusioncatcher data is specified, we do not run it and there are no
    # outputs
    if not config["fusioncatcher_data"]:
        return list()

    fs = expand("{sample}/fusion/{sample}.fusioncatcher", sample=samples)
    fuma = expand("{sample}/fusion/{sample}.fuma", sample=samples)
    is_svg = expand("{sample}/fusion/{sample}.sf-isect-circos/fsnviz.svg", sample=samples)
    is_png = expand("{sample}/fusion/{sample}.sf-isect-circos/fsnviz.png", sample=samples)
    fc_svg = expand("{sample}/fusion/{sample}.fusioncatcher-circos/fsnviz.svg", sample=samples)
    fc_png = expand("{sample}/fusion/{sample}.fusioncatcher-circos/fsnviz.png", sample=samples)
    combined_svg = expand("{sample}/fusion/{sample}.fusions-combined.svg", sample=samples)
    subset = expand("{sample}/fusion/{sample}.sf-isect", sample=samples)
    return fs + fuma + is_svg + is_png + fc_svg + fc_png + combined_svg + subset

def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]

def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]
