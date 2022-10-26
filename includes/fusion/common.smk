from types import SimpleNamespace


pepfile: config["pepfile"]


containers = {
    "debian": "docker://debian:buster-slim",
    "fsnviz": "docker://quay.io/biocontainers/fsnviz:0.3.0--py_4",
    "fuma": "docker://quay.io/biocontainers/fuma:4.0.0--pyhb7b1952_0",
    "fusioncatcher": "docker://quay.io/biocontainers/fusioncatcher:1.20--2",
    "star-fusion": "docker://quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1",
    "crimson": "docker://quay.io/biocontainers/crimson:1.1.0--pyh5e36f6f_0",
}


def set_default(key, value):
    """Set default value for settings"""
    if key not in config:
        config[key] = value


set_default("sf_subset_script", srcdir("scripts/subset_sf.py"))
set_default("plot_combined_script", srcdir("scripts/combine_svgs.py"))
set_default("sf_rewrite_script", srcdir("scripts/rewrite_star_fusion_header.py"))
set_default("fusioncatcher_data", False)


def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


## Functions for module outputs ##
def star_fusion(wildcards):
    return f"{wildcards.sample}/fusion/{wildcards.sample}.star-fusion"


def star_fusion_fig(wildcards, ext):
    return (
        f"{wildcards.sample}/fusion/{wildcards.sample}.star-fusion-circos/fsnviz.{ext}"
    )


def json(wildcards):
    return f"{wildcards.sample}/fusion/fusion-output.json"


def fusion_catcher(wildcards):
    if config["fusioncatcher_data"]:
        return f"{wildcards.sample}/fusion/{wildcards.sample}.fusioncatcher"


def fusion_catcher_fig(wildcards, ext):
    if config["fusioncatcher_data"]:
        return f"{wildcards.sample}/fusion/{wildcards.sample}.fusioncatcher.{ext}"
    else:
        return list()


def intersect(wildcards):
    if config["fusioncatcher_data"]:
        return f"{wildcards.sample}/fusion/{wildcards.sample}.fuma"


def intersect_fig(wildcards):
    if config["fusioncatcher_data"]:
        return (
            f"{wildcards.sample}/fusion/{wildcards.sample}.sf-isect-circos/fsnviz.png"
        )


def subset_predictions(wildcards):
    if config["fusioncatcher_data"]:
        return f"{wildcards.sample}/fusion/{wildcards.sample}.sf-isect"


# Optional module outputs
optional = SimpleNamespace(
    fusion_catcher=fusion_catcher if config["fusioncatcher_data"] else lambda x: [],
    intersect=intersect if config["fusioncatcher_data"] else lambda x: [],
    intersect_fig=intersect_fig if config["fusioncatcher_data"] else lambda x: [],
    subset_predictions=subset_predictions
    if config["fusioncatcher_data"]
    else lambda x: [],
)
module_output = SimpleNamespace(
    star_fusion=star_fusion,
    star_fusion_fig=star_fusion_fig,
    json=json,
    fusion_catcher_fig=fusion_catcher_fig,
    optional=optional,
)
