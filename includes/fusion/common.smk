from types import SimpleNamespace


pepfile: config["pepfile"]


containers = {
    "fsnviz": "docker://quay.io/biocontainers/fsnviz:0.3.0--pyhdfd78af_6",
    "fuma": "docker://quay.io/biocontainers/fuma:4.0.0--pyhb7b1952_0",
    "fusioncatcher": "docker://quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_4",
    "star-fusion": "docker://quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1",
    "crimson": "docker://quay.io/biocontainers/crimson:1.1.0--pyh5e36f6f_0",
}


def set_default(key, value):
    """Set default value for settings"""
    if key not in config:
        config[key] = value


set_default("sf_subset_script", srcdir("scripts/subset_sf.py"))
set_default("fusioncatcher_data", False)


def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R2"]


## Functions for module outputs ##
def star_fusion(wildcards):
    return f"{wildcards.sample}/fusion/{wildcards.sample}.star-fusion"


def star_fusion_fig(wildcards):
    return f"{wildcards.sample}/fusion/{wildcards.sample}.star-fusion-circos/fsnviz.png"


def json(wildcards):
    return f"{wildcards.sample}/fusion/fusion-output.json"


def fusion_catcher(wildcards):
    if config["fusioncatcher_data"]:
        return f"{wildcards.sample}/fusion/{wildcards.sample}.fusioncatcher"


def fusioncatcher_fig(wildcards):
    if config["fusioncatcher_data"]:
        sample = wildcards.sample
        return f"{sample}/fusion/{sample}.fusioncatcher-circos/fsnviz.png"


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


# Optional module outputs if they exist
if config["fusioncatcher_data"]:
    optional = SimpleNamespace(
        fusion_catcher=fusion_catcher,
        intersect=intersect,
        intersect_fig=intersect_fig,
        fusioncatcher_fig=fusioncatcher_fig,
        subset_predictions=subset_predictions,
    )
# If the optional outputs do not exist, they should be a function that returns
# an empty list. This is needed for compatibility with Snakemake
else:
    optional = SimpleNamespace(
        fusion_catcher=lambda x: list(),
        intersect=lambda x: list(),
        intersect_fig=lambda x: list(),
        fusioncatcher_fig=lambda x: list(),
        subset_predictions=lambda x: list(),
    )
module_output = SimpleNamespace(
    star_fusion=star_fusion,
    star_fusion_fig=star_fusion_fig,
    json=json,
    optional=optional,
)
