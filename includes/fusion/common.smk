from types import SimpleNamespace


pepfile: config["pepfile"]


# If we run with the full hamlet configuration, subset the configuration
if "fusion" in config:
    config = config["fusion"]


# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


containers = {
    "arriba": "docker://quay.io/biocontainers/arriba:2.4.0--h0033a41_2",
    "poppler": "docker://quay.io/biocontainers/keggcharter:0.6.0--hdfd78af_0",
}


def get_bam(wildcards):
    return pep.sample_table.loc[wildcards.sample, "bam"]


def get_bai(wildcards):
    return f"{get_bam(wildcards)}.bai"


## Functions for module outputs ##
def json(wildcards):
    return f"{wildcards.sample}/fusion/fusion-output.json"


def arriba(wildcards):
    return f"{wildcards.sample}/fusion/arriba/fusions.tsv"


module_output = SimpleNamespace(arriba=arriba, json=json)
