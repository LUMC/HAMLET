from types import SimpleNamespace


pepfile: config["pepfile"]


# If we run with the full hamlet configuration, subset the configuration
if "itd" in config:
    config = config["itd"]

# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


containers = {
    "bwa-0.7.17-samtools-1.3.1-picard-2.9.2": "docker://quay.io/biocontainers/mulled-v2-1c6be8ad49e4dfe8ab70558e8fb200d7b2fd7509:5900b4e68c4051137fffd99165b00e98f810acae-0",
    "rose": "docker://quay.io/redmar_van_den_berg/rose-dt:0.4",
}


def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R2"]


## Functions for module outputs ##
def get_csv(wildcards, gene):
    return f"{wildcards.sample}/itd/{wildcards.sample}.{gene}.csv"


def get_plot(wildcards, gene):
    return f"{wildcards.sample}/itd/{wildcards.sample}.{gene}.png"


def get_kmt2a(wildcards):
    return get_csv(wildcards, "kmt2a")


def get_kmt2a_plot(wildcards):
    return get_plot(wildcards, "kmt2a")


def get_flt3(wildcards):
    return get_csv(wildcards, "flt3")


def get_flt3_plot(wildcards):
    return get_plot(wildcards, "flt3")


def get_json(wildcards):
    return f"{wildcards.sample}/itd/itd-output.json"


module_output = SimpleNamespace(
    kmt2a_csv=get_kmt2a,
    kmt2a_plot=get_kmt2a_plot,
    flt3_csv=get_flt3,
    flt3_plot=get_flt3_plot,
    json=get_json,
)
