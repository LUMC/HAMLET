from types import SimpleNamespace

containers = {
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:5.0--py39hbcbf7aa_0",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.22.1--pyhdfd78af_0",
    "sequali": "docker://quay.io/biocontainers/sequali:0.12.0--py311haab0aaa_1",
}


pepfile: config["pepfile"]


# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


def get_input_fastq(sample, pair):
    """Get all fastq files for the specified sample"""
    fastq = pep.sample_table.loc[sample, pair]
    # If a single fastq file is specified, we put it in a list
    if isinstance(fastq, str):
        fastq = [fastq]
    return fastq


def get_forward_input(wildcards):
    """Get the input forward fastq file"""
    return get_input_fastq(wildcards.sample, "R1")


def get_reverse_input(wildcards):
    """Get the input reverse fastq file"""
    return get_input_fastq(wildcards.sample, "R2")


def get_forward_output(wildcards):
    return get_output_fastq(wildcards, "R1")


def get_reverse_output(wildcards):
    return get_output_fastq(wildcards, "R2")


def get_output_fastq(wildcards, pair):
    return f"{wildcards.sample}/qc-seq/{wildcards.sample}-{pair}.fq.gz"


def multiqc_files():
    cutadapt = [
        f"{wildcards.sample}/qc-seq/{wildcards.sample}.cutadapt.json"
        for wildcards in samples
    ]
    sequali = [
        f"{wildcards.sample}/qc-seq/sequali/{wildcards.sample}.json"
        for wildcards in samples
    ]
    return cutadapt + sequali


def multiqc_modules():
    """Define which MultiQC modules to run here"""
    modules = ["cutadapt", "sequali"]
    return [f" --module {module}" for module in modules]


module_output = SimpleNamespace(
    forward=get_forward_output,
    reverse=get_reverse_output,
    multiqc_files=multiqc_files(),
)
