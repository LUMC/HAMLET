from types import SimpleNamespace

containers = {
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:4.6--py39hf95cd2a_1",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0",
}


pepfile: config["pepfile"]


# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]


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


module_output = SimpleNamespace(
    forward=get_forward_output,
    reverse=get_reverse_output,
)
