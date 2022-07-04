from os import path

containers = {
    "crimson": "docker://quay.io/biocontainers/crimson:0.3.0--py27_1",
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:1.14--py36_0",
    "debian": "docker://debian:buster-slim",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.11.8--1"
}

config["rg_stats_script"] = srcdir(path.join("scripts", "gather_rg_stats.py"))
config["sample_stats_script"] = srcdir(path.join("scripts", "gather_sample_stats.py"))

def get_forward(wildcards):
    return get_fastq(wildcards, "R1")


def get_reverse(wildcards):
    return get_fastq(wildcards, "R2")


def get_fastq(wildcards, pair):
    fastq = pep.sample_table.loc[wildcards.sample, pair]

    # If a single fastq file is specified, we put it in a list
    if isinstance(fastq, str):
        fastq = [fastq]
        nr_fastq = 1
    # If multiple fastq files were specified, it is already a list
    else:
        nr_fastq = len(fastq)

    # Here, we use the readgroup wildcard to pick the correct fastq file
    read_groups = {f"rg_{rg+1}": fastq for rg, fastq in zip(range(len(fastq)), fastq)}
    return read_groups[wildcards.read_group]

def get_all_raw_fastq(wildcards):
    return pep.sample_table.loc[wildcards.sample, wildcards.pair]

def get_all_trimmed_fastq(wildcards, pair):
    trimmed_files = list()
    for rg, sample in get_readgroup_per_sample():
        if sample == wildcards.sample:
            trimmed_files.append(f"{sample}/qc-seq/{rg}/{sample}-{rg}-{pair}.fq.gz")
    return trimmed_files

def get_all_trimmed_forward(wildcards):
    return get_all_trimmed_fastq(wildcards, "R1")

def get_all_trimmed_reverse(wildcards):
    return get_all_trimmed_fastq(wildcards, "R2")

def get_r1(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]

def get_r2(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R2"]

def get_readgroup_per_sample():
    for sample in pep.sample_table["sample_name"]:
        fastq = pep.sample_table.loc[sample, "R1"]
        # If there is only a single fastq, put it in a list
        if isinstance(fastq, str):
            fastq = [fastq]
        nr_readgroups = len(fastq)
        for i in range(nr_readgroups):
            yield f"rg_{i+1}", sample

def get_readgroup(wildcards):
    stat_files = list()
    for read_group, sample in get_readgroup_per_sample():
        if sample == wildcards.sample:
            stat_files.append(f"{sample}/qc-seq/{read_group}/stats.json")
    return stat_files
