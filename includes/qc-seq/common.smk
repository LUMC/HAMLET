from os import path

containers = {
    "crimson": "docker://quay.io/biocontainers/crimson:0.3.0--py27_1",
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:1.14--py36_0",
    "debian": "docker://debian:buster-slim",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.11.8--1"
}

config["rg_stats_script"] = srcdir(path.join("scripts", "gather_rg_stats.py"))
config["sample_stats_script"] = srcdir(path.join("scripts", "gather_sample_stats.py"))

def get_r(strand, wildcards):
    """Get fastq files on a single strand for a sample"""
    s = config["samples"].get(wildcards.sample)
    rs = []
    for rg in sorted(s["read_groups"].keys()):
        rs.append(s["read_groups"][rg][strand])
    return rs

get_r1 = partial(get_r, "R1")
get_r2 = partial(get_r, "R2")
