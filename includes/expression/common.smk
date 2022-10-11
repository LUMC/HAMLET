from types import SimpleNamespace

containers = {
    "bedtools-2.17-python-2.7": "docker://quay.io/biocontainers/mulled-v2-a9ddcbd438a66450297b5e0b61ac390ee9bfdb61:e60f3cfda0dfcf4a72f2091c6fa1ebe5a5400220-0",
    "htseq": "docker://quay.io/biocontainers/htseq:0.11.2--py27h637b7d7_1",
    "picard": "docker://quay.io/biocontainers/picard:2.20.5--0",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.2",
}


def set_default(key, value):
    """Set default value for settings"""
    if key not in config:
        config[key] = value


set_default("base_count_script", srcdir("scripts/hist2count.py"))
set_default("aggr_base_count_script", srcdir("scripts/aggr_base_count.R"))
set_default("calc_ratio_script", srcdir("scripts/calc_ratio.py"))
set_default("relative_gene_name", "HMBS")

if "exon_names" not in config:
    raise ValueError("No exon names for exon ratio calculation defined")


def get_bamfile(wildcards):
    return pep.sample_table.loc[wildcards.sample, "bam"]


## Functions for module outputs ##


def output_file(wildcards, extension):
    return f"{wildcards.sample}/expression/{wildcards.sample}.{extension}"


output_files = [
    "bases_per_exon",
    "bases_per_gene",
    "exon_ratios",
    "fragments_per_gene",
    "raw_base",
]

output_functions = {
    out_file: partial(output_file, extension=out_file) for out_file in output_files
}

module_output = SimpleNamespace(**output_functions)
