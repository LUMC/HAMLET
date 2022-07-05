containers = {
    "bedtools-2.17-python-2.7": "docker://quay.io/biocontainers/mulled-v2-a9ddcbd438a66450297b5e0b61ac390ee9bfdb61:e60f3cfda0dfcf4a72f2091c6fa1ebe5a5400220-0",
    "htseq": "docker://quay.io/biocontainers/htseq:0.11.2--py27h637b7d7_1",
    "picard": "docker://quay.io/biocontainers/picard:2.20.5--0",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.2"
}

settings=config["settings"]

# Set the default settings
def set_default(key, value):
    """ Set default value for settings """
    if key not in settings:
        settings[key] = value

set_default("base_count_script", srcdir("scripts/hist2count.py"))
set_default("aggr_base_count_script", srcdir("scripts/aggr_base_count.R"))
set_default("calc_ratio_script", srcdir("scripts/calc_ratio.py"))
set_default("relative_gene_name", "HMBS")

if "exon_names" not in settings:
    raise ValueError("No exon names for exon ratio calculation defined")

def get_bamfile(wildcards):
    print(pep.sample_table)
    print(wildcards)
    return pep.sample_table.loc[wildcards.sample, "bam"]
