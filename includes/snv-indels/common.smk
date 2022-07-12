pepfile: config["pepfile"]

samples = pep.sample_table["sample_name"]

containers = {
    "bedtools-2.27-grep-2.14-gawk-5.0-click-7-python-3.7": "docker://quay.io/biocontainers/mulled-v2-a4b89e0b16b1d7db92e5a069e5c40405b3b53aab:98c4ac2f0e27869be58f6a4d8bb7ae3bc02a3a70-0",
    "debian": "docker://debian:buster-slim",
    "gsnap": "docker://quay.io/biocontainers/gmap:2020.06.30--pl526h2f06484_0",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.2",
    "picard": "docker://quay.io/biocontainers/picard:2.20.5--0",
    "python3": "docker://python:3.7.4-slim-stretch",
    "samtools": "docker://quay.io/biocontainers/samtools:1.6--h244ad75_4",
    "varscan-2.4.2-samtools-1.3.1-tabix-0.2.6-grep-2.14":
    "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:101.0--pl526hecda079_1"
}

def set_default(key, value):
    """ Set default value for settings """
    if key not in config:
        config[key] = value

set_default("genome_dict", config["genome_fasta"].rsplit(".", 1)[0] + ".dict")
set_default("genome_fai", config["genome_fasta"] + ".fai")
set_default("exon_cov_script", srcdir("scripts/aggr_exon_cov.py"))
set_default("extract_script", srcdir("scripts/vcf2json.py"))
set_default("csv_script", srcdir("scripts/json2csv.py"))
set_default("plot_script", srcdir("scripts/plotVariants.R"))

def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]

def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]
