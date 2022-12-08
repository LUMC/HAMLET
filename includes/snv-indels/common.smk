from types import SimpleNamespace


pepfile: config["pepfile"]


containers = {
    "bedtools-2.27-grep-2.14-gawk-5.0-python-3.7": "docker://quay.io/biocontainers/mulled-v2-a4b89e0b16b1d7db92e5a069e5c40405b3b53aab:98c4ac2f0e27869be58f6a4d8bb7ae3bc02a3a70-0",
    "gsnap": "docker://quay.io/biocontainers/gmap:2021.08.25--pl5321h67092d7_2",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.2",
    "picard": "docker://quay.io/biocontainers/picard:2.27.4--hdfd78af_0",
    "varscan-2.4.2-samtools-1.3.1-tabix-0.2.6-grep-2.14": "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0",
    "crimson": "docker://quay.io/biocontainers/crimson:1.1.0--pyh5e36f6f_0",
}


def set_default(key, value):
    """Set default value for settings"""
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


## Functions for module outputs ##
def get_bam_output(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.snv-indel.bam"


def get_variant_plot_dir(wildcards):
    return f"{wildcards.sample}/snv-indels/variant_plots"


def get_all_csv(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.variants_all.csv"


def get_high_csv(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.variants_hi.csv"


def get_json(wildcards):
    return f"{wildcards.sample}/snv-indels/snv-indels-output.json"


module_output = SimpleNamespace(
    bam=get_bam_output,
    variant_plot_dir=get_variant_plot_dir,
    var_all=get_all_csv,
    var_hi=get_high_csv,
    json=get_json,
)
