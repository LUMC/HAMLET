from types import SimpleNamespace


pepfile: config["pepfile"]


containers = {
    "bedtools-2.27-grep-2.14-gawk-5.0-python-3.7": "docker://quay.io/biocontainers/mulled-v2-a4b89e0b16b1d7db92e5a069e5c40405b3b53aab:98c4ac2f0e27869be58f6a4d8bb7ae3bc02a3a70-0",
    "picard": "docker://quay.io/biocontainers/picard:2.27.4--hdfd78af_0",
    "varscan-2.4.2-samtools-1.3.1-tabix-0.2.6-grep-2.14": "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
    "vardict": "docker://quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0",
    "star": "docker://quay.io/biocontainers/star:2.7.10b--h9ee0642_0",
    "crimson": "docker://quay.io/biocontainers/crimson:1.1.0--pyh5e36f6f_0",
}


def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R2"]


## Functions for module outputs ##
def get_bam_output(wildcards):
    return f"{wildcards.sample}/snv-indels/Aligned.sortedByCoord.out.bam"


def get_bai_output(wildcards):
    return get_bam_output(wildcards) + ".bai"


def get_vep_high(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.vep.high.txt.gz"


def get_vep_target(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.vep.target.txt.gz"


def get_json(wildcards):
    return f"{wildcards.sample}/snv-indels/snv-indels-output.json"


def get_hotspot(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.hotspot.vcf"


module_output = SimpleNamespace(
    bam=get_bam_output,
    bai=get_bai_output,
    vep_high=get_vep_high,
    vep_target=get_vep_target,
    json=get_json,
    hotspot=get_hotspot,
)
