from types import SimpleNamespace
from scripts import gtf

containers = {
    "pysam": "docker://quay.io/biocontainers/pysam:0.22.1--py39h61809e1_2",
}


pepfile: config["pepfile"]


# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


def get_bam(wildcards):
    return pep.sample_table.loc[wildcards.sample, "bam"]


def check_housekeeping():
    """Check if we can find each housekeeping gene in the GTF file"""
    # Read the mapping from ENSG to gene name
    with open(config["gtf"]) as fin:
        ensg_to_name = gtf.gene_id_name(fin)

    # Create mapping from name to ENSG
    name_to_ensg = {v: k for k, v in ensg_to_name.items()}

    # Raise an error on unknown genes
    for gene in config["housekeeping"]:
        if gene not in name_to_ensg:
            msg = f"Unknown housekeeping gene: {gene}"
            raise RuntimeError(msg)


## Functions for module outputs ##
def coverage(wildcards):
    return f"{wildcards.sample}/expression/coverage.csv"


check_housekeeping()

module_output = SimpleNamespace(coverage=coverage)
