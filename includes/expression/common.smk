from types import SimpleNamespace
from scripts import gtf

containers = {}


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


def check_housekeeping():
    """Check if can find each housekeeping gene in the GTF file"""
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


check_housekeeping()

module_output = SimpleNamespace()
