from types import SimpleNamespace
from scripts import gtf

containers = {
    "pysam": "docker://quay.io/biocontainers/pysam:0.22.1--py39h61809e1_2",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0",
    "seamless": "docker://quay.io/biocontainers/r-seamless:0.1.1--r44h3121a25_1",
    "amlmapr": "docker://quay.io/biocontainers/r-amlmapr:0.1.0--r44hdfd78af_0",
}


pepfile: config["pepfile"]


# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


def str_to_list(value):
    """Convert a space separated str to a list"""
    if value is None:
        return list()
    elif isinstance(value, str):
        return value.split(" ")
    elif isinstance(value, list):
        return value
    else:
        raise RuntimeError


def set_genes_of_interest():
    """Set the genes of interest to a list, if it is a string"""
    # Genes of interest is either a list of str, or a string with spaces
    goi = config.get("genes_of_interest")
    config["genes_of_interest"] = str_to_list(goi)


def set_housekeeping_genes():
    housekeeping = config.get("housekeeping")
    config["housekeeping"] = str_to_list(housekeeping)


def set_pdf_report_genes():
    report_genes = config.get("report")
    config["report"] = str_to_list(report_genes)


## Input functions ##
def get_bam(wildcards):
    return pep.sample_table.loc[wildcards.sample, "bam"]


def get_counts(wildcards):
    return pep.sample_table.loc[wildcards.sample, "count"]


def get_strand(wildcards):
    try:
        strandedness = pep.sample_table.loc[wildcards.sample, "strandedness"]
    except KeyError:
        strandedness = "unstranded"
    if strandedness is None:
        raise ValueError("Please specify a strandedness for every sample")
    return strandedness


## Check specified settings ##
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


def _bed_names():
    """Get gene names from the Bed file"""
    if "bed" not in config:
        return list()

    names = list()
    with open(config["bed"]) as fin:
        for line in fin:
            spline = line.strip("\n").split("\t")
            names.append(spline[3])
    return names


def check_bed_genes_of_interest():
    """Check for duplicates genes bed file and genes_of_interest"""
    # Get the names from the bed file
    names = _bed_names()

    for gene in config["genes_of_interest"]:
        if gene in names:
            msg = f"{gene} is specified twice: in the bed file and genes_of_interest"
            raise RuntimeError(msg)


def check_genes_for_report():
    quantified_genes = _bed_names() + config["genes_of_interest"]
    for gene in config["report"]:
        if gene not in quantified_genes:
            msg = f"{gene} is specified for the report, but not in genes_of_interest or the bed file"
            raise RuntimeError(msg)


## Functions for module outputs ##
def coverage(wildcards):
    return f"{wildcards.sample}/expression/coverage.csv"


def normalized(wildcards):
    return f"{wildcards.sample}/expression/coverage.normalized.csv"


def multiqc_files():
    unstranded = ("merged_expression_unstranded_mqc.tsv",)
    stranded = ("merged_expression_stranded_mqc.tsv",)
    cell_types = ("merged_expression_cell_types_mqc.tsv",)

    return unstranded + stranded + cell_types


def get_json(wildcards):
    return f"{wildcards.sample}/expression/expression-output.json"


# Set default values, parse specified values
set_genes_of_interest()
set_housekeeping_genes()
set_pdf_report_genes()

# Check the sanity of the reported settings
check_housekeeping()
check_bed_genes_of_interest()
check_genes_for_report()

module_output = SimpleNamespace(
    coverage=coverage,
    normalized_expression=normalized,
    multiqc_files=multiqc_files(),
    json=get_json,
)
