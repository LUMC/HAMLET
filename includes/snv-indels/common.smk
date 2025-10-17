from types import SimpleNamespace


pepfile: config["pepfile"]


containers = {
    "bedtools-2.27-grep-2.14-gawk-5.0-python-3.7": "docker://quay.io/biocontainers/mulled-v2-a4b89e0b16b1d7db92e5a069e5c40405b3b53aab:98c4ac2f0e27869be58f6a4d8bb7ae3bc02a3a70-0",
    "picard": "docker://quay.io/biocontainers/picard:3.3.0--hdfd78af_0",
    "varscan-2.4.2-samtools-1.3.1-tabix-0.2.6-grep-2.14": "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
    "vardict": "docker://quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1",
    "star": "docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_5",
    "crimson": "docker://quay.io/biocontainers/crimson:1.1.0--pyh5e36f6f_0",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0",
    "mutalyzer": "docker://quay.io/biocontainers/mutalyzer_hgvs_parser:0.3.8--pyh7e72e81_0",
}

# If we run with the full hamlet configuration, subset the configuration
if "snv-indels" in config:
    config = config["snv-indels"]

# Old features that are no longer supported
if "bed_variant_call_regions" in config:
    msg = """
    'bed_variant_call_regions' is no longer supported, regions of interest are determined automatically from the GTF file.
    """
    raise DeprecationWarning(msg)

if "ref_id_mapping" in config:
    msg = """
    'ref_id_mapping' is no longer supported, the mapping is determined automatically from the GTF file.
    """
    raise DeprecationWarning(msg)

if "blacklist" in config:
    msg = """
    'blacklist' is no longer supported, please specify artifacts in the `known_variants` file.
    """
    raise DeprecationWarning(msg)

if "vep_include_consequence" in config:
    msg = """
    'vep_include_consequence' is no longer supported.
    """
    raise DeprecationWarning(msg)

# Put each sample name in a SimpleNamespace to mimic Snakemake wildcard usage
# (e.g {wildcards.sample}). This is only used in the 'all' rule.
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]

for s in samples:
    if " " in s.sample:
        raise RuntimeError(f'Spaces in samples are not supported ("{s.sample}")')


def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]


def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R2"]


def get_strand(wildcards):
    try:
        strandedness = pep.sample_table.loc[wildcards.sample, "strandedness"]
    except KeyError:
        strandedness = "unstranded"
    if strandedness is None:
        raise ValueError("Please specify a strandedness for every sample")
    return strandedness


def get_strand_picard(wildcards):
    strand = get_strand(wildcards)

    mapping = {
        "unstranded": "NONE",
        "forward": "FIRST_READ_TRANSCRIPTION_STRAND",
        "reverse": "SECOND_READ_TRANSCRIPTION_STRAND",
    }
    return mapping[strand]


## Functions for module outputs ##
def get_bam_output(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.bam"


def get_bai_output(wildcards):
    return get_bam_output(wildcards) + ".bai"


def get_filter_vep(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.vep.annotated.txt.gz"


def get_json(wildcards):
    return f"{wildcards.sample}/snv-indels/snv-indels-output.json"


def get_star_count(wildcards):
    return f"{wildcards.sample}/snv-indels/{wildcards.sample}.ReadsPerGene.out.tab"


def multiqc_files():
    star_count = [get_star_count(wildcards) for wildcards in samples]

    star_log = [f"{wildcards.sample}/snv-indels/Log.final.out" for wildcards in samples]

    picard_stats = list()
    for tool in ["rna_metrics", "alignment_summary_metrics", "insert_size_metrics"]:
        picard_stats += [
            f"{wildcards.sample}/snv-indels/{wildcards.sample}.{tool}"
            for wildcards in samples
        ]

    vep_stats = [
        f"{wildcards.sample}/snv-indels/{wildcards.sample}.vep_stats.txt"
        for wildcards in samples
    ]

    return star_count + star_log + picard_stats + vep_stats


def multiqc_modules():
    """Define which MultiQC modules to run here"""
    modules = ["star", "picard", "vep"]
    return [f" --module {module}" for module in modules]


module_output = SimpleNamespace(
    bam=get_bam_output,
    bai=get_bai_output,
    counts=get_star_count,
    filter_vep=get_filter_vep,
    id_mapping="id_mapping.txt",
    json=get_json,
    multiqc_files=multiqc_files(),
    multiqc_parquet="multiqc_snv_indels_data/multiqc.parquet",
)
