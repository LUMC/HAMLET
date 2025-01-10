from uuid import uuid4


pepfile: config["pepfile"]


samples = pep.sample_table["sample_name"]

for s in samples:
    if " " in s:
        raise RuntimeError(f'Spaces in samples are not supported ("{s}")')


containers = {
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.3",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.22.1--pyhdfd78af_0",
}

# The version of HAMLET
PIPELINE_VERSION = "v2.2.2-dev"
