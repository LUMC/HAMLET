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


def list_files(folder) -> list[str]:
    """Recursively list all files in folder"""
    files = list()

    for item in os.listdir(folder):
        path = os.path.join(folder, item)
        if os.path.isfile(path):
            files.append(path)
        elif os.path.isdir(path):
            files += list_files(path)
        else:
            msg = f"Unknown path: {path=}"
            raise RuntimeError(msg)
    return files


def report_files(wildcards):
    return list_files("report")


# The version of HAMLET
PIPELINE_VERSION = "v2.3.1-dev"
