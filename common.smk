from functools import partial
from os.path import dirname
from uuid import uuid4

pepfile: config["pepfile"]

samples = pep.sample_table["sample_name"]

containers = {
    "debian": "docker://debian:buster-slim",
    "crimson": "docker://quay.io/biocontainers/crimson:1.0.0--pyh5e36f6f_0",
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.3",
    "zip": "docker://quay.io/redmar_van_den_berg/zip:3.0"
}

settings=config["settings"]

# The version of HAMLET
PIPELINE_VERSION = "v1.0.2-dev-1"
RUN_NAME = settings.get("run_name") or f"hamlet-{uuid4().hex[:8]}"

def make_pattern(extension, dirname):
    """Helper function to create a wildcard-containing path for output files."""
    return f"{{sample}}/{dirname}/{{sample}}{extension}"

seqqc_output = partial(make_pattern, dirname="qc-seq")
var_output = partial(make_pattern, dirname="snv-indels")
fusion_output = partial(make_pattern, dirname="fusion")
expr_output = partial(make_pattern, dirname="expression")
itd_output = partial(make_pattern, dirname="itd")
