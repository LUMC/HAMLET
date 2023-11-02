from uuid import uuid4


pepfile: config["pepfile"]


samples = pep.sample_table["sample_name"]

containers = {
    "hamlet-scripts": "docker://quay.io/redmar_van_den_berg/hamlet-scripts:0.3",
}

# The version of HAMLET
PIPELINE_VERSION = "v2.0.1"
