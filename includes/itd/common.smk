containers = {
    "bwa-0.7.17-samtools-1.3.1-picard-2.9.2": "docker://quay.io/biocontainers/mulled-v2-1c6be8ad49e4dfe8ab70558e8fb200d7b2fd7509:5900b4e68c4051137fffd99165b00e98f810acae-0",
    "rose": "docker://quay.io/redmar_van_den_berg/rose-dt:0.4"
}

settings = config ["settings"]

pepfile: config["pepfile"]

samples = pep.sample_table["sample_name"]


# Set the default settings
def set_default(key, value):
    """ Set default value for settings """
    if key not in settings:
        settings[key] = value

# Region of interest for kmt2a ~ exon 2-10 in transcript coordinates.
set_default("kmt2a_name", "KMT2A-213")
set_default("kmt2a_start", 406)
set_default("kmt2a_end", 4769)

# Region of interest for flt3 ~ exon 14-15 in transcript coordinates.
set_default("flt3_name", "FLT3-001")
set_default("flt3_start", 1787)
set_default("flt3_end", 2024)

def get_forward(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]

def get_reverse(wildcards):
    return pep.sample_table.loc[wildcards.sample, "R1"]
