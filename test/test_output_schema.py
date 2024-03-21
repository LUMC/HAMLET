import json
import jsonschema
import pathlib
import pytest

def validate_files(output_file, schema_file):
    # Load the output file
    with open(output_file) as fin:
        hamlet_output = json.load(fin)

    # Load the schema
    with open(schema_file) as fin:
        schema = json.load(fin)

    jsonschema.validate(instance=hamlet_output, schema=schema)

@pytest.mark.workflow('test-hamlet-chrM')
def test_output_against_schema(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/{sample}.summary.json")
    schema_file = pathlib.Path(workflow_dir, "utilities/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-hamlet-targetted-RNA')
def test_output_targetted_seq_against_schema(workflow_dir):
    sample = "MO1-RNAseq-1-16714"
    output_file = pathlib.Path(workflow_dir, f"{sample}/{sample}.summary.json")
    schema_file = pathlib.Path(workflow_dir, "utilities/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-fusion-chrM')
def test_fusion_schema(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/fusion/fusion-output.json")
    schema_file = pathlib.Path(workflow_dir, "includes/fusion/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-snv-indels-chrM')
def test_snv_indel_schema(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/snv-indels/snv-indels-output.json")
    schema_file = pathlib.Path(workflow_dir, "includes/snv-indels/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-itd')
def test_itd_schema(workflow_dir):
    sample = "SRR8616218"
    output_file = pathlib.Path(workflow_dir, f"{sample}/itd//itd-output.json")
    schema_file = pathlib.Path(workflow_dir, "includes/itd/output-schema.json")
    validate_files(output_file, schema_file)


if __name__ == '__main__':
    import sys
    instance = sys.argv[1]
    schema = sys.argv[2]

    validate_files(instance, schema)
