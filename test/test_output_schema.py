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
