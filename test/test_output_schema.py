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

@pytest.mark.workflow('test-qc-trio')
def test_qc_seq_schema(workflow_dir):
    sample = "TestSample3"
    output_file = pathlib.Path(workflow_dir, f"{sample}/qc-seq/{sample}.seq_stats.json")
    schema_file = pathlib.Path(workflow_dir, "includes/qc-seq/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-expression-chrM')
def test_expression_schema(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/expression/{sample}.exon_ratios.json")
    schema_file = pathlib.Path(workflow_dir, "includes/expression/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-fusion-chrM')
def test_fusion_schema(workflow_dir):
    sample = "SRR8615409"
    output_file = pathlib.Path(workflow_dir, f"{sample}/fusion/fusion-output.json")
    schema_file = pathlib.Path(workflow_dir, "includes/fusion/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-fusion-output')
def test_example_fusion_output_schema(workflow_dir):
    output_file = pathlib.Path(workflow_dir, "fusion-output.json")
    schema_file = pathlib.Path(workflow_dir, "includes/fusion/output-schema.json")
    validate_files(output_file, schema_file)

@pytest.mark.workflow('test-fusion-output')
def test_example_fusion_output_content(workflow_dir):
    output_file = pathlib.Path(workflow_dir, "fusion-output.json")
    with open(output_file) as fin:
        js = json.load(fin)

    # Test if we recognised the data is intersected
    assert js["fusion"]["intersected"]

    # Test if we have a table for each tool
    tools = ["fusioncatcher", "star-fusion", "intersection"]
    for tool in tools:
        assert tool in js["fusion"]["tables"]

    # Test first and last result for intersection
    results = js["fusion"]["tables"]["intersection"]["top20"]
    first = results[0]
    last = results[-1]

    assert first["name"] == "PLAA--MIR31HG"
    assert first["sf_count"] == 101

    assert last["type"] == "INCL_NON_REF_SPLICE"
    assert last["jr_count"] == 3

    # Test first result for star-fusion
    results = js["fusion"]["tables"]["intersection"]["top20"]
    first = results[0]
    last = results[-1]

    assert first["name"] == "PLAA--MIR31HG"
    assert first["sf_count"] == 101

@pytest.mark.workflow('test-fusion-output')
def test_example_fusion_output_pngs(workflow_dir):
    """ Test if we have a png plot for each tool """
    output_file = pathlib.Path(workflow_dir, "fusion-output.json")
    with open(output_file) as fin:
        js = json.load(fin)

    # Get the plots for each tool
    fc = js["fusion"]["plots"]["fusioncatcher"]
    isect = js["fusion"]["plots"]["intersection"]
    star = js["fusion"]["plots"]["star-fusion"]

    # Check the content
    assert fc.endswith("sample.fusioncatcher-circos/fsnviz.png")
    assert isect.endswith("sample.sf-isect-circos/fsnviz.png")
    assert star.endswith("sample.star-fusion-circos/fsnviz.png")
