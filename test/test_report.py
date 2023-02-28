import pytest
import bs4

@pytest.mark.workflow('test-report')
def test_variant_overview(workflow_dir):
    """ Test the content of the variant overview

    The structure of the data that gets put in the variant overview is quite
    complex.
    - Every gene can have zero or more variants
    - Every variant can overlap zero or more transcripts of interest

    This test ensures that the the gene name spans the apropriate number of
    columns in the variant table.
    """
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find('table', id='var-overview')

    # Extract the gene name columns
    gene1, gene2 = variant_table.findAll(attrs={"class": "hl"})
    assert gene1["rowspan"] == "3"
    assert gene1.text == "MT-ATP6"

    assert gene2["rowspan"] == "1"
    assert gene2.text == "MT-ATP8"
