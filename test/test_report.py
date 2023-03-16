import pytest
import bs4

def get_rows(table):
    """Get rows from table, supports rowspan in the first column only"""
    # To deal with rowspan
    rowspan = None
    gene = None

    rows = table.find("tbody").find_all("tr")
    for row in rows:
        cells = row.find_all("td")
        text = [x.get_text() for x in cells]

        # We need to do some special magic for rowspan in the first column
        if cells[0].get("rowspan"):
            gene = cells[0].get_text()
            rowspan = int(cells[0]["rowspan"])
        elif rowspan > 0:
            # Prepend the gene name
            text.insert(0, gene)
            rowspan -= 1

        yield text

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


@pytest.mark.workflow('test-report')
def test_is_in_hotspot(workflow_dir):
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find('table', id='var-overview')
    expected_values = ['no', 'yes', 'yes', 'yes']
    hotspot_column = 4

    for row, expected in zip(get_rows(variant_table), expected_values):
        assert row[hotspot_column] == expected


@pytest.mark.workflow('test-full-report')
def test_full_variant_overview(workflow_dir):
    """ Test the content of the variant overview

    The structure of the data that gets put in the variant overview is quite
    complex.
    - Every gene can have zero or more variants
    - Every variant can overlap zero or more transcripts of interest

    This test uses a modified summary.json file where the first variant is made
    to influence two transcripts of interest (by simply duplicating one variant
    of interest from the original file).
    """
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find('table', id='var-overview')

    # Extract the gene name columns
    gene1, gene2 = variant_table.findAll(attrs={"class": "hl"})

    # gene1 (MT-ATP6) has 3 variants, but the first variant overlaps two
    # transcripts of interest. Therefore, it should span 4 rows
    assert gene1["rowspan"] == "4"
    assert gene1.text == "MT-ATP6"

    assert gene2["rowspan"] == "1"
    assert gene2.text == "MT-ATP8"
