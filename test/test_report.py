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
        elif rowspan and rowspan > 0:
            # Prepend the gene name
            text.insert(0, gene)
            rowspan -= 1

        yield text

def parse_table(table):
    """Convert a table to a list of dicts"""
    data = list()
    # Get the headers
    headers = [header.get_text() for header in table.find_all('th')]

    # Get the rows
    for row in table.find("tbody").find_all("tr"):
        d = {k: v.get_text() for k, v in zip(headers, row.find_all("td"))}
        # Convert the price field to float
        for k, v in d.items():
            if "€" in v:
                d[k] = price_to_float(v)
        data.append(d)
    return data

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
def test_fusion_overview(workflow_dir):
    """ Test the content of the fusion overview

    """
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the fusion table
    fusion_table = parse_table(soup.find('table', id='fusion-overview'))

    # Check the headers
    assert "Split reads 1" in fusion_table[0]
    assert "Split reads 2" in fusion_table[0]
    assert "Discordant mates" in fusion_table[0]
    assert "Confidence" in fusion_table[0]


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


@pytest.mark.workflow('test-report')
def test_database_identifiers(workflow_dir):
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find('table', id='var-overview')

    expected_values = [
            "rs3020563",
            "COSV104419767, rs2000975",
            "rs2001031",
            ""
    ]

    database_column = 2

    for row, expected in zip(get_rows(variant_table), expected_values):
        assert row[database_column] == expected


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

@pytest.mark.workflow('test-report-vardict')
def test_full_variant_overview_vardict(workflow_dir):
    """ Test the content of the variant overview from vardict
    """
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find('table', id='var-overview')

    # Check the first row
    row = parse_table(variant_table)[0]
    assert row["Ref"] == '0'
    assert row["Alt"] == '27'
    assert row["Total"] == '27'

    # The allele frequency should be given in percentage: 100%, not 1
    assert row["Allele frequency"] == '100.0%'
