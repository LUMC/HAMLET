import pytest
import bs4
import base64


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
    headers = [header.get_text() for header in table.find_all("th")]

    # Get the rows
    for row in table.find("tbody").find_all("tr"):
        d = {k: v.get_text() for k, v in zip(headers, row.find_all("td"))}
        data.append(d)
    return data


def convert_image_to_base64(img_path):
    with open(img_path, "rb") as f:
        b64_bytes = base64.b64encode(f.read())
    return b64_bytes.decode("utf-8")


@pytest.mark.workflow("test-report")
def test_variant_overview(workflow_dir):
    """Test the content of the variant overview

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
    variant_table = soup.find("table", id="var-overview")

    # Extract the gene name columns
    gene1, gene2 = variant_table.findAll(attrs={"class": "hl"})
    assert gene1["rowspan"] == "3"
    assert gene1.text == "MT-ATP6"

    assert gene2["rowspan"] == "1"
    assert gene2.text == "MT-ATP8"


@pytest.mark.workflow("test-report")
def test_genes_in_order(workflow_dir):
    """Test if genes in table 2 are in order"""
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract table 2, genes of interest
    genes_of_interest = parse_table(soup.find("table", id="genes-of-interest"))
    assert "MT-ATP6" in genes_of_interest[0].values()
    assert "MT-ATP8" in genes_of_interest[1].values()


@pytest.mark.workflow("test-report")
def test_chr_location_exon(workflow_dir):
    """Test if the genomic HGVS and exons are in the table"""
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Get rows from the table, accounting for rowspan
    rows = get_rows(soup.find("table", id="var-overview"))

    # The genomic HGVS is in the second row
    next(rows)
    row = next(rows)

    # Convert into a dictionary
    header = "gene hgvs database VAF_exon annotation ref alt total".split()
    data = {k: v for k, v in zip(header, row)}

    # Test that the genomic HGVS is in the table
    assert "chrM:g.8701A>G" in data["hgvs"]
    # Test that the exon number is in the table in the VAF/Exon column
    assert data["VAF_exon"] == "95.83%1/1"
    assert data["annotation"] == "Hotspot"


@pytest.mark.workflow("test-report")
def test_fusion_overview(workflow_dir):
    """Test the content of the fusion overview"""
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the fusion table
    fusion_table = parse_table(soup.find("table", id="fusion-overview"))
    row = fusion_table[0]

    # Check the headers
    assert "Split reads 1" in row
    assert "Split reads 2" in row
    assert "Discordant mates" in row
    assert "Confidence" in row

    # Check that the breakpoints are in the table
    assert "(chr22, chr9)" in row["Fusion name"]


@pytest.mark.workflow("test-report")
def test_is_in_hotspot(workflow_dir):
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find("table", id="var-overview")
    expected_values = ["", "Hotspot", "Hotspot", "Hotspot"]
    hotspot_column = 4

    for row, expected in zip(get_rows(variant_table), expected_values):
        assert expected in row[hotspot_column]


@pytest.mark.workflow("test-report")
def test_database_identifiers(workflow_dir):
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find("table", id="var-overview")

    expected_values = ["rs3020563", "COSV104419767 rs2000975", "rs2001031", ""]

    database_column = 2

    for row, expected in zip(get_rows(variant_table), expected_values):
        assert row[database_column] == expected


@pytest.mark.workflow("test-full-report")
def test_full_variant_overview(workflow_dir):
    """Test the content of the variant overview

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
    variant_table = soup.find("table", id="var-overview")

    # Extract the gene name columns
    gene1, gene2 = variant_table.findAll(attrs={"class": "hl"})

    # gene1 (MT-ATP6) has 3 variants, but the first variant overlaps two
    # transcripts of interest. Therefore, it should span 4 rows
    assert gene1["rowspan"] == "4"
    assert gene1.text == "MT-ATP6"

    assert gene2["rowspan"] == "1"
    assert gene2.text == "MT-ATP8"


@pytest.mark.workflow("test-report-vardict")
def test_full_variant_overview_vardict(workflow_dir):
    """Test the content of the variant overview from vardict"""
    report = f"{workflow_dir}/report.html"
    with open(report) as fin:
        soup = bs4.BeautifulSoup(fin, features="html.parser")

    # Extract the variant table
    variant_table = soup.find("table", id="var-overview")

    # Check the first row
    row = parse_table(variant_table)[0]
    assert row["Ref/Alt(Total)"] == "0/27(27)"

    # The allele frequency should be given in percentage: 100%, not 1
    assert row["VAF/Exon"] == "100.0%1/1"


@pytest.mark.workflow("Test report expression genes")
def test_variant_overview_expression(workflow_dir):
    """Test the content of the variant overview

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
    expression_table = soup.find("table", id="gene-expression")
    table = parse_table(expression_table)

    # First row
    expected = {"Gene": "MT-ND4", "Raw count": "314", "Normalized expression": "1.11"}
    assert table[0] == expected

    # Second row
    expected = {"Gene": "MT-TH", "Raw count": "3", "Normalized expression": "0.01"}
    assert table[1] == expected


@pytest.mark.workflow("Test cell type and fusion images are embedded as base64 string")
def test_base64_image_in_report(workflow_dir):
    html_path = f"{workflow_dir}/report.html"
    imgs_path = [
        f"{workflow_dir}/SRR8615409/expression/seAMLess/cell-types.png",
        f"{workflow_dir}/SRR8615409/fusion/arriba/plots/fusion-1.png",
    ]

    with open(html_path, "r", encoding="utf-8") as f:
        html_content = f.read()

    for img_path in imgs_path:
        b64_str = convert_image_to_base64(img_path)
        assert (
            b64_str in html_content
        ), f"Base64 string for {img_path} not found in HTML report"
