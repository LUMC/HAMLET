# Configuration file for the Sphinx documentation builder.

# -- Project information

project = "HAMLET"
copyright = "2018, LUMC"

author = " Redmar van den Berg, Wibowo Arindrarto, Xiaoyun Liu, Anne van der Grinten, Ramin Shirali"

# Fudge the latex author list, insert \and where you want the line break to appear
latex_elements = {
    "maketitle": r"""\author{Redmar van den Berg, Wibowo Arindrarto, Xiaoyun Liu, \and Anne van der Grinten, Ramin Shirali}
\sphinxmaketitle
"""
}

release = "v2.5.2-dev"
version = ".".join(release.split(".")[0:2])

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
master_doc = "index"

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"

linkcheck_retries = 3
linkecheck_timeout = 60
