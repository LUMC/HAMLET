# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'HAMLET'
copyright = '2018, LUMC'
author = 'Wibowo Arindrarto, Redmar van den Berg'

release = '2.0'
version = '2.0.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
master_doc = 'index'

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
