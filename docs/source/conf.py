# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Sam54000-eeg_modules'
copyright = '2023, Samuel Louviot'
author = 'Samuel Louviot'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'mne': ('https://mne.tools/stable/', None),
    'mne_bids': ('https://mne.tools/mne-bids/stable/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

napoleon_google_docstring = True
napoleon_numpy_docstring = True

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
