# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------

project = 'OpenPNM'
copyright = '2021, PMEAL'
author = 'OpenPNM Dev Team'

# The full version, including alpha/beta/rc tags
from openpnm import __version__
release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.duration',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinx_panels',
    'nbsphinx',
    'nbsphinx_link',
    'numpydoc',
    'matplotlib.sphinxext.plot_directive'
]

autosummary_imported_members = True

panels_add_bootstrap_css = False  # to fix narrow width

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

html_js_files = ['js/custom.js']
html_css_files = ['css/custom.css']

html_logo = '_static/images/openpnm_logo.png'

html_theme_options = {
    "logo_link": "https://www.openpnm.org",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/PMEAL/OpenPNM",
            "icon": "fab fa-github-square",
        },
        {
            "name": "Substack",
            "url": "https://openpnm.substack.com/",
            "icon": "fas fa-envelope-square",
        },
        {
            "name": "Twitter",
            "url": "https://twitter.com/OpenPNM",
            "icon": "fab fa-twitter-square",
        },
    ],
    "external_links": [
        {
            "name": "Issue Tracker", "url": "https://github.com/PMEAL/OpenPNM/issues"
        },
        {
            "name": "Get Help", "url": "https://github.com/PMEAL/OpenPNM/discussions"
        },
    ],
    "navigation_with_keys": False,
    "show_prev_next": False,
    "icon_links_label": "Quick Links",
    "use_edit_page_button": False,
    "search_bar_position": "sidebar",
    "navbar_align": "left",
}

html_sidebars = {
    "contributing": ["sidebar-search-bs.html"],
    "changelog": [],
    "examples/*": []
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
