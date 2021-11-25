#------------------------------------------------------------------------#
# Path setup                                                             #
#------------------------------------------------------------------------#
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))

#------------------------------------------------------------------------#
# Project info                                                           #
#------------------------------------------------------------------------#

project = 'OpenPNM'
copyright = '2021, PMEAL'
author = 'OpenPNM Dev Team'

# The full version, including alpha/beta/rc tags
from openpnm import __version__
release = __version__

#------------------------------------------------------------------------#
# General config                                                         #
#------------------------------------------------------------------------#

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # 'sphinx.ext.autodoc',
    # 'sphinx.ext.autosummary',
    # 'sphinx.ext.doctest',
    # 'sphinx.ext.duration',
    # 'sphinx.ext.mathjax',
    # 'sphinx.ext.napoleon',
    # 'sphinx.ext.viewcode',
    # 'sphinx_copybutton',
    # 'sphinx_panels',
    # 'nbsphinx',
    # 'nbsphinx_link',
    # 'numpydoc',
    # 'matplotlib.sphinxext.plot_directive'
]

autosummary_imported_members = True

panels_add_bootstrap_css = False  # to fix narrow width

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']


#------------------------------------------------------------------------#
# Matplotlib plot_directive options                                      #
#------------------------------------------------------------------------#

plot_pre_code = """
import numpy as np
np.random.seed(123)
"""
plot_include_source = False
plot_formats = [('svg', 200), 'pdf']
plot_html_show_formats = False
plot_html_show_source_link = False

import math
phi = (math.sqrt(5) + 1)/2

font_size = 13*72/96.0  # 13 px

plot_rcparams = {
    'font.size': font_size,
    'axes.titlesize': font_size,
    'axes.labelsize': font_size,
    'xtick.labelsize': font_size,
    'ytick.labelsize': font_size,
    'legend.fontsize': font_size,
    'figure.autolayout': True,
    'figure.figsize': (3*phi, 3),
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.85,
    'figure.subplot.wspace': 0.4,
    'text.usetex': False,
}

import matplotlib.pyplot as plt
plt.ioff()


#------------------------------------------------------------------------#
# HTML/theme options                                                     #
#------------------------------------------------------------------------#

html_theme = 'pydata_sphinx_theme!'

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
