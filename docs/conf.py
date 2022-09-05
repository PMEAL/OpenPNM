#------------------------------------------------------------------------#
# Path setup                                                             #
#------------------------------------------------------------------------#
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))

#------------------------------------------------------------------------#
# Project info                                                           #
#------------------------------------------------------------------------#

project = 'OpenPNM'
copyright = f'{datetime.now().year}, PMEAL'
author = 'OpenPNM Dev Team'

# The full version, including alpha/beta/rc tags
from openpnm import __version__
release = __version__

# Copy examples folder from OpenPNM root to docs folder
import shutil
shutil.copytree('../examples', '_examples', dirs_exist_ok=True)

#------------------------------------------------------------------------#
# General config                                                         #
#------------------------------------------------------------------------#

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.duration',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinx_design',
    'sphinx.ext.ifconfig',
    'matplotlib.sphinxext.plot_directive',
    'myst_nb',
    'numpydoc',
]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]

autosummary_imported_members = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']

numpydoc_show_inherited_class_members = True

members_to_exclude = [
    # NumPy methods we don't want to document
    'append', 'clear', 'copy', 'count', 'extend', 'fromkeys', 'get',
    'index', 'insert', 'items', 'keys', 'pop', 'popitem', 'remove',
    'reverse', 'setdefault', 'sort', 'update', 'values',
    'all', 'any', 'argmax', 'argmin', 'argpartition', 'argsort',
    'astype', 'base', 'byteswap', 'choose', 'clip', 'compress', 'conj',
    'conjugate', 'copy', 'ctypes', 'cumprod', 'cumsum', 'data',
    'diagonal', 'dot', 'dtype', 'dump', 'dumps', 'fill', 'flags', 'flat',
    'flatten', 'getfield', 'imag', 'item', 'itemset', 'itemsize', 'max',
    'mean', 'min', 'nbytes', 'ndim', 'newbyteorder', 'nonzero',
    'partition', 'prod', 'ptp', 'put', 'ravel', 'real', 'repeat',
    'reshape', 'resize', 'round', 'searchsorted', 'setfield', 'setflags',
    'shape', 'size', 'sort', 'squeeze', 'std', 'strides', 'sum',
    'swapaxes', 'take', 'tobytes', 'tofile', 'tolist', 'tostring',
    'trace', 'transpose', 'var', 'view', '__call__',
]

autodoc_default_options = {
    'exclude-members': ", ".join(members_to_exclude)
}

#------------------------------------------------------------------------#
# Matplotlib plot_directive options                                      #
#------------------------------------------------------------------------#

plot_pre_code = """
import numpy as np
np.random.seed(123)
"""
plot_include_source = True
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

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_js_files = ['js/custom.js']
html_css_files = ['css/custom.css']
html_logo = '_static/images/openpnm_logo.jpg'

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/PMEAL/OpenPNM",
            "icon": "fab fa-github-square",
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
    "navbar_align": "left",
}

html_sidebars = {
    # "examples_index": [],
}

html_js_files = ['js/custom.js']

nbsphinx_execute = 'always'
nbsphinx_prompt_width = "0"
nbsphinx_allow_errors = True

exclude_patterns = ['_build', '_templates']

add_module_names = False

autosummary_generate = True

globaltoc_maxdepth = 3

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ['openpnm']

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False
