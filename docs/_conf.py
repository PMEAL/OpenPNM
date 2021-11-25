import os
import sys
import datetime

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))

#------------------------------------------------------------------------------
# General configuration
#------------------------------------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.autosummary',
              'sphinx.ext.ifconfig',
              'sphinx.ext.viewcode',
              'sphinx.ext.mathjax',
              'sphinx_copybutton',
              'nbsphinx',
              'nbsphinx_link',
              #'numpydoc',
              'sphinx_panels',
              'matplotlib.sphinxext.plot_directive']

#------------------------------------------------------------------------------
# Matplotlib plot_directive options
#------------------------------------------------------------------------------

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

# So that 'sphinx-copybutton' only copies the actual code, not the prompt
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# html_js_files = ['js/custom.js']

# html_css_files = ['css/custom.css']

nbsphinx_prompt_width = "0"

panels_add_bootstrap_css = False  # to fix narrow width

exclude_patterns = ['_build']

add_module_names = False

autosummary_generate = True

autosummary_imported_members = True

globaltoc_maxdepth = 2

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# General information about the project.

project = 'OpenPNM'
year = datetime.datetime.now().year
copyright = '%d OpenPNM Team' % year
author = 'PMEAL Team'

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

html_theme = 'furo'

html_logo = '_static/images/openpnm_logo.png'

html_static_path = ['_static']
