# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'MolR'
copyright = '2025, Abhishek Tiwari'
author = 'Abhishek Tiwari'

# The full version, including alpha/beta/rc tags
import molr
release = molr.__version__
version = molr.__version__

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_logo = '../../molr-logo-animated.svg'

# Custom CSS files
html_css_files = [
    'custom.css',
]
# sphinx_book_theme options
# HTML theme options
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'google_analytics_id': 'G-PLZNCJY1B7',
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#343131',  # Dark header background
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    "repository_url": "https://github.com/abhishektiwari/molr",
    "repository_provider": "github",
    "repository_branch": "main",
    "path_to_docs": "docs/source",
    "use_issues_button": True,
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "use_fullscreen_button": True,
    "use_search_button": True,
    "use_sidenotes": True,
    "icon_links_label": "Quick Links",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/abhishektiwari/molr",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
        {
            "name": "Abhishek Tiwari",
            "url": "https://www.abhishek-tiwari.com",
            "icon": "https://www.abhishek-tiwari.com/images/logo.svg",
            "type": "local",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/molr/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
   ]
}

# -- Extension configuration -------------------------------------------------

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Type hints
always_document_param_types = True
typehints_defaults = 'comma'

# Handle forward references
autodoc_type_aliases = {
    'BondList': 'molr.BondList',
}