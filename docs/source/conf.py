# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import coderdata 

project = 'CoderData'
copyright = '2025, Sara Gosline'
author = 'Sara Gosline'
release = '1.0.0'

import os
import sys
sys.path.insert(0, os.path.abspath('../../coderdata')) #need to add /coderdata to get build functions


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", # used to pull documentation
              "sphinx.ext.coverage",
              "sphinx.ext.extlinks",
              "sphinx.ext.napoleon",
              "nbsphinx", # conversion with jupyter notebook
              "myst_parser", # allows us to include md files
              "sphinx.ext.doctest", # reates test documentation
              "sphinx_tabs.tabs", # allow for tabs in web struct THIS IS NOT WHAT I THINK IT SMH. BUT COULD BE USED FOR WINDOWS VS MAC EXAMPLES 
              "sphinx.ext.autosectionlabel",
 #             "sphinx_design",
 #             "sphinx.ext.apidoc" # to document a whole package
] 


# this is also to allow for md files 
source_suffix = ['.rst','.md']

autosummary_generate = True
autodoc_typehints = ["none"]

# myst_enable_extensions = ["colon_fence", "linkify", "substitution"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

master_doc = 'index'
templates_path = ['_templates']
html_static_path = ['_static']

html_theme = 'sphinxdoc'
html_theme_path = ['_static'] 
html_css_files = ["custom.css"] #modify in this location: docs/source/_static
                                #this is how to override the builtin sphinx theme 'sphinxdoc'
                                

html_title = 'CoderData'
html_logo = "_static/coderdata_logo3.png"

html_show_sourcelink = True #hide source for the html 

html_sidebars = {
    '**': ['my_custom_sidebar.html', 'localtoc.html', 'searchbox.html'],
}
# my_custom_sidebar = Useful links 
# localtoc = Table of Contents
# searchbox = Quick search


#ignore header warnings and non-referenced documents that does not interefere with the build
suppress_warnings = ["myst.header", "myst.reference", "toc.not_readable"] 
nitpicky= False
# exclude_pattterns =['library/xml' – ignores the library/xml directory]


 








