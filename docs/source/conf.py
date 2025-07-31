# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CoderData'
copyright = '2025, Sara Gosline'
author = 'Sara Gosline'
release = '1.0.0'

import os
import sys
sys.path.insert(0, os.path.abspath('../../coderdata')) 

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", # used to pull documentation
              "sphinx.ext.coverage",
              "sphinx.ext.extlinks",
              "sphinx.ext.napoleon",
              "sphinx.ext.doctest", # reads test documentation
              "sphinx.ext.autosectionlabel",
              # Above extensions are all built-in Sphinx
              "nbsphinx", # conversion with jupyter notebook 
                            #install nbsphinx
              "myst_parser", # allows us to include md files
                            #install myst-parser
              "sphinx_tabs.tabs", # allow for tabs in web struct  
                            #install sphinx-tabs
              # "sphinx.ext.apidoc" # to document a whole package
                            ##install sphinx-design
] 


exclude_patterns= ['_notebooks']
nbsphinx_allow_errors= True

# this allows for md files 
source_suffix = {'.rst': 'restructuredtext',
                 '.md': 'markdown',
}

# docstrings changes
autosummary_generate = True
autodoc_typehints = ["none"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

master_doc = 'index'
templates_path = ['_templates']
html_static_path = ['_static']

html_theme = 'sphinxdoc'
html_theme_path = ['_static'] 
html_css_files = ["custom.css"] #modify in this location: docs/source/_static
                                #this is how to override the builtin sphinx theme 'sphinxdoc'
                                
html_baseurl = 'https://pnnl-compbio.github.io/coderdata/'

html_title = 'CoderData'
html_logo = "_static/coderdata_logo3.jpg"

html_show_sourcelink = True #hide source for the html 

html_sidebars = {
    '**': ['my_custom_sidebar.html', 'localtoc.html', 'searchbox.html'],
}
# my_custom_sidebar = Useful links 
# localtoc = Table of Contents
# searchbox = Quick search


#ignore header warnings and non-referenced documents that does not interfer with the build process.
suppress_warnings = ["myst.header", "myst.reference", "toc.not_readable","autosectionlabel","docutils"] 
nitpicky= False


# -----------process to convert notebook files to html-----------------

import os
import subprocess

# Directory where the Jupyter notebooks are located
notebook_dir = os.path.abspath('./_notebooks')  # Adjust path as needed

# Directory to store the converted HTML files
output_dir = os.path.abspath('../build/html')  # Adjust path as needed

def convert_notebooks_to_html():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for notebook in os.listdir(notebook_dir):
        if notebook.endswith('.ipynb'):
            input_path = os.path.join(notebook_dir, notebook)
            output_path = os.path.join(output_dir, notebook.replace('.ipynb', '.html'))
            cmd = [
                'jupyter', 'nbconvert',
                '--to', 'html',
                '--output', output_path,
                input_path
            ]
            subprocess.run(cmd, check=True)

# Run the conversion before Sphinx builds the documentation
convert_notebooks_to_html()



 








