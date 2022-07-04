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
import shutil
from pathlib import Path
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'MRI_DistortionQA'
copyright = '2021, Brendan Whelan(s)'
author = 'Brendan Whelan(s)'



# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    # 'recommonmark',
    'sphinx_markdown_tables',
    'sphinxcontrib.mermaid',
    'myst_parser',
    'sphinxemoji']

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

autodoc_mock_imports = ["utilities"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# copy report demo -----------------------------------------------------------
def copy_and_overwrite(from_path, to_path):
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree(from_path, to_path)



this_dir = Path(__file__).parent
print(f'{this_dir} exists? {os.path.isdir(this_dir)}')
print(f'this_dir / _static? {os.path.isdir(this_dir / "_static")}')
print(f'MR_QA_report_20_05_2022.html exists? {os.path.isfile(this_dir / "_static" / "MR_QA_report_20_05_2022.html")}')
print(f'this_dir / _build? {os.path.isdir(this_dir / "_build")}')
print(f'this_dir / docs? {os.path.isdir(this_dir / "docs")}')
print(f'this_dir / .. / docs? {os.path.isdir(this_dir.parent / "docs")}')


shutil.copy(this_dir / '_static' / 'MR_QA_report_20_05_2022.html',
            this_dir / '_build' / 'html' / '_static' / 'MR_QA_report_20_05_2022.html')
copy_and_overwrite(this_dir / '_static' / 'plots',
            this_dir / '_build' / 'html' / '_static' / 'plots')
copy_and_overwrite(this_dir / '_static' / 'themes',
            this_dir / '_build' / 'html' / '_static' / 'themes')

