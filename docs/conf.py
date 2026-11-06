
import sys, os.path as op

sys.path.insert(0, op.join(op.dirname(__file__), "..", "python"))
sys.path.insert(0, op.join(op.dirname(__file__), "_ext"))

try:
    import moist
except ImportError:
    pass


project = "moist"
author = "Lukas Wittmann"
copyright = f"2026 {author}"

extensions = [
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "solvents",
]

# Status icons shared by the support tables
rst_prolog = """
.. |yes| replace:: :octicon:`check-circle-fill;1em;sd-text-success`
.. |light| replace:: :octicon:`check-circle-fill;1em;sd-text-success tick-light`
.. |part| replace:: :octicon:`check-circle-fill;1em;sd-text-warning`
.. |no| replace:: :octicon:`x-circle;1em;sd-text-muted`
.. |soon| replace:: :bdg-info:`soon`
"""

# Global bibliography (sphinxcontrib-bibtex)
bibtex_bibfiles = ["_static/references.bib"]

html_theme = "sphinx_book_theme"
html_title = project
html_logo = "_static/logos/moist-light.svg"
html_favicon = html_logo

html_theme_options = {
    "logo": {
        "image_light": "_static/logos/moist-light.svg",
        "image_dark": "_static/logos/moist-dark.svg",
        "alt_text": "MOIST",
    },
    "repository_url": "https://github.com/lukaswittmann/moist",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "path_to_docs": "docs",
    "show_navbar_depth": 3
}

html_sidebars = {}

html_static_path = ["_static"]
html_css_files = ["custom.css"]
templates_path = ["_templates"]
exclude_patterns = ["_build"]

autodoc_typehints = "none"
autodoc_mock_imports = ["moist.library", "numpy", "ase", "pyscf"]
autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

master_doc = "index"
