
import sys, os.path as op

sys.path.insert(0, op.join(op.dirname(__file__), "..", "python"))

try:
    import moist  # noqa: F401
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
]

# Global bibliography (sphinxcontrib-bibtex). Cite from any page with the
# :cite:t: / :cite:p: roles; the entries are defined in references.bib and the
# collected list is rendered once on the References page.
bibtex_bibfiles = ["_static/references.bib"]

html_theme = "sphinx_book_theme"
html_title = project
html_logo = "_static/moist.svg"
html_favicon = html_logo

html_theme_options = {
    "repository_url": "https://github.com/lukaswittmann/moist",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "path_to_docs": "doc",
    "show_navbar_depth": 3
}

html_sidebars = {}

html_css_files = [
    "css/custom.css",
]
html_static_path = ["_static"]
templates_path = ["_templates"]

autodoc_typehints = "none"
autodoc_mock_imports = ["moist.library", "numpy", "ase", "pyscf"]
autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

master_doc = "index"
