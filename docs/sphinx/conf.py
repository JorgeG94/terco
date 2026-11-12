# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "terco"
copyright = "2026, Jorge Luis Galvez Vallejo"
author = "Jorge Luis Galvez Vallejo"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx_fortran_domain",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# MyST: let a fenced ```{directive}``` block span a colon-fence, which is how
# the f: domain's directives (f:module, f:function, ...) are written in the
# Markdown pages below.
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- sphinx-fortran-domain ---------------------------------------------------
# https://github.com/fortran-lang/sphinx-fortran-domain
#
# terco documents with FORD's own convention (AGENTS.md): a `!!` block goes
# AFTER the entity it documents -- after a declaration, or as the first lines
# of a procedure's body, before any declaration. That is FORD's native
# docmark, so the "ford" lexer (backed by an actual FORD parse, not a regex
# approximation) is what attaches these correctly; "regex" does not know the
# after-entity convention.
#
# Only the public API is parsed here, never src/generated/ (generator
# output, AGENTS.md: "never edit by hand", and not what a consumer reads
# against). trc_api.F90 is the Fortran facade; trc_capi.F90 is the `bind(c)`
# ABI and carries almost all of the FORD `!!` doc comments in this API, since
# trc_api.F90 itself only re-exports names from lower-layer modules (see
# api.md for what that means for the rendered page).
fortran_sources = [
    "../../src/api/trc_api.F90",
    "../../src/api/trc_capi.F90",
]
fortran_lexer = "ford"
# Single character after "!": terco's doc marker is "!!", never "!>".
fortran_doc_chars = ["!"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
