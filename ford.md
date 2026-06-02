project: MOIST
summary: >
  Modular and Open-source Implicit Solvation Toolkit providing reusable Fortran
  components, CLI tools, and bindings for implicit solvation models.
project_url: https://github.com/lukaswittmann/moist
github: lukaswittmann/moist
license: LGPL-3.0-or-later
authors:
  - Lukas Wittmann
output_dir: doc/ford
src_dir:
  - src
  - app
  - test/unit
docmark: "!>"
display: public
graph: true
graphviz_dot: dot
project_download: https://github.com/lukaswittmann/moist/releases
source: true
sort: alpha
search: true
exclude_dir:
  - assets
  - build
  - config
  - doc
  - include
  - man
  - python
  - subprojects
  - test/api
  - tools

---

# MOIST developer documentation

FORD consumes this file to produce browsable API pages for the Fortran
components. The defaults above expose the source under `src/`, the CLI entry
points in `app/`, and the Fortran unit tests in `test/unit/`. The generated HTML
is written to `doc/ford`.

## Quick start

```sh
pip install ford # or use your package manager
sudo apt-get install graphviz # or equivalent, needed for graphs
ford ford.md
```

Open `doc/ford/index.html` in your browser to view the rendered documentation.

## Customizing

- Add Markdown pages (e.g., design notes) next to this file and list them under
  `extras` to make them appear in the navigation tree.
- Adjust `src_dir` if you introduce additional Fortran trees (for example, a new
  module under `python/` that needs documentation).
- Update the author list, summary, or license fields to match release metadata
  when the project description evolves.

See the FORD manual for more configuration options:
https://github.com/Fortran-FOSS-Programmers/ford.
