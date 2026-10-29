"""``solvent-table`` directive: the solvent descriptor table, read from the Fortran data

Parses ``src/moist/data/solvents/descriptors.inc`` at build time, so the table
always matches the compiled library.
"""

import re
from pathlib import Path

from docutils import nodes
from docutils.statemachine import StringList
from sphinx.util.docutils import SphinxDirective

DESCRIPTORS = Path(__file__).parents[2] / "src" / "moist" / "data" / "solvents" / "descriptors.inc"
CAS = re.compile(r"\d{2,7}-\d{2}-\d")
FIELDS = ("eps", "refr", "A", "B", "g", "rho")


def _number(text, digits):
    """Fortran real literal rounded to ``digits`` decimals, trailing zeros removed"""
    value = round(float(text.replace("_wp", "")), digits)
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def read_solvents(path=DESCRIPTORS):
    """One dict per solvent block: id, name, descriptors, aliases, CAS number and notes"""
    text = path.read_text()
    solvents = []
    for block in re.split(r"\n   !-{10,}\n", text)[1:]:
        entry = {"notes": []}
        for line in block.splitlines():
            if line.startswith("   ! "):
                entry["notes"].append(line[5:].strip())
            elif match := re.match(r"\s*data (\w+)\(\d+\) / (.+) /", line):
                key, value = match.groups()
                entry[key] = value.strip('"')
        aliases = re.findall(r'"([^"]*)"', block[block.index("alias_list"):])
        aliases = [a for a in aliases if a]
        entry["cas"] = next((a for a in aliases if CAS.fullmatch(a)), "")
        entry["aliases"] = [a for a in aliases if a != entry["cas"] and a != entry["name_list"].lower()]
        solvents.append(entry)
    return solvents


class SolventTable(SphinxDirective):
    """List table of every solvent, with footnotes for deviating reference temperatures"""

    def run(self):
        self.env.note_dependency(str(DESCRIPTORS))
        solvents = read_solvents()
        rst = [
            ".. list-table::",
            "   :header-rows: 1",
            "",
            "   * - ID",
            "     - Name",
            "     - CAS",
            "     - :math:`\\varepsilon`",
            "     - :math:`n_\\mathrm{D}`",
            "     - :math:`\\alpha`",
            "     - :math:`\\beta`",
            "     - :math:`\\gamma` / mN m\\ :sup:`-1`",
            "     - :math:`\\rho` / kg m\\ :sup:`-3`",
            "     - Aliases",
        ]
        notes = []
        for s in solvents:
            mark = ""
            if s["notes"]:
                label = f"solvent-{s['id_list']}"
                notes.append(f".. [#{label}] Reference temperatures: {' '.join(s['notes'])}")
                mark = f" [#{label}]_"
            rst += [
                f"   * - {s['id_list']}",
                f"     - {s['name_list'].lower()}{mark}",
                f"     - {s['cas']}",
                f"     - {_number(s['eps'], 4)}",
                f"     - {_number(s['refr'], 4)}",
                f"     - {_number(s['A'], 3)}",
                f"     - {_number(s['B'], 3)}",
                f"     - {_number(s['g'], 2)}",
                f"     - {_number(s['rho'], 1)}",
                f"     - {', '.join(s['aliases'])}",
            ]
        if notes:
            rst += [""] + notes
        container = nodes.container()
        self.state.nested_parse(StringList(rst, source=str(DESCRIPTORS)), self.content_offset, container)
        return container.children


def setup(app):
    app.add_directive("solvent-table", SolventTable)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
