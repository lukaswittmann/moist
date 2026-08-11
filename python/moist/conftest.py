"""Pytest configuration shared by the moist test modules.

Markers are registered here rather than in ``python/pyproject.toml`` because the
tests are normally run against the *built* or *installed* package -- e.g.
``pytest build/python/moist/test_pyscf.py`` -- where ``python/pyproject.toml``
is not pytest's rootdir and its ``[tool.pytest.ini_options]`` is never read.  A
conftest next to the tests travels with them and applies in every layout, so
the markers stay registered and ``PytestUnknownMarkWarning`` stays quiet.
"""

#: Marker name -> description. Each is also a meson test target; see
#: python/moist/meson.build.
MARKERS = {
    # cost
    "slow": "larger solutes and longer sweeps, kept out of the default run",
    # moist <-> PySCF coupling, by the layer a failure points at
    "host": "the PySCF host adapter's own derivatives, independent of moist",
    "vdw": "solute-vdW cavity: electrostatic conventions, no level-set response",
    "isodensity": "isodensity cavity: the level-set chain rule",
    "conventions": "channel semantics and negative controls",
    "scf": "self-consistent solvated SCF",
    # solvation components, crossed with `isodensity`
    "cpcm": "CPCM component only",
    "pv": "PV component only",
    "cpcm_pv": "CPCM and PV together",
    # GOSTSHYP pressure wall
    "gostshyp_integrals": "libcint fakemol constants and cartesian orders",
    "gostshyp_params": "surface-parameter derivatives at a frozen cavity",
    "gostshyp_fock": "energy and dE/dP on the moving surface",
    "gostshyp_gradient": "the nuclear-gradient channels",
    "gostshyp_conventions": "channel semantics and negative controls",
}


def pytest_configure(config):
    for name, description in MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {description}")
