"""Pytest configuration shared by the moist test modules
"""

#: Marker name -> description. Each is also a meson test target; see
#: python/moist/meson.build.
MARKERS = {
    # cost
    "slow": "larger solutes and longer sweeps, kept out of the default run",
    # moist <-> PySCF (QM) coupling
    "host": "the PySCF host adapter's own derivatives, independent of moist",
    "vdw": "solute-vdW cavity: electrostatic conventions, no level-set response",
    "isodensity": "isodensity cavity: the level-set chain rule",
    "conventions": "channel semantics and negative controls",
    "scf": "self-consistent solvated SCF",
    # Solvation model components with isodensity caivty
    "cpcm": "CPCM component only",
    "pv": "PV component only",
    "cpcm_pv": "CPCM and PV together",
    # GOSTSHYP pressure model
    "gostshyp": "umbrella marker applied to every test in test_gostshyp.py",
    "gostshyp_integrals": "libcint fakemol constants and cartesian orders",
    "gostshyp_params": "surface-parameter derivatives at a frozen cavity",
    "gostshyp_fock": "energy and dE/dP on the moving surface",
    "gostshyp_gradient": "the nuclear-gradient channels",
    "gostshyp_conventions": "channel semantics and negative controls",
}


def pytest_configure(config):
    for name, description in MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {description}")
