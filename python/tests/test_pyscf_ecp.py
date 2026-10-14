"""ECP nuclei in isodensity cavities.

An ECP pseudo-density can fall below ``rho_iso`` at its own nucleus, which
gives the isosurface a spurious inner sheet around it. moist fills every
nucleus whose nuclear charge, the host's ``nuclear_charges`` answer, lies below
its atomic number inside a ball that stays clear of the outer surface
(``core_fill``):

- a pseudo-atom whose density vanishes as ``r**4`` at its nucleus shows the
  sheet and its removal without any SCF
- AuH, 60 of gold's electrons behind the def2 ECP, reproduces the cavity of
  its all-electron ANO-RCC (X2C) density

Each test fails without an active fill:

- the outer surface lies far outside the fill and must stay bitwise unchanged;
  that holds trivially without a fill, so these tests also raise the isovalue
  until the surface reaches into the fill ball, where an active fill must
  refuse the build
- every projector but one keeps to the outer surface, the inner sheet being
  some 4 bohr from every anchor; the Newton-deflation projector
  (``proj_level=6``) converges onto it for some anchors when the fill is
  missing, so with the fill it must return the default projector's cavity
"""

import os

import numpy as np
import pytest

try:
    from pyscf import gto, scf
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from moist import Context, DensityMatrixRequest, GaussianBasis, NuclearChargesRequest, Structure
from moist.parameters import DROPParameters, IsodensityParameters
from moist.pyscf import DROP, Isodensity, PySCFHost

#: Run context shared by every cavity in this module
CONTEXT = Context()
#: Projection level whose deflated Newton search reaches the inner sheet
DEFLATION = 6
#: Distance from an ECP nucleus, in bohr, inside which a surface point lies on the inner sheet
INNER = 1.5


def _cavity(answers, structure, *, core_fill=True, rho_iso=None, **drop):
    """Internal isodensity cavity over ``answers = (basis, dcart, z)``, built as a host would."""
    basis, dcart, z = answers
    iso = {} if rho_iso is None else dict(rho_iso=rho_iso, scale=1.0 / rho_iso)
    lsf = Isodensity(parameters=IsodensityParameters(core_fill=core_fill, **iso))
    cavity = DROP(lsf=lsf, parameters=DROPParameters(**drop)).build(basis=basis, context=CONTEXT)
    cavity.update(structure)
    coupling = cavity.new_coupling()
    cavity.prepare_build(coupling)
    for request in coupling:
        if isinstance(request, DensityMatrixRequest):
            coupling.answer(dcart=dcart)
        elif isinstance(request, NuclearChargesRequest):
            coupling.answer(z=z)
        else:
            raise AssertionError(f"the build phase asked for {request.name!r}")
    cavity.build(coupling)
    return cavity


def _assert_fill_active(source, structure, rho_iso):
    """At an isovalue whose outer surface lies inside the 2 bohr fill, only an active fill refuses the build."""
    _cavity(source, structure, core_fill=False, rho_iso=rho_iso)
    with pytest.raises(RuntimeError, match="ECP core fill"):
        _cavity(source, structure, rho_iso=rho_iso)


def _assert_matches_default_projector(source, structure, **drop):
    """The deflated projector, which reaches an unfilled inner sheet, returns the default cavity."""
    deflated = _cavity(source, structure, proj_level=DEFLATION, **drop)
    reference = _cavity(source, structure, **drop)
    assert len(deflated.xyz) == len(reference.xyz)
    # Measured 9e-13: the two projectors converge onto the same points
    assert deflated.area == pytest.approx(reference.area, rel=1e-9)
    return deflated


def _inner_points(cavity, center):
    return int(np.sum(np.linalg.norm(cavity.xyz - center, axis=1) < INNER))


def _assert_same_surface(left, right):
    np.testing.assert_array_equal(left.xyz, right.xyz)
    np.testing.assert_array_equal(left.a, right.a)


# -----------------------------------------------------------------------------
# Pseudo-atom
# -----------------------------------------------------------------------------

_PSEUDO_BASIS = GaussianBasis(shell_atom=[0], shell_l=[2], shell_nprim=[1],
                              exponents=[0.5], coefficients=[1.0])


def _pseudo_atom(charge):
    """Xe nucleus carrying a pure d-shell density, which vanishes as ``r**4``.

    ``charge`` is the host's nuclear charge; below 54 it marks an ECP nucleus
    """
    answers = (_PSEUDO_BASIS, np.eye(_PSEUDO_BASIS.ncart), [charge])
    return answers, Structure([54], [[0.0, 0.0, 0.0]])


def test_pseudo_atom_outer_surface_is_untouched():
    ecp, structure = _pseudo_atom(8)
    all_electron, _ = _pseudo_atom(54)
    filled = _cavity(ecp, structure)
    _assert_same_surface(filled, _cavity(ecp, structure, core_fill=False))
    _assert_same_surface(filled, _cavity(all_electron, structure))
    # At rho_iso = 0.4 the outer surface lies at 1.6 bohr; a charge of Z is never filled
    _assert_fill_active(ecp, structure, 0.4)
    _cavity(all_electron, structure, rho_iso=0.4)


def test_pseudo_atom_inner_sheet_is_filled():
    # The d-shell density vanishes at the nucleus, so the sheet exists without the fill
    source, structure = _pseudo_atom(8)
    deflated = _assert_matches_default_projector(source, structure, nleb=110)
    assert _inner_points(deflated, 0.0) == 0


# -----------------------------------------------------------------------------
# AuH: def2-SVP with a 60-electron ECP against all-electron ANO-RCC
# -----------------------------------------------------------------------------

AUH = "Au 0 0 0; H 0 0 1.524"


def _host(mol, mf):
    mf.kernel()
    assert mf.converged
    host = PySCFHost(mol)
    host.dm = mf.make_rdm1()
    return host


def _answers(host):
    """The host's build-phase answers: its basis, ``dcart`` and ``z``."""
    return host.basis, host.cartesian_density(), host.nuclear_charges


def _structure(mol):
    # PySCFHost.structure() still refuses ECPs, so the elements are spelled out;
    # the nuclear charges reach moist in the build walk
    numbers = [gto.charge(mol.atom_pure_symbol(i)) for i in range(mol.natm)]
    return Structure(numbers, mol.atom_coords())


@pytest.fixture(scope="module")
def auh_ecp():
    mol = gto.M(atom=AUH, basis="def2-svp", ecp="def2-svp", verbose=0)
    assert mol.atom_nelec_core(0) == 60
    return mol, _host(mol, scf.RHF(mol)), _structure(mol)


@pytest.fixture(scope="module")
def auh_all_electron():
    # ANO-RCC ships with PySCF; x2c-SVPall would need basis-set-exchange
    mol = gto.M(atom=AUH, basis={"Au": "ano@7s6p4d2f", "H": "ano@2s1p"}, verbose=0)
    assert not mol.has_ecp()
    return mol, _host(mol, scf.RHF(mol).x2c()), _structure(mol)


def test_auh_ecp_cavity_matches_all_electron(auh_ecp, auh_all_electron):
    mol, host, structure = auh_ecp
    _, ae_host, ae_structure = auh_all_electron
    ecp = _cavity(_answers(host), structure)
    all_electron = _cavity(_answers(ae_host), ae_structure)
    gold = mol.atom_coord(0)

    # Measured 0.3 % in area and 0.03 bohr in closest approach (def2-SVP vs
    # ANO-RCC); the ECP only replaces the core, so the outer density agrees
    assert ecp.area == pytest.approx(all_electron.area, rel=0.01)
    closest = [np.linalg.norm(c.xyz - gold, axis=1).min() for c in (ecp, all_electron)]
    assert closest[0] == pytest.approx(closest[1], abs=0.1)


def test_auh_core_fill_leaves_outer_surface_untouched(auh_ecp):
    _, host, structure = auh_ecp
    source = _answers(host)
    _assert_same_surface(_cavity(source, structure), _cavity(source, structure, core_fill=False))
    # The ECP density stays below 0.2 beyond 1.8 bohr from gold
    _assert_fill_active(source, structure, 0.25)


def test_auh_inner_sheet_is_filled(auh_ecp):
    mol, host, structure = auh_ecp
    gold = mol.atom_coord(0)
    # Without the fill the sheet exists: the ECP density dips below rho_iso at the nucleus
    assert host.density(gold, 1)[0] < IsodensityParameters().rho_iso
    deflated = _assert_matches_default_projector(_answers(host), structure, nleb=26)
    assert _inner_points(deflated, gold) == 0
