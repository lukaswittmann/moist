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
- at the default ``branch_dpdist_max`` the branch-search window stays near the
  outer surface, the inner sheet being some 4 bohr from every anchor; widening
  ``branch_dpdist_max`` at ``proj_level=3`` lets the coarse-octree search reach
  the unfilled inner cavity
- the sheet's normals point back at the nucleus, so the default horizon switch
  closes it; the search runs with the horizon switch off and a reach band past
  the widened window (``INNER_GATE``)
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
#: Coarse-octree projection level used to search for an unfilled inner sheet
GLOBAL_SEARCH = 3
#: Distance window (bohr) opening the level-3 search far enough to reach the nucleus
INNER_BRANCH_DPDIST_MAX = 3.6
#: Distance from an ECP nucleus, in bohr, inside which a surface point lies on the inner sheet
INNER = 1.5
#: Gate admitting the unfilled inner sheet, 3.4 to 4.0 bohr from its anchors:
#:
#: - horizon switch off: it closes every point whose normal points back into the
#:   anchor's sphere, as on the inner sheet
#: - reach open past the widened branch window (at most 5.5 bohr here), so the
#:   reach gate is one on every point; the default band, 3.78 to 5.67 bohr,
#:   would thin the pseudo-atom's sheet
INNER_GATE = dict(w_0hor_on=0.0, branch_reach_on=6.0, branch_reach_off=7.0)
#: Default ``reach_on`` in bohr, from which the default reach gate falls
DEFAULT_REACH_ON = 3.78
#: Lebedev points per atom of the inner-sheet search
NLEB = 302
#: Settings of the wide search for the inner sheet
SEARCH = dict(proj_level=GLOBAL_SEARCH, branch_dpdist_max=INNER_BRANCH_DPDIST_MAX, **INNER_GATE)


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


def _assert_inner_sheet_is_filled(source, structure, center, **drop):
    """The wide level-3 search reaches the unfilled inner sheet, and returns the default cavity with the fill.

    Returns the unfilled cavity
    """
    unfilled = _cavity(source, structure, core_fill=False, **SEARCH, **drop)
    filled = _cavity(source, structure, **SEARCH, **drop)
    reference = _cavity(source, structure, **drop)

    _assert_shows_inner_sheet(unfilled, center)
    assert len(unfilled.xyz) > len(reference.xyz), "the unfilled sheet adds no points to the default cavity"

    assert _inner_points(filled, center) == 0, "the core fill leaves points on the inner sheet"
    assert len(filled.xyz) == len(reference.xyz), "with the core fill the wide search changes the point count"
    np.testing.assert_allclose(filled.xyz, reference.xyz, rtol=0.0, atol=1e-10)
    assert filled.area == pytest.approx(reference.area, rel=1e-9)
    return unfilled


def _assert_shows_inner_sheet(cavity, center):
    """The wide search shows the inner sheet of an unfilled density, whose points the default gate closes."""
    inner = _inner_mask(cavity, center)
    n_inner = int(inner.sum())
    assert n_inner > 0, "the unfilled density shows no inner sheet: the wide search or its gate misses it"
    # The sheet's normals point back at the nucleus: the default horizon switch closes all of them
    inward = np.sum(cavity.normal0 * (cavity.xyz - center), axis=1) < 0.0
    n_inward = int(np.sum(inward & inner))
    assert n_inward == n_inner, "an inner-sheet normal points away from the nucleus"
    # Reach gate one on every point: only the core fill and the search decide what the cavity shows
    farthest = float(cavity.pdist.max())
    assert farthest < INNER_GATE["branch_reach_on"], "a point lies in the reach band of INNER_GATE"


def _inner_mask(cavity, center):
    return np.linalg.norm(cavity.xyz - center, axis=1) < INNER


def _inner_points(cavity, center):
    return int(np.sum(_inner_mask(cavity, center)))


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
    ecp, structure = _pseudo_atom(8)
    all_electron, _ = _pseudo_atom(54)
    unfilled = _assert_inner_sheet_is_filled(ecp, structure, 0.0, nleb=NLEB)
    # The d-shell sheet lies past the default reach_on, where the default band would thin it
    inner = _inner_mask(unfilled, 0.0)
    assert np.any(unfilled.pdist[inner] > DEFAULT_REACH_ON), "the inner sheet lies inside the default reach"
    # A charge of Z is never filled, whatever the density
    unfilled_ae = _cavity(all_electron, structure, **SEARCH, nleb=NLEB)
    _assert_shows_inner_sheet(unfilled_ae, 0.0)


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
    _assert_inner_sheet_is_filled(_answers(host), structure, gold, nleb=NLEB)
