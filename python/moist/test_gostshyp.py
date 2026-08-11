"""End-to-end tests of the GOSTSHYP pressure wall on a moist isodensity cavity.

GOSTSHYP is almost entirely QM-side integral work; what moist contributes is
the response of the cavity itself.  These tests close that loop: every analytic
quantity is checked against a finite difference of the energy the model
actually reports, at 50 GPa.

The suite is layered so a failure localises:

``integrals``
    The libcint fakemol constants, the cartesian d/f component orders, and the
    identity ``f == n . grad_r g``.  moist is used only to place the  grid points --
    a failure here is integral bookkeeping, not the cavity.
``params``
    Surface-parameter derivatives of ``gtilde``/``ftilde`` at a *frozen* cavity,
    against their own finite-difference reference.  These are the moments every
    surface weight is built from.
``fock``
    Energy and ``dE/dP`` with the surface following the density -- the level-set
    chain rule.
``gradient``
    ``dE/dR`` at fixed density, split into its integral, field and anchor
    channels and then assembled.
``conventions``
    Channel semantics and the negative controls.

Each layer carries a negative control, because a finite-difference test whose
extra term is numerically negligible passes while testing nothing.  GOSTSHYP is
unusually sharp here: it has no electrostatic component at all, so there is no
dominant term for a broken level-set route to hide behind.  The cavity response
is ~84% of ``dE/dP`` rather than a correction, and the anchor channel is
comparable in size to the other two throughout the nuclear gradient.
"""

import functools
import math
from dataclasses import dataclass

import numpy as np
import pytest

try:
    # Import what is actually used rather than the top-level package alone, and
    # skip on any ImportError: pyscf imports fine while pyscf.scf fails when a
    # broken h5py shadows the HDF5 runtime, and importorskip re-raises for an
    # ImportError raised by a transitive dependency rather than the named module.
    from pyscf import dft, gto, scf
except ImportError as exc:  # pragma: no cover - optional dependency
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from . import gostshyp as gostshyp_module  # noqa: E402
from .gostshyp import (  # noqa: E402
    _D_CART_ORDER,
    _F_RHO2_FIRST_MOMENT,
    _S_NORM,
    _S_OVER_D_NORM,
    _S_OVER_F_NORM,
    _S_OVER_P_NORM,
    GPA_TO_AU,
    GostshypWall,
    _int3c1e,
)
from .interface import Structure  # noqa: E402
from .pyscf import PySCFIsodensityHost  # noqa: E402

#: Every test in this module is a GOSTSHYP test; the second marker selects the
#: layer, and meson turns each into its own target.
pytestmark = pytest.mark.gostshyp

#: 50 GPa in Hartree / bohr^3.
PRESSURE = 50.0 * GPA_TO_AU
#: Lebedev order
NLEB = 26
#: Cavity projection tolerance
PROJ_TOL = 1e-12
#: FD step on the density matrix, in units of a unit-Frobenius-norm direction.
STEP_DM = 1e-4
#: FD step on nuclear coordinates in bohr; matches test_pyscf.py's tuned value.
STEP_R = 2.5e-4
#: FD step on a grid point centre, in bohr.
STEP_C = 1e-4
#: *Relative* FD step on the Gaussian exponent.  omega = pi ln2 / a spans orders
#: of magnitude across the  grid points, so one absolute step cannot serve them all.
STEP_W = 1e-4

#: Tolerances, combined as test-drive's ``check(..., thr_abs=, thr_rel=)`` does:
#: a value passes when ``|diff| <= max(ABS_THR, REL_THR * |reference|)``.
#:
#: This tier covers every finite difference that rebuilds the cavity, which is
#: the dominant noise source -- each stencil point re-runs the projection Newton
#: solve.  Measured worst deviation over the full matrix including the slow row:
#: 0.7% of tolerance for the Fock (fluoroacetate/STO-3G) and 0.3% for the
#: nuclear gradient (water/def2-SVP), i.e. ~140x and ~300x of headroom.
#:
#: That margin is larger than the measurement alone would justify and is
#: deliberate: the numbers come from one BLAS on one platform, and the limiting
#: noise is the cavity rather than the difference -- the FD values agree with
#: each other to ~1e-12 across steps from 5e-5 to 4e-4, and tightening PROJ_TOL
#: from 1e-12 to 1e-14 does not move the deviation at all.  Sensitivity does not
#: suffer for it: injecting a 0.0025% error into any single factor of the
#: analytic derivatives fails ten of these tests, and 0.05% into the area fold
#: fails six.
REL_THR = 1e-9
ABS_THR = REL_THR / 10.0

#: Frozen-cavity finite differences: no projection solve is repeated, so these
#: are limited only by FD truncation on smooth integrals.
FROZEN_REL_THR = 1e-7
FROZEN_ABS_THR = FROZEN_REL_THR / 10.0

#: The surface-parameter derivatives use a plain 2-point central difference, so
#: they carry an O(h^2) truncation error that no tightening of the cavity
#: removes.  Measured worst deviation 2.7% of tolerance (fluoroacetate/STO-3G,
#: ``dfdw``), the exponent route being the least well conditioned of the four.
PARAM_REL_THR = 1e-6

#: Independent-quadrature comparison for the fakemol constants.  A level-9
#: unpruned Becke grid integrates these Gaussian products to ~1e-13.
QUAD_ATOL = 1e-10
QUAD_RTOL = 1e-8

#: A negative control must miss by at least this many tolerances
VACUITY_FACTOR = 100.0
#: An analytic quantity must be at least this large before agreeing with a
#: finite difference says anything at all.
MIN_SIGNAL = 1e-12


@dataclass(frozen=True)
class System:
    """A test solute.  Geometries are in Angstrom."""

    atom: str
    charge: int = 0


#: Test solutes.  These are *fixtures*, not reference structures: bond lengths
#: are idealised and the geometries are not optimised.  What the tests need is a
#: fixed, well-conditioned geometry with a converging SCF.
SYSTEMS = {
    # Slightly asymmetric water, so no cartesian component is degenerate.
    "water": System(
        """O  0.0000  0.0000 -0.3893
           H  0.7629  0.0000  0.1947
           H -0.7991  0.0953  0.2223"""
    ),
    "fluoroacetate": System(
        """C  0.000  0.000  0.000
           C  1.540  0.000  0.000
           O  2.150  1.070  0.000
           O  2.150 -1.070  0.000
           F -0.760  0.000  1.180
           H -0.360  0.900 -0.520
           H -0.360 -0.900 -0.520""",
        charge=-1,
    ),
}

#: Active (system, basis) cases.  GOSTSHYP retains dense ``(nao, nao, ngrid)``
#: tensors and, for the response moments, a ``(ncart, ncart, ngrid, 10)``
#: f-fakemol block, so cost grows as nao^2.
#:
#: Basis choice is not arbitrary: the *grid point* leg carries d and f fakemols for
#: every basis, but the AO leg's cartesian-to-spherical transform is the
#: identity without l >= 2 AOs.  STO-3G has no d functions on O/F, so a
#: cart/sph mix-up in ``_to_spherical`` is invisible to it -- keep a def2-SVP row.
CASES = [
    ("water", "sto-3g"),
    ("water", "def2-svp"),
]
#: Larger rows, kept out of the default sweep.
SLOW_CASES = [
    ("fluoroacetate", "sto-3g"),
]
#: The case used by tests that pin a convention once rather than sweeping.
PRIMARY_CASE = CASES[0]


def _case_params():
    """Parametrisation entries, marking everything outside CASES slow."""
    entries = []
    for system, basis in CASES + SLOW_CASES:
        marks = () if (system, basis) in CASES else (pytest.mark.slow,)
        entries.append(pytest.param(system, basis, id=f"{system}-{basis}", marks=marks))
    return entries


CASE_PARAMS = _case_params()


def deviation(actual, reference, *, thr_abs=None, thr_rel=None) -> float:
    """Deviation measured in units of the tolerance: ``<= 1`` passes.

    Combines the two thresholds the way test-drive's
    ``check(..., thr_abs=, thr_rel=)`` does -- ``max(thr_abs, thr_rel*|ref|)`` --
    so the absolute floor covers references near zero while the relative bound
    scales with the magnitude.  Reporting the ratio rather than the raw
    difference makes a failure message say how many tolerances were missed.
    """
    thr_abs = ABS_THR if thr_abs is None else thr_abs
    thr_rel = REL_THR if thr_rel is None else thr_rel
    return abs(actual - reference) / max(thr_abs, thr_rel * abs(reference))


def array_deviation(actual, reference, *, thr_rel, thr_abs=None) -> float:
    """Worst elementwise deviation, scaled by the *array's* magnitude.

    Per-element relative errors are meaningless for arrays whose entries span
    orders of magnitude -- the small entries are differences of large ones -- so
    the tolerance is set from the largest entry of the reference.
    """
    actual = np.asarray(actual)
    reference = np.asarray(reference)
    scale = float(np.max(np.abs(reference), initial=0.0))
    threshold = max(thr_abs if thr_abs is not None else 0.0, thr_rel * scale)
    return float(np.max(np.abs(actual - reference), initial=0.0)) / threshold


FD4_OFFSETS = (2, 1, -1, -2)


def fd4(values, step: float) -> float:
    """4-point central difference from samples at ``(+2, +1, -1, -2) * step``."""
    fpp, fp, fm, fmm = values
    return (-fpp + 8.0 * fp - 8.0 * fm + fmm) / (12.0 * step)


def fd2(plus, minus, step):
    """2-point central difference, for the array-valued references."""
    return (np.asarray(plus) - np.asarray(minus)) / (2.0 * step)


# ----------------------------------------------------------------------
# system construction (cached: every parametrised test reuses them)
# ----------------------------------------------------------------------


@functools.lru_cache(maxsize=None)
def molecule(system: str, basis: str):
    spec = SYSTEMS[system]
    return gto.M(
        atom=spec.atom, basis=basis, charge=spec.charge, unit="Angstrom", verbose=0
    )


@functools.lru_cache(maxsize=None)
def reference_density(system: str, basis: str):
    """Converged gas-phase RHF density, then held fixed as a free parameter."""
    mean_field = scf.RHF(molecule(system, basis))
    mean_field.conv_tol = 1e-12
    mean_field.kernel()
    assert mean_field.converged, f"gas-phase SCF failed for {system}/{basis}"
    return mean_field.make_rdm1()


def make_wall(mol, positions=None, *, dm, pressure=PRESSURE):
    """Host plus an updated GOSTSHYP wall at ``dm``, optionally displaced."""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    host = PySCFIsodensityHost(mol)
    host.dm = dm
    wall = GostshypWall(host, pressure, nleb=NLEB, tolerance=PROJ_TOL)
    wall.update(dm)
    return host, wall


def fd_density(mol, dm, direction, *, pressure=PRESSURE):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        _host, wall = make_wall(mol, dm=dm + offset * STEP_DM * direction, pressure=pressure)
        samples.append(wall.energy)
        grids.append(wall.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, pressure=PRESSURE):
    """dE/dR along one cartesian coordinate at fixed density matrix."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        displaced = positions.copy()
        displaced.flat[index] += offset * STEP_R
        _host, wall = make_wall(mol, displaced, dm=dm, pressure=pressure)
        samples.append(wall.energy)
        grids.append(wall.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_R)


def symmetric_directions(nao, count, seed=11):
    """Unit-norm symmetric density-matrix perturbations."""
    rng = np.random.default_rng(seed)
    for _ in range(count):
        direction = rng.standard_normal((nao, nao))
        direction = 0.5 * (direction + direction.T)
        yield direction / np.linalg.norm(direction)


def sampled_coordinates(natm):
    """Cartesian coordinates to difference: all of them for a small solute."""
    total = 3 * natm
    if total <= 9:
        return list(range(total))
    return [int(i) for i in np.linspace(0, total - 1, 3, dtype=int)]


def frozen_surface_energy(mol, dm, centers, areas, omega, normals, pressure=PRESSURE):
    """GOSTSHYP energy for an explicit surface, independent of any wall state.

    Deliberately re-derived from :func:`~moist.gostshyp._int3c1e` rather than
    calling into the model, so the finite differences that use it cannot inherit
    a bug from the code under test.
    """
    ngrid = int(np.asarray(centers).shape[0])
    g_cart = _int3c1e(mol, centers, omega, 0)
    p_cart = _int3c1e(mol, centers, omega, 1)
    ncart = g_cart.shape[0]
    p_cart = p_cart.reshape(ncart, ncart, ngrid, 3)

    c2s = np.eye(mol.nao_nr()) if mol.cart else np.asarray(mol.cart2sph_coeff(normalized="sp"))
    # BLAS leaves dirty FP status flags for these shapes; see _density_matrix_cart.
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        dm_cart = c2s @ np.asarray(dm) @ c2s.T

    gtilde = np.einsum("pqj,pq->j", g_cart, dm_cart, optimize=True)
    moment = _S_OVER_P_NORM * np.einsum("pqja,pq->ja", p_cart, dm_cart, optimize=True)
    ftilde = -2.0 * np.asarray(omega) * np.einsum(
        "ja,ja->j", np.asarray(normals), moment, optimize=True
    )

    floor = gostshyp_module._OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = pressure * np.asarray(areas) / ftilde
    alpha = np.where(np.isfinite(alpha) & (np.abs(ftilde) > floor), alpha, 0.0)
    return float(alpha @ gtilde)


# ----------------------------------------------------------------------
# integrals -- libcint conventions, moist only places the  grid points
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_integrals
def test_integrals_fakemol_normalization_constants():
    """The four angular constants and the cartesian d/f orders, vs quadrature.

    libcint attaches a different normalization to a coefficient-1 shell of each
    angular momentum.  Restoring the ratios ``N_s/N_l`` is what puts every
    moment in the same units as ``g`` and makes ``f = n . grad g`` exact rather
    than exact-up-to-a-factor.  An independent Becke-grid quadrature of the same
    monomial moments pins the constants *and* the component orders together: a
    transposed ``_D_CART_ORDER`` entry or a mis-assigned
    ``_F_RHO2_FIRST_MOMENT`` triple is invisible to any self-consistency check.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # A handful of well-conditioned  grid points is enough and keeps the grid cheap.
    picks = np.argsort(-np.abs(wall._state.ftilde))[:4]
    centers = wall.centers[picks]
    omega = wall.omega[picks]

    grids = dft.gen_grid.Grids(mol)
    grids.level = 9
    grids.prune = None
    grids.build()
    ao = mol.eval_gto("GTOval_cart", grids.coords)
    weighted = ao * grids.weights[:, None]

    for slot in range(len(picks)):
        delta = grids.coords - centers[slot]
        rho2 = np.einsum("ga,ga->g", delta, delta)
        gauss = np.exp(-omega[slot] * rho2)

        def quad(factor):
            return np.einsum("gu,gv,g->uv", weighted, ao, gauss * factor, optimize=True)

        reference_s = _S_NORM * quad(np.ones_like(rho2))
        actual_s = _int3c1e(mol, centers, omega, 0)[:, :, slot]
        np.testing.assert_allclose(actual_s, reference_s, rtol=QUAD_RTOL, atol=QUAD_ATOL)

        p_block = _int3c1e(mol, centers, omega, 1)
        ncart = p_block.shape[0]
        p_block = p_block.reshape(ncart, ncart, len(picks), 3)
        for axis in range(3):
            actual = _S_OVER_P_NORM * p_block[:, :, slot, axis]
            np.testing.assert_allclose(
                actual, _S_NORM * quad(delta[:, axis]), rtol=QUAD_RTOL, atol=QUAD_ATOL
            )

        d_block = _int3c1e(mol, centers, omega, 2).reshape(ncart, ncart, len(picks), 6)
        for component, (a, b) in enumerate(_D_CART_ORDER):
            actual = _S_OVER_D_NORM * d_block[:, :, slot, component]
            np.testing.assert_allclose(
                actual,
                _S_NORM * quad(delta[:, a] * delta[:, b]),
                rtol=QUAD_RTOL,
                atol=QUAD_ATOL,
            )

        f_block = _int3c1e(mol, centers, omega, 3).reshape(ncart, ncart, len(picks), 10)
        for axis, triple in enumerate(_F_RHO2_FIRST_MOMENT):
            actual = _S_OVER_F_NORM * f_block[:, :, slot, list(triple)].sum(axis=-1)
            np.testing.assert_allclose(
                actual,
                _S_NORM * quad(delta[:, axis] * rho2),
                rtol=QUAD_RTOL,
                atol=QUAD_ATOL,
            )


@pytest.mark.gostshyp_integrals
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_integrals_f_is_the_field_point_normal_gradient_of_g(system, basis):
    """``ftilde == n . grad_r gtilde``, i.e. minus the grid point-centre gradient.

    Displacing the field point is the opposite of displacing the Gaussian
    centre, which is the whole content of the sign in ``_build_integrals``.  Get
    it wrong and the wall pushes outward instead of inward.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    dm_cart = wall._density_matrix_cart(dm)
    ftilde = wall._state.ftilde

    samples = []
    for offset in FD4_OFFSETS:
        centers = wall.centers + offset * STEP_C * wall.normals
        samples.append(wall._gf_tilde(dm_cart, centers, wall.omega, wall.normals)[0])
    centre_gradient = fd4(samples, STEP_C)

    assert array_deviation(ftilde, -centre_gradient, thr_rel=FROZEN_REL_THR) <= 1.0


@pytest.mark.gostshyp_integrals
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_integrals_ftilde_is_the_normal_projection_of_the_f_vector(system, basis):
    """``ftilde_j == n_j . gvfield_j`` exactly.

    ``gvfield`` is the unprojected f-vector and doubles as ``dftilde/dn``, so
    the two must agree to round-off or the normal weight is built from a
    different quantity than the one it differentiates.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gvfield = np.einsum("uvja,uv->ja", wall._Fvec, dm, optimize=True)
    projected = np.einsum("ja,ja->j", wall.normals, gvfield, optimize=True)
    np.testing.assert_allclose(projected, wall._state.ftilde, rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# params -- surface-parameter derivatives at a frozen cavity
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_derivatives_match_fd(system, basis):
    """All four moment-built derivatives of ``gtilde``/``ftilde``.

    These are the raw material of every surface weight: the position
    derivatives come from the p and d moments, the exponent derivatives from the
    d and f moments.  Checking them here means a later failure in ``fock`` or
    ``gradient`` is a chain-rule error, not a bad moment.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    dm_cart = wall._density_matrix_cart(dm)

    analytic = wall._param_derivatives(dm_cart)
    numerical = wall._param_derivatives_fd(dm_cart, h_r=STEP_C, h_w=STEP_W)

    for name, exact, approx in zip(("dgdr", "dfdr", "dgdw", "dfdw"), analytic, numerical):
        assert array_deviation(exact, approx, thr_rel=PARAM_REL_THR) <= 1.0, name


@pytest.mark.gostshyp_params
def test_params_second_moment_matches_a_second_difference():
    """``Mt_ab`` against a second difference of ``gtilde`` in the grid point centre.

    With ``G = exp(-omega rho^2)``, ``d^2 gtilde/dC_a dC_b = 4 omega^2 Mt_ab -
    2 omega delta_ab gtilde``.  This isolates the d-fakemol constant and the
    ``_D_CART_ORDER`` mapping, which ``test_params_derivatives_match_fd`` only
    sees through the contracted combination ``n . Mt``.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    dm_cart = wall._density_matrix_cart(dm)

    _gt, _pt, moment, _rt = wall._surface_moments(dm_cart)
    gtilde = wall._state.gtilde
    omega = wall.omega
    step = 2.0e-3

    def gt_at(shift):
        centers = wall.centers + shift
        return wall._gf_tilde(dm_cart, centers, omega, wall.normals)[0]

    for a in range(3):
        for b in range(3):
            ea = np.zeros(3)
            ea[a] = step
            eb = np.zeros(3)
            eb[b] = step
            second = (gt_at(ea + eb) - gt_at(ea - eb) - gt_at(-ea + eb) + gt_at(-ea - eb)) / (
                4.0 * step * step
            )
            expected = 4.0 * omega**2 * moment[:, a, b] - 2.0 * omega * gtilde * (a == b)
            assert array_deviation(expected, second, thr_rel=1e-4) <= 1.0, (a, b)


# ----------------------------------------------------------------------
# fock -- the level-set chain rule
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_fock
def test_fock_energy_is_linear_in_pressure_and_sets_the_effective_volume():
    """``E = p V_eff`` with ``V_eff`` independent of ``p``, and ``E(0) = 0``.

    Pure identities on one cavity: the surface does not depend on the pressure,
    and the amplitudes are proportional to it.  A mis-masked ``active`` would
    break the volume identity without breaking linearity.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)

    _host, wall = make_wall(mol, dm=dm)
    _host2, doubled = make_wall(mol, dm=dm, pressure=2.0 * PRESSURE)
    _host3, vacuum = make_wall(mol, dm=dm, pressure=0.0)

    assert wall.energy > 0.0
    assert doubled.energy == pytest.approx(2.0 * wall.energy, rel=1e-12)
    assert vacuum.energy == 0.0
    assert wall.effective_volume() == pytest.approx(wall.energy / PRESSURE, rel=1e-12)
    # Not the cavity volume: the wall is weighted by the local density overlap.
    assert wall.effective_volume() != pytest.approx(wall.cavity.cavity.volume, rel=1e-3)

    # The three walls share a density and therefore a surface, so the volume is
    # common to them.  At p = 0 the energy vanishes with the pressure while the
    # volume does not, which is why it cannot be recovered as E/p there.
    assert vacuum.effective_volume() == pytest.approx(wall.effective_volume(), rel=1e-12)
    assert doubled.effective_volume() == pytest.approx(wall.effective_volume(), rel=1e-12)


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_fock_frozen_surface_matches_fd(system, basis):
    """The eq-16 Fock is ``dE/dP`` with the  grid points held fixed.

    Isolates the explicit density dependence from the cavity response, so a
    failure in the headline test below can be attributed to one or the other.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    for direction in symmetric_directions(dm.shape[0], 2):
        samples = [
            frozen_surface_energy(
                mol,
                dm + offset * STEP_DM * direction,
                wall.centers,
                wall.areas,
                wall.omega,
                wall.normals,
            )
            for offset in FD4_OFFSETS
        ]
        numerical = fd4(samples, STEP_DM)
        analytic = float(np.einsum("uv,uv->", wall._state.fock, direction))
        assert deviation(
            analytic, numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR
        ) <= 1.0


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_fock_matches_fd(system, basis):
    """dE/dP with the surface following the density -- the headline Fock test.

    Closes the whole chain: the surface weights, moist's contraction of them
    into level-set adjoints, and the host's contraction of those with dS/dP.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    fock = wall.fock(dm)

    for direction in symmetric_directions(dm.shape[0], 2):
        numerical = fd_density(mol, dm, direction)
        analytic = float(np.einsum("uv,uv->", fock, direction))
        assert deviation(analytic, numerical) <= 1.0


# ----------------------------------------------------------------------
# gradient -- three nuclear routes
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_gradient
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_gradient_integral_channel_matches_frozen_surface_fd(system, basis):
    """AO centres move,  grid points frozen: the ``int3c1e_ip1`` term alone.

    Pins the sign of the integral derivative, the both-legs factor of two, and
    the mapping from cartesian AO rows onto atoms.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)
    analytic = wall._integral_nuclear_gradient(dm).flatten(order="F")

    for index in sampled_coordinates(mol.natm):
        samples = []
        for offset in FD4_OFFSETS:
            displaced = positions.copy()
            displaced.flat[index] += offset * STEP_R
            moved = mol.set_geom_(displaced, unit="Bohr", inplace=False)
            samples.append(
                frozen_surface_energy(
                    moved, dm, wall.centers, wall.areas, wall.omega, wall.normals
                )
            )
        numerical = fd4(samples, STEP_R)
        assert deviation(
            analytic[index], numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR
        ) <= 1.0


@pytest.mark.gostshyp_gradient
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_gradient_anchor_channel_matches_fd(system, basis):
    """Rigid owner-atom motion of the  grid points at a frozen level-set field.

    Displacing only the structure handed to the cavity -- never the molecule the
    level set is built from -- drags the atom-anchored reference grid while the
    density stays put.  That is exactly the channel moist's anchor pass reports,
    and it validates both ``xyz1_rA`` (including the normal's point-motion fold)
    and ``a_i1_rA`` (the switching-inclusive area route) at once.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host, wall = make_wall(mol, dm=dm)
    analytic = wall._anchor_nuclear_gradient(dm).flatten(order="F")
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)

    def anchored_energy(displaced):
        # The cavity moves; host.mol -- and hence the level set -- does not.
        cavity = host.make_cavity(nleb=NLEB, tolerance=PROJ_TOL)
        cavity.update(Structure(numbers, displaced))
        result = cavity.cavity
        centers = np.ascontiguousarray(result.xyz.T)
        areas = np.ascontiguousarray(result.a)
        normals = np.ascontiguousarray(result.normal0.T)
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            omega = np.pi * math.log(2.0) / areas
        omega[~np.isfinite(omega)] = 0.0
        return frozen_surface_energy(mol, dm, centers, areas, omega, normals), result.ngrid

    for index in sampled_coordinates(mol.natm):
        samples, grids = [], []
        for offset in FD4_OFFSETS:
            displaced = positions.copy()
            displaced.flat[index] += offset * STEP_R
            energy, ngrid = anchored_energy(displaced)
            samples.append(energy)
            grids.append(ngrid)
        assert len(set(grids)) == 1, f"grid point count drifted: {grids}"
        numerical = fd4(samples, STEP_R)
        assert deviation(analytic[index], numerical) <= 1.0


@pytest.mark.gostshyp_gradient
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_gradient_matches_fd(system, basis):
    """Total dE/dR at fixed density -- the headline nuclear-gradient test."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)
    analytic = wall.nuclear_gradient(dm).flatten(order="F")

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, positions, dm, index)
        assert deviation(analytic[index], numerical) <= 1.0


@pytest.mark.gostshyp_gradient
def test_gradient_is_translationally_invariant():
    """The net force vanishes.

    A rigid translation moves the AOs, the density and the cavity together, so
    the energy cannot change.  This is a global identity that knows nothing
    about the channel split, which is what makes it able to catch a sign error
    or a mis-binned AO-to-atom map that each per-channel test tolerates.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    gradient = wall.nuclear_gradient(dm)

    assert np.abs(gradient).max() > MIN_SIGNAL
    np.testing.assert_allclose(gradient.sum(axis=1), 0.0, atol=1e-10)


@pytest.mark.gostshyp_gradient
def test_gradient_is_exactly_the_sum_of_its_channels():
    """No channel is silently dropped inside :meth:`nuclear_gradient`.

    Also the correctness half of the shared-adjoint optimisation: the channels
    here build their own surface adjoints, so bit-level agreement proves that
    reusing one instance across the field and anchor routes changes nothing.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    channels = (
        wall._integral_nuclear_gradient(dm)
        + wall._field_nuclear_gradient(dm)
        + wall._anchor_nuclear_gradient(dm)
    )
    np.testing.assert_allclose(wall.nuclear_gradient(dm), channels, rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# conventions -- channel semantics and the negative controls
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_conventions
def test_conventions_surface_weights_leave_w_f_zero():
    """GOSTSHYP has no dependence on the DROP focusing criterion.

    The switching function is an anchor-only iSwiG overlap, so ``df/dS = 0`` for
    a callback level set and the whole area route runs through the Gaussian
    width.  Exactly zero, not merely small.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    adjoints = wall.surface_adjoints(dm)

    np.testing.assert_array_equal(adjoints.w_f, np.zeros_like(adjoints.w_f))
    # The routes that must *not* be zero, or the controls below are vacuous.
    assert np.abs(adjoints.w_xi).max() > MIN_SIGNAL
    assert np.abs(adjoints.w_xyz).max() > MIN_SIGNAL
    assert np.abs(adjoints.w_n).max() > MIN_SIGNAL


@pytest.mark.gostshyp_conventions
def test_conventions_fock_requires_cavity_response():
    """The eq-16 Fock alone is not ``dE/dP`` for a surface that follows the density.

    With no electrostatic component there is nothing to dilute the error: the
    cavity response is the majority of the derivative, not a correction.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    response = wall.fock(dm) - wall.fock(dm, include_cavity_response=False)
    assert np.abs(response).max() > MIN_SIGNAL

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction)
    starved = float(np.einsum("uv,uv->", wall.fock(dm, include_cavity_response=False), direction))
    assert deviation(starved, numerical) > VACUITY_FACTOR


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize("channel", ["w_n", "w_xi"])
def test_conventions_fock_requires_every_surface_route(channel):
    """Dropping either the normal or the area route breaks the Fock outright.

    ``w_n`` is what moist folds into both the direct ``P_tan/|grad S|`` term and
    the Hessian point-motion term; ``w_xi`` is the host-side stand-in for the
    ``w_a`` channel the C API does not carry.  Neither is a small correction.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host, wall = make_wall(mol, dm=dm)

    adjoints = wall.surface_adjoints(dm)
    w_xi = np.zeros_like(adjoints.w_xi) if channel == "w_xi" else adjoints.w_xi
    w_n = None if channel == "w_n" else np.asfortranarray(adjoints.w_n.T)
    w_lsf = wall.cavity.contract_surface_lsf_weights(
        w_xi, adjoints.w_f, np.asfortranarray(adjoints.w_xyz.T), w_n=w_n
    )

    weights = gostshyp_module._LsfWeights(*w_lsf)
    starved = wall._state.fock + host._fock_lsf(wall.centers, weights)

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction)
    analytic = float(np.einsum("uv,uv->", starved, direction))
    assert deviation(analytic, numerical) > VACUITY_FACTOR


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize("dropped", ["integral", "field", "anchor"])
def test_conventions_gradient_requires_every_channel(dropped):
    """Each of the three nuclear routes carries a non-negligible share."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)

    channels = {
        "integral": wall._integral_nuclear_gradient(dm),
        "field": wall._field_nuclear_gradient(dm),
        "anchor": wall._anchor_nuclear_gradient(dm),
    }
    assert np.abs(channels[dropped]).max() > MIN_SIGNAL
    starved = sum(value for name, value in channels.items() if name != dropped)

    index = int(np.argmax(np.abs(channels[dropped].flatten(order="F"))))
    numerical = fd_position(mol, positions, dm, index)
    assert deviation(starved.flatten(order="F")[index], numerical) > VACUITY_FACTOR


@pytest.mark.gostshyp_conventions
def test_conventions_anchor_area_route_needs_the_true_area_derivative():
    """``a_i1_rA`` is not recoverable from the Gaussian-width proxy.

    The grid point area is ``a_i ~ f_i / xi_i**2``, and under *nuclear* motion the
    switching factor moves too -- a route ``xi1_rA`` cannot see.  Reusing the
    Fock-side ``dE_da -> w_xi`` fold here is the natural mistake, so it is
    pinned as a control rather than left to a comment.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host, wall = make_wall(mol, dm=dm)

    adjoints = wall.surface_adjoints(dm)
    grad_norm, hess = wall._lsf_jet_at_grid_points()
    nwn = np.einsum("ia,ia->i", wall.normals, adjoints.w_n, optimize=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        normal_grad = (adjoints.w_n - wall.normals * nwn[:, None]) / grad_norm[:, None]
    normal_grad[~np.isfinite(normal_grad)] = 0.0
    w_xyz_total = adjoints.w_xyz + np.einsum("iab,ib->ia", hess, normal_grad, optimize=True)

    wall.cavity.compute_anchor_gradient()
    anchor = wall.cavity.get_anchor_gradient()
    proxy = np.einsum("ic,caAi->aA", w_xyz_total, anchor.xyz1_rA, optimize=True)
    # The wrong area route: fold through xi instead of using da_i/dR_A.
    proxy += np.einsum("i,aAi->aA", adjoints.w_xi, anchor.xi1_rA, optimize=True)

    starved = (
        wall._integral_nuclear_gradient(dm) + wall._field_nuclear_gradient(dm) + proxy
    ).flatten(order="F")

    misses = 0
    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, positions, dm, index)
        if deviation(starved[index], numerical) > VACUITY_FACTOR:
            misses += 1
    assert misses > 0


@pytest.mark.gostshyp_conventions
def test_conventions_overlap_floor_keeps_the_model_differentiable():
    """The activity threshold is applied to the energy, not just its derivative.

    Every GOSTSHYP quantity is a ratio ``gtilde_j/ftilde_j`` of two exponentially
    small numbers, so  grid points that have left the density carry no information.
    Dropping them is necessary, but *where* the cut falls is a modelling choice
    -- so the requirement is not that the value be insensitive to it (it is not;
    a looser floor discards real energy) but that the model stay exactly
    differentiable wherever it is placed.  That holds only because the threshold
    is applied once, to the amplitudes: gating only ``dE/da`` would leave the
    energy depending on  grid points whose sensitivity is defined to be zero, and
    the finite difference below would fail at every floor rather than none.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    direction = next(iter(symmetric_directions(dm.shape[0], 1)))

    original = gostshyp_module._OVERLAP_FLOOR
    energies = []
    try:
        for floor in (1e-6, 1e-9, 1e-12):
            gostshyp_module._OVERLAP_FLOOR = floor
            _host, wall = make_wall(mol, dm=dm)
            energies.append(wall.energy)
            analytic = float(np.einsum("uv,uv->", wall.fock(dm), direction))
            numerical = fd_density(mol, dm, direction)
            assert deviation(analytic, numerical) <= 1.0, f"floor={floor:.0e}"
    finally:
        gostshyp_module._OVERLAP_FLOOR = original

    # Tightening the floor readmits energy monotonically: a looser cut can only
    # discard  grid points, and every grid point's contribution p_j gtilde_j is
    # positive for a compressive wall.
    assert energies[0] < energies[1] < energies[2]
    # ...and the default discards under a part in 10^6 (measured 1.7e-7), while
    # the loosest floor tried here gives up an order of magnitude more.
    assert (energies[2] - energies[1]) / energies[2] < 1e-6
    assert (energies[2] - energies[0]) / energies[2] > 1e-6


@pytest.mark.gostshyp_conventions
def test_conventions_anchor_area_derivatives_sum_to_the_total():
    """``sum_i a_i1_rA == A_tot1_rA``, pinning the anchor buffer layout.

    The binding hands moist a rank-4 Fortran buffer and a capacity pair whose
    order is the reverse of the array indexing.  This identity is independent of
    any GOSTSHYP physics, so it bisects a layout error away from a weight error.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    wall.cavity.compute_anchor_gradient()
    anchor = wall.cavity.get_anchor_gradient()
    np.testing.assert_allclose(
        anchor.a_i1_rA.sum(axis=2), anchor.A_tot1_rA, rtol=0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        anchor.v_i1_rA.sum(axis=2), anchor.V_tot1_rA, rtol=0.0, atol=1e-12
    )


@pytest.mark.gostshyp_conventions
def test_conventions_failed_update_leaves_no_stale_results(monkeypatch):
    """A rebuild that raises must not leave the previous geometry readable.

    :meth:`update` moves the host density and the cavity before it can fail --
    a density with no isodensity surface, a projection that will not converge --
    so results cached from the previous call describe a surface that no longer
    exists.  The failure is forced here rather than hunted for: what is being
    pinned is the invalidation, not any particular way of tripping it.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # Vacuous unless there is a real result to go stale.
    assert wall.effective_volume() > 0.0
    assert wall.amplitudes is not None

    def boom(*_args, **_kwargs):
        raise RuntimeError("integral build failed")

    monkeypatch.setattr(wall, "_build_integrals", boom)
    with pytest.raises(RuntimeError, match="integral build failed"):
        wall.update(dm)

    assert wall.amplitudes is None
    assert wall.energy == 0.0
    for read in (wall.effective_volume, lambda: wall.fock(dm), lambda: wall.nuclear_gradient(dm)):
        with pytest.raises(RuntimeError, match="update"):
            read()


@pytest.mark.gostshyp_conventions
def test_conventions_gradient_builds_the_surface_adjoints_once(monkeypatch):
    """The field and anchor routes share one :meth:`surface_adjoints` result.

    Building them runs :meth:`_surface_moments`, which is four dense
    three-centre integral blocks -- the d and f shells among them -- and is the
    dominant cost of the whole gradient.  Computing them per route is correct
    but doubles that cost, and nothing about the returned numbers would show it,
    so the call count is pinned here.  Bit-level agreement between the total and
    its channels is checked by ``test_gradient_is_exactly_the_sum_of_its_channels``.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    original = wall.surface_adjoints
    calls = []

    def counted(dm_arg):
        calls.append(dm_arg)
        return original(dm_arg)

    monkeypatch.setattr(wall, "surface_adjoints", counted)
    gradient = wall.nuclear_gradient(dm)

    assert len(calls) == 1
    assert np.abs(gradient).max() > MIN_SIGNAL
