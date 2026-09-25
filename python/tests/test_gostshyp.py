"""End-to-end tests of the GOSTSHYP pressure model on a moist isodensity cavity.

GOSTSHYP is almost entirely QM-side integral work; what moist contributes is
the response of the cavity itself.  These tests close that loop: every analytic
quantity is checked against a finite difference of the energy the model
actually reports, at 50 GPa.

The physics lives in moist's ``gostshyp`` component, so this file tests the
*host* half and the assembled result.  The component's own surface weights are
finite-differenced channel by channel in
``test/unit/test_model/component/gostshyp.f90``, against a model density whose
Gaussian moments are available in closed form.  Where a control moved there,
the comment at its old site says so.

The suite is layered so a failure localises:

``integrals``
    The libcint fakemol constants, the cartesian d/f component orders, and the
    identity ``f == n . grad_r g``.  moist is used only to place the grid points;
    a failure here is integral bookkeeping, not the cavity.  These constants
    are the one convention that stayed host-side, so this layer is load-bearing.
``params``
    The Gaussian moments this module hands the component, and the hand-off
    itself: a transposed moment is the failure mode the language boundary
    introduces.
``fock``
    Energy and ``dE/dP`` with the surface following the density -- the level-set
    chain rule.
``gradient``
    ``dE/dR`` at fixed density, split into its integral, field and surface
    channels and then assembled.  The first two are the host's; the third is
    moist contracting the component's weights in reverse mode.
``conventions``
    Channel semantics, the negative controls, and the frozen reference values.

Each layer carries a negative control, because a finite-difference test whose
extra term is numerically negligible passes while testing nothing.  GOSTSHYP is
unusually sharp here: it has no electrostatic component at all, so there is no
dominant term for a broken level-set route to hide behind.  The cavity response
is ~84% of ``dE/dP`` rather than a correction, and all three gradient channels
are comparable in size (1.3e-2, 1.4e-2, 1.5e-2 for water/STO-3G at 50 GPa).

One caveat the layering cannot fix: the finite differences are self-consistent,
comparing the model's derivative against the model's own energy.  A convention
error applied consistently to both would leave every one of them passing.  That
is what ``GOLDEN`` and ``test_conventions_matches_the_frozen_reference`` are
for, and why they cover the energy, amplitudes, adjoints and gradient rather
than the energy alone.
"""

import functools
import math
import os
from dataclasses import dataclass

import numpy as np
import pytest

try:
    from pyscf import dft, gto, scf
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from .interface import (
    DensityResponse,
    GostshypAmplitudeResponse,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    Structure,
)
from .parameters import IsodensityParameters
from .pyscf import (
    _D_CART_ORDER,
    _F_RHO2_FIRST_MOMENT,
    _S_NORM,
    _S_OVER_D_NORM,
    _S_OVER_F_NORM,
    _S_OVER_P_NORM,
    DROP,
    GPA_TO_AU,
    GaussianMoments,
    Isodensity,
    PySCFHost,
    PySCFSolvation,
    _int3c1e,
)

#: Every test in this module is a GOSTSHYP test; the second marker selects the
#: layer, and meson turns each into its own target.
pytestmark = pytest.mark.gostshyp

#: 50 GPa in Hartree / bohr^3
PRESSURE = 50.0 * GPA_TO_AU
#: Dielectric constant used by the component-composition check.
EPSILON = 80.0
#: Lebedev order
NLEB = 26
#: Cavity projection tolerance
PROJ_TOL = 1e-12

#: FD step on the density matrix, in units of a unit-Frobenius-norm direction
STEP_DM = 1e-4
#: FD step on nuclear coordinates in bohr; matches test_pyscf.py's tuned value
STEP_R = 2.5e-4
#: FD step on a grid point center, in bohr
STEP_C = 1e-4

#: Tolerances
REL_THR = 1e-9
ABS_THR = REL_THR / 10.0
#: Frozen-cavity tolerances
FROZEN_REL_THR = 1e-9
FROZEN_ABS_THR = FROZEN_REL_THR / 10.0
#: Surface-parameter derivatives tolerances
PARAM_REL_THR = 1e-9
PARAM_ABS_THR = PARAM_REL_THR / 10.0
#: Independent-quadrature comparison
QUAD_ATOL = 1e-10
QUAD_RTOL = 1e-9

#: Mirrors ``overlap_floor`` in ``src/moist/model/component/gostshyp.f90``
OVERLAP_FLOOR = 1.0e-9
#: A negative control must miss by at least this many tolerances
VACUITY_FACTOR = 100.0
#: An analytic quantity must be at least this large
MIN_SIGNAL = 1e-12


@dataclass(frozen=True)
class System:
    """A test solute.  Geometries are in Angstrom."""

    atom: str
    charge: int = 0


#: Test solutes
SYSTEMS = {
    # Slightly asymmetric water
    "water": System(
        """O  0.0000  0.0000 -0.3893
           H  0.7629  0.0000  0.1947
           H -0.7991  0.0953  0.2223"""
    ),
    # Glyciine zwitterion
    "glycine_zwitterion": System(
        """C  0.000  0.000  0.000
           C  1.540  0.000  0.000
           O  2.150  1.070  0.000
           O  2.150 -1.070  0.000
           N -0.760  0.000  1.280
           H -0.470  0.850 -0.490
           H -0.470 -0.850 -0.490
           H -1.760  0.000  1.100
           H -0.500  0.830  1.830
           H -0.500 -0.830  1.830"""
    ),
    # Anion
    "fluoroacetate": System(
        """C  0.000  0.000  0.000
           C  1.550  0.000  0.000
           F -0.640  1.240  0.000
           O  2.160  1.080  0.000
           O  2.160 -1.080  0.000
           H -0.390 -0.540  0.860
           H -0.390 -0.540 -0.860""",
        charge=-1,
    ),
}

#: Active (system, basis) cases
CASES = [
    ("water", "sto-3g"),                  # ~6 s   the workhorse
    ("water", "def2-svp"),                # ~6 s   covers l >= 2
    # ("fluoroacetate", "sto-3g"),         # ~21 s  anion
    ("fluoroacetate", "def2-svp"),         # ~24 s
    # ("glycine_zwitterion", "sto-3g"),   # ~21 s  neutral, charge-separated
    # ("glycine_zwitterion", "def2-svp"), # ~35 s
]
#: The case used by tests that pin a convention once rather than sweeping.
PRIMARY_CASE = CASES[0]

#: Reference valuesfrom the pure-Python implementation before creating the
#: standalone ``gostshyp`` component
GOLDEN = {
    ("water", "sto-3g"): {
        "energy": 0.1340185971628404,
        "ngrid": 71,
        "alpha_sum": 703.3145198328141,
        "alpha_norm": 301.9375922565026,
        "alpha_min": 0.0,
        "alpha_max": 246.68894895104054,
        "fock_frozen_sum": 0.002772468173207137,
        "fock_frozen_norm": 0.007890395554882459,
        "fock_frozen_min": -0.0013445039833400265,
        "fock_frozen_max": 0.0050487842334908314,
        "w_lsf0_sum": -0.07225582537940567,
        "w_lsf0_norm": 0.011984317074386375,
        "w_lsf0_min": -0.0035446644483471557,
        "w_lsf0_max": 0.0,
        "w_lsf1_sum": -0.00373430836922124,
        "w_lsf1_norm": 0.0060409051959967,
        "w_lsf1_min": -0.0014134523386366649,
        "w_lsf1_max": 0.0015364658891712886,
        "w_lsf2_sum": 0.0726577022287183,
        "w_lsf2_norm": 0.010664055525717617,
        "w_lsf2_min": -0.000934026258796762,
        "w_lsf2_max": 0.0021791011349693848,
        "gradient_sum": -5.70463815074973e-15,
        "gradient_norm": 0.027279652884437636,
        "gradient_min": -0.018309048197592433,
        "gradient_max": 0.016562186729215386,
    },
    ("water", "def2-svp"): {
        "energy": 0.15034462492176853,
        "ngrid": 71,
        "alpha_sum": 679.4410221485954,
        "alpha_norm": 302.2629419120036,
        "alpha_min": 0.0,
        "alpha_max": 249.91722445491104,
        "fock_frozen_sum": 0.20115384963534755,
        "fock_frozen_norm": 0.07009732976838826,
        "fock_frozen_min": -0.0059689710567790745,
        "fock_frozen_max": 0.0424777918638201,
        "w_lsf0_sum": -0.0692949109781968,
        "w_lsf0_norm": 0.012383979219126862,
        "w_lsf0_min": -0.004503088163928888,
        "w_lsf0_max": 0.0,
        "w_lsf1_sum": 0.0038287978871702553,
        "w_lsf1_norm": 0.006232115054651985,
        "w_lsf1_min": -0.002468816383539373,
        "w_lsf1_max": 0.002537238688237087,
        "w_lsf2_sum": 0.11295177621586769,
        "w_lsf2_norm": 0.0148415051860289,
        "w_lsf2_min": -0.0012398761281829382,
        "w_lsf2_max": 0.004416488076594353,
        "gradient_sum": 7.13925446538255e-15,
        "gradient_norm": 0.02023312701021516,
        "gradient_min": -0.012902438573779188,
        "gradient_max": 0.011492781673742245,
    },
}

#: The golden values are a different platform's arithmetic away from exact, but
#: the cavity projection is the only real variability and it converges to
#: PROJ_TOL.  Anything looser than this would stop catching a convention error.
GOLDEN_RTOL = 1e-9


CASE_PARAMS = [
    pytest.param(system, basis, id=f"{system}-{basis}") for system, basis in CASES
]


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


def cavity_config(rho_iso=None):
    """The isodensity DROP configuration the wall is built on."""
    lsf = Isodensity() if rho_iso is None else Isodensity(rho_iso=rho_iso)
    return DROP(lsf=lsf, nleb=NLEB, tolerance=PROJ_TOL)


def make_wall(mol, positions=None, *, dm, pressure=PRESSURE, rho_iso=None):
    """Host plus an evaluated GOSTSHYP driver ("wall") at ``dm``, optionally displaced."""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    wall = PySCFSolvation(mol, cavity_config(rho_iso), [ModelComponentGOSTSHYP(pressure)])
    wall.evaluate(dm)
    return wall.host, wall


def item_of(response, kind):
    """The item of class ``kind`` the walk over ``response`` meets, or None."""
    found = [item for item in response if isinstance(item, kind)]
    assert len(found) <= 1, found
    return found[0] if found else None


def fd_density(mol, dm, direction, *, pressure=PRESSURE):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        _host, wall = make_wall(mol, dm=dm + offset * STEP_DM * direction, pressure=pressure)
        samples.append(wall.energy)
        grids.append(wall.model.cavity.ngrid)
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
        grids.append(wall.model.cavity.ngrid)
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

    floor = OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
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
    picks = np.argsort(-np.abs(wall.moments.traces(dm)[1]))[:4]
    centers = wall.moments.centers[picks]
    omega = wall.moments.omega[picks]

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
    """``ftilde == n . grad_r gtilde``, i.e. minus the grid point-center gradient.

    Displacing the field point is the opposite of displacing the Gaussian
    center, which is the whole content of the sign in ``_build_integrals``.  Get
    it wrong and the wall pushes outward instead of inward.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    ftilde = wall.moments.traces(dm)[1]

    samples = []
    for offset in FD4_OFFSETS:
        centers = wall.moments.centers + offset * STEP_C * wall.moments.normals
        samples.append(wall.moments.traces(dm, centers=centers)[0])
    center_gradient = fd4(samples, STEP_C)

    assert array_deviation(ftilde, -center_gradient, thr_rel=FROZEN_REL_THR) <= 1.0


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

    gvfield = np.einsum("uvja,uv->ja", wall.moments.f_vector(), dm, optimize=True)
    projected = np.einsum("ja,ja->j", wall.moments.normals, gvfield, optimize=True)
    np.testing.assert_allclose(projected, wall.moments.traces(dm)[1], rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# params -- surface-parameter derivatives at a frozen cavity
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_first_moment_is_the_center_gradient(system, basis):
    """``dgtilde/dC_a == 2 omega Pt_a``, pinning the p moment component by component.

    ``Pt`` is the only moment the energy itself depends on -- through
    ``ftilde = -2 omega (n . Pt)`` -- so an error here moves every downstream
    quantity at once.  The Becke quadrature in the integrals layer pins the
    fakemol block; this pins the contraction that turns it into the moment the
    component is handed.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    _gt, pt, _mt, _rt = wall.moments._surface_moments(wall.moments._density_matrix_cart(dm))

    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = STEP_C
        samples = [
            wall.moments.traces(dm, centers=wall.moments.centers + offset * shift)[0]
            for offset in FD4_OFFSETS
        ]
        derivative = fd4(samples, STEP_C)
        expected = 2.0 * wall.moments.omega * pt[:, axis]
        assert array_deviation(expected, derivative, thr_rel=PARAM_REL_THR, thr_abs=PARAM_ABS_THR) <= 1.0, axis


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_moments_reach_the_component_intact(system, basis):
    """The energy moist reports is the one the host's own moments imply.

    The moments cross a language boundary and change layout on the way -- numpy
    ``(ngrid, 3)`` becomes Fortran ``(3, ngrid)``, ``(ngrid, 3, 3)`` becomes
    ``(3, 3, ngrid)``.  A transposed hand-off is the failure mode this refactor
    introduces and nothing else in this layer would see, so the whole exchange
    is closed here against the documented formula
    ``E = sum_j p a_j gtilde_j / ftilde_j``.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gt, ftilde = wall.moments.traces(dm)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = PRESSURE * wall.moments.areas * gt / ftilde
    # Grid points the component switched off carry a zero amplitude.
    active = item_of(wall.response, GostshypAmplitudeResponse).w_overlap != 0.0
    expected = float(np.sum(np.where(active, terms, 0.0)))

    assert wall.energy == pytest.approx(expected, rel=1e-12)


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_energy_is_insensitive_to_the_activity_floor(system, basis):
    """Moving the cut six decades must not move the energy.

    The floor exists for the *derivatives* -- ``beta ~ gtilde^2/ftilde^2``
    diverges as the overlap vanishes -- and it is only defensible if it costs
    the energy nothing, because it does not fall in a gap.  Water gives it about
    a decade of clearance, but that is the easy case: measured once on
    fluoroacetate/STO-3G, the smallest kept point sat at ``|ftilde|/max = 2.0e-9``
    and the largest dropped at ``6.7e-10`` -- a factor of three apart, 24 of 159
    points discarded.  A cut through a continuum is arbitrary, and what makes an
    arbitrary cut acceptable is exactly the insensitivity checked here.  Enable a
    larger row in ``CASES`` to run this against that harder case.

    So this is the test that fails if the floor ever becomes load-bearing.  If
    the energy starts depending on where the cut falls, the cut is discarding
    physics -- most likely a point approaching the genuine pole at
    ``ftilde = 0``, where the wall turns locally attractive -- and the model
    needs a treatment of that neighbourhood rather than a threshold.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gt, ftilde = wall.moments.live_traces
    relative = np.abs(ftilde) / np.abs(ftilde).max()
    # Masked explicitly rather than via nansum: a NaN among the *kept* points is
    # a real failure and must not be swallowed with the discarded ones.
    with np.errstate(divide="ignore", invalid="ignore"):
        contribution = PRESSURE * wall.moments.areas * gt / ftilde

    def energy_at(floor):
        keep = relative > floor
        assert np.all(np.isfinite(contribution[keep])), f"non-finite kept term at {floor:.0e}"
        return float(contribution[keep].sum())

    # The default floor reproduces what the component reported.
    assert energy_at(OVERLAP_FLOOR) == pytest.approx(wall.energy, rel=1e-12)

    # ...and six decades either side of it change nothing that matters.
    energies = [energy_at(floor) for floor in (1e-6, 1e-8, 1e-10, 1e-12)]
    for floor, energy in zip((1e-6, 1e-8, 1e-10, 1e-12), energies):
        assert energy == pytest.approx(wall.energy, rel=1e-4), f"floor={floor:.0e}"

    # Vacuous unless the floor is actually cutting something.
    assert wall.moments.inactive_count > 0
    assert wall.moments.inactive_count < wall.model.cavity.ngrid // 2


@pytest.mark.gostshyp_params
def test_params_inactive_count_does_not_follow_the_pressure():
    """The dropped-point count describes the surface, not the pressure.

    The component short-circuits at zero pressure and hands back amplitudes that
    are identically zero, so a count read off them calls the whole grid dropped
    -- the trap :meth:`effective_volume` already documents one screen above.
    This number exists to be watched across a trajectory for a moment supply
    going wrong, which it cannot do if it also tracks the applied pressure.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, live = make_wall(mol, dm=dm)
    _host, quiet = make_wall(mol, dm=dm, pressure=0.0)

    # The cavity is isodensity, so dropping the pressure must not move a point.
    assert quiet.model.cavity.ngrid == live.model.cavity.ngrid

    # Vacuous unless the floor is actually cutting something.
    assert 0 < live.moments.inactive_count < live.model.cavity.ngrid // 2
    assert quiet.moments.inactive_count == live.moments.inactive_count


@pytest.mark.gostshyp_params
def test_params_update_drops_the_previous_surfaces_moments():
    """A cavity update invalidates the moments, and reading without new ones fails.

    moist cannot rebuild them itself -- they are AO-basis integrals only the host
    can form -- and their shapes cannot betray a stale set, because a geometry
    step normally preserves the point count.  So the update drops them, turning
    "forgot to re-supply" into the same error as "never supplied" instead of a
    plausible energy for a surface that has moved.

    :meth:`PySCFSolvation.evaluate` always re-supplies in the same breath, so
    this only bites a host driving the model by hand.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    before = wall.energy
    original_xyz = wall.model.cavity.xyz.copy()
    original_ngrid = wall.model.cavity.ngrid

    displaced = mol.atom_coords()
    displaced[0, 2] += 0.05
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)
    wall.model.update(Structure(numbers, displaced))

    # Vacuous unless the surface really moved while keeping its point count --
    # exactly the case a shape check cannot distinguish.
    assert wall.model.cavity.ngrid == original_ngrid
    assert not np.allclose(wall.model.cavity.xyz, original_xyz)

    coupling = wall.model.new_coupling()
    wall.model.prepare_energy(coupling)
    energy = np.array(0.0)
    with pytest.raises(RuntimeError, match="missing required outputs"):
        wall.model.get_energy(coupling, energy)

    # Going through the wrapper re-supplies, and the energy tracks the geometry.
    _host2, moved = make_wall(mol, displaced, dm=dm)
    assert moved.energy != pytest.approx(before, rel=1e-9)


@pytest.mark.gostshyp_params
def test_params_second_moment_matches_a_second_difference():
    """``Mt_ab`` against a second difference of ``gtilde`` in the grid point center.

    With ``G = exp(-omega rho^2)``, ``d^2 gtilde/dC_a dC_b = 4 omega^2 Mt_ab -
    2 omega delta_ab gtilde``.  This isolates the d-fakemol constant and the
    ``_D_CART_ORDER`` mapping, which ``test_params_derivatives_match_fd`` only
    sees through the contracted combination ``n . Mt``.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    dm_cart = wall.moments._density_matrix_cart(dm)

    _gt, _pt, moment, _rt = wall.moments._surface_moments(dm_cart)
    gtilde = wall.moments.traces(dm)[0]
    omega = wall.moments.omega
    step = 2.0e-3

    def gt_at(shift):
        return wall.moments.traces(dm, centers=wall.moments.centers + shift)[0]

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
def test_gostshyp_is_a_composable_model_component():
    """One PySCF coupling completes GOSTSHYP alone or alongside CPCM."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)

    def evaluate(components):
        solvation = PySCFSolvation(mol, cavity_config(), components)
        solvation.evaluate(dm)
        return solvation

    gostshyp = evaluate([ModelComponentGOSTSHYP(PRESSURE)])
    cpcm = evaluate([ModelComponentCPCM(EPSILON)])
    combined = evaluate([ModelComponentCPCM(EPSILON), ModelComponentGOSTSHYP(PRESSURE)])

    assert combined.energy == pytest.approx(cpcm.energy + gostshyp.energy, rel=1e-12)
    np.testing.assert_allclose(combined.fock, cpcm.fock + gostshyp.fock, rtol=1e-11)
    np.testing.assert_allclose(
        combined.gradient(dm),
        cpcm.gradient(dm) + gostshyp.gradient(dm),
        rtol=1e-10,
        atol=1e-11,
    )


@pytest.mark.gostshyp_fock
def test_results_are_frozen_and_pinned_to_the_evaluated_density():
    """Results are immutable values, and a mutated host density cannot leak in."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host, wall = make_wall(mol, dm=dm)

    result = wall.result
    expected_gradient = wall.gradient(dm)
    host.dm = np.zeros_like(dm)

    assert result.energy == wall.energy
    assert wall.evaluate(dm) is result
    np.testing.assert_allclose(wall.gradient(dm), expected_gradient)
    assert not result.fock.flags.writeable


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
    assert wall.moments.effective_volume() == pytest.approx(wall.energy / PRESSURE, rel=1e-12)
    # Not the cavity volume: the wall is weighted by the local density overlap.
    assert wall.moments.effective_volume() != pytest.approx(wall.model.cavity.volume, rel=1e-3)

    # The three walls share a density and therefore a surface, so the volume is
    # common to them.  At p = 0 the energy vanishes with the pressure while the
    # volume does not, which is why it cannot be recovered as E/p there.
    assert vacuum.moments.effective_volume() == pytest.approx(wall.moments.effective_volume(), rel=1e-12)
    assert doubled.moments.effective_volume() == pytest.approx(wall.moments.effective_volume(), rel=1e-12)


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
                wall.moments.centers,
                wall.moments.areas,
                wall.moments.omega,
                wall.moments.normals,
            )
            for offset in FD4_OFFSETS
        ]
        numerical = fd4(samples, STEP_DM)
        frozen = wall.frozen_fock()
        analytic = float(np.einsum("uv,uv->", frozen, direction))
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
    fock = wall.fock

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
    """AO centers move,  grid points frozen: the ``int3c1e_ip1`` term alone.

    Pins the sign of the integral derivative, the both-legs factor of two, and
    the mapping from cartesian AO rows onto atoms.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)
    analytic = wall.gradient_channels(dm)["moments"].flatten(order="C")

    for index in sampled_coordinates(mol.natm):
        samples = []
        for offset in FD4_OFFSETS:
            displaced = positions.copy()
            displaced.flat[index] += offset * STEP_R
            moved = mol.set_geom_(displaced, unit="Bohr", inplace=False)
            samples.append(
                frozen_surface_energy(
                    moved, dm, wall.moments.centers, wall.moments.areas, wall.moments.omega, wall.moments.normals
                )
            )
        numerical = fd4(samples, STEP_R)
        assert deviation(
            analytic[index], numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR
        ) <= 1.0


@pytest.mark.gostshyp_gradient
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_gradient_surface_channel_matches_fd(system, basis):
    """The cavity's own response at a frozen level-set field.

    Displacing only the structure handed to the cavity -- never the molecule the
    level set is built from -- drags the atom-anchored reference grid while the
    density stays put.  That is exactly the channel moist contracts in reverse
    mode, so this closes the component's ``w_a``/``w_xyz``/``w_normal`` against a
    finite difference of the energy through moist's own cavity derivatives:
    the area route including its switching factor, and the position route
    including the normal's point-motion fold.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host, wall = make_wall(mol, dm=dm)
    analytic = wall.gradient_channels(dm)["model"].flatten(order="C")
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)

    def anchored_energy(displaced):
        # The cavity moves; host.mol -- and hence the level set -- does not.
        cavity = cavity_config().build(source=host)
        cavity.update(Structure(numbers, displaced))
        result = cavity.cavity
        centers = np.ascontiguousarray(result.xyz)
        areas = np.ascontiguousarray(result.a)
        normals = np.array(result.normal0, order="C", copy=True)
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
    analytic = wall.gradient(dm).flatten(order="C")

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
    gradient = wall.gradient(dm)

    assert np.abs(gradient).max() > MIN_SIGNAL
    np.testing.assert_allclose(gradient.sum(axis=0), 0.0, atol=1e-10)


@pytest.mark.gostshyp_gradient
def test_gradient_is_exactly_the_sum_of_its_channels():
    """No channel is silently dropped inside :meth:`nuclear_gradient`."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    channels = (
        wall.gradient_channels(dm)["moments"]
        + wall.gradient_channels(dm)["density"]
        + wall.gradient_channels(dm)["model"]
    )
    np.testing.assert_allclose(wall.gradient(dm), channels, rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# conventions -- channel semantics and the negative controls
# ----------------------------------------------------------------------


# The surface weights themselves are no longer visible from Python: they are
# built and consumed inside the moist ``gostshyp`` component.  What used to be
# checked here -- ``w_f == 0`` exactly, and a per-channel starve control for
# ``w_normal``/``w_a`` -- now lives in test/unit/test_model/component/gostshyp.f90,
# where ``check_surface_weights`` finite-differences each channel separately
# rather than only detecting that some channel is missing.


@pytest.mark.gostshyp_conventions
def test_conventions_fock_requires_cavity_response():
    """The eq-16 Fock alone is not ``dE/dP`` for a surface that follows the density.

    With no electrostatic component there is nothing to dilute the error: the
    cavity response is the majority of the derivative, not a correction.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    response = wall.fock - wall.frozen_fock()
    assert np.abs(response).max() > MIN_SIGNAL

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction)
    starved = float(np.einsum("uv,uv->", wall.frozen_fock(), direction))
    assert deviation(starved, numerical) > VACUITY_FACTOR


@pytest.mark.gostshyp_conventions
def test_conventions_amplitudes_are_both_needed():
    """Both host amplitudes carry a non-negligible share of the frozen Fock.

    ``w_overlap`` and ``w_normal_deriv`` arrive as a pair with their signs already
    folded in, so the natural failure is using one and dropping the other, or
    flipping the fold.  Either would leave a Fock that still looks plausible.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    amplitude = item_of(wall.response, GostshypAmplitudeResponse)

    assert np.abs(amplitude.w_overlap).max() > MIN_SIGNAL
    assert np.abs(amplitude.w_normal_deriv).max() > MIN_SIGNAL
    # The fold is a sign, not an absolute value: the two channels oppose.
    assert amplitude.w_overlap.max() > 0.0
    assert amplitude.w_normal_deriv.min() < 0.0

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd4(
        [
            frozen_surface_energy(
                mol,
                dm + offset * STEP_DM * direction,
                wall.moments.centers,
                wall.moments.areas,
                wall.moments.omega,
                wall.moments.normals,
            )
            for offset in FD4_OFFSETS
        ],
        STEP_DM,
    )

    g_only = np.einsum(
        "j,uvj->uv", amplitude.w_overlap, wall.moments._G, optimize=True
    )
    starved = float(np.einsum("uv,uv->", 0.5 * (g_only + g_only.T), direction))
    assert deviation(starved, numerical, thr_rel=FROZEN_REL_THR) > VACUITY_FACTOR
    # ...and the pair together is the quantity that does close.
    full = float(np.einsum("uv,uv->", wall.frozen_fock(), direction))
    assert deviation(full, numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR) <= 1.0


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize("dropped", ["integral", "field", "surface"])
def test_conventions_gradient_requires_every_channel(dropped):
    """Each of the three nuclear routes carries a non-negligible share.

    Two are the host's and one is moist's, so this also pins the split itself:
    dropping the moist term must break the gradient, or the component is not
    actually contributing what its surface weights claim.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)

    channels = {
        "integral": wall.gradient_channels(dm)["moments"],
        "field": wall.gradient_channels(dm)["density"],
        "surface": wall.gradient_channels(dm)["model"],
    }
    assert np.abs(channels[dropped]).max() > MIN_SIGNAL
    starved = sum(value for name, value in channels.items() if name != dropped)

    index = int(np.argmax(np.abs(channels[dropped].flatten(order="C"))))
    numerical = fd_position(mol, positions, dm, index)
    assert deviation(starved.flatten(order="C")[index], numerical) > VACUITY_FACTOR


# Two controls retired with the port, both because the quantity they starved is
# no longer reachable from Python:
#
#   * the area route must use the true ``da_i/dR_A`` rather than the Gaussian-
#     width proxy -- now moist's internal choice, and covered end to end by
#     ``test_gradient_surface_channel_matches_fd``;
#   * the activity floor must gate the amplitudes rather than only ``dE/da``.
#     That is now structural: the component computes the mask once, beside the
#     amplitudes, and ``test_gostshyp_energy_matches_amplitudes`` in
#     test/unit/test_model/component/gostshyp.f90 checks the energy and the host
#     amplitudes agree, which is exactly the property a split mask breaks.


@pytest.mark.gostshyp_conventions
def test_conventions_anchor_area_derivatives_sum_to_the_total():
    """``sum_i a_i1_rA == A_tot1_rA``, pinning the anchor buffer layout.

    The binding exposes C arrays with the native axes reversed.  This identity is independent of
    any GOSTSHYP physics, so it bisects a layout error away from a weight error.

    Driven on a standalone cavity: the wall's cavity now belongs to the model,
    and this checks the binding rather than the wall.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = PySCFHost(mol)
    host.dm = dm
    cavity = cavity_config().build(source=host)
    cavity.update(host.structure())

    cavity.compute_anchor_gradient()
    anchor = cavity.get_anchor_gradient()
    np.testing.assert_allclose(
        anchor.a_i1_rA.sum(axis=0), anchor.A_tot1_rA, rtol=0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        anchor.v_i1_rA.sum(axis=0), anchor.V_tot1_rA, rtol=0.0, atol=1e-12
    )


@pytest.mark.gostshyp_conventions
def test_conventions_failed_evaluation_leaves_no_stale_results(monkeypatch):
    """An evaluation that raises must not leave the previous density readable.

    :meth:`evaluate` moves the host density and the cavity before it can fail --
    a density with no isodensity surface, a projection that will not converge --
    so results cached from the previous call describe a surface that no longer
    exists.  The failure is forced here rather than hunted for: what is being
    pinned is the invalidation, not any particular way of tripping it.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # Vacuous unless there is a real result to go stale.
    assert wall.moments.effective_volume() > 0.0
    assert item_of(wall.response, GostshypAmplitudeResponse) is not None

    def boom(*_args, **_kwargs):
        raise RuntimeError("integral build failed")

    monkeypatch.setattr(wall.moments, "_build_integrals", boom)
    with pytest.raises(RuntimeError, match="integral build failed"):
        wall.evaluate(1.001 * dm)

    assert wall.result is None
    for read in (lambda: wall.energy, lambda: wall.fock, lambda: wall.response):
        with pytest.raises(RuntimeError, match="evaluate"):
            read()
    # ...and the gradient re-evaluates rather than reading anything stale.
    with pytest.raises(RuntimeError, match="integral build failed"):
        wall.gradient(dm)


@pytest.mark.gostshyp_conventions
def test_conventions_failed_evaluation_publishes_nothing_partial(monkeypatch):
    """A failure in the *last* step of :meth:`evaluate` publishes nothing either.

    The cavity and the moments are built by then, so the energy waiting behind
    that call is a perfectly good number -- which is why this is worth pinning
    rather than obvious: it is the tempting one to keep.  :meth:`evaluate`
    raised, so no caller was ever handed the state it belongs to, and a
    half-published result is one no accessor can tell from a whole one.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    def boom(*_args, **_kwargs):
        raise RuntimeError("response assembly failed")

    monkeypatch.setattr(wall.model, "get_response", boom)
    with pytest.raises(RuntimeError, match="response assembly failed"):
        wall.evaluate(1.001 * dm)

    assert wall.result is None
    for read in (lambda: wall.energy, lambda: wall.fock, lambda: wall.response):
        with pytest.raises(RuntimeError, match="evaluate"):
            read()


@pytest.mark.gostshyp_conventions
def test_conventions_moments_are_built_once_per_evaluation(monkeypatch):
    """Build each required moment block once, with no rebuild on result reads."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)

    original = GaussianMoments._surface_moments
    calls = []

    def counted(self, *args, **kwargs):
        calls.append(frozenset(kwargs["required"]))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(GaussianMoments, "_surface_moments", counted)

    _host, wall = make_wall(mol, dm=dm)
    assert calls == [frozenset({"gt", "pt"}), frozenset({"mt", "rt"})]

    fock = wall.fock
    gradient = wall.gradient(dm)
    volume = wall.moments.effective_volume()
    assert len(calls) == 2, "a result read rebuilt the moments"
    assert np.abs(fock).max() > MIN_SIGNAL
    assert np.abs(gradient).max() > MIN_SIGNAL
    assert volume > 0.0

    # The cached traces are the ones a fresh build would produce.
    gt, ftilde = wall.moments.live_traces
    rebuilt_gt, rebuilt_ftilde = wall.moments.traces(dm)
    np.testing.assert_allclose(gt, rebuilt_gt, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(ftilde, rebuilt_ftilde, rtol=0.0, atol=0.0)


def _golden_summary(name, array):
    """Reduce an array the way ``scratchpad/capture_golden.py`` did."""
    flat = np.asarray(array, dtype=np.float64).ravel()
    return {
        f"{name}_sum": float(flat.sum()),
        f"{name}_norm": float(np.linalg.norm(flat)),
        f"{name}_min": float(flat.min()),
        f"{name}_max": float(flat.max()),
    }


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize(
    "system,basis",
    [pytest.param(*case, id=f"{case[0]}-{case[1]}") for case in GOLDEN],
)
def test_conventions_matches_the_frozen_reference(system, basis):
    """Every reported quantity still equals its pre-port value.

    The finite differences elsewhere in this file are self-consistent: they
    check the model's derivative against the model's own energy, so a
    convention error applied consistently to both -- a wrong angular constant,
    a mis-ordered cartesian component -- leaves all of them passing.  This is
    the only test that would notice, which is why it covers the energy, the
    amplitudes, the level-set adjoints and the gradient rather than just the
    energy.

    Parametrised over ``GOLDEN``, *not* ``CASES``: these values were captured
    from the pure-Python implementation that existed before the physics moved
    into the Fortran component, so they exist only for the cases that were
    enabled at that time.  Enabling a new row in ``CASES`` must not fail here.

    Nor should a new row simply be added to ``GOLDEN``.  That pre-port
    implementation is gone, so values captured for a new case would come from
    the code under test -- a regression pin against future drift, which is
    useful, but not the independent cross-check the water rows carry.  If you
    add one, say in a comment which of the two it is.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    # The independent reference was recorded with the former PySCF contour.
    # Keep that physical input explicit; do not regenerate the reference values.
    _host, wall = make_wall(mol, dm=dm, rho_iso=4e-4)

    amplitude = item_of(wall.response, GostshypAmplitudeResponse)
    weights = item_of(wall.response, DensityResponse)
    actual = {"energy": float(wall.energy), "ngrid": int(wall.model.cavity.ngrid)}
    actual.update(_golden_summary("alpha", amplitude.w_overlap))
    actual.update(
        _golden_summary("fock_frozen", wall.frozen_fock())
    )
    # The frozen reference predates the `dS/drho` fold, so it holds adjoints of
    # the level set.  moist now returns adjoints of the density, larger by
    # exactly `-scale`; dividing it back out compares the same quantity the
    # pre-port implementation reported rather than renumbering the reference.
    def level_set(array):
        return np.asarray(array) / -IsodensityParameters().scale

    actual.update(_golden_summary("w_lsf0", level_set(weights.w_rho)))
    actual.update(_golden_summary("w_lsf1", level_set(weights.w_grad_rho)))
    actual.update(_golden_summary("w_lsf2", level_set(weights.w_hess_rho)))
    actual.update(_golden_summary("gradient", wall.gradient(dm)))

    expected = GOLDEN[(system, basis)]
    assert actual.keys() == expected.keys()
    for key, reference in expected.items():
        if key == "ngrid":
            assert actual[key] == reference, key
        elif abs(reference) < MIN_SIGNAL:
            # Translational invariance of the gradient, and the sign-definite
            # channels whose extremum is exactly zero.
            assert abs(actual[key]) < MIN_SIGNAL, key
        else:
            assert actual[key] == pytest.approx(reference, rel=GOLDEN_RTOL), key
