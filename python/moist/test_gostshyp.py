"""End-to-end tests of the GOSTSHYP pressure model on a moist isodensity cavity

GOSTSHYP is mostly QM-side integral work; moist contributes the cavity response
- every analytic quantity is checked against a finite difference of the energy
  the model reports, at 50 GPa

Physics lives in moist's ``gostshyp`` component; this file tests the host half
and the assembled result
- component surface weights are finite-differenced channel by channel in
  ``test/unit/test_model/component/gostshyp.f90``, against a model density with
  closed-form Gaussian moments
- where a control moved there, the comment at its old site says so

Layers, so a failure localises:

``integrals``
    libcint fakemol constants, cartesian d/f component orders, identity
    ``f == n . grad_r g``
    - moist only places the grid points; a failure is integral bookkeeping
    - the constants are the one convention that stayed host-side: load-bearing
``params``
    Gaussian moments this module hands the component, and the hand-off itself
    - a transposed moment is the failure mode the language boundary introduces
``fock``
    Energy and ``dE/dP`` with the surface following the density
    - the level-set chain rule
``gradient``
    ``dE/dR`` at fixed density, split into integral, field and surface channels,
    then assembled
    - integral and field are the host's; surface is moist contracting the
      component's weights in reverse mode
``conventions``
    Channel semantics, negative controls, frozen reference values

Each layer carries a negative control: an FD test whose extra term is
numerically negligible passes while testing nothing
- GOSTSHYP has no electrostatic component, so no dominant term hides a broken
  level-set route
- cavity response is ~84% of ``dE/dP``, not a correction
- the three gradient channels are comparable (1.3e-2, 1.4e-2, 1.5e-2 for
  water/STO-3G at 50 GPa)

Caveat the layering cannot fix: the FDs are self-consistent (model derivative
vs the model's own energy)
- a convention error applied to both leaves every one of them passing
- hence ``GOLDEN`` and ``test_conventions_matches_the_frozen_reference``,
  covering energy, amplitudes, adjoints and gradient, not the energy alone
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

from .gostshyp import (
    _D_CART_ORDER,
    _F_RHO2_FIRST_MOMENT,
    _S_NORM,
    _S_OVER_D_NORM,
    _S_OVER_F_NORM,
    _S_OVER_P_NORM,
    GPA_TO_AU,
    GostshypModel,
    GostshypWall,
    _int3c1e,
)
from .interface import (
    CavityDROPIsodensityCallback,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    SolvationModel,
    Structure,
)
from .pyscf import PySCFHost

#: Every test is a GOSTSHYP test
#: - second marker selects the layer, meson makes each a target
pytestmark = pytest.mark.gostshyp

#: 50 GPa in Hartree / bohr^3
PRESSURE = 50.0 * GPA_TO_AU
#: Dielectric constant of the component-composition check
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
    """Test solute, geometry in Angstrom"""

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
    # Glycine zwitterion
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
#: Case for tests that pin a convention once, no sweep
PRIMARY_CASE = CASES[0]

#: Reference values from the pure-Python implementation, before the standalone
#: ``gostshyp`` component existed
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

#: Golden relative tolerance
#: - a different platform's arithmetic from exact; the cavity projection is the
#:   only real variability and converges to PROJ_TOL
#: - looser would stop catching a convention error
GOLDEN_RTOL = 1e-9


CASE_PARAMS = [
    pytest.param(system, basis, id=f"{system}-{basis}") for system, basis in CASES
]


def deviation(actual, reference, *, thr_abs=None, thr_rel=None) -> float:
    """Deviation in units of the tolerance, ``<= 1`` passes

    Thresholds combine as test-drive's ``check(..., thr_abs=, thr_rel=)`` does:
    ``max(thr_abs, thr_rel*|ref|)``
    - absolute floor covers references near zero, relative bound scales with magnitude
    - a ratio makes a failure say how many tolerances were missed
    """
    thr_abs = ABS_THR if thr_abs is None else thr_abs
    thr_rel = REL_THR if thr_rel is None else thr_rel
    return abs(actual - reference) / max(thr_abs, thr_rel * abs(reference))


def array_deviation(actual, reference, *, thr_rel, thr_abs=None) -> float:
    """Worst elementwise deviation, scaled by the *array's* magnitude

    - per-element relative errors are meaningless when entries span orders of
      magnitude (small ones are differences of large ones)
    - tolerance is set from the largest entry of the reference
    """
    actual = np.asarray(actual)
    reference = np.asarray(reference)
    scale = float(np.max(np.abs(reference), initial=0.0))
    threshold = max(thr_abs if thr_abs is not None else 0.0, thr_rel * scale)
    return float(np.max(np.abs(actual - reference), initial=0.0)) / threshold


FD4_OFFSETS = (2, 1, -1, -2)


def fd4(values, step: float) -> float:
    """4-point central difference from samples at ``(+2, +1, -1, -2) * step``"""
    fpp, fp, fm, fmm = values
    return (-fpp + 8.0 * fp - 8.0 * fm + fmm) / (12.0 * step)


def fd2(plus, minus, step):
    """2-point central difference for the array-valued references"""
    return (np.asarray(plus) - np.asarray(minus)) / (2.0 * step)


# ----------------------------------------------------------------------
# System construction (cached: every parametrised test reuses them)


@functools.lru_cache(maxsize=None)
def molecule(system: str, basis: str):
    spec = SYSTEMS[system]
    return gto.M(
        atom=spec.atom, basis=basis, charge=spec.charge, unit="Angstrom", verbose=0
    )


@functools.lru_cache(maxsize=None)
def reference_density(system: str, basis: str):
    """Converged gas-phase RHF density, held fixed as a free parameter"""
    mean_field = scf.RHF(molecule(system, basis))
    mean_field.conv_tol = 1e-12
    mean_field.kernel()
    assert mean_field.converged, f"gas-phase SCF failed for {system}/{basis}"
    return mean_field.make_rdm1()


def make_wall(mol, positions=None, *, dm, pressure=PRESSURE):
    """Build a host and an updated GOSTSHYP wall at ``dm``, optionally displaced"""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    host = PySCFHost(mol)
    host.dm = dm
    with pytest.deprecated_call(match="Gostshyp"):
        wall = GostshypModel(host, pressure, nleb=NLEB, tolerance=PROJ_TOL)
    wall.update(dm)
    return host, wall


def fd_density(mol, dm, direction, *, pressure=PRESSURE):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample"""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        _host, wall = make_wall(mol, dm=dm + offset * STEP_DM * direction, pressure=pressure)
        samples.append(wall.energy)
        grids.append(wall.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, pressure=PRESSURE):
    """dE/dR along one cartesian coordinate at fixed density matrix"""
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
    """Unit-norm symmetric density-matrix perturbations"""
    rng = np.random.default_rng(seed)
    for _ in range(count):
        direction = rng.standard_normal((nao, nao))
        direction = 0.5 * (direction + direction.T)
        yield direction / np.linalg.norm(direction)


def sampled_coordinates(natm):
    """Cartesian coordinates to difference, all for a small solute"""
    total = 3 * natm
    if total <= 9:
        return list(range(total))
    return [int(i) for i in np.linspace(0, total - 1, 3, dtype=int)]


def frozen_surface_energy(mol, dm, centers, areas, omega, normals, pressure=PRESSURE):
    """GOSTSHYP energy for an explicit surface, independent of any wall state

    Re-derived from :func:`~moist.gostshyp._int3c1e`, not the model, so the finite
    differences using it cannot inherit a bug from the code under test
    """
    ngrid = int(np.asarray(centers).shape[0])
    g_cart = _int3c1e(mol, centers, omega, 0)
    p_cart = _int3c1e(mol, centers, omega, 1)
    ncart = g_cart.shape[0]
    p_cart = p_cart.reshape(ncart, ncart, ngrid, 3)

    c2s = np.eye(mol.nao_nr()) if mol.cart else np.asarray(mol.cart2sph_coeff(normalized="sp"))
    # BLAS leaves dirty FP status flags for these shapes; see _density_matrix_cart
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
# integrals -- libcint conventions, moist only places the grid points


@pytest.mark.gostshyp_integrals
def test_integrals_fakemol_normalization_constants():
    """Four angular constants and cartesian d/f orders, vs quadrature

    libcint gives a coefficient-1 shell of each angular momentum its own normalization
    - restoring the ratios ``N_s/N_l`` puts every moment in the units of ``g`` and
      makes ``f = n . grad g`` exact, not exact up to a factor
    - independent Becke-grid quadrature of the same monomial moments pins constants
      and component orders together
    - a transposed ``_D_CART_ORDER`` entry or mis-assigned ``_F_RHO2_FIRST_MOMENT``
      triple is invisible to any self-consistency check
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # Few well-conditioned grid points, cheap grid
    picks = np.argsort(-np.abs(wall.traces(dm)[1]))[:4]
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
    """``ftilde == n . grad_r gtilde``, i.e. minus the grid point-center gradient

    - displacing the field point is the opposite of displacing the Gaussian center:
      the sign in ``_build_integrals``
    - wrong sign: the wall pushes outward instead of inward
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    ftilde = wall.traces(dm)[1]

    samples = []
    for offset in FD4_OFFSETS:
        centers = wall.centers + offset * STEP_C * wall.normals
        samples.append(wall.traces(dm, centers=centers)[0])
    center_gradient = fd4(samples, STEP_C)

    assert array_deviation(ftilde, -center_gradient, thr_rel=FROZEN_REL_THR) <= 1.0


@pytest.mark.gostshyp_integrals
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_integrals_ftilde_is_the_normal_projection_of_the_f_vector(system, basis):
    """``ftilde_j == n_j . gvfield_j`` exactly

    - ``gvfield`` is the unprojected f-vector and doubles as ``dftilde/dn``
    - must agree to round-off, else the normal weight is built from a different
      quantity than the one it differentiates
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gvfield = np.einsum("uvja,uv->ja", wall.f_vector(), dm, optimize=True)
    projected = np.einsum("ja,ja->j", wall.normals, gvfield, optimize=True)
    np.testing.assert_allclose(projected, wall.traces(dm)[1], rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# params -- surface-parameter derivatives at a frozen cavity
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_first_moment_is_the_center_gradient(system, basis):
    """``dgtilde/dC_a == 2 omega Pt_a``, pinning the p moment component by component

    - ``Pt`` is the only moment the energy depends on, via
      ``ftilde = -2 omega (n . Pt)``: an error here moves every downstream quantity
    - Becke quadrature in the integrals layer pins the fakemol block; this pins the
      contraction into the moment handed to the component
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    _gt, pt, _mt, _rt = wall._surface_moments(wall._density_matrix_cart(dm))

    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = STEP_C
        samples = [
            wall.traces(dm, centers=wall.centers + offset * shift)[0]
            for offset in FD4_OFFSETS
        ]
        derivative = fd4(samples, STEP_C)
        expected = 2.0 * wall.omega * pt[:, axis]
        assert array_deviation(expected, derivative, thr_rel=PARAM_REL_THR, thr_abs=PARAM_ABS_THR) <= 1.0, axis


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_moments_reach_the_component_intact(system, basis):
    """Reported energy equals the one the host's own moments imply

    - moments cross a language boundary and change layout: numpy ``(ngrid, 3)``
      becomes Fortran ``(3, ngrid)``, ``(ngrid, 3, 3)`` becomes ``(3, 3, ngrid)``
    - a transposed hand-off is this refactor's failure mode; nothing else in this
      layer would see it
    - whole exchange closed against the documented formula
      ``E = sum_j p a_j gtilde_j / ftilde_j``
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gt, ftilde = wall.traces(dm)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = PRESSURE * wall.areas * gt / ftilde
    # Switched-off grid points carry zero amplitude
    active = wall.amplitudes != 0.0
    expected = float(np.sum(np.where(active, terms, 0.0)))

    assert wall.energy == pytest.approx(expected, rel=1e-12)


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_energy_is_insensitive_to_the_activity_floor(system, basis):
    """Moving the cut six decades must not move the energy

    Floor exists for the *derivatives*: ``beta ~ gtilde^2/ftilde^2`` diverges as the
    overlap vanishes
    - defensible only if it costs the energy nothing, as it does not fall in a gap
    - water: about a decade of clearance, but the easy case
    - fluoroacetate/STO-3G, measured once: smallest kept point
      ``|ftilde|/max = 2.0e-9``, largest dropped ``6.7e-10``, factor three apart,
      24 of 159 points discarded
    - a cut through a continuum is arbitrary; the insensitivity checked here is what
      makes it acceptable
    - enable a larger row in ``CASES`` to run the harder case

    Fails if the floor ever becomes load-bearing
    - then the cut discards physics, most likely a point nearing the genuine pole at
      ``ftilde = 0``, where the wall turns locally attractive
    - the model then needs a treatment of that neighbourhood, not a threshold
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gt, ftilde = wall.live_traces
    relative = np.abs(ftilde) / np.abs(ftilde).max()
    # Explicit mask, not nansum: a NaN among the *kept* points is a real failure
    with np.errstate(divide="ignore", invalid="ignore"):
        contribution = PRESSURE * wall.areas * gt / ftilde

    def energy_at(floor):
        keep = relative > floor
        assert np.all(np.isfinite(contribution[keep])), f"non-finite kept term at {floor:.0e}"
        return float(contribution[keep].sum())

    # Default floor reproduces the component's energy
    assert energy_at(OVERLAP_FLOOR) == pytest.approx(wall.energy, rel=1e-12)

    # Six decades either side change nothing that matters
    energies = [energy_at(floor) for floor in (1e-6, 1e-8, 1e-10, 1e-12)]
    for floor, energy in zip((1e-6, 1e-8, 1e-10, 1e-12), energies):
        assert energy == pytest.approx(wall.energy, rel=1e-4), f"floor={floor:.0e}"

    # Vacuous unless the floor cuts something
    assert wall.inactive_count > 0
    assert wall.inactive_count < wall.ngrid // 2


@pytest.mark.gostshyp_params
def test_params_inactive_count_does_not_follow_the_pressure():
    """Dropped-point count describes the surface, not the pressure

    - component short-circuits at zero pressure and returns zero amplitudes, so a
      count read off them calls the whole grid dropped (trap
      :meth:`effective_volume` documents one screen above)
    - the number is watched across a trajectory for a moment supply going wrong; it
      cannot be if it also tracks the pressure
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, live = make_wall(mol, dm=dm)
    _host, quiet = make_wall(mol, dm=dm, pressure=0.0)

    # Isodensity cavity: zero pressure moves no point
    assert quiet.ngrid == live.ngrid

    # Vacuous unless the floor cuts something
    assert 0 < live.inactive_count < live.ngrid // 2
    assert quiet.inactive_count == live.inactive_count


@pytest.mark.gostshyp_params
def test_params_supply_rejects_mis_shaped_moments():
    """Wrong-shaped moments are refused, not reinterpreted

    - binding hands moist four raw pointers; an accepted wrong shape would be read
      as whatever the grid size implies
    - cheap to pin, the only guard between a host bug and a plausible wrong number
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    gt, pt, mt, rt = wall._surface_moments(wall._density_matrix_cart(dm))
    good = (gt, np.asfortranarray(pt.T), np.asfortranarray(mt.transpose(1, 2, 0)),
            np.asfortranarray(rt.T))

    # Untransposed vector moments: the natural mistake
    with pytest.raises(ValueError, match="pt"):
        wall.model.supply_gostshyp(gt, pt, good[2], good[3])
    with pytest.raises(ValueError, match="mt"):
        wall.model.supply_gostshyp(gt, good[1], mt, good[3])
    with pytest.raises(ValueError, match="rt"):
        wall.model.supply_gostshyp(gt, good[1], good[2], rt)
    # Grid size no longer matching the live cavity
    with pytest.raises(RuntimeError):
        wall.model.supply_gostshyp(gt[:-1], good[1][:, :-1], good[2][:, :, :-1], good[3][:, :-1])


@pytest.mark.gostshyp_params
def test_params_update_drops_the_previous_surfaces_moments():
    """Cavity update invalidates the moments; reading without new ones fails

    - moist cannot rebuild them: AO-basis integrals only the host can form
    - shapes cannot betray a stale set, a geometry step normally preserves the
      point count
    - dropping them makes "forgot to re-supply" the same error as "never supplied",
      not a plausible energy for a moved surface
    - :meth:`GostshypModel.update` always re-supplies, so this only bites a host
      driving the model by hand
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    before = wall.energy
    original_xyz = wall.model.cavity.xyz.copy()

    displaced = mol.atom_coords()
    displaced[0, 2] += 0.05
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)
    wall.model.update(Structure(numbers, displaced))

    # Vacuous unless the surface moved with its point count kept
    # - the case a shape check cannot distinguish
    assert wall.model.ngrid == wall.ngrid
    assert not np.allclose(wall.model.cavity.xyz, original_xyz)

    with pytest.raises(RuntimeError, match="not supplied by the host"):
        wall.model.get_energy()

    # Wrapper re-supplies, energy tracks the geometry
    _host2, moved = make_wall(mol, displaced, dm=dm)
    assert moved.energy != pytest.approx(before, rel=1e-9)


@pytest.mark.gostshyp_params
def test_params_second_moment_matches_a_second_difference():
    """``Mt_ab`` against a second difference of ``gtilde`` in the grid point center

    With ``G = exp(-omega rho^2)``,
    ``d^2 gtilde/dC_a dC_b = 4 omega^2 Mt_ab - 2 omega delta_ab gtilde``
    - isolates the d-fakemol constant and the ``_D_CART_ORDER`` mapping
    - ``test_params_derivatives_match_fd`` sees them only through the contracted
      ``n . Mt``
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    dm_cart = wall._density_matrix_cart(dm)

    _gt, _pt, moment, _rt = wall._surface_moments(dm_cart)
    gtilde = wall.traces(dm)[0]
    omega = wall.omega
    step = 2.0e-3

    def gt_at(shift):
        return wall.traces(dm, centers=wall.centers + shift)[0]

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
    """One PySCF coupling completes GOSTSHYP alone or alongside CPCM"""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)

    def evaluate(components):
        host = PySCFHost(mol)
        cavity = CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL)
        model = SolvationModel(cavity=cavity, components=components)
        return model.evaluate(coupling=host.coupling(dm))

    gostshyp = evaluate([ModelComponentGOSTSHYP(PRESSURE)])
    cpcm = evaluate([ModelComponentCPCM(EPSILON)])
    combined = evaluate([ModelComponentCPCM(EPSILON), ModelComponentGOSTSHYP(PRESSURE)])

    assert combined.energy == pytest.approx(cpcm.energy + gostshyp.energy, rel=1e-12)
    np.testing.assert_allclose(combined.fock, cpcm.fock + gostshyp.fock, rtol=1e-11)
    np.testing.assert_allclose(
        combined.gradient,
        cpcm.gradient + gostshyp.gradient,
        rtol=1e-10,
        atol=1e-11,
    )


@pytest.mark.gostshyp_fock
def test_evaluation_owns_complete_gostshyp_results():
    """Generic result includes the wall's reverse-mode host response"""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = PySCFHost(mol)
    with pytest.deprecated_call(match="Gostshyp"):
        wall = GostshypModel(host, PRESSURE, nleb=NLEB, tolerance=PROJ_TOL)

    result = wall.evaluate(dm)
    expected_gradient = wall.nuclear_gradient(dm)
    host.dm = np.zeros_like(dm)

    assert result.energy == wall.energy
    np.testing.assert_allclose(result.fock, wall.fock(dm))
    np.testing.assert_allclose(result.gradient, expected_gradient)
    assert not result.fock.flags.writeable
    assert not result.gradient.flags.writeable


@pytest.mark.gostshyp_conventions
def test_gostshyp_wall_name_remains_compatible():
    mol = molecule(*PRIMARY_CASE)
    host = PySCFHost(mol)

    with pytest.deprecated_call(match="GostshypModel"):
        wall = GostshypWall(host, PRESSURE, nleb=NLEB, tolerance=PROJ_TOL)

    assert isinstance(wall, GostshypModel)


@pytest.mark.gostshyp_fock
def test_fock_energy_is_linear_in_pressure_and_sets_the_effective_volume():
    """``E = p V_eff`` with ``V_eff`` independent of ``p``, and ``E(0) = 0``

    Pure identities on one cavity
    - surface does not depend on the pressure; amplitudes are proportional to it
    - a mis-masked ``active`` breaks the volume identity without breaking linearity
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
    # Not the cavity volume: wall is weighted by the local density overlap
    assert wall.effective_volume() != pytest.approx(wall.model.cavity.volume, rel=1e-3)

    # Three walls share a density, hence a surface, hence the volume
    # - at p = 0 the energy vanishes but the volume does not: not recoverable as E/p
    assert vacuum.effective_volume() == pytest.approx(wall.effective_volume(), rel=1e-12)
    assert doubled.effective_volume() == pytest.approx(wall.effective_volume(), rel=1e-12)


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_fock_frozen_surface_matches_fd(system, basis):
    """Eq-16 Fock is ``dE/dP`` with the grid points held fixed

    - isolates the explicit density dependence from the cavity response, so a
      headline-test failure is attributable to one or the other
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
        frozen = wall.fock(dm, include_cavity_response=False)
        analytic = float(np.einsum("uv,uv->", frozen, direction))
        assert deviation(
            analytic, numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR
        ) <= 1.0


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_fock_matches_fd(system, basis):
    """dE/dP with the surface following the density, the headline Fock test

    Closes the whole chain
    - surface weights
    - moist's contraction of them into level-set adjoints
    - host's contraction of those with dS/dP
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
    """AO centers move, grid points frozen: the ``int3c1e_ip1`` term alone

    - pins the sign of the integral derivative, the both-legs factor of two, and
      the mapping of cartesian AO rows onto atoms
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
def test_gradient_surface_channel_matches_fd(system, basis):
    """Cavity's own response at a frozen level-set field

    - displacing only the structure handed to the cavity, never the molecule the
      level set is built from, drags the atom-anchored reference grid while the
      density stays put
    - that is the channel moist contracts in reverse mode: closes the component's
      ``w_a``/``w_xyz``/``w_n`` against an FD of the energy through moist's own
      cavity derivatives
    - area route, including its switching factor
    - position route, including the normal's point-motion fold
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host, wall = make_wall(mol, dm=dm)
    analytic = wall._surface_nuclear_gradient().flatten(order="F")
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)

    def anchored_energy(displaced):
        # Cavity moves; host.mol, hence the level set, does not
        cavity = CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL)
        cavity.update(Structure(numbers, displaced))
        result = cavity.cavity
        centers = np.ascontiguousarray(result.xyz.T)
        areas = np.ascontiguousarray(result.a)
        normals = np.array(result.normal0.T, order="C", copy=True)
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
    """Total dE/dR at fixed density, the headline nuclear-gradient test"""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)
    analytic = wall.nuclear_gradient(dm).flatten(order="F")

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, positions, dm, index)
        assert deviation(analytic[index], numerical) <= 1.0


@pytest.mark.gostshyp_gradient
def test_gradient_is_translationally_invariant():
    """Net force vanishes

    - a rigid translation moves AOs, density and cavity together: energy cannot change
    - global identity blind to the channel split, so it catches a sign error or
      mis-binned AO-to-atom map that each per-channel test tolerates
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    gradient = wall.nuclear_gradient(dm)

    assert np.abs(gradient).max() > MIN_SIGNAL
    np.testing.assert_allclose(gradient.sum(axis=1), 0.0, atol=1e-10)


@pytest.mark.gostshyp_gradient
def test_gradient_is_exactly_the_sum_of_its_channels():
    """No channel silently dropped inside :meth:`nuclear_gradient`"""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    channels = (
        wall._integral_nuclear_gradient(dm)
        + wall._field_nuclear_gradient(dm)
        + wall._surface_nuclear_gradient()
    )
    np.testing.assert_allclose(wall.nuclear_gradient(dm), channels, rtol=0.0, atol=1e-14)


# ----------------------------------------------------------------------
# conventions -- channel semantics and the negative controls
# ----------------------------------------------------------------------


# Surface weights are not visible from Python: built and consumed inside the moist
# ``gostshyp`` component
# Moved to test/unit/test_model/component/gostshyp.f90
# - ``w_f == 0`` exactly
# - per-channel starve control for ``w_n``/``w_a``
# - ``check_surface_weights`` finite-differences each channel separately, not
#   just detecting that some channel is missing


@pytest.mark.gostshyp_conventions
def test_conventions_fock_requires_cavity_response():
    """Eq-16 Fock alone is not ``dE/dP`` for a surface following the density

    - no electrostatic component to dilute the error
    - cavity response is the majority of the derivative, not a correction
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
def test_conventions_amplitudes_are_both_needed():
    """Both host amplitudes carry a non-negligible share of the frozen Fock

    - ``w_overlap`` and ``w_normal_deriv`` arrive as a pair with signs already folded in
    - natural failures: using one and dropping the other, or flipping the fold
    - either leaves a plausible-looking Fock
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    response = wall.response

    assert np.abs(response.gostshyp.w_overlap).max() > MIN_SIGNAL
    assert np.abs(response.gostshyp.w_normal_deriv).max() > MIN_SIGNAL
    # Fold is a sign, not an absolute value: the two channels oppose
    assert response.gostshyp.w_overlap.max() > 0.0
    assert response.gostshyp.w_normal_deriv.min() < 0.0

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd4(
        [
            frozen_surface_energy(
                mol,
                dm + offset * STEP_DM * direction,
                wall.centers,
                wall.areas,
                wall.omega,
                wall.normals,
            )
            for offset in FD4_OFFSETS
        ],
        STEP_DM,
    )

    g_only = np.einsum(
        "j,uvj->uv", response.gostshyp.w_overlap, wall._G, optimize=True
    )
    starved = float(np.einsum("uv,uv->", 0.5 * (g_only + g_only.T), direction))
    assert deviation(starved, numerical, thr_rel=FROZEN_REL_THR) > VACUITY_FACTOR
    # Pair together is the quantity that closes
    full = float(np.einsum("uv,uv->", wall.fock(dm, include_cavity_response=False), direction))
    assert deviation(full, numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR) <= 1.0


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize("dropped", ["integral", "field", "surface"])
def test_conventions_gradient_requires_every_channel(dropped):
    """Each of the three nuclear routes carries a non-negligible share

    - two are the host's, one is moist's: pins the split itself
    - dropping the moist term must break the gradient, else the component does not
      contribute what its surface weights claim
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    _host, wall = make_wall(mol, dm=dm)

    channels = {
        "integral": wall._integral_nuclear_gradient(dm),
        "field": wall._field_nuclear_gradient(dm),
        "surface": wall._surface_nuclear_gradient(),
    }
    assert np.abs(channels[dropped]).max() > MIN_SIGNAL
    starved = sum(value for name, value in channels.items() if name != dropped)

    index = int(np.argmax(np.abs(channels[dropped].flatten(order="F"))))
    numerical = fd_position(mol, positions, dm, index)
    assert deviation(starved.flatten(order="F")[index], numerical) > VACUITY_FACTOR


# Two controls retired with the port: the starved quantity is no longer reachable
# from Python
# - area route must use the true ``da_i/dR_A``, not the Gaussian-width proxy
#   - now moist's internal choice, covered end to end by
#     ``test_gradient_surface_channel_matches_fd``
# - activity floor must gate the amplitudes, not only ``dE/da``
#   - now structural: the component computes the mask once, beside the amplitudes
#   - ``test_gostshyp_energy_matches_amplitudes`` in
#     test/unit/test_model/component/gostshyp.f90 checks energy and host
#     amplitudes agree, the property a split mask breaks


@pytest.mark.gostshyp_conventions
def test_conventions_anchor_area_derivatives_sum_to_the_total():
    """``sum_i a_i1_rA == A_tot1_rA``, pinning the anchor buffer layout

    - binding hands moist a rank-4 Fortran buffer and a capacity pair ordered
      opposite to the array indexing
    - independent of GOSTSHYP physics: bisects a layout error from a weight error
    - driven on a standalone cavity, as the wall's cavity now belongs to the model;
      checks the binding, not the wall
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = PySCFHost(mol)
    host.dm = dm
    cavity = CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL)
    cavity.update(host.structure())

    cavity.compute_anchor_gradient()
    anchor = cavity.get_anchor_gradient()
    np.testing.assert_allclose(
        anchor.a_i1_rA.sum(axis=2), anchor.A_tot1_rA, rtol=0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        anchor.v_i1_rA.sum(axis=2), anchor.V_tot1_rA, rtol=0.0, atol=1e-12
    )


@pytest.mark.gostshyp_conventions
def test_conventions_failed_update_leaves_no_stale_results(monkeypatch):
    """A rebuild that raises must not leave the previous geometry readable

    - :meth:`update` moves the host density and the cavity before it can fail
      (density with no isodensity surface, projection that will not converge)
    - results cached from the previous call then describe a surface that no longer
      exists
    - failure is forced, not hunted: the invalidation is pinned, not any particular
      way of tripping it
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # Vacuous unless a real result exists to go stale
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
def test_conventions_failed_update_publishes_nothing_partial(monkeypatch):
    """A failure in the *last* step of :meth:`update` publishes nothing either

    - cavity and moments are built by then; the traces and energy waiting behind
      that call are good numbers, the tempting ones to keep
    - :meth:`update` raised, so no caller was handed the state they belong to
    - a half-published result is one no accessor can tell from a whole one
    - the sibling test injects its failure in ``_build_integrals``, never reaching
      this path
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    def boom(*_args, **_kwargs):
        raise RuntimeError("response assembly failed")

    monkeypatch.setattr(wall.model, "response", boom)
    with pytest.raises(RuntimeError, match="response assembly failed"):
        wall.update(dm)

    assert wall.amplitudes is None
    assert wall.energy == 0.0
    for read in (
        wall.effective_volume,
        lambda: wall.live_traces,
        lambda: wall.inactive_count,
        lambda: wall.fock(dm),
    ):
        with pytest.raises(RuntimeError, match="update"):
            read()


@pytest.mark.gostshyp_conventions
def test_conventions_moments_are_built_once_per_update(monkeypatch):
    """One moment build per :meth:`update`, none per result read

    - :meth:`_surface_moments` is four dense three-center integral blocks (d and f
      shells among them), the most expensive thing in this module
    - needed once, to hand the component its moments; later reads contract against
      cached amplitudes
    - returned numbers would not reveal a rebuild, so the count is pinned
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    original = wall._surface_moments
    calls = []

    def counted(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(wall, "_surface_moments", counted)

    wall.update(dm)
    assert len(calls) == 1

    fock = wall.fock(dm)
    gradient = wall.nuclear_gradient(dm)
    volume = wall.effective_volume()
    assert len(calls) == 1, "a result read rebuilt the moments"
    assert np.abs(fock).max() > MIN_SIGNAL
    assert np.abs(gradient).max() > MIN_SIGNAL
    assert volume > 0.0

    # Cached traces equal a fresh build's
    gt, ftilde = wall.live_traces
    rebuilt_gt, rebuilt_ftilde = wall.traces(dm)
    np.testing.assert_allclose(gt, rebuilt_gt, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(ftilde, rebuilt_ftilde, rtol=0.0, atol=0.0)


def _golden_summary(name, array):
    """Reduce an array the way ``scratchpad/capture_golden.py`` did"""
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
    """Every reported quantity still equals its pre-port value

    FDs elsewhere are self-consistent (model derivative vs the model's own energy)
    - a convention error applied to both (wrong angular constant, mis-ordered
      cartesian component) leaves all of them passing
    - this is the only test that notices, hence it covers energy, amplitudes,
      level-set adjoints and gradient

    Parametrised over ``GOLDEN``, *not* ``CASES``
    - values come from the pure-Python implementation that preceded the Fortran
      component, so they exist only for the cases enabled then
    - enabling a new row in ``CASES`` must not fail here

    Do not simply add a new row to ``GOLDEN``
    - the pre-port implementation is gone: new values would come from the code
      under test
    - that is a regression pin against drift, useful, but not the independent
      cross-check the water rows carry
    - if one is added, say in a comment which of the two it is
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    weights = wall.response
    actual = {"energy": float(wall.energy), "ngrid": int(wall.ngrid)}
    actual.update(_golden_summary("alpha", wall.amplitudes))
    actual.update(
        _golden_summary("fock_frozen", wall.fock(dm, include_cavity_response=False))
    )
    actual.update(_golden_summary("w_lsf0", weights.lsf.w_value))
    actual.update(_golden_summary("w_lsf1", weights.lsf.w_gradient))
    actual.update(_golden_summary("w_lsf2", weights.lsf.w_hessian))
    actual.update(_golden_summary("gradient", wall.nuclear_gradient(dm)))

    expected = GOLDEN[(system, basis)]
    assert actual.keys() == expected.keys()
    for key, reference in expected.items():
        if key == "ngrid":
            assert actual[key] == reference, key
        elif abs(reference) < MIN_SIGNAL:
            # Gradient translational invariance; sign-definite channels with extremum
            # exactly zero
            assert abs(actual[key]) < MIN_SIGNAL, key
        else:
            assert actual[key] == pytest.approx(reference, rel=GOLDEN_RTOL), key
