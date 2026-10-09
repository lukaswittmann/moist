"""End-to-end tests of the GOSTSHYP pressure model on a moist isodensity cavity.

Every analytic quantity is checked against a finite difference of the energy
the model reports, at 50 GPa.

This file tests the *host* half and the assembled result. The component's
own surface weights are finite-differenced channel by channel in
``test/unit/test_model_component_gostshyp.f90``, against a model density
whose Gaussian moments are available in closed form.

The suite is organised into layers:

``integrals``
    The libcint fakemol constants, the cartesian d/f component orders, and
    the identity ``f == n . grad_r g``. moist is used only to place the grid
    points, so a failure here is integral bookkeeping, not the cavity. These
    constants are the one convention that stayed host-side.
``params``
    The Gaussian moments this module hands the component, and the hand-off
    itself: a transposed moment is the failure mode the language boundary
    introduces.
``fock``
    Energy and ``dE/dP`` with the surface following the density -- the
    level-set chain rule.
``gradient``
    ``dE/dR`` at fixed density, split into its integral, field and surface
    channels and then assembled. The first two are the host's; the third is
    moist contracting the component's weights in reverse mode.
``conventions``
    Channel semantics, the negative controls, and the frozen reference values.

Each layer carries a negative control. GOSTSHYP has no electrostatic
component: the cavity response is ~84% of ``dE/dP`` rather than a
correction, and all three gradient channels are comparable in size (1.3e-2,
1.4e-2, 1.5e-2 for water/STO-3G at 50 GPa).

The finite differences here are self-consistent: they compare the model's
derivative against the model's own energy, so a convention error applied
consistently to both would leave every one of them passing. ``GOLDEN`` and
``test_conventions_matches_the_frozen_reference`` check the energy,
amplitudes, adjoints and gradient against frozen reference values instead.
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

from moist.interface import (
    GaussianAmplitudeResponse,
    GaussianMomentRequest,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    SolvationModel,
    Structure,
)
from moist.parameters import DROPParameters, IsodensityParameters, GOSTSHYPParameters
from moist import Context
from moist.pyscf import (
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

#: Run context shared by every cavity and model in this module
CONTEXT = Context()

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
#: Grid-center displacement for the eighth-order second difference (Bohr).
STEP_SECOND_MOMENT = 2e-2

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
    # One atom: no seams, so no narrow points
    "neon": System("Ne 0.0 0.0 0.0"),
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

#: Reference values
GOLDEN = {
    ("water", "sto-3g"): {"energy": 0.1340185971628404, "ngrid": 71},
    ("water", "def2-svp"): {"energy": 0.15034462492176853, "ngrid": 71},
}

#: The golden values are a different platform's arithmetic away from exact,
#: but the cavity projection is the only real variability, and it converges
#: to PROJ_TOL.
GOLDEN_RTOL = 1e-9


CASE_PARAMS = [
    pytest.param(system, basis, id=f"{system}-{basis}") for system, basis in CASES
]


def deviation(actual, reference, *, thr_abs=None, thr_rel=None) -> float:
    """Deviation measured in units of the tolerance: ``<= 1`` passes.

    Combines the two thresholds the way test-drive's
    ``check(..., thr_abs=, thr_rel=)`` does: ``max(thr_abs, thr_rel*|ref|)``.
    The absolute floor covers references near zero; the relative bound
    scales with the magnitude. The return value is the number of tolerances
    ``actual`` missed ``reference`` by.
    """
    thr_abs = ABS_THR if thr_abs is None else thr_abs
    thr_rel = REL_THR if thr_rel is None else thr_rel
    return abs(actual - reference) / max(thr_abs, thr_rel * abs(reference))


def array_deviation(actual, reference, *, thr_rel, thr_abs=None) -> float:
    """Worst elementwise deviation, scaled by the *array's* magnitude.

    Per-element relative errors are meaningless for arrays whose entries span
    orders of magnitude -- the small entries are differences of large ones.
    The tolerance is set from the largest entry of the reference.
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
    lsf = Isodensity() if rho_iso is None else Isodensity(
        parameters=IsodensityParameters(rho_iso=rho_iso))
    return DROP(lsf=lsf, parameters=DROPParameters(nleb=NLEB, tolerance=PROJ_TOL))


def make_wall(mol, positions=None, *, dm, pressure=PRESSURE, rho_iso=None, parameters=None):
    """Host plus an evaluated GOSTSHYP driver ("wall") at ``dm``, optionally displaced."""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    wall = PySCFSolvation(mol, cavity_config(rho_iso), [ModelComponentGOSTSHYP(pressure, parameters=parameters)],
                          context=CONTEXT)
    wall.evaluate(dm)
    return wall.host, wall


def item_of(response, kind):
    """The item of class ``kind`` the walk over ``response`` meets, or None."""
    found = [item for item in response if isinstance(item, kind)]
    assert len(found) <= 1, found
    return found[0] if found else None


def fd_density(mol, dm, direction, *, pressure=PRESSURE, parameters=None):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        _host, wall = make_wall(mol, dm=dm + offset * STEP_DM * direction, pressure=pressure, parameters=parameters)
        samples.append(wall.energy)
        grids.append(wall.model.cavity.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, pressure=PRESSURE, parameters=None):
    """dE/dR along one cartesian coordinate at fixed density matrix."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        displaced = positions.copy()
        displaced.flat[index] += offset * STEP_R
        _host, wall = make_wall(mol, displaced, dm=dm, pressure=pressure, parameters=parameters)
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


def component_energy(model, mol, dm):
    """The model's energy on its current surface, answered by a fresh host for ``mol`` and ``dm``.

    The component applies its own widths, switches and trace window, so none
    of its numerics is restated here. Moving ``mol`` moves the basis centres
    while the grid stays where it is.
    """
    moments = GaussianMoments(PySCFHost(mol))
    moments.bind(model.cavity, dm)
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    for request in coupling:
        assert isinstance(request, GaussianMomentRequest), request.name
        moments.answer(coupling, request)
    energy = np.array(0.0)
    model.get_energy(coupling, energy)
    return float(energy)


# ----------------------------------------------------------------------
# integrals -- libcint conventions, moist only places the grid points
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_integrals
def test_integrals_fakemol_normalization_constants():
    """The four angular constants and the cartesian d/f orders, vs quadrature.

    libcint attaches a different normalization to a coefficient-1 shell of
    each angular momentum. Restoring the ratios ``N_s/N_l`` puts every
    moment in the same units as ``g`` and makes ``f = n . grad g`` exact
    rather than exact-up-to-a-factor. An independent Becke-grid quadrature
    of the same monomial moments checks the constants and the component
    orders together, catching a transposed ``_D_CART_ORDER`` entry or a
    mis-assigned ``_F_RHO2_FIRST_MOMENT`` triple.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # A handful of well-conditioned grid points is enough.
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
        normalized_s = _int3c1e(mol, centers, omega, 0, normalized=True)[:, :, slot]
        np.testing.assert_allclose(normalized_s, (omega[slot]/np.pi)**1.5 * quad(np.ones_like(rho2)),
                                   rtol=QUAD_RTOL, atol=QUAD_ATOL)
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
    center, the sign convention used in ``_build_integrals``. Get it wrong
    and the wall pushes outward instead of inward.
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

    ``gvfield`` is the unprojected f-vector and doubles as ``dftilde/dn``.
    The two must agree to round-off, or the normal weight is built from a
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

    ``Pt`` is the only moment the energy depends on, through
    ``ftilde = -2 omega (n . Pt)``. The Becke quadrature in the integrals
    layer checks the fakemol block; this checks the contraction that turns
    it into the moment the component is handed.
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
def test_params_moments_reach_the_component_intact():
    """The energy moist reports is the one the host's own moments imply.

    The moments cross a language boundary and change layout on the way --
    numpy ``(ngrid, 3)`` becomes Fortran ``(3, ngrid)``. A transposed hand-off
    is a failure mode nothing else in this layer would see. On one atom no
    point is narrow and every trace is large, so the documented formula
    ``E = sum_j p a_j gtilde_j / ftilde_j`` holds with nothing switched.
    """
    mol, dm = molecule("neon", "def2-svp"), reference_density("neon", "def2-svp")
    _host, wall = make_wall(mol, dm=dm)

    gt, ftilde = wall.moments.traces(dm)
    # Far from the narrow-point switch and from the trace window
    assert wall.moments.omega.max() < 10.0
    assert np.abs(ftilde).min() > 1e3 * GOSTSHYPParameters().regularization_end
    expected = PRESSURE * float(np.sum(wall.moments.areas * gt / ftilde))
    assert wall.energy == pytest.approx(expected, rel=1e-12)


@pytest.mark.gostshyp_params
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_params_window_only_moves_points_inside_it(system, basis):
    """Windows below every trace change nothing; a start above every trace turns all off.

    In between each point keeps a share ``S(|f|)`` in ``[0, 1]`` of a positive
    contribution, so the energy lies strictly between zero and the plain
    energy. The switch's exact values are pinned by
    ``gostshyp_regularization_branches`` in the Fortran suite.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    _gt, ftilde = wall.moments.live_traces
    live = wall.moments.omega > 0.0
    magnitude = np.abs(ftilde[live])
    # The ordering below needs positive contributions
    assert np.all(ftilde[live] > 0.0)
    assert magnitude.min() > 1e3 * GOSTSHYPParameters().regularization_end

    def energy(start, end):
        parameters = GOSTSHYPParameters(regularization_start=start, regularization_end=end)
        return make_wall(mol, dm=dm, parameters=parameters)[1].energy

    for start, end in ((1e-14, 1e-12), (0.25 * magnitude.min(), 0.5 * magnitude.min())):
        assert energy(start, end) == pytest.approx(wall.energy, rel=1e-12), (start, end)
    assert energy(magnitude.max(), 2.0 * magnitude.max()) == 0.0

    # Windows centred on the tenth and the half quantile of the points still on
    for quantile in np.quantile(magnitude, (0.1, 0.5)):
        banded = energy(0.5 * quantile, 2.0 * quantile)
        assert 0.0 < banded < wall.energy, quantile
        assert wall.energy - banded > VACUITY_FACTOR * 1e-12 * wall.energy, quantile


@pytest.mark.gostshyp_params
def test_params_update_drops_the_previous_surfaces_moments():
    """A cavity update invalidates the moments, and reading without new ones fails.

    moist cannot rebuild the moments itself: they are AO-basis integrals
    only the host can form, and their shapes cannot distinguish a stale set
    from a current one. The update drops them; reading afterward raises
    rather than returning an energy for a surface that has moved.

    :meth:`PySCFSolvation.evaluate` always re-supplies immediately; only a
    host driving the model by hand hits this.
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
    step = STEP_SECOND_MOMENT

    def gt_at(shift):
        return wall.moments.traces(dm, centers=wall.moments.centers + shift)[0]

    def directional_second(direction):
        # Eighth-order second difference of values, at +/-1..4 steps.
        # Subtract the center before summing to limit cancellation.
        center = gtilde
        second = np.zeros_like(center)
        for offset, weight in enumerate((8 / 5, -1 / 5, 8 / 315, -1 / 560), 1):
            plus = gt_at(offset * step * direction)
            minus = gt_at(-offset * step * direction)
            second += weight * ((plus - center) + (minus - center))
        return second / step**2

    # h/2 and 2h also pass; diagonal cancellation near zero is covered by
    # the absolute floor without discarding any component or grid point.
    references = {}
    for a in range(3):
        for b in range(3):
            key = tuple(sorted((a, b)))
            if key not in references:
                ea, eb = np.eye(3)[a], np.eye(3)[b]
                if a == b:
                    references[key] = directional_second(ea)
                else:
                    references[key] = (
                        directional_second(ea + eb) - directional_second(ea - eb)
                    ) / 4.0
            second = references[key]
            # Reconstruct the moment from the finite-difference Hessian. The
            # reverse subtraction loses the Hessian for omega > 1e12 even if
            # the supplied moment is correctly rounded.
            on = omega > 0.0
            reference = (second[on] + 2.0 * omega[on] * gtilde[on] * (a == b)) / (4.0 * omega[on]**2)
            threshold = np.maximum(PARAM_ABS_THR / (4.0 * omega[on]**2),
                                   PARAM_REL_THR * np.abs(reference))
            assert np.all(np.abs(moment[on, a, b] - reference) <= threshold), (a, b)
            # Switched-off points are requested at zero width and carry nothing
            assert np.all(moment[~on, a, b] == 0.0), (a, b)


# ----------------------------------------------------------------------
# fock -- the level-set chain rule
# ----------------------------------------------------------------------


@pytest.mark.gostshyp_fock
def test_gostshyp_is_a_composable_model_component():
    """One PySCF coupling completes GOSTSHYP alone or alongside CPCM."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)

    def evaluate(components):
        solvation = PySCFSolvation(mol, cavity_config(), components, context=CONTEXT)
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
def test_fock_energy_is_linear_in_pressure():
    """``E`` is odd and linear in ``p``, with and without points in the band.

    The band acts on ``ftilde``, a property of the density alone, so the same
    points sit inside it at every pressure. The isodensity surface does not
    depend on the pressure either.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, reference = make_wall(mol, dm=dm)
    _gt, ftilde = reference.moments.live_traces
    live = ftilde != 0.0
    smallest = float(np.abs(ftilde[live]).min())
    banded = GOSTSHYPParameters(regularization_start=0.5 * smallest, regularization_end=2.0 * smallest)
    assert np.any(live & (np.abs(ftilde) < banded.regularization_end))

    for parameters in (None, banded):
        _host, wall = make_wall(mol, dm=dm, parameters=parameters)
        _host, doubled = make_wall(mol, dm=dm, pressure=2.0 * PRESSURE, parameters=parameters)
        _host, inverted = make_wall(mol, dm=dm, pressure=-PRESSURE, parameters=parameters)
        _host, vacuum = make_wall(mol, dm=dm, pressure=0.0, parameters=parameters)
        assert vacuum.model.cavity.ngrid == wall.model.cavity.ngrid
        assert wall.energy > 0.0
        assert doubled.energy == pytest.approx(2.0 * wall.energy, rel=1e-12)
        assert inverted.energy == pytest.approx(-wall.energy, rel=1e-12)
        assert vacuum.energy == 0.0
        # Not the cavity volume: the wall is weighted by the local density overlap
        assert wall.energy / PRESSURE != pytest.approx(wall.model.cavity.volume, rel=1e-3)


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_fock_frozen_surface_matches_fd(system, basis):
    """The eq-16 Fock is ``dE/dP`` with the grid points held fixed.

    Isolates the explicit density dependence from the cavity response.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    for direction in symmetric_directions(dm.shape[0], 2):
        samples = [
            component_energy(wall.model, mol, dm + offset * STEP_DM * direction)
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
    """AO centers move, grid points frozen: the ``int3c1e_ip1`` term alone.

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
            samples.append(component_energy(wall.model, moved, dm))
        numerical = fd4(samples, STEP_R)
        assert deviation(
            analytic[index], numerical, thr_abs=FROZEN_ABS_THR, thr_rel=FROZEN_REL_THR
        ) <= 1.0


@pytest.mark.gostshyp_gradient
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_gradient_surface_channel_matches_fd(system, basis):
    """The cavity's own response at a frozen level-set field.

    Displacing only the structure handed to the cavity -- never the molecule
    the level set is built from -- drags the atom-anchored reference grid
    while the density stays put: the channel moist contracts in reverse
    mode. This checks the component's ``w_a``/``w_xyz``/``w_normal`` against
    a finite difference of the energy through moist's own cavity
    derivatives: the area route including its switching factor, and the
    position route including the normal's point-motion fold.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host, wall = make_wall(mol, dm=dm)
    analytic = wall.gradient_channels(dm)["model"].flatten(order="C")
    numbers = np.asarray(mol.atom_charges(), dtype=np.int64)

    def anchored_energy(displaced):
        # The cavity moves; host.mol -- and hence the level set -- does not.
        model = SolvationModel(CONTEXT, cavity_config().build(source=host), [ModelComponentGOSTSHYP(PRESSURE)])
        model.update(Structure(numbers, displaced))
        return component_energy(model, mol, dm), model.cavity.ngrid

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

    A rigid translation moves the AOs, the density and the cavity together,
    so the energy cannot change. This is a global identity that knows
    nothing about the channel split, so it catches a sign error or a
    mis-binned AO-to-atom map that each per-channel test tolerates.
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


# The surface weights are built and consumed inside the moist ``gostshyp``
# component and are not visible from Python. ``w_f == 0`` exactly, and a
# per-channel starve control for ``w_normal``/``w_a``, are checked in
# test/unit/test_model_component_gostshyp.f90's ``check_surface_weights``,
# which finite-differences each channel separately.


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

    ``w_overlap`` and ``w_normal_deriv`` arrive as a pair with their signs
    already folded in. A natural failure is using one and dropping the
    other, or flipping the fold; either would leave a Fock that still looks
    plausible.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)
    amplitude = item_of(wall.response, GaussianAmplitudeResponse)

    assert np.abs(amplitude.w_overlap).max() > MIN_SIGNAL
    assert np.abs(amplitude.w_normal_deriv).max() > MIN_SIGNAL
    # The fold is a sign, not an absolute value: the two channels oppose.
    assert amplitude.w_overlap.max() > 0.0
    assert amplitude.w_normal_deriv.min() < 0.0

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd4(
        [
            component_energy(wall.model, mol, dm + offset * STEP_DM * direction)
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

    Two are the host's and one is moist's; this also pins the split itself.
    Dropping the moist term must break the gradient, or the component is not
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


# Two constraints on the amplitudes, checked elsewhere:
#
#   * the area route must use the true ``da_i/dR_A`` rather than a
#     Gaussian-width proxy, covered end to end by
#     ``test_gradient_surface_channel_matches_fd``;
#   * the regularization must act on the amplitudes, not only ``dE/da``. The
#     component evaluates the band once, beside the amplitudes, and
#     ``gostshyp_energy_matches_amplitudes`` in
#     test/unit/test_model_component_gostshyp.f90 checks that the energy and
#     the host amplitudes agree.


@pytest.mark.gostshyp_conventions
def test_conventions_anchor_area_derivatives_sum_to_the_total():
    """``sum_i a_i1_rA == A_tot1_rA``, pinning the anchor buffer layout.

    The binding exposes C arrays with the native axes reversed. This
    identity is independent of any GOSTSHYP physics, isolating a layout
    error from a weight error.

    Driven on a standalone cavity: the wall's cavity belongs to the model,
    so this checks the binding rather than the wall.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = PySCFHost(mol)
    host.dm = dm
    cavity = cavity_config().build(source=host, context=CONTEXT)
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

    :meth:`evaluate` moves the host density and the cavity before it can
    fail -- a density with no isodensity surface, a projection that will
    not converge. Results cached from the previous call would describe a
    surface that no longer exists. The failure is forced here rather than
    hunted for: this pins the invalidation, not any particular way of
    tripping it.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm)

    # Vacuous unless there is a real result to go stale.
    assert wall.energy > 0.0
    assert item_of(wall.response, GaussianAmplitudeResponse) is not None

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

    The cavity and the moments are built by then, and the energy waiting
    behind that call is a perfectly good number. :meth:`evaluate` raised, so
    no caller was ever handed the state it belongs to, and a half-published
    result is one no accessor can tell from a whole one.
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
    gt, ftilde = wall.moments.live_traces
    assert len(calls) == 2, "a result read rebuilt the moments"
    assert np.abs(fock).max() > MIN_SIGNAL
    assert np.abs(gradient).max() > MIN_SIGNAL

    # The cached traces are the ones a fresh build would produce.
    rebuilt_gt, rebuilt_ftilde = wall.moments.traces(dm)
    np.testing.assert_allclose(gt, rebuilt_gt, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(ftilde, rebuilt_ftilde, rtol=0.0, atol=0.0)


@pytest.mark.gostshyp_conventions
@pytest.mark.parametrize(
    "system,basis",
    [pytest.param(*case, id=f"{case[0]}-{case[1]}") for case in GOLDEN],
)
def test_conventions_matches_the_frozen_reference(system, basis):
    """The energy over the reference's point set still equals its pre-port value.

    The finite differences elsewhere in this file are self-consistent: they
    check the model's derivative against the model's own energy, so a
    convention error applied consistently to both -- a wrong angular
    constant, a mis-ordered cartesian component -- leaves all of them
    passing. The pre-port energy does not.

    The reference used raw shells on every point and cut points at a relative
    floor on the raw trace. Its point set is rebuilt from traces at the full
    width ``pi ln2/a``; the normalization cancels in ``a g/f``. Its amplitudes,
    adjoints and gradient include points the narrow-point switch now turns
    off, so only the energy and the grid remain comparable.

    Parametrised over ``GOLDEN``, *not* ``CASES``: these values exist only
    for the cases enabled when they were captured. A value captured for a new
    case would come from the code under test, not the pre-port
    implementation -- say in a comment if a row is such a regression pin.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    _host, wall = make_wall(mol, dm=dm, rho_iso=4e-4)
    assert int(wall.model.cavity.ngrid) == GOLDEN[(system, basis)]["ngrid"]

    areas = wall.moments.areas
    assert np.all(areas > 0.0)
    omega = np.pi * math.log(2.0) / areas
    gt, ftilde = wall.moments.traces(dm, omega=omega)
    raw_f = ftilde / ((omega / np.pi)**1.5 / _S_NORM)
    kept = np.abs(raw_f) > 1e-9 * np.abs(raw_f).max()
    energy = PRESSURE * float(np.sum(areas[kept] * gt[kept] / ftilde[kept]))
    assert energy == pytest.approx(GOLDEN[(system, basis)]["energy"], rel=GOLDEN_RTOL)


@pytest.mark.gostshyp_fock
@pytest.mark.parametrize("suppress", [False, True])
def test_regularized_fock_matches_energy_difference(suppress):
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    parameters = GOSTSHYPParameters(regularization_start=5e-3, regularization_end=2e-2,
                                    suppress_negative_amplitudes=suppress)
    _, wall = make_wall(mol, dm=dm, parameters=parameters)
    magnitude = np.abs(wall.moments.live_traces[1])
    assert np.any((magnitude > parameters.regularization_start) & (magnitude < parameters.regularization_end))
    for direction in symmetric_directions(mol.nao_nr(), 2):
        analytic = np.einsum("uv,uv", wall.fock, direction)
        numerical = fd_density(mol, dm, direction, parameters=parameters)
        assert deviation(analytic, numerical) <= 1.0


@pytest.mark.gostshyp_gradient
def test_regularized_nuclear_gradient_matches_energy_difference():
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    parameters = GOSTSHYPParameters(regularization_start=5e-3, regularization_end=2e-2)
    _, wall = make_wall(mol, dm=dm, parameters=parameters)
    analytic = wall.gradient(dm)
    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, mol.atom_coords(), dm, index, parameters=parameters)
        assert deviation(analytic.flat[index], numerical) <= 1.0, index
