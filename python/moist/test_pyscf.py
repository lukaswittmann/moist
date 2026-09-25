"""End-to-end tests of the moist/host chain rule against PySCF

moist hands out adjoint weights and expects the host to finish the chain rule
with its own density derivatives, so the correctness of an isodensity cavity is
split across a language boundary.  These tests close that loop: every analytic
quantity is checked against a finite difference of the energy that moist itself
returns.

The suite is layered so a failure localises:

``L0``
    Solute-vdW cavity, whose surface does not depend on the density.  The
    level-set response is absent, so these tests pin the electrostatic
    conventions (``phi``, ``w_xyz``, nuclear charges) on their own.
``L1``
    Isodensity cavity at a *fixed* density matrix, over three component sets.
    Only the ``lsf`` routes are new relative to L0.
``L2``
    Self-consistent solvated SCF and its total nuclear gradient.

Each layer also carries a negative control, because an FD test whose extra term
is numerically negligible passes while testing nothing.  The PV component is
what makes the level-set route dominant rather than a 7% correction: with PV
alone the surface charges vanish and the *entire* Fock matrix is the level-set
contraction.

Tests carry one marker per concern, so each is a separate meson target under
the ``moist_pyscf`` suite: ``host``, ``vdw``, ``isodensity`` (crossed with
``cpcm`` / ``pv`` / ``cpcm_pv``), ``conventions`` and ``scf``.  A failure names
the layer that broke.

Mutation testing -- injecting sign flips, dropped terms, factor errors and index
aliases into the analytic derivatives -- confirms this suite catches every such
error down to ~1 part in 10^6 of the level-set term; below that the injected
error falls under the finite-difference noise.
"""

import functools
import os
from dataclasses import dataclass

import numpy as np
import pytest

try:
    from pyscf import grad, gto, scf
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from .interface import (
    DensityResponse,
    GaussianMomentRequest,
    GaussianPotentialRequest,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    ModelComponentPV,
    Response,
)
from .pyscf import (
    CFC, DROP, ISwiG, Isodensity, SvdW,
    GaussianMoments, PySCFHost, PySCFSolvation, moist_for_scf,
)

#: Dielectric constant of water
EPSILON = 80.0
#: 10 GPa (atomic units)
PRESSURE = 1.0e10 / 2.9421015697e13
#: Lebedev order
NLEB = 50
#: Cavity projection tolerance
PROJ_TOL = 1e-13

#: FD step on the density matrix, in units of a unit-Frobenius-norm direction
STEP_DM = 1e-4
#: FD step on nuclear coordinates in bohr; matches test_helpers.f90's tuned value
STEP_R = 2.5e-4

#: Tolerances
REL_THR = 1e-9
ABS_THR = REL_THR / 10.0
#: Converged-SCF gradient tolerances
#:
#: The absolute floor is deliberately not ``SCF_REL_THR / 10``.  These
#: references are differences of total energies near -75 Ha, so each sample
#: carries about an ULP of noise, and ``fd4`` multiplies that by ``18 / (12 h)``
#: -- roughly 6e3 at this step.  The FD side therefore cannot be trusted below
#: ~1e-10 no matter how tightly the SCF converges: it moves by that much merely
#: from changing the OpenMP thread count, and neither a smaller nor a larger
#: step reduces it.  1e-10 would be sitting on that floor, so the check would
#: report the quadrature noise rather than the gradient.
SCF_REL_THR = 1e-10
SCF_ABS_THR = 1e-9

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

#: Component sets driven on the isodensity cavity.  ``pv`` is the sharpest of
#: the three: with no electrostatic component the surface charges are zero, so
#: the whole Fock matrix is the level-set contraction.
COMPONENTS = {
    "cpcm": lambda: [ModelComponentCPCM(EPSILON)],
    "pv": lambda: [ModelComponentPV(PRESSURE)],
    "cpcm+pv": lambda: [ModelComponentCPCM(EPSILON), ModelComponentPV(PRESSURE)],
}


CASE_PARAMS = [
    pytest.param(system, basis, id=f"{system}-{basis}") for system, basis in CASES
]
#: One marker per component set so each can be its own meson target.
COMPONENT_PARAMS = [
    pytest.param("cpcm", id="cpcm", marks=pytest.mark.cpcm),
    pytest.param("pv", id="pv", marks=pytest.mark.pv),
    pytest.param("cpcm+pv", id="cpcm_pv", marks=pytest.mark.cpcm_pv),
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


FD4_OFFSETS = (2, 1, -1, -2)


def answer_with_zeros(coupling, request, ngrid) -> None:
    """Answer the current request of ``coupling`` with zeros, whatever shape it wants.

    ``request`` is the snapshot the loop yielded for it.  Used by the negative
    controls, which need every *other* request satisfied so the failure they
    provoke is unambiguous.
    """
    if isinstance(request, GaussianMomentRequest):
        coupling.answer(
            gt=np.zeros(ngrid),
            pt=np.zeros((ngrid, 3)),
            mt=np.zeros((ngrid, 3, 3)),
            rt=np.zeros((ngrid, 3)),
        )
    elif isinstance(request, GaussianPotentialRequest):
        coupling.answer(**{name: np.zeros((ngrid, 3)) if name == "dphi_dr" else np.zeros(ngrid)
                           for name in request.missing})



def fd4(values, step: float) -> float:
    """4-point central difference from samples at ``(+2, +1, -1, -2) * step``."""
    fpp, fp, fm, fmm = values
    return (-fpp + 8.0 * fp - 8.0 * fm + fmm) / (12.0 * step)


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


def make_host(mol, *, dm):
    """Host bound to ``mol`` at fixed ``dm``, for the host-only tests."""
    host = PySCFHost(mol)
    host.dm = dm
    return host


def cavity_config(isodensity):
    """The DROP configuration of one layer: isodensity or solute-vdW."""
    if isodensity:
        return DROP(lsf=Isodensity(), nleb=NLEB, tolerance=PROJ_TOL)
    return DROP(lsf=SvdW(), nleb=NLEB)


def solve(mol, positions=None, *, dm, isodensity, components="cpcm"):
    """Driver for ``mol`` displaced to ``positions`` (bohr), evaluated at ``dm``."""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    solvation = PySCFSolvation(mol, cavity_config(isodensity), COMPONENTS[components]())
    solvation.evaluate(dm)
    return solvation


def solvated_rhf(mol, epsilon):
    """Converged solvated RHF on the isodensity cavity used by the L2 tier."""
    mean_field = moist_for_scf(
        scf.RHF(mol), cavity=cavity_config(True), components=[ModelComponentCPCM(epsilon)],
    )
    mean_field.conv_tol = 1e-13
    mean_field.conv_tol_grad = 1e-9
    mean_field.kernel()
    return mean_field


@pytest.mark.conventions
def test_pyscf_host_is_an_isodensity_cavity_source():
    """A host is the density source of an isodensity cavity built from it."""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    host = make_host(mol, dm=dm)

    cavity = cavity_config(True).build(source=host)
    cavity.update(host.structure())

    assert cavity.density_dependent
    assert cavity.ngrid > 0


@pytest.mark.isodensity
@pytest.mark.cpcm
def test_results_are_immutable_and_independent_of_the_input_density():
    """The driver's results are frozen values, not views of the caller's arrays."""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    density = np.array(dm, copy=True)
    solvation = PySCFSolvation(mol, cavity_config(True), COMPONENTS["cpcm"]())

    result = solvation.evaluate(density)
    gradient = solvation.gradient(density)
    density.fill(0.0)

    with pytest.raises(ValueError, match="WRITEABLE"):
        result.fock.setflags(write=True)
    with pytest.raises(ValueError, match="WRITEABLE"):
        solvation.density_matrix.setflags(write=True)
    assert solvation.evaluate(dm) is result
    np.testing.assert_allclose(solvation.gradient(dm), gradient)


def fd_density(mol, dm, direction, *, isodensity, components="cpcm"):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        solvation = solve(mol, dm=dm + offset * STEP_DM * direction,
                          isodensity=isodensity, components=components)
        samples.append(solvation.energy)
        grids.append(solvation.model.cavity.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, isodensity, components="cpcm"):
    """dE/dR along one cartesian coordinate at fixed density matrix."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        displaced = positions.copy()
        displaced.flat[index] += offset * STEP_R
        solvation = solve(mol, displaced, dm=dm, isodensity=isodensity, components=components)
        samples.append(solvation.energy)
        grids.append(solvation.model.cavity.ngrid)
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


# ----------------------------------------------------------------------
# L0 -- solute-vdW cavity: electrostatic conventions only
# ----------------------------------------------------------------------


@pytest.mark.vdw
@pytest.mark.parametrize(
    "cavity",
    [DROP(lsf=SvdW(), nleb=26), DROP(lsf=CFC(), nleb=26), ISwiG(nleb=26)],
    ids=("svdw-drop", "cfc-drop", "iswig"),
)
@pytest.mark.parametrize(
    "component_type",
    [ModelComponentCPCM, ModelComponentCOSMO],
    ids=("cpcm", "cosmo"),
)
def test_l0_pcm_components_and_cavity_types_share_the_pyscf_driver(cavity, component_type):
    """Every PCM/cavity combination runs through the same driver."""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    solvation = PySCFSolvation(mol, cavity, [component_type(EPSILON)])

    result = solvation.evaluate(dm)

    assert np.isfinite(result.energy)
    assert result.fock.shape == dm.shape
    assert solvation.gradient(dm).shape == (mol.natm, 3)
    assert solvation.model.cavity.ngrid > 0


@pytest.mark.vdw
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l0_fock_matches_fd(system, basis):
    """dE/dP through the surface potential alone, with a fixed surface."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    fock = solve(mol, dm=dm, isodensity=False).fock

    for direction in symmetric_directions(dm.shape[0], 2):
        numerical = fd_density(mol, dm, direction, isodensity=False)
        analytic = float(np.einsum("uv,uv->", fock, direction))
        assert deviation(analytic, numerical) <= 1.0


@pytest.mark.vdw
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l0_gradient_matches_fd(system, basis):
    """dE/dR at fixed P: moist's geometry terms plus the host's AO derivative."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    gradient = solve(mol, dm=dm, isodensity=False).gradient(dm)

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, positions, dm, index, isodensity=False)
        assert deviation(gradient.flatten(order="C")[index], numerical) <= 1.0


@pytest.mark.conventions
def test_l0_gradient_requires_position_weight():
    """Without ``w_xyz`` the gradient is refused rather than silently wrong.

    The total surface-position weight carries the whole surface-motion term
    of the gradient; read as zero it would leave a plausible, wrong gradient
    (only the direct nuclear term).  The request is mandatory, so the
    gradient is refused by name.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    model = solve(mol, dm=dm, isodensity=False).model

    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    # Answer everything except the position weight.
    for request in coupling:
        if not (isinstance(request, GaussianPotentialRequest) and "dphi_dr" in request.missing):
            answer_with_zeros(coupling, request, model.cavity.ngrid)
    with pytest.raises(RuntimeError, match="gaussian_potential"):
        _gradient(model, coupling)


# ----------------------------------------------------------------------
# L1 -- isodensity cavity at fixed density matrix
# ----------------------------------------------------------------------


@pytest.mark.host
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_density_callback_derivatives(system, basis):
    """The callback's own derivative orders are mutually consistent.

    Validates the Leibniz expansion and the PySCF derivative-component ordering
    independently of moist, so a failure here cannot be blamed on the cavity.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    point = mol.atom_coords()[0] + np.array([0.9, 0.4, 1.3])
    step = 1e-4
    _, drho, d2rho, d3rho = host.density(point, 3)

    def stencil(func, axis):
        samples = []
        for offset in FD4_OFFSETS:
            shifted = point.copy()
            shifted[axis] += offset * step
            samples.append(func(shifted))
        return fd4(samples, step)

    numerical_grad = np.array([stencil(lambda p: host.density(p, 1)[0], k) for k in range(3)])
    numerical_hess = np.array([stencil(lambda p: host.density(p, 1)[1], k) for k in range(3)])
    numerical_third = np.array([stencil(lambda p: host.density(p, 2)[2], k) for k in range(3)])

    assert np.abs(numerical_grad - drho).max() < 1e-9
    assert np.abs(numerical_hess.T - d2rho).max() < 1e-9
    assert np.abs(np.moveaxis(numerical_third, 0, -1) - d3rho).max() < 1e-8
    assert np.abs(d2rho - d2rho.T).max() < 1e-14


@pytest.mark.host
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_density_callback_is_finite(system, basis):
    """Pins the suppressed BLAS status flags in the callback as false positives.

    The products there raise spurious divide-by-zero and overflow, so the flags
    are ignored; this asserts the values really are finite and identical to a
    BLAS-free reference, both at a normal point and far into the tail where the
    density underflows.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    center = mol.atom_coords()[0]
    for offset in ([0.9, 0.4, 1.3], [40.0, 0.0, 0.0], [0.0, -60.0, 25.0]):
        point = center + np.array(offset)
        rho, drho, d2rho, d3rho = host.density(point, 3)
        assert np.isfinite(rho)
        for tensor in (drho, d2rho, d3rho):
            assert np.isfinite(tensor).all()

        ao = host._ao(point.reshape(1, 3), 3)[:, 0, :]
        reference = np.einsum("cu,uv->cv", ao, dm, optimize=False)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            product = ao @ dm
        # Not bitwise: BLAS and the einsum loop sum in different orders. Equal to
        # rounding is all that is needed to show the flags carry no information.
        assert np.isfinite(product).all()
        np.testing.assert_allclose(product, reference, rtol=5e-12, atol=1e-16)


@pytest.mark.isodensity
@pytest.mark.parametrize("components", COMPONENT_PARAMS)
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_fock_matches_fd(system, basis, components):
    """dE/dP with the surface following the density -- the headline Fock test."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    fock = solve(mol, dm=dm, isodensity=True, components=components).fock

    for direction in symmetric_directions(dm.shape[0], 2):
        numerical = fd_density(
            mol, dm, direction, isodensity=True, components=components
        )
        analytic = float(np.einsum("uv,uv->", fock, direction))
        assert deviation(analytic, numerical) <= 1.0


@pytest.mark.isodensity
@pytest.mark.parametrize("components", COMPONENT_PARAMS)
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_gradient_matches_fd(system, basis, components):
    """dE/dR at fixed P, including the level set's own basis-center derivative."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    gradient = solve(mol, dm=dm, isodensity=True, components=components).gradient(dm)

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(
            mol, positions, dm, index, isodensity=True, components=components
        )
        assert deviation(gradient.flatten(order="C")[index], numerical) <= 1.0


@pytest.mark.conventions
def test_the_driver_stops_on_an_item_it_cannot_contract():
    """An item the driver does not know raises instead of dropping its term.

    Unlike an unanswered request, which ``get_*`` reports, a skipped response
    item fails nowhere: the Fock matrix and the gradient would silently lack
    its contribution.
    """

    class UnknownItem:
        name = "unknown"

    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    solvation = solve(mol, dm=dm, isodensity=False)
    response = Response([*solvation.response, UnknownItem()])
    with pytest.raises(NotImplementedError, match="Unsupported response item 'unknown'"):
        solvation._fock(response)
    # The gradient walks the response of the cached evaluation at this density.
    solvation._response = response
    with pytest.raises(NotImplementedError, match="Unsupported response item 'unknown'"):
        solvation.gradient_channels(dm)


@pytest.mark.conventions
def test_pv_alone_makes_the_fock_purely_level_set():
    """With no electrostatic component the whole Fock is the ``lsf`` contraction.

    A pv-only model produces no electrostatic channel at all, so nothing
    survives except the cavity-shape response -- the sharpest available
    isolation of the weights.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    solvation = solve(mol, dm=dm, isodensity=True, components="pv")
    # The walk meets the density item alone: no potential adjoint, no amplitudes.
    (density,) = list(solvation.response)
    assert isinstance(density, DensityResponse)
    full = solvation.fock
    np.testing.assert_allclose(
        full, solvation.host._fock_lsf(solvation.model.cavity.xyz, density), atol=1e-14
    )

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction, isodensity=True, components="pv")
    analytic = float(np.einsum("uv,uv->", full, direction))
    assert abs(analytic) > MIN_SIGNAL
    assert deviation(analytic, numerical) <= 1.0


@pytest.mark.conventions
@pytest.mark.parametrize("components", ["cpcm", "cpcm+pv"])
def test_l1_fock_requires_lsf_term(components):
    """Dropping the level-set response must break the Fock test outright."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    without = solve(mol, dm=dm, isodensity=True, components=components).frozen_fock()

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction, isodensity=True, components=components)
    analytic = float(np.einsum("uv,uv->", without, direction))
    assert deviation(analytic, numerical) > VACUITY_FACTOR


@pytest.mark.conventions
@pytest.mark.parametrize("components", ["cpcm", "cpcm+pv"])
def test_l1_gradient_requires_lsf_term(components):
    """The isodensity level set reports zero nuclear partials by construction.

    moist therefore returns a gradient that is missing the density's own
    dependence on the nuclei; the host has to add it.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    channels = solve(mol, dm=dm, isodensity=True, components=components).gradient_channels(dm)
    assert np.abs(channels["density"]).max() > MIN_SIGNAL
    starved = channels["model"] + channels["potential"]

    numerical = fd_position(
        mol, positions, dm, 2, isodensity=True, components=components
    )
    assert deviation(starved.flatten(order="C")[2], numerical) > VACUITY_FACTOR


@pytest.mark.conventions
def test_l1_potential_requires_point_potentials():
    """``w_xyz`` carries the dominant part of the cavity response.

    When the density changes the  grid points move and ``phi(r_i)`` moves with them.
    moist cannot see that route -- ``phi`` is the host's function -- so it has to
    arrive as ``w_xyz`` before the potential is read.  The request is
    mandatory for a density-dependent cavity: omitting it is refused by name
    rather than returning ``lsf`` weights that look healthy and are wrong.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    model = solve(mol, dm=dm, isodensity=True).model

    coupling = model.new_coupling()
    model.prepare_response(coupling)
    for request in coupling:
        if not (isinstance(request, GaussianPotentialRequest) and "dphi_dr" in request.missing):
            answer_with_zeros(coupling, request, model.cavity.ngrid)
    with pytest.raises(RuntimeError, match="gaussian_potential"):
        model.get_response(coupling)


@pytest.mark.conventions
def test_the_host_is_asked_for_the_potential_once_per_evaluation(monkeypatch):
    """``phi`` is not charge-dependent, so it is missing in the energy phase only.

    The pre-protocol exchange supplied it twice -- once bare to obtain the
    charges, then again alongside the position weights -- because the host had
    to drive moist's ordering itself.  The phase loop owns that ordering now,
    so the potential is asked for exactly once and only ``w_xyz`` comes back
    for the response and the gradient.
    """
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    solvation = PySCFSolvation(mol, cavity_config(True), COMPONENTS["cpcm"]())
    host = solvation.host

    potential_calls = []
    weight_calls = []
    surface_potential = host.surface_potential
    spatial_derivative = host._grad_phi_elec

    def counted_potential(coords, xi=None):
        potential_calls.append(1)
        return surface_potential(coords, xi)

    def counted_weights(coords, xi=None):
        weight_calls.append(1)
        return spatial_derivative(coords, xi)

    monkeypatch.setattr(host, "surface_potential", counted_potential)
    monkeypatch.setattr(host, "_grad_phi_elec", counted_weights)

    solvation.evaluate(dm)

    assert potential_calls == [1]
    assert weight_calls == [1]
    solvation.gradient(dm)
    # Raw derivatives remain valid when the gradient phase follows response.
    assert potential_calls == [1]
    assert len(weight_calls) == 1


@pytest.mark.conventions
def test_a_model_without_a_moment_request_builds_no_gaussian_integrals(monkeypatch):
    """A model that does not ask for moments never makes the host form them.

    The Gaussian moments are dense three-centre AO integrals.  A host cannot
    know whether they are wanted, so it used to be told by a hand-written table
    on the Python side; now the model's own declaration decides, and a CPCM
    model simply never produces the request that would trigger them.
    """
    builds = []
    build_integrals = GaussianMoments._build_integrals

    def counted(self):
        builds.append(1)
        return build_integrals(self)

    monkeypatch.setattr(GaussianMoments, "_build_integrals", counted)

    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)

    electrostatic = PySCFSolvation(mol, cavity_config(True), COMPONENTS["cpcm"]())
    electrostatic.evaluate(dm)
    assert builds == []
    assert electrostatic.moments is None

    # Vacuous unless the same driver does build them when a component asks.
    pressurised = PySCFSolvation(mol, cavity_config(True), [ModelComponentGOSTSHYP(PRESSURE)])
    pressurised.evaluate(dm)
    assert builds == [1]


@pytest.mark.conventions
def test_gradient_path_reads_host_surface_weights():
    """The gradient contracts the same total ``w_xyz`` the response uses.

    There is one surface-position weight, ``q_i grad phi_total(r_i)``, read by
    both paths; moist adds no nuclear field of its own to it.  Scaling it
    therefore scales the surface-motion term of the gradient: the difference
    between the two gradients below is linear in the factor and nonzero.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    solvation = solve(mol, dm=dm, isodensity=True)
    host, model = solvation.host, solvation.model
    coords, xi = model.cavity.xyz, model.cavity.xi0

    coupling = model.new_coupling()
    gradients = []
    for factor in (1.0, 2.0, 3.0):
        model.prepare_energy(coupling)
        for request in coupling:
            if isinstance(request, GaussianPotentialRequest):
                coupling.answer(phi=host.surface_potential(coords, xi))
        _energy(model, coupling)
        model.prepare_gradient(coupling)
        for request in coupling:
            if isinstance(request, GaussianPotentialRequest) and "dphi_dr" in request.missing:
                coupling.answer(
                    dphi_dr=factor * host.surface_potential_gradient(coords, xi),
                    dphi_dxi=host._dphi_dxi(coords, xi),
                )
        gradients.append(_gradient(model, coupling)[0])
    step = gradients[1] - gradients[0]
    assert np.max(np.abs(step)) > 1.0e-6
    np.testing.assert_allclose(
        gradients[2] - gradients[1], step, rtol=1.0e-9, atol=1.0e-12
    )


# ----------------------------------------------------------------------
# L2 -- self-consistent solvated SCF (slow tier)
# ----------------------------------------------------------------------


@pytest.mark.scf
def test_l2_scf_converges_and_stabilises():
    """A solvated SCF converges, is stabilising, and reduces to gas phase."""
    mol = molecule(*PRIMARY_CASE)
    gas = scf.RHF(mol).run()
    solvated = solvated_rhf(mol, EPSILON)
    assert solvated.converged
    assert solvated.e_tot < gas.e_tot

    vacuum = solvated_rhf(mol, 1.0)
    assert vacuum.e_tot == pytest.approx(gas.e_tot, abs=1e-9)


@pytest.mark.scf
def test_l2_total_gradient_matches_fd():
    """Total solvated SCF energy gradient against FD of the converged energy.

    The density response drops out at convergence, so the analytic gradient is
    the ordinary RHF gradient built from the solvated orbitals plus the explicit
    solvation terms evaluated at the converged density.
    """
    mol = molecule(*PRIMARY_CASE)
    positions = mol.atom_coords()
    solvated = solvated_rhf(mol, EPSILON)
    dm = solvated.make_rdm1()

    analytic = grad.RHF(solvated).kernel() + solve(mol, dm=dm, isodensity=True).gradient(dm)

    for index in sampled_coordinates(mol.natm)[:3]:
        samples = []
        for offset in FD4_OFFSETS:
            displaced = positions.copy()
            displaced.flat[index] += offset * STEP_R
            moved = mol.set_geom_(displaced, unit="Bohr", inplace=False)
            samples.append(solvated_rhf(moved, EPSILON).e_tot)
        numerical = fd4(samples, STEP_R)
        assert deviation(analytic.flatten(order="C")[index], numerical,
                 thr_abs=SCF_ABS_THR, thr_rel=SCF_REL_THR) <= 1.0

    # The solvated gradient must actually differ from the gas-phase one, or the
    # solvation terms above are not being exercised.
    gas_gradient = grad.RHF(scf.RHF(mol).run()).kernel()
    assert np.abs(analytic - gas_gradient).max() > MIN_SIGNAL


@pytest.mark.host
@pytest.mark.parametrize("cart", [False, True])
def test_gaussian_raw_derivatives(cart):
    """Spatial and inverse-length derivatives include the coincident-point limit."""
    mol = molecule("water", "def2-svp").copy()
    mol.cart = cart
    host = make_host(mol, dm=scf.RHF(mol).run(verbose=0).make_rdm1())
    coords = np.vstack([mol.atom_coords()[0], [0.2, -0.3, 0.1], [2.0, 1.0, 1.0]])
    xi = np.array([0.5, 1.2, 2.0])
    step = 1e-5
    numerical = (host.surface_potential(coords, xi + step)
                 - host.surface_potential(coords, xi - step)) / (2 * step)
    np.testing.assert_allclose(host._dphi_dxi(coords, xi), numerical, atol=1e-8, rtol=1e-8)
    spatial = host._grad_phi_nuc(coords, xi) + host._grad_phi_elec(coords, xi)
    for axis in range(3):
        shift = np.zeros_like(coords)
        shift[:, axis] = step
        numerical = (host.surface_potential(coords + shift, xi)
                     - host.surface_potential(coords - shift, xi)) / (2 * step)
        np.testing.assert_allclose(spatial[axis], numerical, atol=1e-8, rtol=1e-8)


@pytest.mark.conventions
def test_gaussian_pcm_matches_pyscf_on_the_same_cavity():
    """Independent PCM matrix, RHS, energy and Fock on identical surface inputs."""
    from pyscf.solvent import pcm
    from . import library

    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    solvation = PySCFSolvation(mol, DROP(lsf=SvdW(), nleb=50), [ModelComponentCPCM(32.0)])
    result = solvation.evaluate(dm)
    host, model = solvation.host, solvation.model
    coords = model.cavity.xyz
    native_matrix, xi = library.assemble_drop_amat(model.cavity._handle)
    _, switch = library.get_cavity_gaussian(model.cavity._handle)
    reference = pcm.PCM(mol)
    reference.surface = dict(grid_coords=coords, charge_exp=xi, switch_fun=switch,
                             norm_vec=np.zeros_like(coords), R_vdw=np.ones(len(xi)))
    _, matrix = pcm.get_D_S(reference.surface, with_D=False)
    np.testing.assert_allclose(native_matrix, matrix, atol=1e-12, rtol=1e-12)
    auxiliary = gto.fakemol_for_charges(coords, expnt=xi**2)
    nuclei = gto.fakemol_for_charges(mol.atom_coords())
    reference.v_grids_n = mol.atom_charges() @ gto.mole.intor_cross(
        mol._add_suffix("int2c2e"), nuclei, auxiliary)
    reference._intermediates = {"K": matrix, "R": -(31.0 / 32.0) * np.eye(len(xi))}
    energy, fock = reference._get_vind(dm)
    np.testing.assert_allclose(result.energy, energy, atol=1e-12, rtol=1e-10)
    np.testing.assert_allclose(result.fock, fock, atol=1e-12, rtol=1e-10)
    phi = reference.v_grids_n - reference._get_v(dm[None, :, :])[0]
    np.testing.assert_allclose(host.surface_potential(coords, xi), phi, atol=1e-10, rtol=1e-10)


@pytest.mark.conventions
def test_point_gaussian_mismatch_decreases_with_grid_order():
    """The operator difference is bounded and decreases with discretization error."""
    from . import library

    mol, dm = molecule("water", "def2-svp"), reference_density("water", "def2-svp")
    differences = []
    for nleb in (50, 194, 770):
        solvation = PySCFSolvation(mol, DROP(lsf=SvdW(), nleb=nleb), [ModelComponentCPCM(78.3553)])
        result = solvation.evaluate(dm)
        matrix, _ = library.assemble_drop_amat(solvation.model.cavity._handle)
        point_phi = solvation.host.surface_potential(solvation.model.cavity.xyz)
        q = -(1.0 - 1.0 / 78.3553) * np.linalg.solve(matrix, point_phi)
        point_energy = 0.5 * q @ point_phi
        differences.append(abs(result.energy - point_energy) * 627.509474)
    assert differences[0] > differences[1] > differences[2] > 0
    assert differences[1] < 0.05  # kcal/mol, default 194-point grid
    # A broad bound on inverse-grid-order convergence, independent of last bits.
    assert 3.0 < differences[0] / differences[2] < 40.0


@pytest.mark.host
def test_ecp_molecule_is_rejected_before_element_lookup():
    mol = gto.M(atom="I 0 0 0", basis="def2-svp", ecp="def2-svp", spin=1, verbose=0)
    host = PySCFHost(mol)
    with pytest.raises(ValueError, match="effective core potentials"):
        host.structure()


def _energy(model, coupling):
    energy = np.array(0.0)
    model.get_energy(coupling, energy)
    return float(energy)


def _gradient(model, coupling):
    gradient = np.zeros((model._natoms, 3))
    response = model.get_gradient(coupling, gradient)
    return gradient, response
