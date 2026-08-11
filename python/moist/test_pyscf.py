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
    conventions (``phi``, ``qefield``, nuclear charges) on their own.
``L1``
    Isodensity cavity at a *fixed* density matrix, over three component sets.
    Only the ``w_lsf`` routes are new relative to L0.
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
from dataclasses import dataclass

import numpy as np
import pytest

try:
    # Import what is actually used rather than the top-level package alone, and
    # skip on any ImportError: pyscf imports fine while pyscf.scf fails when a
    # broken h5py shadows the HDF5 runtime, and importorskip re-raises for an
    # ImportError raised by a transitive dependency rather than the named module.
    from pyscf import grad, gto, scf
except ImportError as exc:  # pragma: no cover - optional dependency
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from .interface import CPCM, PV, DROPCavity, GeneralSolvationModel, IsodensityDROPCavity
from .pyscf import PySCFHost, PySCFIsodensityHost, solvated_rhf

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
SCF_REL_THR = 1e-9
SCF_ABS_THR = SCF_REL_THR / 10.0

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
    "cpcm": lambda: [CPCM(EPSILON)],
    "pv": lambda: [PV(PRESSURE)],
    "cpcm+pv": lambda: [CPCM(EPSILON), PV(PRESSURE)],
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


def make_host(mol, positions=None, *, dm):
    """Host bound to ``mol`` displaced to ``positions`` (bohr), at fixed ``dm``."""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    host = PySCFHost(mol)
    host.dm = dm
    return host


@pytest.mark.host
def test_pyscf_host_is_an_isodensity_cavity_source():
    mol = molecule(*PRIMARY_CASE)
    host = PySCFHost(mol)

    cavity = IsodensityDROPCavity(host, nleb=NLEB, tolerance=PROJ_TOL)

    assert isinstance(cavity, IsodensityDROPCavity)
    with pytest.deprecated_call(match="PySCFHost"):
        legacy = PySCFIsodensityHost(mol)
    with pytest.deprecated_call(match="IsodensityDROPCavity"):
        compatibility_cavity = host.make_cavity(nleb=NLEB)
    assert isinstance(legacy, PySCFHost)
    assert isinstance(compatibility_cavity, IsodensityDROPCavity)


def solve(host, *, isodensity, components="cpcm"):
    """Build a cavity plus components and evaluate one coherent coupling."""
    if isodensity:
        cavity = IsodensityDROPCavity(host, nleb=NLEB, tolerance=PROJ_TOL)
    else:
        cavity = DROPCavity(nleb=NLEB)
    model = GeneralSolvationModel(cavity, COMPONENTS[components]())
    result = model.evaluate(coupling=host.coupling(host.dm))
    return result.energy, result.potential, result.cavity.xyz.T, model


@pytest.mark.host
def test_evaluation_exposes_complete_pyscf_results():
    """The evaluation owns the host Fock and complete nuclear gradient."""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    host = make_host(mol, dm=dm)
    model = GeneralSolvationModel(
        IsodensityDROPCavity(host, nleb=NLEB, tolerance=PROJ_TOL),
        COMPONENTS["cpcm"](),
    )

    density = np.array(dm, copy=True)
    coupling = host.coupling(density)
    result = model.evaluate(coupling=coupling)
    coords = result.cavity.xyz.T

    np.testing.assert_allclose(
        result.fock,
        host.fock(coords, result.potential, include_lsf=True),
    )
    with pytest.raises(ValueError, match="WRITEABLE"):
        coupling.density_matrix.setflags(write=True)
    with pytest.raises(AttributeError):
        coupling.density_matrix = np.zeros_like(density)
    density.fill(0.0)
    np.testing.assert_allclose(
        result.gradient,
        model.gradient() + host.gradient(coords, result.potential, include_lsf=True),
    )


def fd_density(mol, dm, direction, *, isodensity, components="cpcm"):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        host = make_host(mol, dm=dm + offset * STEP_DM * direction)
        energy, _, _, model = solve(host, isodensity=isodensity, components=components)
        samples.append(energy)
        grids.append(model.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, isodensity, components="cpcm"):
    """dE/dR along one cartesian coordinate at fixed density matrix."""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        displaced = positions.copy()
        displaced.flat[index] += offset * STEP_R
        host = make_host(mol, displaced, dm=dm)
        energy, _, _, model = solve(host, isodensity=isodensity, components=components)
        samples.append(energy)
        grids.append(model.ngrid)
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
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l0_fock_matches_fd(system, basis):
    """dE/dP through the surface potential alone, with a fixed surface."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, potential, coords, _ = solve(host, isodensity=False)
    fock = host.fock(coords, potential, include_lsf=False)

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
    host = make_host(mol, dm=dm)
    _, potential, coords, model = solve(host, isodensity=False)
    gradient = model.get_gradient(mol.natm) + host.gradient(
        coords, potential, include_lsf=False
    )

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(mol, positions, dm, index, isodensity=False)
        assert deviation(gradient.flatten(order="F")[index], numerical) <= 1.0


@pytest.mark.conventions
def test_l0_gradient_requires_qefield():
    """Without ``qefield`` moist keeps only the nuclear half of the surface motion.

    The omission is silent -- an unsupplied channel is treated as zero -- so this
    guards the failure mode rather than the formula.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    positions = mol.atom_coords()
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=False)

    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)  # no qefield
    starved = model.get_gradient(mol.natm) + host.gradient(
        coords, model.get_potential(), include_lsf=False
    )

    numerical = fd_position(mol, positions, dm, 2, isodensity=False)
    assert deviation(starved.flatten(order="F")[2], numerical) > VACUITY_FACTOR


# ----------------------------------------------------------------------
# L1 -- isodensity cavity at fixed density matrix
# ----------------------------------------------------------------------


@pytest.mark.host
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_lsf_callback_derivatives(system, basis):
    """The callback's own derivative orders are mutually consistent.

    Validates the Leibniz expansion and the PySCF derivative-component ordering
    independently of moist, so a failure here cannot be blamed on the cavity.
    """
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    point = mol.atom_coords()[0] + np.array([0.9, 0.4, 1.3])
    step = 1e-4
    _, grad_, hess, third = host.lsf(point, 3)

    def stencil(func, axis):
        samples = []
        for offset in FD4_OFFSETS:
            shifted = point.copy()
            shifted[axis] += offset * step
            samples.append(func(shifted))
        return fd4(samples, step)

    numerical_grad = np.array([stencil(lambda p: host.lsf(p, 1)[0], k) for k in range(3)])
    numerical_hess = np.array([stencil(lambda p: host.lsf(p, 1)[1], k) for k in range(3)])
    numerical_third = np.array([stencil(lambda p: host.lsf(p, 2)[2], k) for k in range(3)])

    assert np.abs(numerical_grad - grad_).max() < 1e-9
    assert np.abs(numerical_hess.T - hess).max() < 1e-9
    assert np.abs(np.moveaxis(numerical_third, 0, -1) - third).max() < 1e-8
    assert np.abs(hess - hess.T).max() < 1e-14


@pytest.mark.host
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_lsf_callback_is_finite(system, basis):
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
        value, grad_, hess, third = host.lsf(point, 3)
        assert np.isfinite(value)
        for tensor in (grad_, hess, third):
            assert np.isfinite(tensor).all()

        ao = host._ao(point.reshape(1, 3), 3)[:, 0, :]
        reference = np.einsum("cu,uv->cv", ao, dm, optimize=False)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            product = ao @ dm
        # Not bitwise: BLAS and the einsum loop sum in different orders. Equal to
        # rounding is all that is needed to show the flags carry no information.
        assert np.isfinite(product).all()
        np.testing.assert_allclose(product, reference, rtol=1e-12, atol=0.0)


@pytest.mark.isodensity
@pytest.mark.parametrize("components", COMPONENT_PARAMS)
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_fock_matches_fd(system, basis, components):
    """dE/dP with the surface following the density -- the headline Fock test."""
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, potential, coords, _ = solve(host, isodensity=True, components=components)
    fock = host.fock(coords, potential, include_lsf=True)

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
    host = make_host(mol, dm=dm)
    _, potential, coords, model = solve(host, isodensity=True, components=components)
    gradient = model.get_gradient(mol.natm) + host.gradient(
        coords, potential, include_lsf=True
    )

    for index in sampled_coordinates(mol.natm):
        numerical = fd_position(
            mol, positions, dm, index, isodensity=True, components=components
        )
        assert deviation(gradient.flatten(order="F")[index], numerical) <= 1.0


@pytest.mark.conventions
def test_pv_alone_makes_the_fock_purely_level_set():
    """With no electrostatic component the whole Fock is the ``w_lsf`` contraction.

    The surface charges are identically zero, so nothing survives except the
    cavity-shape response -- the sharpest available isolation of the weights.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, potential, coords, _ = solve(host, isodensity=True, components="pv")

    np.testing.assert_array_equal(potential.w_umol, np.zeros_like(potential.w_umol))
    full = host.fock(coords, potential, include_lsf=True)
    np.testing.assert_allclose(full, host._fock_lsf(coords, potential), atol=1e-14)

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
    host = make_host(mol, dm=dm)
    _, potential, coords, _ = solve(host, isodensity=True, components=components)
    without = host.fock(coords, potential, include_lsf=False)

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
    host = make_host(mol, dm=dm)
    _, potential, coords, model = solve(host, isodensity=True, components=components)
    starved = model.get_gradient(mol.natm) + host.gradient(
        coords, potential, include_lsf=False
    )

    numerical = fd_position(
        mol, positions, dm, 2, isodensity=True, components=components
    )
    assert deviation(starved.flatten(order="F")[2], numerical) > VACUITY_FACTOR


@pytest.mark.conventions
def test_l1_potential_requires_surface_position_weights():
    """``w_xyz`` carries the dominant part of the cavity response.

    When the density changes the  grid points move and ``phi(r_i)`` moves with them.
    moist cannot see that route -- ``phi`` is the host's function -- so it has to
    arrive as ``w_xyz`` before the potential is read.  Omitting it returns
    ``w_lsf`` weights that look perfectly healthy and are badly wrong.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=True)

    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)  # no w_xyz
    starved = host.fock(coords, model.get_potential(), include_lsf=True)

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction, isodensity=True)
    analytic = float(np.einsum("uv,uv->", starved, direction))
    assert deviation(analytic, numerical) > VACUITY_FACTOR


@pytest.mark.conventions
def test_cpcm_leaves_qmol_untouched():
    """CPCM writes no normal-derivative trace, so ``w_qmol`` must stay zero."""
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, potential, _, _ = solve(host, isodensity=True)
    np.testing.assert_array_equal(potential.w_qmol, np.zeros_like(potential.w_qmol))


@pytest.mark.conventions
def test_gradient_path_ignores_host_surface_weights():
    """Documents an asymmetry between the potential and gradient paths.

    ``w_xyz`` is read when the potential is assembled but dropped when the
    gradient is, so scaling it by a thousand moves the gradient by exactly zero.
    That is what makes it safe for :meth:`PySCFHost.solve` to supply
    ``w_xyz`` and ``qefield`` together: were the gradient path to start reading
    ``w_xyz``, the surface-motion term would be counted twice and
    ``test_l1_gradient_matches_fd`` would begin to fail.
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=True)
    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)
    charges, _ = model.get_trace_potential()

    gradients = []
    for factor in (0.0, 1000.0):
        model.supply_electrostatics(
            phi,
            w_xyz=factor * host.surface_position_weights(coords, charges),
            qefield=host.qefield(coords, charges),
        )
        gradients.append(model.get_gradient(mol.natm))
    np.testing.assert_array_equal(gradients[0], gradients[1])


# ----------------------------------------------------------------------
# L2 -- self-consistent solvated SCF (slow tier)
# ----------------------------------------------------------------------


@pytest.mark.scf
def test_l2_scf_converges_and_stabilises():
    """A solvated SCF converges, is stabilising, and reduces to gas phase."""
    mol = molecule(*PRIMARY_CASE)
    gas = scf.RHF(mol).run()
    solvated = solvated_rhf(mol, EPSILON, nleb=NLEB, tolerance=PROJ_TOL)
    assert solvated.converged
    assert solvated.e_tot < gas.e_tot

    vacuum = solvated_rhf(mol, 1.0, nleb=NLEB, tolerance=PROJ_TOL)
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
    solvated = solvated_rhf(mol, EPSILON, nleb=NLEB, tolerance=PROJ_TOL)
    dm = solvated.make_rdm1()

    host = make_host(mol, dm=dm)
    _, potential, coords, model = solve(host, isodensity=True)
    analytic = (
        grad.RHF(solvated).kernel().T
        + model.get_gradient(mol.natm)
        + host.gradient(coords, potential, include_lsf=True)
    )

    for index in sampled_coordinates(mol.natm)[:3]:
        samples = []
        for offset in FD4_OFFSETS:
            displaced = positions.copy()
            displaced.flat[index] += offset * STEP_R
            moved = mol.set_geom_(displaced, unit="Bohr", inplace=False)
            samples.append(
                solvated_rhf(moved, EPSILON, nleb=NLEB, tolerance=PROJ_TOL).e_tot
            )
        numerical = fd4(samples, STEP_R)
        assert deviation(analytic.flatten(order="F")[index], numerical,
                 thr_abs=SCF_ABS_THR, thr_rel=SCF_REL_THR) <= 1.0

    # The solvated gradient must actually differ from the gas-phase one, or the
    # solvation terms above are not being exercised.
    gas_gradient = grad.RHF(scf.RHF(mol).run()).kernel().T
    assert np.abs(analytic - gas_gradient).max() > MIN_SIGNAL
