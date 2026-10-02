"""End-to-end tests of the moist/host chain rule against PySCF

moist hands out adjoint weights, the host finishes the chain rule with its own
density derivatives: isodensity-cavity correctness spans a language boundary.
Every analytic quantity is checked against a finite difference of the energy
moist itself returns.

Layers, so a failure localises:

``L0``
    Solute-vdW cavity, surface independent of the density
    - no level-set response: pins the electrostatic conventions (``phi``,
      ``qefield``, nuclear charges) alone
``L1``
    Isodensity cavity at a *fixed* density matrix, over three component sets
    - only the ``lsf`` routes are new relative to L0
``L2``
    Self-consistent solvated SCF and its total nuclear gradient

Each layer carries a negative control: an FD test whose extra term is
numerically negligible passes while testing nothing.
- PV makes the level-set route dominant, not a 7% correction
- with PV alone the surface charges vanish and the *entire* Fock matrix is the
  level-set contraction

One marker per concern, each a separate meson target of the ``moist_pyscf``
suite: ``host``, ``vdw``, ``isodensity`` (crossed with ``cpcm`` / ``pv`` /
``cpcm_pv``), ``conventions``, ``scf``. A failure names the layer that broke.

Mutation testing (sign flips, dropped terms, factor errors, index aliases
injected into the analytic derivatives) is caught down to ~1 part in 10^6 of
the level-set term
- below that the injected error falls under the FD noise
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
    CavityDROP,
    CavityDROPCFC,
    CavityDROPIsodensityCallback,
    CavityDROPIsodensityInternal,
    CavityDROPSvdW,
    CavityISwiG,
    GeneralSolvationModel,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentPV,
)
from .pyscf import PySCFHost, PySCFIsodensityHost, solvated_rhf

#: Dielectric constant of water
EPSILON = 80.0
#: 10 GPa (atomic units)
PRESSURE = 1.0e10 / 2.9421015697e13
#: Lebedev order
NLEB = 50
#: Cavity projection tolerance
PROJ_TOL = 1e-13

#: FD step on the density matrix, along a unit-Frobenius-norm direction
STEP_DM = 1e-4
#: FD step on nuclear coordinates, bohr; tuned value of test_helpers.f90
STEP_R = 2.5e-4

#: Tolerances
REL_THR = 1e-9
ABS_THR = REL_THR / 10.0
#: Converged-SCF gradient tolerances
#:
#: Absolute floor deliberately not ``SCF_REL_THR / 10``
#: - references are total-energy differences near -75 Ha, ~1 ULP noise per sample
#: - ``fd4`` multiplies it by ``18 / (12 h)``, ~6e3 at this step
#: - FD side untrustworthy below ~1e-10 however tight the SCF
#: - it moves by that much with the OpenMP thread count alone
#: - neither a smaller nor a larger step reduces it
#: - 1e-10 would sit on that floor and report quadrature noise, not the gradient
SCF_REL_THR = 1e-10
SCF_ABS_THR = 1e-9

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
#: Case of tests that pin a convention once, not sweep
PRIMARY_CASE = CASES[0]

#: Component sets driven on the isodensity cavity
#: - ``pv`` sharpest: no electrostatic component, so the surface charges are zero
#:   and the whole Fock matrix is the level-set contraction
COMPONENTS = {
    "cpcm": lambda: [ModelComponentCPCM(EPSILON)],
    "pv": lambda: [ModelComponentPV(PRESSURE)],
    "cpcm+pv": lambda: [ModelComponentCPCM(EPSILON), ModelComponentPV(PRESSURE)],
}


CASE_PARAMS = [
    pytest.param(system, basis, id=f"{system}-{basis}") for system, basis in CASES
]
#: One marker per component set, each its own meson target
COMPONENT_PARAMS = [
    pytest.param("cpcm", id="cpcm", marks=pytest.mark.cpcm),
    pytest.param("pv", id="pv", marks=pytest.mark.pv),
    pytest.param("cpcm+pv", id="cpcm_pv", marks=pytest.mark.cpcm_pv),
]


def deviation(actual, reference, *, thr_abs=None, thr_rel=None) -> float:
    """Deviation in units of the tolerance, ``<= 1`` passes

    Thresholds combine as in test-drive's ``check(..., thr_abs=, thr_rel=)``:
    ``max(thr_abs, thr_rel*|ref|)``
    - absolute floor covers references near zero, relative bound scales with size
    - a ratio, so a failure message says how many tolerances were missed
    """
    thr_abs = ABS_THR if thr_abs is None else thr_abs
    thr_rel = REL_THR if thr_rel is None else thr_rel
    return abs(actual - reference) / max(thr_abs, thr_rel * abs(reference))


FD4_OFFSETS = (2, 1, -1, -2)


def fd4(values, step: float) -> float:
    """4-point central difference from samples at ``(+2, +1, -1, -2) * step``"""
    fpp, fp, fm, fmm = values
    return (-fpp + 8.0 * fp - 8.0 * fm + fmm) / (12.0 * step)


# ----------------------------------------------------------------------
# System construction (cached, reused by every parametrised test)
# ----------------------------------------------------------------------


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


def make_host(mol, positions=None, *, dm):
    """Host bound to ``mol`` displaced to ``positions`` (bohr), at fixed ``dm``"""
    if positions is not None:
        mol = mol.set_geom_(positions, unit="Bohr", inplace=False)
    host = PySCFHost(mol)
    host.dm = dm
    return host


@pytest.mark.conventions
def test_pyscf_host_is_an_isodensity_cavity_source():
    mol = molecule(*PRIMARY_CASE)
    host = PySCFHost(mol)

    cavity = CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL)

    assert isinstance(cavity, CavityDROPIsodensityCallback)
    with pytest.deprecated_call(match="PySCFHost"):
        legacy = PySCFIsodensityHost(mol)
    with pytest.deprecated_call(match="CavityDROPIsodensityCallback"):
        compatibility_cavity = host.make_cavity(nleb=NLEB)
    assert isinstance(legacy, PySCFHost)
    assert isinstance(compatibility_cavity, CavityDROPIsodensityCallback)


def solve(host, *, isodensity, components="cpcm"):
    """Build a cavity plus components and evaluate one coherent coupling"""
    if isodensity:
        cavity = CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL)
    else:
        cavity = CavityDROP(nleb=NLEB)
    model = GeneralSolvationModel(cavity, COMPONENTS[components]())
    result = model.evaluate(coupling=host.coupling(host.dm))
    return result.energy, result.response, result.cavity.xyz.T, model


@pytest.mark.isodensity
@pytest.mark.cpcm
def test_evaluation_exposes_complete_pyscf_results():
    """Evaluation carries the host Fock and complete nuclear gradient"""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    host = make_host(mol, dm=dm)
    model = GeneralSolvationModel(
        CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL),
        COMPONENTS["cpcm"](),
    )

    density = np.array(dm, copy=True)
    coupling = host.coupling(density)
    result = model.evaluate(coupling=coupling)
    coords = result.cavity.xyz.T

    np.testing.assert_allclose(
        result.fock,
        host.fock(coords, result.response, include_lsf=True),
    )
    with pytest.raises(ValueError, match="WRITEABLE"):
        coupling.density_matrix.setflags(write=True)
    with pytest.raises(AttributeError):
        coupling.density_matrix = np.zeros_like(density)
    density.fill(0.0)
    np.testing.assert_allclose(
        result.gradient,
        model.gradient() + host.gradient(coords, result.response, include_lsf=True),
    )


def fd_density(mol, dm, direction, *, isodensity, components="cpcm"):
    """dE/dt along ``dm + t * direction``, rebuilding the cavity each sample"""
    samples, grids = [], []
    for offset in FD4_OFFSETS:
        host = make_host(mol, dm=dm + offset * STEP_DM * direction)
        energy, _, _, model = solve(host, isodensity=isodensity, components=components)
        samples.append(energy)
        grids.append(model.ngrid)
    assert len(set(grids)) == 1, f"grid point count drifted across the stencil: {grids}"
    return fd4(samples, STEP_DM)


def fd_position(mol, positions, dm, index, *, isodensity, components="cpcm"):
    """dE/dR along one cartesian coordinate at fixed density matrix"""
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
    """Unit-norm symmetric density-matrix perturbations"""
    rng = np.random.default_rng(seed)
    for _ in range(count):
        direction = rng.standard_normal((nao, nao))
        direction = 0.5 * (direction + direction.T)
        yield direction / np.linalg.norm(direction)


def sampled_coordinates(natm):
    """Cartesian coordinates to difference: all for a small solute"""
    total = 3 * natm
    if total <= 9:
        return list(range(total))
    return [int(i) for i in np.linspace(0, total - 1, 3, dtype=int)]


# ----------------------------------------------------------------------
# L0 -- solute-vdW cavity: electrostatic conventions only
# ----------------------------------------------------------------------


@pytest.mark.vdw
@pytest.mark.parametrize(
    "cavity_type",
    [CavityDROPSvdW, CavityDROPCFC, CavityISwiG],
    ids=("svdw-drop", "cfc-drop", "iswig"),
)
@pytest.mark.parametrize(
    "component_type",
    [ModelComponentCPCM, ModelComponentCOSMO],
    ids=("cpcm", "cosmo"),
)
def test_l0_pcm_components_and_cavity_types_share_the_pyscf_coupling(
    cavity_type,
    component_type,
):
    """Every PCM/cavity combination uses the same PySCF coupling"""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    host = PySCFHost(mol)
    model = GeneralSolvationModel(cavity_type(nleb=26), [component_type(EPSILON)])

    result = model.evaluate(coupling=host.coupling(dm))

    assert np.isfinite(result.energy)
    assert result.fock.shape == dm.shape
    assert result.gradient.shape == (3, mol.natm)
    assert result.cavity.ngrid > 0


@pytest.mark.vdw
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l0_fock_matches_fd(system, basis):
    """dE/dP through the surface potential alone, with a fixed surface"""
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
    """dE/dR at fixed P: moist's geometry terms plus the host's AO derivative"""
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
    """Without ``qefield`` the gradient is refused, not silently wrong

    An unsupplied channel used to read as zero: only the nuclear half of the
    surface motion was kept, giving a plausible wrong gradient. The
    external-potential gradient path now requires the channel
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=False)

    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)  # no qefield
    with pytest.raises(RuntimeError, match="electrostatics%qefield"):
        model.get_gradient(mol.natm)


# ----------------------------------------------------------------------
# L1 -- isodensity cavity at fixed density matrix
# ----------------------------------------------------------------------


@pytest.mark.host
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_density_callback_derivatives(system, basis):
    """Callback derivative orders are mutually consistent

    Validates the Leibniz expansion and the PySCF derivative-component ordering
    independently of moist, so a failure is not the cavity's
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
    """Suppressed BLAS status flags in the callback are false positives

    The products raise spurious divide-by-zero and overflow, so the flags are
    ignored; values asserted finite and equal to a BLAS-free reference
    - at a normal point and far in the tail, where the density underflows
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
        # Not bitwise: BLAS and the einsum loop sum in different orders
        # Equal to rounding suffices: the flags carry no information
        assert np.isfinite(product).all()
        np.testing.assert_allclose(product, reference, rtol=5e-12, atol=1e-16)


@pytest.mark.host
def test_grid_integrals_are_contracted_in_place():
    """Grid-integral contractions read PySCF's buffers in place

    - at a few thousand surface points the tensors run to gigabytes, copying
      first costs more than the integrals
    - the copy-free view depends on PySCF's memory layout: pinned, with the pair
      index the contractions assume
    """
    from .pyscf import _pair_view

    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    nao = mol.nao
    coords = mol.atom_coords() + np.array([1.1, -0.7, 0.9])
    for name, lead in (("int1e_grids", ()), ("int1e_grids_ip", (3,))):
        ints = mol.intor(name, grids=coords)
        view = _pair_view(ints)
        assert np.shares_memory(view, ints), f"{name} is copied"
        assert view.shape == lead + (nao * nao, len(coords))
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            contracted = dm.T.ravel() @ view
        assert np.isfinite(contracted).all()
        np.testing.assert_allclose(
            contracted,
            np.einsum("...iuv,uv->...i", ints, dm, optimize=False),
            rtol=1e-13,
            atol=1e-15,
        )


@pytest.mark.host
@pytest.mark.parametrize("cart", [False, True], ids=["spherical", "cartesian"])
def test_internal_cavity_density_is_the_host_density(cart):
    """Density the host installs on an internal cavity is its own

    moist's internal evaluator sums the installed matrix over bare cartesian
    monomials on the contracted radials it was given
    - rebuilt here in numpy, compared with PySCF's density at random points
    - cc-pVTZ brings f functions and general contractions
    - cartesian basis takes another transform than spherical
    - internal vs callback cavities compared end to end in
      ``test_hessian_directional.py``
    """
    mol = gto.M(atom=SYSTEMS["water"].atom, basis="cc-pvtz", cart=cart, unit="Angstrom", verbose=0)
    dm = scf.RHF(mol).get_init_guess()
    host = make_host(mol, dm=dm)
    cavity = host.internal_cavity(nleb=NLEB)
    layout, basis = cavity.layout, host._internal_basis.arrays
    assert max(basis["shell_l"]) == 3

    points = np.random.default_rng(5).normal(scale=1.5, size=(64, 3))
    monomials = np.empty((len(points), layout.ncart))
    primitive = np.concatenate([[0], np.cumsum(basis["shell_nprim"])])
    for s, atom in enumerate(basis["shell_atom"]):
        offset = points - mol.atom_coord(atom)
        exps = basis["exps"][primitive[s]:primitive[s + 1]]
        coeffs = basis["coeffs"][primitive[s]:primitive[s + 1]]
        radial = np.exp(-np.square(offset).sum(axis=1)[:, None] * exps) @ coeffs
        for c in range(layout.shell_offset[s], layout.shell_offset[s + 1]):
            monomials[:, c] = np.prod(offset ** layout.powers[c], axis=1) * radial

    dcart = host._internal_basis.density(dm)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        internal = np.einsum("ic,cd,id->i", monomials, dcart, monomials)
    reference = [host.density(point, 1)[0] for point in points]
    np.testing.assert_allclose(internal, reference, rtol=1e-12, atol=1e-16)


@pytest.mark.conventions
def test_internal_cavity_is_fed_only_by_its_host():
    """Internal cavity not built by the host is refused, not silently stale"""
    mol, dm = molecule(*PRIMARY_CASE), reference_density(*PRIMARY_CASE)
    host = make_host(mol, dm=dm)

    fed = host.internal_cavity(nleb=NLEB, tolerance=PROJ_TOL)
    model = GeneralSolvationModel(fed, COMPONENTS["cpcm"]())
    reference = GeneralSolvationModel(
        CavityDROPIsodensityCallback(host, nleb=NLEB, tolerance=PROJ_TOL), COMPONENTS["cpcm"]()
    )
    energy = model.evaluate(coupling=host.coupling(dm)).energy
    assert energy == pytest.approx(reference.evaluate(coupling=host.coupling(dm)).energy, rel=1e-12)

    # Same basis, built by hand: the host never installs a density
    stray = CavityDROPIsodensityInternal(
        **host._internal_basis.arrays, rho_iso=host.rho_iso, scale=host.scale, nleb=NLEB
    )
    stray.set_density(host._internal_basis.density(dm))
    unfed = GeneralSolvationModel(stray, COMPONENTS["cpcm"]())
    with pytest.raises(ValueError, match="host.internal_cavity"):
        unfed.evaluate(coupling=host.coupling(dm))

    from .hessian import rhf_hessian

    with pytest.raises(ValueError, match="host.internal_cavity"):
        rhf_hessian(scf.RHF(mol), unfed, host)


@pytest.mark.conventions
def test_solvated_rhf_selects_the_isodensity_backend():
    mol = molecule(*PRIMARY_CASE)
    with pytest.raises(ValueError, match="isodensity backend"):
        solvated_rhf(mol, EPSILON, isodensity="grid")
    with pytest.raises(TypeError, match="model_factory owns"):
        solvated_rhf(mol, model_factory=lambda host: None, isodensity="internal")


@pytest.mark.isodensity
@pytest.mark.parametrize("components", COMPONENT_PARAMS)
@pytest.mark.parametrize("system,basis", CASE_PARAMS)
def test_l1_fock_matches_fd(system, basis, components):
    """dE/dP with the surface following the density, the headline Fock test"""
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
    """dE/dR at fixed P, including the level set's own basis-center derivative"""
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
    """With no electrostatic component the whole Fock is the ``lsf`` contraction

    A pv-only model has no electrostatic channel, so only the cavity-shape
    response survives: the sharpest isolation of the weights
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, potential, coords, _ = solve(host, isodensity=True, components="pv")

    assert potential.electrostatics is None
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
    """Dropping the level-set response breaks the Fock test outright"""
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
    """Isodensity level set reports zero nuclear partials by construction

    moist's gradient therefore lacks the density's own dependence on the nuclei;
    the host must add it
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
    """``w_xyz`` carries the dominant part of the cavity response

    - a density change moves the grid points, and ``phi(r_i)`` with them
    - moist cannot see that route (``phi`` is the host's function), so it must
      arrive as ``w_xyz`` before the potential is read
    - omitted, ``lsf`` weights look healthy and are badly wrong
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=True)

    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)  # no w_xyz
    starved = host.fock(coords, model.response(), include_lsf=True)

    direction = next(iter(symmetric_directions(dm.shape[0], 1)))
    numerical = fd_density(mol, dm, direction, isodensity=True)
    analytic = float(np.einsum("uv,uv->", starved, direction))
    assert deviation(analytic, numerical) > VACUITY_FACTOR


@pytest.mark.conventions
def test_gradient_path_ignores_host_surface_weights():
    """Gradient path ignores host surface weights, unlike the potential path

    ``w_xyz`` is read when the potential is assembled, dropped for the gradient:
    scaling it by a thousand leaves the gradient unchanged
    - makes it safe for :meth:`PySCFHost.solve` to supply ``w_xyz`` and
      ``qefield`` together
    - were the gradient path to read ``w_xyz``, the surface-motion term would
      count twice and ``test_l1_gradient_matches_fd`` would fail
    - bound tight, not exact
    """
    system, basis = PRIMARY_CASE
    mol, dm = molecule(system, basis), reference_density(system, basis)
    host = make_host(mol, dm=dm)
    _, _, coords, model = solve(host, isodensity=True)
    phi = host.surface_potential(coords)
    model.supply_electrostatics(phi)
    charges = model.trace_response().electrostatics.surface_charge

    gradients = []
    for factor in (0.0, 1000.0):
        model.supply_electrostatics(
            phi,
            w_xyz=factor * host.surface_position_weights(coords, charges),
            qefield=host.qefield(coords, charges),
        )
        gradients.append(model.get_gradient(mol.natm))
    np.testing.assert_allclose(gradients[0], gradients[1], rtol=1.0e-11, atol=1.0e-14)


# ----------------------------------------------------------------------
# L2 -- self-consistent solvated SCF (slow tier)
# ----------------------------------------------------------------------


@pytest.mark.scf
def test_l2_scf_converges_and_stabilises():
    """Solvated SCF converges, is stabilising, and reduces to gas phase"""
    mol = molecule(*PRIMARY_CASE)
    gas = scf.RHF(mol).run()
    solvated = solvated_rhf(mol, EPSILON, nleb=NLEB, tolerance=PROJ_TOL)
    assert solvated.converged
    assert solvated.e_tot < gas.e_tot

    vacuum = solvated_rhf(mol, 1.0, nleb=NLEB, tolerance=PROJ_TOL)
    assert vacuum.e_tot == pytest.approx(gas.e_tot, abs=1e-9)


@pytest.mark.scf
def test_l2_total_gradient_matches_fd():
    """Total solvated SCF gradient against FD of the converged energy

    Density response drops out at convergence: the analytic gradient is the
    ordinary RHF gradient of the solvated orbitals plus the explicit solvation
    terms at the converged density
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

    # Solvated gradient must differ from gas phase, else the solvation terms
    # above are not exercised
    gas_gradient = grad.RHF(scf.RHF(mol).run()).kernel().T
    assert np.abs(analytic - gas_gradient).max() > MIN_SIGNAL
