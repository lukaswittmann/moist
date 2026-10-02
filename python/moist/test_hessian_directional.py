"""Directional second-order host protocol and the SCF Hessian built on it

Companion of ``test_hessian.py``, which pins the dense reference backend
(solvent second derivatives in the full host parameter space). This module
pins the *protocol* shared by the Fortran, C and Python layers:
- moist differentiates along a batch of directions
- the host returns tangents of the host data entering surface map and
  electrostatics
- moist returns a response tangent per direction; the host completes the
  Hessian columns and Fock tangents from it

A direction has a nuclear part moist sees and a host-private part (here a
density direction) it never sees, so one call covers the RR, RP, PR, PP blocks.
Pinned, in rising order of assembly:

  * ``protocol_matches_fd``: completed nuclear columns and Fock tangents vs
    central differences of the *analytic* gradient and Fock matrix
    - mixed nuclear/density direction, two steps, rebuilt models
    - whole chain in one assertion; the only test that notices a missing term
  * ``branched_cavity_matches_fd``: same assertion on a cavity with
    multi-branch anchor groups (softmax branch weights couple grid points)
    - plus the dense backend's second derivatives contracted with the direction
    - only test of the branch-weight terms on an isodensity cavity, either backend
  * ``internal_backend_matches_callback``: internal isodensity cavity (density
    from basis and transformed density matrix) vs the callback one
    - surface, first derivatives, Hessian columns, to rounding
  * ``response_tangent_channels_match_fd``: returned tangents (surface points,
    surface charges, level-set weights) vs differences of the same quantities
    - hosts read these directly, so pinned independently of the completion
  * ``matches_dense``: total SCF Hessian, directional vs dense backend
    (``test_hessian.py`` validates dense against FD of SCF gradients)
    - no shared code below ``Evaluation``: machine-precision agreement checks both
  * ``scf_hessian_matches_fd``: three-atom systems vs central differences of
    reconverged analytic SCF gradients, anchoring the chain outside moist
  * ``scf_hvp_matches_fd``: same anchor for fluoroacetate, contracted with one
    direction instead of 6N columns
    - two reconverged SCF runs make a seven-atom anion affordable
  * the guards: host callback that raises returns its own traceback; unknown
    method and superseded evaluation are refused

Systems are :data:`SYSTEMS`:
- water at STO-3G everywhere, plus def2-SVP (only d functions in the suite)
- fluoroacetate anion (only polyatomic, charged, strongly asymmetric cavity)
- water with a neon atom (only cavity with multi-branch anchor groups)

Tolerances of the two reference kinds, six orders apart
-------------------------------------------------------
References are central differences, four-point except the two-point ones of
:data:`FD_TOL`; four-point truncation falls as ``h**4``. What differs is the
noise floor of the differenced quantity.

``protocol_matches_fd``: analytic gradient and Fock matrix, no SCF in the loop
- noise is round-off only: 1.4e-12 over every case at ``h = 1e-3``
- bounded at :data:`PROTOCOL_TOL`
- says the moist half (cavity plus model) is right; tighten this one first

The SCF references cannot follow, for reasons outside moist
- PySCF's own gas-phase analytic Hessian (no cavity, no moist) vs four-point
  difference of its own analytic gradient: 7.1e-9 water/STO-3G, 5.8e-8
  water/def2-SVP
- its analytic Hessian is asymmetric by 1.6e-8 (water/def2-SVP), 6.0e-7
  (gas-phase fluoroacetate), 9.7e-16 (water/STO-3G); a difference of a true
  gradient is symmetric to its noise, so the asymmetry is on the analytic side
- both insensitive to:
  - ``conv_tol_grad`` (floors at 1e-10 for the solvated SCF; 1e-11 does not converge)
  - ``conv_tol_cpscf`` over 1e-8..1e-14 (Hessian bit-identical)
  - ``direct_scf_tol`` over 1e-13..1e-20
  - cavity projection tolerance
  - ``h`` over a decade
- ``conv_tol`` cannot tighten: 1e-14 is below one ulp of a -76 Ha energy, SCF
  stops converging

The solvent moves none of it
- solvated water/STO-3G 7.3e-9 vs gas phase 7.1e-9
- solvated fluoroacetate *less* asymmetric (2.6e-7) than gas phase (6.0e-7)
- so :data:`SCF_TOL` is set from the measured floor per system
- below ~1e-8 would mean validating PySCF's Hessian first, not moist's job
- still the only anchor entirely outside moist
"""

import os

import numpy as np
import pytest

try:
    from pyscf import gto, lib, scf
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from .hessian import (
    DirectionalHost,
    finite_difference_gradient_tangent,
    finite_difference_hessian,
    solvated_rhf_gradient,
    solvent_hvp,
)
from .interface import (
    CavityDROPIsodensityCallback,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentPV,
    SolvationModel,
)
from .pyscf import PySCFHost, solvated_rhf


pytestmark = pytest.mark.scf

#: Cavity settings shared by every model
#: - small enough to rebuild many times
#: - tight enough that the projection does not move under an FD step
CAVITY = dict(nleb=50, tolerance=1.0e-13)

#: Nondefault level-set scale, so a dropped or squared ``scale`` shows
SCALE = 2000.0

#: Two-point central-difference steps
#: - two, so a value agreeing at one step only still fails
STEPS = (1.0e-4, 5.0e-5)

#: Agreement bound of the two-point references
#: - worst measured over every configuration and both steps: 5.5e-10, values ~0.3
FD_TOL = 1.0e-8

#: Central four-point rule, offsets in steps with weights over a common denominator
#: - the rule of ``finite_difference_hessian(stencil=4)``
#: - written out: several quantities of different shapes differenced in one pass
FOUR_POINT = ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0))
FOUR_POINT_DENOMINATOR = 12.0

#: Steps of the four-point references
#: - truncation falls as ``h**4``; binding wall is round-off, *growing* as h shrinks
#: - worst over every protocol case: 1.4e-12, 6.1e-12, 1.5e-11, 7.1e-11 at
#:   h = 1e-3, 3e-4, 1e-4, 3e-5
#: - hence steps far above the two-point ones; smaller is worse
PROTOCOL_STEPS = (1.0e-3, 3.0e-4)

#: Agreement bound of the four-point references
#: - ~8x the worst measured deviation (6.1e-12, at the smaller step)
PROTOCOL_TOL = 5.0e-11

#: Step of the four-point SCF references
#: - truncation ``h**4``, round-off ``eps max|g| / h``: both far below the floor hit
#: - so a *large* step is right
SCF_STEP = 1.0e-2

#: Per-system agreement bounds of the SCF FD references
#: - each ~5x the worst deviation measured here
#: - set by PySCF's own gradient/Hessian consistency, not moist (see module
#:   docstring); measured gas-phase floor in the line comment
SCF_TOL = {
    "water/sto-3g": 5.0e-8,          # measured 7.3e-9; gas phase alone 7.1e-9
    "water/def2-svp": 7.0e-7,        # measured 1.4e-7; gas phase alone 5.8e-8
    "fluoroacetate/sto-3g": 2.5e-7,  # measured 4.3e-8
}

#: Reference below this carries no information
VACUITY = 1.0e-4

COMPONENT_SETS = {
    "cpcm": [(ModelComponentCPCM, 80.0)],
    "cosmo": [(ModelComponentCOSMO, 12.0)],
    "cpcm+cosmo": [(ModelComponentCPCM, 4.0), (ModelComponentCOSMO, 12.0)],
    "pv": [(ModelComponentPV, 1.0e-3)],
}


#: Asymmetric geometry, exercises off-diagonal atom and Cartesian blocks
WATER = "O 0 0 -0.3893; H 0.7629 0 0.1947; H -0.7991 0.0953 0.2223"

#: Fluoroacetate, FCH2-COO(-), asymmetric non-stationary geometry
#: - no test differentiates at a minimum
#: - symmetric cavity would risk multi-branch projections, which
#:   :data:`WATER_NEON` pins on purpose
FLUOROACETATE = (
    "C 0.000 0.000 0.000; F 1.180 0.560 0.420; "
    "H -0.420 -0.900 0.470; H -0.120 0.180 -1.060; "
    "C -0.870 1.060 0.680; O -1.980 1.300 0.130; O -0.480 1.600 1.760"
)

#: Water with a neon atom 3.2 A below the oxygen
#: - concave isodensity seam between the two; four points in multi-branch anchor groups
#: - one pair with softmax weights ~0.42/0.58, where weights move fastest
#: - one pair at 0/1, where the weight guards act
#: - branching depends on distance (3.0 A, 3.4 A give other patterns), so the
#:   test asserts the fixture still branches
WATER_NEON = WATER + "; Ne 0.1 0.2 -3.2"

#: Systems the isodensity derivatives are pinned on
#: - water/STO-3G: fast default, every component set
#: - water/def2-SVP: only d functions, so only test of the AO derivative table
#:   beyond s and p
#: - fluoroacetate: only polyatomic, charged, strongly asymmetric case
#:   - anion SCF needs a raised iteration cap to reach these tolerances
#: - water-neon: only multi-branch cavity
SYSTEMS = {
    "water/sto-3g": dict(atom=WATER, basis="sto-3g", charge=0, max_cycle=None),
    "water/def2-svp": dict(atom=WATER, basis="def2-svp", charge=0, max_cycle=None),
    "fluoroacetate/sto-3g": dict(atom=FLUOROACETATE, basis="sto-3g", charge=-1, max_cycle=300),
    "water-neon/sto-3g": dict(atom=WATER_NEON, basis="sto-3g", charge=0, max_cycle=None),
}


def molecule(name):
    """Build a PySCF molecule and its SCF cycle cap from a :data:`SYSTEMS` entry"""
    spec = SYSTEMS[name]
    mol = gto.M(
        atom=spec["atom"],
        basis=spec["basis"],
        charge=spec["charge"],
        unit="Angstrom",
        verbose=0,
    )
    return mol, spec["max_cycle"]


@pytest.fixture
def water():
    return molecule("water/sto-3g")[0]


def converged_rhf(mol, max_cycle=None):
    mean_field = scf.RHF(mol)
    mean_field.conv_tol = 1.0e-13
    mean_field.conv_tol_grad = 1.0e-10
    if max_cycle is not None:
        mean_field.max_cycle = max_cycle
    mean_field.kernel()
    assert mean_field.converged, "SCF must converge before differentiating"
    return mean_field


def build(mol, specs, *, internal=False):
    """Build a host and model of one component set on the isodensity cavity

    ``internal`` selects the cavity evaluating the density from the basis, not
    through the host callback
    """
    host = PySCFHost(mol, scale=SCALE)
    components = [component(value) for component, value in specs]
    if internal:
        cavity = host.internal_cavity(**CAVITY)
    else:
        cavity = CavityDROPIsodensityCallback(host, **CAVITY)
    return host, SolvationModel(cavity, components)


def mixed_direction(mol, dm, seed=11):
    """Build one normalised direction with a nuclear and a density part"""
    rng = np.random.default_rng(seed)
    nuclear = rng.normal(size=(3, mol.natm))
    nuclear /= np.linalg.norm(nuclear)
    density = rng.normal(size=dm.shape)
    density = 0.5 * (density + density.T)
    density /= np.linalg.norm(density)
    return np.asfortranarray(nuclear), density


def displaced(mol, dm, nuclear, density, step):
    """Move molecule and density along a mixed direction by ``step``"""
    moved = mol.set_geom_(mol.atom_coords() + step * nuclear.T, unit="Bohr", inplace=False)
    return moved, dm + step * density


def worst(value, reference):
    """Worst absolute deviation and its scale"""
    value = np.asarray(value)
    reference = np.asarray(reference)
    return np.abs(value - reference).max(), np.abs(reference).max()


#: `(system, component set)` pairs of the protocol test
#: - water/STO-3G runs every component set
#: - the two costlier systems run one electrostatic, one non-electrostatic set
#: - together they touch every exchange channel without the full product
PROTOCOL_CASES = [("water/sto-3g", name) for name in COMPONENT_SETS] + [
    (system, name)
    for system in ("water/def2-svp", "fluoroacetate/sto-3g")
    for name in ("cpcm", "pv")
]


@pytest.mark.parametrize("system,name", PROTOCOL_CASES)
def test_protocol_matches_fd(system, name):
    """Whole exchange matches differences of the analytic first derivatives

    Sharpest check of the *moist* half (cavity plus model)
    - analytic gradient and Fock matrix differenced along a mixed
      nuclear/density direction, no SCF in the loop
    - a wrong cavity or component term shows undiluted by electronic relaxation
    """
    specs = COMPONENT_SETS[name]
    with lib.with_omp_threads(1):
        mol, max_cycle = molecule(system)
        dm = converged_rhf(mol, max_cycle).make_rdm1()
        nuclear, density = mixed_direction(mol, dm)

        host, model = build(mol, specs)
        evaluation = model.evaluate(coupling=host.coupling(dm))
        columns, fock = solvent_hvp(evaluation, host, nuclear[:, :, None], density[None])

        def first_derivatives(offset, step):
            moved, moved_dm = displaced(mol, dm, nuclear, density, offset * step)
            moved_host, moved_model = build(moved, specs)
            result = moved_model.evaluate(coupling=moved_host.coupling(moved_dm))
            points = moved_model.cavity.snapshot().xyz
            return (np.asarray(result.gradient).copy(),
                    np.asarray(result.fock).copy(),
                    points.shape)

        for step in PROTOCOL_STEPS:
            sampled = {o: first_derivatives(o, step) for o, _ in FOUR_POINT}
            shapes = {value[2] for value in sampled.values()}
            assert len(shapes) == 1, "the projected grid changed under the step"

            def combine(channel):
                total = sum(weight * sampled[o][channel] for o, weight in FOUR_POINT)
                return total / (FOUR_POINT_DENOMINATOR * step)

            reference = combine(0)
            deviation, scale = worst(columns[:, :, 0], reference)
            assert scale > VACUITY, "the differenced gradient carries no information"
            assert deviation < PROTOCOL_TOL, (
                f"nuclear columns miss the differenced gradient by {deviation:.3e} "
                f"against {scale:.3e} at step {step:g}"
            )

            reference = combine(1)
            deviation, scale = worst(fock[0], reference)
            assert scale > VACUITY, "the differenced Fock matrix carries no information"
            assert deviation < PROTOCOL_TOL, (
                f"Fock tangents miss the differenced Fock matrix by {deviation:.3e} "
                f"against {scale:.3e} at step {step:g}"
            )


#: Dense vs directional agreement on the branched cavity
#: - measured 2.9e-16 (nuclear), 2.2e-16 (Fock) on values ~0.05, ~0.2
#: - without the dense branch-weight pass: 2.1e-3, 4.1e-4
#: - nuclear deviation with each term of that pass removed:
#:   - owner shift of the branch objective 4.0e-3
#:   - first-order weight motion 1.1e-3
#:   - quadratic weight term 1.1e-3
#:   - second derivative of the weights 2.0e-4
#:   - cross term with the frozen width derivative 2.3e-5
#:   - quadratic term of the objective's second derivative 1.4e-5
#: - measured 2026-10-01
BRANCHED_DENSE_TOL = 1.0e-12

#: Steps and agreement bound of the four-point reference on the branched cavity
#: - softmax enlarges high derivatives, so truncation shows at ``h = 1e-3``
#:   (PV: 8.7e-11, sixteenfold less per halving); steps below :data:`PROTOCOL_STEPS`
#: - worst measured at these steps: 4.6e-11 (PV), 6.4e-13 (CPCM), values ~0.03, ~0.05
#: - bound ~6x that
BRANCHED_STEPS = (5.0e-4, 3.0e-4)
BRANCHED_FD_TOL = 3.0e-10


@pytest.mark.parametrize("name", ["cpcm", "pv"])
def test_branched_cavity_matches_fd(name):
    """Both backends carry the branch weights of a multi-branch cavity

    :func:`test_protocol_matches_fd` on :data:`WATER_NEON`, whose cavity has
    multi-branch anchor groups
    - softmax branch weights scale the Lebedev weights, hence the Gaussian
      widths; their derivatives couple the points of a group
    - dense second derivatives contracted with the same direction, held to the
      directional result; PV has no dense counterpart
    """
    specs = COMPONENT_SETS[name]
    with lib.with_omp_threads(1):
        mol, max_cycle = molecule("water-neon/sto-3g")
        dm = converged_rhf(mol, max_cycle).make_rdm1()
        nuclear, density = mixed_direction(mol, dm)

        host, model = build(mol, specs)
        evaluation = model.evaluate(coupling=host.coupling(dm))
        snapshot = model.cavity.snapshot()
        weights = snapshot.wbranch[snapshot.branch_count > 1]
        assert weights.size > 0, "fixture stopped producing branches"
        assert ((weights > 0.1) & (weights < 0.9)).any(), "no branch weight left to move"
        columns, fock = solvent_hvp(evaluation, host, nuclear[:, :, None], density[None])

        def first_derivatives(offset, step):
            moved, moved_dm = displaced(mol, dm, nuclear, density, offset * step)
            moved_host, moved_model = build(moved, specs)
            result = moved_model.evaluate(coupling=moved_host.coupling(moved_dm))
            moved_snapshot = moved_model.cavity.snapshot()
            # Point entering/leaving an anchor group keeps the grid size but makes
            # the stencil discontinuous: groups must stay put
            topology = tuple(
                np.asarray(getattr(moved_snapshot, field)).copy()
                for field in ("anchor_id", "branch", "branch_count")
            )
            return (np.asarray(result.gradient).copy(),
                    np.asarray(result.fock).copy(),
                    topology)

        for step in BRANCHED_STEPS:
            sampled = {o: first_derivatives(o, step) for o, _ in FOUR_POINT}
            topologies = [value[2] for value in sampled.values()]
            for topology in topologies[1:]:
                assert len(topology[0]) == len(topologies[0][0]), "the projected grid changed under the step"
                assert all(np.array_equal(a, b) for a, b in zip(topology, topologies[0])), (
                    "the branch topology changed under the step"
                )
            for channel, value, label in ((0, columns[:, :, 0], "nuclear columns"),
                                          (1, fock[0], "Fock tangents")):
                reference = sum(weight * sampled[o][channel] for o, weight in FOUR_POINT)
                reference /= FOUR_POINT_DENOMINATOR * step
                deviation, scale = worst(value, reference)
                assert scale > VACUITY, f"the differenced {label} carry no information"
                assert deviation < BRANCHED_FD_TOL, (
                    f"{label} miss their finite difference by {deviation:.3e} "
                    f"against {scale:.3e} at step {step:g}"
                )

        if name == "pv":
            return
        dense_host, dense_model = build(mol, specs)
        coupling = dense_host.coupling(dm)
        parameters = coupling.second_order()
        dense = dense_model.evaluate(coupling=coupling).linearize(parameters)
        rows, cols = parameters.rows, parameters.cols
        packed = density[rows, cols]
        dense_columns = dense.nuclear @ nuclear.T.ravel() + dense.mixed @ packed
        dense_fock = dense.mixed.T @ nuclear.T.ravel() + dense.density_response(packed)

    deviation, scale = worst(dense_columns, columns[:, :, 0].T.ravel())
    assert scale > VACUITY
    assert deviation < BRANCHED_DENSE_TOL, (
        f"dense nuclear columns miss the directional ones by {deviation:.3e} against {scale:.3e}"
    )
    packed_fock = fock[0][rows, cols] * np.where(rows == cols, 1, 2)
    deviation, scale = worst(dense_fock, packed_fock)
    assert scale > VACUITY
    assert deviation < BRANCHED_DENSE_TOL, (
        f"dense Fock tangents miss the directional ones by {deviation:.3e} against {scale:.3e}"
    )


#: Internal vs callback isodensity agreement, relative to the largest entry
#: - same density, so agreement to rounding
#: - worst measured 7.7e-14 (fluoroacetate Hessian columns), rest below 7e-15
#: - measured 2026-10-02, also at cc-pVTZ for f functions
BACKEND_TOL = 1.0e-12


@pytest.mark.parametrize("system", ["water/sto-3g", "water/def2-svp", "fluoroacetate/sto-3g"])
def test_internal_backend_matches_callback(system):
    """Internal isodensity cavity reproduces the callback one

    - callback: density through PySCF's AOs point by point
    - internal: density from the basis and the density matrix the host
      transforms into moist's cartesian monomials
    - on one density: same surface, first derivatives and, along a mixed
      direction, Hessian columns and Fock tangents
    - wrong transform, uninstalled or stale density moves all far beyond rounding
    - density is PySCF's initial guess, no SCF
    """
    specs = COMPONENT_SETS["cpcm"] + COMPONENT_SETS["pv"]
    with lib.with_omp_threads(1):
        mol, _ = molecule(system)
        dm = scf.RHF(mol).get_init_guess()
        nuclear, density = mixed_direction(mol, dm)

        def quantities(internal):
            host, model = build(mol, specs, internal=internal)
            evaluation = model.evaluate(coupling=host.coupling(dm))
            snapshot = model.cavity.snapshot()
            columns, fock = solvent_hvp(evaluation, host, nuclear[:, :, None], density[None])
            return {
                "surface points": snapshot.xyz,
                "area": snapshot.area,
                "volume": snapshot.volume,
                "energy": evaluation.energy,
                "Fock matrix": np.asarray(evaluation.fock),
                "gradient": np.asarray(evaluation.gradient),
                "Hessian columns": columns,
                "Fock tangents": fock,
            }

        internal, callback = quantities(True), quantities(False)

    assert internal["surface points"].shape == callback["surface points"].shape
    for label, reference in callback.items():
        deviation, scale = worst(internal[label], reference)
        assert scale > VACUITY, f"the {label} carry no information"
        assert deviation <= BACKEND_TOL * scale, (
            f"internal {label} miss the callback ones by {deviation:.3e} against {scale:.3e}"
        )


def response_channels(mol, dm, specs):
    """Read the base quantities of the response tangent through one call"""
    host, model = build(mol, specs)
    evaluation = model.evaluate(coupling=host.coupling(dm))
    points = model.cavity.snapshot().xyz
    probe = np.zeros((3, mol.natm, 1), order="F")
    probe[0, 0, 0] = 1.0
    tangent_host = DirectionalHost(
        host, dm, points.T, probe, potential=evaluation.response.electrostatics is not None
    )
    _, tangent = evaluation.hvp_coupled(probe, tangent_host, level_set=True)
    channels = {
        "xyz": points,
        "w_value": tangent.w_value,
        "w_gradient": tangent.w_gradient,
        "w_hessian": tangent.w_hessian,
    }
    if evaluation.response.electrostatics is not None:
        channels["surface_charge"] = evaluation.charges
    return channels


@pytest.mark.parametrize("name", ["cpcm", "cpcm+cosmo", "pv"])
def test_response_tangent_channels_match_fd(water, name):
    """Every returned channel is the derivative of what it claims"""
    specs = COMPONENT_SETS[name]
    with lib.with_omp_threads(1):
        dm = converged_rhf(water).make_rdm1()
        nuclear, density = mixed_direction(water, dm, seed=5)

        host, model = build(water, specs)
        evaluation = model.evaluate(coupling=host.coupling(dm))
        points = model.cavity.snapshot().xyz
        tangent_host = DirectionalHost(
            host,
            dm,
            points.T,
            nuclear[:, :, None],
            density[None],
            potential=evaluation.response.electrostatics is not None,
        )
        _, tangent = evaluation.hvp_coupled(nuclear[:, :, None], tangent_host, level_set=True)

        analytic = {
            "xyz": tangent.xyz[:, :, 0],
            "w_value": tangent.dw_value[:, 0],
            "w_gradient": tangent.dw_gradient[:, :, 0],
            "w_hessian": tangent.dw_hessian[:, :, :, 0],
        }
        if evaluation.response.electrostatics is not None:
            analytic["surface_charge"] = tangent.surface_charge[:, 0]

        for step in STEPS:
            plus = response_channels(*displaced(water, dm, nuclear, density, step), specs)
            minus = response_channels(*displaced(water, dm, nuclear, density, -step), specs)
            assert plus["xyz"].shape == minus["xyz"].shape, "the grid changed under the step"
            for channel, value in analytic.items():
                reference = (plus[channel] - minus[channel]) / (2 * step)
                deviation, scale = worst(value, reference)
                assert scale > 0.0, f"the differenced {channel} carries no information"
                assert deviation < FD_TOL * max(1.0, scale), (
                    f"{channel} tangent misses its difference by {deviation:.3e} "
                    f"against {scale:.3e} at step {step:g}"
                )


# PV has no dense counterpart (no second-order channel, pinned in
# ``test_hessian.py``); the directional path covers it in the two FD tests above
@pytest.mark.parametrize("name", ["cpcm", "cosmo", "cpcm+cosmo"])
def test_scf_hessian_matches_dense(water, name):
    """Directional SCF Hessian reproduces the dense one"""
    specs = COMPONENT_SETS[name]

    def model_factory(host):
        return SolvationModel(
            CavityDROPIsodensityCallback(host, **CAVITY),
            [component(value) for component, value in specs],
        )

    with lib.with_omp_threads(1):
        mean_field = solvated_rhf(water, model_factory=model_factory, conv_tol_grad=1.0e-10)
        assert mean_field.converged
        directional = mean_field.Hessian("directional").kernel()
        dense = mean_field.Hessian("dense").kernel()

    assert np.abs(dense).max() > VACUITY
    np.testing.assert_allclose(directional, dense, rtol=0, atol=1.0e-10)
    np.testing.assert_allclose(
        directional, directional.transpose(1, 0, 3, 2), rtol=0, atol=1.0e-10
    )


def solvated_gradient_at(mol, specs, max_cycle):
    """Build a callable giving the reconverged solvated analytic gradient"""

    def model_factory(host):
        return SolvationModel(
            CavityDROPIsodensityCallback(host, **CAVITY),
            [component(value) for component, value in specs],
        )

    def gradient(positions):
        moved = mol.set_geom_(positions, unit="Bohr", inplace=False)
        mean_field = solvated_rhf(
            moved, model_factory=model_factory, conv_tol_grad=1.0e-10, max_cycle=max_cycle
        )
        host = PySCFHost(moved, scale=1000.0)
        return solvated_rhf_gradient(mean_field, model_factory(host), host)

    return model_factory, gradient


@pytest.mark.parametrize("system", ["water/sto-3g", "water/def2-svp"])
def test_scf_hessian_matches_fd(system):
    """Full Hessian matches finite differences of SCF gradients

    Directional path end to end
    - solvated SCF reconverged at every displaced geometry, gradient differenced:
      electronic relaxation in reference and Hessian alike
    - three-atom systems only: 6N reconverged solvated SCF runs;
      fluoroacetate uses :func:`test_scf_hvp_matches_fd`
    """
    mol, max_cycle = molecule(system)
    specs = COMPONENT_SETS["cpcm"]
    model_factory, gradient = solvated_gradient_at(mol, specs, max_cycle)

    with lib.with_omp_threads(1):
        mean_field = solvated_rhf(
            mol, model_factory=model_factory, conv_tol_grad=1.0e-10, max_cycle=max_cycle
        )
        analytic = mean_field.Hessian("directional").kernel()
        numerical = finite_difference_hessian(
            gradient, mol.atom_coords(), step=SCF_STEP, stencil=4
        )

    assert np.isfinite(numerical).all()
    np.testing.assert_allclose(analytic, numerical, rtol=0, atol=SCF_TOL[system])


def test_scf_hvp_matches_fd():
    """Fluoroacetate end to end, one direction instead of the whole block

    :func:`test_scf_hessian_matches_fd` contracted with a single nuclear direction
    - analytic Hessian vs differences of reconverged solvated SCF gradients
      (coupled-perturbed solve in the reference too)
    - two SCF runs instead of 6N: makes the seven-atom charged anion affordable
    - a direction with a component on every atom still touches every column
    """
    mol, max_cycle = molecule("fluoroacetate/sto-3g")
    specs = COMPONENT_SETS["cpcm"]
    model_factory, gradient = solvated_gradient_at(mol, specs, max_cycle)

    rng = np.random.default_rng(7)
    direction = rng.normal(size=(mol.natm, 3))
    direction /= np.linalg.norm(direction)

    with lib.with_omp_threads(1):
        mean_field = solvated_rhf(
            mol, model_factory=model_factory, conv_tol_grad=1.0e-10, max_cycle=max_cycle
        )
        analytic = mean_field.Hessian("directional").kernel()
        numerical = finite_difference_gradient_tangent(
            gradient, mol.atom_coords(), direction, step=SCF_STEP, stencil=4
        )

    # H[A, B, a, b] contracted over the (B, b) half
    reference = np.einsum("ABab,Bb->Aa", analytic, direction)

    assert np.isfinite(numerical).all()
    assert np.abs(numerical).max() > VACUITY, "the differenced gradient carries no information"
    np.testing.assert_allclose(
        reference, numerical, rtol=0, atol=SCF_TOL["fluoroacetate/sto-3g"]
    )


def test_host_callback_failure_propagates(water):
    """Host raising inside the exchange comes back with its own error"""

    class Failing(DirectionalHost):
        def level_set_tangent(self, first, dirs, xyz):
            raise ZeroDivisionError("the host's own failure")

    with lib.with_omp_threads(1):
        dm = converged_rhf(water).make_rdm1()
        host, model = build(water, COMPONENT_SETS["cpcm"])
        evaluation = model.evaluate(coupling=host.coupling(dm))
        points = model.cavity.snapshot().xyz
        probe = np.zeros((3, water.natm, 1), order="F")
        probe[0, 0, 0] = 1.0
        failing = Failing(host, dm, points.T, probe)
        with pytest.raises(ZeroDivisionError, match="the host's own failure"):
            evaluation.hvp_coupled(probe, failing, level_set=True)


def test_reductions_rebuild_when_their_operand_changes(water):
    """Cached reductions are keyed on operand value, not identity

    ``_PotentialIntegrals`` outlives one Hessian-vector call (kept in the
    solve's ``shared`` dict)
    - cached reductions read by the nuclear block and every coupled-perturbed
      iteration; charges and density do not change within a solve
    - a cache answering for a changed operand returns a wrong Hessian silently,
      so the miss is pinned as well as the hit
    """
    from .hessian import _contract, _PotentialIntegrals

    rng = np.random.default_rng(17)
    coords = rng.normal(scale=2.0, size=(23, 3))
    ints = _PotentialIntegrals(water, coords)

    first = rng.normal(size=(water.nao, water.nao))
    first = first + first.T
    second = rng.normal(size=(water.nao, water.nao))
    second = second + second.T
    charges = rng.normal(size=len(coords))
    other = rng.normal(size=len(coords))

    # Hit is the same object: cache is load-bearing, not merely correct
    assert ints.density_contracted(first) is ints.density_contracted(first.copy())
    assert ints.charge_weighted(charges) is ints.charge_weighted(charges.copy())

    # Every miss, incl. return to an operand seen before, rebuilds
    for dm in (first, second, first):
        assert np.array_equal(
            ints.density_contracted(dm).field, _contract("cuvi,uv->ci", ints.ip2, dm)
        )
    for q in (charges, other, charges):
        assert np.array_equal(
            ints.charge_weighted(q).bra, _contract("i,cuvi->cuv", q, ints.ip1)
        )

    # Caller mutating its own array afterwards must not be answered from cache
    held = first.copy()
    ints.density_contracted(held)
    held += 1.0
    assert np.array_equal(
        ints.density_contracted(held).field, _contract("cuvi,uv->ci", ints.ip2, held)
    )


def test_direction_threads_do_not_change_the_result(water):
    """Threaded per-direction loops leave the result unchanged

    Both host callbacks and the completion loop over a pool sized by PySCF's
    thread count
    - directions share read-only operands, write disjoint slices
    - results must be equal *bit for bit*: a tolerance would hide the race this
      catches; the one comparison in the suite made with ``array_equal``

    Only the Python layer switches between the two runs
    - tangent and response taken once under one thread, three loops replayed
    - rerunning the Fortran under a second thread count would test its
      reductions, not these loops
    """
    with lib.with_omp_threads(1):
        dm = converged_rhf(water).make_rdm1()
        host, model = build(water, COMPONENT_SETS["cpcm"])
        evaluation = model.evaluate(coupling=host.coupling(dm))
        points = model.cavity.snapshot().xyz

        ndir = 3 * water.natm
        rng = np.random.default_rng(5)
        dirs = np.asfortranarray(rng.normal(size=(3, water.natm, ndir)))
        dm_dirs = rng.normal(size=(ndir,) + dm.shape)
        dm_dirs = 0.5 * (dm_dirs + dm_dirs.transpose(0, 2, 1))

        exchange = DirectionalHost(host, dm, points.T, dirs, dm_dirs)
        hvp, tangent = evaluation.hvp_coupled(dirs, exchange, level_set=True)
        q = evaluation.charges

    def under(nthread):
        with lib.with_omp_threads(nthread):
            assert lib.num_threads() == nthread, "the thread count did not take"
            jets = exchange.level_set_tangent(0, dirs, points)
            efield, dphi, defield = exchange.field_tangent(0, dirs, points, tangent.xyz)
            nuclear, fock = exchange.complete(q, hvp, tangent)
        return jets, efield, dphi, defield, nuclear, fock

    assert ndir > 1, "one direction would take the serial path in both runs"
    for serial, threaded in zip(under(1), under(4)):
        assert np.array_equal(serial, threaded)


def test_guards(water):
    """Missing half of the exchange, stale evaluation and bad method are refused"""
    with lib.with_omp_threads(1):
        dm = converged_rhf(water).make_rdm1()
        probe = np.zeros((3, water.natm, 1), order="F")
        probe[0, 0, 0] = 1.0

        # No host object: level set is the host's, all its tangents would be
        # an unphysical zero
        host, model = build(water, COMPONENT_SETS["cpcm"])
        evaluation = model.evaluate(coupling=host.coupling(dm))
        with pytest.raises(RuntimeError, match="level set is the host's"):
            evaluation.hvp_coupled(probe)

        # Host answering for the level set but not its own potential, which an
        # electrostatic component on a host potential needs
        points = model.cavity.snapshot().xyz
        complete = DirectionalHost(host, dm, points.T, probe)

        class LevelSetOnly:
            def level_set_tangent(self, first, dirs, xyz):
                return complete.level_set_tangent(first, dirs, xyz)

        with pytest.raises(RuntimeError, match="field tangent callback"):
            evaluation.hvp_coupled(probe, LevelSetOnly(), level_set=True)

        # Model reading no host potential needs no field callback, produces no
        # charges to differentiate
        pv_host, pv_model = build(water, COMPONENT_SETS["pv"])
        pv = pv_model.evaluate(coupling=pv_host.coupling(dm))
        pv_points = pv_model.cavity.snapshot().xyz
        pv_tangent_host = DirectionalHost(pv_host, dm, pv_points.T, probe, potential=False)
        columns, tangent = pv.hvp_coupled(probe, pv_tangent_host, level_set=False)
        assert columns.shape == probe.shape
        assert tangent.w_value is None, "level-set channels appeared unasked"
        assert np.abs(tangent.xyz).max() > 0.0
        assert np.abs(tangent.surface_charge).max() == 0.0

        # On by default for a density-defined cavity, the only kind whose
        # level-set weights a host uses
        _, tangent = pv.hvp_coupled(probe, pv_tangent_host)
        assert tangent.w_value is not None

        model.evaluate(coupling=host.coupling(dm))
        with pytest.raises(RuntimeError, match="superseded"):
            evaluation.hvp_coupled(probe, complete, level_set=True)

    mean_field = converged_rhf(water)
    with pytest.raises(ValueError, match="unknown Hessian method"):
        from .hessian import rhf_hessian

        rhf_hessian(mean_field, model, host, method="bogus")
