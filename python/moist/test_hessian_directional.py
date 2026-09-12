"""The directional second-order host protocol, and the SCF Hessian built on it.

The companion of ``test_hessian.py``. That module pins the dense reference
backend, which assembles the solvent second derivatives in the full host
parameter space. This one pins the *protocol* the Fortran, C and Python layers
share: moist differentiates along a batch of directions, asks the host for the
tangents of whatever host data enter the surface map and the electrostatics,
and hands back a response tangent per direction; the host completes the
Hessian columns and the Fock-matrix tangents from it.

A direction has a nuclear part moist sees and a host-private part -- here a
density direction -- that it never does, so one call covers the RR, RP, PR and
PP blocks alike. What is pinned, in rising order of assembly:

  * ``protocol_matches_fd``: the completed nuclear columns and Fock tangents
    against central differences of the *analytic* gradient and Fock matrix
    along a mixed nuclear/density direction, at two steps, over rebuilt
    models. This is the whole chain in one assertion and the only test that
    would notice a term missing from either side of the exchange.
  * ``response_tangent_channels_match_fd``: the tangents moist returns --
    surface points, surface charges, level-set weights -- against differences
    of the same quantities. A host reads these directly, so they are pinned
    independently of what the completion does with them.
  * ``matches_dense``: the total SCF Hessian through the directional path
    against the dense backend, which ``test_hessian.py`` validates against
    finite differences of SCF gradients. The two share no code below
    ``Evaluation``, so agreement to machine precision is a strong check of
    both, and a cheap one.
  * ``scf_hessian_matches_fd``: the three-atom systems carried all the way to
    central differences of reconverged analytic SCF gradients, so the chain is
    anchored to something outside moist as well.
  * ``scf_hvp_matches_fd``: the same anchor for fluoroacetate, contracted with
    one direction rather than assembled column by column -- two reconverged
    SCF runs instead of 6N, which is what makes a seven-atom anion affordable.
  * the guards: a host callback that raises comes back with its own traceback,
    an unknown method and a superseded evaluation are refused.

The systems are :data:`SYSTEMS`: water at STO-3G everywhere, plus def2-SVP for
the only d functions in the suite and the fluoroacetate anion for the only
polyatomic, charged, strongly asymmetric cavity.

Why the two kinds of reference have tolerances six orders apart
---------------------------------------------------------------
Every reference here is a four-point central difference, whose truncation falls
as ``h**4``; what differs is the noise floor of the quantity being differenced.

``protocol_matches_fd`` differences an *analytic* gradient and Fock matrix with
no SCF in the loop, so the only noise is round-off: it lands at 1.4e-12 over
every case at ``h = 1e-3`` and is bounded at :data:`PROTOCOL_TOL`. That is the
test that says the moist half -- cavity plus model -- is right, and it is the
one to tighten further if any of this ever needs to be sharper.

The SCF references cannot follow it there, for reasons outside moist. Measured
on this machine, PySCF's *own* gas-phase analytic Hessian -- no cavity, no
moist -- differs from a four-point difference of its own analytic gradient by
7.1e-9 for water/STO-3G and 5.8e-8 for water/def2-SVP, and its analytic Hessian
is asymmetric by 1.6e-8 (water/def2-SVP) and 6.0e-7 (gas-phase fluoroacetate),
against 9.7e-16 for water/STO-3G. A finite difference of a true gradient is
symmetric to its own noise, so an asymmetry that large lives on the analytic
side. Both numbers are insensitive to ``conv_tol_grad`` (which floors at 1e-10
for the solvated SCF, and at 1e-11 does not converge at all), to
``conv_tol_cpscf`` over 1e-8..1e-14 (the Hessian comes back bit-identical), to
``direct_scf_tol`` over 1e-13..1e-20, to the cavity projection tolerance, and
to ``h`` over a decade. ``conv_tol`` cannot be tightened either: at 1e-14 it is
below one ulp of a -76 Ha energy and the SCF stops converging.

Adding the solvent moves none of it: solvated water/STO-3G lands at 7.3e-9
against the gas-phase 7.1e-9, and solvated fluoroacetate is *less* asymmetric
(2.6e-7) than the same molecule in gas phase (6.0e-7). So :data:`SCF_TOL` is
set from the measured floor per system, and pushing these particular tests
below ~1e-8 would mean validating PySCF's Hessian first, which is not moist's
job. They remain worth having as the only anchor outside moist entirely.
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
    CavityDROPIsodensity,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentPV,
    SolvationModel,
)
from .pyscf import PySCFHost, solvated_rhf


pytestmark = pytest.mark.scf

#: Cavity settings shared by every model here: small enough to rebuild many
#: times, tight enough that the projection does not move under an FD step.
CAVITY = dict(nleb=50, tolerance=1.0e-13)

#: A nondefault level-set scale, so that a dropped or squared ``scale`` shows.
SCALE = 2000.0

#: Central-difference steps; two, so a value that agrees at one step only
#: still fails. Used by the two-point references.
STEPS = (1.0e-4, 5.0e-5)

#: Agreement bound of the two-point references. Worst measured over every
#: configuration and both steps is 5.5e-10 on quantities of order 0.3.
FD_TOL = 1.0e-8

#: The central four-point rule, as offsets in units of the step with their
#: weights over the common denominator -- the rule
#: ``finite_difference_hessian(stencil=4)`` applies, written out because these
#: tests difference several quantities of different shapes in one pass.
FOUR_POINT = ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0))
FOUR_POINT_DENOMINATOR = 12.0

#: Steps for the four-point references. Truncation falls as ``h**4`` here, so
#: the binding wall is round-off, which *grows* as the step shrinks: measured
#: over every protocol case the worst deviation runs 1.4e-12, 6.1e-12, 1.5e-11,
#: 7.1e-11 at h = 1e-3, 3e-4, 1e-4, 3e-5. Hence steps far larger than the
#: two-point ones -- going smaller makes this reference worse, not better.
PROTOCOL_STEPS = (1.0e-3, 3.0e-4)

#: Agreement bound of the four-point references, about eight times the worst
#: measured deviation (6.1e-12, at the smaller of the two steps).
PROTOCOL_TOL = 5.0e-11

#: Step for the four-point SCF references. Truncation falls as ``h**4`` and
#: round-off as ``eps max|g| / h``, and at this step both sit far below the
#: floor these tests actually hit, so a *large* step is the right one.
SCF_STEP = 1.0e-2

#: Per-system agreement bounds for the SCF finite-difference references, each
#: about five times the worst deviation measured here. They are set by PySCF's
#: own gradient/Hessian consistency rather than by moist -- see the note at the
#: end of this docstring -- with the measured gas-phase floor in the comment.
SCF_TOL = {
    "water/sto-3g": 5.0e-8,          # measured 7.3e-9; gas phase alone 7.1e-9
    "water/def2-svp": 7.0e-7,        # measured 1.4e-7; gas phase alone 5.8e-8
    "fluoroacetate/sto-3g": 2.5e-7,  # measured 4.3e-8
}

#: Below this a reference carries no information.
VACUITY = 1.0e-4

COMPONENT_SETS = {
    "cpcm": [(ModelComponentCPCM, 80.0)],
    "cosmo": [(ModelComponentCOSMO, 12.0)],
    "cpcm+cosmo": [(ModelComponentCPCM, 4.0), (ModelComponentCOSMO, 12.0)],
    "pv": [(ModelComponentPV, 1.0e-3)],
}


#: An asymmetric geometry exercises off-diagonal atom and Cartesian blocks.
WATER = "O 0 0 -0.3893; H 0.7629 0 0.1947; H -0.7991 0.0953 0.2223"

#: Fluoroacetate, FCH2-COO(-), in a deliberately asymmetric non-stationary
#: geometry: no test here differentiates at a minimum, and a symmetric cavity
#: would risk the multi-branch projections the Hessian refuses.
FLUOROACETATE = (
    "C 0.000 0.000 0.000; F 1.180 0.560 0.420; "
    "H -0.420 -0.900 0.470; H -0.120 0.180 -1.060; "
    "C -0.870 1.060 0.680; O -1.980 1.300 0.130; O -0.480 1.600 1.760"
)

#: The systems the isodensity derivatives are pinned on. Water at STO-3G is the
#: fast default every component set runs through. The other two each add one
#: thing nothing else in the suite covers: def2-SVP is the only basis here with
#: d functions, so it is the only test of the AO derivative table beyond s and
#: p; fluoroacetate is the only polyatomic, charged and strongly asymmetric
#: case, and its minimal-basis anion SCF needs a raised iteration cap to reach
#: these tolerances at all.
SYSTEMS = {
    "water/sto-3g": dict(atom=WATER, basis="sto-3g", charge=0, max_cycle=None),
    "water/def2-svp": dict(atom=WATER, basis="def2-svp", charge=0, max_cycle=None),
    "fluoroacetate/sto-3g": dict(atom=FLUOROACETATE, basis="sto-3g", charge=-1, max_cycle=300),
}


def molecule(name):
    """One entry of :data:`SYSTEMS` as a PySCF molecule and its SCF cycle cap."""
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


def build(mol, specs):
    """A host and a model of one component set on the isodensity cavity."""
    host = PySCFHost(mol, scale=SCALE)
    components = [component(value) for component, value in specs]
    return host, SolvationModel(CavityDROPIsodensity(host, **CAVITY), components)


def mixed_direction(mol, dm, seed=11):
    """One normalised direction with a nuclear and a density part."""
    rng = np.random.default_rng(seed)
    nuclear = rng.normal(size=(3, mol.natm))
    nuclear /= np.linalg.norm(nuclear)
    density = rng.normal(size=dm.shape)
    density = 0.5 * (density + density.T)
    density /= np.linalg.norm(density)
    return np.asfortranarray(nuclear), density


def displaced(mol, dm, nuclear, density, step):
    """The molecule and density moved along a mixed direction by ``step``."""
    moved = mol.set_geom_(mol.atom_coords() + step * nuclear.T, unit="Bohr", inplace=False)
    return moved, dm + step * density


def worst(value, reference):
    """Worst absolute deviation, and the scale it is measured against."""
    value = np.asarray(value)
    reference = np.asarray(reference)
    return np.abs(value - reference).max(), np.abs(reference).max()


#: `(system, component set)` pairs for the protocol test. Water at STO-3G runs
#: every component set; the two costlier systems run one electrostatic and one
#: non-electrostatic set, which between them touch every channel of the
#: exchange without paying for the full product.
PROTOCOL_CASES = [("water/sto-3g", name) for name in COMPONENT_SETS] + [
    (system, name)
    for system in ("water/def2-svp", "fluoroacetate/sto-3g")
    for name in ("cpcm", "pv")
]


@pytest.mark.parametrize("system,name", PROTOCOL_CASES)
def test_protocol_matches_fd(system, name):
    """The whole exchange against differences of the analytic first derivatives.

    The sharpest statement that the *moist* half -- cavity plus model -- is
    right: the analytic gradient and Fock matrix are differenced along a mixed
    nuclear/density direction with no SCF in the loop, so a wrong cavity or
    component term shows up undiluted by electronic relaxation.
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


def response_channels(mol, dm, specs):
    """The base quantities of the response tangent, read back through one call."""
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
    """Every channel moist hands back is the derivative of what it says it is."""
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


# PV has no dense counterpart: the dense backend has no second-order channel
# for it, which ``test_hessian.py`` pins. The directional path covers it in the
# two finite-difference tests above.
@pytest.mark.parametrize("name", ["cpcm", "cosmo", "cpcm+cosmo"])
def test_scf_hessian_matches_dense(water, name):
    """The directional SCF Hessian reproduces the dense backend's."""
    specs = COMPONENT_SETS[name]

    def model_factory(host):
        return SolvationModel(
            CavityDROPIsodensity(host, **CAVITY),
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
    """A callable giving the reconverged solvated analytic gradient at a geometry."""

    def model_factory(host):
        return SolvationModel(
            CavityDROPIsodensity(host, **CAVITY),
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
    """The full Hessian carried to finite differences of SCF gradients.

    The directional path end to end: the solvated SCF is reconverged at every
    displaced geometry and its analytic gradient differenced, so electronic
    relaxation is in the reference as well as in the Hessian. Restricted to the
    three-atom systems -- the cost is 6N reconverged solvated SCF runs, which
    fluoroacetate answers through :func:`test_scf_hvp_matches_fd` instead.
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
    """Fluoroacetate end to end, one direction instead of the whole block.

    The same statement as :func:`test_scf_hessian_matches_fd` -- an analytic
    Hessian against differences of reconverged solvated SCF gradients, so the
    coupled-perturbed solve is in the reference too -- contracted with a single
    nuclear direction. That is two reconverged SCF runs rather than 6N, which
    is what makes a seven-atom, charged, minimal-basis anion affordable here at
    all; a direction with a component on every atom still touches every column
    of the block.
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
    """A host that raises inside the exchange comes back with its own error."""

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
    """The cached reductions are keyed on their operand's value, not its identity.

    ``_PotentialIntegrals`` outlives a single Hessian-vector call -- it is kept
    in the solve's ``shared`` dict -- so the reductions it caches are read by the
    nuclear block and again by every coupled-perturbed iteration. Neither the
    converged charges nor the converged density change within a solve, which is
    what makes the cache worth having; but a cache that answered for an operand
    that *had* changed would return a wrong Hessian silently rather than fail,
    so the miss is worth pinning as well as the hit.
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

    # A hit is the same object: the cache is load-bearing, not merely correct
    assert ints.density_contracted(first) is ints.density_contracted(first.copy())
    assert ints.charge_weighted(charges) is ints.charge_weighted(charges.copy())

    # and every miss, including a return to an operand seen before, rebuilds
    for dm in (first, second, first):
        assert np.array_equal(
            ints.density_contracted(dm).field, _contract("cuvi,uv->ci", ints.ip2, dm)
        )
    for q in (charges, other, charges):
        assert np.array_equal(
            ints.charge_weighted(q).bra, _contract("i,cuvi->cuv", q, ints.ip1)
        )

    # a caller mutating its own array afterwards must not be answered from cache
    held = first.copy()
    ints.density_contracted(held)
    held += 1.0
    assert np.array_equal(
        ints.density_contracted(held).field, _contract("cuvi,uv->ci", ints.ip2, held)
    )


def test_direction_threads_do_not_change_the_result(water):
    """The per-direction loops are threaded, and the answer does not know it.

    Both host callbacks and the completion run their direction loop over a pool
    sized by PySCF's thread count. The directions share only read-only operands
    and write disjoint slices, so the schedule cannot reach the arithmetic and
    the results must be equal *bit for bit*, not merely close -- a tolerance
    here would hide exactly the race it is meant to catch, which is why this is
    the one comparison in the suite made with ``array_equal``.

    Only the Python layer is switched between the two runs: the tangent and the
    response are taken once, under one thread, and the three loops are then
    replayed against them. Re-running the Fortran under a second thread count
    would test its reductions rather than these loops.
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
    """What a missing half of the exchange, a stale evaluation and a bad method do."""
    with lib.with_omp_threads(1):
        dm = converged_rhf(water).make_rdm1()
        probe = np.zeros((3, water.natm, 1), order="F")
        probe[0, 0, 0] = 1.0

        # No host object at all: the level set is the host's, so every one of
        # its tangents would be an unphysical zero
        host, model = build(water, COMPONENT_SETS["cpcm"])
        evaluation = model.evaluate(coupling=host.coupling(dm))
        with pytest.raises(RuntimeError, match="level set is the host's"):
            evaluation.hvp_coupled(probe)

        # A host that answers for the level set but not for its own potential,
        # which an electrostatic component on a host potential has to have
        points = model.cavity.snapshot().xyz
        complete = DirectionalHost(host, dm, points.T, probe)

        class LevelSetOnly:
            def level_set_tangent(self, first, dirs, xyz):
                return complete.level_set_tangent(first, dirs, xyz)

        with pytest.raises(RuntimeError, match="field tangent callback"):
            evaluation.hvp_coupled(probe, LevelSetOnly(), level_set=True)

        # A model that reads no host potential needs no field callback, and
        # produces no charges to differentiate
        pv_host, pv_model = build(water, COMPONENT_SETS["pv"])
        pv = pv_model.evaluate(coupling=pv_host.coupling(dm))
        pv_points = pv_model.cavity.snapshot().xyz
        pv_tangent_host = DirectionalHost(pv_host, dm, pv_points.T, probe, potential=False)
        columns, tangent = pv.hvp_coupled(probe, pv_tangent_host, level_set=False)
        assert columns.shape == probe.shape
        assert tangent.w_value is None, "level-set channels appeared unasked"
        assert np.abs(tangent.xyz).max() > 0.0
        assert np.abs(tangent.surface_charge).max() == 0.0

        # They are on by default for a density-defined cavity, which is the
        # only kind whose level-set weights a host has any use for
        _, tangent = pv.hvp_coupled(probe, pv_tangent_host)
        assert tangent.w_value is not None

        model.evaluate(coupling=host.coupling(dm))
        with pytest.raises(RuntimeError, match="superseded"):
            evaluation.hvp_coupled(probe, complete, level_set=True)

    mean_field = converged_rhf(water)
    with pytest.raises(ValueError, match="unknown Hessian method"):
        from .hessian import rhf_hessian

        rhf_hessian(mean_field, model, host, method="bogus")
