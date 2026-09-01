from typing import Callable

import numpy as np
import pytest
from pytest import approx, raises

import moist
from moist import library
from moist.interface import (
    ArrayCoupling,
    Cavity,
    CavityDROP,
    CavityDROPCFC,
    CavityDROPSvdW,
    CavityISwiG,
    CavitySnapshot,
    CavitySnapshotDROP,
    CouplingChannel,
    CouplingTransaction,
    Electrostatics,
    Evaluation,
    GeneralSolvationModel,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentPV,
    PCMSolver,
    SolvationCoupling,
    SolvationModel,
    Structure,
)


@pytest.fixture
def numbers() -> np.ndarray:
    return np.array([8, 1, 1])


@pytest.fixture
def positions() -> np.ndarray:
    return np.array(
        [
            [ 0.00000000000000, 0.00000000000000, -0.73578586109551],
            [ 1.44183152868459, 0.00000000000000,  0.36789293054775],
            [-1.44183152868459, 0.00000000000000,  0.36789293054775],
        ]
    )


@pytest.fixture
def diatomic() -> Callable[..., Structure]:
    """Factory for an H2 structure with the atoms ``bond`` bohr apart along z.

    The model and coupling tests only need *a* valid structure, and two atoms
    keep the surface small enough that building several models per test stays
    cheap -- the ``numbers``/``positions`` fixtures above are water, which those
    tests would pay for without testing anything more.
    """

    def build(bond: float = 1.4) -> Structure:
        half = 0.5 * bond
        return Structure(
            np.array([1, 1], dtype=np.int32),
            np.array([[0.0, 0.0, -half], [0.0, 0.0, half]]),
        )

    return build


def test_default_drop_surface_matches_the_svdw_defaults(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """The default cavity built through this API is the library's default surface.

    The C entry point takes the SvdW shape parameters as nullable pointers and
    has to supply its own fallbacks.  Those fallbacks once drifted from the
    defaults declared on ``moist_cavity_drop_lsf_svdw_param_type`` (k 2.0 vs
    5.5, 2b 1.0 vs 0.0, 3b 1.0 vs 3.0), so a Fortran caller and a Python caller
    asking for "the default cavity" got measurably different surfaces -- 23%
    apart in total area for this molecule.

    Nothing caught it: every other cavity test is a finite difference or an
    internal consistency check, and those follow whichever surface is built.
    Only an absolute value can see a changed default, which is what this is.
    Regenerate deliberately if the default surface is ever meant to change.
    """
    cavity = CavityDROP(nleb=26)
    cavity.update(Structure(numbers, positions))
    result = cavity.cavity

    assert len(result.a) == 70
    assert np.asarray(result.a).sum() == approx(151.6278477636, rel=1e-10)
    assert result.volume == approx(173.4128767106, rel=1e-10)


@pytest.mark.parametrize(
    "cavity,snapshot_type",
    [
        (
            CavityDROPSvdW(
                nleb=26,
                blend_k=5.5,
                blend_1b=1.0,
                blend_2b=0.0,
                blend_3b=3.0,
                debug=False,
                verbosity=0,
                do_fine=True,
                tolerance=1.0e-10,
                proj_maxiter=200,
                proj_level=2,
                branch_weight_s=0.08,
                rho_grid_h=0.8,
                wleb_prune_level=1,
            ),
            CavitySnapshotDROP,
        ),
        (
            CavityDROPCFC(
                nleb=26,
                a1=-15.0,
                a2=-9.0,
                c=5.0,
                m=4,
                debug=False,
                verbosity=0,
                do_fine=True,
                tolerance=1.0e-10,
                proj_maxiter=200,
                proj_level=2,
                branch_weight_s=0.08,
                rho_grid_h=0.8,
                wleb_prune_level=1,
            ),
            CavitySnapshotDROP,
        ),
        (
            CavityISwiG(
                nleb=26,
                cut_a=0.0,
                cut_f=1.0e-10,
                debug=False,
                verbosity=0,
            ),
            CavitySnapshot,
        ),
    ],
    ids=("svdw-drop", "cfc-drop", "iswig"),
)
def test_public_cavity_types_accept_their_optional_arguments(
    cavity: Cavity,
    snapshot_type,
    numbers: np.ndarray,
    positions: np.ndarray,
) -> None:
    cavity.update(Structure(numbers, positions))
    result = cavity.snapshot()

    assert isinstance(result, snapshot_type)
    assert result.ngrid > 0
    assert result.nsph == len(numbers)
    assert np.isfinite(result.area)
    assert np.isfinite(result.volume)


@pytest.fixture
def branching_cross() -> Structure:
    """Five-carbon cross whose concave seams branch at ``proj_level=7``.

    The same fixture the Fortran suite uses (``get_test_cross``); the default
    projection level finds no second solution on it, which is why the branching
    tests below raise the level rather than the geometry.
    """

    aatoau = 1.8897261246257702
    return Structure(
        np.array([6, 6, 6, 6, 6]),
        np.array(
            [
                [0.00, 4.21, 0.00],
                [0.00, 0.00, 4.22],
                [0.00, -4.18, 0.00],
                [0.00, 0.00, -4.15],
                [0.02, 0.10, -0.20],
            ]
        )
        * aatoau,
    )


def test_cavity_declares_its_own_results(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """Every declared field describes itself well enough to be read blind."""
    cavity = CavityDROPSvdW()
    cavity.update(Structure(numbers, positions))

    fields = cavity.fields()
    assert fields

    results = cavity.results()
    assert set(results) == {field.name for field in fields}

    for field in fields:
        value = results[field.name]
        assert np.shape(value) == field.shape
        assert np.asarray(value).dtype == field.dtype
        assert cavity.describe(field.name)


def test_named_results_agree_with_the_snapshot(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """The typed snapshot is a view of the same declarations, not a second read."""
    cavity = CavityDROPSvdW()
    cavity.update(Structure(numbers, positions))

    snapshot = cavity.snapshot()
    results = cavity.results()

    assert results["ngrid"] == snapshot.ngrid
    assert results["nsph"] == snapshot.nsph
    assert results["area"] == approx(snapshot.area)
    assert np.array_equal(results["owner"], snapshot.owner)
    assert np.array_equal(results["xyz"], snapshot.xyz)
    assert results["xyz"].shape == (3, snapshot.ngrid)
    assert np.array_equal(results["numbering"], snapshot.numbering)
    assert np.array_equal(results["branch_count"], snapshot.branch_count)


def test_uncomputed_results_are_absent_rather_than_zero(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """A property that was never requested must not read back as zeros."""
    structure = Structure(numbers, positions)

    plain = CavityDROPSvdW()
    plain.update(structure)
    assert "k1" not in {field.name for field in plain.fields()}
    with raises(KeyError, match="k1"):
        plain.get("k1")

    fine = CavityDROPSvdW(do_fine=True)
    fine.update(structure)
    curvature = fine.get("k1")
    assert curvature.shape == (fine.ngrid,)
    assert np.isfinite(curvature).all()

    with raises(KeyError, match="not_a_field"):
        fine.get("not_a_field")


def test_branching_is_read_from_the_cavity(branching_cross: Structure) -> None:
    """Branch data comes from moist's own arrays, not from unpacking an id."""
    cavity = CavityDROPSvdW(proj_level=7)
    cavity.update(branching_cross)
    snapshot = cavity.snapshot()

    assert snapshot.branch.min() == 1
    assert snapshot.branch.max() > 1, "fixture stopped producing branches"
    assert snapshot.branched.sum() > 0
    assert np.array_equal(snapshot.branched, snapshot.branch_count > 1)

    # numbering is the packing of the two, so the arrays and the id agree.
    base = snapshot.nsph * cavity.get("num_leb")
    assert np.array_equal(
        snapshot.numbering, snapshot.anchor_id + base * (snapshot.branch - 1)
    )

    # Every point in a branched group reports the same group size.
    for anchor in np.unique(snapshot.anchor_id[snapshot.branched]):
        group = snapshot.anchor_id == anchor
        assert snapshot.branch_count[group].min() == group.sum()
        assert np.array_equal(np.unique(snapshot.branch[group]), snapshot.branch[group])


def test_named_results_reach_every_cavity_type(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """The field API is a cavity feature, not a DROP one."""
    cavity = CavityISwiG()
    cavity.update(Structure(numbers, positions))

    names = {field.name for field in cavity.fields()}
    assert {"xyz", "a", "owner", "numbering"} <= names
    assert np.array_equal(cavity.get("owner"), cavity.snapshot().owner)


def test_named_results_need_a_built_cavity() -> None:
    cavity = CavityDROPSvdW()
    with raises(RuntimeError, match="not been successfully updated"):
        cavity.fields()
    with raises(RuntimeError, match="not been successfully updated"):
        cavity.get("xyz")


def test_cavity_drop_is_the_svdw_surface() -> None:
    """The unqualified DROP name is the SvdW cavity, not a distinct surface."""
    assert CavityDROP is CavityDROPSvdW
    assert moist.CavityDROP is moist.CavityDROPSvdW


@pytest.mark.parametrize("cavity_type", [CavityDROPSvdW, CavityDROPCFC])
def test_drop_constructor_controls_are_validated_natively(cavity_type) -> None:
    with raises(RuntimeError, match="wleb_prune_level.*0-6"):
        cavity_type(wleb_prune_level=7)


def test_cavity_specific_options_reach_the_native_implementations(
    numbers: np.ndarray,
    positions: np.ndarray,
) -> None:
    structure = Structure(numbers, positions)

    def surface(cavity: Cavity) -> CavitySnapshot:
        cavity.update(structure)
        return cavity.snapshot()

    svdw_default = surface(CavityDROPSvdW(nleb=26))
    svdw_custom = surface(CavityDROPSvdW(nleb=26, blend_k=4.0))
    assert svdw_custom.area != approx(svdw_default.area, rel=1.0e-6)

    cfc_default = surface(CavityDROPCFC(nleb=26))
    cfc_custom = surface(
        CavityDROPCFC(
            nleb=26,
            a1=-12.0,
            a2=-8.0,
            c=4.0,
            m=4,
        )
    )
    assert cfc_custom.area != approx(cfc_default.area, rel=1.0e-6)

    iswig_default = surface(CavityISwiG(nleb=26))
    iswig_custom = surface(CavityISwiG(nleb=26, cut_f=1.0e-2))
    assert iswig_custom.ngrid < iswig_default.ngrid


def test_structure(numbers: np.ndarray, positions: np.ndarray) -> None:
    with raises(ValueError, match="positions must have shape"):
        Structure(np.array([1, 1]), positions)

    with raises(ValueError, match="positions must have shape"):
        Structure(numbers, np.random.default_rng().random(7))

    structure = Structure(numbers, positions)

    with raises(ValueError, match="positions must have shape"):
        structure.update(np.random.default_rng().random(7))

    with raises(ValueError, match="lattice must have shape"):
        structure.update(positions, np.random.default_rng().random(7))

    with raises(ValueError, match="positions must have shape"):
        Structure(np.array([1, 1]), np.zeros((3, 2)))

    with raises(ValueError, match="numbers must have shape"):
        Structure(numbers.reshape(1, -1), positions)


def test_general_model_iterates_cpcm_and_pv_components(diatomic) -> None:
    """A heterogeneous component list shares one authoritative live cavity."""

    structure = diatomic()
    pressure = 2.5e-4
    model = GeneralSolvationModel(
        CavityDROP(nleb=26),
        [ModelComponentCPCM(32.0), ModelComponentPV(pressure)],
    )
    model.update(structure)

    energy, charges = model.solve(np.zeros(model.ngrid))
    response = model.response()

    assert energy == approx(pressure * model.cavity.volume, abs=2.0e-13)
    np.testing.assert_array_equal(charges, np.zeros(model.ngrid))
    assert response.electrostatics.surface_charge.shape == (model.ngrid,)
    assert response.lsf.w_value.shape == (model.ngrid,)
    assert response.lsf.w_gradient.shape == (3, model.ngrid)
    assert response.lsf.w_hessian.shape == (3, 3, model.ngrid)


def test_cosmo_is_a_standalone_pcm_component() -> None:
    cpcm = ModelComponentCPCM(32.0, solver="lu")
    cosmo = ModelComponentCOSMO(32.0, solver="lu")

    assert cpcm.epsilon == cosmo.epsilon == 32.0
    assert cpcm.solver is cosmo.solver is PCMSolver.LU
    assert cpcm.coupling_channels == cosmo.coupling_channels
    assert moist.CPCMSolver is moist.PCMSolver


def test_cosmo_rejects_an_invalid_dielectric() -> None:
    with raises(RuntimeError, match="Dielectric constant must be >= 1"):
        ModelComponentCOSMO(0.5)


def test_cosmo_uses_its_own_dielectric_scaling(
    numbers: np.ndarray,
    positions: np.ndarray,
) -> None:
    """COSMO composes like CPCM but applies its distinct screening factor."""
    structure = Structure(numbers, positions)
    epsilon = 32.0
    results = {}

    for component_type in (ModelComponentCPCM, ModelComponentCOSMO):
        model = SolvationModel(CavityDROPSvdW(nleb=26), [component_type(epsilon)])
        model.update(structure)
        phi = np.linspace(-0.2, 0.3, model.ngrid)
        results[component_type] = model.solve(phi)

    cpcm_energy, cpcm_charges = results[ModelComponentCPCM]
    cosmo_energy, cosmo_charges = results[ModelComponentCOSMO]
    factor_ratio = ((epsilon - 1.0) / (epsilon + 0.5)) / (
        (epsilon - 1.0) / epsilon
    )

    assert cosmo_energy == approx(factor_ratio * cpcm_energy, rel=1.0e-13)
    np.testing.assert_allclose(
        cosmo_charges,
        factor_ratio * cpcm_charges,
        rtol=1.0e-13,
        atol=1.0e-14,
    )


def test_general_model_evaluates_a_complete_array_coupling(diatomic) -> None:
    structure = diatomic()
    pressure = 2.5e-4
    model = GeneralSolvationModel(
        CavityDROP(nleb=26),
        [ModelComponentCPCM(32.0), ModelComponentPV(pressure)],
    )
    coupling = ArrayCoupling(
        structure,
        # qefield is required by the external-potential gradient path; with
        # phi = 0 the correct value is an explicit zero field.
        electrostatics=lambda cavity, _trace: Electrostatics(
            np.zeros(cavity.ngrid),
            qefield=np.zeros((3, cavity.ngrid), order="F"),
        ),
    )

    result = model.evaluate(coupling=coupling)

    assert isinstance(result, Evaluation)
    assert isinstance(result.cavity, CavitySnapshotDROP)
    assert isinstance(model.cavity, Cavity)
    assert result.energy == approx(pressure * result.cavity.volume, abs=2.0e-13)
    # CPCM is present and was handed phi = 0, so zero charges is a genuine
    # result rather than an absent channel.
    np.testing.assert_array_equal(result.charges, np.zeros(result.cavity.ngrid))
    assert result.fock is None
    assert result.gradient.shape == (3, len(structure))
    assert not result.cavity.xyz.flags.writeable
    assert not result.response.lsf.w_value.flags.writeable
    assert not result.gradient.flags.writeable


def test_solvation_model_is_the_canonical_constructor(diatomic) -> None:
    structure = diatomic()
    coupling = ArrayCoupling(structure)
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentPV(2.5e-4)])

    result = model.evaluate(structure, coupling=coupling)

    assert GeneralSolvationModel is SolvationModel
    assert result.energy == approx(2.5e-4 * result.cavity.volume, abs=2.0e-13)


def test_evaluation_rejects_a_structure_from_another_coupling(diatomic) -> None:
    structure = diatomic()
    other = diatomic(1.6)
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentPV(2.5e-4)])

    with raises(ValueError, match="does not match"):
        model.evaluate(structure, coupling=ArrayCoupling(other))


def test_evaluation_results_use_irreversibly_read_only_buffers(diatomic) -> None:
    structure = diatomic()
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentPV(2.5e-4)])
    result = model.evaluate(structure)

    with raises(ValueError, match="WRITEABLE"):
        result.cavity.xyz.setflags(write=True)
    with raises(ValueError, match="WRITEABLE"):
        result.response.lsf.w_value.setflags(write=True)


def test_failed_coupling_preparation_invalidates_the_model(diatomic) -> None:
    structure = diatomic()

    class FailingCoupling(SolvationCoupling):
        channels = frozenset({CouplingChannel.ELECTROSTATICS})

        @property
        def structure(self) -> Structure:
            return structure

        def prepare(self, transaction: CouplingTransaction) -> None:
            def fail_after_first_supply(cavity, trace):
                if trace is not None:
                    raise RuntimeError("host response failed")
                return Electrostatics(np.zeros(cavity.ngrid))

            transaction.exchange_electrostatics(fail_after_first_supply)

    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])

    with raises(RuntimeError, match="host response failed"):
        model.evaluate(coupling=FailingCoupling())
    with raises(RuntimeError, match="successfully updated"):
        _ = model.energy
    with raises(RuntimeError, match="successfully updated"):
        model.response()
    with raises(RuntimeError, match="successfully updated"):
        model.cavity.snapshot()


def test_general_model_rejects_missing_coupling_channels(diatomic) -> None:
    structure = diatomic()
    model = GeneralSolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])

    with raises(ValueError, match="electrostatics"):
        model.evaluate(structure)


def test_evaluation_gradient_rejects_a_superseded_model_state(diatomic) -> None:
    structure = diatomic()
    model = GeneralSolvationModel(CavityDROP(nleb=26), [ModelComponentPV(2.5e-4)])
    first = model.evaluate(structure)
    model.evaluate(structure)

    with raises(RuntimeError, match="superseded"):
        _ = first.gradient


def test_general_model_update_invalidates_supplied_electrostatics(diatomic) -> None:
    """Every geometry update requires a fresh external potential trace."""

    structure = diatomic()
    model = GeneralSolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    model.solve(np.zeros(model.ngrid))

    previous_ngrid = model.ngrid
    model.update(structure)
    assert model.ngrid == previous_ngrid

    with raises(RuntimeError, match="phi not supplied"):
        model.get_energy()
    with raises(RuntimeError, match="phi not supplied"):
        model.response()
    with raises(RuntimeError, match="phi not supplied"):
        model.get_gradient(len(structure))


def test_borrowed_model_cavity_rejects_standalone_updates(diatomic) -> None:
    """A model-owned cavity may be inspected but not rebuilt out of band."""

    structure = diatomic()
    model = GeneralSolvationModel(CavityDROP(nleb=26), [ModelComponentPV(1.0e-4)])
    model.update(structure)
    original_area = model.cavity.area

    with raises(RuntimeError, match="model-owned cavity"):
        model.cavity.update(structure)
    with pytest.deprecated_call():
        cavity_handle = model.cavity_handle
    with raises(RuntimeError, match="borrowed cavity"):
        library.update_cavity(cavity_handle, structure._mol)
    with pytest.deprecated_call():
        cavity_handle = model.cavity_handle
    with raises(RuntimeError, match="borrowed cavity"):
        library.error_check(library.lib.moist_update_drop_cavity)(
            cavity_handle.handle,
            structure._mol.handle,
            library.ffi.NULL,
        )

    borrowed = library.error_check(library.lib.moist_get_solvation_model_cavity)(
        model._model.handle
    )
    borrowed_ref = library.ffi.new("moist_cavity *")
    borrowed_ref[0] = borrowed
    library.lib.moist_delete_cavity(borrowed_ref)

    assert borrowed_ref[0] == library.ffi.NULL
    assert model.cavity.area == original_area
