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
