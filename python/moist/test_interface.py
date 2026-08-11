import numpy as np
import pytest
from pytest import approx, raises

from moist import library
from moist.interface import (
    ArrayCoupling,
    CPCM,
    Cavity,
    CouplingTransaction,
    DROPCavitySnapshot,
    DROPCavity,
    Electrostatics,
    Evaluation,
    GeneralSolvationModel,
    PV,
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
    cavity = DROPCavity(nleb=26)
    cavity.update(Structure(numbers, positions))
    result = cavity.cavity

    assert len(result.a) == 70
    assert np.asarray(result.a).sum() == approx(151.6278477636, rel=1e-10)
    assert result.volume == approx(173.4128767106, rel=1e-10)


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


def test_general_model_iterates_cpcm_and_pv_components() -> None:
    """A heterogeneous component list shares one authoritative live cavity."""

    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    pressure = 2.5e-4
    model = GeneralSolvationModel(
        DROPCavity(nleb=26),
        [CPCM(32.0), PV(pressure)],
    )
    model.update(structure)

    energy, charges = model.solve(np.zeros(model.ngrid))
    potential = model.get_potential()

    assert energy == approx(pressure * model.cavity.volume, abs=2.0e-13)
    np.testing.assert_array_equal(charges, np.zeros(model.ngrid))
    assert potential.w_umol.shape == (model.ngrid,)
    assert potential.w_lsf0.shape == (model.ngrid,)
    assert potential.w_lsf1.shape == (3, model.ngrid)
    assert potential.w_lsf2.shape == (3, 3, model.ngrid)


def test_general_model_evaluates_a_complete_array_coupling() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    pressure = 2.5e-4
    model = GeneralSolvationModel(
        DROPCavity(nleb=26),
        [CPCM(32.0), PV(pressure)],
    )
    coupling = ArrayCoupling(
        structure,
        electrostatics=lambda cavity, _trace: Electrostatics(
            np.zeros(cavity.ngrid)
        ),
    )

    result = model.evaluate(coupling=coupling)

    assert isinstance(result, Evaluation)
    assert isinstance(result.cavity, DROPCavitySnapshot)
    assert isinstance(model.cavity, Cavity)
    assert result.energy == approx(pressure * result.cavity.volume, abs=2.0e-13)
    np.testing.assert_array_equal(result.charges, np.zeros(result.cavity.ngrid))
    assert result.fock is None
    assert result.gradient.shape == (3, len(structure))
    assert result.cavity.grid_points is result.cavity.xyz
    assert result.cavity.grid_areas is result.cavity.a
    assert result.potential.level_set_value_weights is result.potential.w_lsf0
    assert not result.cavity.xyz.flags.writeable
    assert not result.potential.w_lsf0.flags.writeable
    assert not result.gradient.flags.writeable


def test_solvation_model_is_the_canonical_constructor() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    coupling = ArrayCoupling(structure)
    model = SolvationModel(DROPCavity(nleb=26), [PV(2.5e-4)])

    result = model.evaluate(structure, coupling=coupling)

    assert GeneralSolvationModel is SolvationModel
    assert result.energy == approx(2.5e-4 * result.cavity.volume, abs=2.0e-13)


def test_evaluation_rejects_a_structure_from_another_coupling() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    other = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.8], [0.0, 0.0, 0.8]]),
    )
    model = SolvationModel(DROPCavity(nleb=26), [PV(2.5e-4)])

    with raises(ValueError, match="does not match"):
        model.evaluate(structure, coupling=ArrayCoupling(other))


def test_evaluation_results_use_irreversibly_read_only_buffers() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    result = SolvationModel(DROPCavity(nleb=26), [PV(2.5e-4)]).evaluate(structure)

    with raises(ValueError, match="WRITEABLE"):
        result.cavity.xyz.setflags(write=True)
    with raises(ValueError, match="WRITEABLE"):
        result.potential.w_lsf0.setflags(write=True)


def test_failed_coupling_preparation_invalidates_the_model() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )

    class FailingCoupling(SolvationCoupling):
        channels = frozenset({"electrostatics"})

        @property
        def structure(self) -> Structure:
            return structure

        def prepare(self, transaction: CouplingTransaction) -> None:
            def fail_after_first_supply(cavity, trace):
                if trace is not None:
                    raise RuntimeError("host response failed")
                return Electrostatics(np.zeros(cavity.ngrid))

            transaction.exchange_electrostatics(fail_after_first_supply)

    model = SolvationModel(DROPCavity(nleb=26), [CPCM(32.0)])

    with raises(RuntimeError, match="host response failed"):
        model.evaluate(coupling=FailingCoupling())
    with raises(RuntimeError, match="successfully updated"):
        _ = model.energy
    with raises(RuntimeError, match="successfully updated"):
        model.potential()
    with raises(RuntimeError, match="successfully updated"):
        model.cavity.snapshot()


def test_general_model_rejects_missing_coupling_channels() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    model = GeneralSolvationModel(DROPCavity(nleb=26), [CPCM(32.0)])

    with raises(ValueError, match="electrostatics"):
        model.evaluate(structure)


def test_evaluation_gradient_rejects_a_superseded_model_state() -> None:
    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    model = GeneralSolvationModel(DROPCavity(nleb=26), [PV(2.5e-4)])
    first = model.evaluate(structure)
    model.evaluate(structure)

    with raises(RuntimeError, match="superseded"):
        _ = first.gradient


def test_general_model_update_invalidates_supplied_electrostatics() -> None:
    """Every geometry update requires a fresh external potential trace."""

    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    model = GeneralSolvationModel(DROPCavity(nleb=26), [CPCM(32.0)])
    model.update(structure)
    model.solve(np.zeros(model.ngrid))

    previous_ngrid = model.ngrid
    model.update(structure)
    assert model.ngrid == previous_ngrid

    with raises(RuntimeError, match="phi not supplied"):
        model.get_energy()
    with raises(RuntimeError, match="phi not supplied"):
        model.get_potential()
    with raises(RuntimeError, match="phi not supplied"):
        model.get_gradient(len(structure))


def test_borrowed_model_cavity_rejects_standalone_updates() -> None:
    """A model-owned cavity may be inspected but not rebuilt out of band."""

    structure = Structure(
        np.array([1, 1], dtype=np.int32),
        np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]]),
    )
    model = GeneralSolvationModel(DROPCavity(nleb=26), [PV(1.0e-4)])
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
