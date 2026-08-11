import numpy as np
import pytest
from pytest import approx, raises

from moist import library
from moist.interface import (
    CPCM,
    DROPCavity,
    GeneralSolvationModel,
    PV,
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
    with raises(ValueError, match="Dimension missmatch"):
        Structure(np.array([1, 1]), positions)

    with raises(ValueError, match="Expected tripels"):
        Structure(numbers, np.random.default_rng().random(7))

    structure = Structure(numbers, positions)

    with raises(ValueError, match="Dimension missmatch for positions"):
        structure.update(np.random.default_rng().random(7))

    with raises(ValueError, match="Invalid lattice provided"):
        structure.update(positions, np.random.default_rng().random(7))


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

    with raises(RuntimeError, match="borrowed cavity"):
        library.update_cavity(model.cavity_handle, structure._mol)
    with raises(RuntimeError, match="borrowed cavity"):
        library.error_check(library.lib.moist_update_drop_cavity)(
            model.cavity_handle.handle,
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
