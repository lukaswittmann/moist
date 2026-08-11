import numpy as np
import pytest
from pytest import approx, raises

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
