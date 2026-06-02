import numpy as np
import pytest
from pytest import approx, raises

from moist.interface import SolvationModel, Structure


@pytest.fixture
def numbers() -> np.ndarray:
    return np.array([8, 1, 1])


@pytest.fixture
def positions() -> np.ndarray:
    return np.array(
        [
            [0.00000000000000, 0.00000000000000, -0.73578586109551],
            [1.44183152868459, 0.00000000000000, 0.36789293054775],
            [-1.44183152868459, 0.00000000000000, 0.36789293054775],
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


# def test_gems_model_invalid_solvent() -> None:
#     with raises(RuntimeError):
#         SolvationModelGEMS("definitely-not-a-solvent")


# def test_gems_model_water(numbers: np.ndarray, positions: np.ndarray) -> None:
#     structure = Structure(numbers, positions)
#     model = SolvationModelGEMS("water")

#     with raises(RuntimeError, match="updated before requesting the energy"):
#         model.get_energy()

#     with raises(RuntimeError, match="updated before requesting the cavity"):
#         _ = model.cavity

#     model.update(structure)

#     energy = model.get_energy()
#     cavity = model.cavity

#     assert approx(energy, abs=1.0e-12) == 0.0064864129243683965
#     assert cavity.ngrid > 0
#     assert cavity.nsph == len(numbers)
#     assert cavity.xyz.shape == (3, cavity.ngrid)
#     assert cavity.a.shape == (cavity.ngrid,)
#     assert cavity.owner.shape == (cavity.ngrid,)
#     assert cavity.converged.shape == (cavity.ngrid,)
#     assert cavity.radii.shape == (cavity.nsph,)
#     assert cavity.asph.shape == (cavity.nsph,)
#     assert cavity.normal0.shape == (3, cavity.ngrid)
#     assert cavity.wleb.shape == (cavity.ngrid,)
#     assert cavity.r_iI0.shape == (cavity.ngrid,)
#     assert cavity.f.shape == (cavity.ngrid,)
#     assert cavity.rho.shape == (cavity.ngrid,)
#     assert cavity.area > 0.0
#     assert cavity.volume > 0.0
#     assert np.all(cavity.owner >= 0)
#     assert np.all(cavity.owner < cavity.nsph)


# def test_gems_model_base_type(numbers: np.ndarray, positions: np.ndarray) -> None:
#     structure = Structure(numbers, positions)
#     model = SolvationModelGEMS("water")

#     assert isinstance(model, SolvationModel)

#     model.update(structure)
#     assert isinstance(model.get_energy(), float)
