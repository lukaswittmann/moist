"""The parameter_adjoint response item: presence, option gating, and dE/dp by finite differences."""

import math

import numpy as np
import pytest
from pytest import approx, raises

import moist
from moist import (
    CavityDROP, DROPParameters, ModelComponentCPCM, ModelComponentLPB,
    ModelComponentNonlocal, ModelComponentSLIC, ParameterAdjointResponse,
    PCMParameters, GaussianPotentialRequest, SolvationModel, Structure,
)

CONTEXT = moist.Context()
_Z = np.array([1.0, 1.0])
_R = np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]])

#: Component factory taking the parameters by name, and the adjoint fields the model has
_MODELS = {
    "cpcm": (
        lambda p, o: ModelComponentCPCM(p["eps"], parameters=o),
        ("dg_deps",),
    ),
    "lpb": (
        lambda p, o: ModelComponentLPB(p["eps"], p["kappa"], parameters=o),
        ("dg_deps", "dg_dkappa"),
    ),
    "nonlocal": (
        lambda p, o: ModelComponentNonlocal(p["eps"], p["eps_inf"], p["lambda"], p["eps_solute"], parameters=o),
        ("dg_deps", "dg_deps_inf", "dg_dlambda", "dg_deps_solute"),
    ),
    "slic": (
        lambda p, o: ModelComponentSLIC(p["eps"], p["alpha"], p["beta"], p["gamma"], p["mu"],
                                        p["phi_static"], parameters=o),
        ("dg_deps", "dg_dalpha", "dg_dbeta", "dg_dgamma", "dg_dmu", "dg_dphi_static"),
    ),
}
_VALUES = {
    "eps": 30.0, "kappa": 0.1, "eps_inf": 2.0, "lambda": 3.0, "eps_solute": 1.5,
    "alpha": 0.3, "beta": 1.2, "gamma": 0.4, "mu": 0.1, "phi_static": 0.02,
}
#: Fraction of the value used as the central-difference step
_STEP = 1.0e-4
#: Relative tolerance on dE/dp against the central difference; measured worst case 1e-8
_RTOL = 1.0e-6


def _potential(point, xi):
    """Potential of the nuclei on unit Gaussians of exponent ``xi`` at ``point`` (ngrid, 3), with derivatives."""
    d = point[:, None, :] - _R[None, :, :]
    r = np.linalg.norm(d, axis=2)
    x = xi[:, None]
    root = np.sqrt(x)
    erf = np.vectorize(math.erf)(root * r)
    gauss = 2 * root / np.sqrt(np.pi) * np.exp(-x * r**2)
    f = erf / r
    f1 = gauss / r - erf / r**2
    f2 = -2 * x * gauss - 2 * gauss / r**2 + 2 * erf / r**3
    dd = np.einsum("iac,iad->iacd", d, d)
    h = (f2 - f1 / r)[..., None, None] * dd / r[..., None, None]**2 + (f1 / r)[..., None, None] * np.eye(3)
    dxi = np.exp(-x * r**2) / np.sqrt(np.pi * x)
    return {
        "phi": np.einsum("a,ia->i", _Z, f),
        "dphi_dr": np.einsum("a,iac->ic", _Z, (f1 / r)[..., None] * d),
        "d2phi_dr2": np.einsum("a,iacd->icd", _Z, h),
        "dphi_dxi": np.einsum("a,ia->i", _Z, dxi),
        "d2phi_dr_dxi": np.einsum("a,iac->ic", _Z, -2 * x[..., None] * dxi[..., None] * d),
    }


def _run(component):
    """Return (energy, response) of one fixed-geometry solve with the nuclear potential on the Gaussian grid."""
    model = SolvationModel(CONTEXT, CavityDROP(parameters=DROPParameters(nleb=26), context=CONTEXT),
                           [component])
    model.update(Structure(np.array([1, 1], dtype=np.int32), _R, charge=1.0))
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    energy = np.array(0.0)

    def answer():
        for request in coupling:
            assert isinstance(request, GaussianPotentialRequest), request.name
            outputs = _potential(model.cavity.xyz, model.cavity.xi0)
            coupling.answer(**{name: outputs[name] for name in request.missing})

    answer()
    model.get_energy(coupling, energy)
    model.prepare_response(coupling)
    answer()
    return float(energy), model.get_response(coupling)


def _item(response):
    found = [item for item in response if isinstance(item, ParameterAdjointResponse)]
    assert len(found) <= 1
    return found[0] if found else None


@pytest.mark.parametrize("model", _MODELS)
def test_parameter_adjoint_matches_central_differences(model) -> None:
    """Exactly the model's parameters are present, each dE/dp equal to a central difference."""
    factory, names = _MODELS[model]
    on = PCMParameters(parameter_derivatives=True)
    _, response = _run(factory(_VALUES, on))
    item = _item(response)
    assert item is not None and item.name == "parameter_adjoint"
    assert set(item.as_dict()) == set(names)
    for name in names:
        key = name[len("dg_d"):]
        h = _STEP * abs(_VALUES[key])
        energies = [_run(factory({**_VALUES, key: _VALUES[key] + s * h}, PCMParameters()))[0]
                    for s in (+1, -1)]
        fd = (energies[0] - energies[1]) / (2 * h)
        value = getattr(item, name)
        assert isinstance(value, float)
        assert value == approx(fd, rel=_RTOL, abs=1e-12), (model, name, value, fd)


@pytest.mark.parametrize("model", _MODELS)
def test_parameter_adjoint_is_absent_by_default(model) -> None:
    factory, _ = _MODELS[model]
    _, response = _run(factory(_VALUES, PCMParameters()))
    assert _item(response) is None
    assert PCMParameters().parameter_derivatives is False


def test_parameter_adjoint_fields_absent_are_none() -> None:
    item = ParameterAdjointResponse(dg_deps=np.array(2.0))
    assert item.dg_deps == 2.0 and item.dg_dkappa is None
    assert item.as_dict() == {"dg_deps": 2.0}
    with raises(AttributeError):
        item.dg_deps = 1.0
