"""Shared configuration, native defaults and input ownership contracts."""

from dataclasses import FrozenInstanceError, asdict, fields, replace
import gc
import pickle
import weakref

import numpy as np
import pytest

import moist
from . import library


@pytest.mark.parametrize("kind", [
    moist.DROPParameters, moist.ISwiGParameters, moist.SvdWParameters,
    moist.CFCParameters, moist.IsodensityParameters, moist.PCMParameters,
    moist.ModelParameters,
])
def test_parameters_cover_native_options_and_defaults(kind):
    parameters = kind()
    native = library._options(parameters._kind)
    members = {name for name, _ in library.ffi.typeof(native[0]).fields}
    members -= {"struct_size", "reserved0"}
    assert {item.name for item in fields(parameters)} == members
    assert asdict(parameters) == {name: getattr(native, name) for name in members}
    assert pickle.loads(pickle.dumps(parameters)) == parameters
    name = fields(parameters)[0].name
    with pytest.raises(FrozenInstanceError):
        setattr(parameters, name, getattr(parameters, name))
    with pytest.raises(TypeError):
        kind(1)


@pytest.mark.parametrize("make", [
    lambda: moist.DROPParameters(nlebb=26),
    lambda: moist.CavityDROP(parameters=moist.ISwiGParameters()),
    lambda: moist.CavityDROP(parameters=moist.DROPParameters(), nleb=26),
    lambda: moist.SvdW(parameters=moist.SvdWParameters(), blend_k=4),
    lambda: moist.ModelComponentCPCM(32, "lu", parameters=moist.PCMParameters()),
    lambda: moist.SolvationModel(moist.CavityDROP(), [moist.ModelComponentPV(1e-5)],
                                verbosity=0, parameters=moist.ModelParameters()),
])
def test_configuration_errors_are_explicit(make):
    with pytest.raises(TypeError):
        make()


def test_parameter_types_and_nonfinite_values():
    with pytest.raises(TypeError):
        moist.DROPParameters(nleb=26.5)
    with pytest.raises(ValueError, match="finite"):
        moist.DROPParameters(tolerance=float("nan"))
    assert moist.PCMParameters(solver="lu").solver is moist.PCMSolver.LU
    with pytest.raises(ValueError):
        moist.PCMParameters(solver=2.5)
    with pytest.raises(ValueError, match="finite"):
        moist.PCMParameters(solver_tol=float("nan"))
    with pytest.raises(TypeError):
        moist.PCMParameters(solver_maxiter=20.5)


@pytest.mark.parametrize("component", [moist.ModelComponentCPCM, moist.ModelComponentCOSMO])
def test_pcm_iterative_controls_round_trip(component):
    defaults = moist.PCMParameters()
    assert defaults.solver_tol > 0.0 and defaults.solver_maxiter >= 1
    parameters = moist.PCMParameters(
        solver="iterative", solver_tol=1e-8, solver_maxiter=200)
    term = component(32, parameters=parameters)
    assert term.parameters.solver is moist.PCMSolver.ITERATIVE
    assert term.parameters.solver_tol == 1e-8
    assert term.parameters.solver_maxiter == 200
    with pytest.raises(RuntimeError):
        component(32, parameters=moist.PCMParameters(solver_maxiter=0))


@pytest.mark.parametrize("lsf", [moist.SvdW(), moist.CFC(), moist.Isodensity()])
def test_all_drop_surfaces_share_parameters(lsf, gaussian_density):
    parameters = moist.DROPParameters(
        nleb=26, tolerance=2e-9, proj_maxiter=200, proj_level=2,
        branch_weight_s=.08, rho_grid_h=.8, wleb_prune_level=1,
    )
    config = moist.DROP(lsf=lsf, parameters=parameters)
    source = gaussian_density if config.density_dependent else None
    cavity = config.build(source=source)
    cavity.update(moist.Structure([1], [[0., 0., 0.]]))
    assert cavity.snapshot().ngrid > 0
    assert cavity.parameters == parameters
    tolerance = library.ffi.new("double *")
    library.error_check(library.lib.moist_get_drop_cavity_tolerance)(
        cavity._as_handle().handle, tolerance
    )
    assert tolerance[0] == parameters.tolerance
    # Invalid settings must reach the native validator for every surface.
    invalid = replace(config, parameters=replace(parameters, wleb_prune_level=7))
    with pytest.raises(RuntimeError, match="wleb_prune_level"):
        invalid.build(source=source)


@pytest.mark.parametrize("factory", [moist.CavityDROP, moist.CavityISwiG])
def test_custom_radii_are_copied_and_control_native_surface(factory):
    values = np.array([2.5])
    radii = moist.CustomRadii(values)
    values[:] = 7
    by_element = moist.CustomRadii([2.5], numbers=[1])
    structure = moist.Structure([1], [[0., 0., 0.]])
    cavities = [factory(nleb=26, radii=item) for item in (radii, by_element)]
    for cavity in cavities:
        cavity.update(structure)
        np.testing.assert_allclose(cavity.radii, [2.5])
    np.testing.assert_allclose(cavities[0].xyz, cavities[1].xyz)
    assert cavities[0].radius_model is radii
    assert pickle.loads(pickle.dumps(radii)) == radii


@pytest.mark.parametrize("radii", [
    moist.CPCMRadii(), moist.SMDRadii(), moist.D3Radii(),
    moist.COSMORadii(), moist.BondiRadii(),
])
def test_builtin_radius_models(radii):
    cavity = moist.CavityISwiG(nleb=26, radii=radii)
    cavity.update(moist.Structure([1], [[0., 0., 0.]]))
    assert np.isfinite(cavity.area) and cavity.area > 0


def test_model_copy_retains_configuration_and_callback(gaussian_density):
    cavity = moist.CavityDROP(lsf=moist.Isodensity(), source=gaussian_density,
                              parameters=moist.DROPParameters(nleb=26))
    model = moist.SolvationModel(cavity, [moist.ModelComponentPV(1e-4)],
                                 parameters=moist.ModelParameters(verbosity=0))
    del cavity
    model.update(moist.Structure([1], [[0., 0., 0.]]))
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    energy = np.array(0.)
    model.get_energy(coupling, energy)
    assert np.isfinite(energy)
    model.prepare_gradient(coupling)
    gradient = np.zeros((1, 3))
    model.get_gradient(coupling, gradient)
    assert np.isfinite(gradient).all()
    assert model.cavity.parameters.nleb == 26
    assert model.cavity.lsf == moist.Isodensity()
    assert model.parameters == moist.ModelParameters(verbosity=0)


def test_low_level_model_retains_callback_owner(gaussian_density):
    cavity = library.new_drop_cavity_isodensity_callback(
        gaussian_density, moist.DROPParameters(nleb=26), moist.IsodensityParameters()
    )
    owner = weakref.ref(cavity)
    component = library.new_pv_component(1e-4)
    model = library.new_general_model(cavity, [component], moist.ModelParameters())
    del cavity
    gc.collect()
    assert owner() is not None
    structure = moist.Structure([1], [[0., 0., 0.]])
    library.update_model(model, structure._as_handle())
    assert library.get_cavity_sizes(library.get_model_cavity(model))[0] > 0
    del model
    gc.collect()
    assert owner() is None


def test_parameter_replacement_creates_independent_live_state():
    config = moist.DROP(lsf=moist.SvdW(), parameters=moist.DROPParameters(nleb=26))
    finer = replace(config, parameters=replace(config.parameters, nleb=50))
    assert config.parameters.nleb == 26
    assert finer.parameters.nleb == 50
    first, second = config.build(), config.build()
    first.update(moist.Structure([1], [[0., 0., 0.]]))
    with pytest.raises(RuntimeError, match="not been successfully updated"):
        second.snapshot()
    assert pickle.loads(pickle.dumps(config)) == config


@pytest.mark.parametrize("component", [moist.ModelComponentCPCM, moist.ModelComponentCOSMO])
def test_pcm_parameter_construction_matches_legacy_solver(component):
    structure = moist.Structure([1], [[0., 0., 0.]])
    energies = []
    for term in (component(32, solver="lu"),
                 component(32, parameters=moist.PCMParameters(solver=moist.PCMSolver.LU))):
        model = moist.SolvationModel(moist.CavityISwiG(nleb=26), [term])
        model.update(structure)
        coupling = model.new_coupling()
        model.prepare_energy(coupling)
        for request in coupling:
            coupling.answer(phi=np.ones(model.cavity.ngrid))
        energy = np.array(0.)
        model.get_energy(coupling, energy)
        energies.append(float(energy))
        assert term.parameters.solver is moist.PCMSolver.LU
    assert energies[0] == pytest.approx(energies[1], abs=1e-14)


def test_structure_owns_input_buffers_on_construction_and_update():
    numbers = np.array([1], dtype=np.int32)
    positions = np.array([[0., 0., 0.]])
    lattice = np.eye(3) * 20
    periodic = np.zeros(3, dtype=np.bool_)
    structure = moist.Structure(numbers, positions, lattice, periodic)
    numbers[:] = 8
    positions[:] = 10
    lattice[:] = 40
    periodic[:] = True
    np.testing.assert_array_equal(structure.numbers, [1])
    np.testing.assert_array_equal(structure.positions, [[0., 0., 0.]])
    np.testing.assert_array_equal(structure.lattice, np.eye(3) * 20)
    np.testing.assert_array_equal(structure.periodic, [False, False, False])
    cavity = moist.CavityISwiG(nleb=26)
    cavity.update(structure)
    np.testing.assert_allclose(cavity.xyz.mean(axis=0), structure.positions[0], atol=1e-14)
    updated = np.array([[2., 0., 0.]])
    new_lattice = np.eye(3) * 30
    structure.update(updated, new_lattice)
    updated[:] = 100
    new_lattice[:] = 100
    cavity.update(structure)
    np.testing.assert_allclose(cavity.xyz.mean(axis=0), structure.positions[0], atol=1e-14)
    np.testing.assert_array_equal(structure.positions, [[2., 0., 0.]])
    np.testing.assert_array_equal(structure.lattice, np.eye(3) * 30)
