"""Shared configuration, native defaults and input ownership contracts."""

from dataclasses import FrozenInstanceError, asdict, fields, replace
import gc
import pickle
import weakref

import numpy as np
import pytest

import moist
from moist import library

#: Run context shared by cavities and models that do not test contexts
CONTEXT = moist.Context()


@pytest.mark.parametrize("kind", [
    moist.DROPParameters, moist.ISwiGParameters, moist.SvdWParameters,
    moist.CFCParameters, moist.IsodensityParameters, moist.PCMParameters, moist.GOSTSHYPParameters,
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
    lambda: moist.CavityDROP(parameters=moist.ISwiGParameters(), context=CONTEXT),
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


def _runtime_threads():
    """Thread count a fresh unset context takes from the OpenMP runtime."""
    return moist.Context(nthreads=0).nthreads


@pytest.fixture
def baseline():
    """This thread's OpenMP setting."""
    return _runtime_threads()


def test_context_thread_count_is_fixed(baseline):
    for nthreads in (2, 1, 0):
        context = moist.Context(nthreads=nthreads)
        expected = nthreads if nthreads > 0 else _runtime_threads()
        assert context.nthreads == expected
        assert _runtime_threads() == baseline
    logged = moist.Context(nthreads=0, verbosity=2, debug=True)
    assert logged.nthreads == baseline
    assert moist.Context(nthreads=0, verbosity=1, debug=False).nthreads == baseline
    assert moist.Context().nthreads == baseline
    with pytest.raises(TypeError):
        moist.Context(nthreads=2.5)
    with pytest.raises(TypeError):
        moist.Context(nthreads=True)
    with pytest.raises(TypeError):
        moist.Context(nthreads=0, debug=1)
    with pytest.raises(RuntimeError, match="Thread count must not be negative"):
        moist.Context(nthreads=-1)
    # Threads and logging belong to the context alone
    for parameters in (moist.DROPParameters, moist.ISwiGParameters):
        for name in ("nthreads", "verbosity", "debug"):
            with pytest.raises(TypeError):
                parameters(**{name: 1})


def test_context_leaves_the_runtime_alone(baseline):
    first = moist.Context(nthreads=baseline + 1)
    second = moist.Context(nthreads=baseline + 2)
    assert _runtime_threads() == baseline
    assert first.nthreads == baseline + 1
    assert second.nthreads == baseline + 2


@pytest.mark.parametrize("config", [
    moist.DROP(lsf=moist.SvdW(), parameters=moist.DROPParameters(nleb=26)),
    moist.ISwiG(parameters=moist.ISwiGParameters(nleb=26)),
])
def test_shared_context_lifetime(config, baseline):
    context = moist.Context(nthreads=2)
    owner = weakref.ref(context)
    cavity = config.build(context=context)
    model = moist.SolvationModel(context, cavity, [moist.ModelComponentPV(1e-4)])
    assert cavity.context is context
    assert model.context is context
    assert model.cavity.context is context
    del context
    gc.collect()
    assert owner() is not None
    model.update(moist.Structure([1], [[0., 0., 0.]]))
    assert model.cavity.snapshot().ngrid > 0
    assert model.context.nthreads == 2
    del model
    gc.collect()
    del cavity
    gc.collect()
    assert owner() is None


def test_parts_keep_their_own_context(baseline):
    context = moist.Context(nthreads=2)
    cavity = moist.CavityISwiG(context=context)
    other = moist.Context(nthreads=1)
    component = moist.ModelComponentPV(1e-4, context=context)
    model = moist.SolvationModel(other, cavity, [component])
    assert cavity.context is context and component.context is context
    assert model.context is other
    assert model.cavity.context is context
    del other
    gc.collect()
    model.update(moist.Structure([1], [[0., 0., 0.]]))
    assert model.context.nthreads == 1
    assert context.nthreads == 2
    del model, cavity, component, context
    gc.collect()


@pytest.mark.parametrize("config", [
    moist.DROP(lsf=moist.SvdW(), parameters=moist.DROPParameters(nleb=26)),
    moist.ISwiG(parameters=moist.ISwiGParameters(nleb=26)),
])
def test_parts_without_context_run_on_the_model(config):
    cavity = config.build()
    component = moist.ModelComponentPV(1e-4)
    assert cavity.context is None and component.context is None
    structure = moist.Structure([1], [[0., 0., 0.]])
    with pytest.raises(RuntimeError, match="Cavity has no context"):
        cavity.update(structure)
    model = moist.SolvationModel(CONTEXT, cavity, [component])
    assert model.cavity.context is CONTEXT
    model.update(structure)
    assert model.cavity.snapshot().ngrid > 0


@pytest.mark.parametrize("make", [
    lambda: moist.CavityDROP(context=object()),
    lambda: moist.CavityISwiG(context=object()),
    lambda: moist.ModelComponentCPCM(78.4, context=object()),
    lambda: moist.ModelComponentPV(1e-4, context=object()),
])
def test_context_type_is_checked(make):
    with pytest.raises(TypeError, match="context must be a Context"):
        make()


@pytest.mark.parametrize("make", [
    lambda: moist.SolvationModel(moist.CavityISwiG(context=CONTEXT),
                                 [moist.ModelComponentPV(1e-4)]),
    lambda: moist.SolvationModel(None, moist.CavityISwiG(context=CONTEXT), [moist.ModelComponentPV(1e-4)]),
])
def test_model_requires_context(make):
    with pytest.raises(TypeError, match="context|required positional"):
        make()


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
    cavity = config.build(source=source, context=CONTEXT)
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
        invalid.build(source=source, context=CONTEXT)


@pytest.mark.parametrize("factory,params", [
    (moist.CavityDROP, moist.DROPParameters(nleb=26)),
    (moist.CavityISwiG, moist.ISwiGParameters(nleb=26)),
])
def test_custom_radii_are_copied_and_control_native_surface(factory, params):
    values = np.array([2.5])
    radii = moist.CustomRadii(values)
    values[:] = 7
    by_element = moist.CustomRadii([2.5], numbers=[1])
    structure = moist.Structure([1], [[0., 0., 0.]])
    cavities = [factory(parameters=params,
                        radii=item, context=CONTEXT) for item in (radii, by_element)]
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
    cavity = moist.CavityISwiG(parameters=moist.ISwiGParameters(nleb=26),
                               radii=radii, context=CONTEXT)
    cavity.update(moist.Structure([1], [[0., 0., 0.]]))
    assert np.isfinite(cavity.area) and cavity.area > 0


def test_model_copy_retains_configuration_and_callback(gaussian_density):
    cavity = moist.CavityDROP(lsf=moist.Isodensity(), source=gaussian_density,
                              parameters=moist.DROPParameters(nleb=26), context=CONTEXT)
    model = moist.SolvationModel(CONTEXT, cavity, [moist.ModelComponentPV(1e-4)])
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
    assert model.context is CONTEXT


def test_low_level_model_retains_callback_owner(gaussian_density):
    cavity = library.new_drop_cavity_isodensity_callback(
        gaussian_density, moist.DROPParameters(nleb=26), moist.IsodensityParameters(),
        context=CONTEXT._as_handle(),
    )
    owner = weakref.ref(cavity)
    component = library.new_pv_component(1e-4)
    model = library.new_general_model(CONTEXT._as_handle(), cavity, [component])
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
    first, second = config.build(context=CONTEXT), config.build(context=CONTEXT)
    first.update(moist.Structure([1], [[0., 0., 0.]]))
    with pytest.raises(RuntimeError, match="not been successfully updated"):
        second.snapshot()
    assert pickle.loads(pickle.dumps(config)) == config


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
    cavity = moist.CavityISwiG(parameters=moist.ISwiGParameters(nleb=26), context=CONTEXT)
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


@pytest.mark.parametrize("factory,parameters", [
    (moist.SvdW, moist.SvdWParameters(blend_k=7.0)),
    (moist.CFC, moist.CFCParameters(a1=0.7)),
    (moist.Isodensity, moist.IsodensityParameters(rho_iso=0.003)),
    (moist.ISwiG, moist.ISwiGParameters(nleb=26)),
])
def test_surface_configuration_preserves_explicit_parameters(factory, parameters):
    config = factory(parameters=parameters)
    assert config.parameters == parameters
    assert pickle.loads(pickle.dumps(config)) == config
    with pytest.raises(TypeError, match="parameters must"):
        factory(parameters=moist.PCMParameters())


@pytest.mark.parametrize("surface", [moist.SvdW(), moist.CFC()])
@pytest.mark.parametrize("kwargs", [{"source": object()}, {"pass_order": True}])
def test_geometric_configuration_rejects_density_inputs(surface, kwargs):
    with pytest.raises(TypeError, match="geometric LSF"):
        moist.DROP(lsf=surface).build(context=CONTEXT, **kwargs)


def test_configuration_rejects_wrong_lsf_and_radius_models():
    with pytest.raises(TypeError, match="LSF"):
        moist.DROP(lsf=object())
    for factory in (moist.ISwiG, lambda **kw: moist.DROP(lsf=moist.SvdW(), **kw)):
        with pytest.raises(TypeError, match="radius model"):
            factory(radii=object())
        with pytest.raises(TypeError, match="radius model"):
            factory(radii=moist.Radii())
    with pytest.raises(TypeError, match="does not accept"):
        moist.ISwiG().build(source=object(), context=CONTEXT)


@pytest.mark.parametrize("factory,parameters", [
    (moist.ISwiG, moist.ISwiGParameters(nleb=26)),
    (lambda **kw: moist.DROP(lsf=moist.SvdW(), **kw), moist.DROPParameters(nleb=26)),
])
def test_configuration_build_preserves_radii_and_parameters(factory, parameters):
    radii = moist.CustomRadii([2.5])
    config = factory(parameters=parameters, radii=radii)
    cavity = config.build(context=CONTEXT)
    assert cavity.parameters == parameters
    assert cavity.radius_model is radii
    cavity.update(moist.Structure([1], [[0., 0., 0.]]))
    np.testing.assert_array_equal(cavity.radii, [2.5])


def test_isodensity_configuration_callback_object_and_order(monkeypatch, gaussian_density):
    class Source:
        density = staticmethod(gaussian_density)

    source = Source()
    parameters = moist.DROPParameters(nleb=26)
    config = moist.DROP(lsf=moist.Isodensity(), parameters=parameters)
    observed = {}
    original = library.new_drop_cavity_isodensity_callback

    def capture(callback, drop, lsf, radii, *, pass_order, context):
        observed.update(callback=callback, pass_order=pass_order)
        return original(callback, drop, lsf, radii, pass_order=pass_order, context=context)

    monkeypatch.setattr(library, "new_drop_cavity_isodensity_callback", capture)
    cavity = config.build(source=source, pass_order=True, context=CONTEXT)
    assert observed == {"callback": gaussian_density, "pass_order": True}
    cavity.update(moist.Structure([1], [[0., 0., 0.]]))
    assert cavity.snapshot().ngrid > 0
    with pytest.raises(TypeError, match="requires a callable"):
        config.build(source=object(), context=CONTEXT)


def test_internal_isodensity_rejects_callback_order():
    basis = moist.GaussianBasis(shell_atom=[0], shell_l=[0], shell_nprim=[1],
                                exponents=[1.], coefficients=[1.])
    source = moist.InternalDensity(basis, [[1.]])
    with pytest.raises(TypeError, match="only to callbacks"):
        moist.DROP(lsf=moist.Isodensity()).build(source=source, pass_order=True, context=CONTEXT)


@pytest.mark.parametrize("values,kwargs", [
    ([], {}), ([[2.5]], {}), ([float("nan")], {}), ([0.], {}), ([-1.], {}),
    ([2.5], {"numbers": [1, 8]}), ([2.5], {"numbers": [1.5]}),
    ([2.5], {"numbers": [0]}), ([2.5], {"numbers": [119]}),
    ([2.5, 3.], {"numbers": [1, 1]}),
])
def test_custom_radii_reject_invalid_data(values, kwargs):
    with pytest.raises(ValueError):
        moist.CustomRadii(values, **kwargs)


def test_custom_radii_element_mapping_handles_reordered_atoms():
    radii = moist.CustomRadii([2.5, 3.], numbers=[1, 8])
    assert radii.numbers == (1, 8)
    config = moist.ISwiG(parameters=moist.ISwiGParameters(nleb=26), radii=radii)
    cavity = config.build(context=CONTEXT)
    cavity.update(moist.Structure([8, 1], [[0., 0., 0.], [7., 0., 0.]]))
    np.testing.assert_array_equal(cavity.radii, [3., 2.5])


@pytest.mark.parametrize("radii,kind", [
    (moist.CPCMRadii(), "cpcm"), (moist.SMDRadii(), "smd"),
    (moist.D3Radii(), "d3"), (moist.COSMORadii(), "cosmo"),
    (moist.BondiRadii(), "bondi"),
])
def test_builtin_radii_select_their_native_radius_set(radii, kind):
    parameters = moist.ISwiGParameters(nleb=26)
    structure = moist.Structure([1, 3, 8], [[0., 0., 0.], [7., 0., 0.], [14., 0., 0.]])
    actual = moist.CavityISwiG(parameters=parameters, radii=radii, context=CONTEXT)
    actual.update(structure)
    expected = library.new_iswig_cavity(parameters,
                                        library.new_radii(kind), context=CONTEXT._as_handle())
    library.update_cavity(expected, structure._as_handle())
    snapshot = library.get_cavity_results(expected)
    np.testing.assert_array_equal(actual.radii, snapshot["radii"])


@pytest.mark.parametrize("config", [
    moist.SvdW(), moist.CFC(), moist.Isodensity(),
    moist.DROP(lsf=moist.SvdW()), moist.ISwiG(),
])
def test_surface_configuration_is_immutable(config):
    with pytest.raises(FrozenInstanceError):
        config.parameters = config.parameters


def test_gostshyp_settings_and_validation():
    parameters = moist.GOSTSHYPParameters(regularization_start=2e-9, regularization_end=2e-8,
                                          suppress_negative_amplitudes=True)
    component = moist.ModelComponentGOSTSHYP(1e-3, parameters=parameters)
    assert component.parameters == parameters
    assert not moist.GOSTSHYPParameters().suppress_negative_amplitudes
    with pytest.raises(RuntimeError, match="start must be nonnegative"):
        moist.ModelComponentGOSTSHYP(1e-3, parameters=replace(parameters, regularization_start=-1e-9))
    for end in (2e-9, 1e-9):
        with pytest.raises(RuntimeError, match="end must exceed start"):
            moist.ModelComponentGOSTSHYP(1e-3, parameters=replace(parameters, regularization_end=end))
    for name in ("regularization_start", "regularization_end"):
        with pytest.raises(ValueError, match="finite"):
            moist.GOSTSHYPParameters(**{name: float("nan")})
    with pytest.raises(TypeError, match="GOSTSHYPParameters"):
        moist.ModelComponentGOSTSHYP(1e-3, parameters=moist.PCMParameters())
    with pytest.raises(RuntimeError, match="pressure and scale must be finite"):
        moist.ModelComponentGOSTSHYP(float("inf"))


def test_gostshyp_registered_parameter_printout():
    model = moist.SolvationModel(CONTEXT, moist.CavityISwiG(),
                                [moist.ModelComponentGOSTSHYP(1e-3)])
    lines = [line.split() for line in model.parameters_text().splitlines()]
    rows = {name: [line for line in lines if line[:1] == [name]]
            for name in ("regularization_start", "regularization_end", "suppress_negative_amplitudes")}
    assert all(len(found) == 1 for found in rows.values())
    assert float(rows["regularization_start"][0][-1]) == 1e-12
    assert float(rows["regularization_end"][0][-1]) == 1e-10
    assert rows["suppress_negative_amplitudes"][0][-1] == "F"
