"""Cross-language contracts exercised through public Python and C entry points."""

import gc
import math
import numpy as np
import pytest
import moist
from moist import library

#: Run context shared by every cavity and model in this module
CONTEXT = moist.Context()


FD_ABS_TOL = 1e-10
FD_REL_TOL = 1e-9
FD_POSITION_STEP = 5e-4  # bohr, for nuclear and surface point positions
FD_SURFACE_RELATIVE_STEP = 5e-4  # fraction of xi (bohr^-1) or f (dimensionless)
FD_PROJECTION_TOL = 1e-14  # DROP level-set residual for displaced nuclear references


def _assert_fd_close(actual, reference):
    assert np.all(np.isfinite(actual))
    assert np.all(np.isfinite(reference))
    bound = np.maximum(FD_ABS_TOL, FD_REL_TOL * np.abs(reference))
    normalized_error = np.abs(actual - reference) / bound
    assert np.all(normalized_error <= 1), f"max FD error / bound: {np.max(normalized_error)}"


def _model(components=None):
    model = moist.SolvationModel(CONTEXT, moist.CavityISwiG(parameters=moist.ISwiGParameters(nleb=26),
                                                   context=CONTEXT),
                                components or [moist.ModelComponentPV(1e-4)])
    structure = moist.Structure([1, 1], [[0., 0., 0.], [2., 1., .2]], z_eff=[1, 1])
    model.update(structure)
    return model, structure


def _energy(model):
    """Energy of a model whose components ask the host for nothing."""
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    assert list(coupling) == []
    energy = np.array(0.)
    model.get_energy(coupling, energy)
    return float(energy)


def _gradient(model, natoms):
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    assert list(coupling) == []
    gradient = np.zeros((natoms, 3))
    model.get_gradient(coupling, gradient)
    return gradient


def test_model_accumulators():
    model, structure = _model()
    coupling = model.new_coupling()
    energy = np.array(7.)
    with pytest.raises(RuntimeError, match="staged"):
        model.get_energy(coupling, energy)
    assert energy == 7
    model.prepare_energy(coupling)
    expected = _energy(model)
    model.get_energy(coupling, energy)
    model.get_energy(coupling, energy)
    assert energy == pytest.approx(7 + 2 * expected)
    gradient = np.full((len(structure), 3), 3.)
    with pytest.raises(RuntimeError, match="staged"):
        model.get_gradient(coupling, gradient)
    np.testing.assert_array_equal(gradient, 3.)
    model.prepare_gradient(coupling)
    expected = _gradient(model, len(structure))
    assert np.linalg.norm(expected) > 1e-9
    model.get_gradient(coupling, gradient)
    model.get_gradient(coupling, gradient)
    np.testing.assert_allclose(gradient, 3 + 2 * expected)
    with pytest.raises(TypeError, match="float64"):
        model.get_energy(coupling, 0.)
    with pytest.raises(ValueError, match="C-contiguous"):
        model.get_gradient(coupling, np.zeros((2, 6))[:, ::2])
    with pytest.raises(TypeError, match="float64"):
        model.get_gradient(coupling, np.zeros((2, 3), dtype=np.float32))


def test_cyclic_garbage_deletes_coupling_before_model(monkeypatch):
    """Deleting a coupling reaches into its model, whatever order GC finds them in."""
    order = []
    for kind in (library.ModelHandle, library.CouplingHandle):
        delete = kind._delete
        monkeypatch.setattr(kind, "_delete", staticmethod(
            lambda handle, kind=kind, delete=delete: (order.append(kind), delete(handle))))
    gc.disable()
    try:
        for _ in range(8):
            model, _ = _model()
            coupling = model.new_coupling()
            box = {"model": model, "coupling": coupling}
            box["box"] = box
            del model, coupling, box
        gc.collect()
    finally:
        gc.enable()
    assert order.count(library.CouplingHandle) == 8
    assert order.index(library.ModelHandle) > max(
        i for i, kind in enumerate(order) if kind is library.CouplingHandle)
    assert not library._COUPLING_PARENTS


def test_independent_answers():
    """Each output is submitted on its own; a bad one stops at itself."""
    model, _ = _model([moist.ModelComponentCPCM(32)])
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    ngrid = model.cavity.ngrid
    visited = []
    for request in coupling:
        visited.append(request)
        assert request.missing == {"phi", "dphi_dr", "dphi_dxi"}
        with pytest.raises(ValueError):
            coupling.answer(phi=np.zeros(ngrid), dphi_dr=["bad"], dphi_dxi=np.zeros(ngrid))
    # phi was stored before the bad output stopped the call; a second pass
    # asks for the rest.
    for request in coupling:
        visited.append(request)
        assert request.missing == {"dphi_dr", "dphi_dxi"}
        coupling.answer(dphi_dr=np.zeros((ngrid, 3)), dphi_dxi=np.zeros(ngrid))
    assert len(visited) == 2
    assert list(coupling) == []
    gradient = np.zeros((2, 3))
    model.get_gradient(coupling, gradient)
    # Staging again ends the walk: answering now raises, and the earlier
    # snapshots stay what they were.
    model.prepare_energy(coupling)
    with pytest.raises(RuntimeError, match="No current coupling request"):
        coupling.answer(phi=np.zeros(ngrid))
    assert visited[1].missing == {"dphi_dr", "dphi_dxi"}


@pytest.mark.parametrize("numbers", [[1.2], [1.], [True], [2**32+1], [0], [119]])
def test_atomic_numbers_not_coerced(numbers):
    with pytest.raises(ValueError):
        moist.Structure(numbers, [[0., 0., 0.]], z_eff=numbers)


def test_response_without_pcm_has_no_potential_adjoint():
    model, _ = _model()
    coupling = model.new_coupling()
    model.prepare_response(coupling)
    response = model.get_response(coupling)
    assert list(response) == []


def test_native_strings():
    assert moist.get_version_string().startswith(library.get_api_version())
    for style in ("full", "short", "ascii", "build"):
        assert moist.get_banner(style)
    with pytest.raises(ValueError):
        moist.get_banner("invalid")


def test_internal_density_matches_callback_and_updates_model(gaussian_density):
    basis = moist.GaussianBasis(shell_atom=[0], shell_l=[0], shell_nprim=[1],
                               exponents=[.5], coefficients=[1.])
    density = np.ones((1, 1))
    source = moist.InternalDensity(basis, density)
    density[:] = 9  # The source owns its input.
    config = moist.DROP(lsf=moist.Isodensity(), parameters=moist.DROPParameters(nleb=26))
    cavity = config.build(source=source, context=CONTEXT)
    offsets, powers = cavity.isodensity_layout()
    np.testing.assert_array_equal(offsets, [0, 1])
    np.testing.assert_array_equal(powers, [[0, 0, 0]])
    structure = moist.Structure([1], [[0., 0., 0.]], z_eff=[1])
    reference = config.build(source=gaussian_density, context=CONTEXT)
    cavity.update(structure)
    reference.update(structure)
    np.testing.assert_allclose(cavity.xyz, reference.xyz, atol=1e-10)
    model = moist.SolvationModel(CONTEXT, cavity, [moist.ModelComponentPV(1e-4)])
    model.update(structure)
    initial = _energy(model)
    weights = (np.zeros(cavity.ngrid), np.zeros(cavity.ngrid), np.ones((cavity.ngrid, 3)))
    standalone_weights = cavity.contract_surface_lsf_weights(*weights)
    # Read the C boundary independently to pin the returned Python channels.
    raw_weights = (np.zeros(cavity.ngrid), np.zeros((cavity.ngrid, 3)),
                   np.zeros((cavity.ngrid, 3, 3)))
    library.error_check(library.lib.moist_contract_surface_lsf_weights_extended)(
        cavity._as_handle().handle,
        *(library._cast("double*", value) for value in (*weights, *raw_weights)),
        library.ffi.NULL, library.ffi.NULL, library.ffi.NULL)
    assert np.linalg.norm(raw_weights[0]) > 1e-6
    assert np.linalg.norm(raw_weights[1]) > 1e-6
    for actual, expected in zip(standalone_weights, raw_weights):
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)
    source.density_matrix = [[2.]]
    model.update(structure)
    assert _energy(model) > initial
    for before, after in zip(standalone_weights, cavity.contract_surface_lsf_weights(*weights)):
        np.testing.assert_array_equal(before, after)
    # The model updates its own copied cavity; no density aliases escape.
    np.testing.assert_array_equal(source.density_matrix, [[2.]])
    matrix = source.density_matrix
    matrix[:] = 10
    np.testing.assert_array_equal(source.density_matrix, [[2.]])
    with pytest.raises(ValueError, match="shape"):
        source.density_matrix = np.eye(2)


def test_c_field_shape_uses_row_major_order():
    model, _ = _model()
    fields = {field.name: field for field in model.cavity.fields()}
    assert fields["xyz"].shape == (model.cavity.ngrid, 3)
    xyz = model.cavity.get("xyz")
    assert xyz.flags.c_contiguous
    np.testing.assert_array_equal(xyz, model.cavity.xyz)


def test_internal_basis_layout():
    basis = moist.GaussianBasis(shell_atom=[0, 0], shell_l=[0, 1], shell_nprim=[1, 1],
                               exponents=[.5, .5], coefficients=[1., 1.])
    source = moist.InternalDensity(basis, np.diag([1., .1, .2, .3]))
    config = moist.DROP(lsf=moist.Isodensity(), parameters=moist.DROPParameters(nleb=26))
    cavity = config.build(source=source, context=CONTEXT)
    offsets, powers = cavity.isodensity_layout()
    np.testing.assert_array_equal(offsets, [0, 1, 4])
    assert {tuple(row) for row in powers} == {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)}


def test_diagnostic_contractions_match_explicit_anchor_tensors():
    cavity = moist.CavityDROP(parameters=moist.DROPParameters(nleb=26), context=CONTEXT)
    structure = moist.Structure([1, 1], [[0., 0., 0.], [2., 1., .2]], z_eff=[1, 1])
    cavity.update(structure)
    cavity.compute_cavity_gradient()
    anchor = cavity.get_anchor_gradient()
    weights = np.arange(cavity.ngrid * 3).reshape(-1, 3) * .001
    actual = cavity.contract_pcm_nuclear_gradient(np.zeros(cavity.ngrid), weights,
                                                  np.ones(2))
    expected = np.einsum("iAkj,ij->Ak", anchor.xyz1_rA, weights)
    np.testing.assert_allclose(actual, expected, atol=1e-12)


@pytest.mark.parametrize("kind", ["energy", "gradient"])
@pytest.mark.parametrize("invalid", ["shape", "writeability", "alignment"])
def test_accumulator_buffer_guards(kind, invalid, monkeypatch):
    model, _ = _model()
    coupling = model.new_coupling()
    getattr(model, f"prepare_{kind}")(coupling)
    shape = () if kind == "energy" else (2, 3)
    if invalid == "shape":
        value = np.zeros((1,) if kind == "energy" else (2, 2))
        match = "shape"
    elif invalid == "writeability":
        value = np.zeros(shape)
        value.flags.writeable = False
        match = "writable"
    else:
        storage = np.zeros((1 if not shape else np.prod(shape)) * 8 + 1, dtype=np.uint8)
        value = np.ndarray(shape, dtype=np.float64, buffer=storage, offset=1)
        assert not value.flags.aligned
        match = "aligned"

    def reject_pointer(*args):
        raise AssertionError("Invalid accumulator reached the native pointer cast")

    # Guard failures must happen before a native call can read an invalid buffer.
    monkeypatch.setattr(library, "_cast", reject_pointer)
    with pytest.raises(ValueError, match=match):
        getattr(model, f"get_{kind}")(coupling, value)


@pytest.mark.parametrize("invalid", ["nan", "infinity", "asymmetry"])
def test_internal_density_rejects_invalid_values(invalid):
    basis = moist.GaussianBasis(shell_atom=[0], shell_l=[1], shell_nprim=[1],
                               exponents=[.5], coefficients=[1.])
    source = moist.InternalDensity(basis, np.eye(3))
    value = np.eye(3)
    if invalid == "asymmetry":
        value[0, 1] = .5
        match = "symmetric"
    else:
        value[0, 0] = np.nan if invalid == "nan" else np.inf
        match = "finite"
    with pytest.raises(ValueError, match=match):
        source.density_matrix = value
    np.testing.assert_array_equal(source.density_matrix, np.eye(3))


@pytest.mark.parametrize("style, native_style", [
    ("full", "moist_banner_full"), ("short", "moist_banner_short"),
    ("ascii", "moist_banner_ascii"), ("build", "moist_banner_build"),
])
def test_banner_style_reaches_native_entry(style, native_style, monkeypatch):
    calls = []

    def text(entry, *args):
        calls.append((entry, args))
        return "banner text"

    monkeypatch.setattr(library, "_native_text", text)
    assert moist.get_banner(style) == "banner text"
    assert calls == [(library.lib.moist_get_banner,
                      (getattr(library.lib, native_style),))]
