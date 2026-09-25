"""Cross-language contracts exercised through public Python and C entry points."""

import numpy as np
import pytest
import moist
from moist import library


def _model(components=None):
    model = moist.SolvationModel(moist.CavityISwiG(parameters=moist.ISwiGParameters(nleb=26)),
                                components or [moist.ModelComponentPV(1e-4)])
    structure = moist.Structure([1, 1], [[0., 0., 0.], [2., 1., .2]])
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
        moist.Structure(numbers, [[0., 0., 0.]])


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
    cavity = config.build(source=source)
    offsets, powers = cavity.isodensity_layout()
    np.testing.assert_array_equal(offsets, [0, 1])
    np.testing.assert_array_equal(powers, [[0, 0, 0]])
    structure = moist.Structure([1], [[0., 0., 0.]])
    reference = config.build(source=gaussian_density)
    cavity.update(structure)
    reference.update(structure)
    np.testing.assert_allclose(cavity.xyz, reference.xyz, atol=1e-10)
    model = moist.SolvationModel(cavity, [moist.ModelComponentPV(1e-4)])
    model.update(structure)
    initial = _energy(model)
    weights = (np.zeros(cavity.ngrid), np.zeros(cavity.ngrid), np.ones((cavity.ngrid, 3)))
    standalone_weights = cavity.contract_surface_lsf_weights(*weights)
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
    cavity = moist.DROP(lsf=moist.Isodensity(), parameters=moist.DROPParameters(nleb=26)).build(source=source)
    offsets, powers = cavity.isodensity_layout()
    np.testing.assert_array_equal(offsets, [0, 1, 4])
    assert {tuple(row) for row in powers} == {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)}


def test_diagnostic_contractions_match_explicit_anchor_tensors():
    cavity = moist.CavityDROP(parameters=moist.DROPParameters(nleb=26))
    structure = moist.Structure([1, 1], [[0., 0., 0.], [2., 1., .2]])
    cavity.update(structure)
    cavity.assemble_amat()
    cavity.compute_cavity_gradient()
    anchor = cavity.get_anchor_gradient()
    q = np.linspace(.1, .2, cavity.ngrid) * cavity.get("f")
    gradient = cavity.contract_amat_nuclear_gradient(q, q)
    # A finite-difference contraction independently pins nuclear-axis ordering.
    h = 1e-5
    energies = []
    for shift in (-h, h):
        xyz = structure.positions
        xyz[1, 2] += shift
        moved = moist.CavityDROP(parameters=moist.DROPParameters(nleb=26))
        moved.update(moist.Structure(structure.numbers, xyz))
        matrix, _ = moved.assemble_amat()
        energies.append(np.einsum("i,ij,j->", q, matrix, q))
    assert gradient[1, 2] == pytest.approx((energies[1] - energies[0])/(2*h), abs=2e-7)
    weights = np.arange(cavity.ngrid * 3).reshape(-1, 3) * .001
    actual = cavity.contract_pcm_nuclear_gradient(np.zeros(cavity.ngrid), weights,
                                                  np.ones(2))
    expected = np.einsum("iAkj,ij->Ak", anchor.xyz1_rA, weights)
    np.testing.assert_allclose(actual, expected, atol=1e-12)
