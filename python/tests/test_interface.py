from typing import Callable

import numpy as np
import pytest
from pytest import approx, raises

import moist
from moist import library
from moist import CFC, Isodensity, SvdW
from moist.interface import (
    Cavity,
    CavityDROP,
    CavityISwiG,
    CavitySnapshot,
    CavitySnapshotDROP,
    Coupling,
    GaussianMomentRequest,
    GaussianPotentialRequest,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    ModelComponentPV,
    DensityResponse,
    PCMSolver,
    PotentialAdjointResponse,
    Response,
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
            CavityDROP(
                lsf=SvdW(blend_k=5.5, blend_1b=1.0, blend_2b=0.0, blend_3b=3.0),
                nleb=26,
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
            CavityDROP(
                lsf=CFC(a1=-15.0, a2=-9.0, c=5.0, m=4),
                nleb=26,
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


@pytest.fixture
def branching_cross() -> Structure:
    """Five-carbon cross whose concave seams branch at ``proj_level=7``.

    The same fixture the Fortran suite uses (``get_test_cross``); the default
    projection level finds no second solution on it, which is why the branching
    tests below raise the level rather than the geometry.
    """

    aatoau = 1.8897261246257702
    return Structure(
        np.array([6, 6, 6, 6, 6]),
        np.array(
            [
                [0.00, 4.21, 0.00],
                [0.00, 0.00, 4.22],
                [0.00, -4.18, 0.00],
                [0.00, 0.00, -4.15],
                [0.02, 0.10, -0.20],
            ]
        )
        * aatoau,
    )


def test_cavity_declares_its_own_results(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """Every declared field describes itself well enough to be read blind."""
    cavity = CavityDROP()
    cavity.update(Structure(numbers, positions))

    fields = cavity.fields()
    assert fields

    results = cavity.results()
    assert set(results) == {field.name for field in fields}

    for field in fields:
        value = results[field.name]
        assert np.shape(value) == field.shape
        assert np.asarray(value).dtype == field.dtype
        assert cavity.describe(field.name)


def test_named_results_agree_with_the_snapshot(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """The typed snapshot is a view of the same declarations, not a second read."""
    cavity = CavityDROP()
    cavity.update(Structure(numbers, positions))

    snapshot = cavity.snapshot()
    results = cavity.results()

    assert results["ngrid"] == snapshot.ngrid
    assert results["nsph"] == snapshot.nsph
    assert results["area"] == approx(snapshot.area)
    assert np.array_equal(results["owner"], snapshot.owner)
    assert np.array_equal(results["xyz"], snapshot.xyz)
    assert results["xyz"].shape == (snapshot.ngrid, 3)
    assert np.array_equal(results["numbering"], snapshot.numbering)
    assert np.array_equal(results["branch_count"], snapshot.branch_count)


def test_uncomputed_results_are_absent_rather_than_zero(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """A property that was never requested must not read back as zeros."""
    structure = Structure(numbers, positions)

    plain = CavityDROP()
    plain.update(structure)
    assert "k1" not in {field.name for field in plain.fields()}
    with raises(KeyError, match="k1"):
        plain.get("k1")

    fine = CavityDROP(do_fine=True)
    fine.update(structure)
    curvature = fine.get("k1")
    assert curvature.shape == (fine.ngrid,)
    assert np.isfinite(curvature).all()

    with raises(KeyError, match="not_a_field"):
        fine.get("not_a_field")


def test_branching_is_read_from_the_cavity(branching_cross: Structure) -> None:
    """Branch data comes from moist's own arrays, not from unpacking an id."""
    cavity = CavityDROP(proj_level=7)
    cavity.update(branching_cross)
    snapshot = cavity.snapshot()

    assert snapshot.branch.min() == 1
    assert snapshot.branch.max() > 1, "fixture stopped producing branches"
    assert snapshot.branched.sum() > 0
    assert np.array_equal(snapshot.branched, snapshot.branch_count > 1)

    # numbering is the packing of the two, so the arrays and the id agree.
    base = snapshot.nsph * cavity.get("num_leb")
    assert np.array_equal(
        snapshot.numbering, snapshot.anchor_id + base * (snapshot.branch - 1)
    )

    # Every point in a branched group reports the same group size.
    for anchor in np.unique(snapshot.anchor_id[snapshot.branched]):
        group = snapshot.anchor_id == anchor
        assert snapshot.branch_count[group].min() == group.sum()
        assert np.array_equal(np.unique(snapshot.branch[group]), snapshot.branch[group])


def test_named_results_reach_every_cavity_type(
    numbers: np.ndarray, positions: np.ndarray
) -> None:
    """The field API is a cavity feature, not a DROP one."""
    cavity = CavityISwiG()
    cavity.update(Structure(numbers, positions))

    names = {field.name for field in cavity.fields()}
    assert {"xyz", "a", "owner", "numbering"} <= names
    assert np.array_equal(cavity.get("owner"), cavity.snapshot().owner)


def test_named_results_need_a_built_cavity() -> None:
    cavity = CavityDROP()
    with raises(RuntimeError, match="not been successfully updated"):
        cavity.fields()
    with raises(RuntimeError, match="not been successfully updated"):
        cavity.get("xyz")


def test_cavity_drop_is_the_svdw_surface() -> None:
    """Unqualified DROP retains its default surface while allowing composition."""
    assert CavityDROP().configuration == CavityDROP(lsf=SvdW()).configuration
    assert isinstance(moist.CavityDROP().lsf, moist.SvdW)


@pytest.mark.parametrize("lsf", [SvdW(), CFC()])
def test_drop_constructor_controls_are_validated_natively(lsf) -> None:
    with raises(RuntimeError, match="wleb_prune_level.*0-6"):
        CavityDROP(lsf=lsf, wleb_prune_level=7)


def test_cavity_specific_options_reach_the_native_implementations(
    numbers: np.ndarray,
    positions: np.ndarray,
) -> None:
    structure = Structure(numbers, positions)

    def surface(cavity: Cavity) -> CavitySnapshot:
        cavity.update(structure)
        return cavity.snapshot()

    svdw_default = surface(CavityDROP(lsf=SvdW(), nleb=26))
    svdw_custom = surface(CavityDROP(lsf=SvdW(blend_k=4.0), nleb=26))
    assert svdw_custom.area != approx(svdw_default.area, rel=1.0e-6)

    cfc_default = surface(CavityDROP(lsf=CFC(), nleb=26))
    cfc_custom = surface(CavityDROP(lsf=CFC(a1=-12.0, a2=-8.0, c=4.0, m=4), nleb=26))
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


# -----------------------------------------------------------------------------
# Host coupling protocol
# -----------------------------------------------------------------------------
#
# The tests below drive the same six steps a Fortran host does: create a
# coupling, prepare a phase, walk the requests with ``for request in coupling``,
# answer the missing outputs of each with ``coupling.answer``, read the result,
# read the response.  Grid inputs come from ``model.cavity``.


def _energy(model, coupling):
    energy = np.array(0.0)
    model.get_energy(coupling, energy)
    return float(energy)


def _gradient(model, coupling):
    gradient = np.zeros((model._natoms, 3))
    response = model.get_gradient(coupling, gradient)
    return gradient, response


def _answer_potential(coupling, phi):
    """Answer every potential request whose only missing output is ``phi``."""
    for request in coupling:
        assert request.missing <= {"phi"}, request.missing
        coupling.answer(phi=phi)


def _solve(model, phi):
    """Energy and potential adjoint of a potential-only model, the Fortran way."""
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    _answer_potential(coupling, phi)
    energy = _energy(model, coupling)
    model.prepare_response(coupling)
    _answer_potential(coupling, phi)
    response = model.get_response(coupling)
    return energy, _item(response, PotentialAdjointResponse).w_phi


def _item(response, kind):
    """The item of class ``kind`` the walk over ``response`` meets, or None."""
    found = [item for item in response if isinstance(item, kind)]
    assert len(found) <= 1, found
    return found[0] if found else None


def _answer_zeros(model, coupling):
    """Answer every missing output with zeros of the declared shape."""
    ngrid = model.cavity.ngrid
    for request in coupling:
        coupling.answer(**{name: np.zeros((ngrid, *request._outputs[name]))
                           for name in request.missing})


def _still_missing(coupling, request):
    """Live native requirement state of the current request, unlike its snapshot."""
    return {name for name in request._outputs
            if library.get_coupling_request_missing(coupling._handle, name)}


def test_general_model_iterates_cpcm_and_pv_components(diatomic) -> None:
    """A heterogeneous component list shares one authoritative live cavity."""

    structure = diatomic()
    pressure = 2.5e-4
    model = SolvationModel(
        CavityDROP(nleb=26),
        [ModelComponentCPCM(32.0), ModelComponentPV(pressure)],
    )
    model.update(structure)
    ngrid = model.cavity.ngrid

    energy, charges = _solve(model, np.zeros(ngrid))

    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    _answer_potential(coupling, np.zeros(ngrid))
    assert _energy(model, coupling) == approx(energy)
    model.prepare_response(coupling)
    # The energy answer is retained, so the response phase owes nothing on a
    # cavity whose geometry is independent of the host.
    assert list(coupling) == []
    response = model.get_response(coupling)

    assert energy == approx(pressure * model.cavity.volume, abs=2.0e-13)
    np.testing.assert_array_equal(charges, np.zeros(ngrid))
    assert isinstance(response, Response)
    # SvdW is a geometric level set, so the walk meets no density item:
    # nothing the host owns moves this surface.
    (adjoint,) = list(response)
    assert isinstance(adjoint, PotentialAdjointResponse)
    assert adjoint.name == "potential_adjoint"
    assert adjoint.w_phi.shape == (ngrid,)
    # A response is a plain value: iterating it again yields the same items.
    assert list(response) == [adjoint]


def test_cosmo_is_a_standalone_pcm_component() -> None:
    cpcm = ModelComponentCPCM(32.0, solver="lu")
    cosmo = ModelComponentCOSMO(32.0, solver="lu")

    assert cpcm.epsilon == cosmo.epsilon == 32.0
    assert cpcm.solver is cosmo.solver is PCMSolver.LU


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
        model = SolvationModel(CavityDROP(nleb=26), [component_type(epsilon)])
        model.update(structure)
        phi = np.linspace(-0.2, 0.3, model.cavity.ngrid)
        results[component_type] = _solve(model, phi)

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


def test_model_drives_the_three_phases(diatomic) -> None:
    """Energy, response and gradient through one coupling, with read-only results."""
    structure = diatomic()
    pressure = 2.5e-4
    model = SolvationModel(
        CavityDROP(nleb=26),
        [ModelComponentCPCM(32.0), ModelComponentPV(pressure)],
    )
    model.update(structure)
    coupling = model.new_coupling()
    ngrid = model.cavity.ngrid

    model.prepare_energy(coupling)
    _answer_zeros(model, coupling)
    energy = _energy(model, coupling)
    model.prepare_response(coupling)
    _answer_zeros(model, coupling)
    response = model.get_response(coupling)
    model.prepare_gradient(coupling)
    _answer_zeros(model, coupling)
    gradient, gradient_response = _gradient(model, coupling)

    assert isinstance(model.cavity, Cavity)
    assert isinstance(model.cavity.snapshot(), CavitySnapshotDROP)
    assert energy == approx(pressure * model.cavity.volume, abs=2.0e-13)
    # CPCM is present and was handed phi = 0, so zero charges is a genuine
    # result rather than an absent item.
    adjoint = _item(response, PotentialAdjointResponse)
    np.testing.assert_array_equal(adjoint.w_phi, np.zeros(ngrid))
    np.testing.assert_array_equal(
        _item(gradient_response, PotentialAdjointResponse).w_phi, np.zeros(ngrid))
    assert gradient.shape == (len(structure), 3)
    assert np.isfinite(gradient).all()
    assert not model.cavity.snapshot().xyz.flags.writeable
    assert not adjoint.w_phi.flags.writeable
    with raises(ValueError, match="WRITEABLE"):
        adjoint.w_phi.setflags(write=True)
    with raises(ValueError, match="WRITEABLE"):
        model.cavity.snapshot().xyz.setflags(write=True)


def test_general_model_names_the_request_no_host_answered(diatomic) -> None:
    """What a model needs is its own declaration, not a capability list.

    The model states its requests, the host answers what it can, and the
    required output left unanswered is reported by name.
    """
    structure = diatomic()
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    coupling = model.new_coupling()
    model.prepare_energy(coupling)

    with raises(RuntimeError, match="gaussian_potential.*missing required outputs: phi"):
        _energy(model, coupling)


def test_staging_a_phase_drops_the_previous_answers(diatomic) -> None:
    """Staging clears the answers; nothing survives into the next phase or update."""

    structure = diatomic()
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)

    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    _answer_potential(coupling, np.zeros(model.cavity.ngrid))
    assert _energy(model, coupling) == approx(0.0)

    previous_ngrid = model.cavity.ngrid
    model.update(structure)
    assert model.cavity.ngrid == previous_ngrid

    # The coupling survives the update -- it belongs to the model -- but every
    # answer in it does not, so the required request is owed again by name.
    model.prepare_energy(coupling)
    assert [request.name for request in coupling] == ["gaussian_potential"]
    with raises(RuntimeError, match="potential"):
        _energy(model, coupling)


def test_staging_ends_the_current_request(diatomic) -> None:
    """A request snapshot held across a staging cannot answer into it.

    Answers go to the request the cursor stands at, and every ``prepare_*``
    restarts the walk.  An answer given after re-staging -- say for a
    snapshot kept from the previous phase -- is refused instead of landing in
    the new staging without a word.  The snapshot itself stays readable.
    """

    structure = diatomic()
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    ngrid = model.cavity.ngrid

    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    for stale in coupling:
        break
    assert isinstance(stale, GaussianPotentialRequest)
    assert stale.missing == {"phi"}

    model.prepare_energy(coupling)

    with raises(RuntimeError, match="No current coupling request"):
        coupling.answer(phi=np.zeros(ngrid))
    # The snapshot is a value, untouched by the staging; repr stays a
    # debugging aid rather than a second failure.
    assert stale.missing == {"phi"}
    assert repr(stale) == "<GaussianPotentialRequest missing=['phi']>"

    # A fresh pass reaches the request again.
    for fresh in coupling:
        coupling.answer(phi=np.zeros(ngrid))
    assert fresh is not stale
    assert _energy(model, coupling) == approx(0.0)


def test_gaussian_width_is_scoped_to_its_request(diatomic) -> None:
    """The width is the input of one request kind, never a cavity field."""

    structure = diatomic()
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    coupling = model.new_coupling()

    assert not hasattr(coupling, "width")
    assert hasattr(GaussianMomentRequest, "width")
    assert not hasattr(GaussianPotentialRequest, "width")
    model.prepare_energy(coupling)
    visited = []
    for request in coupling:
        visited.append(request.name)
        assert not hasattr(request, "width")
        with raises(RuntimeError, match="gaussian_potential has no input 'width'"):
            library.get_coupling_request_width(coupling._handle, np.empty(model.cavity.ngrid))
    assert visited == ["gaussian_potential"]


def test_request_outputs_are_named_and_ordered(diatomic) -> None:
    """Outputs are the contract; answer takes them by keyword only."""

    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    ngrid = model.cavity.ngrid

    visited = []
    for request in coupling:
        visited.append(request.name)
        assert request.name == "gaussian_potential"
        assert tuple(request._outputs) == ("phi", "dphi_dr", "dphi_dxi")
        assert request.missing == set(request._outputs)
        with raises(TypeError, match="has no output 'gt'"):
            coupling.answer(gt=np.zeros(ngrid))
        with raises(TypeError):
            coupling.answer(np.zeros(ngrid))
        # Neither refusal stored anything.
        assert _still_missing(coupling, request) == set(request._outputs)

        coupling.answer(phi=np.zeros(ngrid), dphi_dr=np.zeros((ngrid, 3)), dphi_dxi=np.zeros(ngrid))
        assert _still_missing(coupling, request) == set()
        # ``None`` skips an output rather than un-answering it.
        coupling.answer(phi=None)
        assert _still_missing(coupling, request) == set()
    assert visited == ["gaussian_potential"]
    # Nothing is left for another pass.
    assert list(coupling) == []


@pytest.mark.parametrize("cavity_type", [CavityISwiG, CavityDROP])
def test_grid_inputs_come_from_the_cavity(diatomic, cavity_type) -> None:
    """The coupling carries no grid; the model's cavity is the source of every input.

    ``xyz``, ``xi0``, ``normal0``, ``a`` and ``f`` are properties of every
    cavity type, as ``model%cavity`` holds them in Fortran, and they read the
    same arrays as the named-field catalogue.
    """

    model = SolvationModel(cavity_type(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    (request,) = list(coupling)
    assert isinstance(request, GaussianPotentialRequest)

    for name in ("grid", "grid_fields", "ngrid", "phase", "model", "requests", "pending"):
        assert not hasattr(coupling, name)
    for name in ("xyz", "xi", "ngrid", "set", "required", "handle", "n_missing", "answer"):
        assert not hasattr(request, name)

    cavity = model.cavity
    ngrid = cavity.ngrid
    shapes = {"xyz": (ngrid, 3), "xi0": (ngrid,), "normal0": (ngrid, 3), "a": (ngrid,), "f": (ngrid,)}
    for name, shape in shapes.items():
        value = getattr(cavity, name)
        assert value.shape == shape, name
        assert not value.flags.writeable, name
        np.testing.assert_array_equal(value, cavity.get(name))
        np.testing.assert_array_equal(value, getattr(cavity.snapshot(), name))
    np.testing.assert_allclose(np.linalg.norm(cavity.normal0, axis=1), 1.0, rtol=1e-12)


def test_native_coupling_retains_model_until_release(diatomic) -> None:
    """Native deletion must run before the owning model is collected."""
    import gc
    import weakref

    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    owner = weakref.ref(model._model)
    coupling = library.new_coupling(model._model)
    library.prepare_model_energy(model._model, coupling)
    del model
    gc.collect()
    assert owner() is not None
    # The coupling still walks the requests of its model.
    assert library.next_coupling_request(coupling)
    assert library.get_coupling_request_name(coupling) == "gaussian_potential"
    del coupling
    gc.collect()
    assert owner() is None


def test_getters_require_the_phase_they_were_staged_for(diatomic) -> None:
    """A ``get_*`` reads its own phase, never whichever answers happen to be fresh.

    ``prepare_response`` and ``prepare_gradient`` leave overlapping answers
    behind, so reading the wrong accessor used to succeed by accident.  The
    guard names both the phase the accessor needs and the one the coupling
    actually carries.
    """

    structure = diatomic()
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)

    unstaged = model.new_coupling()
    with raises(RuntimeError, match="is not staged"):
        _energy(model, unstaged)

    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    _answer_potential(coupling, np.zeros(model.cavity.ngrid))
    _energy(model, coupling)

    with raises(RuntimeError, match="staged for the energy phase"):
        _gradient(model, coupling)
    with raises(RuntimeError, match="staged for the energy phase"):
        model.get_response(coupling)

    # The same coupling, staged for the phase it is read in, is fine.
    model.prepare_response(coupling)
    model.get_response(coupling)
    with raises(RuntimeError, match="staged for the response phase"):
        _energy(model, coupling)


def test_each_phase_walks_what_it_still_owes(diatomic) -> None:
    """A walk visits a request only while the staged phase lacks one of its outputs.

    Nothing is owed before staging.  The energy phase owes the potential; the
    response phase keeps the valid energy answer and owes nothing more on a
    cavity that ignores the host; the gradient phase owes the two raw
    derivatives.
    """

    structure = diatomic()
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    ngrid = model.cavity.ngrid

    coupling = model.new_coupling()
    assert isinstance(coupling, Coupling)
    # An unstaged coupling has nothing to visit, as in Fortran
    assert list(coupling) == []

    model.prepare_energy(coupling)
    walked = []
    for request in coupling:
        walked.append((request.name, request.missing))
        coupling.answer(phi=np.zeros(ngrid))
    assert walked == [("gaussian_potential", {"phi"})]
    assert list(coupling) == []
    _energy(model, coupling)

    model.prepare_response(coupling)
    assert list(coupling) == []

    model.prepare_gradient(coupling)
    (request,) = list(coupling)
    assert request.missing == {"dphi_dr", "dphi_dxi"}


def test_native_cursor_queries_are_checked(diatomic) -> None:
    """The C request queries need a current request, and name what is wrong.

    Asking whether an undeclared output is missing is not an error -- it is
    simply not missing -- while answering one is.
    """

    structure = diatomic()
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(structure)
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    handle = coupling._handle
    no_current = r"No current coupling request - call next\(\) first"

    for query in (
        lambda: library.get_coupling_request_name(handle),
        lambda: library.get_coupling_request_missing(handle, "phi"),
        lambda: library.answer_coupling_request(handle, "phi", np.zeros(ngrid)),
    ):
        with raises(RuntimeError, match=no_current):
            query()

    model.prepare_energy(coupling)
    assert library.next_coupling_request(handle) is True
    assert library.get_coupling_request_name(handle) == "gaussian_potential"
    assert library.get_coupling_request_missing(handle, "phi") is True
    # Declared but not required in this phase, and not declared by this kind.
    assert library.get_coupling_request_missing(handle, "dphi_dr") is False
    assert library.get_coupling_request_missing(handle, "gt") is False
    with raises(RuntimeError, match="gaussian_potential has no output 'gt'"):
        library.answer_coupling_request(handle, "gt", np.zeros(ngrid))
    library.answer_coupling_request(handle, "phi", np.zeros(ngrid))
    assert library.get_coupling_request_missing(handle, "phi") is False

    # That was the only request: the pass ends and rewinds, and a new pass
    # finds nothing missing.
    assert library.next_coupling_request(handle) is False
    with raises(RuntimeError, match=no_current):
        library.get_coupling_request_name(handle)
    assert library.next_coupling_request(handle) is False
    # A failure is raised, never read as the end of a pass.
    with raises(RuntimeError, match="next_coupling_request"):
        library.next_coupling_request(library.CouplingHandle.null())


def test_an_empty_response_ends_the_walk_at_once(diatomic) -> None:
    """Absence is physics: a model without an item to hand back yields none."""

    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentPV(1.0e-4)])
    model.update(diatomic())
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    model.prepare_response(coupling)
    response = model.get_response(coupling)

    assert list(response) == []
    handle = model._response_handle()
    assert library.next_response_item(handle) is False
    # Nothing is current, so there is nothing to name or to read.
    with raises(RuntimeError, match="No current response item"):
        library.get_response_item_name(handle)
    with raises(RuntimeError, match="No current response item"):
        library.get_response_array(handle, "w_phi", np.empty(ngrid))
    # A failure is raised, never read as the end of a pass.
    with raises(RuntimeError, match="next_response_item"):
        library.next_response_item(library.ResponseHandle.null())


def test_response_arrays_are_read_by_name(gaussian_density) -> None:
    """The cursor walks the items; each array of the current one is read by name, grid axis first."""

    model = SolvationModel(
        CavityDROP(lsf=Isodensity(), nleb=26, source=gaussian_density),
        [ModelComponentCPCM(32.0), ModelComponentPV(1.0e-4)],
    )
    model.update(Structure(np.array([1]), np.zeros((1, 3))))
    ngrid = model.cavity.ngrid
    rng = np.random.default_rng(3)
    coupling = model.new_coupling()
    for prepare in (model.prepare_energy, model.prepare_response):
        prepare(coupling)
        for request in coupling:
            coupling.answer(**{name: rng.standard_normal((ngrid, *request._outputs[name]))
                               for name in request.missing})
    response = model.get_response(coupling)
    handle = model._response_handle()

    expected = {
        "potential_adjoint": {"w_phi": (ngrid,)},
        "density": {"w_rho": (ngrid,), "w_grad_rho": (ngrid, 3), "w_hess_rho": (ngrid, 3, 3)},
    }
    # Native order: the components' items, then the cavity's density item.
    assert [item.name for item in response] == list(expected)
    assert [type(item) for item in response] == [PotentialAdjointResponse, DensityResponse]
    # Building the value walked the native cursor to the end, which rewound it,
    # so this pass meets the same items again.
    for item in response:
        assert library.next_response_item(handle) is True
        assert library.get_response_item_name(handle) == item.name
        for array, shape in expected[item.name].items():
            value = getattr(item, array)
            assert value.shape == shape, array
            assert value.flags.c_contiguous, array
            copy = np.empty_like(value)
            library.get_response_array(handle, array, copy)
            np.testing.assert_array_equal(value, copy)
    assert library.next_response_item(handle) is False
    assert np.abs(_item(response, DensityResponse).w_hess_rho).max() > 0.0

    # An array the current item does not have is refused by name before
    # anything is written, whether or not another item has it.
    assert library.next_response_item(handle) is True
    buffer = np.empty(ngrid)
    with raises(RuntimeError, match="potential_adjoint has no array 'w_rho'"):
        library.get_response_array(handle, "w_rho", buffer)
    with raises(RuntimeError, match="potential_adjoint has no array 'q'"):
        library.get_response_array(handle, "q", buffer)
    assert library.next_response_item(handle) is True
    assert library.get_response_item_name(handle) == "density"
    with raises(RuntimeError, match="density has no array 'w_phi'"):
        library.get_response_array(handle, "w_phi", buffer)

    # Filling the response again starts a new walk, even in the middle of one.
    library.get_model_response(model._model, coupling._handle, handle)
    assert library.next_response_item(handle) is True
    assert library.get_response_item_name(handle) == "potential_adjoint"


def test_borrowed_model_cavity_rejects_standalone_updates(diatomic) -> None:
    """A model-owned cavity may be inspected but not rebuilt out of band."""

    structure = diatomic()
    model = SolvationModel(CavityDROP(nleb=26), [ModelComponentPV(1.0e-4)])
    model.update(structure)
    original_area = model.cavity.area

    with raises(RuntimeError, match="model-owned cavity"):
        model.cavity.update(structure)
    with raises(RuntimeError, match="borrowed cavity"):
        library.update_cavity(model.cavity._as_handle(), structure._mol)

    borrowed = library.error_check(library.lib.moist_get_model_cavity)(
        model._model.handle
    )
    borrowed_ref = library.ffi.new("moist_cavity *")
    borrowed_ref[0] = borrowed
    library.lib.moist_delete_cavity(borrowed_ref)

    assert borrowed_ref[0] == library.ffi.NULL
    assert model.cavity.area == original_area


def _gaussian_answers(ngrid):
    return {
        "phi": np.zeros(ngrid),
        "dphi_dr": np.zeros((ngrid, 3)),
        "dphi_dxi": np.zeros(ngrid),
    }


@pytest.mark.parametrize("output", ["dphi_dr", "dphi_dxi"])
def test_python_side_rejection_leaves_the_stored_answer(diatomic, output):
    """A value that is no number never reaches the library, so nothing stored changes."""
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    valid = _gaussian_answers(model.cavity.ngrid)
    for request in coupling:
        coupling.answer(**valid)
        # ``get_*`` leaves the cursor alone, so the request stays current.
        expected, _ = _gradient(model, coupling)
        with pytest.raises(ValueError):
            coupling.answer(**{output: ["invalid"]})
    assert list(coupling) == []
    actual, _ = _gradient(model, coupling)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("output", ["dphi_dr", "dphi_dxi"])
def test_native_rejection_leaves_the_output_missing(diatomic, output):
    """A native rejection un-answers exactly that output; the others survive."""
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    valid = _gaussian_answers(model.cavity.ngrid)
    for request in coupling:
        coupling.answer(**valid)
        expected, _ = _gradient(model, coupling)
        with pytest.raises(RuntimeError, match="finite"):
            coupling.answer(**{output: np.full(valid[output].shape, np.nan)})
    with pytest.raises(RuntimeError, match="missing required outputs: " + output):
        _gradient(model, coupling)

    # A second pass retries exactly the rejected output.
    retried = []
    for retry in coupling:
        retried.append(retry.missing)
        coupling.answer(**{output: valid[output]})
    assert retried == [{output}]
    actual, _ = _gradient(model, coupling)
    np.testing.assert_array_equal(actual, expected)


def test_rejected_replacement_is_missing_and_retryable(diatomic):
    model = SolvationModel(CavityISwiG(nleb=26), [ModelComponentCPCM(32.0)])
    model.update(diatomic())
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    ngrid = model.cavity.ngrid
    for request in coupling:
        # A valid answer to an answered output replaces the stored value.
        coupling.answer(phi=np.ones(ngrid))
        assert _energy(model, coupling) != 0.0
        coupling.answer(phi=np.zeros(ngrid))
        assert _energy(model, coupling) == 0.0
        # Python-side conversion: the stored answer is untouched.
        with pytest.raises(ValueError):
            coupling.answer(phi=["invalid"] * ngrid)
        assert _still_missing(coupling, request) == set()
        assert _energy(model, coupling) == 0.0
        # Native rejection: the output is missing until a valid retry.
        with pytest.raises(RuntimeError, match="finite"):
            coupling.answer(phi=np.full(ngrid, np.nan))
        assert _still_missing(coupling, request) == {"phi"}
        with pytest.raises(RuntimeError, match="missing"):
            _energy(model, coupling)
    retried = []
    for retry in coupling:
        retried.append(retry.missing)
        coupling.answer(phi=np.zeros(ngrid))
    assert retried == [{"phi"}]
    assert _energy(model, coupling) == 0.0


@pytest.mark.parametrize("cavity_type", [CavityISwiG, CavityDROP])
@pytest.mark.parametrize("component", ["pcm", "gostshyp", "combined", "disabled"])
def test_builtin_components_declare_their_requests(diatomic, cavity_type, component):
    """The walk is the model's declaration in component order; the moment width is GOSTSHYP's."""
    components = []
    expected = []
    if component in ("pcm", "combined"):
        components.append(ModelComponentCPCM(32.0))
        expected.append("gaussian_potential")
    if component in ("gostshyp", "combined", "disabled"):
        components.append(ModelComponentGOSTSHYP(0.0 if component == "disabled" else 1e-5))
        if component != "disabled":
            expected.append("gaussian_moments")
    model = SolvationModel(cavity_type(nleb=26), components)
    model.update(diatomic())
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    assert list(coupling) == []

    model.prepare_energy(coupling)
    walked = []
    for request in coupling:
        walked.append(request.name)
        if not isinstance(request, GaussianMomentRequest):
            continue
        assert tuple(request._outputs) == ("gt", "pt", "mt", "rt")
        assert request.missing == {"gt", "pt"}
        width = request.width
        a = model.cavity.a
        assert width.shape == (ngrid,)
        assert not width.flags.writeable
        active = a > 0.0
        np.testing.assert_allclose(width[active], np.pi * np.log(2.0) / a[active], rtol=1e-12)
        np.testing.assert_array_equal(width[~active], 0.0)
        # The snapshot is a copy of the request's own input.
        copy = np.empty(ngrid)
        library.get_coupling_request_width(coupling._handle, copy)
        np.testing.assert_array_equal(width, copy)
    assert walked == expected


def test_removed_protocol_names_stay_removed() -> None:
    """Requests have no handles, no count and no indexed list; responses no lookup by kind."""
    assert not hasattr(moist, "CPCMSolver")
    assert not hasattr(SolvationModel, "ngrid")
    for name in ("requests", "pending", "model"):
        assert not hasattr(Coupling, name)
    for name in ("answer", "handle", "n_missing", "outputs", "shapes"):
        assert not hasattr(moist.CouplingRequest, name)
    for name in ("potential_adjoint", "density", "gostshyp_amplitude", "items", "__len__"):
        assert not hasattr(Response, name)
    for name in ("coupling_n_requests", "coupling_request_handle", "coupling_request_missing_count",
                 "response_has", "response_get_potential_adjoint", "response_get_density",
                 "response_get_gostshyp_amplitude",
                 "coupling_next", "coupling_request_name", "coupling_request_missing",
                 "coupling_get_moment_width", "coupling_answer", "response_next",
                 "response_item_name", "response_get", "general_model_get_energy",
                 "general_model_prepare_energy"):
        assert not hasattr(library, name)
    for name in ("get_coupling_request_count", "get_coupling_request_handle",
                 "get_coupling_request_missing_count", "has_response",
                 "get_response_potential_adjoint", "get_response_density",
                 "get_response_gostshyp_amplitude"):
        assert not hasattr(library.lib, "moist_" + name)


# -----------------------------------------------------------------------------
# The cursor
# -----------------------------------------------------------------------------
#
# ``for request in coupling`` drives the native cursor, so a loop has the
# Fortran and C semantics: one visit per request and pass, in declaration
# order; the end of a loop rewinds; leaving a loop early resumes the pass;
# staging and a model update restart the walk.  A model with a potential and a
# moment request makes resuming and restarting tell apart.

#: The energy-phase walk of the ``two_requests`` model, in declaration order
BOTH = ["gaussian_potential", "gaussian_moments"]


@pytest.fixture
def two_requests(diatomic) -> SolvationModel:
    """ISwiG model whose phases walk a potential request, then a moment request."""
    model = SolvationModel(
        CavityISwiG(nleb=26), [ModelComponentCPCM(32.0), ModelComponentGOSTSHYP(1.0e-5)]
    )
    model.update(diatomic())
    return model


def _walk(coupling):
    """Names of the requests one loop over ``coupling`` visits, answering none."""
    return [request.name for request in coupling]


def test_cursor_visits_each_request_once_per_pass_in_declaration_order(two_requests):
    model = two_requests
    coupling = model.new_coupling()
    assert _walk(coupling) == []
    model.prepare_energy(coupling)
    # Nothing is answered, yet each request is visited once and the loop ends.
    assert _walk(coupling) == BOTH
    # The end of the pass rewound the cursor, so the next loop is a new pass.
    assert _walk(coupling) == BOTH


def test_cursor_second_pass_retries_a_rejected_output(two_requests):
    """A second loop visits only what is still missing, and nothing once all is answered."""
    model = two_requests
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    for request in coupling:
        if isinstance(request, GaussianPotentialRequest):
            with raises(RuntimeError, match="gaussian_potential: phi must be finite"):
                coupling.answer(phi=np.full(ngrid, np.nan))
        else:
            coupling.answer(gt=np.ones(ngrid), pt=np.ones((ngrid, 3)))

    retried = []
    for request in coupling:
        retried.append((request.name, request.missing))
        coupling.answer(phi=np.zeros(ngrid))
    assert retried == [("gaussian_potential", {"phi"})]
    assert _walk(coupling) == []
    assert np.isfinite(_energy(model, coupling))


@pytest.mark.parametrize("leave", ["break", "raise"])
def test_cursor_leaving_a_loop_early_resumes_the_pass(two_requests, leave):
    model = two_requests
    coupling = model.new_coupling()
    model.prepare_energy(coupling)

    try:
        for request in coupling:
            if leave == "break":
                break
            raise KeyError(request.name)
    except KeyError:
        pass
    assert request.name == "gaussian_potential"
    # The loop left the cursor on that request, which is still current...
    assert library.get_coupling_request_name(coupling._handle) == "gaussian_potential"
    # ...and the next loop resumes after it instead of starting over.
    assert _walk(coupling) == ["gaussian_moments"]
    # That loop finished the pass; the one after starts a new pass.
    assert _walk(coupling) == BOTH


@pytest.mark.parametrize("phase", ["energy", "response", "gradient"])
def test_cursor_staging_restarts_a_half_walked_pass(two_requests, phase):
    model = two_requests
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    for request in coupling:
        break
    getattr(model, f"prepare_{phase}")(coupling)
    with raises(RuntimeError, match="No current coupling request"):
        coupling.answer(phi=np.zeros(model.cavity.ngrid))
    # The walk starts over at the first request instead of resuming after it.
    assert _walk(coupling) == BOTH


def test_cursor_model_update_ends_the_walk(two_requests, diatomic):
    """An update unstages every coupling: no current request and nothing to walk."""
    model = two_requests
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    for request in coupling:
        break
    model.update(diatomic())
    with raises(RuntimeError, match="No current coupling request"):
        coupling.answer(phi=np.zeros(model.cavity.ngrid))
    assert _walk(coupling) == []
    with raises(RuntimeError, match="requires a coupling staged by prepare_energy"):
        _energy(model, coupling)
    model.prepare_energy(coupling)
    assert _walk(coupling) == BOTH


def test_cursor_answer_needs_a_current_request(two_requests):
    """``answer`` acts on the current request only, and says so when there is none."""
    model = two_requests
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    phi = np.zeros(ngrid)
    no_current = r"No current coupling request - call next\(\) first"

    with raises(RuntimeError, match=no_current):  # never staged
        coupling.answer(phi=phi)
    model.prepare_energy(coupling)
    with raises(RuntimeError, match=no_current):  # staged, no loop yet
        coupling.answer(phi=phi)
    # Even with nothing to submit, the call needs a request to address.
    with raises(RuntimeError, match=no_current):
        coupling.answer()

    for request in coupling:
        break
    # A loop left early keeps its request current, so it can be answered after.
    coupling.answer(phi=phi)
    assert _still_missing(coupling, request) == set()
    assert _walk(coupling) == ["gaussian_moments"]
    with raises(RuntimeError, match=no_current):  # the pass has ended
        coupling.answer(gt=np.ones(ngrid))


def test_answer_validates_keywords(two_requests):
    """Unknown keywords and values that are no numbers are refused before moist sees them."""
    model = two_requests
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    visited = []
    for request in coupling:
        visited.append(request.name)
        if isinstance(request, GaussianPotentialRequest):
            with raises(TypeError, match="gaussian_potential has no output 'gt'"):
                coupling.answer(gt=np.zeros(ngrid))
        else:
            with raises(TypeError, match="gaussian_moments has no output 'phi'"):
                coupling.answer(phi=np.zeros(ngrid))
        with raises(ValueError):
            coupling.answer(**{next(iter(request._outputs)): ["invalid"] * ngrid})
        # None of the refusals stored or un-answered anything.
        assert _still_missing(coupling, request) == request.missing
    assert visited == BOTH


def test_answer_rejects_a_wrong_shape(two_requests):
    """An ndarray carries its shape, so a wrong one is refused before it reaches moist."""
    model = two_requests
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    model.prepare_gradient(coupling)
    visited = []
    for request in coupling:
        visited.append(request.name)
        for name, extents in request._outputs.items():
            good = (ngrid, *extents)
            for bad in ((ngrid - 1, *extents), (ngrid, *extents, 3), (ngrid,) if extents else (ngrid, 3)):
                with raises(ValueError, match=rf"{request.name}: {name} must have shape"):
                    coupling.answer(**{name: np.zeros(bad)})
        # A refused shape stores nothing; the request is still fully missing.
        assert _still_missing(coupling, request) == request.missing
        # The outputs before a refusal are kept, the rest of the call is skipped.
        first, second = list(request._outputs)[:2]
        with raises(ValueError, match=second):
            coupling.answer(**{first: np.zeros((ngrid, *request._outputs[first])),
                               second: np.zeros((ngrid - 1, *request._outputs[second]))})
        assert _still_missing(coupling, request) == request.missing - {first}
    assert visited == BOTH


def test_request_snapshots_are_values(two_requests, diatomic):
    """A snapshot keeps what the cursor saw: answers and updates leave it alone."""
    model = two_requests
    ngrid = model.cavity.ngrid
    coupling = model.new_coupling()
    model.prepare_energy(coupling)
    snapshots = []
    for request in coupling:
        snapshots.append(request)
        coupling.answer(**{name: np.ones((ngrid, *request._outputs[name]))
                           for name in request.missing})
        # The coupling holds the answers now; the snapshot still shows what
        # was missing when the cursor got here.
        assert _still_missing(coupling, request) == set()
    potential, moments = snapshots
    assert potential.missing == {"phi"}
    assert moments.missing == {"gt", "pt"}
    with raises(AttributeError):
        potential.missing = frozenset()
    width = moments.width.copy()
    with raises(ValueError, match="WRITEABLE"):
        moments.width.setflags(write=True)

    # A new surface gives the request new widths; the snapshot keeps its copy.
    model.update(diatomic(1.6))
    model.prepare_energy(coupling)
    (fresh,) = [request for request in coupling if isinstance(request, GaussianMomentRequest)]
    assert not np.array_equal(fresh.width, width)
    np.testing.assert_array_equal(moments.width, width)
    assert moments.missing == {"gt", "pt"}
