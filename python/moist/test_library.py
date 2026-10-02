import functools
from types import SimpleNamespace

import numpy as np
import pytest
from pytest import approx, raises

from moist.interface import (
    CavityDROPIsodensityCallback,
    CavityDROPIsodensityInternal,
    GeneralSolvationModel,
    ModelComponentPV,
    Structure,
)
from moist.library import _callback_takes_order, get_api_version


def test_api_version_format() -> None:
    version = get_api_version()
    parts = version.split(".")
    assert len(parts) == 3
    assert all(part.isdigit() for part in parts)


def test_callback_arity_detection() -> None:
    """Both callback forms are recognised without being called"""

    def legacy(point):
        return 0.0, np.zeros(3)

    def with_order(point, order):
        return 0.0, np.zeros(3)

    def keyword_order(point, *, order=1):
        return 0.0, np.zeros(3)

    def varargs(*args):
        return 0.0, np.zeros(3)

    assert not _callback_takes_order(legacy)
    assert _callback_takes_order(with_order)
    # `order` keyword-only: the two-positional call would fail, so legacy form
    assert not _callback_takes_order(keyword_order)
    # *args could accept anything; safe assumption is the legacy form
    assert not _callback_takes_order(varargs)
    # Bound partial that already consumed `point` still exposes `order`
    assert _callback_takes_order(functools.partial(with_order))
    # Builtins are not introspectable: fall back, do not raise
    assert not _callback_takes_order(len)


#: Density contour of the test surfaces, in electrons/Bohr^3
#: - owned by moist, so it reaches the cavity constructor, not the callback
_RHO_ISO = 1.0e-3


class _GaussianDensity:
    """rho(r) = sum_A c exp(-a |r - R_A|^2), returned to moist

    - moist forms the level set S = scale (rho_iso - rho) from it
    - records the derivative orders asked for, so a test can tell that moist
      skips the expensive orders during projection
    """

    def __init__(self, centers: np.ndarray, alpha: float = 0.3):
        # Diagonal single-s-per-atom density of test/api/example.c, known to
        # give a well-behaved surface
        coeff = (2.0 * alpha / np.pi) ** 0.75
        self.centers = np.asarray(centers, dtype=np.float64)
        self.c = 2.0 * coeff**2
        self.a = 2.0 * alpha
        self.rho_iso = _RHO_ISO
        self.orders: set[int] = set()
        self.calls = 0

    def _derivatives(self, point, order):
        self.calls += 1
        self.orders.add(order)

        point = np.asarray(point, dtype=np.float64)
        d = point[None, :] - self.centers
        g = self.c * np.exp(-self.a * np.einsum("ai,ai->a", d, d))

        rho = g.sum()
        drho = np.einsum("a,ai->i", -2.0 * self.a * g, d)

        if order < 2:
            return rho, drho

        eye = np.eye(3)
        d2rho = np.einsum("a,ai,aj->ij", 4.0 * self.a**2 * g, d, d) - 2.0 * self.a * g.sum() * eye
        if order < 3:
            return rho, drho, d2rho

        d3rho = np.einsum("a,ai,aj,ak->ijk", -8.0 * self.a**3 * g, d, d, d)
        d3rho += 4.0 * self.a**2 * (
            np.einsum("a,ai,jk->ijk", g, d, eye)
            + np.einsum("a,aj,ik->ijk", g, d, eye)
            + np.einsum("a,ak,ij->ijk", g, d, eye)
        )
        return rho, drho, d2rho, d3rho

    def with_order(self, point, order):
        return self._derivatives(point, order)

    def legacy(self, point):
        """Pre-order form, always computes everything"""
        return self._derivatives(point, 3)


class _GaussianSource(_GaussianDensity):
    """Density provider matching the public isodensity source interface"""

    scale = 750.0

    def density(self, point, order):
        return self._derivatives(point, order)


@pytest.fixture
def water() -> tuple[np.ndarray, np.ndarray]:
    numbers = np.array([8, 1, 1])
    positions = np.array(
        [
            [0.0000, 0.0000, 0.1173],
            [0.0000, 1.4309, -0.9370],
            [0.0000, -1.4309, -0.9370],
        ]
    )
    return numbers, positions


@pytest.fixture
def flaky_lsf(water) -> SimpleNamespace:
    """Density callback that can be switched to fail partway into a build

    - ``arm(True)`` raises once moist is 50 grid points in, far enough that the
      abort unwinds a partly built surface
    - ``arm`` also resets the call counter: trip point is measured from the
      start of the next build, not the fixture's lifetime
    """
    _, positions = water
    lsf = _GaussianDensity(positions)
    state = SimpleNamespace(lsf=lsf, rho_iso=lsf.rho_iso, failing=False)

    def callback(point, order):
        if state.failing and lsf.calls >= 50:
            raise RuntimeError("no density here")
        return lsf._derivatives(point, order)

    def arm(fail: bool = True) -> None:
        state.failing = fail
        lsf.calls = 0

    state.callback = callback
    state.arm = arm
    return state


def _build(callback, water, rho_iso=_RHO_ISO, **kwargs):
    numbers, positions = water
    structure = Structure(numbers, positions)
    cavity = CavityDROPIsodensityCallback(callback, rho_iso=rho_iso, nleb=26, **kwargs)
    cavity.update(structure)
    return cavity.cavity


def test_isodensity_cavity_accepts_a_density_source(water) -> None:
    """Cavity constructor owns source adaptation and parameter matching"""
    numbers, positions = water
    source = _GaussianSource(positions)
    cavity = CavityDROPIsodensityCallback(source, nleb=26)

    cavity.update(Structure(numbers, positions))

    assert cavity.snapshot().ngrid > 0
    assert source.calls > 0
    with raises(ValueError, match="scale must match"):
        CavityDROPIsodensityCallback(source, nleb=26, scale=1000.0)
    with raises(ValueError, match="rho_iso must match"):
        CavityDROPIsodensityCallback(source, nleb=26, rho_iso=2.0e-3)
    with raises(TypeError, match="rho_iso must be given"):
        CavityDROPIsodensityCallback(source.with_order, nleb=26)


def test_isodensity_cavity_retains_callback_keyword_compatibility(water) -> None:
    _, positions = water
    source = _GaussianDensity(positions)

    with pytest.deprecated_call(match="source"):
        cavity = CavityDROPIsodensityCallback(callback=source.with_order, rho_iso=source.rho_iso, nleb=26)

    assert isinstance(cavity, CavityDROPIsodensityCallback)


def test_callback_both_forms_agree(water) -> None:
    """Legacy one-argument callback gives the same cavity

    The forms differ only in whether moist can tell the callback to skip
    unneeded derivatives, so the surface is the same
    """
    _, positions = water

    new_lsf = _GaussianDensity(positions)
    new_cavity = _build(new_lsf.with_order, water)

    old_lsf = _GaussianDensity(positions)
    old_cavity = _build(old_lsf.legacy, water)

    assert new_cavity.ngrid == old_cavity.ngrid
    assert new_cavity.area == approx(old_cavity.area, rel=1.0e-13)
    assert new_cavity.volume == approx(old_cavity.volume, rel=1.0e-13)

    # Order-aware form was spared some work
    assert new_lsf.orders == {1, 2}
    assert old_lsf.orders == {3}


def test_callback_order_override(water) -> None:
    """`pass_order` overrides introspection in both directions"""
    _, positions = water

    # Two-argument callback forced to the legacy form would raise TypeError
    # inside the CFFI trampoline, so it is rejected up front
    lsf = _GaussianDensity(positions)
    with raises(TypeError):
        _build(lsf.with_order, water, pass_order=False)

    # Forcing legacy form on a legacy callback: no-op
    lsf = _GaussianDensity(positions)
    cavity = _build(lsf.legacy, water, pass_order=False)
    assert cavity.ngrid > 0
    assert lsf.orders == {3}

    # Forcing order form on a callback that accepts it: works
    lsf = _GaussianDensity(positions)
    cavity = _build(lsf.with_order, water, pass_order=True)
    assert cavity.ngrid > 0
    assert lsf.orders == {1, 2}


def test_callback_missing_derivative_is_reported(water) -> None:
    """Callback omitting a requested derivative is reported clearly"""
    _, positions = water
    lsf = _GaussianDensity(positions)

    def truncated(point, order):
        # Never returns a Hessian, whatever moist asks
        value, grad = lsf._derivatives(point, 1)
        return value, grad

    with raises(ValueError, match="did not return the requested Hessian"):
        _build(truncated, water)


def test_callback_failure_aborts_build(water) -> None:
    """Raising callback aborts the build and surfaces the real exception

    - raised on the 50th evaluation, not the first
    - moist evaluates the level set inside OpenMP parallel loops: a mid-loop
      abort exercises the failure channel
    - a first-call failure would pass even with broken parallel handling
    """
    _, positions = water
    lsf = _GaussianDensity(positions)

    class Boom(RuntimeError):
        pass

    def flaky(point, order):
        if lsf.calls >= 50:
            raise Boom("the host density is unavailable here")
        return lsf._derivatives(point, order)

    with raises(Boom, match="the host density is unavailable here"):
        _build(flaky, water)

    # Exception arrives with its own traceback, not a moist error rewrapped
    # around a downstream symptom
    try:
        lsf.calls = 0
        _build(flaky, water)
    except Boom as exc:
        frames = []
        tb = exc.__traceback__
        while tb is not None:
            frames.append(tb.tb_frame.f_code.co_name)
            tb = tb.tb_next
        assert "flaky" in frames
    else:  # pragma: no cover - the call above must raise
        raise AssertionError("a raising callback still produced a cavity")


def test_callback_failure_stops_further_calls(water) -> None:
    """Failed callback is not called again"""
    _, positions = water
    lsf = _GaussianDensity(positions)
    seen = []

    def flaky(point, order):
        seen.append(order)
        if len(seen) > 50:
            raise RuntimeError("no density here")
        return lsf._derivatives(point, order)

    with raises(RuntimeError, match="no density here"):
        _build(flaky, water)

    # A build ignoring the failure would keep calling for the whole grid;
    # the abort unwinds promptly instead
    assert 50 < len(seen) < 500


def test_callback_failure_is_not_sticky(water, flaky_lsf) -> None:
    """Cavity whose callback failed once rebuilds cleanly"""
    structure = Structure(*water)
    cavity = CavityDROPIsodensityCallback(flaky_lsf.callback, rho_iso=flaky_lsf.rho_iso, nleb=26)
    flaky_lsf.arm()

    with raises(RuntimeError, match="no density here"):
        cavity.update(structure)

    flaky_lsf.arm(False)
    cavity.update(structure)
    assert cavity.cavity.ngrid > 0


def test_failed_rebuild_invalidates_previous_cavity_results(water, flaky_lsf) -> None:
    """Failed second build does not leave the first surface readable as current"""
    structure = Structure(*water)
    cavity = CavityDROPIsodensityCallback(flaky_lsf.callback, rho_iso=flaky_lsf.rho_iso, nleb=26)
    cavity.update(structure)
    assert cavity.snapshot().ngrid > 0

    flaky_lsf.arm()
    with raises(RuntimeError, match="no density here"):
        cavity.update(structure)

    with raises(RuntimeError, match="successfully updated"):
        cavity.snapshot()


def test_failed_model_rebuild_invalidates_its_cavity_view(water, flaky_lsf) -> None:
    """Model updates propagate callback failures and invalidate the live view"""
    structure = Structure(*water)
    model = GeneralSolvationModel(
        CavityDROPIsodensityCallback(flaky_lsf.callback, rho_iso=flaky_lsf.rho_iso, nleb=26), [ModelComponentPV(1.0e-4)]
    )
    model.update(structure)
    assert model.cavity.snapshot().ngrid > 0

    flaky_lsf.arm()
    with raises(RuntimeError, match="no density here"):
        model.update(structure)

    with raises(RuntimeError, match="successfully updated"):
        model.cavity.snapshot()


# -- internal isodensity backend ----------------------------------------------

_ALPHA = 0.3


def _s_basis(natoms: int, alpha: float = _ALPHA) -> dict:
    """One normalized s primitive per atom

    With the diagonal density ``2 * eye``: exactly the density
    :class:`_GaussianDensity` returns, as in test/api/example.c
    """
    return dict(
        shell_atom=np.arange(natoms),
        shell_l=np.zeros(natoms, dtype=int),
        shell_nprim=np.ones(natoms, dtype=int),
        exps=np.full(natoms, alpha),
        coeffs=np.full(natoms, (2.0 * alpha / np.pi) ** 0.75),
    )


def _internal(water, **kwargs) -> CavityDROPIsodensityInternal:
    numbers, _ = water
    return CavityDROPIsodensityInternal(
        **_s_basis(len(numbers)), rho_iso=_RHO_ISO, nleb=26, **kwargs
    )


def test_internal_isodensity_layout_matches_the_basis(water) -> None:
    layout = _internal(water).layout
    assert (layout.ncart, layout.nshell) == (3, 3)
    assert layout.shell_offset.tolist() == [0, 1, 2, 3]
    assert not layout.powers.any()

    # p and d shell: lx descending, then ly descending
    mixed = CavityDROPIsodensityInternal(
        shell_atom=[0, 0],
        shell_l=[1, 2],
        shell_nprim=[1, 2],
        exps=[0.5, 1.0, 0.2],
        coeffs=[1.0, 0.6, 0.4],
        rho_iso=_RHO_ISO,
    )
    assert mixed.ncart == 9
    assert mixed.layout.shell_offset.tolist() == [0, 3, 9]
    assert mixed.layout.powers.tolist() == [
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2],
    ]
    with raises(ValueError):
        mixed.layout.powers[0, 0] = 3


def test_internal_isodensity_layouts_compare_by_value(water) -> None:
    """Layouts compare by value: equal on one basis, unequal on another"""
    first, second = _internal(water).layout, _internal(water).layout
    assert first is not second
    assert first == second
    assert hash(first) == hash(second)
    assert len({first, second}) == 1

    wider = CavityDROPIsodensityInternal(**_s_basis(4), rho_iso=_RHO_ISO).layout
    assert first != wider
    assert first != "not a layout"


def test_internal_isodensity_matches_the_callback_backend(water) -> None:
    """Same density through both backends gives the same surface"""
    numbers, positions = water
    structure = Structure(numbers, positions)

    internal = _internal(water)
    internal.set_density(2.0 * np.eye(internal.ncart))
    internal.update(structure)

    reference = CavityDROPIsodensityCallback(
        _GaussianDensity(positions).with_order, rho_iso=_RHO_ISO, nleb=26
    )
    reference.update(structure)

    got, ref = internal.snapshot(), reference.snapshot()
    assert got.ngrid == ref.ngrid
    assert got.area == approx(ref.area, abs=1.0e-8)
    assert got.volume == approx(ref.volume, abs=1.0e-8)
    assert np.allclose(got.xyz, ref.xyz, rtol=0.0, atol=1.0e-8)


def test_internal_isodensity_refuses_to_build_without_a_density(water) -> None:
    """No density is a clean error; installing one afterwards recovers"""
    structure = Structure(*water)
    cavity = _internal(water)
    with raises(RuntimeError, match="density matrix has not been installed"):
        cavity.update(structure)

    cavity.set_density(2.0 * np.eye(cavity.ncart))
    cavity.update(structure)
    assert cavity.snapshot().ngrid > 0


def test_internal_isodensity_rejects_inconsistent_input(water) -> None:
    cavity = _internal(water)
    with raises(RuntimeError, match="does not match the basis"):
        cavity.set_density(np.eye(cavity.ncart + 1))
    with raises(ValueError, match="square"):
        cavity.set_density(np.ones((3, 2)))

    basis = _s_basis(3)
    with raises(ValueError, match="sum\\(shell_nprim\\)"):
        CavityDROPIsodensityInternal(**{**basis, "exps": basis["exps"][:2]}, rho_iso=_RHO_ISO)
    with raises(ValueError, match="same length"):
        CavityDROPIsodensityInternal(**{**basis, "shell_l": [0, 0]}, rho_iso=_RHO_ISO)
    with raises(ValueError, match="non-negative"):
        CavityDROPIsodensityInternal(**{**basis, "shell_atom": [-1, 0, 1]}, rho_iso=_RHO_ISO)
    with raises(RuntimeError, match="angular momentum"):
        CavityDROPIsodensityInternal(**{**basis, "shell_l": [0, 0, 9]}, rho_iso=_RHO_ISO)

    # Shells on a missing atom would be read out of bounds; moist itself
    # refuses, so the C API is covered too
    wide = CavityDROPIsodensityInternal(**_s_basis(4), rho_iso=_RHO_ISO, nleb=26)
    wide.set_density(2.0 * np.eye(wide.ncart))
    with raises(RuntimeError, match="needs 4 atoms, but the structure has only 3"):
        wide.update(Structure(*water))
    with raises(RuntimeError, match="needs 4 atoms, but the structure has only 3"):
        GeneralSolvationModel(wide, [ModelComponentPV(1.0e-4)]).update(Structure(*water))


def test_internal_isodensity_model_follows_the_latest_density(water) -> None:
    """Model cavity copy tracks the density, whichever object received it"""
    structure = Structure(*water)
    pressure = 1.0e-4
    eye = np.eye(3)

    def standalone_volume(dcart):
        cavity = _internal(water)
        cavity.set_density(dcart)
        cavity.update(structure)
        return cavity.volume

    volume_2 = standalone_volume(2.0 * eye)
    volume_3 = standalone_volume(3.0 * eye)
    assert volume_3 > volume_2

    # Density installed before the model is built: carried by the copy
    cavity = _internal(water)
    cavity.set_density(2.0 * eye)
    model = GeneralSolvationModel(cavity, [ModelComponentPV(pressure)])
    assert model.cavity.density_dependent
    assert model.cavity.layout == cavity.layout
    model.update(structure)
    assert model.cavity.volume == approx(volume_2, rel=1.0e-12)
    assert model.energy == approx(pressure * volume_2, rel=1.0e-12)

    # Set after the model was built, on the host-held object
    cavity.set_density(3.0 * eye)
    model.update(structure)
    assert model.cavity.volume == approx(volume_3, rel=1.0e-12)

    # Set through the model's own view
    model.cavity.set_density(2.0 * eye)
    model.update(structure)
    assert model.cavity.volume == approx(volume_2, rel=1.0e-12)

    # Model built before any density: refuses, then recovers
    fresh = GeneralSolvationModel(_internal(water), [ModelComponentPV(pressure)])
    with raises(RuntimeError, match="density matrix has not been installed"):
        fresh.update(structure)
    with raises(RuntimeError, match="successfully updated"):
        fresh.cavity.snapshot()
    fresh.cavity.set_density(2.0 * eye)
    fresh.update(structure)
    assert fresh.cavity.volume == approx(volume_2, rel=1.0e-12)

    # Model built from another model's view: follows the same density
    nested = GeneralSolvationModel(model.cavity, [ModelComponentPV(pressure)])
    assert nested.cavity.density_owner is cavity
    nested.cavity.set_density(3.0 * eye)
    nested.update(structure)
    assert nested.cavity.volume == approx(volume_3, rel=1.0e-12)
    model.update(structure)
    assert model.cavity.volume == approx(volume_3, rel=1.0e-12)
