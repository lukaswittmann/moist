import functools

import numpy as np
import pytest
from pytest import raises

from moist.interface import IsodensityDROPCavity, Structure
from moist.library import _callback_takes_order, get_api_version


def test_api_version_format() -> None:
    version = get_api_version()
    parts = version.split(".")
    assert len(parts) == 3
    assert all(part.isdigit() for part in parts)


def test_callback_arity_detection() -> None:
    """Both callback forms must be recognised without being called."""

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
    # `order` is keyword-only here, so the two-positional-argument call would
    # fail -- it must be treated as the legacy form
    assert not _callback_takes_order(keyword_order)
    # *args could accept anything; the safe assumption is the legacy form
    assert not _callback_takes_order(varargs)
    # A bound partial that has already consumed `point` still exposes `order`
    assert _callback_takes_order(functools.partial(with_order))
    # Builtins are not introspectable and must fall back, not raise
    assert not _callback_takes_order(len)


class _GaussianLSF:
    """rho(r) = sum_A c exp(-a |r - R_A|^2), level set S = rho_iso - rho.

    Records which derivative orders it was asked for, so a test can tell that
    moist really does skip the expensive orders during projection.
    """

    def __init__(self, centers: np.ndarray, alpha: float = 0.3):
        # Mirrors the diagonal single-s-per-atom density used by the C example
        # in test/api/example.c, which is known to give a well-behaved surface.
        coeff = (2.0 * alpha / np.pi) ** 0.75
        self.centers = np.asarray(centers, dtype=np.float64)
        self.c = 2.0 * coeff**2
        self.a = 2.0 * alpha
        self.rho_iso = 1.0e-3
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

        value = self.rho_iso - rho
        grad = -drho
        if order < 2:
            return value, grad

        eye = np.eye(3)
        d2rho = np.einsum("a,ai,aj->ij", 4.0 * self.a**2 * g, d, d) - 2.0 * self.a * g.sum() * eye
        if order < 3:
            return value, grad, -d2rho

        d3rho = np.einsum("a,ai,aj,ak->ijk", -8.0 * self.a**3 * g, d, d, d)
        d3rho += 4.0 * self.a**2 * (
            np.einsum("a,ai,jk->ijk", g, d, eye)
            + np.einsum("a,aj,ik->ijk", g, d, eye)
            + np.einsum("a,ak,ij->ijk", g, d, eye)
        )
        return value, grad, -d2rho, -d3rho

    def with_order(self, point, order):
        return self._derivatives(point, order)

    def legacy(self, point):
        """The pre-order form: always computes everything."""
        return self._derivatives(point, 3)


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


def _build(callback, water, **kwargs):
    numbers, positions = water
    structure = Structure(numbers, positions)
    cavity = IsodensityDROPCavity(callback, nleb=26, **kwargs)
    cavity.update(structure)
    return cavity.cavity


def test_callback_both_forms_agree(water) -> None:
    """A legacy one-argument callback must still work and give the same cavity.

    The two forms differ only in whether moist can tell the callback to skip
    computing derivatives it does not need, so the resulting surface must be
    identical -- not merely close.
    """
    _, positions = water

    new_lsf = _GaussianLSF(positions)
    new_cavity = _build(new_lsf.with_order, water)

    old_lsf = _GaussianLSF(positions)
    old_cavity = _build(old_lsf.legacy, water)

    assert new_cavity.ngrid == old_cavity.ngrid
    assert new_cavity.area == old_cavity.area
    assert new_cavity.volume == old_cavity.volume

    # The order-aware form must actually have been spared some work
    assert new_lsf.orders == {1, 2}
    assert old_lsf.orders == {3}


def test_callback_order_override(water) -> None:
    """`pass_order` overrides introspection in both directions."""
    _, positions = water

    # A two-argument callback forced into the legacy call form would raise
    # TypeError inside the CFFI trampoline, so this must be rejected up front
    lsf = _GaussianLSF(positions)
    with raises(TypeError):
        _build(lsf.with_order, water, pass_order=False)

    # Forcing the legacy form on a callback that really is legacy is a no-op
    lsf = _GaussianLSF(positions)
    cavity = _build(lsf.legacy, water, pass_order=False)
    assert cavity.ngrid > 0
    assert lsf.orders == {3}

    # Forcing the order form on a callback that accepts it works too
    lsf = _GaussianLSF(positions)
    cavity = _build(lsf.with_order, water, pass_order=True)
    assert cavity.ngrid > 0
    assert lsf.orders == {1, 2}


def test_callback_missing_derivative_is_reported(water) -> None:
    """A callback that omits a requested derivative must say so clearly."""
    _, positions = water
    lsf = _GaussianLSF(positions)

    def truncated(point, order):
        # Never returns a Hessian, whatever moist asks for
        value, grad = lsf._derivatives(point, 1)
        return value, grad

    with raises(ValueError, match="did not return the requested Hessian"):
        _build(truncated, water)


def test_callback_failure_aborts_build(water) -> None:
    """A raising callback must abort the build and surface the real exception.

    The failure is raised on the 50th evaluation rather than the first: moist
    evaluates the level set from inside OpenMP parallel loops, so a mid-loop
    abort is what actually exercises the failure channel. A first-call failure
    would pass even if the parallel handling were broken.
    """
    _, positions = water
    lsf = _GaussianLSF(positions)

    class Boom(RuntimeError):
        pass

    def flaky(point, order):
        if lsf.calls >= 50:
            raise Boom("the host density is unavailable here")
        return lsf._derivatives(point, order)

    with raises(Boom, match="the host density is unavailable here"):
        _build(flaky, water)

    # The exception must arrive with its own traceback, not as a moist error
    # rewrapped around a downstream symptom.
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
    """Once the callback has failed, moist must stop calling it."""
    _, positions = water
    lsf = _GaussianLSF(positions)
    seen = []

    def flaky(point, order):
        seen.append(order)
        if len(seen) > 50:
            raise RuntimeError("no density here")
        return lsf._derivatives(point, order)

    with raises(RuntimeError, match="no density here"):
        _build(flaky, water)

    # A build that ignored the failure would keep calling for the whole grid;
    # the abort has to unwind promptly instead.
    assert 50 < len(seen) < 500


def test_callback_failure_is_not_sticky(water) -> None:
    """A cavity whose callback failed once must rebuild cleanly afterwards."""
    numbers, positions = water
    lsf = _GaussianLSF(positions)
    fail = [True]

    def flaky(point, order):
        if fail[0] and lsf.calls >= 50:
            raise RuntimeError("no density here")
        return lsf._derivatives(point, order)

    structure = Structure(numbers, positions)
    cavity = IsodensityDROPCavity(flaky, nleb=26)

    with raises(RuntimeError, match="no density here"):
        cavity.update(structure)

    fail[0] = False
    cavity.update(structure)
    assert cavity.cavity.ngrid > 0
