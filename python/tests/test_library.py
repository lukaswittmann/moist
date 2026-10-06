import functools
import importlib
from pathlib import Path
import runpy
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from pytest import approx, raises

from moist import DROPParameters, Isodensity, IsodensityParameters
from moist.interface import (
    CavityDROP,
    ModelComponentPV,
    SolvationModel,
    Structure,
)
from moist.library import _callback_takes_order, get_api_version
from moist import Context

#: Run context shared by every cavity and model in this module
CONTEXT = Context()

REQUIREMENTS_SCRIPT = Path(__file__).parents[1] / "check_requirements.py"


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
    # A keyword-only `order` is treated as the legacy form.
    assert not _callback_takes_order(keyword_order)
    # An arbitrary *args signature is treated as the legacy form.
    assert not _callback_takes_order(varargs)
    # A partial with `point` already bound still exposes `order`.
    assert _callback_takes_order(functools.partial(with_order))
    # Builtins are not introspectable and must fall back, not raise.
    assert not _callback_takes_order(len)


#: Density contour the test surfaces follow, in electrons/Bohr^3.
_RHO_ISO = 1.0e-3


class _GaussianDensity:
    """rho(r) = sum_A c exp(-a |r - R_A|^2), returned to moist.

    moist forms the level set S = scale (rho_iso - rho) from it.

    Records which derivative orders were requested; used to verify that
    moist skips the expensive orders during projection.
    """

    def __init__(self, centers: np.ndarray, alpha: float = 0.3):
        # Diagonal single-s-per-atom density, matching the C example in
        # test/api/example.c.
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
        """The pre-order form: always computes everything."""
        return self._derivatives(point, 3)


class _GaussianSource(_GaussianDensity):
    """Density provider matching the public isodensity source interface."""

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
    """A density callback that can be switched to fail partway into a build.

    ``arm(True)`` makes the callback raise once moist is 50 grid points in,
    and resets the call counter.
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


def _isodensity(source, rho_iso=_RHO_ISO, **kwargs):
    """A callback-backed DROP cavity; the contour is an LSF setting, not the source's."""
    return CavityDROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=rho_iso)),
                      parameters=DROPParameters(nleb=26), source=source, context=CONTEXT,
                      **kwargs)


def _build(callback, water, rho_iso=_RHO_ISO, **kwargs):
    numbers, positions = water
    structure = Structure(numbers, positions)
    cavity = _isodensity(callback, rho_iso, **kwargs)
    cavity.update(structure)
    return cavity.cavity


def test_isodensity_cavity_accepts_a_density_source(water) -> None:
    """A source is a callable or an object with ``density(point, order)``."""
    numbers, positions = water
    source = _GaussianSource(positions)
    cavity = _isodensity(source)

    cavity.update(Structure(numbers, positions))

    assert cavity.snapshot().ngrid > 0
    assert source.calls > 0
    assert cavity.lsf == Isodensity(parameters=IsodensityParameters(rho_iso=_RHO_ISO))
    with raises(TypeError, match="callable source"):
        _isodensity(object())


def test_callback_both_forms_agree(water) -> None:
    """A legacy one-argument callback must still work and give the same cavity.

    The two forms differ only in whether moist can tell the callback to skip
    computing derivatives it does not need.
    """
    _, positions = water

    new_lsf = _GaussianDensity(positions)
    new_cavity = _build(new_lsf.with_order, water)

    old_lsf = _GaussianDensity(positions)
    old_cavity = _build(old_lsf.legacy, water)

    assert new_cavity.ngrid == old_cavity.ngrid
    assert new_cavity.area == approx(old_cavity.area, rel=1.0e-13)
    assert new_cavity.volume == approx(old_cavity.volume, rel=1.0e-13)

    # The order-aware form must actually have been spared some work
    assert new_lsf.orders == {1, 2}
    assert old_lsf.orders == {3}


def test_callback_order_override(water) -> None:
    """`pass_order` overrides introspection in both directions."""
    _, positions = water

    # A two-argument callback forced into the legacy call form raises
    # TypeError inside the CFFI trampoline.
    lsf = _GaussianDensity(positions)
    with raises(TypeError):
        _build(lsf.with_order, water, pass_order=False)

    # Forcing the legacy form on a callback that really is legacy is a no-op
    lsf = _GaussianDensity(positions)
    cavity = _build(lsf.legacy, water, pass_order=False)
    assert cavity.ngrid > 0
    assert lsf.orders == {3}

    # Forcing the order form on a callback that accepts it works too
    lsf = _GaussianDensity(positions)
    cavity = _build(lsf.with_order, water, pass_order=True)
    assert cavity.ngrid > 0
    assert lsf.orders == {1, 2}


def test_callback_missing_derivative_is_reported(water) -> None:
    """A callback that omits a requested derivative must say so clearly."""
    _, positions = water
    lsf = _GaussianDensity(positions)

    def truncated(point, order):
        # Never returns a Hessian, whatever moist asks for
        value, grad = lsf._derivatives(point, 1)
        return value, grad

    with raises(ValueError, match="did not return the requested Hessian"):
        _build(truncated, water)


def test_callback_failure_aborts_build(water) -> None:
    """A raising callback must abort the build and surface the real exception.

    The failure is raised on the 50th evaluation rather than the first.
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

    # The exception must arrive with its own traceback.
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
    lsf = _GaussianDensity(positions)
    seen = []

    def flaky(point, order):
        seen.append(order)
        if len(seen) > 50:
            raise RuntimeError("no density here")
        return lsf._derivatives(point, order)

    with raises(RuntimeError, match="no density here"):
        _build(flaky, water)

    # The abort unwinds well before the whole grid is evaluated.
    assert 50 < len(seen) < 500


def test_callback_failure_is_not_sticky(water, flaky_lsf) -> None:
    """A cavity whose callback failed once must rebuild cleanly afterwards."""
    structure = Structure(*water)
    cavity = _isodensity(flaky_lsf.callback, flaky_lsf.rho_iso)
    flaky_lsf.arm()

    with raises(RuntimeError, match="no density here"):
        cavity.update(structure)

    flaky_lsf.arm(False)
    cavity.update(structure)
    assert cavity.cavity.ngrid > 0


def test_failed_rebuild_invalidates_previous_cavity_results(water, flaky_lsf) -> None:
    """A failed second build must not leave the first surface readable as current."""
    structure = Structure(*water)
    cavity = _isodensity(flaky_lsf.callback, flaky_lsf.rho_iso)
    cavity.update(structure)
    assert cavity.snapshot().ngrid > 0

    flaky_lsf.arm()
    with raises(RuntimeError, match="no density here"):
        cavity.update(structure)

    with raises(RuntimeError, match="successfully updated"):
        cavity.snapshot()


def test_failed_model_rebuild_invalidates_its_cavity_view(water, flaky_lsf) -> None:
    """Model updates propagate callback failures and invalidate their live view."""
    structure = Structure(*water)
    model = SolvationModel(
        CONTEXT, _isodensity(flaky_lsf.callback, flaky_lsf.rho_iso), [ModelComponentPV(1.0e-4)]
    )
    model.update(structure)
    assert model.cavity.snapshot().ngrid > 0

    flaky_lsf.arm()
    with raises(RuntimeError, match="no density here"):
        model.update(structure)

    with raises(RuntimeError, match="successfully updated"):
        model.cavity.snapshot()


@pytest.mark.parametrize("name,args,ctype,initial", [
    ("get_coupling_request_missing", (b"phi",), "bool", True),
])
def test_scalar_query_errors_preserve_outputs(name, args, ctype, initial):
    """C query failures use the error handle and never overwrite the result."""
    from moist.library import ffi, lib

    error = lib.moist_new_error()
    value = ffi.new(f"{ctype} *", initial)
    query = getattr(lib, "moist_" + name)
    try:
        query(error, ffi.NULL, *args, value)
        assert lib.moist_check_error(error) == lib.moist_failure
        assert value[0] == initial
        query(error, ffi.NULL, *args, ffi.NULL)
        assert lib.moist_check_error(error) == lib.moist_failure
        message = ffi.new("char[512]")
        size = ffi.new("int *", 512)
        lib.moist_get_error(error, message, size)
        assert b"Output pointer is missing" in ffi.string(message)
        query(ffi.NULL, ffi.NULL, *args, value)
        assert lib.moist_check_error(ffi.NULL) == lib.moist_invalid_error
        assert value[0] == initial
    finally:
        lib.moist_delete_error(ffi.new("moist_error *", error))


def test_next_coupling_request_fails_closed():
    """The cursor entry returns false on every failure.

    False is also the ordinary end of a pass; the wrapper checks the error
    handle after every call and raises instead.
    """
    from moist.library import CouplingHandle, next_coupling_request, ffi, lib

    error = lib.moist_new_error()
    try:
        assert not lib.moist_next_coupling_request(error, ffi.NULL)
        assert lib.moist_check_error(error) == lib.moist_failure
        assert not lib.moist_next_coupling_request(ffi.NULL, ffi.NULL)
    finally:
        lib.moist_delete_error(ffi.new("moist_error *", error))
    with raises(RuntimeError, match="next_coupling_request"):
        next_coupling_request(CouplingHandle.null())


def test_next_response_item_fails_closed():
    """The response cursor fails closed like the coupling's: false on every failure."""
    from moist.library import ResponseHandle, ffi, lib, next_response_item

    error = lib.moist_new_error()
    try:
        assert not lib.moist_next_response_item(error, ffi.NULL)
        assert lib.moist_check_error(error) == lib.moist_failure
        assert not lib.moist_next_response_item(ffi.NULL, ffi.NULL)
    finally:
        lib.moist_delete_error(ffi.new("moist_error *", error))
    with raises(RuntimeError, match="next_response_item"):
        next_response_item(ResponseHandle.null())


def test_callback_uninspectable_falls_back():
    class Uninspectable:
        @property
        def __signature__(self):
            raise ValueError("no signature")

        def __call__(self, point):
            raise AssertionError("introspection must not call the callback")

    assert not _callback_takes_order(Uninspectable())


@pytest.mark.parametrize("order", [1, 2, 3])
@pytest.mark.parametrize("attributes", [False, True])
def test_callback_transports_requested_density_derivatives(order, attributes):
    from moist.library import ffi

    point = np.array([0.2, -0.7, 1.3])
    gradient = np.array([2.0, -3.0, 5.0])
    hessian = np.arange(9, dtype=float).reshape(3, 3) + 7.0
    third = np.arange(27, dtype=float).reshape(3, 3, 3) - 11.0
    seen = []

    def callback(actual_point, actual_order):
        seen.append((actual_point.copy(), actual_order))
        values = (4.25, gradient, hessian, third)[:actual_order + 1]
        if attributes:
            return SimpleNamespace(**dict(zip(("rho", "drho", "d2rho", "d3rho"), values)))
        return values

    handle = _isodensity(callback)._as_handle()
    native_point = ffi.new("double[3]", point.tolist())
    rho = ffi.new("double *")
    drho = ffi.new("double[3]")
    d2rho = ffi.new("double[9]") if order >= 2 else ffi.NULL
    d3rho = ffi.new("double[27]") if order == 3 else ffi.NULL
    assert handle.callback_ref(ffi.NULL, native_point, rho, drho, d2rho, d3rho) == 0
    assert len(seen) == 1
    np.testing.assert_array_equal(seen[0][0], point)
    assert seen[0][1] == order
    assert rho[0] == 4.25
    np.testing.assert_array_equal(list(drho), gradient)
    if order >= 2:
        np.testing.assert_array_equal(list(d2rho), hessian.ravel(order="C"))
    if order == 3:
        np.testing.assert_array_equal(list(d3rho), third.ravel(order="C"))


@pytest.mark.parametrize("values,order,message", [
    ((1.0, np.zeros((1, 3))), 1, "gradient must have shape"),
    ((1.0, np.zeros(3), np.zeros(9)), 2, "Hessian must have shape"),
    ((1.0, np.zeros(3), np.zeros((3, 3)), np.zeros(27)), 3, "third derivative must have shape"),
    ((1.0, np.zeros(3)), 2, "did not return the requested Hessian"),
    ((1.0, np.zeros(3), np.zeros((3, 3))), 3, "did not return the requested third derivative"),
])
def test_callback_rejects_malformed_derivatives(values, order, message):
    from moist.library import ffi

    handle = _isodensity(lambda point, order: values)._as_handle()
    status = handle.callback_ref(
        ffi.NULL, ffi.new("double[3]"), ffi.new("double *"), ffi.new("double[3]"),
        ffi.new("double[9]") if order >= 2 else ffi.NULL,
        ffi.new("double[27]") if order == 3 else ffi.NULL,
    )
    assert status != 0
    with raises(ValueError, match=message):
        handle.callback_state.raise_if_failed()


def test_callback_preserves_base_exception():
    from moist.library import ffi

    failure = KeyboardInterrupt("density interrupted")

    def callback(point, order):
        raise failure

    handle = _isodensity(callback)._as_handle()
    assert handle.callback_ref(
        ffi.NULL, ffi.new("double[3]"), ffi.new("double *"), ffi.new("double[3]"),
        ffi.NULL, ffi.NULL,
    ) != 0
    with raises(KeyboardInterrupt) as caught:
        handle.callback_state.raise_if_failed()
    assert caught.value is failure


def test_callback_state_first_failure_consumption_and_reset():
    from moist.library import CallbackState

    state = CallbackState()
    first = RuntimeError("first")
    state.record(first)
    state.record(ValueError("later"))
    with raises(RuntimeError) as caught:
        state.raise_if_failed()
    assert caught.value is first
    state.raise_if_failed()
    state.record(ValueError("reset me"))
    state.reset()
    state.raise_if_failed()


@pytest.mark.parametrize("missing", [None, "cffi", "numpy", "pytest", "pyscf",
                                    "pyscf.dft", "pyscf.grad", "pyscf.gto", "pyscf.scf"])
def test_python_requirements_cli_checks_every_dependency(monkeypatch, capsys, missing):
    """The registered gate must fail for each unusable root or PySCF submodule."""
    names = ["cffi", "numpy", "pytest", "pyscf", "pyscf.dft", "pyscf.grad",
             "pyscf.gto", "pyscf.scf"]
    called = []

    def import_dependency(name):
        called.append(name)
        if name == missing:
            raise ImportError("dependency's internal import failed")
        return object()

    monkeypatch.setattr(importlib, "import_module", import_dependency)
    monkeypatch.setattr(sys, "argv", [str(REQUIREMENTS_SCRIPT), "cffi", "numpy", "pytest",
                                      "pyscf:dft,grad,gto,scf"])
    if missing:
        with raises(SystemExit) as caught:
            runpy.run_path(str(REQUIREMENTS_SCRIPT), run_name="__main__")
        message = caught.value.code
        assert isinstance(message, str)
        assert f"  {missing}: dependency's internal import failed" in message
        assert "Install them" in message and "-Dpython=false" in message
        assert "all Python test requirements are importable" not in capsys.readouterr().out
    else:
        runpy.run_path(str(REQUIREMENTS_SCRIPT), run_name="__main__")
        assert capsys.readouterr().out.strip() == "all Python test requirements are importable"
    expected = names[:4] if missing == "pyscf" else names
    assert set(called) == set(expected)


def test_python_requirements_aggregates_and_parses_specs(monkeypatch):
    main = runpy.run_path(str(REQUIREMENTS_SCRIPT))["main"]
    called = []

    def unavailable(name):
        called.append(name)
        if name in {"first", "second.a", "second.b"}:
            raise ImportError(f"broken {name}")
        return object()

    monkeypatch.setattr(importlib, "import_module", unavailable)
    with raises(SystemExit) as caught:
        main(["first:a,b", "second:a,,b,", "third:"])
    assert set(called) == {"first", "second", "second.a", "second.b", "third"}
    message = caught.value.code
    for name in ("first", "second.a", "second.b"):
        assert f"  {name}: broken {name}" in message


def test_python_requirements_fails_closed_on_unexpected_import_error(monkeypatch):
    main = runpy.run_path(str(REQUIREMENTS_SCRIPT))["main"]
    failure = RuntimeError("dependency initialization failed")

    def broken(name):
        raise failure

    monkeypatch.setattr(importlib, "import_module", broken)
    with raises((RuntimeError, SystemExit)) as caught:
        main(["numpy"])
    if isinstance(caught.value, SystemExit):
        assert caught.value.code not in (None, 0)
