"""Thin wrapper around the moist CFFI extension."""

import functools
import inspect
from typing import Optional

import numpy as np

try:
    from ._libmoist import ffi, lib
except ImportError as exc:
    raise ImportError("moist C extension unimportable, cannot use C-API") from exc


def get_api_version() -> str:
    """Return the current API version from moist."""
    api_version = lib.moist_get_version()
    return "{}.{}.{}".format(
        api_version // 10000,
        api_version % 10000 // 100,
        api_version % 100,
    )


class Handle:
    """Base wrapper for opaque C handles."""

    def __init__(self, handle):
        self.handle = handle

    @classmethod
    def with_gc(cls, handle):
        return cls(ffi.gc(handle, cls._delete))

    @classmethod
    def null(cls):
        return cls(ffi.NULL)

    @staticmethod
    def _delete(handle):
        raise NotImplementedError("Delete function not implemented")


class StructureHandle(Handle):
    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_structure *")
        ptr[0] = handle
        lib.moist_delete_structure(ptr)


class ModelHandle(Handle):
    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_model *")
        ptr[0] = handle
        lib.moist_delete_solvation_model(ptr)


class ComponentHandle(Handle):
    """Owning handle for a standalone general-model component."""

    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_component *")
        ptr[0] = handle
        lib.moist_delete_solvation_component(ptr)


class CavityHandle(Handle):
    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_cavity *")
        ptr[0] = handle
        lib.moist_delete_cavity(ptr)


def _delete_error(error):
    ptr = ffi.new("moist_error *")
    ptr[0] = error
    lib.moist_delete_error(ptr)


def new_error():
    return ffi.gc(lib.moist_new_error(), _delete_error)


def error_check(func):
    """Handle errors for moist library functions."""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        error = new_error()
        value = func(error, *args, **kwargs)
        if lib.moist_check_error(error):
            buffer_size = ffi.new("int *", 512)
            message = ffi.new("char[]", buffer_size[0])
            lib.moist_get_error(error, message, buffer_size)
            raise RuntimeError(ffi.string(message).decode())
        return value

    return handle_error


def new_structure(
    natoms: int,
    numbers: np.ndarray,
    positions: np.ndarray,
    lattice: Optional[np.ndarray],
    periodic: Optional[np.ndarray],
) -> StructureHandle:
    return StructureHandle.with_gc(
        error_check(lib.moist_new_structure)(
            natoms,
            _cast("int*", numbers),
            _cast("double*", positions),
            _cast("double*", lattice),
            _cast("bool*", periodic),
        )
    )


def update_structure(
    mol: StructureHandle,
    positions: np.ndarray,
    lattice: Optional[np.ndarray],
) -> None:
    return error_check(lib.moist_update_structure)(
        mol.handle,
        _cast("double*", positions),
        _cast("double*", lattice),
    )


def new_drop_cavity(
    nleb: Optional[int] = None,
    debug: bool = False,
    verbosity: int = 0,
    blend_k: Optional[float] = None,
    blend_1b: Optional[float] = None,
    blend_2b: Optional[float] = None,
    blend_3b: Optional[float] = None,
    do_fine: bool = False,
    tolerance: Optional[float] = None,
    proj_maxiter: Optional[int] = None,
    proj_level: Optional[int] = None,
    branch_weight_s: Optional[float] = None,
    rho_grid_h: Optional[float] = None,
    wleb_prune_level: Optional[int] = None,
) -> CavityHandle:
    """
    Create a standard solute-vdW (SvdW) DROP cavity with default CPCM radii.

    ``tolerance`` overrides the master numerical tolerance (``None`` keeps the
    compiled DROP default).
    """

    return CavityHandle.with_gc(
        error_check(lib.moist_new_drop_cavity)(
            _ref("int", nleb),
            _ref("bool", debug),
            _ref("int", verbosity),
            _ref("double", blend_k),
            _ref("double", blend_1b),
            _ref("double", blend_2b),
            _ref("double", blend_3b),
            _ref("bool", do_fine),
            _ref("double", tolerance),
            _ref("int", proj_maxiter),
            _ref("int", proj_level),
            _ref("double", branch_weight_s),
            _ref("double", rho_grid_h),
            _ref("int", wleb_prune_level),
        )
    )


def new_cfc_drop_cavity(
    nleb: Optional[int] = None,
    debug: bool = False,
    verbosity: int = 0,
    a1: Optional[float] = None,
    a2: Optional[float] = None,
    c: Optional[float] = None,
    m: Optional[int] = None,
    do_fine: bool = False,
    tolerance: Optional[float] = None,
    proj_maxiter: Optional[int] = None,
    proj_level: Optional[int] = None,
    branch_weight_s: Optional[float] = None,
    rho_grid_h: Optional[float] = None,
    wleb_prune_level: Optional[int] = None,
) -> CavityHandle:
    """Create a CFC-DROP cavity with default CPCM radii."""

    return CavityHandle.with_gc(
        error_check(lib.moist_new_cfc_drop_cavity)(
            _ref("int", nleb),
            _ref("bool", debug),
            _ref("int", verbosity),
            _ref("double", a1),
            _ref("double", a2),
            _ref("double", c),
            _ref("int", m),
            _ref("bool", do_fine),
            _ref("double", tolerance),
            _ref("int", proj_maxiter),
            _ref("int", proj_level),
            _ref("double", branch_weight_s),
            _ref("double", rho_grid_h),
            _ref("int", wleb_prune_level),
        )
    )


def new_iswig_cavity(
    nleb: Optional[int] = None,
    debug: bool = False,
    verbosity: int = 0,
    cut_a: Optional[float] = None,
    cut_f: Optional[float] = None,
) -> CavityHandle:
    """Create an iSwiG cavity with default CPCM radii."""

    return CavityHandle.with_gc(
        error_check(lib.moist_new_iswig_cavity)(
            _ref("int", nleb),
            _ref("bool", debug),
            _ref("int", verbosity),
            _ref("double", cut_a),
            _ref("double", cut_f),
        )
    )


class CallbackState:
    """Carries the *exception object* out of a CFFI callback frame.

    The callback ABI has a failure channel -- a nonzero return aborts the cavity
    build with a proper moist error -- but a return code only says *that* the
    callback failed. It cannot carry the Python exception, and CFFI will not let
    one cross the C frame either. So the wrapper returns nonzero to stop moist
    and stashes the exception here, and the wrapper that drove the C call
    re-raises it afterwards, with its original traceback, in place of moist's
    (correct but generic) "external LSF evaluation failed" error.
    """

    def __init__(self):
        self.exception: Optional[BaseException] = None

    def record(self, exc: BaseException) -> None:
        """Keep the first failure; later ones are usually knock-on effects."""
        if self.exception is None:
            self.exception = exc

    def reset(self) -> None:
        self.exception = None

    def raise_if_failed(self) -> None:
        """Re-raise a recorded callback failure, once."""
        exc = self.exception
        if exc is not None:
            self.exception = None
            raise exc


def _callback_takes_order(callback) -> bool:
    """Whether an isodensity callback accepts the derivative-order argument.

    Callbacks written before the order argument existed take ``point`` alone, so
    calling them with two arguments raises ``TypeError`` *inside* the CFFI
    trampoline, where it cannot become a moist error -- the Fortran caller would
    simply read unwritten buffers. Deciding the arity once, here, keeps that
    failure out of the hot path entirely.

    Anything introspection cannot resolve (builtins, ``*args``, C callables) is
    treated as the one-argument form: that form always works, it just forgoes
    the skip-computation speedup.
    """
    try:
        params = inspect.signature(callback).parameters.values()
    except (TypeError, ValueError):
        return False

    positional = 0
    for param in params:
        if param.kind in (param.POSITIONAL_ONLY, param.POSITIONAL_OR_KEYWORD):
            positional += 1
        elif param.kind is param.VAR_POSITIONAL:
            return False
    return positional >= 2


def new_drop_cavity_isodensity_callback(
    callback,
    nleb: Optional[int] = None,
    scale: float = 1000.0,
    debug: bool = False,
    verbosity: int = 0,
    do_fine: bool = False,
    wleb_prune_level: Optional[int] = None,
    tolerance: Optional[float] = None,
    pass_order: Optional[bool] = None,
) -> tuple[CavityHandle, object]:
    """Create a DROP cavity backed by a Python isodensity LSF callback.

    Two callback forms are accepted:

    * ``callback(point, order)`` -- ``point`` is a single point in Bohr and
      ``order`` is the highest derivative moist needs (1, 2 or 3). moist passes
      NULL Hessian/third-derivative buffers when it does not need them, so a
      callback in this form can skip *computing* the expensive high-order
      derivatives during the value+gradient-only projection phase.
    * ``callback(point)`` -- the original form. It always computes everything,
      which is correct but slower.

    The form is detected from the callback's signature; pass ``pass_order``
    explicitly to override that when introspection cannot decide (builtins,
    ``*args``, :func:`functools.partial` over a C callable).

    Either form must return ``(value, grad[, hess[, third]])`` or an object with
    ``value``, ``grad`` and optional ``hess``/``third`` attributes. The returned
    CFFI callback must be kept alive by the caller for at least as long as the
    cavity handle.

    ``tolerance`` overrides the master numerical tolerance (``None`` keeps the
    compiled DROP default).

    An exception raised by the callback aborts the cavity build: the wrapper
    reports failure through the callback's return code, which moist turns into a
    normal error on the update call. The exception object itself cannot cross
    the C frame, so it is recorded on ``handle.callback_state``; whoever drives
    the update calls :meth:`CallbackState.raise_if_failed` afterwards to re-raise
    the real cause with its own traceback.
    """

    if pass_order is None:
        pass_order = _callback_takes_order(callback)

    state = CallbackState()

    @ffi.callback("moist_isodensity_lsf_callback")
    def c_callback(context, point_ptr, value_ptr, grad_ptr, hess_ptr, third_ptr):
        want_hess = hess_ptr != ffi.NULL
        want_third = third_ptr != ffi.NULL
        try:
            point = np.frombuffer(ffi.buffer(point_ptr, 24), dtype=np.float64)
            if pass_order:
                order = 3 if want_third else (2 if want_hess else 1)
                result = callback(point, order)
            else:
                result = callback(point)

            if hasattr(result, "value"):
                value = result.value
                grad = result.grad
                hess = getattr(result, "hess", None)
                third = getattr(result, "third", None)
            else:
                value, grad = result[0], result[1]
                hess = result[2] if len(result) > 2 else None
                third = result[3] if len(result) > 3 else None

            grad = np.asarray(grad, dtype=np.float64)
            if grad.shape != (3,):
                raise ValueError("Isodensity callback gradient must have shape (3,)")

            value_ptr[0] = float(value)
            np.frombuffer(ffi.buffer(grad_ptr, 24), dtype=np.float64)[:] = grad
            if want_hess:
                if hess is None:
                    raise ValueError("Isodensity callback did not return the requested Hessian")
                hess = np.asarray(hess, dtype=np.float64)
                if hess.shape != (3, 3):
                    raise ValueError("Isodensity callback Hessian must have shape (3, 3)")
                np.frombuffer(ffi.buffer(hess_ptr, 72), dtype=np.float64)[:] = hess.ravel(order="F")
            if want_third:
                if third is None:
                    raise ValueError(
                        "Isodensity callback did not return the requested third derivative"
                    )
                third = np.asarray(third, dtype=np.float64)
                if third.shape != (3, 3, 3):
                    raise ValueError(
                        "Isodensity callback third derivative must have shape (3, 3, 3)"
                    )
                np.frombuffer(ffi.buffer(third_ptr, 216), dtype=np.float64)[:] = third.ravel(
                    order="F"
                )
        except BaseException as exc:
            # CFFI cannot propagate an exception through the C frame, but the
            # callback ABI has a failure channel: returning nonzero aborts the
            # build with a moist error. Keep the exception so the caller can
            # re-raise the real cause (see CallbackState) and leave the output
            # buffers alone -- moist ignores them once the status is nonzero.
            state.record(exc)
            return 1
        return 0

    handle = CavityHandle.with_gc(
        error_check(lib.moist_new_drop_cavity_isodensity_callback)(
            c_callback,
            ffi.NULL,
            _ref("double", scale),
            _ref("int", nleb),
            _ref("bool", debug),
            _ref("int", verbosity),
            _ref("bool", do_fine),
            _ref("int", wleb_prune_level),
            _ref("double", tolerance),
        )
    )
    #: Attached rather than returned so the (handle, callback) result stays a
    #: two-tuple; callers that drive an update read it back off the handle.
    handle.callback_state = state
    return handle, c_callback


def update_model(model: ModelHandle, structure: StructureHandle) -> None:
    return error_check(lib.moist_update_solvation_model)(
        model.handle,
        structure.handle,
    )


def get_model_energy(model: ModelHandle) -> float:
    energy = np.array(0.0, dtype=np.float64)
    error_check(lib.moist_get_solvation_model_energy)(
        model.handle,
        _cast("double*", energy),
    )
    return float(energy)


def get_model_cavity(model: ModelHandle) -> CavityHandle:
    return CavityHandle.with_gc(
        error_check(lib.moist_get_solvation_model_cavity)(model.handle)
    )


def new_cpcm_component(epsilon: float, solver: int) -> ComponentHandle:
    """Create a CPCM component for a general solvation model."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_cpcm_component)(float(epsilon), int(solver))
    )


def new_cosmo_component(epsilon: float, solver: int) -> ComponentHandle:
    """Create a COSMO component for a general solvation model."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_cosmo_component)(float(epsilon), int(solver))
    )


def new_pv_component(pressure: float) -> ComponentHandle:
    """Create a pressure-volume component whose energy is pressure times volume."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_pv_component)(float(pressure))
    )


def new_gostshyp_component(pressure: float) -> ComponentHandle:
    """Create a GOSTSHYP hydrostatic-pressure component.

    ``pressure`` is in Hartree/bohr^3.  The component needs Gaussian density
    moments supplied through :func:`general_model_supply_gostshyp` after every
    cavity update; it cannot form them itself.
    """

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_gostshyp_component)(float(pressure))
    )


def new_general_model(
    cavity: CavityHandle,
    components: list[ComponentHandle],
    debug: bool = False,
    verbosity: int = 0,
) -> ModelHandle:
    """Create a general model and append copies of the requested components."""

    model = ModelHandle.with_gc(
        error_check(lib.moist_new_general_solvation_model)(
            cavity.handle,
            bool(debug),
            int(verbosity),
        )
    )
    for component in components:
        error_check(lib.moist_general_model_add_component)(
            model.handle, component.handle
        )
    return model


def general_model_supply_electrostatics(
    model: ModelHandle,
    phi: np.ndarray,
    w_xi: Optional[np.ndarray] = None,
    w_f: Optional[np.ndarray] = None,
    w_xyz: Optional[np.ndarray] = None,
    w_n: Optional[np.ndarray] = None,
    qefield: Optional[np.ndarray] = None,
) -> None:
    """Supply CPCM traces and optional direct geometry-response arrays."""

    _phi = np.ascontiguousarray(phi, dtype=np.float64).reshape(-1)
    ngrid = int(_phi.size)

    def scalar(name, value):
        if value is None:
            return None
        array = np.ascontiguousarray(value, dtype=np.float64)
        if array.shape != (ngrid,):
            raise ValueError(f"{name} must have shape (ngrid,)")
        return array

    def vector(name, value):
        if value is None:
            return None
        array = np.asarray(value, dtype=np.float64, order="F")
        if array.shape != (3, ngrid):
            raise ValueError(f"{name} must have shape (3, ngrid)")
        return array

    _w_xi = scalar("w_xi", w_xi)
    _w_f = scalar("w_f", w_f)
    _w_xyz = vector("w_xyz", w_xyz)
    _w_n = vector("w_n", w_n)
    _qefield = vector("qefield", qefield)
    error_check(lib.moist_general_model_supply_electrostatics)(
        model.handle,
        ngrid,
        _cast("double*", _phi),
        _cast("double*", _w_xi),
        _cast("double*", _w_f),
        _cast("double*", _w_xyz),
        _cast("double*", _w_n),
        _cast("double*", _qefield),
    )


def general_model_supply_gostshyp(
    model: ModelHandle,
    gt: np.ndarray,
    pt: np.ndarray,
    mt: np.ndarray,
    rt: np.ndarray,
) -> None:
    """Supply the Gaussian density moments the GOSTSHYP component consumes.

    Moments of the solute density against the unnormalized Gaussian
    ``exp(-w_i |r - r_i|^2)`` on each grid point, in native cavity order:
    ``gt = <G>``, ``pt = <(r-r_i) G>``, ``mt = <(r-r_i)(r-r_i) G>`` and
    ``rt = <(r-r_i) |r-r_i|^2 G>``.  All four are required.
    """

    _gt = np.ascontiguousarray(gt, dtype=np.float64).reshape(-1)
    ngrid = int(_gt.size)

    def vector(name, value):
        array = np.asarray(value, dtype=np.float64, order="F")
        if array.shape != (3, ngrid):
            raise ValueError(f"{name} must have shape (3, ngrid)")
        return array

    _pt = vector("pt", pt)
    _rt = vector("rt", rt)
    _mt = np.asarray(mt, dtype=np.float64, order="F")
    if _mt.shape != (3, 3, ngrid):
        raise ValueError("mt must have shape (3, 3, ngrid)")

    error_check(lib.moist_general_model_supply_gostshyp)(
        model.handle,
        ngrid,
        _cast("double*", _gt),
        _cast("double*", _pt),
        _cast("double*", _mt),
        _cast("double*", _rt),
    )


def general_model_get_trace_response(
    model: ModelHandle, ngrid: int
) -> np.ndarray:
    """Return the accumulated surface charges from all model components.

    The surface charge ``q_i`` equals ``dE/dphi_i`` by stationarity; the host
    contracts it as ``F += sum_i q_i V(r_i)``.  Requesting it from a model that
    produces no surface charges raises rather than returning zeros.
    """

    surface_charge = np.zeros(ngrid, dtype=np.float64)
    error_check(lib.moist_general_model_get_trace_response)(
        model.handle,
        int(ngrid),
        _cast("double*", surface_charge),
    )
    return surface_charge


def general_model_get_response(
    model: ModelHandle,
    ngrid: int,
    *,
    electrostatics: bool = True,
    lsf: bool = True,
    gostshyp: bool = False,
) -> dict:
    """Return the requested response channels from a general model.

    Each keyword selects one channel group.  A group that is *not* requested is
    skipped entirely (a NULL pointer at the C boundary) and comes back ``None``.
    A group that *is* requested must be produced by the model configuration, or
    the call raises: absence is a legitimate physical answer here -- a cavity
    with field-independent geometry has no level-set response -- so zeros could
    not be told apart from a genuine result and are never returned silently.

    Requesting ``gostshyp`` uses the extended entry point, which reports the
    amplitudes conjugate to the host's Gaussian integral blocks; the host
    completes its Fock contribution as
    ``F += sum_i [w_overlap[i] g[..., i] + w_normal_deriv[i] f[..., i]]``.

    Prefer one call with every group you need: assembling a response contracts
    the cavity surface adjoints once, and splitting the read pays that twice.
    """

    surface_charge = np.zeros(ngrid, dtype=np.float64) if electrostatics else None
    w_value = np.zeros(ngrid, dtype=np.float64) if lsf else None
    w_gradient = (
        np.zeros((3, ngrid), dtype=np.float64, order="F") if lsf else None
    )
    w_hessian = (
        np.zeros((3, 3, ngrid), dtype=np.float64, order="F") if lsf else None
    )
    w_overlap = np.zeros(ngrid, dtype=np.float64) if gostshyp else None
    w_normal_deriv = np.zeros(ngrid, dtype=np.float64) if gostshyp else None

    if gostshyp:
        error_check(lib.moist_general_model_get_response_extended)(
            model.handle,
            int(ngrid),
            _cast("double*", surface_charge),
            _cast("double*", w_value),
            _cast("double*", w_gradient),
            _cast("double*", w_hessian),
            _cast("double*", w_overlap),
            _cast("double*", w_normal_deriv),
        )
    else:
        error_check(lib.moist_general_model_get_response)(
            model.handle,
            int(ngrid),
            _cast("double*", surface_charge),
            _cast("double*", w_value),
            _cast("double*", w_gradient),
            _cast("double*", w_hessian),
        )

    return {
        "electrostatics": (
            {"surface_charge": surface_charge} if electrostatics else None
        ),
        "lsf": (
            {
                "w_value": w_value,
                "w_gradient": w_gradient,
                "w_hessian": w_hessian,
            }
            if lsf
            else None
        ),
        "gostshyp": (
            {"w_overlap": w_overlap, "w_normal_deriv": w_normal_deriv}
            if gostshyp
            else None
        ),
    }


def general_model_get_gradient(model: ModelHandle, natoms: int) -> np.ndarray:
    """Return the accumulated nuclear gradient with shape (3, natoms)."""

    gradient = np.zeros((3, natoms), dtype=np.float64, order="F")
    error_check(lib.moist_general_model_get_gradient)(
        model.handle,
        int(natoms),
        _cast("double*", gradient),
    )
    return gradient


def update_cavity(cavity: CavityHandle, structure: StructureHandle) -> None:
    return error_check(lib.moist_update_cavity)(
        cavity.handle,
        structure.handle,
    )


def get_cavity_sizes(cavity: CavityHandle) -> tuple[int, int]:
    ngrid = ffi.new("int *")
    nsph = ffi.new("int *")
    error_check(lib.moist_get_cavity_sizes)(cavity.handle, ngrid, nsph)
    return int(ngrid[0]), int(nsph[0])


def get_cavity_results(cavity: CavityHandle) -> dict:
    """Return the generic cavity results.

    The buffers are allocated from the cavity's current sizes and those sizes
    are handed to the C entry point as the array capacities, so a cavity
    rebuilt between the two calls raises a clean API error instead of writing
    past the buffers.
    """

    ngrid, nsph = get_cavity_sizes(cavity)

    area = np.array(0.0, dtype=np.float64)
    volume = np.array(0.0, dtype=np.float64)
    out_ngrid = ffi.new("int *")
    out_nsph = ffi.new("int *")
    xyz = np.zeros((3, ngrid), dtype=np.float64, order="F")
    weights = np.zeros(ngrid, dtype=np.float64)
    owner = np.zeros(ngrid, dtype=np.int32)
    converged = np.zeros(ngrid, dtype=np.bool_)
    radii = np.zeros(nsph, dtype=np.float64)
    asph = np.zeros(nsph, dtype=np.float64)

    error_check(lib.moist_get_cavity_results)(
        cavity.handle,
        ngrid,
        nsph,
        _cast("double*", area),
        _cast("double*", volume),
        out_ngrid,
        out_nsph,
        _cast("double*", xyz),
        _cast("double*", weights),
        _cast("int*", owner),
        _cast("bool*", converged),
        _cast("double*", radii),
        _cast("double*", asph),
    )

    return {
        "area": float(area),
        "volume": float(volume),
        "ngrid": int(out_ngrid[0]),
        "nsph": int(out_nsph[0]),
        "xyz": xyz,
        "a": weights,
        "owner": owner,
        "converged": converged,
        "radii": radii,
        "asph": asph,
    }


def get_drop_specific(cavity: CavityHandle, ngrid: Optional[int] = None) -> dict:
    """Return the DROP-specific cavity fields.

    ``ngrid`` sizes the buffers and is passed on as the array capacity; it has
    to be at least the cavity's own ngrid, otherwise the C entry point reports
    an error and writes nothing. Pass ``None`` to read the size from the cavity.
    """

    if ngrid is None:
        ngrid, _ = get_cavity_sizes(cavity)

    nmax = ffi.new("int *")
    normal0 = np.zeros((3, ngrid), dtype=np.float64, order="F")
    wleb = np.zeros(ngrid, dtype=np.float64)
    r_iI0 = np.zeros(ngrid, dtype=np.float64)
    switch_f = np.zeros(ngrid, dtype=np.float64)
    rho = np.zeros(ngrid, dtype=np.float64)

    error_check(lib.moist_get_drop_specific)(
        cavity.handle,
        ngrid,
        nmax,
        _cast("double*", normal0),
        _cast("double*", wleb),
        _cast("double*", r_iI0),
        _cast("double*", switch_f),
        _cast("double*", rho),
    )

    return {
        "nmax": int(nmax[0]),
        "normal0": normal0,
        "wleb": wleb,
        "r_iI0": r_iI0,
        "f": switch_f,
        "rho": rho,
    }


def get_drop_numbering(cavity: CavityHandle, ngrid: Optional[int] = None) -> np.ndarray:
    """Return the stable per-grid-point numbering of a DROP cavity.

    The numbering packs the anchor and the branch index into one integer as
    ``anchor_id + nmax*(branch - 1)``, so it identifies the same physical
    surface point across rebuilds and lets callers recover which points are
    branches of a common anchor.

    ``ngrid`` sizes the buffer and is passed on as the array capacity; pass
    ``None`` to read the size from the cavity.
    """

    if ngrid is None:
        ngrid, _ = get_cavity_sizes(cavity)

    numbering = np.zeros(ngrid, dtype=np.int32)

    error_check(lib.moist_get_drop_numbering)(
        cavity.handle,
        ngrid,
        _cast("int*", numbering),
    )

    return numbering


def assemble_drop_amat(cavity: CavityHandle) -> tuple[np.ndarray, np.ndarray]:
    """Assemble the Gaussian CPCM A-matrix and return it with xi values.

    No longer DROP-specific: the underlying C entry point now works for every
    Gaussian-discretized cavity.
    """

    ngrid, _ = get_cavity_sizes(cavity)
    amat = np.zeros((ngrid, ngrid), dtype=np.float64, order="F")
    xi = np.zeros(ngrid, dtype=np.float64)

    error_check(lib.moist_assemble_amat)(
        cavity.handle,
        ngrid,
        _cast("double*", amat),
        _cast("double*", xi),
    )
    return amat, xi


def get_cavity_gaussian(cavity: CavityHandle) -> tuple[np.ndarray, np.ndarray]:
    """Return Gaussian widths and switching factors in native cavity order."""

    ngrid, _ = get_cavity_sizes(cavity)
    xi = np.zeros(ngrid, dtype=np.float64)
    switch = np.zeros(ngrid, dtype=np.float64)
    error_check(lib.moist_get_cavity_gaussian)(
        cavity.handle,
        ngrid,
        _cast("double*", xi),
        _cast("double*", switch),
    )
    return xi, switch


def compute_anchor_gradient(cavity: CavityHandle) -> None:
    """Compute the anchor-only nuclear derivatives of a DROP cavity.

    Restricts every grid point's nuclear coupling to its owner atom's rigid
    anchor motion. For a callback level set the field's own nuclear partials
    are zero, so this is the entire nuclear route moist can see; the host adds
    the field route by contracting the LSF adjoints with its own ``dS/dR``.
    """

    error_check(lib.moist_compute_anchor_gradient)(cavity.handle)


def get_anchor_gradient(cavity: CavityHandle) -> dict:
    """Return the anchor-channel nuclear derivatives in native cavity order.

    Requires a preceding :func:`compute_anchor_gradient`. The buffers are sized
    from the cavity's current sizes and those same sizes are handed over as the
    array capacities, so a cavity rebuilt in between fails with a clean API
    error instead of writing past the buffers.

    The per-point area derivative ``a_i1_rA`` is the one a geometric surface
    functional needs: the grid point area carries a switching-function dependence
    (``a_i ~ f_i / xi_i**2``), so it is not recoverable from ``xi1_rA`` alone.
    """

    ngrid, nsph = get_cavity_sizes(cavity)

    xyz1_rA = np.zeros((3, 3, nsph, ngrid), dtype=np.float64, order="F")
    xi1_rA = np.zeros((3, nsph, ngrid), dtype=np.float64, order="F")
    a_i1_rA = np.zeros((3, nsph, ngrid), dtype=np.float64, order="F")
    v_i1_rA = np.zeros((3, nsph, ngrid), dtype=np.float64, order="F")
    A_tot1_rA = np.zeros((3, nsph), dtype=np.float64, order="F")
    V_tot1_rA = np.zeros((3, nsph), dtype=np.float64, order="F")

    # The capacity pair is (nsph, ngrid) -- the reverse of the order the grid
    # index appears in the array shapes above.
    error_check(lib.moist_get_anchor_gradient)(
        cavity.handle,
        int(nsph),
        int(ngrid),
        _cast("double*", xyz1_rA),
        _cast("double*", xi1_rA),
        _cast("double*", a_i1_rA),
        _cast("double*", v_i1_rA),
        _cast("double*", A_tot1_rA),
        _cast("double*", V_tot1_rA),
    )
    return {
        "xyz1_rA": xyz1_rA,
        "xi1_rA": xi1_rA,
        "a_i1_rA": a_i1_rA,
        "v_i1_rA": v_i1_rA,
        "A_tot1_rA": A_tot1_rA,
        "V_tot1_rA": V_tot1_rA,
    }


def contract_amat1_q1q2_surface_weights(
    cavity: CavityHandle,
    q1: np.ndarray,
    q2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contract Gaussian PCM A-matrix derivatives to per-grid surface weights."""

    ngrid, _ = get_cavity_sizes(cavity)
    _q1 = np.ascontiguousarray(q1, dtype=np.float64)
    _q2 = np.ascontiguousarray(q2, dtype=np.float64)
    if _q1.shape != (ngrid,) or _q2.shape != (ngrid,):
        raise ValueError("q1 and q2 must have shape (ngrid,)")

    w_xi = np.zeros(ngrid, dtype=np.float64)
    w_f = np.zeros(ngrid, dtype=np.float64)
    w_xyz = np.zeros((3, ngrid), dtype=np.float64, order="F")

    error_check(lib.moist_contract_amat1_q1q2_surface_weights)(
        cavity.handle,
        _cast("double*", _q1),
        _cast("double*", _q2),
        _cast("double*", w_xi),
        _cast("double*", w_f),
        _cast("double*", w_xyz),
    )
    return w_xi, w_f, w_xyz


def contract_surface_lsf_weights(
    cavity: CavityHandle,
    w_xi: np.ndarray,
    w_f: np.ndarray,
    w_xyz: np.ndarray,
    w_n: Optional[np.ndarray] = None,
    w_k1: Optional[np.ndarray] = None,
    w_k2: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contract DROP surface weights to LSF adjoint weights.

    The outward-normal (``w_n``, shape ``(3, ngrid)``) and principal-curvature
    (``w_k1``/``w_k2``, shape ``(ngrid,)``) channels are optional; ``None``
    skips the channel entirely.
    """

    ngrid, _ = get_cavity_sizes(cavity)
    _w_xi = np.ascontiguousarray(w_xi, dtype=np.float64)
    _w_f = np.ascontiguousarray(w_f, dtype=np.float64)
    _w_xyz = np.asarray(w_xyz, dtype=np.float64, order="F")
    if _w_xi.shape != (ngrid,) or _w_f.shape != (ngrid,) or _w_xyz.shape != (3, ngrid):
        raise ValueError("w_xi/w_f must have shape (ngrid,), w_xyz must have shape (3, ngrid)")

    _w_n = None if w_n is None else np.asarray(w_n, dtype=np.float64, order="F")
    _w_k1 = None if w_k1 is None else np.ascontiguousarray(w_k1, dtype=np.float64)
    _w_k2 = None if w_k2 is None else np.ascontiguousarray(w_k2, dtype=np.float64)
    if _w_n is not None and _w_n.shape != (3, ngrid):
        raise ValueError("w_n must have shape (3, ngrid)")
    if _w_k1 is not None and _w_k1.shape != (ngrid,):
        raise ValueError("w_k1 must have shape (ngrid,)")
    if _w_k2 is not None and _w_k2.shape != (ngrid,):
        raise ValueError("w_k2 must have shape (ngrid,)")

    w_lsf0 = np.zeros(ngrid, dtype=np.float64)
    w_lsf1 = np.zeros((3, ngrid), dtype=np.float64, order="F")
    w_lsf2 = np.zeros((3, 3, ngrid), dtype=np.float64, order="F")

    error_check(lib.moist_contract_surface_lsf_weights_extended)(
        cavity.handle,
        _cast("double*", _w_xi),
        _cast("double*", _w_f),
        _cast("double*", _w_xyz),
        _cast("double*", w_lsf0),
        _cast("double*", w_lsf1),
        _cast("double*", w_lsf2),
        _cast("double*", _w_n),
        _cast("double*", _w_k1),
        _cast("double*", _w_k2),
    )
    return w_lsf0, w_lsf1, w_lsf2


def _char(value: Optional[str]):
    return ffi.new("char[]", value.encode()) if value is not None else ffi.NULL


def _ref(ctype: str, value):
    if value is None:
        return ffi.NULL
    ref = ffi.new(ctype + " *")
    ref[0] = value
    return ref


def _cast(ctype: str, array):
    return ffi.cast(ctype, array.ctypes.data) if array is not None else ffi.NULL
