"""Thin wrapper around the moist CFFI extension."""

import functools
import inspect
from dataclasses import dataclass
from typing import Iterable, Optional

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


def _native_text(entry, *args):
    length = ffi.new("size_t *")
    call = error_check(entry)
    call(*args, ffi.NULL, 0, length)
    buffer = ffi.new("char[]", length[0] + 1)
    call(*args, buffer, length[0] + 1, length)
    return ffi.string(buffer).decode("utf-8")


def get_version_string():
    """Full native version, including prerelease information."""
    return _native_text(lib.moist_get_version_string)


def get_banner(style="full"):
    """Return native banner text; style is full, short, ascii or build."""
    if style not in ("full", "short", "ascii", "build"):
        raise ValueError("banner style must be full, short, ascii or build")
    return _native_text(lib.moist_get_banner, getattr(lib, "moist_banner_" + style))


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
        lib.moist_delete_model(ptr)


class ComponentHandle(Handle):
    """Owning handle for a standalone general-model component."""

    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_component *")
        ptr[0] = handle
        lib.moist_delete_component(ptr)


class CavityHandle(Handle):
    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_cavity *")
        ptr[0] = handle
        lib.moist_delete_cavity(ptr)


class RadiiHandle(Handle):
    @staticmethod
    def _delete(handle):
        lib.moist_delete_radii(ffi.new("moist_radii *", handle))


class CouplingHandle(Handle):
    """Owning handle for a model's host coupling: requests, answers and cursor."""

    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_coupling *")
        ptr[0] = handle
        lib.moist_delete_coupling(ptr)


class ResponseHandle(Handle):
    """Owning handle for the item list a phase hands back to the host."""

    @staticmethod
    def _delete(handle):
        ptr = ffi.new("moist_response *")
        ptr[0] = handle
        lib.moist_delete_response(ptr)


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


def _options(kind, **values):
    """Initialize a typed options struct and apply explicit overrides."""
    options = ffi.new(f"moist_{kind}_options *")
    error_check(getattr(lib, f"moist_init_{kind}_options"))(
        options, ffi.sizeof(options[0])
    )
    for name, value in values.items():
        if value is not None:
            setattr(options, name, value)
    return options


def _drop_from_lsf(lsf, options, radii=None):
    """Copy an LSF into a cavity and release the temporary LSF handle."""
    try:
        return CavityHandle.with_gc(
            error_check(lib.moist_new_drop_cavity)(
                lsf, ffi.NULL if radii is None else radii.handle, options
            )
        )
    finally:
        ptr = ffi.new("moist_lsf *", lsf)
        lib.moist_delete_lsf(ptr)


def new_drop_cavity(parameters, lsf_parameters, radii=None) -> CavityHandle:
    """Copy a geometric LSF and radii into a DROP cavity."""
    options = parameters._as_options()
    lsf_options = lsf_parameters._as_options()
    constructor = getattr(lib, f"moist_new_{lsf_parameters._kind}_lsf")
    return _drop_from_lsf(error_check(constructor)(lsf_options), options, radii)


def new_internal_isodensity_cavity(source, parameters, lsf_parameters, radii=None):
    basis = source.basis
    arrays = [np.asarray(getattr(basis, name), dtype=dtype) for name, dtype in (
        ("shell_atom", np.int32), ("shell_l", np.int32), ("shell_nprim", np.int32),
        ("exponents", np.float64), ("coefficients", np.float64))]
    lsf = error_check(lib.moist_new_isodensity_lsf)(
        len(basis.shell_atom), *[_cast("int*" if i < 3 else "double*", a)
                                 for i, a in enumerate(arrays)], lsf_parameters._as_options())
    cavity = _drop_from_lsf(lsf, parameters._as_options(), radii)
    set_isodensity_density(cavity, source.density_matrix)
    return cavity


def isodensity_layout(cavity):
    ncart, nshell = ffi.new("int *"), ffi.new("int *")
    entry = error_check(lib.moist_get_isodensity_cart_layout)
    entry(cavity.handle, ncart, nshell, ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL)
    offsets = np.empty(nshell[0] + 1, dtype=np.int32)
    powers = [np.empty(ncart[0], dtype=np.int32) for _ in range(3)]
    entry(cavity.handle, ncart, nshell, _cast("int*", offsets),
          *[_cast("int*", a) for a in powers])
    return offsets, np.stack(powers, axis=1)


def set_isodensity_density(handle, density, *, model=False):
    array = np.asarray(density, dtype=np.float64, order="C")
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("density must be a square matrix")
    entry = lib.moist_set_model_isodensity_density if model else lib.moist_set_isodensity_density
    error_check(entry)(handle.handle, len(array), _cast("double*", array))


def new_iswig_cavity(parameters, radii=None) -> CavityHandle:
    return CavityHandle.with_gc(
        error_check(lib.moist_new_iswig_cavity)(
            ffi.NULL if radii is None else radii.handle, parameters._as_options()
        )
    )


def new_radii(kind: str) -> RadiiHandle:
    return RadiiHandle.with_gc(error_check(getattr(lib, f"moist_new_{kind}_radii"))())


def set_custom_radii_atoms(radii, values):
    error_check(lib.moist_set_custom_radii_atoms)(
        radii.handle, len(values), _cast("double*", values)
    )


def set_custom_radii_elements(radii, numbers, values):
    error_check(lib.moist_set_custom_radii_elements)(
        radii.handle, len(values), _cast("int*", numbers), _cast("double*", values)
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
    callback, parameters, lsf_parameters, radii=None, *, pass_order=None,
) -> CavityHandle:
    """Build a callback cavity, retaining its callback and exception state.

    Callbacks accept ``(point, order)`` or ``(point)`` and return bare density
    derivatives. The highest requested order is 1, 2 or 3; optional native
    buffers are only written when requested. ``pass_order`` overrides signature
    detection. Exceptions are retained for the high-level update guard.
    """

    if pass_order is None:
        pass_order = _callback_takes_order(callback)

    state = CallbackState()

    @ffi.callback("moist_isodensity_lsf_callback")
    def c_callback(context, point_ptr, rho_ptr, drho_ptr, d2rho_ptr, d3rho_ptr):
        want_hess = d2rho_ptr != ffi.NULL
        want_third = d3rho_ptr != ffi.NULL
        try:
            point = np.frombuffer(ffi.buffer(point_ptr, 24), dtype=np.float64)
            if pass_order:
                order = 3 if want_third else (2 if want_hess else 1)
                result = callback(point, order)
            else:
                result = callback(point)

            if hasattr(result, "rho"):
                rho = result.rho
                drho = result.drho
                d2rho = getattr(result, "d2rho", None)
                d3rho = getattr(result, "d3rho", None)
            else:
                rho, drho = result[0], result[1]
                d2rho = result[2] if len(result) > 2 else None
                d3rho = result[3] if len(result) > 3 else None

            drho = np.asarray(drho, dtype=np.float64)
            if drho.shape != (3,):
                raise ValueError("Isodensity callback density gradient must have shape (3,)")

            rho_ptr[0] = float(rho)
            np.frombuffer(ffi.buffer(drho_ptr, 24), dtype=np.float64)[:] = drho
            if want_hess:
                if d2rho is None:
                    raise ValueError("Isodensity callback did not return the requested Hessian")
                d2rho = np.asarray(d2rho, dtype=np.float64)
                if d2rho.shape != (3, 3):
                    raise ValueError("Isodensity callback density Hessian must have shape (3, 3)")
                np.frombuffer(ffi.buffer(d2rho_ptr, 72), dtype=np.float64)[:] = d2rho.ravel(
                    order="C"
                )
            if want_third:
                if d3rho is None:
                    raise ValueError(
                        "Isodensity callback did not return the requested third derivative"
                    )
                d3rho = np.asarray(d3rho, dtype=np.float64)
                if d3rho.shape != (3, 3, 3):
                    raise ValueError(
                        "Isodensity callback third derivative must have shape (3, 3, 3)"
                    )
                np.frombuffer(ffi.buffer(d3rho_ptr, 216), dtype=np.float64)[:] = d3rho.ravel(
                    order="C"
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

    options = parameters._as_options()
    lsf_options = lsf_parameters._as_options()
    lsf = error_check(lib.moist_new_isodensity_callback_lsf)(
        c_callback, ffi.NULL, lsf_options
    )
    handle = _drop_from_lsf(lsf, options, radii)
    handle.callback_state = state
    handle.callback_ref = c_callback
    return handle


def update_model(model: ModelHandle, structure: StructureHandle) -> None:
    return error_check(lib.moist_update_model)(
        model.handle,
        structure.handle,
    )


def get_model_cavity(model: ModelHandle) -> CavityHandle:
    handle = CavityHandle.with_gc(
        error_check(lib.moist_get_model_cavity)(model.handle)
    )
    handle._owner = model
    return handle


def new_cpcm_component(epsilon: float, parameters) -> ComponentHandle:
    """Create a CPCM component for a general solvation model."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_cpcm_component)(float(epsilon), parameters._as_options())
    )


def new_cosmo_component(epsilon: float, parameters) -> ComponentHandle:
    """Create a COSMO component for a general solvation model."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_cosmo_component)(float(epsilon), parameters._as_options())
    )


def new_pv_component(pressure: float) -> ComponentHandle:
    """Create a pressure-volume component whose energy is pressure times volume."""

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_pv_component)(float(pressure))
    )


def new_gostshyp_component(pressure: float) -> ComponentHandle:
    """Create a GOSTSHYP hydrostatic-pressure component.

    ``pressure`` is in Hartree/bohr^3.  The component declares a Gaussian
    moment request on every coupling and cannot form the moments itself; the
    host answers it with :func:`answer_coupling_request` and reads the amplitudes back
    with :func:`get_response_array`.
    """

    return ComponentHandle.with_gc(
        error_check(lib.moist_new_gostshyp_component)(float(pressure))
    )


def new_general_model(
    cavity: CavityHandle,
    components: list[ComponentHandle],
    parameters,
) -> ModelHandle:
    """Create a general model and append copies of the requested components."""

    model = ModelHandle.with_gc(
        error_check(lib.moist_new_model)(
            cavity.handle,
            parameters._as_options(),
        )
    )
    # Native copies borrow Python callbacks. Retain their owning handle even
    # when callers use this low-level constructor without a SolvationModel.
    model._source_cavity = cavity
    for component in components:
        error_check(lib.moist_add_model_component)(
            model.handle, component.handle
        )
    return model


# -----------------------------------------------------------------------------
# Host coupling protocol
# -----------------------------------------------------------------------------
#
# A coupling carries the model's requests and a cursor over them; one response
# handle receives what a phase hands back and carries a cursor over its items.
# next_coupling_request() and next_response_item() move the cursors, and the
# other request and response functions act on the request or item they stopped
# at, failing by name when none is current.  Requests and response items are
# identified by their scientific names, never by an index or a token.  Grid
# inputs are read from the model's cavity.  Like the C entries they bind, these
# functions take no sizes: moist reads and writes exactly the documented
# shape, so an array passed in has to have it.  Only the element type and the
# layout are checked here, which the pointer cast relies on.


def new_coupling(model: ModelHandle) -> CouplingHandle:
    """Declare the host coupling of an updated general model."""

    handle = error_check(lib.moist_new_coupling)(model.handle)

    def release(handle, owner=model):
        # Keep the parent alive through native deletion, including cyclic GC.
        CouplingHandle._delete(handle)

    return CouplingHandle(ffi.gc(handle, release))


def new_response() -> ResponseHandle:
    """Create an empty response handle, reusable across phases and iterations."""

    return ResponseHandle.with_gc(error_check(lib.moist_new_response)())


def prepare_model_energy(
    model: ModelHandle, coupling: CouplingHandle
) -> None:
    error_check(lib.moist_prepare_model_energy)(
        model.handle, coupling.handle
    )


def prepare_model_response(
    model: ModelHandle, coupling: CouplingHandle
) -> None:
    error_check(lib.moist_prepare_model_response)(
        model.handle, coupling.handle
    )


def prepare_model_gradient(
    model: ModelHandle, coupling: CouplingHandle
) -> None:
    error_check(lib.moist_prepare_model_gradient)(
        model.handle, coupling.handle
    )


def _accumulator(value, shape, name):
    """Require writable native buffers; never silently copy an accumulator."""
    if not isinstance(value, np.ndarray) or value.dtype != np.dtype(np.float64):
        raise TypeError(f"{name} must be a float64 NumPy array")
    if value.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not value.flags.c_contiguous or not value.flags.writeable or not value.flags.aligned:
        raise ValueError(f"{name} must be writable, aligned and C-contiguous")
    return value


def _native_array(value, name, writable=False):
    """Require a float64, aligned, C-contiguous array; its shape is trusted."""
    if not isinstance(value, np.ndarray) or value.dtype != np.dtype(np.float64):
        raise TypeError(f"{name} must be a float64 NumPy array")
    if not value.flags.c_contiguous or not value.flags.aligned:
        raise ValueError(f"{name} must be aligned and C-contiguous")
    if writable and not value.flags.writeable:
        raise ValueError(f"{name} must be writable")
    return value


def get_model_energy(model, coupling, energy) -> None:
    """Accumulate staged energy into a writable float64 scalar array."""
    _accumulator(energy, (), "energy")
    error_check(lib.moist_get_model_energy)(
        model.handle, coupling.handle, _cast("double*", energy)
    )


def get_model_response(
    model: ModelHandle, coupling: CouplingHandle, response: ResponseHandle
) -> None:
    """Fill ``response`` with the host part of a staged response phase."""

    error_check(lib.moist_get_model_response)(
        model.handle, coupling.handle, response.handle
    )


def get_model_gradient(model, coupling, response, natoms, gradient) -> None:
    """Accumulate the native contribution into ``gradient[natoms,3]``.

    ``response`` is cleared and refilled with the host part of that phase.
    """
    _accumulator(gradient, (int(natoms), 3), "gradient")
    error_check(lib.moist_get_model_gradient)(
        model.handle, coupling.handle, response.handle, int(natoms),
        _cast("double*", gradient),
    )


def next_coupling_request(coupling: CouplingHandle) -> bool:
    """Advance the cursor to the next request with a missing output of the staged phase.

    Each request is visited at most once per pass.  False ends the pass and
    rewinds the cursor, so the next call starts a new pass; an unstaged
    coupling has nothing to visit.  Otherwise only staging and a model update
    move the cursor back.  The native entry also returns false on a failure,
    which is raised here instead.
    """

    return bool(error_check(lib.moist_next_coupling_request)(coupling.handle))


def get_coupling_request_name(coupling: CouplingHandle) -> str:
    """Scientific name of the current request, e.g. ``"gaussian_potential"``."""

    buffer = ffi.new(f"char[{lib.MOIST_NAME_MAX + 1}]")
    error_check(lib.moist_get_coupling_request_name)(coupling.handle, buffer)
    return ffi.string(buffer).decode()


def get_coupling_request_missing(coupling: CouplingHandle, name: str) -> bool:
    """Whether an output of the current request is required by the staged phase and unanswered.

    A name the request does not declare is not missing; ``answer_coupling_request``
    and the ``get_*`` accessors report a misspelt output by name.
    """

    value = ffi.new("bool *")
    error_check(lib.moist_get_coupling_request_missing)(
        coupling.handle, name.encode(), value)
    return bool(value[0])


def get_coupling_request_width(coupling: CouplingHandle, width: np.ndarray) -> None:
    """Copy the Gaussian moment exponents of the current request into ``width``.

    ``width`` is a float64 array of the cavity's grid size, filled in bohr**-2.
    """

    _native_array(width, "width", writable=True)
    error_check(lib.moist_get_coupling_request_width)(
        coupling.handle, _cast("double*", width))


def answer_coupling_request(coupling: CouplingHandle, name: str, values: np.ndarray) -> None:
    """Submit one named output of the current request.

    ``values`` is a C-contiguous float64 array ``(ngrid, ...)`` on the cavity
    grid with the output's trailing extents; moist reads exactly that many
    values.  A rejected output stays missing until a valid retry; the others
    are unaffected.
    """

    _native_array(values, name)
    error_check(lib.moist_answer_coupling_request)(
        coupling.handle, name.encode(), _cast("double*", values))


def next_response_item(response: ResponseHandle) -> bool:
    """Advance the cursor to the next item of the response.

    Every item the model produced is visited once per pass, in native order.
    False ends the pass and rewinds the cursor, so the next call starts a new
    pass; an empty response gives false at once.  Filling the response again
    restarts the walk.  The native entry also returns false on a failure,
    which is raised here instead.
    """

    return bool(error_check(lib.moist_next_response_item)(response.handle))


def get_response_item_name(response: ResponseHandle) -> str:
    """Scientific name of the current item, e.g. ``"potential_adjoint"``."""

    buffer = ffi.new(f"char[{lib.MOIST_NAME_MAX + 1}]")
    error_check(lib.moist_get_response_item_name)(response.handle, buffer)
    return ffi.string(buffer).decode()


def get_response_array(response: ResponseHandle, array: str, values: np.ndarray) -> None:
    """Copy one named array of the current item into ``values``.

    ``values`` is a C-contiguous float64 array ``(ngrid, ...)``: ``w_phi``,
    ``w_rho``, ``w_overlap`` and ``w_normal_deriv`` are ``(ngrid,)``,
    ``w_grad_rho`` is ``(ngrid, 3)`` and ``w_hess_rho`` is ``(ngrid, 3, 3)``
    with indices ``[point, b, a]`` for native ``(a, b, point)``; moist writes
    exactly that many values.  No current item and an array the current item
    does not have are refused by name before anything is written.
    """

    _native_array(values, array, writable=True)
    error_check(lib.moist_get_response_array)(
        response.handle, array.encode(), _cast("double*", values)
    )


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
    xyz = np.zeros((ngrid, 3), dtype=np.float64, order="C")
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


# Element type tags from moist.h; a field is read with the accessor matching its
# tag. Macro values are supplied by the compiled cffi extension.
FIELD_REAL = lib.MOIST_FIELD_REAL
FIELD_INT = lib.MOIST_FIELD_INT
FIELD_BOOL = lib.MOIST_FIELD_BOOL

_FIELD_READER = {
    FIELD_REAL: ("moist_get_cavity_field_real", np.float64, "double*"),
    FIELD_INT: ("moist_get_cavity_field_int", np.int32, "int*"),
    FIELD_BOOL: ("moist_get_cavity_field_bool", np.bool_, "bool*"),
}

_FIELD_MAX_RANK = lib.MOIST_FIELD_MAX_RANK

_FIELD_NAME_CAP = lib.MOIST_FIELD_NAME_MAX + 1


@dataclass(frozen=True)
class CavityField:
    """Shape and type of one readable cavity field.

    ``shape`` is empty for a scalar and otherwise carries the extents in
    C order, slowest-varying first. Grid vectors have shape ``(ngrid, 3)``.
    """

    name: str
    dtype: np.dtype
    shape: tuple[int, ...]
    count: int


def get_cavity_field_count(cavity: CavityHandle) -> int:
    """Return how many named result fields the cavity currently holds."""

    nfield = ffi.new("int *")
    error_check(lib.moist_get_cavity_field_count)(cavity.handle, nfield)
    return int(nfield[0])


def get_cavity_field_info(cavity: CavityHandle, index: int) -> CavityField:
    """Describe the field at ``index``, counting from zero."""

    name = ffi.new(f"char[{_FIELD_NAME_CAP}]")
    dtype = ffi.new("int *")
    rank = ffi.new("int *")
    dims = ffi.new(f"int[{_FIELD_MAX_RANK}]")
    count = ffi.new("int *")

    error_check(lib.moist_get_cavity_field_info)(
        cavity.handle,
        int(index),
        name,
        dtype,
        rank,
        dims,
        count,
    )

    tag = int(dtype[0])
    if tag not in _FIELD_READER:
        raise ValueError(f"moist reported an unknown field type tag {tag}")

    return CavityField(
        name=ffi.string(name).decode(),
        dtype=np.dtype(_FIELD_READER[tag][1]),
        shape=tuple(int(dims[i]) for i in range(int(rank[0]))),
        count=int(count[0]),
    )


def list_cavity_fields(cavity: CavityHandle) -> tuple[CavityField, ...]:
    """Describe every result the cavity currently holds.

    The list is what the cavity itself declares, so it grows with the cavity
    type and omits optional properties that were never computed. Read it again
    after an update: a rebuild changes the extents.
    """

    return tuple(
        get_cavity_field_info(cavity, index)
        for index in range(get_cavity_field_count(cavity))
    )


def get_cavity_field_about(cavity: CavityHandle, name: str) -> str:
    """Return the one-line description moist attaches to a field."""

    length = ffi.new("size_t *")
    get_about = error_check(lib.moist_get_cavity_field_about)
    get_about(cavity.handle, _char(name), ffi.NULL, 0, length)
    capacity = length[0] + 1
    about = ffi.new("char[]", capacity)
    get_about(cavity.handle, _char(name), about, capacity, length)
    return ffi.string(about).decode()


def get_cavity_field(
    cavity: CavityHandle,
    name: str,
    info: Optional[CavityField] = None,
) -> np.ndarray:
    """Return one named cavity result.

    Fields come back C-contiguous with the reported C shape. A name the
    cavity does not currently hold -- unknown, or an optional property that was
    not requested -- raises rather than returning zeros.

    ``info`` skips the shape lookup when the caller already has a descriptor
    from :func:`list_cavity_fields`; it has to describe the same cavity state.
    """

    if info is None:
        info = _find_cavity_field(cavity, name)

    reader, dtype, ctype = _FIELD_READER[_tag_of(info.dtype)]
    values = np.zeros(info.count, dtype=dtype)
    error_check(getattr(lib, reader))(
        cavity.handle,
        _char(name),
        _cast(ctype, values),
    )

    if not info.shape:
        return values[0]
    return values.reshape(info.shape, order="C")


def get_cavity_fields(
    cavity: CavityHandle,
    names: Optional[Iterable[str]] = None,
) -> dict:
    """Return the named cavity results as a ``{name: value}`` mapping.

    With ``names`` omitted this reads everything the cavity holds; the fields
    are enumerated once and the descriptors reused, so the mapping is
    consistent with a single cavity state.
    """

    fields = list_cavity_fields(cavity)
    if names is not None:
        wanted = list(names)
        known = {field.name: field for field in fields}
        missing = [name for name in wanted if name not in known]
        if missing:
            raise KeyError(
                "cavity does not hold the field(s) "
                + ", ".join(repr(name) for name in missing)
                + "; available: "
                + ", ".join(field.name for field in fields)
            )
        fields = [known[name] for name in wanted]

    return {
        field.name: get_cavity_field(cavity, field.name, info=field)
        for field in fields
    }


def _find_cavity_field(cavity: CavityHandle, name: str) -> CavityField:
    for field in list_cavity_fields(cavity):
        if field.name == name:
            return field
    raise KeyError(f"cavity does not hold a field named {name!r}")


def _tag_of(dtype: np.dtype) -> int:
    for tag, (_, candidate, _) in _FIELD_READER.items():
        if dtype == np.dtype(candidate):
            return tag
    raise ValueError(f"no moist field accessor for dtype {dtype}")


def assemble_drop_amat(cavity: CavityHandle) -> tuple[np.ndarray, np.ndarray]:
    """Assemble the Gaussian CPCM A-matrix and return it with xi values.

    No longer DROP-specific: the underlying C entry point now works for every
    Gaussian-discretized cavity.
    """

    ngrid, _ = get_cavity_sizes(cavity)
    amat = np.zeros((ngrid, ngrid), dtype=np.float64, order="C")
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


def compute_cavity_gradient(cavity):
    """Build diagnostic forward derivatives; normal model gradients avoid this."""
    error_check(lib.moist_compute_cavity_gradient)(cavity.handle)


def contract_amat_nuclear_gradient(cavity, q1, q2):
    """Diagnostic forward contraction, after compute_cavity_gradient."""
    ngrid, nsph = get_cavity_sizes(cavity)
    arrays = [np.asarray(q, dtype=np.float64, order="C") for q in (q1, q2)]
    if any(a.shape != (ngrid,) for a in arrays):
        raise ValueError("q1 and q2 must have shape (ngrid,)")
    gradient = np.zeros((nsph, 3))
    error_check(lib.moist_contract_amat1_q1q2_rA)(cavity.handle,
        *[_cast("double*", a) for a in arrays], _cast("double*", gradient))
    return gradient


def contract_pcm_nuclear_gradient(cavity, w_phi, w_xyz, charges):
    """Diagnostic electrostatic contraction, after compute_cavity_gradient.

    ``w_phi`` is the potential adjoint ``dE/dphi`` (the surface charge of a
    stationary PCM), ``charges`` the nuclear charges.
    """
    ngrid, nsph = get_cavity_sizes(cavity)
    arrays = [np.asarray(a, dtype=np.float64, order="C")
              for a in (w_phi, w_xyz, charges)]
    for a, shape in zip(arrays, ((ngrid,), (ngrid, 3), (nsph,))):
        if a.shape != shape:
            raise ValueError(f"contraction input must have shape {shape}")
    gradient = np.zeros((nsph, 3))
    error_check(lib.moist_contract_pcm_nuclear_gradient)(cavity.handle,
        *[_cast("double*", a) for a in arrays], _cast("double*", gradient))
    return gradient


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

    xyz1_rA = np.zeros((ngrid, nsph, 3, 3), dtype=np.float64, order="C")
    xi1_rA = np.zeros((ngrid, nsph, 3), dtype=np.float64, order="C")
    a_i1_rA = np.zeros((ngrid, nsph, 3), dtype=np.float64, order="C")
    v_i1_rA = np.zeros((ngrid, nsph, 3), dtype=np.float64, order="C")
    A_tot1_rA = np.zeros((nsph, 3), dtype=np.float64, order="C")
    V_tot1_rA = np.zeros((nsph, 3), dtype=np.float64, order="C")

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
    w_xyz = np.zeros((ngrid, 3), dtype=np.float64, order="C")

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
    w_normal: Optional[np.ndarray] = None,
    w_k1: Optional[np.ndarray] = None,
    w_k2: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contract DROP surface weights to LSF adjoint weights.

    The outward-normal (``w_normal``, shape ``(ngrid, 3)``) and principal-curvature
    (``w_k1``/``w_k2``, shape ``(ngrid,)``) channels are optional; ``None``
    skips the channel entirely.
    """

    ngrid, _ = get_cavity_sizes(cavity)
    _w_xi = np.ascontiguousarray(w_xi, dtype=np.float64)
    _w_f = np.ascontiguousarray(w_f, dtype=np.float64)
    _w_xyz = np.asarray(w_xyz, dtype=np.float64, order="C")
    if _w_xi.shape != (ngrid,) or _w_f.shape != (ngrid,) or _w_xyz.shape != (ngrid, 3):
        raise ValueError("w_xi/w_f must have shape (ngrid,), w_xyz must have shape (ngrid, 3)")

    _w_normal = (
        None if w_normal is None else np.asarray(w_normal, dtype=np.float64, order="C")
    )
    _w_k1 = None if w_k1 is None else np.ascontiguousarray(w_k1, dtype=np.float64)
    _w_k2 = None if w_k2 is None else np.ascontiguousarray(w_k2, dtype=np.float64)
    if _w_normal is not None and _w_normal.shape != (ngrid, 3):
        raise ValueError("w_normal must have shape (ngrid, 3)")
    if _w_k1 is not None and _w_k1.shape != (ngrid,):
        raise ValueError("w_k1 must have shape (ngrid,)")
    if _w_k2 is not None and _w_k2.shape != (ngrid,):
        raise ValueError("w_k2 must have shape (ngrid,)")

    w_lsf0 = np.zeros(ngrid, dtype=np.float64)
    w_lsf1 = np.zeros((ngrid, 3), dtype=np.float64, order="C")
    w_lsf2 = np.zeros((ngrid, 3, 3), dtype=np.float64, order="C")

    error_check(lib.moist_contract_surface_lsf_weights_extended)(
        cavity.handle,
        _cast("double*", _w_xi),
        _cast("double*", _w_f),
        _cast("double*", _w_xyz),
        _cast("double*", w_lsf0),
        _cast("double*", w_lsf1),
        _cast("double*", w_lsf2),
        _cast("double*", _w_normal),
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
