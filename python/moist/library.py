"""Thin wrapper around the moist CFFI extension."""

import functools
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


# def new_gems_model(
#     solvent: str,
#     debug: bool = False,
#     verbosity: int = 0,
#     parameter_file: Optional[str] = None,
# ) -> ModelHandle:
#     read_parameters = parameter_file is not None
#     return ModelHandle.with_gc(
#         error_check(lib.moist_new_gems_solvation_model)(
#             _char(solvent),
#             _ref("bool", debug),
#             _ref("int", verbosity),
#             _ref("bool", read_parameters),
#             _char(parameter_file),
#         )
#     )


def new_drop_cavity(
    nleb: Optional[int] = None,
    debug: bool = False,
    verbosity: int = 0,
    do_fine: bool = False,
) -> CavityHandle:
    """
    Create a standard solute-vdW (SvdW) DROP cavity with default CPCM radii.
    """

    return CavityHandle.with_gc(
        error_check(lib.moist_new_drop_cavity)(
            _ref("int", nleb),
            _ref("bool", debug),
            _ref("int", verbosity),
            _ref("double", None),
            _ref("double", None),
            _ref("double", None),
            _ref("double", None),
            _ref("bool", do_fine),
        )
    )


def new_drop_cavity_isodensity_callback(
    callback,
    nleb: Optional[int] = None,
    scale: float = 1000.0,
    debug: bool = False,
    verbosity: int = 0,
    do_fine: bool = False,
    wleb_prune_level: Optional[int] = None,
) -> tuple[CavityHandle, object]:
    """Create a DROP cavity backed by a Python isodensity LSF callback.

    The Python callback receives a single point in Bohr and must return either
    ``(value, grad, hess, third)`` or an object with ``value``, ``grad``,
    ``hess``, and ``third`` attributes. The returned CFFI callback must be kept alive by the caller for
    at least as long as the cavity handle.
    """

    @ffi.callback("moist_isodensity_lsf_callback")
    def c_callback(context, point_ptr, value_ptr, grad_ptr, hess_ptr, third_ptr):
        point = np.array([point_ptr[i] for i in range(3)], dtype=np.float64)
        result = callback(point)

        if hasattr(result, "value"):
            value = result.value
            grad = result.grad
            hess = result.hess
            third = result.third
        else:
            value, grad, hess, third = result

        grad = np.asarray(grad, dtype=np.float64)
        hess = np.asarray(hess, dtype=np.float64)
        third = np.asarray(third, dtype=np.float64)
        if grad.shape != (3,):
            raise ValueError("Isodensity callback gradient must have shape (3,)")
        if hess.shape != (3, 3):
            raise ValueError("Isodensity callback Hessian must have shape (3, 3)")
        if third.shape != (3, 3, 3):
            raise ValueError("Isodensity callback third derivative must have shape (3, 3, 3)")

        value_ptr[0] = float(value)
        for i in range(3):
            grad_ptr[i] = float(grad[i])
        for j in range(3):
            for i in range(3):
                hess_ptr[i + 3 * j] = float(hess[i, j])
        for k in range(3):
            for j in range(3):
                for i in range(3):
                    third_ptr[i + 3 * j + 9 * k] = float(third[i, j, k])

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
        )
    )
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


def assemble_drop_amat(cavity: CavityHandle) -> tuple[np.ndarray, np.ndarray]:
    """Assemble the DROP CPCM A-matrix and return it with xi values."""

    ngrid, _ = get_cavity_sizes(cavity)
    amat = np.zeros((ngrid, ngrid), dtype=np.float64, order="F")
    xi = np.zeros(ngrid, dtype=np.float64)

    error_check(lib.moist_assemble_amat)(
        cavity.handle,
        _ref("int", ngrid),
        _cast("double*", amat),
        _cast("double*", xi),
    )
    return amat, xi


def contract_amat1_q1q2_surface_weights(
    cavity: CavityHandle,
    q1: np.ndarray,
    q2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contract DROP CPCM A-matrix derivatives to per-grid surface weights."""

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
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contract DROP surface weights to LSF adjoint weights."""

    ngrid, _ = get_cavity_sizes(cavity)
    _w_xi = np.ascontiguousarray(w_xi, dtype=np.float64)
    _w_f = np.ascontiguousarray(w_f, dtype=np.float64)
    _w_xyz = np.asarray(w_xyz, dtype=np.float64, order="F")
    if _w_xi.shape != (ngrid,) or _w_f.shape != (ngrid,) or _w_xyz.shape != (3, ngrid):
        raise ValueError("w_xi/w_f must have shape (ngrid,), w_xyz must have shape (3, ngrid)")

    w_lsf0 = np.zeros(ngrid, dtype=np.float64)
    w_lsf1 = np.zeros((3, ngrid), dtype=np.float64, order="F")
    w_lsf2 = np.zeros((3, 3, ngrid), dtype=np.float64, order="F")

    error_check(lib.moist_contract_surface_lsf_weights)(
        cavity.handle,
        _cast("double*", _w_xi),
        _cast("double*", _w_f),
        _cast("double*", _w_xyz),
        _cast("double*", w_lsf0),
        _cast("double*", w_lsf1),
        _cast("double*", w_lsf2),
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
