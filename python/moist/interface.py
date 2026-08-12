"""Object-oriented Python interface for moist solvation models.

The native C interface is intentionally procedural.  This module puts the
Python seam around a complete model evaluation instead: live cavity objects own
their behaviour, snapshots are explicit values, coupling adapters hide host
exchange ordering, and :class:`Evaluation` represents one coherent model state.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from enum import Enum, IntEnum
from typing import Callable, Iterator, Optional, Protocol, Union
import warnings

import numpy as np

from . import library


# -----------------------------------------------------------------------------
# Result values
# -----------------------------------------------------------------------------


def _immutable_array(array: np.ndarray) -> np.ndarray:
    """Copy an array onto a buffer whose write protection cannot be reversed."""
    order = "F" if array.flags.f_contiguous and not array.flags.c_contiguous else "C"
    buffer = array.tobytes(order=order)
    return np.frombuffer(buffer, dtype=array.dtype).reshape(array.shape, order=order)


def _freeze_result_arrays(value: _ImmutableArrayValue) -> None:
    """Replace ndarray fields with immutable-buffer copies."""
    for field in fields(value):
        array = getattr(value, field.name)
        if isinstance(array, np.ndarray):
            object.__setattr__(value, field.name, _immutable_array(array))


class _ImmutableArrayValue:
    """Dataclass mixin that gives ndarray fields immutable backing buffers."""

    def __post_init__(self) -> None:
        _freeze_result_arrays(self)


@dataclass(frozen=True)
class CavitySnapshot(_ImmutableArrayValue):
    """Generic cavity data copied from one successful update.

    Field names are moist's own, so a value read here can be matched against the
    native arrays without a translation step:

    ``xyz``
        ``(3, ngrid)`` grid-point coordinates in bohr.
    ``a``
        ``(ngrid,)`` grid-point areas.
    ``owner``
        ``(ngrid,)`` index of the sphere each grid point belongs to.
    ``converged``
        ``(ngrid,)`` per-point projection success flags.
    ``radii``, ``asph``
        ``(nsph,)`` sphere radii and per-sphere surface areas.
    """

    area: float
    volume: float
    ngrid: int
    nsph: int
    xyz: np.ndarray
    a: np.ndarray
    owner: np.ndarray
    converged: np.ndarray
    radii: np.ndarray
    asph: np.ndarray


@dataclass(frozen=True)
class CavitySnapshotDROP(CavitySnapshot):
    """DROP-specific cavity data copied from one successful update.

    ``nmax``
        Grid points allocated per sphere before pruning.
    ``normal0``
        ``(3, ngrid)`` initial (pre-projection) surface normals.
    ``wleb``
        ``(ngrid,)`` Lebedev quadrature weights.
    ``r_iI0``
        ``(3, ngrid)`` displacement of each grid point from its anchor atom.
    ``f``
        ``(ngrid,)`` switching-function values.
    ``rho``
        ``(ngrid,)`` level-set (density) values at the projected points.
    """

    nmax: int
    normal0: np.ndarray
    wleb: np.ndarray
    r_iI0: np.ndarray
    f: np.ndarray
    rho: np.ndarray


@dataclass(frozen=True)
class AnchorGradient(_ImmutableArrayValue):
    """Anchor-channel nuclear derivatives of a DROP cavity, native grid order.

    Every array is Fortran-ordered, matching moist's own layout:

    ``xyz1_rA``
        ``(3, 3, nsph, ngrid)`` indexed ``(j, alpha, A, i)`` --
        ``d(r_i)_j / d(R_A)_alpha``.
    ``xi1_rA``, ``a_i1_rA``, ``v_i1_rA``
        ``(3, nsph, ngrid)`` indexed ``(alpha, A, i)``.
    ``A_tot1_rA``, ``V_tot1_rA``
        ``(3, nsph)`` -- the grid sums of ``a_i1_rA``/``v_i1_rA``.
    """

    xyz1_rA: np.ndarray
    xi1_rA: np.ndarray
    a_i1_rA: np.ndarray
    v_i1_rA: np.ndarray
    A_tot1_rA: np.ndarray
    V_tot1_rA: np.ndarray


@dataclass(frozen=True)
class TracePotential(_ImmutableArrayValue):
    """Direct host-trace adjoints used while preparing a coupling."""

    molecular: np.ndarray
    normal: np.ndarray

    def __iter__(self) -> Iterator[np.ndarray]:
        """Preserve tuple unpacking used by the pre-evaluation interface."""
        yield self.molecular
        yield self.normal


@dataclass(frozen=True)
class GeneralPotential(_ImmutableArrayValue):
    """Every adjoint channel returned by a general solvation model.

    Each ``w_*`` is the derivative of the model energy with respect to the host
    quantity it names, so the host completes the chain rule by contracting it
    with its own derivative of that quantity:

    ``w_umol``, ``w_qmol``
        ``(ngrid,)`` adjoints of the molecular potential and of the surface
        charges.
    ``w_lsf0``, ``w_lsf1``, ``w_lsf2``
        ``(ngrid,)``, ``(3, ngrid)`` and ``(3, 3, ngrid)`` adjoints of the
        **scaled** level set and of its gradient and Hessian.
    ``w_gauss_g``, ``w_gauss_f``
        ``(ngrid,)`` GOSTSHYP amplitudes for the Gaussian value and normal
        derivative; ``None`` unless a GOSTSHYP component is present.
    """

    w_umol: np.ndarray
    w_qmol: np.ndarray
    w_lsf0: np.ndarray
    w_lsf1: np.ndarray
    w_lsf2: np.ndarray
    w_gauss_g: Optional[np.ndarray] = None
    w_gauss_f: Optional[np.ndarray] = None


# -----------------------------------------------------------------------------
# Host input values
# -----------------------------------------------------------------------------
#
# These travel the other way: a coupling adapter fills them from host data and
# hands them to the model.  They deliberately do not use _ImmutableArrayValue,
# because the buffers belong to the host and copying them on every exchange
# would cost more than the protection is worth here.


@dataclass(frozen=True)
class Electrostatics:
    """Host electrostatic traces supplied for one cavity surface.

    ``phi`` is the only universally required field.  The response arrays are
    optional at this low-level value-object seam; coupling adapters decide which
    ones are required for a complete evaluation and supply them atomically.
    """

    phi: np.ndarray
    w_xi: Optional[np.ndarray] = None
    w_f: Optional[np.ndarray] = None
    w_xyz: Optional[np.ndarray] = None
    w_n: Optional[np.ndarray] = None
    qefield: Optional[np.ndarray] = None


@dataclass(frozen=True)
class GostshypMoments:
    """Gaussian density moments consumed by a GOSTSHYP model component."""

    gt: np.ndarray
    pt: np.ndarray
    mt: np.ndarray
    rt: np.ndarray


# -----------------------------------------------------------------------------
# Molecular structure
# -----------------------------------------------------------------------------


class Structure:
    """Validated molecular structure owning a native moist handle.

    Coordinates and lattice vectors are in Bohr.  Positions use the natural
    NumPy shape ``(natoms, 3)``; the wrapper converts that row-major memory to
    the native Fortran ``(3, natoms)`` view without an intermediate transpose.
    """

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ) -> None:
        _numbers = np.asarray(numbers)
        if _numbers.ndim != 1:
            raise ValueError("numbers must have shape (natoms,)")
        natoms = int(_numbers.size)

        _positions = self._positions_array(positions, natoms)
        _lattice = self._lattice_array(lattice)

        if periodic is None:
            _periodic = None
        else:
            _periodic = np.asarray(periodic)
            if _periodic.shape != (3,):
                raise ValueError("periodic must have shape (3,)")

        self._natoms = natoms
        self._numbers = np.ascontiguousarray(_numbers, dtype=np.int32)
        self._positions = _positions
        self._lattice = _lattice
        self._periodic = (
            None if _periodic is None else np.ascontiguousarray(_periodic, dtype=np.bool_)
        )
        self._handle = library.new_structure(
            self._natoms,
            self._numbers,
            self._positions,
            self._lattice,
            self._periodic,
        )

    @staticmethod
    def _positions_array(positions: np.ndarray, natoms: int) -> np.ndarray:
        array = np.asarray(positions)
        if array.shape != (natoms, 3):
            raise ValueError(f"positions must have shape ({natoms}, 3)")
        return np.ascontiguousarray(array, dtype=np.float64)

    @staticmethod
    def _lattice_array(lattice: Optional[np.ndarray]) -> Optional[np.ndarray]:
        if lattice is None:
            return None
        array = np.asarray(lattice)
        if array.shape != (3, 3):
            raise ValueError("lattice must have shape (3, 3)")
        return np.ascontiguousarray(array, dtype=np.float64)

    def _same_geometry(self, other: Structure) -> bool:
        """Whether two structures describe exactly the same native update."""
        return (
            np.array_equal(self._numbers, other._numbers)
            and np.array_equal(self._positions, other._positions)
            and self._optional_array_equal(self._lattice, other._lattice)
            and self._optional_array_equal(self._periodic, other._periodic)
        )

    @staticmethod
    def _optional_array_equal(
        left: Optional[np.ndarray], right: Optional[np.ndarray]
    ) -> bool:
        return (left is None and right is None) or (
            left is not None and right is not None and np.array_equal(left, right)
        )

    def __len__(self) -> int:
        return self._natoms

    @property
    def natoms(self) -> int:
        return self._natoms

    @property
    def numbers(self) -> np.ndarray:
        return self._numbers.copy()

    @property
    def positions(self) -> np.ndarray:
        return self._positions.copy()

    @property
    def lattice(self) -> Optional[np.ndarray]:
        return None if self._lattice is None else self._lattice.copy()

    @property
    def periodic(self) -> Optional[np.ndarray]:
        return None if self._periodic is None else self._periodic.copy()

    @property
    def _mol(self) -> library.StructureHandle:
        """Compatibility name for package-internal pre-refactor callers."""
        return self._handle

    def _as_handle(self) -> library.StructureHandle:
        return self._handle

    def update(self, positions: np.ndarray, lattice: Optional[np.ndarray] = None) -> None:
        _positions = self._positions_array(positions, self._natoms)
        _lattice = self._lattice_array(lattice)

        library.update_structure(self._handle, _positions, _lattice)
        self._positions = _positions
        if _lattice is not None:
            self._lattice = _lattice


# -----------------------------------------------------------------------------
# Live cavities
# -----------------------------------------------------------------------------


def _guarded_native_update(cavity: Cavity, native: Callable[[], None]) -> None:
    """Run one native update, surfacing a Python callback failure either way.

    A callback-backed cavity records an exception raised inside its level-set
    callback rather than letting it cross the native frames, so the recorded
    failure has to be re-raised on both paths: after a native error, where it is
    the more specific cause, and after an apparently successful call, which the
    native side can report when the callback was the one that failed.  Both
    :meth:`Cavity.update` and :meth:`SolvationModel.update` drive the same
    protocol; only the surrounding bookkeeping differs.
    """
    cavity._before_native_update()
    try:
        native()
    except Exception:
        cavity._raise_callback_failure()
        raise
    cavity._raise_callback_failure()


class Cavity(ABC):
    """Live cavity object.

    A cavity owns behaviour and native state.  :meth:`snapshot` returns an
    explicit copied value for a particular successful update.  Model-owned
    cavity views expose the same read/derivative behaviour but cannot be rebuilt
    independently of their model.
    """

    density_dependent = False

    def __init__(self, handle: library.CavityHandle, *, owned: bool = True) -> None:
        self._handle = handle
        self._owned = owned
        self._updated = False
        self._snapshot_cache: Optional[CavitySnapshot] = None

    def _as_handle(self) -> library.CavityHandle:
        return self._handle

    def _invalidate(self) -> None:
        self._updated = False
        self._snapshot_cache = None

    def _mark_updated(self) -> None:
        self._updated = True
        self._snapshot_cache = None

    def _before_native_update(self) -> None:
        """Hook for callback-backed cavities."""

    def _raise_callback_failure(self) -> None:
        """Hook for callback-backed cavities."""

    def update(self, structure: Structure) -> None:
        if not self._owned:
            raise RuntimeError("A model-owned cavity must be updated through its model")
        self._invalidate()
        _guarded_native_update(
            self, lambda: library.update_cavity(self._handle, structure._as_handle())
        )
        self._mark_updated()

    def _require_updated(self) -> None:
        if not self._updated:
            raise RuntimeError("Cavity has not been successfully updated")

    def snapshot(self) -> CavitySnapshot:
        """Copy the current native cavity results into an immutable value object."""
        self._require_updated()
        if self._snapshot_cache is None:
            self._snapshot_cache = self._read_snapshot()
        return self._snapshot_cache

    @abstractmethod
    def _read_snapshot(self) -> CavitySnapshot:
        """Read a concrete snapshot from the native handle."""

    @abstractmethod
    def _model_view(self, handle: library.CavityHandle) -> Cavity:
        """Wrap the model-owned native copy without taking update ownership."""

    @property
    def cavity(self) -> CavitySnapshot:
        """Compatibility alias for :meth:`snapshot`."""
        return self.snapshot()

    @property
    def area(self) -> float:
        return self.snapshot().area

    @property
    def volume(self) -> float:
        return self.snapshot().volume

    @property
    def ngrid(self) -> int:
        return self.snapshot().ngrid

    @property
    def nsph(self) -> int:
        return self.snapshot().nsph

    @property
    def xyz(self) -> np.ndarray:
        return self.snapshot().xyz

    @property
    def a(self) -> np.ndarray:
        return self.snapshot().a

    @property
    def owner(self) -> np.ndarray:
        return self.snapshot().owner

    @property
    def converged(self) -> np.ndarray:
        return self.snapshot().converged

    @property
    def radii(self) -> np.ndarray:
        return self.snapshot().radii

    @property
    def asph(self) -> np.ndarray:
        return self.snapshot().asph


class _CavityGenericBase(Cavity):
    """Shared implementation for non-DROP native cavities."""

    def _read_snapshot(self) -> CavitySnapshot:
        return CavitySnapshot(**library.get_cavity_results(self._handle))

    def _model_view(self, handle: library.CavityHandle) -> Cavity:
        return _CavityGenericBorrowed(handle, self)


class _CavityGenericBorrowed(_CavityGenericBase):
    """High-level view of a non-DROP cavity owned by a model."""

    def __init__(self, handle: library.CavityHandle, source: Cavity) -> None:
        super().__init__(handle, owned=False)
        self._source = source

    @property
    def density_dependent(self) -> bool:
        return self._source.density_dependent


class CavityISwiG(_CavityGenericBase):
    """iSwiG switching-Gaussian cavity with default CPCM radii.

    ``nleb`` controls the Lebedev grid.  ``cut_a`` selects an area cutoff when
    positive; otherwise ``cut_f`` is the switching-function cutoff.
    """

    def __init__(
        self,
        nleb: Optional[int] = None,
        cut_a: Optional[float] = None,
        cut_f: Optional[float] = None,
        debug: bool = False,
        verbosity: int = 0,
    ) -> None:
        super().__init__(
            library.new_iswig_cavity(
                nleb=nleb,
                cut_a=cut_a,
                cut_f=cut_f,
                debug=debug,
                verbosity=verbosity,
            )
        )


class _CavityDROPBase(Cavity):
    """Shared behaviour for standalone and model-owned DROP cavities."""

    def _read_snapshot(self) -> CavitySnapshotDROP:
        generic = library.get_cavity_results(self._handle)
        drop = library.get_drop_specific(self._handle, ngrid=generic["ngrid"])
        return CavitySnapshotDROP(**generic, **drop)

    def _model_view(self, handle: library.CavityHandle) -> Cavity:
        return _CavityDROPBorrowed(handle, self)

    def assemble_amat(self) -> tuple[np.ndarray, np.ndarray]:
        self._require_updated()
        return library.assemble_drop_amat(self._handle)

    def get_gaussian(self) -> tuple[np.ndarray, np.ndarray]:
        """Return Gaussian widths and switching factors without assembling A."""
        self._require_updated()
        return library.get_cavity_gaussian(self._handle)

    def compute_anchor_gradient(self) -> None:
        """Compute the anchor-only nuclear derivatives."""
        self._require_updated()
        library.compute_anchor_gradient(self._handle)

    def get_anchor_gradient(self) -> AnchorGradient:
        """Return the anchor-channel nuclear derivatives in native grid order."""
        self._require_updated()
        return AnchorGradient(**library.get_anchor_gradient(self._handle))

    def contract_amat_surface_weights(
        self,
        q1: np.ndarray,
        q2: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._require_updated()
        return library.contract_amat1_q1q2_surface_weights(self._handle, q1, q2)

    def contract_surface_lsf_weights(
        self,
        w_xi: np.ndarray,
        w_f: np.ndarray,
        w_xyz: np.ndarray,
        w_n: Optional[np.ndarray] = None,
        w_k1: Optional[np.ndarray] = None,
        w_k2: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._require_updated()
        return library.contract_surface_lsf_weights(
            self._handle, w_xi, w_f, w_xyz, w_n, w_k1, w_k2
        )

    @property
    def nmax(self) -> int:
        return self.snapshot().nmax

    @property
    def normal0(self) -> np.ndarray:
        return self.snapshot().normal0

    @property
    def wleb(self) -> np.ndarray:
        return self.snapshot().wleb

    @property
    def r_iI0(self) -> np.ndarray:
        return self.snapshot().r_iI0

    @property
    def f(self) -> np.ndarray:
        return self.snapshot().f

    @property
    def rho(self) -> np.ndarray:
        return self.snapshot().rho


class _CavityDROPBorrowed(_CavityDROPBase):
    """High-level view of the authoritative cavity owned by a model."""

    def __init__(self, handle: library.CavityHandle, source: Cavity) -> None:
        super().__init__(handle, owned=False)
        self._source = source  # also keeps callback-backed source objects alive

    @property
    def density_dependent(self) -> bool:
        return self._source.density_dependent


class CavityDROPSvdW(_CavityDROPBase):
    """Smooth-van-der-Waals DROP cavity with default CPCM radii."""

    def __init__(
        self,
        nleb: Optional[int] = None,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
        tolerance: Optional[float] = None,
        blend_k: Optional[float] = None,
        blend_1b: Optional[float] = None,
        blend_2b: Optional[float] = None,
        blend_3b: Optional[float] = None,
        proj_maxiter: Optional[int] = None,
        proj_level: Optional[int] = None,
        branch_weight_s: Optional[float] = None,
        rho_grid_h: Optional[float] = None,
        wleb_prune_level: Optional[int] = None,
    ) -> None:
        super().__init__(
            library.new_drop_cavity(
                nleb=nleb,
                debug=debug,
                verbosity=verbosity,
                blend_k=blend_k,
                blend_1b=blend_1b,
                blend_2b=blend_2b,
                blend_3b=blend_3b,
                do_fine=do_fine,
                tolerance=tolerance,
                proj_maxiter=proj_maxiter,
                proj_level=proj_level,
                branch_weight_s=branch_weight_s,
                rho_grid_h=rho_grid_h,
                wleb_prune_level=wleb_prune_level,
            )
        )


class CavityDROPCFC(_CavityDROPBase):
    """COSMO Fine Cavity discretized with DROP."""

    def __init__(
        self,
        nleb: Optional[int] = None,
        a1: Optional[float] = None,
        a2: Optional[float] = None,
        c: Optional[float] = None,
        m: Optional[int] = None,
        screen_k: Optional[float] = None,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
        tolerance: Optional[float] = None,
        proj_maxiter: Optional[int] = None,
        proj_level: Optional[int] = None,
        branch_weight_s: Optional[float] = None,
        rho_grid_h: Optional[float] = None,
        wleb_prune_level: Optional[int] = None,
    ) -> None:
        super().__init__(
            library.new_cfc_drop_cavity(
                nleb=nleb,
                a1=a1,
                a2=a2,
                c=c,
                m=m,
                screen_k=screen_k,
                debug=debug,
                verbosity=verbosity,
                do_fine=do_fine,
                tolerance=tolerance,
                proj_maxiter=proj_maxiter,
                proj_level=proj_level,
                branch_weight_s=branch_weight_s,
                rho_grid_h=rho_grid_h,
                wleb_prune_level=wleb_prune_level,
            )
        )


# Short name for the default DROP cavity: SvdW is the surface a caller who does
# not name a level set means.
CavityDROP = CavityDROPSvdW


#: ``(value, grad)``, ``(value, grad, hess)`` or ``(value, grad, hess, third)``:
#: the tuple grows with the requested order so a caller never pays for a
#: derivative it did not ask for.
LevelSetDerivatives = tuple[Union[float, np.ndarray], ...]

#: A raw level-set callback.  It takes ``(point, order)``, or just ``(point)``
#: when the cavity was built with ``pass_order=False``.
IsodensityCallback = Callable[..., LevelSetDerivatives]


class IsodensitySource(Protocol):
    """Provider of an unscaled isodensity level set and its MOIST scale."""

    scale: float

    def lsf(self, point: np.ndarray, order: int) -> LevelSetDerivatives:
        """Return the level-set value and spatial derivatives through ``order``."""
        ...


class CavityDROPIsodensity(_CavityDROPBase):
    """DROP cavity driven by an isodensity source or a raw Python callback.

    A source exposes ``lsf(point, order)`` and ``scale``.  Passing a source is
    the canonical interface because it keeps the callback and its scaling
    invariant together.  Raw callbacks remain supported for non-host callers.
    """

    density_dependent = True

    def __init__(
        self,
        source: Optional[Union[IsodensitySource, IsodensityCallback]] = None,
        nleb: Optional[int] = None,
        scale: Optional[float] = None,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
        wleb_prune_level: Optional[int] = None,
        tolerance: Optional[float] = None,
        pass_order: Optional[bool] = None,
        callback: Optional[IsodensityCallback] = None,
    ) -> None:
        if callback is not None:
            if source is not None:
                raise TypeError("Pass source or callback, not both")
            warnings.warn(
                "callback= is deprecated; pass the callback as source instead",
                DeprecationWarning,
                stacklevel=2,
            )
            source = callback
        if source is None:
            raise TypeError("CavityDROPIsodensity requires a level-set source")

        provider_callback = getattr(source, "lsf", None)
        if callable(provider_callback):
            if not hasattr(source, "scale"):
                raise TypeError("An isodensity source must expose a scale")
            source_scale = float(source.scale)
            if scale is not None and float(scale) != source_scale:
                raise ValueError("cavity scale must match the isodensity source scale")
            resolved_callback = provider_callback
            resolved_scale = source_scale
        elif callable(source):
            resolved_callback = source
            resolved_scale = 1000.0 if scale is None else float(scale)
        else:
            raise TypeError(
                "source must be callable or expose callable lsf(point, order)"
            )

        handle, callback_ref = library.new_drop_cavity_isodensity_callback(
            callback=resolved_callback,
            nleb=nleb,
            scale=resolved_scale,
            debug=debug,
            verbosity=verbosity,
            do_fine=do_fine,
            wleb_prune_level=wleb_prune_level,
            tolerance=tolerance,
            pass_order=pass_order,
        )
        super().__init__(handle)
        self._source = source
        self._callback_ref = callback_ref

    def _before_native_update(self) -> None:
        self._handle.callback_state.reset()

    def _raise_callback_failure(self) -> None:
        self._handle.callback_state.raise_if_failed()


# -----------------------------------------------------------------------------
# Solvation component configurations
# -----------------------------------------------------------------------------


class CouplingChannel(str, Enum):
    """Host-data capability required by a solvation component."""

    ELECTROSTATICS = "electrostatics"
    GOSTSHYP = "gostshyp"


class SolvationModelComponent:
    """Immutable model-component configuration backed by a native constructor."""

    coupling_channels: frozenset[CouplingChannel] = frozenset()

    def __init__(self, handle: library.ComponentHandle) -> None:
        self._handle = handle

    def _as_handle(self) -> library.ComponentHandle:
        return self._handle


class PCMSolver(IntEnum):
    """Linear solver shared by PCM-family components."""

    INVERSION = 1
    LU = 2
    CHOLESKY = 3
    ITERATIVE = 4


# Historical name retained for callers that imported it directly.
CPCMSolver = PCMSolver


class _ModelComponentPCMBase(SolvationModelComponent):
    """Shared immutable configuration for PCM-family components."""

    coupling_channels = frozenset({CouplingChannel.ELECTROSTATICS})
    _SOLVERS = {solver.name.lower(): solver for solver in PCMSolver}

    def __init__(
        self,
        epsilon: float,
        solver: str | int | PCMSolver,
        constructor: Callable[[float, int], library.ComponentHandle],
    ) -> None:
        if isinstance(solver, str):
            try:
                solver_value = self._SOLVERS[solver.lower()]
            except KeyError as exc:
                choices = ", ".join(self._SOLVERS)
                raise ValueError(f"Unknown PCM solver {solver!r}; choose {choices}") from exc
        else:
            try:
                solver_value = PCMSolver(int(solver))
            except (TypeError, ValueError) as exc:
                raise ValueError("PCM solver enumeration must be between 1 and 4") from exc
        self._epsilon = float(epsilon)
        self._solver = solver_value
        super().__init__(constructor(self._epsilon, int(self._solver)))

    @property
    def epsilon(self) -> float:
        return self._epsilon

    @property
    def solver(self) -> PCMSolver:
        return self._solver


class ModelComponentCPCM(_ModelComponentPCMBase):
    """Conductor-like polarizable continuum component."""

    def __init__(
        self,
        epsilon: float,
        solver: str | int | PCMSolver = PCMSolver.CHOLESKY,
    ) -> None:
        super().__init__(epsilon, solver, library.new_cpcm_component)


class ModelComponentCOSMO(_ModelComponentPCMBase):
    """Conductor-like screening-model component."""

    def __init__(
        self,
        epsilon: float,
        solver: str | int | PCMSolver = PCMSolver.CHOLESKY,
    ) -> None:
        super().__init__(epsilon, solver, library.new_cosmo_component)


class ModelComponentPV(SolvationModelComponent):
    """Pressure-volume energy component ``pressure * cavity volume``."""

    def __init__(self, pressure: float) -> None:
        self._pressure = float(pressure)
        super().__init__(library.new_pv_component(self._pressure))

    @property
    def pressure(self) -> float:
        return self._pressure


class ModelComponentGOSTSHYP(SolvationModelComponent):
    """GOSTSHYP hydrostatic-pressure component."""

    coupling_channels = frozenset({CouplingChannel.GOSTSHYP})

    def __init__(self, pressure: float) -> None:
        self._pressure = float(pressure)
        super().__init__(library.new_gostshyp_component(self._pressure))

    @property
    def pressure(self) -> float:
        return self._pressure


# -----------------------------------------------------------------------------
# Coupling adapters and evaluations
# -----------------------------------------------------------------------------


class SolvationCoupling(ABC):
    """Adapter between a host representation and one moist evaluation."""

    channels: frozenset[CouplingChannel] = frozenset()

    @property
    @abstractmethod
    def structure(self) -> Structure:
        """Structure associated with this coupling."""

    def activate(self) -> None:
        """Publish adapter state needed by callback-backed cavity construction."""

    @abstractmethod
    def prepare(self, transaction: CouplingTransaction) -> None:
        """Supply host data through one model-owned coupling transaction."""

    def fock(
        self,
        cavity: CavitySnapshot,
        potential: GeneralPotential,
    ) -> Optional[np.ndarray]:
        """Return a host Fock contribution, or ``None`` when unavailable."""
        return None

    def gradient(
        self,
        cavity: CavitySnapshot,
        potential: GeneralPotential,
        model_gradient: Callable[[], np.ndarray],
    ) -> np.ndarray:
        """Return the complete nuclear gradient for this coupling."""
        return model_gradient()


ElectrostaticsProvider = Callable[
    [CavitySnapshot, Optional[TracePotential]], Electrostatics
]


class CouplingTransaction:
    """Restricted interface for supplying host data during one evaluation.

    Custom coupling adapters receive this object instead of the model itself.
    It owns multi-pass native ordering and exposes only the operations valid
    between a successful cavity update and result assembly.
    """

    __slots__ = ("_model", "_cavity")

    def __init__(self, model: SolvationModel) -> None:
        self._model = model
        self._cavity = model.cavity.snapshot()

    @property
    def cavity(self) -> CavitySnapshot:
        return self._cavity

    @property
    def density_dependent(self) -> bool:
        return self._model.cavity.density_dependent

    def requires(self, channel: CouplingChannel) -> bool:
        return channel in self._model.required_coupling_channels

    def exchange_electrostatics(self, provider: ElectrostaticsProvider) -> None:
        """Complete moist's two-pass electrostatic host exchange."""
        self._model._supply_electrostatics(provider(self._cavity, None))
        trace = self._model.trace_potential()
        self._model._supply_electrostatics(provider(self._cavity, trace))

    def supply_gostshyp(self, moments: GostshypMoments) -> None:
        """Supply Gaussian moments for the transaction's current surface."""
        self._model._supply_gostshyp(moments)


class ArrayCoupling(SolvationCoupling):
    """Low-level adapter for hosts that already own the required arrays.

    ``electrostatics`` may be a fixed :class:`Electrostatics` value or a
    callable.  A callable is invoked first with ``trace=None`` and then with the
    model's direct trace potential, allowing charge-dependent response arrays to
    be constructed without exposing the two-pass ordering to the caller.
    """

    def __init__(
        self,
        structure: Structure,
        *,
        electrostatics: Optional[Electrostatics | ElectrostaticsProvider] = None,
        gostshyp: Optional[GostshypMoments] = None,
    ) -> None:
        self._structure = structure
        self._electrostatics = electrostatics
        self._gostshyp = gostshyp
        channels = set()
        if electrostatics is not None:
            channels.add(CouplingChannel.ELECTROSTATICS)
        if gostshyp is not None:
            channels.add(CouplingChannel.GOSTSHYP)
        self.channels = frozenset(channels)

    @property
    def structure(self) -> Structure:
        return self._structure

    def _electrostatic_data(
        self,
        cavity: CavitySnapshot,
        trace: Optional[TracePotential],
    ) -> Electrostatics:
        provider = self._electrostatics
        if provider is None:
            raise RuntimeError("No electrostatics were configured")
        return provider(cavity, trace) if callable(provider) else provider

    def prepare(self, transaction: CouplingTransaction) -> None:
        if self._electrostatics is not None:
            transaction.exchange_electrostatics(self._electrostatic_data)
        if self._gostshyp is not None:
            transaction.supply_gostshyp(self._gostshyp)


class Evaluation:
    """Results from one coherent model/coupling evaluation.

    Energy, potential, cavity and Fock data are captured before the evaluation
    is returned.  The usually more expensive gradient is lazy and may only be
    requested while this remains the model's current evaluation; this prevents
    a later model update from being mixed with an older potential or cavity.
    """

    __slots__ = (
        "_model",
        "_epoch",
        "_coupling",
        "_gradient",
        "_cavity_result",
        "_energy",
        "_potential_result",
        "_fock_result",
    )

    def __init__(
        self,
        *,
        model: SolvationModel,
        epoch: int,
        coupling: SolvationCoupling,
        cavity: CavitySnapshot,
        energy: float,
        potential: GeneralPotential,
        fock: Optional[np.ndarray],
    ) -> None:
        self._model = model
        self._epoch = epoch
        self._coupling = coupling
        self._gradient: Optional[np.ndarray] = None
        self._cavity_result = cavity
        self._energy = float(energy)
        self._potential_result = potential
        self._fock_result = None if fock is None else _immutable_array(fock)

    @property
    def cavity(self) -> CavitySnapshot:
        return self._cavity_result

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def potential(self) -> GeneralPotential:
        return self._potential_result

    @property
    def fock(self) -> Optional[np.ndarray]:
        return self._fock_result

    @property
    def charges(self) -> np.ndarray:
        """Direct molecular-potential adjoints (surface charges for CPCM)."""
        return self.potential.w_umol

    @property
    def gradient(self) -> np.ndarray:
        if self._gradient is None:
            if self._model.epoch != self._epoch:
                raise RuntimeError(
                    "This evaluation was superseded; request its gradient before "
                    "evaluating the model again"
                )
            self._gradient = _immutable_array(
                self._coupling.gradient(
                    self.cavity,
                    self.potential,
                    self._model.gradient,
                )
            )
        return self._gradient


# -----------------------------------------------------------------------------
# Solvation models
# -----------------------------------------------------------------------------


class SolvationModel:
    """Compose one live cavity with an ordered set of solvation components."""

    def __init__(
        self,
        cavity: Cavity,
        components: list[SolvationModelComponent] | tuple[SolvationModelComponent, ...],
        debug: bool = False,
        verbosity: int = 0,
    ) -> None:
        if not isinstance(cavity, Cavity):
            raise TypeError("cavity must be a moist Cavity object")
        items = tuple(components)
        if not items:
            raise ValueError("A solvation model requires at least one component")
        if any(not isinstance(item, SolvationModelComponent) for item in items):
            raise TypeError("components must contain only SolvationModelComponent objects")

        self._updated = False
        self._natoms: Optional[int] = None
        self._epoch = 0
        self._source_cavity = cavity
        self._components = items
        self._required_coupling_channels = frozenset().union(
            *(item.coupling_channels for item in items)
        )
        self._model = library.new_general_model(
            cavity._as_handle(),
            [item._as_handle() for item in items],
            debug=debug,
            verbosity=verbosity,
        )
        borrowed = library.get_model_cavity(self._model)
        self._cavity = cavity._model_view(borrowed)

    @property
    def epoch(self) -> int:
        return self._epoch

    @property
    def components(self) -> tuple[SolvationModelComponent, ...]:
        return self._components

    @property
    def required_coupling_channels(self) -> frozenset[CouplingChannel]:
        return self._required_coupling_channels

    def _invalidate(self) -> None:
        self._epoch += 1
        self._updated = False
        self._natoms = None
        self._cavity._invalidate()

    def update(self, structure: Structure) -> None:
        self._invalidate()
        _guarded_native_update(
            self._source_cavity,
            lambda: library.update_model(self._model, structure._as_handle()),
        )
        self._natoms = len(structure)
        self._updated = True
        self._cavity._mark_updated()

    def _require_updated(self) -> None:
        if not self._updated:
            raise RuntimeError("Model has not been successfully updated")

    @property
    def energy(self) -> float:
        self._require_updated()
        return library.get_model_energy(self._model)

    def get_energy(self) -> float:
        """Compatibility method for :attr:`energy`."""
        return self.energy

    def gradient(self) -> np.ndarray:
        """Return the native model gradient for the last updated structure."""
        self._require_updated()
        if self._natoms is None:
            raise RuntimeError("Model has no updated structure to differentiate")
        return library.general_model_get_gradient(self._model, self._natoms)

    def get_gradient(self, natoms: Optional[int] = None) -> np.ndarray:
        """Compatibility method; ``natoms`` is now inferred from ``update``."""
        if natoms is not None and self._natoms is not None and int(natoms) != self._natoms:
            raise ValueError(
                f"natoms={natoms} does not match the updated structure ({self._natoms})"
            )
        return self.gradient()

    @property
    def cavity(self) -> Cavity:
        """The authoritative model-owned live cavity."""
        return self._cavity

    @property
    def cavity_handle(self) -> library.CavityHandle:
        """Deprecated low-level escape hatch; use :attr:`cavity` instead."""
        warnings.warn(
            "cavity_handle is deprecated; model.cavity exposes the live cavity object",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._cavity._as_handle()

    @property
    def ngrid(self) -> int:
        return self._cavity.ngrid

    def _supply_electrostatics(self, data: Electrostatics) -> None:
        self._require_updated()
        self._epoch += 1
        library.general_model_supply_electrostatics(
            self._model,
            data.phi,
            data.w_xi,
            data.w_f,
            data.w_xyz,
            data.w_n,
            data.qefield,
        )

    def supply_electrostatics(
        self,
        phi: np.ndarray,
        *,
        w_xi: Optional[np.ndarray] = None,
        w_f: Optional[np.ndarray] = None,
        w_xyz: Optional[np.ndarray] = None,
        w_n: Optional[np.ndarray] = None,
        qefield: Optional[np.ndarray] = None,
    ) -> None:
        """Compatibility shim for manually staged electrostatic coupling."""
        self._supply_electrostatics(
            Electrostatics(phi, w_xi, w_f, w_xyz, w_n, qefield)
        )

    def _supply_gostshyp(self, moments: GostshypMoments) -> None:
        self._require_updated()
        self._epoch += 1
        library.general_model_supply_gostshyp(
            self._model,
            moments.gt,
            moments.pt,
            moments.mt,
            moments.rt,
        )

    def supply_gostshyp(
        self,
        gt: np.ndarray,
        pt: np.ndarray,
        mt: np.ndarray,
        rt: np.ndarray,
    ) -> None:
        """Compatibility shim for manually staged GOSTSHYP moments."""
        self._supply_gostshyp(GostshypMoments(gt, pt, mt, rt))

    def trace_potential(self) -> TracePotential:
        self._require_updated()
        molecular, normal = library.general_model_get_trace_potential(
            self._model, self.ngrid
        )
        return TracePotential(molecular, normal)

    def get_trace_potential(self) -> tuple[np.ndarray, np.ndarray]:
        """Compatibility tuple form of :meth:`trace_potential`."""
        trace = self.trace_potential()
        return trace.molecular, trace.normal

    def get_potential(self) -> GeneralPotential:
        """Compatibility read without optional Gaussian response channels."""
        self._require_updated()
        return GeneralPotential(
            **library.general_model_get_potential(self._model, self.ngrid)
        )

    def potential(self) -> GeneralPotential:
        """Return every response channel in one composed value."""
        self._require_updated()
        return GeneralPotential(
            **library.general_model_get_potential_extended(self._model, self.ngrid)
        )

    def get_potential_extended(self) -> GeneralPotential:
        """Compatibility alias for :meth:`potential`."""
        return self.potential()

    def solve(self, phi: np.ndarray) -> tuple[float, np.ndarray]:
        """Compatibility helper for energy/charge-only electrostatic solves."""
        self.supply_electrostatics(phi)
        return self.energy, self.trace_potential().molecular

    def evaluate(
        self,
        structure: Optional[Structure] = None,
        *,
        coupling: Optional[SolvationCoupling] = None,
    ) -> Evaluation:
        """Evaluate the model and host coupling as one coherent transaction."""
        if coupling is None:
            if structure is None:
                raise TypeError("evaluate requires a structure or coupling")
            coupling = ArrayCoupling(structure)
        elif not isinstance(coupling, SolvationCoupling):
            raise TypeError("coupling must implement SolvationCoupling")

        channels = frozenset(CouplingChannel(channel) for channel in coupling.channels)
        missing = self.required_coupling_channels - channels
        if missing:
            names = ", ".join(sorted(channel.value for channel in missing))
            raise ValueError(f"Coupling does not provide required channels: {names}")

        coupling_structure = coupling.structure
        if structure is None:
            structure = coupling_structure
        elif not structure._same_geometry(coupling_structure):
            raise ValueError("structure does not match the coupling structure")

        try:
            coupling.activate()
            self.update(structure)
            coupling.prepare(CouplingTransaction(self))

            energy = self.energy
            potential = self.potential()
            cavity = self._cavity.snapshot()
            fock = coupling.fock(cavity, potential)
            return Evaluation(
                model=self,
                epoch=self.epoch,
                coupling=coupling,
                cavity=cavity,
                energy=energy,
                potential=potential,
                fock=fock,
            )
        except Exception:
            self._invalidate()
            raise


# Compatibility name retained for callers of the pre-refactor general model.
GeneralSolvationModel = SolvationModel
