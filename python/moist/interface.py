"""Object-oriented Python interface for moist solvation models.

The native C interface is intentionally procedural.  This module puts the
Python objects around it: live cavity objects own their behaviour, snapshots
are explicit values, and a model drives the same host coupling protocol as the
Fortran interface -- create a coupling, prepare a phase, walk the requests
with ``for request in coupling`` and answer the missing outputs of each with
``coupling.answer(...)``, then read the result and contract every item of the
response with ``for item in response``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from typing import Callable, ClassVar, Iterable, Iterator, Optional, Protocol, Union

import numpy as np

from . import library
from .library import CavityField
from .configuration import CFC, DROP, ISwiG, Isodensity, SvdW, LevelSet
from .radii import Radii
from .density import InternalDensity
from .parameters import (
    DROPParameters, ISwiGParameters, ModelParameters, PCMParameters, PCMSolver, _resolve,
)


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
        ``(ngrid, 3)`` grid-point coordinates in bohr.
    ``a``
        ``(ngrid,)`` grid-point areas.
    ``owner``
        ``(ngrid,)`` index of the sphere each grid point belongs to.
    ``converged``
        ``(ngrid,)`` per-point projection success flags.
    ``radii``, ``asph``
        ``(nsph,)`` sphere radii and per-sphere surface areas.
    ``xi0``
        ``(ngrid,)`` Gaussian width of each grid point.
    ``f``
        ``(ngrid,)`` Gaussian switching factor of each grid point.
    ``normal0``
        ``(ngrid, 3)`` outward unit normals.

    ``xyz``, ``xi0``, ``normal0``, ``a`` and ``f`` are the grid inputs a host
    evaluates the coupling requests on.
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
    xi0: np.ndarray
    f: np.ndarray
    normal0: np.ndarray


@dataclass(frozen=True)
class CavitySnapshotDROP(CavitySnapshot):
    """DROP-specific cavity data copied from one successful update.

    ``nmax``
        Grid points allocated per sphere before pruning.
    ``wleb``
        ``(ngrid,)`` Lebedev quadrature weights.
    ``r_iI0``
        ``(ngrid,)`` distance from each grid point to its owner sphere centre.
    ``rho``
        ``(ngrid,)`` distance each point was projected from its anchor.
    ``numbering``
        ``(ngrid,)`` stable point ids, packed as ``anchor_id + base*(branch-1)``
        with ``base = nsph * num_leb``.  The same physical point keeps its id
        across rebuilds, which is what makes a surface comparable between
        geometry steps.
    ``anchor_id``, ``branch``, ``branch_count``
        ``(ngrid,)`` the same information unpacked, as moist itself tracks it.
        Points sharing an ``anchor_id`` are branches of one anchor, ``branch``
        is one-based, and ``branch_count`` counts the *surviving* branches in
        the point's anchor group -- so ``branch_count > 1`` marks the branched
        points and a point whose siblings were filtered away reports one.
    ``wbranch``
        ``(ngrid,)`` branch weights, one for an unbranched point.

    Everything the cavity holds beyond these -- curvatures, grid densities,
    projection diagnostics -- is reachable by name through
    :meth:`Cavity.results` and :meth:`Cavity.get`.
    """

    nmax: int
    wleb: np.ndarray
    r_iI0: np.ndarray
    rho: np.ndarray
    numbering: np.ndarray
    anchor_id: np.ndarray
    branch: np.ndarray
    branch_count: np.ndarray
    wbranch: np.ndarray

    @property
    def branched(self) -> np.ndarray:
        """``(ngrid,)`` mask of points sharing an anchor with another point."""
        return self.branch_count > 1


@dataclass(frozen=True)
class AnchorGradient(_ImmutableArrayValue):
    """Anchor-channel nuclear derivatives of a DROP cavity, native grid order.

    Every array is C-contiguous, with reversed native Fortran axes:

    ``xyz1_rA``
        ``(ngrid, nsph, 3, 3)`` indexed ``(i, A, alpha, j)`` --
        ``d(r_i)_j / d(R_A)_alpha``.
    ``xi1_rA``, ``a_i1_rA``, ``v_i1_rA``
        ``(ngrid, nsph, 3)`` indexed ``(i, A, alpha)``.
    ``A_tot1_rA``, ``V_tot1_rA``
        ``(nsph, 3)`` -- the grid sums of ``a_i1_rA``/``v_i1_rA``.
    """

    xyz1_rA: np.ndarray
    xi1_rA: np.ndarray
    a_i1_rA: np.ndarray
    v_i1_rA: np.ndarray
    A_tot1_rA: np.ndarray
    V_tot1_rA: np.ndarray


@dataclass(frozen=True)
class PotentialAdjointResponse(_ImmutableArrayValue):
    """Potential adjoint item of a solvation response.

    ``w_phi``
        ``(ngrid,)`` weights conjugate to the host potential on the grid,
        ``dE/dphi_i``.  The host contracts them as ``F += sum_i w_phi[i] V(r_i)``
        and with the basis-center derivative of ``V`` in the gradient.  For a
        stationary PCM (CPCM, COSMO) they are the induced surface charges
        ``q_i``; for a non-symmetric response matrix (IEF-PCM, SS(V)PE) they
        are the symmetrized adjoint, which is not the apparent charge.
    """

    #: Native item name
    name: ClassVar[str] = "potential_adjoint"

    w_phi: np.ndarray


@dataclass(frozen=True)
class DensityResponse(_ImmutableArrayValue):
    """Density item of a solvation response.

    ``w_rho``, ``w_grad_rho``, ``w_hess_rho``
        ``(ngrid,)``, ``(ngrid, 3)`` and ``(ngrid, 3, 3)`` weights conjugate to
        the solute density and its first two spatial derivatives on the grid.
        The host contracts them with its own ``d rho/dP`` (Fock) or
        ``d rho/dR`` (gradient).  moist folds the ``dS/drho`` factor of its
        level set in, so no level-set convention crosses the boundary.

    The item exists **only** for a cavity whose surface follows the density,
    and **only in the response phase**.  A geometric cavity produces no density
    item at all; that absence is the physics, and it is why a walk over such a
    response never meets one rather than an item of zeros.

    The weights are ``dE/drho`` at fixed nuclei, the same object in either
    phase, so moist forms them once and the gradient phase leaves the
    contraction out -- a host that does not need them would otherwise pay for
    it.  On a density-backed cavity, run the response phase before the gradient
    and reuse them: against your own ``d rho/dP`` they complete the Fock
    matrix, against the basis-centre ``d rho/dR`` the nuclear gradient.  Going
    straight to the gradient phase drops that term and nothing reports it.
    Hessian weights use ``[point,b,a]`` for native ``(a,b,point)``.
    They need not be symmetric; preserve both Cartesian axes in contractions.
    """

    #: Native item name
    name: ClassVar[str] = "density"

    w_rho: np.ndarray
    w_grad_rho: np.ndarray
    w_hess_rho: np.ndarray


@dataclass(frozen=True)
class GostshypAmplitudeResponse(_ImmutableArrayValue):
    """GOSTSHYP amplitude item of a solvation response.

    ``w_overlap``, ``w_normal_deriv``
        ``(ngrid,)`` amplitudes for the Gaussian value and its normal
        derivative; the host completes its Fock contribution as
        ``F += sum_i [w_overlap[i] g[..., i] + w_normal_deriv[i] f[..., i]]``.
    """

    #: Native item name
    name: ClassVar[str] = "gostshyp_amplitude"

    w_overlap: np.ndarray
    w_normal_deriv: np.ndarray


#: One item of a solvation response
_ResponseItem = Union[PotentialAdjointResponse, DensityResponse, GostshypAmplitudeResponse]


class Response:
    """Every item a model phase hands back to the host.

    Each item is the derivative of the model energy with respect to the host
    quantities it names, so the host completes the chain rule by contracting it
    with its own derivative of those quantities.

    Iterating a response yields the items the model produced, in native order:
    :class:`PotentialAdjointResponse`, :class:`DensityResponse` and
    :class:`GostshypAmplitudeResponse`, each with its native item name as
    ``name`` and its arrays as attributes.  Select on the item's class,
    contract every item, and raise on one the host does not know: unlike an
    unanswered request, a skipped item fails nowhere.  An item this model, on
    this cavity, in this phase does not produce is not there -- a fixed cavity
    has no density response, a model without GOSTSHYP no amplitudes -- and
    absence is never filled in with zeros.

    A response is a plain value copied out of the native result of one
    ``get_response``/``get_gradient`` call; it holds no reference to the
    coupling or to the native response.
    """

    __slots__ = ("_items",)

    def __init__(self, items: Iterable[_ResponseItem] = ()) -> None:
        object.__setattr__(self, "_items", tuple(items))

    def __setattr__(self, name: str, value) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __iter__(self) -> Iterator[_ResponseItem]:
        """Yield the items the model produced, in native order."""
        return iter(self._items)

    def __repr__(self) -> str:
        return f"<Response {[item.name for item in self._items]}>"

    @classmethod
    def _from_handle(
        cls, handle: library.ResponseHandle, ngrid: int
    ) -> "Response":
        """Walk the native response and copy every item it carries, each array by name."""
        items = []
        while library.next_response_item(handle):
            name = library.get_response_item_name(handle)
            kind = _RESPONSE_ITEMS.get(name)
            if kind is None:
                raise RuntimeError(f"Unknown response item '{name}' from the native library")
            items.append(kind(**{
                array.name: _read_response_array(handle, array.name, ngrid)
                for array in fields(kind)
            }))
        return cls(items)


#: Native response item name -> the value type holding its arrays
_RESPONSE_ITEMS = {
    kind.name: kind
    for kind in (PotentialAdjointResponse, DensityResponse, GostshypAmplitudeResponse)
}

#: Extents after the grid axis of every response array, keyed by array name
#: (unique across items), as moist.h documents them: moist writes exactly this
#: shape, so the buffer is allocated to it
_RESPONSE_ARRAYS = {
    "w_phi": (),
    "w_rho": (),
    "w_grad_rho": (3,),
    "w_hess_rho": (3, 3),
    "w_overlap": (),
    "w_normal_deriv": (),
}


def _read_response_array(handle: library.ResponseHandle, array: str, ngrid: int) -> np.ndarray:
    """Copy one array of the current native item, ``(ngrid, ...)``."""
    values = np.empty((ngrid, *_RESPONSE_ARRAYS[array]), dtype=np.float64)
    library.get_response_array(handle, array, values)
    return values


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
        if _numbers.dtype.kind not in "iu":
            raise ValueError("numbers must be an integer vector")
        if np.any(_numbers < 1) or np.any(_numbers > 118):
            raise ValueError("atomic numbers must be between 1 and 118")
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
        self._numbers = np.array(_numbers, dtype=np.int32, order="C", copy=True)
        self._positions = _positions
        self._lattice = _lattice
        self._periodic = (
            None if _periodic is None else np.array(_periodic, dtype=np.bool_, order="C", copy=True)
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
        return np.array(array, dtype=np.float64, order="C", copy=True)

    @staticmethod
    def _lattice_array(lattice: Optional[np.ndarray]) -> Optional[np.ndarray]:
        if lattice is None:
            return None
        array = np.asarray(lattice)
        if array.shape != (3, 3):
            raise ValueError("lattice must have shape (3, 3)")
        return np.array(array, dtype=np.float64, order="C", copy=True)

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

    @property
    def configuration(self):
        """Immutable recipe used to construct this cavity."""
        return self._configuration

    @property
    def parameters(self):
        return self.configuration.parameters

    @property
    def radius_model(self):
        return self.configuration.radii

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
        source = getattr(self, "_density_source", None)
        if isinstance(source, InternalDensity):
            library.set_isodensity_density(self._handle, source.density_matrix)
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

    def fields(self) -> tuple[library.CavityField, ...]:
        """Describe every result this cavity currently holds.

        The list is the cavity's own declaration, so it covers whichever type
        this is and omits optional properties that were never computed.  Read
        it again after an update; a rebuild changes the extents.
        """
        self._require_updated()
        return library.list_cavity_fields(self._handle)

    def get(self, name: str) -> np.ndarray:
        """Return one named cavity result, using moist's own name for it.

        Names not carried by :meth:`snapshot` are reachable here -- ``k1``,
        ``KM``, ``rho_grid``, ``phi0`` and the rest.  A field the cavity does
        not currently hold raises; an optional property that was not requested
        is absent rather than zero.

        Each call copies afresh out of the native arrays, so the result is
        writable -- unlike the snapshot's, which back a cached value.
        """
        self._require_updated()
        return library.get_cavity_field(self._handle, name)

    def describe(self, name: str) -> str:
        """Return moist's one-line description of a named result."""
        self._require_updated()
        return library.get_cavity_field_about(self._handle, name)

    def results(self, names: Optional[Iterable[str]] = None) -> dict:
        """Return the cavity's named results as a ``{name: value}`` mapping.

        With ``names`` omitted this reads everything the cavity holds, which is
        the general form of :meth:`snapshot`: the snapshot is a typed value for
        the fields most callers want, this is the whole declared set.
        """
        self._require_updated()
        return library.get_cavity_fields(self._handle, names)

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

    @property
    def xi0(self) -> np.ndarray:
        """``(ngrid,)`` Gaussian width of each grid point, an inverse length in bohr**-1."""
        return self.snapshot().xi0

    @property
    def f(self) -> np.ndarray:
        """``(ngrid,)`` Gaussian switching factor of each grid point."""
        return self.snapshot().f

    @property
    def normal0(self) -> np.ndarray:
        """``(ngrid, 3)`` outward unit normal of each grid point."""
        return self.snapshot().normal0


#: Named fields the generic snapshot reads on top of the generic results
_GRID_INPUT_FIELDS = ("xi0", "f", "normal0")


class _CavityGenericBase(Cavity):
    """Shared implementation for non-DROP native cavities."""

    def _read_snapshot(self) -> CavitySnapshot:
        return CavitySnapshot(
            **library.get_cavity_results(self._handle),
            **library.get_cavity_fields(self._handle, _GRID_INPUT_FIELDS),
        )

    def _model_view(self, handle: library.CavityHandle) -> Cavity:
        return _CavityGenericBorrowed(handle, self)


class _CavityGenericBorrowed(_CavityGenericBase):
    """High-level view of a non-DROP cavity owned by a model."""

    def __init__(self, handle: library.CavityHandle, source: Cavity) -> None:
        super().__init__(handle, owned=False)
        self._configuration = source.configuration
        self._source = source

    @property
    def density_dependent(self) -> bool:
        return self._source.density_dependent


class CavityISwiG(_CavityGenericBase):
    """iSwiG cavity built from shared parameters and a radius model."""

    def __init__(self, *, parameters: ISwiGParameters | None = None,
                 radii: Radii | None = None, **settings) -> None:
        self._configuration = ISwiG(parameters=parameters, radii=radii, **settings)
        super().__init__(library.new_iswig_cavity(
            self.parameters, self.radius_model._as_handle()
        ))


#: DROP fields the typed snapshot carries on top of the generic results.
#: Everything else the cavity declares stays reachable through
#: :meth:`Cavity.results`.
_DROP_SNAPSHOT_FIELDS = _GRID_INPUT_FIELDS + (
    "nmax",
    "wleb",
    "r_iI0",
    "rho",
    "numbering",
    "anchor_id",
    "branch",
    "branch_count",
    "wbranch",
)


class _CavityDROPBase(Cavity):
    """Shared behaviour for standalone and model-owned DROP cavities."""

    @property
    def lsf(self):
        return self.configuration.lsf

    def _read_snapshot(self) -> CavitySnapshotDROP:
        generic = library.get_cavity_results(self._handle)
        drop = library.get_cavity_fields(self._handle, _DROP_SNAPSHOT_FIELDS)
        drop["nmax"] = int(drop["nmax"])
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

    def compute_cavity_gradient(self):
        """Build diagnostic forward derivatives; model gradients do not need this."""
        self._require_updated()
        library.compute_cavity_gradient(self._handle)

    def contract_amat_nuclear_gradient(self, q1, q2):
        """Return a diagnostic forward A-matrix contraction, shape (natoms,3)."""
        self._require_updated()
        return library.contract_amat_nuclear_gradient(self._handle, q1, q2)

    def contract_pcm_nuclear_gradient(self, w_phi, w_xyz, charges):
        """Return a diagnostic forward PCM contraction, shape (natoms,3).

        ``w_phi`` is the potential adjoint ``dE/dphi`` (the surface charge of a
        stationary PCM), ``charges`` the nuclear charges.
        """
        self._require_updated()
        return library.contract_pcm_nuclear_gradient(self._handle, w_phi, w_xyz, charges)

    def isodensity_layout(self):
        """Return shell offsets and monomial powers ``(ncart,3)`` for internal density."""
        return library.isodensity_layout(self._handle)

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
        w_normal: Optional[np.ndarray] = None,
        w_k1: Optional[np.ndarray] = None,
        w_k2: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._require_updated()
        return library.contract_surface_lsf_weights(
            self._handle, w_xi, w_f, w_xyz, w_normal, w_k1, w_k2
        )

    @property
    def nmax(self) -> int:
        return self.snapshot().nmax

    @property
    def wleb(self) -> np.ndarray:
        return self.snapshot().wleb

    @property
    def r_iI0(self) -> np.ndarray:
        return self.snapshot().r_iI0

    @property
    def rho(self) -> np.ndarray:
        return self.snapshot().rho


class _CavityDROPBorrowed(_CavityDROPBase):
    """High-level view of the authoritative cavity owned by a model."""

    def __init__(self, handle: library.CavityHandle, source: Cavity) -> None:
        super().__init__(handle, owned=False)
        self._configuration = source.configuration
        self._source = source  # also keeps callback-backed source objects alive

    @property
    def density_dependent(self) -> bool:
        return self._source.density_dependent


class CavityDROP(_CavityDROPBase):
    """DROP cavity composed from an LSF, parameters and radii.

    Omitted ``lsf`` selects SvdW for compatibility. Density-backed surfaces
    accept a live ``source`` separately from their immutable LSF parameters.
    """

    def __init__(self, *, lsf: LevelSet | None = None,
                 parameters: DROPParameters | None = None, radii: Radii | None = None, source=None,
                 pass_order=None, **settings) -> None:
        self._configuration = DROP(
            lsf=SvdW() if lsf is None else lsf,
            parameters=parameters, radii=radii, **settings,
        )
        self._density_source = source
        super().__init__(self.lsf._new_cavity(
            self.parameters, self.radius_model._as_handle(), source, pass_order
        ))

    @property
    def density_dependent(self):
        return self.lsf.density_dependent

    def _before_native_update(self) -> None:
        state = getattr(self._handle, "callback_state", None)
        if state is not None:
            state.reset()

    def _raise_callback_failure(self) -> None:
        state = getattr(self._handle, "callback_state", None)
        if state is not None:
            state.raise_if_failed()


#: ``(rho, drho)``, ``(rho, drho, d2rho)`` or ``(rho, drho, d2rho, d3rho)``:
#: the tuple grows with the requested order so a caller never pays for a
#: derivative it did not ask for.
DensityDerivatives = tuple[Union[float, np.ndarray], ...]

#: A raw density callback.  It takes ``(point, order)``, or just ``(point)``
#: when the cavity was built with ``pass_order=False``.  It returns the bare
#: electron density: moist builds ``S = scale * (rho_iso - rho)`` from it.
IsodensityCallback = Callable[..., DensityDerivatives]


class IsodensitySource(Protocol):
    """Provider of bare density derivatives, independent of LSF settings."""

    def density(self, point: np.ndarray, order: int) -> DensityDerivatives:
        """Return the density and its spatial derivatives through ``order``."""
        ...


# -----------------------------------------------------------------------------
# Host coupling requests
# -----------------------------------------------------------------------------
#
# A model does not see the host's wavefunction.  It declares which raw
# quantities it needs on the cavity grid, and a host answers whichever outputs
# are missing in the phase it is driving.  One class per request kind, named
# and addressed exactly as in the Fortran interface.  The output names and
# their documented shapes are the contract: the native library reports a wrong
# name by name and trusts the shape, as a pointer carries none; an ndarray
# carries its shape, so :meth:`Coupling.answer` checks it before the call.


class CouplingRequest:
    """Snapshot of one host calculation, taken when the coupling's cursor reached it.

    Iterating a :class:`Coupling` yields one snapshot per request that misses
    an output; its class says which calculation it is.  A snapshot is a plain
    value: it holds no reference to the coupling, answering does not change
    it, and it cannot answer itself -- :meth:`Coupling.answer` answers the
    request the cursor stands at.

    :attr:`missing` names the outputs to compute.
    """

    #: Scientific name of the request kind
    name: str = ""
    #: Output names this kind can answer, in native order, each with its
    #: extents after the grid axis
    _outputs: dict[str, tuple[int, ...]] = {}

    __slots__ = ("_missing",)

    def __init__(self, missing: Iterable[str]) -> None:
        self._missing = frozenset(missing)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} missing={sorted(self._missing)}>"

    @property
    def missing(self) -> frozenset[str]:
        """Outputs required by the staged phase and unanswered when the cursor got here.

        Answer exactly these.  The set is fixed at the snapshot; the next pass
        over the coupling shows what is still missing after the answers.
        """
        return self._missing

    @classmethod
    def _capture(cls, coupling: Coupling) -> CouplingRequest:
        """Snapshot the current request of ``coupling``, which is of this kind."""
        handle = coupling._handle
        missing = frozenset(
            output for output in cls._outputs
            if library.get_coupling_request_missing(handle, output)
        )
        return cls(missing, **cls._inputs(coupling))

    @classmethod
    def _inputs(cls, coupling: Coupling) -> dict:
        """Inputs of this kind beyond the grid, copied from the current request."""
        return {}


class PointPotentialRequest(CouplingRequest):
    """Point potential and its raw spatial derivative on the grid."""
    name = "point_potential"
    _outputs = {"phi": (), "dphi_dr": (3,)}
    __slots__ = ()


class GaussianPotentialRequest(CouplingRequest):
    """Gaussian potential, its spatial derivative and its inverse-length derivative."""
    name = "gaussian_potential"
    _outputs = {"phi": (), "dphi_dr": (3,), "dphi_dxi": ()}
    __slots__ = ()


class GaussianMomentRequest(CouplingRequest):
    """Independently requested s, p, d and contracted f density moments.

    The moments are taken with the Gaussian exponents :attr:`width`, the
    request's own input, copied with the snapshot.
    """
    name = "gaussian_moments"
    _outputs = {"gt": (), "pt": (3,), "mt": (3, 3), "rt": (3,)}
    __slots__ = ("_width",)

    def __init__(self, missing: Iterable[str], width: np.ndarray) -> None:
        super().__init__(missing)
        self._width = _immutable_array(np.asarray(width, dtype=np.float64))

    @property
    def width(self) -> np.ndarray:
        """``(ngrid,)`` Gaussian exponents in bohr**-2, chosen by the component; never recompute them."""
        return self._width

    @classmethod
    def _inputs(cls, coupling: Coupling) -> dict:
        width = np.empty(coupling._model.cavity.ngrid, dtype=np.float64)
        library.get_coupling_request_width(coupling._handle, width)
        return {"width": width}


_REQUEST_CLASSES = {
    cls.name: cls for cls in (PointPotentialRequest, GaussianPotentialRequest, GaussianMomentRequest)
}


class Coupling:
    """A model's requests to the host, with a cursor that walks them.

    Create one with :meth:`SolvationModel.new_coupling`; several may coexist
    on one model, each with its own answers and cursor.  One phase of the host
    loop reads::

        model.prepare_energy(coupling)
        for request in coupling:
            if isinstance(request, GaussianPotentialRequest) and "phi" in request.missing:
                coupling.answer(phi=potential(model.cavity.xyz, model.cavity.xi0))
        energy = np.array(0.0)
        model.get_energy(coupling, energy)

    Iterating drives the native cursor: each step advances it to the next
    request, in declaration order, with an output the staged phase requires
    and has no valid answer for, and yields a snapshot of that request
    (:class:`CouplingRequest`).  Each request is visited at most once per
    pass.  The end of a loop ends the pass and rewinds the cursor, so another
    loop starts a new pass: it retries whatever is still missing, a rejected
    answer for instance, and yields nothing once everything is answered.
    Leaving a loop early, by ``break`` or an exception, leaves the cursor
    where it is, and the next loop resumes the pass after that request.
    ``prepare_*`` and a model update restart the walk; a coupling that is not
    staged yields nothing.  The cursor belongs to the coupling, not to the
    loop, so loops over one coupling share it.

    :meth:`answer` answers the request the cursor stands at, the one the loop
    yielded last.  Grid inputs are read from the model's cavity:
    ``model.cavity.xyz``, ``xi0``, ``normal0``, ``a`` and ``f``.  There is no
    separate completeness check: ``get_*`` fails by name on a wrong staging
    and on every required output still missing.
    """

    def __init__(self, model: SolvationModel) -> None:
        # Private: keeps the model alive and supplies the cavity grid size
        self._model = model
        self._handle = library.new_coupling(model._model)

    def __iter__(self) -> Iterator[CouplingRequest]:
        """Walk the current pass, one snapshot per request with a missing output."""
        while library.next_coupling_request(self._handle):
            yield self._current_kind()._capture(self)

    def _current_kind(self) -> type[CouplingRequest]:
        """Class of the current request; raises when no request is current."""
        name = library.get_coupling_request_name(self._handle)
        try:
            return _REQUEST_CLASSES[name]
        except KeyError:
            raise RuntimeError(f"Unknown coupling request kind {name!r}") from None

    def answer(self, **outputs) -> None:
        """Answer outputs of the current request by keyword; ``None`` skips one.

        Each output is an array ``(ngrid, ...)`` shaped as the request kind
        declares it (see :ref:`coupling-requests`), with ``ngrid`` the grid size
        of ``model.cavity``.  The shape is checked against that declaration
        before the value reaches moist, which reads exactly that many values.
        The outputs are stored one at a time in the kind's output order: a
        rejected output stays missing until a valid retry, those before it keep
        their answers, and those after it are not submitted.  Answering an
        output that already has an answer replaces it.

        :raises RuntimeError: when no request is current -- no loop has
            reached one yet, the pass has ended, or a ``prepare_*`` or model
            update restarted the walk -- or when moist rejects a value, such
            as a non-finite one.
        :raises TypeError: for a keyword the current request does not declare.
        :raises ValueError: for an array that does not convert to float64 or
            does not have the output's declared shape.
        """
        kind = self._current_kind()
        unknown = set(outputs) - set(kind._outputs)
        if unknown:
            raise TypeError(f"{kind.name} has no output {sorted(unknown)[0]!r}")
        ngrid = self._model.cavity.ngrid
        for name, extents in kind._outputs.items():
            value = outputs.get(name)
            if value is None:
                continue
            array = np.ascontiguousarray(value, dtype=np.float64)
            if array.shape != (ngrid, *extents):
                raise ValueError(
                    f"{kind.name}: {name} must have shape {(ngrid, *extents)}, got {array.shape}"
                )
            library.answer_coupling_request(self._handle, name, array)


# -----------------------------------------------------------------------------
# Solvation component configurations
# -----------------------------------------------------------------------------


class SolvationModelComponent:
    """Immutable model-component configuration backed by a native constructor."""

    def __init__(self, handle: library.ComponentHandle) -> None:
        self._handle = handle

    def _as_handle(self) -> library.ComponentHandle:
        return self._handle


class _ModelComponentPCMBase(SolvationModelComponent):
    """PCM physical input and immutable numerical parameters."""

    def __init__(self, epsilon, solver, constructor, parameters) -> None:
        self._parameters = _resolve(PCMParameters, parameters, dict(solver=solver))
        self._epsilon = float(epsilon)
        super().__init__(constructor(self._epsilon, self.parameters))

    @property
    def parameters(self) -> PCMParameters:
        return self._parameters

    @property
    def epsilon(self) -> float:
        return self._epsilon

    @property
    def solver(self) -> PCMSolver:
        return self.parameters.solver


class ModelComponentCPCM(_ModelComponentPCMBase):
    """Conductor-like polarizable continuum component."""

    def __init__(self, epsilon: float, solver=None, *, parameters: PCMParameters | None = None) -> None:
        super().__init__(epsilon, solver, library.new_cpcm_component, parameters)


class ModelComponentCOSMO(_ModelComponentPCMBase):
    """Conductor-like screening-model component."""

    def __init__(self, epsilon: float, solver=None, *, parameters: PCMParameters | None = None) -> None:
        super().__init__(epsilon, solver, library.new_cosmo_component, parameters)


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

    def __init__(self, pressure: float) -> None:
        self._pressure = float(pressure)
        super().__init__(library.new_gostshyp_component(self._pressure))

    @property
    def pressure(self) -> float:
        return self._pressure


# -----------------------------------------------------------------------------
# Solvation models
# -----------------------------------------------------------------------------


class SolvationModel:
    """Compose one live cavity with an ordered set of solvation components."""

    def __init__(
        self,
        cavity: Cavity,
        components: list[SolvationModelComponent] | tuple[SolvationModelComponent, ...],
        debug: Optional[bool] = None,
        verbosity: Optional[int] = None,
        *,
        parameters: Optional[ModelParameters] = None,
    ) -> None:
        if not isinstance(cavity, Cavity):
            raise TypeError("cavity must be a moist Cavity object")
        items = tuple(components)
        if not items:
            raise ValueError("A solvation model requires at least one component")
        if any(not isinstance(item, SolvationModelComponent) for item in items):
            raise TypeError("components must contain only SolvationModelComponent objects")

        self._parameters = _resolve(ModelParameters, parameters, dict(debug=debug, verbosity=verbosity))
        self._updated = False
        self._natoms: Optional[int] = None
        self._source_cavity = cavity
        self._components = items
        #: Reusable native response handle
        self._response: Optional[library.ResponseHandle] = None
        self._model = library.new_general_model(
            cavity._as_handle(),
            [item._as_handle() for item in items],
            parameters=self.parameters,
        )
        borrowed = library.get_model_cavity(self._model)
        self._cavity = cavity._model_view(borrowed)

    @property
    def parameters(self) -> ModelParameters:
        return self._parameters

    @property
    def components(self) -> tuple[SolvationModelComponent, ...]:
        return self._components

    @property
    def cavity(self) -> Cavity:
        """The authoritative model-owned live cavity; the source of every grid input."""
        return self._cavity

    def _invalidate(self) -> None:
        self._updated = False
        self._natoms = None
        self._cavity._invalidate()

    def update(self, structure: Structure) -> None:
        """Rebuild the cavity and the components.

        Invalidates every coupling of the model, even at an unchanged grid
        size: its answers go stale and its walk ends, so prepare it again.
        """
        self._invalidate()
        source = getattr(self._source_cavity, "_density_source", None)
        if isinstance(source, InternalDensity):
            library.set_isodensity_density(self._model, source.density_matrix, model=True)
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

    # ------------------------------------------------------------------
    # host coupling protocol
    # ------------------------------------------------------------------

    def new_coupling(self) -> Coupling:
        """Declare a coupling: every component's requests, answered by the host.

        A coupling belongs to the model and survives :meth:`update`; the
        ``prepare_*`` calls refresh it.  Nothing is staged yet, so iterating
        it yields nothing before the first ``prepare_*``.
        """
        self._require_updated()
        return Coupling(self)

    def _response_handle(self) -> library.ResponseHandle:
        if self._response is None:
            self._response = library.new_response()
        return self._response

    def prepare_energy(self, coupling: Coupling) -> None:
        """Stage the energy phase: mark every answer stale and restart the walk."""
        self._require_updated()
        library.prepare_model_energy(self._model, coupling._handle)

    def prepare_response(self, coupling: Coupling) -> None:
        """Stage the response phase, retaining valid answers; restarts the walk."""
        self._require_updated()
        library.prepare_model_response(self._model, coupling._handle)

    def prepare_gradient(self, coupling: Coupling) -> None:
        """Stage the gradient phase, retaining valid answers; restarts the walk."""
        self._require_updated()
        library.prepare_model_gradient(self._model, coupling._handle)

    def get_energy(self, coupling: Coupling, energy: np.ndarray) -> None:
        """Accumulate the staged energy into a writable float64 scalar array.

        Fails by name on a coupling not staged by :meth:`prepare_energy` and on
        every required output the host left unanswered, including one whose
        answer was rejected; the accumulator is then unchanged.
        """
        self._require_updated()
        library.get_model_energy(self._model, coupling._handle, energy)

    def get_response(self, coupling: Coupling) -> Response:
        """Host part of a staged response phase.

        Exactly what this model produces in this phase and nothing else: the
        potential adjoint, the density weights of a density-backed cavity and
        the GOSTSHYP amplitudes.  An item this configuration does not produce
        is absent, never zero.
        """
        self._require_updated()
        handle = self._response_handle()
        library.get_model_response(self._model, coupling._handle, handle)
        return Response._from_handle(handle, self._cavity.ngrid)

    def get_gradient(self, coupling: Coupling, gradient: np.ndarray) -> Response:
        """Nuclear gradient of a staged gradient phase, plus its host part.

        Adds the model contribution to ``gradient[natoms, 3]`` and returns the
        response the host contracts with its own geometry derivatives.  The
        accumulator must be a writable C-contiguous float64 array.
        """
        self._require_updated()
        if self._natoms is None:
            raise RuntimeError("Model has no updated structure to differentiate")
        handle = self._response_handle()
        library.get_model_gradient(
            self._model, coupling._handle, handle, self._natoms, gradient
        )
        return Response._from_handle(handle, self._cavity.ngrid)
