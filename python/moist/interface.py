"""High-level Python interface for moist solvation models."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import library


@dataclass
class Cavity:
    """Snapshot of cavity data after a model update."""

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
    nmax: int
    normal0: np.ndarray
    wleb: np.ndarray
    r_iI0: np.ndarray
    f: np.ndarray
    rho: np.ndarray


@dataclass
class AnchorGradient:
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


class Structure:
    """Wrapped molecular structure object."""

    _mol = library.StructureHandle.null()

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        if positions.size % 3 != 0:
            raise ValueError("Expected tripels of cartesian coordinates")
        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype=np.int32)
        _positions = np.ascontiguousarray(positions, dtype=np.float64)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype=np.float64)
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype=np.bool_)
        else:
            _periodic = None

        self._mol = library.new_structure(
            self._natoms,
            _numbers,
            _positions,
            _lattice,
            _periodic,
        )

    def __len__(self):
        return self._natoms

    def update(self, positions: np.ndarray, lattice: Optional[np.ndarray] = None) -> None:
        if 3 * len(self) != positions.size:
            raise ValueError("Dimension missmatch for positions")

        _positions = np.ascontiguousarray(positions, dtype=np.float64)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype=np.float64)
        else:
            _lattice = None

        library.update_structure(self._mol, _positions, _lattice)


class SolvationModel(ABC):
    """Shared high-level interface for moist solvation models."""

    _model = library.ModelHandle.null()

    def __init__(self):
        self._updated = False

    @classmethod
    @abstractmethod
    def _from_constructor(cls, *args, **kwargs):
        """Create a model instance from a model-specific constructor."""

    def update(self, structure: Structure) -> None:
        library.update_model(self._model, structure._mol)
        self._updated = True

    def get_energy(self) -> float:
        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting the energy")
        return library.get_model_energy(self._model)

    @staticmethod
    def _cavity_from_handle(cavity_handle) -> Cavity:
        """Snapshot a (borrowed) cavity handle into a :class:`Cavity`."""

        generic = library.get_cavity_results(cavity_handle)
        drop = library.get_drop_specific(cavity_handle, ngrid=generic["ngrid"])

        return Cavity(
            area=generic["area"],
            volume=generic["volume"],
            ngrid=generic["ngrid"],
            nsph=generic["nsph"],
            xyz=generic["xyz"],
            a=generic["a"],
            owner=generic["owner"],
            converged=generic["converged"],
            radii=generic["radii"],
            asph=generic["asph"],
            nmax=drop["nmax"],
            normal0=drop["normal0"],
            wleb=drop["wleb"],
            r_iI0=drop["r_iI0"],
            f=drop["f"],
            rho=drop["rho"],
        )

    @property
    def cavity(self) -> Cavity:
        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting the cavity")

        return self._cavity_from_handle(library.get_model_cavity(self._model))


class SolvationComponent(ABC):
    """Base class for components accepted by :class:`GeneralSolvationModel`."""

    _component = library.ComponentHandle.null()


class CPCM(SolvationComponent):
    """Conductor-like polarizable continuum component.

    Parameters
    ----------
    epsilon
        Relative dielectric constant of the solvent.
    solver
        One of ``"inversion"``, ``"lu"``, ``"cholesky"``, or
        ``"iterative"`` (or the corresponding integer enumeration 1--4,
        exported to C as ``moist_pcm_solver``).
    """

    _SOLVERS = {
        "inversion": 1,
        "lu": 2,
        "cholesky": 3,
        "iterative": 4,
    }

    def __init__(self, epsilon: float, solver: str | int = "cholesky") -> None:
        if isinstance(solver, str):
            try:
                solver_id = self._SOLVERS[solver.lower()]
            except KeyError as exc:
                choices = ", ".join(self._SOLVERS)
                raise ValueError(f"Unknown CPCM solver {solver!r}; choose {choices}") from exc
        else:
            solver_id = int(solver)
            if solver_id not in self._SOLVERS.values():
                raise ValueError("CPCM solver enumeration must be between 1 and 4")
        self.epsilon = float(epsilon)
        self.solver = solver_id
        self._component = library.new_cpcm_component(self.epsilon, self.solver)


class PV(SolvationComponent):
    """Pressure-volume energy component ``pressure * cavity volume``."""

    def __init__(self, pressure: float) -> None:
        self.pressure = float(pressure)
        self._component = library.new_pv_component(self.pressure)


@dataclass
class GeneralPotential:
    """Adjoint channels returned by a general solvation model."""

    w_umol: np.ndarray
    w_qmol: np.ndarray
    w_lsf0: np.ndarray
    w_lsf1: np.ndarray
    w_lsf2: np.ndarray


class GeneralSolvationModel(SolvationModel):
    """Compose one live cavity with an ordered list of solvation components."""

    def __init__(
        self,
        cavity,
        components: list[SolvationComponent] | tuple[SolvationComponent, ...],
        debug: bool = False,
        verbosity: int = 0,
    ) -> None:
        super().__init__()
        if not hasattr(cavity, "_cavity"):
            raise TypeError("cavity must be a moist standalone cavity object")
        items = list(components)
        if not items:
            raise ValueError("A general solvation model requires at least one component")
        if any(not isinstance(item, SolvationComponent) for item in items):
            raise TypeError("components must contain only SolvationComponent objects")

        self._source_cavity = cavity
        self.components = tuple(items)
        self._model = library.new_general_model(
            cavity._cavity,
            [item._component for item in items],
            debug=debug,
            verbosity=verbosity,
        )
        # General models own a cavity copy. Keep its borrowed handle so hosts
        # can install mutable isodensity data before the first geometry update.
        self._cavity_handle = library.get_model_cavity(self._model)

    @classmethod
    def _from_constructor(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    @property
    def cavity_handle(self) -> library.CavityHandle:
        """Borrowed handle to the model-owned live cavity."""

        return self._cavity_handle

    @property
    def ngrid(self) -> int:
        """Number of points on the current live cavity."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting its grid size")
        return library.get_cavity_sizes(self._cavity_handle)[0]

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
        """Supply the molecular potential and optional direct response arrays."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before supplying electrostatics")
        library.general_model_supply_electrostatics(
            self._model, phi, w_xi, w_f, w_xyz, w_n, qefield
        )

    def solve(self, phi: np.ndarray) -> tuple[float, np.ndarray]:
        """Solve all electrostatic components and return energy and charges."""

        self.supply_electrostatics(phi)
        charges, _ = library.general_model_get_trace_potential(self._model, self.ngrid)
        return self.get_energy(), charges

    def get_trace_potential(self) -> tuple[np.ndarray, np.ndarray]:
        """Return direct molecular-potential and normal-trace adjoints."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting a potential")
        return library.general_model_get_trace_potential(self._model, self.ngrid)

    def get_potential(self) -> GeneralPotential:
        """Return accumulated direct and cavity-response adjoint channels."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting a potential")
        result = library.general_model_get_potential(self._model, self.ngrid)
        return GeneralPotential(**result)

    def get_gradient(self, natoms: int) -> np.ndarray:
        """Return the accumulated nuclear gradient with shape ``(3, natoms)``."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting a gradient")
        return library.general_model_get_gradient(self._model, int(natoms))

    @property
    def cavity(self) -> Cavity:
        """Snapshot the current model-owned cavity."""

        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting the cavity")
        return self._cavity_from_handle(self._cavity_handle)


class _DROPCavityBase:
    """Shared standalone DROP cavity result handling."""

    def update(self, structure: Structure) -> None:
        library.update_cavity(self._cavity, structure._mol)
        self._updated = True

    @property
    def cavity(self) -> Cavity:
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before requesting results")
        generic = library.get_cavity_results(self._cavity)
        drop = library.get_drop_specific(self._cavity, ngrid=generic["ngrid"])

        return Cavity(
            area=generic["area"],
            volume=generic["volume"],
            ngrid=generic["ngrid"],
            nsph=generic["nsph"],
            xyz=generic["xyz"],
            a=generic["a"],
            owner=generic["owner"],
            converged=generic["converged"],
            radii=generic["radii"],
            asph=generic["asph"],
            nmax=drop["nmax"],
            normal0=drop["normal0"],
            wleb=drop["wleb"],
            r_iI0=drop["r_iI0"],
            f=drop["f"],
            rho=drop["rho"],
        )

    def assemble_amat(self) -> tuple[np.ndarray, np.ndarray]:
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before assembling the A-matrix")
        return library.assemble_drop_amat(self._cavity)

    def get_gaussian(self) -> tuple[np.ndarray, np.ndarray]:
        """Return Gaussian widths and switching factors without assembling A."""
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before reading Gaussian data")
        return library.get_cavity_gaussian(self._cavity)

    def compute_anchor_gradient(self) -> None:
        """Compute the anchor-only nuclear derivatives (DROP cavities only).

        Must precede :meth:`get_anchor_gradient`.  For a callback level set the
        field's nuclear partials are zero, so the anchor pass is the whole
        nuclear route moist can see on its own.
        """
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before computing the anchor gradient")
        library.compute_anchor_gradient(self._cavity)

    def get_anchor_gradient(self) -> AnchorGradient:
        """Return the anchor-channel nuclear derivatives in native grid order.

        Requires a preceding :meth:`compute_anchor_gradient`.  The grid point area
        carries a switching-function dependence (``a_i ~ f_i / xi_i**2``), so
        ``a_i1_rA`` is the area route a geometric surface functional needs and
        is not recoverable from ``xi1_rA``.
        """
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before reading the anchor gradient")
        return AnchorGradient(**library.get_anchor_gradient(self._cavity))

    def contract_amat_surface_weights(
        self,
        q1: np.ndarray,
        q2: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before contracting A-matrix weights")
        return library.contract_amat1_q1q2_surface_weights(self._cavity, q1, q2)

    def contract_surface_lsf_weights(
        self,
        w_xi: np.ndarray,
        w_f: np.ndarray,
        w_xyz: np.ndarray,
        w_n: Optional[np.ndarray] = None,
        w_k1: Optional[np.ndarray] = None,
        w_k2: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return LSF adjoint weights for a contracted drop surface response.

        The outward-normal (``w_n``) and principal-curvature (``w_k1``,
        ``w_k2``) channels are optional and skipped when omitted.
        """
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before contracting LSF weights")
        return library.contract_surface_lsf_weights(
            self._cavity, w_xi, w_f, w_xyz, w_n, w_k1, w_k2)


class DROPCavity(_DROPCavityBase):
    """Standard solute-vdW (SvdW) DROP cavity with a density-independent surface.

    The surface is built from atomic van-der-Waals spheres (default CPCM radii)
    and depends only on the molecular structure, so the cavity geometry is fixed
    across an SCF.  Shares the same result/A-matrix accessors as the isodensity
    cavity through :class:`_DROPCavityBase`.
    """

    def __init__(
        self,
        nleb: Optional[int] = None,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
    ):
        self._updated = False
        self._cavity = library.new_drop_cavity(
            nleb=nleb,
            debug=debug,
            verbosity=verbosity,
            do_fine=do_fine,
        )


class IsodensityDROPCavity(_DROPCavityBase):
    """DROP cavity whose level set is provided by a Python callback.

    The callback is invoked either as ``callback(point, order)`` -- where
    ``order`` is the highest derivative moist currently needs, so the callback
    can skip computing the rest -- or as ``callback(point)``. The form is
    detected from the signature; set ``pass_order`` to force it. See
    :func:`moist.library.new_drop_cavity_isodensity_callback` for the full
    contract.

    A callback that raises aborts the build: :meth:`update` fails with the
    original exception instead of returning a cavity built on substitute values.
    """

    def __init__(
        self,
        callback,
        nleb: Optional[int] = None,
        scale: float = 1000.0,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
        wleb_prune_level: Optional[int] = None,
        tolerance: Optional[float] = None,
        pass_order: Optional[bool] = None,
    ):
        self._updated = False
        self._cavity, self._callback_ref = library.new_drop_cavity_isodensity_callback(
            callback=callback,
            nleb=nleb,
            scale=scale,
            debug=debug,
            verbosity=verbosity,
            do_fine=do_fine,
            wleb_prune_level=wleb_prune_level,
            tolerance=tolerance,
            pass_order=pass_order,
        )

    def update(self, structure: Structure) -> None:
        """Build the cavity, surfacing any failure inside the Python callback.

        An exception in the callback aborts the build: the binding reports it
        through the callback's return code and moist fails this call with an
        "external LSF evaluation failed" error. That error is accurate but says
        nothing about *why* the callback failed, so the exception itself -- which
        cannot travel back through the C frame -- is carried out of band and
        re-raised here, with its original traceback, in place of moist's.
        """
        state = self._cavity.callback_state
        state.reset()
        try:
            super().update(structure)
        except Exception:
            state.raise_if_failed()
            raise
        state.raise_if_failed()
