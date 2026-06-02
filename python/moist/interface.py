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

    @property
    def cavity(self) -> Cavity:
        if not self._updated:
            raise RuntimeError("Model has to be updated before requesting the cavity")

        cavity_handle = library.get_model_cavity(self._model)
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
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return LSF adjoint weights for a contracted drop surface response."""
        if not self._updated:
            raise RuntimeError("Cavity has to be updated before contracting LSF weights")
        return library.contract_surface_lsf_weights(self._cavity, w_xi, w_f, w_xyz)


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
    """DROP cavity whose level set is provided by a Python callback."""

    def __init__(
        self,
        callback,
        nleb: Optional[int] = None,
        scale: float = 1000.0,
        debug: bool = False,
        verbosity: int = 0,
        do_fine: bool = False,
        wleb_prune_level: Optional[int] = None,
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
        )


# class SolvationModelGEMS(SolvationModel):
#     """Minimal Python wrapper for the moist GEMS solvation model."""

#     def __init__(
#         self,
#         solvent: str,
#         debug: bool = False,
#         verbosity: int = 0,
#         parameter_file: Optional[str] = None,
#     ):
#         super().__init__()
#         self._model = library.new_gems_model(
#             solvent=solvent,
#             debug=debug,
#             verbosity=verbosity,
#             parameter_file=parameter_file,
#         )

#     @classmethod
#     def _from_constructor(
#         cls,
#         solvent: str,
#         debug: bool = False,
#         verbosity: int = 0,
#         parameter_file: Optional[str] = None,
#     ):
#         return cls(
#             solvent=solvent,
#             debug=debug,
#             verbosity=verbosity,
#             parameter_file=parameter_file,
#         )
