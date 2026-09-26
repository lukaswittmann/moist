"""Immutable settings shared by the Python API and host integrations.

Defaults come from the linked MOIST library. Use ``dataclasses.replace`` to
derive a new configuration; changing settings never mutates a live object.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from enum import IntEnum
from functools import lru_cache
import math
from typing import ClassVar

from . import library


class PCMSolver(IntEnum):
    """PCM solver values from the C ABI."""

    INVERSION = library.lib.moist_pcm_solver_inversion
    LU = library.lib.moist_pcm_solver_lu
    CHOLESKY = library.lib.moist_pcm_solver_cholesky
    ITERATIVE = library.lib.moist_pcm_solver_iterative


@lru_cache(maxsize=None)
def _default(kind, name):
    return getattr(library._options(kind), name)


def _native_default(kind, name):
    return field(default_factory=lambda: _default(kind, name))


class _Parameters:
    _kind: ClassVar[str]

    def __post_init__(self):
        # Validate the ABI representation only; the native constructor checks
        # scientific constraints.
        options = self._as_options()
        for item in fields(self):
            value = getattr(options, item.name)
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"{item.name} must be finite")
            object.__setattr__(self, item.name, value)

    def _as_options(self):
        return library._options(self._kind, **asdict(self))


@dataclass(frozen=True, kw_only=True)
class DROPParameters(_Parameters):
    """DROP discretization and projection settings, for every level set.

    Lengths are in bohr. ``do_fine`` enables optional surface diagnostics.
    Field names follow ``moist_drop_options``; derived Fortran settings stay native.
    """

    _kind = "drop"
    nleb: int = _native_default("drop", "nleb")
    debug: bool = _native_default("drop", "debug")
    verbosity: int = _native_default("drop", "verbosity")
    do_fine: bool = _native_default("drop", "do_fine")
    tolerance: float = _native_default("drop", "tolerance")
    proj_maxiter: int = _native_default("drop", "proj_maxiter")
    proj_level: int = _native_default("drop", "proj_level")
    branch_weight_s: float = _native_default("drop", "branch_weight_s")
    rho_grid_h: float = _native_default("drop", "rho_grid_h")
    wleb_prune_level: int = _native_default("drop", "wleb_prune_level")


@dataclass(frozen=True, kw_only=True)
class ISwiGParameters(_Parameters):
    """iSwiG grid settings; positive ``cut_a`` selects the area cutoff."""

    _kind = "iswig"
    nleb: int = _native_default("iswig", "nleb")
    debug: bool = _native_default("iswig", "debug")
    verbosity: int = _native_default("iswig", "verbosity")
    cut_a: float = _native_default("iswig", "cut_a")
    cut_f: float = _native_default("iswig", "cut_f")


@dataclass(frozen=True, kw_only=True)
class SvdWParameters(_Parameters):
    """Smooth van der Waals level-set settings."""

    _kind = "svdw"
    blend_k: float = _native_default("svdw", "blend_k")
    blend_1b: float = _native_default("svdw", "blend_1b")
    blend_2b: float = _native_default("svdw", "blend_2b")
    blend_3b: float = _native_default("svdw", "blend_3b")


@dataclass(frozen=True, kw_only=True)
class CFCParameters(_Parameters):
    """COSMO fine cavity level-set settings."""

    _kind = "cfc"
    a1: float = _native_default("cfc", "a1")
    a2: float = _native_default("cfc", "a2")
    c: float = _native_default("cfc", "c")
    m: int = _native_default("cfc", "m")


@dataclass(frozen=True, kw_only=True)
class IsodensityParameters(_Parameters):
    """Level set ``scale * (rho_iso - rho)``; density in electrons/bohr**3."""

    _kind = "isodensity"
    rho_iso: float = _native_default("isodensity", "rho_iso")
    scale: float = _native_default("isodensity", "scale")


@dataclass(frozen=True, kw_only=True)
class PCMParameters(_Parameters):
    """Numerical settings shared by CPCM and COSMO.

    ``solver_tol`` and ``solver_maxiter`` are read only by the iterative solver.
    """

    _kind = "pcm"
    solver: PCMSolver = _native_default("pcm", "solver")
    solver_tol: float = _native_default("pcm", "solver_tol")
    solver_maxiter: int = _native_default("pcm", "solver_maxiter")

    def __post_init__(self):
        value = self.solver
        if isinstance(value, str):
            try:
                value = PCMSolver[value.upper()]
            except KeyError as exc:
                raise ValueError(f"Unknown PCM solver {value!r}") from exc
        object.__setattr__(self, "solver", PCMSolver(value))
        super().__post_init__()
        object.__setattr__(self, "solver", PCMSolver(self.solver))


@dataclass(frozen=True, kw_only=True)
class ModelParameters(_Parameters):
    """Logging settings for a composed solvation model."""

    _kind = "model"
    debug: bool = _native_default("model", "debug")
    verbosity: int = _native_default("model", "verbosity")


def _resolve(parameter_type, parameters):
    """Return ``parameters``, or the defaults of ``parameter_type`` for ``None``."""
    if parameters is None:
        return parameter_type()
    if not isinstance(parameters, parameter_type):
        raise TypeError(f"parameters must be {parameter_type.__name__}")
    return parameters
