"""Host-independent surface configurations shared with PySCF.

A configuration holds immutable settings only, never a host density or
native cavity. Calling ``build`` creates independent live state.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar

from . import library
from .density import GaussianBasis
from .parameters import (
    CFCParameters, DROPParameters, ISwiGParameters, IsodensityParameters,
    SvdWParameters, _resolve,
)
from .radii import Radii, _resolve_radii


class LevelSet(ABC):
    """Immutable surface definition; subclasses construct a native cavity."""

    density_dependent = False

    @abstractmethod
    def _new_cavity(self, parameters, radii, source, pass_order, basis, *, context):
        """Create independently owned native state for this surface."""


class _GeometricLevelSet(LevelSet):
    def _new_cavity(self, parameters, radii, source, pass_order, basis, *, context):
        if source is not None or pass_order is not None or basis is not None:
            raise TypeError("A geometric LSF does not accept a density source or basis")
        return library.new_drop_cavity(parameters, self.parameters, radii, context=context)


@dataclass(frozen=True, init=False)
class SvdW(_GeometricLevelSet):
    """Smooth van der Waals surface."""

    parameters: SvdWParameters
    name: ClassVar[str] = "SvdW"

    def __init__(self, *, parameters: SvdWParameters | None = None):
        object.__setattr__(self, "parameters", _resolve(SvdWParameters, parameters))


@dataclass(frozen=True, init=False)
class CFC(_GeometricLevelSet):
    """COSMO fine cavity surface."""

    parameters: CFCParameters
    name: ClassVar[str] = "CFC"

    def __init__(self, *, parameters: CFCParameters | None = None):
        object.__setattr__(self, "parameters", _resolve(CFCParameters, parameters))


@dataclass(frozen=True, init=False)
class Isodensity(LevelSet):
    """Electron-density surface, evaluated internally or by a callback.

    Built with ``basis=``, a :class:`GaussianBasis`, moist evaluates the
    density itself; the host answers the density matrix ``dcart`` in the build
    phase of every coupling walk. Built with ``source=``, a callable or an
    object supplying ``density(point, order)``, moist pulls the bare density
    point by point during the build, outside the coupling. With the core fill
    on, the build phase also asks for the nuclear charges ``z``. The LSF
    parameters define the contour and scale.
    """

    parameters: IsodensityParameters
    name: ClassVar[str] = "Isodensity"
    density_dependent: ClassVar[bool] = True

    def __init__(self, *, parameters: IsodensityParameters | None = None):
        object.__setattr__(self, "parameters", _resolve(IsodensityParameters, parameters))

    def _new_cavity(self, parameters, radii, source, pass_order, basis, *, context):
        if basis is not None:
            if source is not None:
                raise TypeError("Isodensity takes a basis or a callback source, not both")
            if not isinstance(basis, GaussianBasis):
                raise TypeError("basis must be a GaussianBasis")
            if pass_order is not None:
                raise TypeError("pass_order applies only to callbacks")
            return library.new_internal_isodensity_cavity(basis, parameters, self.parameters, radii, context=context)
        callback = getattr(source, "density", source)
        if not callable(callback):
            raise TypeError("Isodensity requires a GaussianBasis, a callable source or density(point, order)")
        return library.new_drop_cavity_isodensity_callback(
            callback, parameters, self.parameters, radii, pass_order=pass_order, context=context
        )


class CavityConfiguration(ABC):
    """Reusable recipe for independent live cavities."""

    density_dependent = False

    @abstractmethod
    def build(self, *, context=None, source=None, basis=None):
        """Create a live cavity; an isodensity surface takes its basis or callback.

        Without ``context`` the cavity runs on the context of its model.
        """


@dataclass(frozen=True, init=False)
class DROP(CavityConfiguration):
    """DROP discretization of an explicit LSF, with reusable parameters."""

    lsf: LevelSet
    parameters: DROPParameters
    radii: Radii
    name: ClassVar[str] = "DROP"

    def __init__(self, *, lsf: LevelSet, parameters: DROPParameters | None = None,
                 radii: Radii | None = None):
        if not isinstance(lsf, LevelSet):
            raise TypeError("DROP requires an LSF such as SvdW, CFC or Isodensity")
        object.__setattr__(self, "lsf", lsf)
        object.__setattr__(self, "parameters", _resolve(DROPParameters, parameters))
        object.__setattr__(self, "radii", _resolve_radii(radii))

    @property
    def density_dependent(self):
        return self.lsf.density_dependent

    def build(self, *, context=None, source=None, pass_order=None, basis=None):
        from .interface import CavityDROP

        return CavityDROP(lsf=self.lsf, parameters=self.parameters, radii=self.radii,
                          source=source, pass_order=pass_order, basis=basis, context=context)


@dataclass(frozen=True, init=False)
class ISwiG(CavityConfiguration):
    """Improved switching-Gaussian discretization."""

    parameters: ISwiGParameters
    radii: Radii
    name: ClassVar[str] = "ISwiG"
    lsf: ClassVar[None] = None
    density_dependent: ClassVar[bool] = False

    def __init__(self, *, parameters: ISwiGParameters | None = None,
                 radii: Radii | None = None):
        object.__setattr__(self, "parameters", _resolve(ISwiGParameters, parameters))
        object.__setattr__(self, "radii", _resolve_radii(radii))

    def build(self, *, context=None, source=None, basis=None):
        from .interface import CavityISwiG

        if source is not None or basis is not None:
            raise TypeError("ISwiG does not accept a density source or basis")
        return CavityISwiG(parameters=self.parameters, radii=self.radii, context=context)
