"""Reusable radius models; radii are in bohr and copied into each cavity."""

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from . import library


@dataclass(frozen=True)
class Radii:
    """Base configuration for a native radius model."""

    _kind: ClassVar[str]

    def _as_handle(self):
        return library.new_radii(self._kind)


@dataclass(frozen=True)
class CPCMRadii(Radii):
    """Default CPCM radii."""

    _kind = "cpcm"


@dataclass(frozen=True)
class SMDRadii(Radii):
    """SMD radii."""

    _kind = "smd"


@dataclass(frozen=True)
class D3Radii(Radii):
    """D3 radii."""

    _kind = "d3"


@dataclass(frozen=True)
class COSMORadii(Radii):
    """COSMO radii."""

    _kind = "cosmo"


@dataclass(frozen=True)
class BondiRadii(Radii):
    """Bondi radii."""

    _kind = "bondi"


@dataclass(frozen=True, init=False)
class CustomRadii(Radii):
    """Positive radii by atom, or by element when ``numbers`` is supplied.

    ``CustomRadii([r_O, r_H, r_H])`` follows molecular atom order.
    ``CustomRadii([r_H, r_O], numbers=[1, 8])`` assigns by atomic number.
    Input arrays are copied to immutable tuples.
    """

    values: tuple[float, ...]
    numbers: tuple[int, ...] | None
    _kind = "custom"

    def __init__(self, values, *, numbers=None):
        array = np.asarray(values, dtype=np.float64)
        if array.ndim != 1 or not array.size:
            raise ValueError("radii must be a nonempty vector")
        if not np.all(np.isfinite(array)) or np.any(array <= 0):
            raise ValueError("radii must be finite and positive")
        elements = None
        if numbers is not None:
            elements = np.asarray(numbers)
            if elements.shape != array.shape or elements.dtype.kind not in "iu":
                raise ValueError("numbers must be an integer vector matching radii")
            if np.any(elements < 1) or np.any(elements > 118):
                raise ValueError("atomic numbers must be between 1 and 118")
            if np.unique(elements).size != elements.size:
                raise ValueError("element radii must have unique atomic numbers")
            elements = tuple(int(value) for value in elements)
        object.__setattr__(self, "values", tuple(float(value) for value in array))
        object.__setattr__(self, "numbers", elements)

    def _as_handle(self):
        handle = super()._as_handle()
        values = np.array(self.values, dtype=np.float64)
        if self.numbers is None:
            library.set_custom_radii_atoms(handle, values)
        else:
            library.set_custom_radii_elements(
                handle, np.array(self.numbers, dtype=np.int32), values
            )
        return handle


def _resolve_radii(radii):
    if radii is None:
        return CPCMRadii()
    if not isinstance(radii, Radii) or type(radii) is Radii:
        raise TypeError("radii must be a radius model such as CPCMRadii or CustomRadii")
    return radii
