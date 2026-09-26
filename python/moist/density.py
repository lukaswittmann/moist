"""Gaussian basis data and live density for MOIST's internal evaluator."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True, init=False)
class GaussianBasis:
    """Cartesian-monomial contracted Gaussian basis.

    Shell atom indices are zero-based. Primitive exponents are in bohr**-2;
    coefficients include primitive normalization. Multiple contractions are
    represented as separate shells. Data are copied into immutable tuples.
    """

    shell_atom: tuple[int, ...]
    shell_l: tuple[int, ...]
    shell_nprim: tuple[int, ...]
    exponents: tuple[float, ...]
    coefficients: tuple[float, ...]

    def __init__(self, *, shell_atom, shell_l, shell_nprim, exponents, coefficients):
        size = None
        for name, values in (("shell_atom", shell_atom), ("shell_l", shell_l),
                             ("shell_nprim", shell_nprim)):
            array = np.asarray(values)
            if array.ndim != 1 or not array.size or array.dtype.kind not in "iu":
                raise ValueError(f"{name} must be a nonempty integer vector")
            if size is not None and array.size != size:
                raise ValueError("shell arrays must have the same length")
            if np.any(array < (1 if name == "shell_nprim" else 0)) or np.any(array > np.iinfo(np.int32).max):
                raise ValueError(f"{name} contains an invalid index or count")
            size = array.size
            object.__setattr__(self, name, tuple(map(int, array)))
        count = sum(self.shell_nprim)
        if count > np.iinfo(np.int32).max:
            raise ValueError("too many primitives for the native interface")
        for name, values in (("exponents", exponents), ("coefficients", coefficients)):
            array = np.asarray(values, dtype=np.float64)
            if array.shape != (count,) or not np.all(np.isfinite(array)):
                raise ValueError(f"{name} must contain {count} finite primitive values")
            if name == "exponents" and np.any(array <= 0):
                raise ValueError("exponents must be positive")
            object.__setattr__(self, name, tuple(map(float, array)))

    @property
    def ncart(self):
        """Number of Cartesian monomials; query cavity layout for their order."""
        return sum((l + 1) * (l + 2) // 2 for l in self.shell_l)


class InternalDensity:
    """Live density source using MOIST's internal Gaussian evaluator.

    The matrix must be expressed in the basis's Cartesian-monomial layout,
    available as ``cavity.isodensity_layout()``. Assign ``density_matrix`` before
    each cavity/model update. Both assignment and reading copy the matrix;
    changing it takes effect at the next explicit update/evaluation.
    """

    def __init__(self, basis: GaussianBasis, density_matrix):
        if not isinstance(basis, GaussianBasis):
            raise TypeError("basis must be a GaussianBasis")
        self._basis = basis
        self.density_matrix = density_matrix

    @property
    def basis(self):
        return self._basis

    @property
    def density_matrix(self):
        return self._density_matrix.copy()

    @density_matrix.setter
    def density_matrix(self, value):
        array = np.asarray(value, dtype=np.float64)
        if array.shape != (self.basis.ncart, self.basis.ncart):
            raise ValueError("density_matrix must have shape (ncart, ncart)")
        if not np.all(np.isfinite(array)):
            raise ValueError("density_matrix must be finite")
        if not np.allclose(array, array.T, rtol=1e-12, atol=1e-14):
            raise ValueError("density_matrix must be symmetric")
        self._density_matrix = np.array(array, order="C", copy=True)
