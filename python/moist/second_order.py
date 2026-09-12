"""Model-owned coupled second derivatives, independent of the electronic host.

Host parameters put nuclear Cartesian coordinates first, then independent
coordinates of the density. Surface and component implementations own the
chain rule from those parameters to the relaxed solvation energy.
"""

from typing import Protocol

import numpy as np


class HostDerivatives(Protocol):
    """Analytic host data consumed during one model linearization.

    ``dirs`` has shape (3, natoms, n). ``level_set_jets`` yields fixed-point
    level-set spatial jets through fourth order and their parameter partials.
    ``potential`` returns first and second electrostatic potential derivatives,
    including surface motion, in the supplied parameter coordinates.
    """

    n: int
    nr: int
    dirs: np.ndarray

    def level_set_jets(self, coords): ...
    def electrostatic_potential(self, coords): ...
    def potential(self, point, surface1, surface2): ...


class CoupledResponse:
    """Solvent nuclear/mixed derivatives and a density-response operation.

    ``nuclear`` is (nr, nr), ``mixed`` is (nr, n-nr). Density coordinates and
    their conjugate derivatives use the host's convention. ``density_response``
    accepts (..., n-nr), returning the same shape. Internal solvent relaxation
    is included; electronic relaxation belongs to the host's response solve.
    """

    def __init__(self, nr, direct, factors=()):
        # E'' = direct - sum(left.T @ right). Keep charge-response factors
        # separate so CPHF does not require an assembled density-density block.
        self._direct = np.array(direct, copy=True)
        self._factors = tuple((np.array(a, copy=True), np.array(b, copy=True)) for a, b in factors)
        self._nr = nr
        self.nuclear = self._direct[:nr, :nr].copy()
        self.mixed = self._direct[:nr, nr:].copy()
        for left, right in self._factors:
            self.nuclear -= np.einsum('ip,iq->pq', left[:, :nr], right[:, :nr])
            self.mixed -= np.einsum('ip,iq->pq', left[:, :nr], right[:, nr:])
        for array in (self._direct, self.nuclear, self.mixed):
            array.flags.writeable = False

    def density_response(self, direction):
        direction = np.asarray(direction)
        nr = self._nr
        if direction.ndim == 0 or direction.shape[-1] != len(self._direct) - nr:
            raise ValueError('Density direction has the wrong parameter dimension')
        result = np.einsum('pq,...q->...p', self._direct[nr:, nr:], direction)
        for left, right in self._factors:
            response = np.einsum('iq,...q->...i', right[:, nr:], direction)
            result -= np.einsum('ip,...i->...p', left[:, nr:], response)
        return result


class SecondOrderTransaction:
    """Model-owned dispatch and shared intermediates for component responses.

    Components implement ``second_order(transaction)``. Unsupported components
    fail explicitly rather than silently dropping a contribution. This dense
    reference backend stores surface second derivatives; the public response
    contract also permits implementations that compute directional derivatives.
    """

    def __init__(self, model, parameters: HostDerivatives, max_memory):
        self.model = model
        self.parameters = parameters
        self.coords = model.cavity.snapshot().xyz.T
        n = parameters.n
        if not 0 <= parameters.nr <= n:
            raise ValueError('Invalid nuclear parameter count')
        required_mb = 8 * (7 * len(self.coords) * n * n + 12 * n * n) / 1e6
        if required_mb > max_memory:
            raise MemoryError(f'Dense solvent response needs about {required_mb:.0f} MB; increase max_memory')
        self._surface = None
        self._pcm = None

    @property
    def surface(self):
        """Shared first/second surface partials, computed by the cavity."""
        if self._surface is None:
            self._surface = self.model.cavity.parameter_surface_derivatives(self.parameters)
        return self._surface

    def assemble(self):
        n = self.parameters.n
        direct = np.zeros((n, n))
        factors = []
        for component in self.model.components:
            term, response = component.second_order(self)
            direct += term
            factors.extend(response)
        return CoupledResponse(self.parameters.nr, direct, factors)

    def pcm(self, dielectric_factor):
        """PCM-family contribution with its component-specific dielectric factor."""
        from scipy.linalg import cho_factor, cho_solve

        n = self.parameters.n
        if dielectric_factor == 0:
            return np.zeros((n, n)), ()
        if self._pcm is None:
            cavity = self.model.cavity
            first, second = self.surface
            amat, _ = cavity.assemble_amat()
            factor = cho_factor(amat)
            # Unit dielectric response: independent components share geometry
            # and host potential, but each supplies its own dielectric factor.
            q = -cho_solve(factor, self.parameters.electrostatic_potential(self.coords))
            phi1 = np.empty((len(self.coords), n))
            direct = np.zeros((n, n))
            for i, point in enumerate(self.coords):
                phi1[i], phi2 = self.parameters.potential(point, first[:, i], second[:, i])
                direct += q[i] * phi2
            aq1, a2 = cavity.amat_host_derivatives(q, first, second)
            rhs = aq1 + phi1
            self._pcm = direct + 0.5 * a2, rhs, cho_solve(factor, rhs)
        direct, rhs, solved = self._pcm
        return dielectric_factor * direct, ((dielectric_factor * rhs, solved),)
