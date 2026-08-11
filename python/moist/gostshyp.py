"""GOSTSHYP hydrostatic pressure on a moist isodensity cavity.

GOSTSHYP simulates hydrostatic pressure by placing an unnormalized Gaussian
potential on every cavity grid point and fixing its amplitude so the force the
wall exerts on the electron density matches ``p_inp`` times the grid point area:

.. math::

    \\omega_j &= \\pi \\ln 2 / a_j                                        \\\\
    g_{\\mu\\nu,j} &= \\langle \\mu | e^{-\\omega_j |r-C_j|^2} | \\nu \\rangle  \\\\
    f_{\\mu\\nu,j} &= n_j \\cdot \\nabla_r g_{\\mu\\nu,j}                    \\\\
    p_j &= p_\\mathrm{inp} a_j / \\tilde f_j                              \\\\
    E^\\mathrm{GOST} &= \\sum_j p_j \\tilde g_j

with ``gtilde``/``ftilde`` the density contractions of ``g``/``f``.  The Fock
matrix at a *frozen* surface is

.. math::

    V_{\\mu\\nu} = \\sum_j \\big[ \\alpha_j g_{\\mu\\nu,j} - \\beta_j f_{\\mu\\nu,j} \\big],
    \\quad \\alpha_j = p_j, \\; \\beta_j = \\tilde g_j p_j / \\tilde f_j

Almost all of that is QM-side integral work.  What moist owns is the other
half: the derivatives of :math:`E^\\mathrm{GOST}` with respect to the *cavity
parameters* — the grid point area, position and outward normal — and the chain
that turns them into a level-set response and a nuclear gradient.  This module
computes the surface weights and hands them to moist; moist returns the
level-set adjoints, which :mod:`moist.pyscf` contracts with ``dS/dP`` and
``dS/dR``.

Every per-grid point Gaussian carries a normalization constant ``N_j`` that
cancels exactly between energy and Fock, so it is never formed.  Only the
*relative* s/p/d/f angular constants are restored, which is what makes
``f = n . grad g`` hold exactly rather than up to a factor.

Conventions, each pinned by a named test in ``test_gostshyp.py``:

* ``w_f`` is identically zero.  The DROP switching function is an anchor-only
  iSwiG overlap, so ``df/dS = 0`` for a callback level set (see the "No w_f
  term" comment in ``src/moist/cavity/drop/derivatives/seeds.f90``); the whole
  area route runs through the Gaussian width.  This holds for *callback* level
  sets only, which is why the model is built on an isodensity cavity.
* The C contraction API exposes no ``w_a`` channel, so the host folds the area
  route itself as ``w_xi = dE_da * (-2 a / xi)`` — the same identity moist
  applies internally in ``src/moist/cavity/drop/derivatives/weights.f90``.
* ``w_n`` is passed to moist **raw**.  ``contract_surface_lsf_weights`` already
  performs both the direct ``P_tan(w_n)/|grad S|`` term and the ``H @
  normal_grad`` point-motion fold (``drop/derivatives/seeds.f90``); pre-folding
  it here would count the normal response twice.  The *anchor* channel is the
  one place the host must fold it, because ``get_anchor_gradient`` exposes only
  the projected-point sensitivity.
* Grid order is moist native throughout.  There is no owner-sorted
  permutation anywhere in moist's Python layer.

All quantities are atomic units: ``pressure`` in Hartree/bohr^3, positions in
bohr, areas in bohr^2, ``omega`` in bohr^-2, the energy in Hartree.  Converting
from GPa is the caller's job; :data:`GPA_TO_AU` is provided for that.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from pyscf import gto

from .interface import IsodensityDROPCavity
from .pyscf import PySCFIsodensityHost

__all__ = ["GostshypWall", "GPA_TO_AU"]

#: 1 GPa in Hartree / bohr^3.
GPA_TO_AU = 1.0e9 * 5.29177210903e-11**3 / 4.3597447222071e-18

#: libcint angular constants for a coefficient-1 fakemol shell.  Restoring the
#: ratios ``N_s/N_l`` puts every moment in the same units as ``g``, which is
#: what makes ``f == n . grad g`` exact.  s: 1/(2 sqrt(pi)), p: sqrt(3/(4 pi)),
#: d: 1, f: 1.
_S_NORM = 1.0 / (2.0 * math.sqrt(math.pi))
_S_OVER_P_NORM = 1.0 / math.sqrt(3.0)
_S_OVER_D_NORM = _S_NORM
_S_OVER_F_NORM = _S_NORM

#: libcint cartesian d order: xx xy xz yy yz zz.
_D_CART_ORDER = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
#: Cartesian f components summing to ``(r_a - C_a) |r - C|^2`` for a = x, y, z.
#: The f order is xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz, so x picks
#: xxx/xyy/xzz, y picks xxy/yyy/yzz and z picks xxz/yyz/zzz.
_F_RHO2_FIRST_MOMENT = ((0, 3, 5), (1, 6, 8), (2, 7, 9))

#: Relative floor on ``|ftilde_j|`` below which a grid point is inactive.
#:
#: ====== ============== =========== ===============
#: floor   E error (rel)  Fock (rel)  gradient (rel)
#: ====== ============== =========== ===============
#: 1e-6   3.9e-6         2.1e-11     1.4e-11
#: 1e-8   5.0e-7         2.7e-11     9.2e-12
#: 1e-9   1.7e-7         3.7e-13     1.9e-11
#: 1e-10  3.5e-9         2.3e-11     3.9e-10
#: 1e-12  3.4e-10        2.1e-11     9.4e-09
#: 1e-14  0 (reference)  6.0e-10     4.3e-08
#: ====== ============== =========== ===============
#:
_OVERLAP_FLOOR = 1.0e-9


def _fakemol_gaussians(coords: np.ndarray, exponents: np.ndarray, angl: int) -> gto.Mole:
    """One coefficient-1 GTO shell of angular momentum ``angl`` per grid point."""

    coords = np.asarray(coords, dtype=np.float64)
    exponents = np.asarray(exponents, dtype=np.float64)
    nshell = coords.shape[0]

    fakemol = gto.Mole()
    fakemol._atm = np.zeros((nshell, gto.ATM_SLOTS), dtype=np.int32)
    fakemol._bas = np.zeros((nshell, gto.BAS_SLOTS), dtype=np.int32)
    env = [0.0] * gto.PTR_ENV_START

    for ishell in range(nshell):
        fakemol._atm[ishell, gto.PTR_COORD] = len(env)
        env.extend(coords[ishell])
        fakemol._bas[ishell, gto.ATOM_OF] = ishell
        fakemol._bas[ishell, gto.ANG_OF] = angl
        fakemol._bas[ishell, gto.NPRIM_OF] = 1
        fakemol._bas[ishell, gto.NCTR_OF] = 1
        fakemol._bas[ishell, gto.PTR_EXP] = len(env)
        fakemol._bas[ishell, gto.PTR_COEFF] = len(env) + 1
        env.extend((float(exponents[ishell]), 1.0))

    fakemol._env = np.asarray(env, dtype=np.float64)
    fakemol._built = True
    return fakemol


def _int3c1e(mol, centers, omega, angl, intor="int3c1e_cart"):
    """Three-centre one-electron integrals over a Gaussian-per-grid point fakemol."""

    fakemol = _fakemol_gaussians(centers, omega, angl)
    nbas = mol.nbas
    shls_slice = (0, nbas, 0, nbas, nbas, nbas + fakemol.nbas)
    return (mol + fakemol).intor(intor, shls_slice=shls_slice)


@dataclass(frozen=True)
class _GostshypState:
    """Per-density GOSTSHYP intermediates on one frozen cavity."""

    #: Total pressure-wall energy, Hartree.
    energy: float
    #: Fixed-surface Fock matrix (eq 16), ``(nao, nao)``.
    fock: np.ndarray
    #: ``alpha_j = p_j = p_inp a_j / ftilde_j``.
    alpha: np.ndarray
    #: ``beta_j = gtilde_j p_j / ftilde_j``.
    beta: np.ndarray
    #: Density-contracted Gaussian overlap.
    gtilde: np.ndarray
    #: Density-contracted normal-projected Gaussian gradient.
    ftilde: np.ndarray
    #: Grid point carrying a usable density overlap.
    active: np.ndarray


@dataclass(frozen=True)
class _SurfaceAdjoints:
    """Energy sensitivities to the DROP surface, per grid point, native order."""

    #: ``(ngrid,)`` area route mapped onto the Gaussian width.
    w_xi: np.ndarray
    #: ``(ngrid,)`` identically zero; see the module docstring.
    w_f: np.ndarray
    #: ``(ngrid, 3)`` ``dE/dr_j`` at a *fixed* normal.
    w_xyz: np.ndarray
    #: ``(ngrid, 3)`` ``dE/dn_j``.
    w_n: np.ndarray
    #: ``(ngrid,)`` raw ``dE/da_j`` for the anchor area route.
    dE_da: np.ndarray


@dataclass(frozen=True)
class _LsfWeights:
    """Level-set adjoints, shaped like :class:`~moist.interface.GeneralPotential`.

    :meth:`~moist.pyscf.PySCFIsodensityHost._fock_lsf` and ``_gradient_lsf``
    read only these three attributes, so this stands in for a potential object
    without a solvation model being involved at all.
    """

    #: ``(ngrid,)``
    w_lsf0: np.ndarray
    #: Fortran ``(3, ngrid)``
    w_lsf1: np.ndarray
    #: Fortran ``(3, 3, ngrid)``
    w_lsf2: np.ndarray


class GostshypWall:
    """GOSTSHYP pressure wall on a moist isodensity DROP cavity.

    :param host: PySCF host supplying the level set, geometry and density.
    :param pressure: ``p_inp`` in Hartree/bohr^3 (multiply GPa by
        :data:`GPA_TO_AU`).
    :param cavity_kwargs: forwarded to
        :meth:`~moist.pyscf.PySCFIsodensityHost.make_cavity` (``nleb``,
        ``tolerance``, ...).

    The cavity follows the density, so :meth:`update` must run before any
    energy, Fock or gradient is read.
    """

    def __init__(
        self,
        host: PySCFIsodensityHost,
        pressure: float,
        **cavity_kwargs,
    ) -> None:
        self.host = host
        self.mol = host.mol
        self.pressure = float(pressure)
        self.cavity: IsodensityDROPCavity = host.make_cavity(**cavity_kwargs)

        self.energy = 0.0
        self.ngrid = 0
        self._state: Optional[_GostshypState] = None
        self._c2s: Optional[np.ndarray] = None
        self._G: Optional[np.ndarray] = None
        self._F: Optional[np.ndarray] = None
        self._Fvec: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # cavity and grid point geometry
    # ------------------------------------------------------------------

    def _set_grid_points(self) -> None:
        """Snapshot the live cavity into the four arrays GOSTSHYP needs."""

        result = self.cavity.cavity
        self.centers = np.ascontiguousarray(result.xyz.T, dtype=np.float64)
        self.areas = np.ascontiguousarray(result.a, dtype=np.float64)
        self.ngrid = int(self.centers.shape[0])
        self.nsph = int(result.nsph)

        normals = np.ascontiguousarray(result.normal0.T, dtype=np.float64)
        norm = np.linalg.norm(normals, axis=1)
        good = norm > 0.0
        normals[good] /= norm[good, None]
        self.normals = normals

        # omega_j = pi ln2 / a_j; a degenerate zero-area grid point is inert.
        with np.errstate(divide="ignore", invalid="ignore"):
            omega = np.pi * math.log(2.0) / self.areas
        omega[~np.isfinite(omega)] = 0.0
        self.omega = omega

        self.xi, _switch = self.cavity.get_gaussian()

    def update(self, dm: np.ndarray) -> float:
        """Rebuild the cavity and integrals at ``dm`` and return the energy.

        The cached results are dropped *before* the rebuild.  The host density
        and the cavity are already mutated by the time anything downstream can
        fail, so keeping the previous state would leave every accessor
        reporting coherent-looking numbers for a surface that no longer exists.
        """

        self._state = None
        self.energy = 0.0

        self.host.dm = dm
        self.cavity.update(self.host.structure())
        self._set_grid_points()
        self._build_integrals()
        self._state = self._compute(dm)
        self.energy = self._state.energy
        return self.energy

    # ------------------------------------------------------------------
    # integrals (dense)
    # ------------------------------------------------------------------

    @property
    def _cart2sph(self) -> np.ndarray:
        """Cartesian-to-spherical AO transform, ``(nao_cart, nao)``."""

        if self._c2s is None:
            if self.mol.cart:
                self._c2s = np.eye(self.mol.nao_nr())
            else:
                self._c2s = np.asarray(self.mol.cart2sph_coeff(normalized="sp"))
        return self._c2s

    def _to_spherical(self, block: np.ndarray) -> np.ndarray:
        """Transform the two AO legs of a ``(ncart, ncart, ...)`` block."""

        c2s = self._cart2sph
        out = np.tensordot(c2s, block, axes=(0, 0))
        return np.tensordot(c2s, out, axes=(0, 1)).swapaxes(0, 1)

    def _build_integrals(self) -> None:
        """Dense ``g``, ``f`` and the unprojected f-vector on the live surface."""

        g_cart = _int3c1e(self.mol, self.centers, self.omega, 0)
        p_cart = _int3c1e(self.mol, self.centers, self.omega, 1)
        ncart = g_cart.shape[0]
        p_cart = p_cart.reshape(ncart, ncart, self.ngrid, 3)

        # dG/dC_a = 2 omega (r_a - C_a) G, and displacing the field point is the
        # opposite of displacing the centre, so grad_r g = -2 omega <(r-C) G>.
        fvec_cart = -_S_OVER_P_NORM * 2.0 * np.einsum(
            "j,pqja->pqja", self.omega, p_cart, optimize=True
        )
        f_cart = np.einsum("pqja,ja->pqj", fvec_cart, self.normals, optimize=True)

        self._G = self._to_spherical(g_cart)
        self._F = self._to_spherical(f_cart)
        self._Fvec = self._to_spherical(fvec_cart)

    def _density_matrix_cart(self, dm: np.ndarray) -> np.ndarray:
        """Density matrix in the cartesian AO basis the fakemol blocks use.

        BLAS leaves dirty floating-point status flags behind for these shapes --
        the SIMD tail reads padding lanes -- so numpy reports divide-by-zero and
        overflow from a product that performs neither.  Suppressed here for the
        same reason :meth:`~moist.pyscf.PySCFIsodensityHost.lsf` suppresses them.
        """

        c2s = self._cart2sph
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            return c2s @ np.asarray(dm) @ c2s.T

    def _gf_tilde(self, dm_cart, centers, omega, normals):
        """``(gtilde, ftilde)`` for *arbitrary* grid point parameters.

        The finite-difference reference for :meth:`_param_derivatives`; it must
        not read any cached surface state.
        """

        g_cart = _int3c1e(self.mol, centers, omega, 0)
        p_cart = _int3c1e(self.mol, centers, omega, 1)
        ncart = g_cart.shape[0]
        ngrid = int(np.asarray(centers).shape[0])
        p_cart = p_cart.reshape(ncart, ncart, ngrid, 3)

        gtilde = np.einsum("pqj,pq->j", g_cart, dm_cart, optimize=True)
        pt = _S_OVER_P_NORM * np.einsum("pqja,pq->ja", p_cart, dm_cart, optimize=True)
        ftilde = -2.0 * np.asarray(omega) * np.einsum(
            "ja,ja->j", np.asarray(normals), pt, optimize=True
        )
        return gtilde, ftilde

    # ------------------------------------------------------------------
    # energy and fixed-surface Fock
    # ------------------------------------------------------------------

    def _compute(self, dm: np.ndarray) -> _GostshypState:
        """Amplitudes, energy and the eq-16 Fock at the current surface."""

        dm = np.asarray(dm)
        gtilde = np.einsum("uvj,uv->j", self._G, dm, optimize=True)
        ftilde = np.einsum("uvj,uv->j", self._F, dm, optimize=True)

        # p_j = p_inp a_j / ftilde_j.  A grid point that has left the density
        # carries no reliable ratio gtilde_j/ftilde_j at all, so it is dropped
        # here once and stays dropped in every derivative; see _OVERLAP_FLOOR.
        floor = _OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha = self.pressure * self.areas / ftilde
        active = np.isfinite(alpha) & (np.abs(ftilde) > floor)
        alpha = np.where(active, alpha, 0.0)

        energy = float(np.dot(alpha, gtilde))

        with np.errstate(divide="ignore", invalid="ignore"):
            beta = gtilde * alpha / ftilde
        beta = np.where(active, beta, 0.0)

        fock = np.einsum("j,uvj->uv", alpha, self._G, optimize=True)
        fock -= np.einsum("j,uvj->uv", beta, self._F, optimize=True)
        fock = 0.5 * (fock + fock.T)

        return _GostshypState(
            energy=energy,
            fock=fock,
            alpha=alpha,
            beta=beta,
            gtilde=gtilde,
            ftilde=ftilde,
            active=active,
        )

    def _require_state(self) -> _GostshypState:
        if self._state is None:
            raise RuntimeError("call update(dm) successfully before reading GOSTSHYP results")
        return self._state

    @property
    def amplitudes(self) -> Optional[np.ndarray]:
        """The per-grid point amplitudes ``p_j`` of the last :meth:`update`."""

        return None if self._state is None else self._state.alpha

    def effective_volume(self) -> float:
        """``E / p_inp`` (eq 11).  Not the cavity volume.

        Evaluated as ``sum_j a_j gtilde_j / ftilde_j`` over the active  grid points
        rather than by dividing the energy: that ratio is what the quantity
        *is*, and it stays well defined at ``p_inp = 0``, where the energy
        vanishes with the pressure but the volume it reports does not.
        """

        state = self._require_state()
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = self.areas * state.gtilde / state.ftilde
        return float(np.sum(np.where(state.active, ratio, 0.0)))

    def fock(self, dm: np.ndarray, *, include_cavity_response: bool = True) -> np.ndarray:
        """Solvation contribution to the Fock matrix, ``(nao, nao)``.

        :param include_cavity_response: contract the level-set adjoints.  With
            ``False`` this is the eq-16 matrix at a frozen surface, which is
            *not* the derivative of the energy an isodensity cavity reports —
            the surface moves with the density.
        """

        state = self._require_state()
        fock = state.fock
        if include_cavity_response:
            fock = fock + self.host._fock_lsf(self.centers, self.lsf_weights(dm))
        return fock

    # ------------------------------------------------------------------
    # surface moments and parameter derivatives
    # ------------------------------------------------------------------

    @staticmethod
    def _contract_p_moments(dm_cart, p_cart) -> np.ndarray:
        return _S_OVER_P_NORM * np.einsum("pqja,pq->ja", p_cart, dm_cart, optimize=True)

    @staticmethod
    def _contract_d_moments(dm_cart, d_cart) -> np.ndarray:
        raw = _S_OVER_D_NORM * np.einsum("pqjc,pq->jc", d_cart, dm_cart, optimize=True)
        ngrid = raw.shape[0]
        moment = np.empty((ngrid, 3, 3), dtype=np.float64)
        for component, (a, b) in enumerate(_D_CART_ORDER):
            moment[:, a, b] = raw[:, component]
            moment[:, b, a] = raw[:, component]
        return moment

    @staticmethod
    def _contract_f_rho2_moments(dm_cart, f_cart) -> np.ndarray:
        raw = _S_OVER_F_NORM * np.einsum("pqjc,pq->jc", f_cart, dm_cart, optimize=True)
        return np.stack([raw[:, list(idx)].sum(axis=1) for idx in _F_RHO2_FIRST_MOMENT], axis=1)

    def _surface_moments(self, dm_cart: np.ndarray):
        """``(gt, Pt, Mt, Rt)`` — the s/p/d/f Gaussian moments of the density."""

        g_cart = _int3c1e(self.mol, self.centers, self.omega, 0)
        p_cart = _int3c1e(self.mol, self.centers, self.omega, 1)
        d_cart = _int3c1e(self.mol, self.centers, self.omega, 2)
        f_cart = _int3c1e(self.mol, self.centers, self.omega, 3)
        ncart = g_cart.shape[0]

        gt = np.einsum("pqj,pq->j", g_cart, dm_cart, optimize=True)
        pt = self._contract_p_moments(dm_cart, p_cart.reshape(ncart, ncart, self.ngrid, 3))
        mt = self._contract_d_moments(dm_cart, d_cart.reshape(ncart, ncart, self.ngrid, 6))
        rt = self._contract_f_rho2_moments(dm_cart, f_cart.reshape(ncart, ncart, self.ngrid, 10))
        return gt, pt, mt, rt

    def _param_derivatives_from_moments(self, moments):
        """Surface-parameter derivatives of ``gtilde``/``ftilde``.

        With ``G = exp(-omega |r-C|^2)`` and the normal held fixed::

            dgtilde/dC_a = 2 omega Pt_a
            dftilde/dC_b = 2 omega n_b gt - 4 omega^2 (n . Mt)_b
            dgtilde/domega = -tr(Mt)
            dftilde/domega = -2 (n . Pt) + 2 omega (n . Rt)
        """

        gt, pt, mt, rt = moments
        omega = self.omega
        normals = self.normals
        two_omega = 2.0 * omega

        n_pt = np.einsum("ja,ja->j", normals, pt, optimize=True)
        n_mt = np.einsum("ja,jab->jb", normals, mt, optimize=True)
        n_rt = np.einsum("ja,ja->j", normals, rt, optimize=True)

        dgdr = two_omega[:, None] * pt
        dfdr = two_omega[:, None] * normals * gt[:, None] - 4.0 * omega[:, None] ** 2 * n_mt
        dgdw = -np.einsum("jaa->j", mt, optimize=True)
        dfdw = -2.0 * n_pt + two_omega * n_rt
        return dgdr, dfdr, dgdw, dfdw

    def _param_derivatives(self, dm_cart: np.ndarray):
        return self._param_derivatives_from_moments(self._surface_moments(dm_cart))

    def _param_derivatives_fd(self, dm_cart, h_r: float = 1.0e-4, h_w: float = 1.0e-4):
        """Central-difference reference for :meth:`_param_derivatives`.

        The ``omega`` step is *relative*: the exponents span orders of magnitude
        across the  grid points, so one absolute step cannot serve them all.
        """

        centers = self.centers
        omega = self.omega
        normals = self.normals

        dgdr = np.zeros((self.ngrid, 3))
        dfdr = np.zeros((self.ngrid, 3))
        for axis in range(3):
            shifted = centers.copy()
            shifted[:, axis] += h_r
            gp, fp = self._gf_tilde(dm_cart, shifted, omega, normals)
            shifted[:, axis] -= 2.0 * h_r
            gm, fm = self._gf_tilde(dm_cart, shifted, omega, normals)
            dgdr[:, axis] = (gp - gm) / (2.0 * h_r)
            dfdr[:, axis] = (fp - fm) / (2.0 * h_r)

        delta = h_w * omega
        gp, fp = self._gf_tilde(dm_cart, centers, omega + delta, normals)
        gm, fm = self._gf_tilde(dm_cart, centers, omega - delta, normals)
        with np.errstate(divide="ignore", invalid="ignore"):
            dgdw = (gp - gm) / (2.0 * delta)
            dfdw = (fp - fm) / (2.0 * delta)
        dgdw[~np.isfinite(dgdw)] = 0.0
        dfdw[~np.isfinite(dfdw)] = 0.0
        return dgdr, dfdr, dgdw, dfdw

    # ------------------------------------------------------------------
    # surface adjoints -- the quantities moist actually consumes
    # ------------------------------------------------------------------

    def surface_adjoints(self, dm: np.ndarray) -> _SurfaceAdjoints:
        """``dE/d(a_j, r_j, n_j)`` per grid point, in moist native grid order.

        This is what a Fortran ``gostshyp`` component would hand to
        ``cavity_surface_adjoint_type``'s ``w_a``/``w_xyz``/``w_n`` channels.
        The C API carries no ``w_a``, so the area route is folded onto the
        Gaussian width here instead.
        """

        state = self._require_state()
        dm_cart = self._density_matrix_cart(dm)
        dgdr, dfdr, dgdw, dfdw = self._param_derivatives(dm_cart)

        alpha = state.alpha
        beta = state.beta
        # dftilde/dn_j is the unprojected f-vector, which is already stored.
        gvfield = np.einsum("uvja,uv->ja", self._Fvec, np.asarray(dm), optimize=True)

        w_xyz = alpha[:, None] * dgdr - beta[:, None] * dfdr
        # Only ftilde depends on the normal, through ftilde_j = n_j . gvfield_j.
        w_n = -beta[:, None] * gvfield

        # a_j enters twice: explicitly through the amplitude p_j = p a_j/ftilde_j,
        # and through the Gaussian exponent omega_j = pi ln2 / a_j.
        with np.errstate(divide="ignore", invalid="ignore"):
            domega_da = -self.omega / self.areas
            omega_route = (alpha * dgdw - beta * dfdw) * domega_da
            dE_da = self.pressure * state.gtilde / state.ftilde + omega_route
        dE_da = np.where(state.active, dE_da, 0.0)
        dE_da[~np.isfinite(dE_da)] = 0.0

        # a_i ~ xi_i^-2, so the area route reaches the level set through the
        # Gaussian width.  Mirrors moist's own fold in the DROP surface-weight
        # contraction (src/moist/cavity/drop/derivatives/weights.f90), which is
        # the code this stands in for until the C API carries a w_a channel.
        with np.errstate(divide="ignore", invalid="ignore"):
            w_xi = dE_da * (-2.0 * self.areas / self.xi)
        w_xi[~np.isfinite(w_xi)] = 0.0

        return _SurfaceAdjoints(
            w_xi=w_xi,
            w_f=np.zeros_like(w_xi),
            w_xyz=w_xyz,
            w_n=w_n,
            dE_da=dE_da,
        )

    def lsf_weights(
        self, dm: np.ndarray, *, adjoints: Optional[_SurfaceAdjoints] = None
    ) -> _LsfWeights:
        """Contract the surface adjoints into level-set adjoints via moist.

        :param adjoints: an already-computed :meth:`surface_adjoints` result to
            reuse.  Building them costs four dense three-centre integral blocks,
            so a caller that needs them anyway should pass them in.
        """

        if adjoints is None:
            adjoints = self.surface_adjoints(dm)
        w_lsf0, w_lsf1, w_lsf2 = self.cavity.contract_surface_lsf_weights(
            adjoints.w_xi,
            adjoints.w_f,
            np.asfortranarray(adjoints.w_xyz.T),
            # Raw, not pre-folded: moist performs both the direct
            # P_tan(w_n)/|grad S| term and the H @ normal_grad point-motion fold.
            w_n=np.asfortranarray(adjoints.w_n.T),
        )
        return _LsfWeights(w_lsf0, w_lsf1, w_lsf2)

    # ------------------------------------------------------------------
    # nuclear gradient
    # ------------------------------------------------------------------

    def _integral_nuclear_gradient(self, dm: np.ndarray) -> np.ndarray:
        """AO centres move, surface frozen.  Fortran ``(3, natm)``."""

        state = self._require_state()
        dm_cart = self._density_matrix_cart(dm)
        ncart = self._cart2sph.shape[0]

        ip1_g = _int3c1e(self.mol, self.centers, self.omega, 0, "int3c1e_ip1_cart")
        ip1_p = _int3c1e(self.mol, self.centers, self.omega, 1, "int3c1e_ip1_cart")
        ip1_g = ip1_g.reshape(3, ncart, ncart, self.ngrid)
        ip1_p = ip1_p.reshape(3, ncart, ncart, self.ngrid, 3)
        ip1_f = -_S_OVER_P_NORM * 2.0 * np.einsum(
            "j,xpqja,ja->xpqj", self.omega, ip1_p, self.normals, optimize=True
        )

        kernel = np.einsum("j,xpqj->xpq", state.alpha, ip1_g, optimize=True)
        kernel -= np.einsum("j,xpqj->xpq", state.beta, ip1_f, optimize=True)
        row = np.einsum("xpq,pq->xp", kernel, dm_cart, optimize=True)

        # Cartesian AO rows fold onto atoms through the cartesian slices.
        gradient = np.zeros((3, self.mol.natm))
        aoslice = self.mol.aoslice_by_atom(self.mol.ao_loc_nr(cart=True))
        for atom in range(self.mol.natm):
            lo, hi = aoslice[atom, 2:4]
            # dg/dR = -ip1, and both AO legs contribute equally for symmetric P.
            gradient[:, atom] = -2.0 * row[:, lo:hi].sum(axis=1)
        return np.asfortranarray(gradient)

    def _field_nuclear_gradient(
        self, dm: np.ndarray, *, adjoints: Optional[_SurfaceAdjoints] = None
    ) -> np.ndarray:
        """The level set's own dependence on the nuclei, through the AO basis."""

        weights = self.lsf_weights(dm, adjoints=adjoints)
        return self.host._gradient_lsf(self.centers, weights)

    def _lsf_jet_at_grid_points(self):
        """``(|grad S|, H)`` at the  grid points from the host's *unscaled* level set.

        moist scales value, gradient and Hessian by one common factor, so the
        combination the anchor fold needs -- ``H / |grad S|`` -- is
        scale-invariant, and the unscaled jet is exact for it.  The residual
        ``normal_grad`` on its own is *not* scale-invariant, which is why it
        never leaves this routine.
        """

        grad_norm = np.empty(self.ngrid)
        hess = np.empty((self.ngrid, 3, 3))
        for igrid in range(self.ngrid):
            _value, gradient, hessian = self.host.lsf(self.centers[igrid], 2)
            grad_norm[igrid] = np.linalg.norm(gradient)
            hess[igrid] = hessian
        return grad_norm, hess

    def _anchor_nuclear_gradient(
        self, dm: np.ndarray, *, adjoints: Optional[_SurfaceAdjoints] = None
    ) -> np.ndarray:
        """Rigid owner-atom motion of the  grid points at a frozen level-set field.

        The area route uses the *true* per-point ``da_i/dR_A``: the grid point area
        carries a switching-function dependence (``a_i ~ f_i/xi_i^2``) that the
        Gaussian-width proxy misses for nuclear motion.
        """

        if adjoints is None:
            adjoints = self.surface_adjoints(dm)
        grad_norm, hess = self._lsf_jet_at_grid_points()

        nwn = np.einsum("ia,ia->i", self.normals, adjoints.w_n, optimize=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            normal_grad = (adjoints.w_n - self.normals * nwn[:, None]) / grad_norm[:, None]
        normal_grad[~np.isfinite(normal_grad)] = 0.0
        w_xyz_total = adjoints.w_xyz + np.einsum("iab,ib->ia", hess, normal_grad, optimize=True)

        self.cavity.compute_anchor_gradient()
        anchor = self.cavity.get_anchor_gradient()

        gradient = np.einsum("ic,caAi->aA", w_xyz_total, anchor.xyz1_rA, optimize=True)
        gradient += np.einsum("i,aAi->aA", adjoints.dE_da, anchor.a_i1_rA, optimize=True)
        return np.asfortranarray(gradient)

    def nuclear_gradient(self, dm: np.ndarray) -> np.ndarray:
        """Total ``dE^GOST/dR`` at fixed ``dm``, Fortran ``(3, natm)``.

        Three routes a nuclear displacement takes: the AO centres move at a
        frozen surface (integral), the density and hence the level set moves
        (field), and the atom-anchored reference grid is dragged rigidly
        (anchor).  At a converged variational SCF this is also the total
        pressure-wall gradient, with no orbital response.

        The field and anchor routes read the same surface adjoints, which cost
        four dense three-centre integral blocks, so they are built once here.
        """

        adjoints = self.surface_adjoints(dm)
        return (
            self._integral_nuclear_gradient(dm)
            + self._field_nuclear_gradient(dm, adjoints=adjoints)
            + self._anchor_nuclear_gradient(dm, adjoints=adjoints)
        )
