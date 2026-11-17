"""Two-stage second-order PySCF host protocol and analytic SCF Hessian."""

from __future__ import annotations

from contextlib import contextmanager
from itertools import permutations, product
import math
from typing import Optional

import numpy as np

from .configuration import ISwiG
from .interface import (
    DensityTangentRequest,
    DensityTangentResponse,
    GaussianAmplitudeResponse,
    GaussianMomentRequest,
    GaussianPotentialRequest,
    GaussianPotentialTangentRequest,
    NuclearChargesRequest,
    PointPotentialRequest,
    PointPotentialTangentRequest,
    PotentialAdjointResponse,
    PotentialAdjointTangentResponse,
)

#: Packed symmetric (3, 3) component order matching moist: [xx, xy, xz, yy, yz, zz].
_SYM2_IDX = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))

#: Packed symmetric (3, 3, 3) component order matching moist.
_SYM3_IDX = (
    (0, 0, 0),
    (0, 0, 1),
    (0, 0, 2),
    (0, 1, 1),
    (0, 1, 2),
    (0, 2, 2),
    (1, 1, 1),
    (1, 1, 2),
    (1, 2, 2),
    (2, 2, 2),
)

_D_TRACE_IDX = (0, 3, 5)
_S_NORM = 1.0 / (2.0 * math.sqrt(math.pi))
_S_OVER_P_NORM = 1.0 / math.sqrt(3.0)
_SQRT_PI = math.sqrt(math.pi)


@contextmanager
def _quiet_blas():
    """Suppress spurious FP status flags from SIMD padding lanes in BLAS."""
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        yield


def _finite(array: np.ndarray, what: str) -> np.ndarray:
    """Verify every element of ``array`` is finite after BLAS flag suppression."""
    if not np.all(np.isfinite(array)):
        raise FloatingPointError(f"Non-finite values in {what}")
    return array


def _contract(subscripts: str, *operands) -> np.ndarray:
    """``np.einsum(..., optimize=True)`` with BLAS flag suppression and finiteness check."""
    with _quiet_blas():
        return _finite(np.einsum(subscripts, *operands, optimize=True), "einsum contraction")


def _sorted_multi_indices(order: int):
    """Sorted Cartesian multi-indices of ``order`` paired with permutation multiplicity."""
    seen: dict[tuple[int, ...], int] = {}
    for axes in product(range(3), repeat=order):
        key = tuple(sorted(axes))
        seen[key] = seen.get(key, 0) + 1
    return list(seen.items())


def _symmetric_outer(tensor: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Symmetrized outer product of ``(3,)*k + (ngrid,)`` with ``(3, ngrid)``."""
    order = tensor.ndim - 1
    outer = _contract("...i,ci->...ci", tensor, vector)
    if order == 0:
        return outer
    total = np.zeros_like(outer)
    axes = list(range(order + 1))
    for position in range(order + 1):
        perm = axes[:position] + [order] + axes[position:order]
        total += np.transpose(outer, perm + [order + 1])
    return total / (order + 1)


def _unpack_sym2(packed: np.ndarray) -> np.ndarray:
    """Unpack ``(ngrid, 6)`` symmetric matrices to ``(3, 3, ngrid)``."""
    ngrid = packed.shape[0]
    out = np.empty((3, 3, ngrid), dtype=np.float64)
    for idx, (a, b) in enumerate(_SYM2_IDX):
        col = packed[:, idx]
        out[a, b] = col
        out[b, a] = col
    return out


class _DensityTables:
    """Geometry-cached AO jet tables to order 4 for bare density HVP coupling."""

    def __init__(self, host, coords: np.ndarray) -> None:
        self.host = host
        self.mol = host.mol
        self.nao = self.mol.nao
        self.ao = host._ao(coords, 4)
        self.atom_of_ao = np.empty(self.nao, dtype=np.intp)
        for atom, (_, _, lo, hi) in enumerate(host._aoslice):
            self.atom_of_ao[lo:hi] = atom
        self._reduced_dm: Optional[tuple[np.ndarray, np.ndarray]] = None

    def _ao_at(self, axes: tuple[int, ...]) -> np.ndarray:
        from .pyscf import _component_index

        return self.ao[_component_index(axes)]

    @staticmethod
    def _project(left: np.ndarray, dm: np.ndarray) -> np.ndarray:
        with _quiet_blas():
            return _finite(left @ dm, "AO density projection")

    def _right_dm(self, dm: np.ndarray) -> np.ndarray:
        """Cache ``ao @ dm`` across directions at fixed converged density."""
        if self._reduced_dm is not None and np.array_equal(self._reduced_dm[0], dm):
            return self._reduced_dm[1]
        with _quiet_blas():
            proj = _finite(self.ao @ dm, "AO density projection")
        self._reduced_dm = (np.array(dm, dtype=np.float64, copy=True), proj)
        return proj

    def _per_atom(self, per_ao: np.ndarray) -> np.ndarray:
        """Sum per-AO contributions ``(3, nao)`` into ``(natm, 3)``."""
        out = np.zeros((self.mol.natm, 3), dtype=np.float64)
        for atom, (_, _, lo, hi) in enumerate(self.host._aoslice):
            out[atom] = per_ao[:, lo:hi].sum(axis=1)
        return out

    def _terms(self, weights: dict[int, Optional[np.ndarray]]):
        for order, tensor in weights.items():
            if tensor is None:
                continue
            for key, multiplicity in _sorted_multi_indices(order):
                w = multiplicity * tensor[key + (slice(None),)] if order else multiplicity * tensor
                if not np.any(w):
                    continue
                for mask in product((False, True), repeat=order):
                    alpha = tuple(a for a, take in zip(key, mask) if take)
                    beta = tuple(a for a, take in zip(key, mask) if not take)
                    yield w, alpha, beta

    def density_tangent(
        self,
        dm: np.ndarray,
        dirs: np.ndarray,
        dm_dirs: Optional[np.ndarray],
        missing: frozenset[str],
    ) -> dict[str, np.ndarray]:
        """Compute bare ``d_rho``, ``d_grad_rho``, ``d_hess_rho``, ``d_d3_rho`` for ``dirs``."""
        ndir = dirs.shape[0]
        ngrid = self.ao.shape[1]
        max_order = 0
        if "d_grad_rho" in missing:
            max_order = 1
        if "d_hess_rho" in missing:
            max_order = 2
        if "d_d3_rho" in missing:
            max_order = 3

        from .pyscf import _component_index

        right_dm = self._right_dm(dm)
        # Precompute AO-pair component tangents for nuclear and orbital parts.
        tensors = {
            order: np.zeros((ndir, ngrid) + (3,) * order, dtype=np.float64)
            for order in range(max_order + 1)
        }
        for j in range(ndir):
            v_r = dirs[j]  # (natm, 3)
            dm_dir = None if dm_dirs is None else dm_dirs[j]
            has_r = bool(np.any(v_r))
            has_p = dm_dir is not None and bool(np.any(dm_dir))
            if not has_r and not has_p:
                continue
            t = v_r[self.atom_of_ao].T  # (3, nao)
            for order in range(max_order + 1):
                target = tensors[order][j]
                for key, _ in _sorted_multi_indices(order):
                    val = np.zeros(ngrid, dtype=np.float64)
                    for mask in product((False, True), repeat=order):
                        alpha = tuple(a for a, take in zip(key, mask) if take)
                        beta = tuple(a for a, take in zip(key, mask) if not take)
                        left = self._ao_at(alpha)
                        ibeta = _component_index(beta)
                        if has_p:
                            right_p = self._project(self.ao[ibeta], dm_dir)
                            val += _contract("iu,iu->i", left, right_p)
                        if has_r:
                            moved = sum(t[d] * self._ao_at(alpha + (d,)) for d in range(3))
                            val -= 2.0 * _contract("iu,iu->i", moved, right_dm[ibeta])
                    for axes in set(permutations(key)):
                        target[(slice(None),) + axes] = val

        outputs: dict[str, np.ndarray] = {}
        if "d_rho" in missing:
            outputs["d_rho"] = tensors[0]
        if "d_grad_rho" in missing:
            outputs["d_grad_rho"] = tensors[1]
        if "d_hess_rho" in missing:
            d_hess = np.empty((ndir, ngrid, 6), dtype=np.float64)
            for idx, (a, b) in enumerate(_SYM2_IDX):
                d_hess[:, :, idx] = tensors[2][:, :, a, b]
            outputs["d_hess_rho"] = d_hess
        if "d_d3_rho" in missing:
            d_d3 = np.empty((ndir, ngrid, 10), dtype=np.float64)
            for idx, (a, b, c) in enumerate(_SYM3_IDX):
                d_d3[:, :, idx] = tensors[3][:, :, a, b, c]
            outputs["d_d3_rho"] = d_d3
        return outputs

    def fock(self, weights: dict[int, Optional[np.ndarray]]) -> np.ndarray:
        """``sum_k weights[k] : d rho^(k) / d P``, ``(nao, nao)``."""
        out = np.zeros((self.nao, self.nao), dtype=np.float64)
        for w, alpha, beta in self._terms(weights):
            out += _contract("i,iu,iv->uv", w, self._ao_at(alpha), self._ao_at(beta))
        return out

    def gradient(self, weights: dict[int, Optional[np.ndarray]], dm: np.ndarray, *,
                 use_cached_dm: bool = False) -> np.ndarray:
        """``sum_k weights[k] : d rho^(k) / d R``, ``(natm, 3)``."""
        from .pyscf import _component_index

        per_ao = np.zeros((3, self.nao), dtype=np.float64)
        right_cached = self._right_dm(dm) if use_cached_dm else None
        for w, alpha, beta in self._terms(weights):
            right = (
                right_cached[_component_index(beta)]
                if right_cached is not None
                else self._project(self._ao_at(beta), dm)
            )
            for c in range(3):
                per_ao[c] += _contract("i,iu,iu->u", w, self._ao_at(alpha + (c,)), right)
        return self._per_atom(-2.0 * per_ao)

    def fock_centre_tangent(
        self, weights: dict[int, Optional[np.ndarray]], dirs: np.ndarray
    ) -> np.ndarray:
        """Directional derivative of ``fock`` along AO center motion ``dirs (natm, 3)``."""
        if not np.any(dirs):
            return np.zeros((self.nao, self.nao), dtype=np.float64)
        t = dirs[self.atom_of_ao].T
        cross = np.zeros((self.nao, self.nao), dtype=np.float64)
        for w, alpha, beta in self._terms(weights):
            left = sum(t[c] * self._ao_at(alpha + (c,)) for c in range(3))
            cross += _contract("i,iu,iv->uv", w, left, self._ao_at(beta))
        return -(cross + cross.T)

    def gradient_centre_tangent(
        self, weights: dict[int, Optional[np.ndarray]], dm: np.ndarray, dirs: np.ndarray
    ) -> np.ndarray:
        """Directional derivative of ``gradient`` along AO center motion ``dirs (natm, 3)``."""
        if not np.any(dirs):
            return np.zeros((self.mol.natm, 3), dtype=np.float64)
        from .pyscf import _component_index

        t = dirs[self.atom_of_ao].T
        right_cached = self._right_dm(dm)
        per_ao = np.zeros((3, self.nao), dtype=np.float64)
        for w, alpha, beta in self._terms(weights):
            right = right_cached[_component_index(beta)]
            moved_beta = sum(t[d] * self._ao_at(beta + (d,)) for d in range(3))
            right_moved = self._project(moved_beta, dm)
            for c in range(3):
                left_moved = sum(t[d] * self._ao_at(alpha + (c, d)) for d in range(3))
                per_ao[c] += _contract("i,iu,iu->u", w, left_moved, right)
                per_ao[c] += _contract("i,iu,iu->u", w, self._ao_at(alpha + (c,)), right_moved)
        return self._per_atom(2.0 * per_ao)

    def complete_density(
        self,
        dm: np.ndarray,
        dirs: np.ndarray,
        dm_dirs: Optional[np.ndarray],
        item: DensityTangentResponse,
        hvp_nuc: np.ndarray,
        fock_tan: np.ndarray,
    ) -> None:
        """Accumulate host AO density HVP and Fock tangent from ``DensityTangentResponse``."""
        w2_raw = np.asarray(item.w_hess_rho, dtype=np.float64).T
        w2_sym = 0.5 * (w2_raw + w2_raw.transpose(1, 0, 2))
        weights = {
            0: np.asarray(item.w_rho, dtype=np.float64),
            1: np.asarray(item.w_grad_rho, dtype=np.float64).T,
            2: w2_sym,
        }
        ndir = dirs.shape[0]
        for j in range(ndir):
            v_r = dirs[j]
            dm_dir = None if dm_dirs is None else dm_dirs[j]
            d_xyz = item.d_xyz[j].T  # (3, ngrid)
            dweights = {
                0: np.asarray(item.dw_rho[j], dtype=np.float64),
                1: np.asarray(item.dw_grad_rho[j], dtype=np.float64).T,
                2: _unpack_sym2(item.dw_hess_rho[j]),
            }
            moving = {
                1: _symmetric_outer(weights[0], d_xyz),
                2: _symmetric_outer(weights[1], d_xyz),
                3: _symmetric_outer(weights[2], d_xyz),
            }
            hvp_nuc[j] += self.gradient(dweights, dm, use_cached_dm=True)
            hvp_nuc[j] += self.gradient(moving, dm, use_cached_dm=True)
            if dm_dir is not None and np.any(dm_dir):
                hvp_nuc[j] += self.gradient(weights, dm_dir)
            hvp_nuc[j] += self.gradient_centre_tangent(weights, dm, v_r)

            fock_tan[j] += self.fock(dweights)
            fock_tan[j] += self.fock(moving)
            fock_tan[j] += self.fock_centre_tangent(weights, v_r)


class _ChargeWeighted:
    """Three-center integrals contracted with converged potential adjoints ``w_phi``."""

    def __init__(self, tables: _ElectrostaticTables, w_phi: np.ndarray) -> None:
        self.bra = _contract("i,cuvi->cuv", w_phi, tables.ip1)
        self.bra_bra = _contract("i,acuvi->acuv", w_phi, tables.ipip1)
        self.bra_ket = _contract("i,acuvi->acuv", w_phi, tables.ipvip1)


class _DensityContracted:
    """Three-center integrals contracted with converged AO density matrix ``dm``."""

    def __init__(self, tables: _ElectrostaticTables, dm: np.ndarray) -> None:
        self.field = _contract("cuvi,uv->ci", tables.ip2, dm)
        self.bra = _contract("cuvi,uv->cui", tables.ip1, dm)
        self.bra_point = _contract("acuvi,uv->acui", tables.ip1ip2, dm)
        # Spatial Hessian of -(uv|i) assembled by translational invariance from AO centers
        # instead of int3c2e_ipip2, which corrupts point-charge (1e16) diagonals.
        bra_bra = _contract("acuvi,uv->aci", tables.ipip1, dm)
        ket_ket = _contract("acuvi,vu->aci", tables.ipip1, dm)
        bra_ket = _contract("acuvi,uv->aci", tables.ipvip1, dm)
        self.point = bra_bra + ket_ket + bra_ket + bra_ket.transpose(1, 0, 2)
        if tables.xi is not None:
            self.dphi_el_dxi = -4.0 * _contract("uvi,uv->i", tables.g_xi, dm)
            self.d2phi_el_dr_dxi = -8.0 * _S_OVER_P_NORM * (tables.xi**2)[None, :] * _contract(
                "auvi,uv->ai", tables.p_xi, dm
            )
            self.d2phi_el_dxi2 = 8.0 * _S_NORM * tables.xi * _contract(
                "uvi,uv->i", tables.d_xi_trace, dm
            )
            # Contract bra AO index to atom later in d_dphi_dxi or completion.
            self.bra_dxi = 8.0 * _contract("kuvi,uv->kui", tables.ip1_xi, dm)
        else:
            self.dphi_el_dxi = None
            self.d2phi_el_dr_dxi = None
            self.d2phi_el_dxi2 = None
            self.bra_dxi = None


class _ElectrostaticTables:
    """Geometry-cached three-center integrals for point or Gaussian potential HVP coupling."""

    def __init__(self, host, coords: np.ndarray, xi: Optional[np.ndarray]) -> None:
        from pyscf import df, gto
        from .pyscf import _int3c1e

        self.host = host
        self.mol = host.mol
        self.nao = self.mol.nao
        self.coords = np.asarray(coords, dtype=np.float64)
        self.xi = None if xi is None else np.asarray(xi, dtype=np.float64)
        self.aoslice = host._aoslice
        self.atom_of_ao = np.empty(self.nao, dtype=np.intp)
        for atom, (_, _, lo, hi) in enumerate(self.aoslice):
            self.atom_of_ao[lo:hi] = atom

        if self.xi is None:
            aux = gto.fakemol_for_charges(self.coords)
        else:
            aux = gto.fakemol_for_charges(self.coords, expnt=self.xi**2)
        aux.cart = self.mol.cart

        def aux_ints(name: str) -> np.ndarray:
            return df.incore.aux_e2(
                self.mol, aux, intor=self.mol._add_suffix(name), aosym="s1"
            )

        self.i0 = aux_ints("int3c2e")
        self.ip1 = aux_ints("int3c2e_ip1")
        self.ip2 = aux_ints("int3c2e_ip2")
        self.ipip1 = aux_ints("int3c2e_ipip1").reshape(3, 3, self.nao, self.nao, -1)
        self.ipvip1 = aux_ints("int3c2e_ipvip1").reshape(3, 3, self.nao, self.nao, -1)
        self.ip1ip2 = aux_ints("int3c2e_ip1ip2").reshape(3, 3, self.nao, self.nao, -1)

        if self.xi is not None:
            c2s = None if self.mol.cart else np.asarray(self.mol.cart2sph_coeff(normalized="sp"))
            omega = self.xi**2
            ngrid = self.xi.size
            g_cart = _int3c1e(self.mol, self.coords, omega, 0, normalized=False)
            ncart = g_cart.shape[0]
            p_cart = (
                _int3c1e(self.mol, self.coords, omega, 1, normalized=False)
                .reshape(ncart, ncart, ngrid, 3)
                .transpose(3, 0, 1, 2)
            )
            d_cart = (
                _int3c1e(self.mol, self.coords, omega, 2, normalized=False)
                .reshape(ncart, ncart, ngrid, 6)[..., _D_TRACE_IDX]
                .sum(axis=-1)
            )
            self.g_xi = self._cart_to_ao_2d(c2s, g_cart)
            self.p_xi = np.stack([self._cart_to_ao_2d(c2s, p_cart[a]) for a in range(3)], axis=0)
            self.d_xi_trace = self._cart_to_ao_2d(c2s, d_cart)
            self.ip1_xi = _int3c1e(
                self.mol,
                self.coords,
                omega,
                0,
                self.mol._add_suffix("int3c1e_ip1"),
                normalized=False,
            ).reshape(3, self.nao, self.nao, ngrid)
        else:
            self.g_xi = None
            self.p_xi = None
            self.d_xi_trace = None
            self.ip1_xi = None

        self._reduced: dict[str, tuple[np.ndarray, object]] = {}

    @staticmethod
    def _cart_to_ao_2d(c2s: Optional[np.ndarray], block: np.ndarray) -> np.ndarray:
        if c2s is None:
            return block
        with _quiet_blas():
            out = np.tensordot(c2s, block, axes=(0, 0))
            return _finite(np.tensordot(c2s, out, axes=(0, 1)).swapaxes(0, 1), "c2s integral transform")

    def _per_atom(self, per_ao: np.ndarray) -> np.ndarray:
        out = np.zeros((self.mol.natm, 3), dtype=np.float64)
        for atom, (_, _, lo, hi) in enumerate(self.aoslice):
            out[atom] = per_ao[:, lo:hi].sum(axis=1)
        return out

    def _reduce(self, slot: str, key: np.ndarray, build):
        cached = self._reduced.get(slot)
        if cached is not None and np.array_equal(cached[0], key):
            return cached[1]
        value = build()
        self._reduced[slot] = (np.array(key, dtype=np.float64, copy=True), value)
        return value

    def charge_weighted(self, w_phi: np.ndarray) -> _ChargeWeighted:
        return self._reduce("w_phi", w_phi, lambda: _ChargeWeighted(self, w_phi))

    def density_contracted(self, dm: np.ndarray) -> _DensityContracted:
        return self._reduce("dm", dm, lambda: _DensityContracted(self, dm))

    def _nuclear_radial(self):
        """Return ``(delta, c0, c1, r2, exp_x2)`` matching ``pcm_elstat_potential_tangent``."""
        from scipy.special import erf

        delta = self.coords[:, None, :] - self.mol.atom_coords()[None, :, :]  # (ngrid, natm, 3)
        r2 = np.sum(delta * delta, axis=2)
        dist = np.sqrt(r2)
        if self.xi is None:
            inv_r3 = np.divide(1.0, r2 * dist, out=np.zeros_like(r2), where=r2 > 1.0e-30)
            c0 = inv_r3
            c1 = -3.0 * np.divide(inv_r3, r2, out=np.zeros_like(r2), where=r2 > 1.0e-30)
            return delta, c0, c1, r2, None

        width = self.xi[:, None]
        w2 = width * width
        x = width * dist
        x2 = x * x
        exp_x2 = np.exp(-x2)
        coeff = 2.0 * exp_x2 / _SQRT_PI
        with np.errstate(divide="ignore", invalid="ignore"):
            c0_far = (erf(x) - x * coeff) / (r2 * dist)
            c1_far = (2.0 * (w2 * width) * coeff - 3.0 * c0_far) / r2
        small = np.abs(x) < 0.05
        c0_near = (4.0 * w2 * width / (3.0 * _SQRT_PI)) * (
            1.0 - 0.6 * x2 + 3.0 * (x2 * x2) / 14.0
        )
        c1_near = (8.0 * w2 * w2 * width / (5.0 * _SQRT_PI)) * (
            -1.0 + (5.0 / 7.0) * x2 - (5.0 / 18.0) * (x2 * x2)
        )
        c0 = np.where(small, c0_near, c0_far)
        c1 = np.where(small, c1_near, c1_far)
        return delta, c0, c1, r2, exp_x2

    def potential_tangent(
        self,
        dm: np.ndarray,
        dirs: np.ndarray,
        dm_dirs: Optional[np.ndarray],
        d_xyz: np.ndarray,
        d_xi: Optional[np.ndarray],
        missing: frozenset[str],
    ) -> dict[str, np.ndarray]:
        """Answer ``PointPotentialTangentRequest`` or ``GaussianPotentialTangentRequest``."""
        if "d_d2phi_dr2" in missing or "d_d2phi_dr_dxi" in missing:
            raise NotImplementedError(
                "PySCF Hessians support PV and scalar-trace PCMs (CPCM, COSMO, IEFPCM, SSVPE, LPB); "
                "normal-trace PCMs requiring second spatial potential tangents are not supported"
            )

        ndir, ngrid, _ = d_xyz.shape
        density = self.density_contracted(dm)
        delta, c0, c1, r2, exp_x2 = self._nuclear_radial()
        charges = np.asarray(self.mol.atom_charges(), dtype=np.float64)

        d_phi = np.empty((ndir, ngrid), dtype=np.float64) if "d_phi" in missing else None
        d_dphi_dr = np.empty((ndir, ngrid, 3), dtype=np.float64) if "d_dphi_dr" in missing else None
        d_dphi_dxi = np.empty((ndir, ngrid), dtype=np.float64) if "d_dphi_dxi" in missing else None

        if self.xi is not None:
            dpsi_dxi = 2.0 * exp_x2 / _SQRT_PI
            dc0_dxi = 4.0 * (self.xi[:, None] ** 2) * exp_x2 / _SQRT_PI
        else:
            dpsi_dxi = dc0_dxi = None

        for j in range(ndir):
            v_r = dirs[j]  # (natm, 3)
            t = v_r[self.atom_of_ao].T  # (3, nao)
            d_r = d_xyz[j].T  # (3, ngrid)
            dm_dir = None if dm_dirs is None else dm_dirs[j]
            has_p = dm_dir is not None and bool(np.any(dm_dir))
            u = d_xyz[j][:, None, :] - v_r[None, :, :]  # (ngrid, natm, 3)
            ru = _contract("iAk,iAk->iA", delta, u)
            dxi_j = None if d_xi is None else d_xi[j]

            if d_phi is not None:
                val = _contract("A,iA,iA->i", charges, -c0, ru)
                val += 2.0 * _contract("cu,cui->i", t, density.bra)
                val += _contract("ci,ci->i", density.field, d_r)
                if has_p:
                    val -= _contract("uvi,uv->i", self.i0, dm_dir)
                if dxi_j is not None:
                    val += _contract("A,iA->i", charges, dpsi_dxi) * dxi_j
                    val += density.dphi_el_dxi * dxi_j
                d_phi[j] = val

            if d_dphi_dr is not None:
                vec = -_contract("A,iA,iAk->ik", charges, c0, u)
                vec -= _contract("A,iA,iAk->ik", charges, c1 * ru, delta)
                vec += (-2.0 * _contract("cu,caui->ai", t, density.bra_point)).T
                vec += (-_contract("aci,ci->ai", density.point, d_r)).T
                if has_p:
                    vec += _contract("auvi,uv->ia", self.ip2, dm_dir)
                if dxi_j is not None:
                    vec -= _contract("A,iA,iAk->ik", charges, dc0_dxi, delta) * dxi_j[:, None]
                    vec += density.d2phi_el_dr_dxi.T * dxi_j[:, None]
                d_dphi_dr[j] = vec

            if d_dphi_dxi is not None:
                if self.xi is None or dxi_j is None:
                    raise RuntimeError("d_dphi_dxi requested without Gaussian widths")
                wval = -_contract("A,iA,iA->i", charges, dc0_dxi, ru)
                wval -= _contract(
                    "A,iA->i", charges, (4.0 * self.xi[:, None] / _SQRT_PI) * r2 * exp_x2
                ) * dxi_j
                wval += _contract("cu,cui->i", t, density.bra_dxi)
                wval += _contract("ai,ai->i", density.d2phi_el_dr_dxi, d_r)
                wval += density.d2phi_el_dxi2 * dxi_j
                if has_p:
                    wval -= 4.0 * _contract("uvi,uv->i", self.g_xi, dm_dir)
                d_dphi_dxi[j] = wval

        outputs: dict[str, np.ndarray] = {}
        if d_phi is not None:
            outputs["d_phi"] = d_phi
        if d_dphi_dr is not None:
            outputs["d_dphi_dr"] = d_dphi_dr
        if d_dphi_dxi is not None:
            outputs["d_dphi_dxi"] = d_dphi_dxi
        return outputs

    def complete_potential(
        self,
        dm: np.ndarray,
        dirs: np.ndarray,
        dm_dirs: Optional[np.ndarray],
        item: PotentialAdjointTangentResponse,
        hvp_nuc: np.ndarray,
        fock_tan: np.ndarray,
    ) -> None:
        """Accumulate host AO potential HVP and Fock tangent from ``PotentialAdjointTangentResponse``."""
        if item.w_dphi_dr is not None or item.dw_dphi_dr is not None:
            raise NotImplementedError(
                "PySCF Hessians do not support normal-trace PCM response items"
            )

        w_phi = np.asarray(item.w_phi, dtype=np.float64)
        charge = self.charge_weighted(w_phi)
        density = self.density_contracted(dm)
        ndir = dirs.shape[0]

        for j in range(ndir):
            v_r = dirs[j]
            t = v_r[self.atom_of_ao].T
            dm_dir = None if dm_dirs is None else dm_dirs[j]
            has_p = dm_dir is not None and bool(np.any(dm_dir))
            dw_phi = np.asarray(item.dw_phi[j], dtype=np.float64)
            d_r = np.asarray(item.d_xyz[j], dtype=np.float64).T  # (3, ngrid)

            # Nuclear HVP electronic terms
            per_ao = 2.0 * _contract("i,cui->cu", dw_phi, density.bra)
            if has_p:
                per_ao += 2.0 * _contract("cuv,uv->cu", charge.bra, dm_dir)
            if np.any(v_r):
                per_ao -= 2.0 * _contract("acuv,cu,uv->au", charge.bra_bra, t, dm)
                per_ao -= 2.0 * _contract("acuv,cv,uv->au", charge.bra_ket, t, dm)
            per_ao -= 2.0 * _contract("acui,ci->au", density.bra_point, w_phi * d_r)

            fock = -_contract("i,uvi->uv", dw_phi, self.i0)
            if np.any(v_r):
                cross = _contract("cu,cuv->uv", t, charge.bra)
                fock += cross + cross.T
            fock += _contract("ci,cuvi->uv", w_phi * d_r, self.ip2)

            if item.d_xi is not None:
                if self.xi is None:
                    raise RuntimeError("d_xi response received without Gaussian widths")
                w_dxi = w_phi * np.asarray(item.d_xi[j], dtype=np.float64)
                per_ao += _contract("i,kui->ku", w_dxi, density.bra_dxi)
                fock -= 4.0 * _contract("i,uvi->uv", w_dxi, self.g_xi)

            hvp_nuc[j] += self._per_atom(per_ao)
            fock_tan[j] += fock


def _same_xi(a: Optional[np.ndarray], b: Optional[np.ndarray]) -> bool:
    if a is None or b is None:
        return a is b
    return np.array_equal(a, b)


def _shared_tables(shared: Optional[dict], host, coords: np.ndarray, xi: Optional[np.ndarray], *,
                   need_density: bool, need_potential: bool):
    """Build or reuse geometry-only host tables for one Hessian solve."""
    if shared is None:
        shared = {}
    elif (
        shared.get("host") is not host
        or not np.array_equal(shared.get("coords"), coords)
        or not _same_xi(shared.get("xi"), xi)
    ):
        shared.clear()
    shared["host"] = host
    shared["coords"] = np.array(coords, dtype=np.float64, copy=True)
    shared["xi"] = None if xi is None else np.array(xi, dtype=np.float64, copy=True)

    density_tables = None
    potential_tables = None
    if need_density:
        density_tables = shared.get("density")
        if density_tables is None:
            density_tables = _DensityTables(host, coords)
            shared["density"] = density_tables
    if need_potential:
        potential_tables = shared.get("potential")
        if potential_tables is None:
            potential_tables = _ElectrostaticTables(host, coords, xi)
            shared["potential"] = potential_tables
    return density_tables, potential_tables


def _validate_second_order_model(solvation) -> None:
    """Reject components/cavities that PySCF second-order host coupling does not support."""
    if isinstance(solvation.configuration, ISwiG):
        raise NotImplementedError("PySCF analytic Hessians do not support ISwiG cavities")
    for item in solvation.response:
        if isinstance(item, GaussianAmplitudeResponse):
            raise NotImplementedError("PySCF analytic Hessians do not support GOSTSHYP")
        if isinstance(item, PotentialAdjointResponse) and item.w_dphi_dr is not None:
            raise NotImplementedError(
                "PySCF analytic Hessians support PV and scalar-trace PCMs only; "
                "normal-trace PCMs are not supported"
            )


def _normalize_hvp_directions(solvation, dirs, dm_dirs):
    """Validate and normalize nuclear and density perturbation batches."""
    mol = solvation.mol
    natm = mol.natm
    nao = mol.nao_nr()
    if dirs is None and dm_dirs is None:
        raise ValueError("provide at least one of dirs or dm_dirs")

    dm_dirs_arr: Optional[np.ndarray] = None
    if dm_dirs is not None:
        raw = np.asarray(dm_dirs)
        if np.iscomplexobj(raw):
            raise NotImplementedError("MOIST supports real density matrices only")
        raw = np.asarray(raw, dtype=np.float64)
        if raw.ndim == 4 and raw.shape[1:] == (2, nao, nao):
            raw = raw.sum(axis=1)
        if raw.ndim != 3 or raw.shape[1:] != (nao, nao):
            raise ValueError(
                f"Expected dm_dirs shape (ndir, {nao}, {nao}) or (ndir, 2, {nao}, {nao}), got {raw.shape}"
            )
        if not np.all(np.isfinite(raw)):
            raise ValueError("dm_dirs must contain finite values")
        dm_dirs_arr = np.ascontiguousarray(raw)

    if dirs is None:
        ndir = int(dm_dirs_arr.shape[0])
        dirs_arr = np.zeros((ndir, natm, 3), dtype=np.float64)
    else:
        dirs_arr = np.ascontiguousarray(dirs, dtype=np.float64)
        if dirs_arr.ndim != 3 or dirs_arr.shape[1:] != (natm, 3):
            raise ValueError(f"Expected dirs shape (ndir, {natm}, 3), got {dirs_arr.shape}")
        if not np.all(np.isfinite(dirs_arr)):
            raise ValueError("dirs must contain finite values")
        if dm_dirs_arr is not None and dm_dirs_arr.shape[0] != dirs_arr.shape[0]:
            raise ValueError("dirs and dm_dirs must have the same number of directions")

    return dirs_arr, dm_dirs_arr


def solvation_hvp(
    solvation,
    dm,
    dirs=None,
    dm_dirs=None,
    *,
    shared: Optional[dict] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute directional nuclear HVP ``(ndir, natm, 3)`` and Fock tangent ``(ndir, nao, nao)``."""
    solvation.evaluate(dm)
    _validate_second_order_model(solvation)

    dirs_arr, dm_dirs_arr = _normalize_hvp_directions(solvation, dirs, dm_dirs)
    ndir = int(dirs_arr.shape[0])
    natm = solvation.mol.natm
    nao = solvation.mol.nao_nr()
    if ndir == 0:
        return np.zeros((0, natm, 3), dtype=np.float64), np.zeros((0, nao, nao), dtype=np.float64)

    dm0 = solvation._dm
    solvation.host.dm = dm0
    coords = solvation.model.cavity.xyz
    xi = solvation._xi
    need_density = bool(solvation.configuration.density_dependent)
    need_potential = any(
        isinstance(item, PotentialAdjointResponse) for item in solvation.response
    )
    density_tables, potential_tables = _shared_tables(
        shared, solvation.host, coords, xi,
        need_density=need_density,
        need_potential=need_potential,
    )

    cpl = solvation.coupling
    timer = solvation.context.timer
    model = solvation.model

    # Stage 1: density tangents at fixed base surface
    model.prepare_hvp_geometry(cpl, dirs_arr)
    with timer("HVP geometry coupling"):
        for request in cpl:
            with timer(request.name):
                if isinstance(request, DensityTangentRequest):
                    if density_tables is None:
                        density_tables = _DensityTables(solvation.host, coords)
                    cpl.answer(
                        **density_tables.density_tangent(dm0, dirs_arr, dm_dirs_arr, request.missing)
                    )
                else:
                    raise NotImplementedError(f"PySCF Hessian host cannot answer {request.name!r}")

    # Stage 2: surface tangents + electrostatic tangents
    model.prepare_hvp(cpl)
    with timer("HVP coupling"):
        for request in cpl:
            with timer(request.name):
                if isinstance(
                    request,
                    (NuclearChargesRequest, PointPotentialRequest, GaussianPotentialRequest),
                ):
                    solvation._answer(request)
                elif isinstance(request, PointPotentialTangentRequest):
                    if potential_tables is None:
                        potential_tables = _ElectrostaticTables(solvation.host, coords, None)
                    cpl.answer(
                        **potential_tables.potential_tangent(
                            dm0, request.dirs, dm_dirs_arr, request.d_xyz, None, request.missing
                        )
                    )
                elif isinstance(request, GaussianPotentialTangentRequest):
                    if potential_tables is None:
                        potential_tables = _ElectrostaticTables(solvation.host, coords, request.xi)
                    cpl.answer(
                        **potential_tables.potential_tangent(
                            dm0, request.dirs, dm_dirs_arr, request.d_xyz, request.d_xi,
                            request.missing,
                        )
                    )
                elif isinstance(request, GaussianMomentRequest):
                    raise NotImplementedError("PySCF analytic Hessians do not support GOSTSHYP")
                else:
                    raise NotImplementedError(f"PySCF Hessian host cannot answer {request.name!r}")

    hvp_nuc = np.zeros((ndir, natm, 3), dtype=np.float64)
    response = model.get_hvp(cpl, hvp_nuc)
    fock_tan = np.zeros((ndir, nao, nao), dtype=np.float64)

    with timer("HVP contraction"):
        for item in response:
            with timer(item.name):
                if isinstance(item, PotentialAdjointTangentResponse):
                    if potential_tables is None:
                        potential_tables = _ElectrostaticTables(solvation.host, coords, xi)
                    potential_tables.complete_potential(
                        dm0, dirs_arr, dm_dirs_arr, item, hvp_nuc, fock_tan
                    )
                elif isinstance(item, DensityTangentResponse):
                    if density_tables is None:
                        density_tables = _DensityTables(solvation.host, coords)
                    density_tables.complete_density(
                        dm0, dirs_arr, dm_dirs_arr, item, hvp_nuc, fock_tan
                    )
                else:
                    raise NotImplementedError(f"Unsupported HVP response item {item.name!r}")

    return hvp_nuc, fock_tan


def solvation_hessian(solvation, dm) -> np.ndarray:
    """Assemble fixed-density solvation nuclear Hessian ``(natm, natm, 3, 3)`` from unit HVPs."""
    natm = solvation.mol.natm
    unit_dirs = np.zeros((3 * natm, natm, 3), dtype=np.float64)
    for atom in range(natm):
        for axis in range(3):
            unit_dirs[3 * atom + axis, atom, axis] = 1.0
    hvp_nuc, _ = solvation_hvp(solvation, dm, unit_dirs)
    return hvp_nuc.reshape(natm, 3, natm, 3).transpose(0, 2, 1, 3)


def _has_df_hessian(mean_field) -> bool:
    """Whether PySCF implements second derivatives for this fitted reference."""
    from pyscf.df.df_jk import _DFHF

    fitting = getattr(mean_field, "with_df", None)
    if fitting is None or not isinstance(mean_field, _DFHF):
        return False
    return not type(fitting).__module__.startswith("pyscf.sgx")


def _validate_scf_reference(mean_field) -> bool:
    """Check SCF convergence and supported electronic reference; return whether unrestricted."""
    from pyscf.scf.hf import RHF
    from pyscf.scf.uhf import UHF

    if not mean_field.converged:
        raise ValueError("SCF must converge before evaluating its Hessian")
    if mean_field.mol.has_ecp():
        raise NotImplementedError("Effective core potentials are not supported in MOIST Hessians")
    if np.iscomplexobj(mean_field.mo_coeff):
        raise NotImplementedError("Complex orbitals are not supported in MOIST Hessians")
    if getattr(mean_field, "with_x2c", None) is not None:
        raise NotImplementedError("X2C relativistic references are not supported in MOIST Hessians")
    fitting = getattr(mean_field, "with_df", None)
    if fitting is not None and not _has_df_hessian(mean_field):
        raise NotImplementedError(
            f"{type(fitting).__name__} has no PySCF second-derivative implementation"
        )

    if isinstance(mean_field, UHF):
        for occ in mean_field.mo_occ:
            if np.any((occ != 0) & (occ != 1)):
                raise NotImplementedError("Fractional occupations are not supported in MOIST Hessians")
        return True
    if isinstance(mean_field, RHF):
        if np.any((mean_field.mo_occ != 0) & (mean_field.mo_occ != 2)):
            raise NotImplementedError(
                "Only closed-shell restricted or integer-occupied unrestricted references are supported"
            )
        return False
    raise NotImplementedError(f"Unsupported mean-field class {type(mean_field).__name__!r}")


def _make_moist_hessian(mean_field):
    """Build a PySCF Hessian object coupled to ``mean_field.with_moist``."""
    from pyscf import lib
    from pyscf.hessian import rhf as rhf_hess, uhf as uhf_hess

    unrestricted = _validate_scf_reference(mean_field)
    base_hess = mean_field.undo_moist().Hessian()
    base_hess.base = mean_field
    nat = mean_field.mol.natm

    class _MoistHessian:
        kernel_callback = None

        def _solvent_terms(self, mo_coeff=None, mo_occ=None):
            mf = self.base
            if mo_coeff is None:
                mo_coeff = mf.mo_coeff
            if mo_occ is None:
                mo_occ = mf.mo_occ
            if mf.with_moist.mol is not self.mol:
                mf.with_moist.reset(self.mol)
            dm0 = mf.make_rdm1(mo_coeff, mo_occ)
            cache = getattr(self, "_moist_cache", None)
            if cache is not None and np.array_equal(cache[0], dm0):
                return cache[1], cache[2], cache[3]

            unit = np.zeros((3 * nat, nat, 3), dtype=np.float64)
            for atom in range(nat):
                for axis in range(3):
                    unit[3 * atom + axis, atom, axis] = 1.0
            shared: dict = {}
            solv = mf.with_moist._driver()
            nuclear, fock = solvation_hvp(solv, dm0, unit, shared=shared)
            hrr = nuclear.reshape(nat, 3, nat, 3).transpose(0, 2, 1, 3)
            hrp = fock.reshape(nat, 3, mf.mol.nao_nr(), mf.mol.nao_nr())
            self._moist_cache = (np.array(dm0, dtype=np.float64, copy=True), hrr, hrp, shared)
            return hrr, hrp, shared

        def partial_hess_elec(
            self,
            mo_energy=None,
            mo_coeff=None,
            mo_occ=None,
            atmlst=None,
            max_memory=4000,
            verbose=None,
        ):
            atoms = list(range(nat)) if atmlst is None else list(atmlst)
            hrr, _, _ = self._solvent_terms(mo_coeff, mo_occ)
            elec = super().partial_hess_elec(
                mo_energy, mo_coeff, mo_occ, atoms, max_memory, verbose
            )
            return elec + hrr[np.ix_(atoms, atoms)]

        def make_h1(self, mo_coeff, mo_occ, chkfile=None, atmlst=None, verbose=None):
            h1 = super().make_h1(mo_coeff, mo_occ, chkfile, atmlst, verbose)
            _, hrp, _ = self._solvent_terms(mo_coeff, mo_occ)
            atoms = range(nat) if atmlst is None else atmlst
            if unrestricted:
                h1a, h1b = h1
                for atom in atoms:
                    h1a[atom] += hrp[atom]
                    h1b[atom] += hrp[atom]
                return h1a, h1b
            for atom in atoms:
                h1[atom] += hrp[atom]
            return h1

        def solve_mo1(
            self,
            mo_energy,
            mo_coeff,
            mo_occ,
            h1ao,
            fx=None,
            atmlst=None,
            max_memory=4000,
            verbose=None,
        ):
            if fx is None:
                mf = self.base
                gas_mf = mf.undo_moist()
                _, _, shared = self._solvent_terms(mo_coeff, mo_occ)
                solv = mf.with_moist._driver()
                dm0 = mf.make_rdm1(mo_coeff, mo_occ)

                if unrestricted:
                    gas_response = uhf_hess.gen_vind(gas_mf, mo_coeff, mo_occ)
                    mocca = mo_coeff[0][:, mo_occ[0] > 0]
                    moccb = mo_coeff[1][:, mo_occ[1] > 0]
                    na_rot = mo_coeff[0].shape[1] * mocca.shape[1]
                    nb_rot = mo_coeff[1].shape[1] * moccb.shape[1]

                    def fx(mo1):
                        mo1_2d = np.asarray(mo1).reshape(-1, na_rot + nb_rot)
                        nset = mo1_2d.shape[0]
                        if self.kernel_callback is not None:
                            self.kernel_callback(nset)
                        xa = mo1_2d[:, :na_rot].reshape(nset, mo_coeff[0].shape[1], mocca.shape[1])
                        xb = mo1_2d[:, na_rot:].reshape(nset, mo_coeff[1].shape[1], moccb.shape[1])
                        dma = _contract("um,pmi,vi->puv", mo_coeff[0], xa, mocca)
                        dmb = _contract("um,pmi,vi->puv", mo_coeff[1], xb, moccb)
                        dm1 = (dma + dma.transpose(0, 2, 1)) + (dmb + dmb.transpose(0, 2, 1))
                        _, f1 = solvation_hvp(solv, dm0, dm_dirs=dm1, shared=shared)
                        va = _contract("um,puv,vi->pmi", mo_coeff[0], f1, mocca).reshape(nset, na_rot)
                        vb = _contract("um,puv,vi->pmi", mo_coeff[1], f1, moccb).reshape(nset, nb_rot)
                        return gas_response(mo1) + np.hstack((va, vb)).reshape(np.asarray(mo1).shape)
                else:
                    gas_response = rhf_hess.gen_vind(gas_mf, mo_coeff, mo_occ)
                    mocc = mo_coeff[:, mo_occ > 0]

                    def fx(mo1):
                        rot = np.asarray(mo1).reshape(-1, mo_coeff.shape[1], mocc.shape[1])
                        if self.kernel_callback is not None:
                            self.kernel_callback(len(rot))
                        dm1 = 2.0 * _contract("um,pmi,vi->puv", mo_coeff, rot, mocc)
                        dm1 += dm1.transpose(0, 2, 1)
                        _, f1 = solvation_hvp(solv, dm0, dm_dirs=dm1, shared=shared)
                        solv_vo = _contract("um,puv,vi->pmi", mo_coeff, f1, mocc)
                        return gas_response(mo1) + solv_vo.reshape(np.asarray(mo1).shape)

            return super().solve_mo1(
                mo_energy, mo_coeff, mo_occ, h1ao, fx, atmlst, max_memory, verbose
            )

        def gen_hop(self, mo_energy=None, mo_coeff=None, mo_occ=None, verbose=None):
            dense = self.kernel(mo_energy, mo_coeff, mo_occ, atmlst=range(nat))
            matrix = dense.transpose(0, 2, 1, 3).reshape(3 * nat, 3 * nat)
            return matrix.dot, matrix.diagonal().copy()

    return lib.set_class(base_hess, (_MoistHessian, base_hess.__class__), "MOIST" + base_hess.__class__.__name__)
