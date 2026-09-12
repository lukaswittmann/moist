"""Nuclear Hessian utilities, using PySCF's atom/atom/Cartesian/Cartesian order."""

from concurrent.futures import ThreadPoolExecutor
from itertools import permutations, product
from typing import Callable

import numpy as np


def _contract(subscripts, *operands):
    """``np.einsum`` with the pairwise path, and BLAS's dirty flags suppressed.

    Two things that have to travel together. Without ``optimize`` numpy
    evaluates a multi-operand expression as one naive nested loop over every
    index, which for the second-order host tensors -- shaped
    ``(3, 3, nao, nao, ngrid)`` -- is the whole cost of a solvent Hessian.
    Measured at nao=37, ngrid~1390: ``i,acuvi,cu,uv->au`` 43.6 -> 2.1 ms,
    ``i,iu,iv->uv`` 1.85 -> 0.06 ms, ``um,pmi,vi->puv`` over 39 directions
    41.4 -> 0.07 ms; at nao=119 the last is 1335 -> 0.40 ms.

    Taking the pairwise path routes the products through BLAS, which brings the
    padding-lane artifact :meth:`_LevelSetFunctional._project` documents: the
    SIMD tail reads lanes that raise divide-by-zero, overflow and invalid from a
    product that performs none of them. The result is unaffected, so the flags
    are silenced here exactly as they are there.

    Not every contraction wants this. A diagonal one (``i,iu,iu->u``) has no
    pairwise factorization to find, and the path search then costs more than it
    saves -- those call sites deliberately still use ``np.einsum``.
    """

    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        return np.einsum(subscripts, *operands, optimize=True)


def _map_directions(body, n):
    """Apply ``body(j)`` to each of ``n`` directions, over threads where they pay.

    The three per-direction loops of the host exchange -- the two callbacks and
    the completion -- are numpy reductions over the host tensors, which release
    the GIL for the whole of their work, so threads here buy real cores rather
    than interleaving. The directions share only read-only operands and write
    disjoint slices of the output, so nothing about the result depends on the
    schedule: it is bit-for-bit what the serial loop produces, which is worth
    more here than the speed, since a Hessian that changed with the thread count
    would be impossible to validate against finite differences.

    The pool is sized by PySCF's thread count -- ``omp_get_max_threads``, the
    same one libmoist reads -- so that ``-t N`` governs this layer too instead of
    the process quietly running more threads than were asked for. The pool is
    built per call rather than kept: these loops run at most a few hundred times
    in a Hessian, against which a pool's lifetime is not worth owning.
    """
    from pyscf import lib

    nthread = min(n, max(1, lib.num_threads()))
    if nthread <= 1:
        for j in range(n):
            body(j)
        return
    with ThreadPoolExecutor(nthread) as pool:
        # Consumed inside the pool's context, so that a direction that raises
        # does so here instead of being discarded at shutdown.
        list(pool.map(body, range(n)))


class _DensityDerivatives:
    """Partial density jets in nuclear and independent symmetric AO-DM coordinates.

    This dense reference implementation trades memory for a direct chain rule.
    Nuclear coordinates come first, followed by the upper triangle of P; an
    off-diagonal parameter changes both symmetric entries of P by one.
    """

    def __init__(self, host, dm):
        self.host = host
        self.dm = dm
        self.nr = 3 * host.mol.natm
        self.rows, self.cols = np.triu_indices(host.mol.nao)
        self.n = self.nr + len(self.rows)
        self.dp = np.zeros((self.n, host.mol.nao, host.mol.nao))
        for p, (u, v) in enumerate(zip(self.rows, self.cols), self.nr):
            self.dp[p, u, v] = self.dp[p, v, u] = 1
        self.dirs = np.zeros((3, host.mol.natm, self.n), order="F")
        for p in range(self.nr):
            self.dirs[p % 3, p // 3, p] = 1

    def level_set_jets(self, coords):
        """Fixed-point level-set jets requested by a density-defined cavity."""
        ao = self.host._ao(coords, 4)
        for i in range(len(coords)):
            yield self.point(ao[:, i])

    def electrostatic_potential(self, coords):
        return self.host.surface_potential(coords)

    def point(self, ao):
        from .pyscf import _component_index

        cache = {}

        def basis(axes, order):
            key = (tuple(sorted(axes)), order)
            if key in cache:
                return cache[key]
            val = ao[_component_index(axes)]
            first = np.zeros((self.n, len(val)))
            second = np.zeros((self.n, self.n, len(val))) if order == 2 else None
            if order:
                for atom, (_, _, lo, hi) in enumerate(self.host._aoslice):
                    for a in range(3):
                        p = 3 * atom + a
                        first[p, lo:hi] = -ao[_component_index(axes + (a,))][lo:hi]
                        if order == 2:
                            for b in range(3):
                                second[p, 3 * atom + b, lo:hi] = ao[
                                    _component_index(axes + (a, b))
                                ][lo:hi]
            cache[key] = val, first, second
            return cache[key]

        def density(axes, order):
            value = 0.0
            first = np.zeros(self.n)
            second = np.zeros((self.n, self.n))
            # Leibniz's rule includes multiplicity through the labelled subsets.
            for mask in product((False, True), repeat=len(axes)):
                left = tuple(a for a, select in zip(axes, mask) if select)
                right = tuple(a for a, select in zip(axes, mask) if not select)
                l0, l1, l2 = basis(left, order)
                r0, r1, r2 = basis(right, order)
                value += l0 @ self.dm @ r0
                if order:
                    first += l1 @ self.dm @ r0 + r1 @ self.dm @ l0
                    first += np.einsum("u,puv,v->p", l0, self.dp, r0)
                if order == 2:
                    second += np.einsum("pqu,u->pq", l2, self.dm @ r0)
                    second += np.einsum("pqu,u->pq", r2, self.dm @ l0)
                    cross = np.einsum("pu,uv,qv->pq", l1, self.dm, r1)
                    cross += np.einsum("pu,quv,v->pq", l1, self.dp, r0)
                    cross += np.einsum("u,puv,qv->pq", l0, self.dp, r1)
                    second += cross + cross.T
            return value, first, second

        jets = [[], [], []]
        density_cache = {}
        for spatial_order in range(5):
            # Symmetric spatial tensors have identical C and Fortran packing.
            for axes in product(range(3), repeat=spatial_order):
                key = tuple(sorted(axes))
                if key not in density_cache:
                    density_cache[key] = density(axes, min(2, 4 - spatial_order))
                val, d1, d2 = density_cache[key]
                jets[0].append(val)
                if spatial_order <= 3:
                    jets[1].append(d1)
                if spatial_order <= 2:
                    jets[2].append(d2)
        jet, jet1, jet2 = [-self.host.scale * np.asarray(x) for x in jets]
        jet[0] += self.host.scale * self.host.rho_iso
        return jet, jet1, jet2

    def potential(self, point, surface1, surface2):
        """Molecular potential derivatives, including moving operator and AO centers."""
        mol = self.host.mol
        xyz1 = surface1[:3].T
        xyz2 = surface2[:3].transpose(1, 2, 0)
        delta = point - mol.atom_coords()
        r = np.linalg.norm(delta, axis=1)
        motion = xyz1[:, None, :] - self.dirs.transpose(2, 1, 0)
        dot = np.einsum("ak,pak->pa", delta, motion)
        z = mol.atom_charges()
        phi1 = -np.einsum("a,pa,a->p", z, dot, r**-3)
        phi2 = 3 * np.einsum("a,pa,qa,a->pq", z, dot, dot, r**-5)
        phi2 -= np.einsum("a,pak,qak,a->pq", z, motion, motion, r**-3)
        phi2 -= np.einsum("a,ak,pqk,a->pq", z, delta, xyz2, r**-3)

        # Translational invariance converts operator motion into the sum of
        # bra and ket spatial derivatives. libcint supplies exact AO integrals.
        with mol.with_rinv_origin(point):
            v = mol.intor("int1e_rinv")
            ip = mol.intor("int1e_iprinv", comp=3)
            aa = mol.intor("int1e_ipiprinv", comp=9).reshape(3, 3, mol.nao, mol.nao)
            ab = mol.intor("int1e_iprinvip", comp=9).reshape(3, 3, mol.nao, mol.nao)
        t = np.empty((self.n, 3, mol.nao))
        for atom, (_, _, lo, hi) in enumerate(self.host._aoslice):
            t[:, :, lo:hi] = motion[:, atom, :, None]
        v1 = np.einsum("pku,kuv->puv", t, ip)
        v1 += v1.transpose(0, 2, 1).copy()
        phi1 -= np.einsum("puv,uv->p", v1, self.dm)
        phi1 -= np.einsum("uv,puv->p", v, self.dp)
        phi2 -= 2 * np.einsum("pku,qlu,kluv,uv->pq", t, t, aa, self.dm, optimize=True)
        phi2 -= 2 * np.einsum("pku,qlv,kluv,uv->pq", t, t, ab, self.dm, optimize=True)
        phi2 -= 2 * np.einsum("pqk,kuv,uv->pq", xyz2, ip, self.dm)
        cross = np.einsum("puv,quv->pq", v1, self.dp)
        phi2 -= cross + cross.T
        return phi1, phi2

    def unpack_fock(self, derivatives):
        """Convert derivatives in symmetric P coordinates to an AO Fock matrix."""
        out = np.zeros(derivatives.shape[:-1] + self.dm.shape)
        values = derivatives / np.where(self.rows == self.cols, 1, 2)
        out[..., self.rows, self.cols] = values
        out[..., self.cols, self.rows] = values
        return out


# -----------------------------------------------------------------------------
# The directional path: one protocol for Fortran, C and Python
# -----------------------------------------------------------------------------


def _sorted_multi_indices(order):
    """Sorted Cartesian multi-indices of one order with their multiplicities."""
    seen = {}
    for axes in product(range(3), repeat=order):
        key = tuple(sorted(axes))
        seen[key] = seen.get(key, 0) + 1
    return list(seen.items())


def _symmetric_outer(tensor, vector):
    """Symmetrised outer product of a symmetric ``(3,)*k + (ngrid,)`` tensor and a ``(3, ngrid)`` vector."""
    k = tensor.ndim - 1
    outer = np.einsum("...i,ci->...ci", tensor, vector)
    if k == 0:
        return outer
    total = np.zeros_like(outer)
    axes = list(range(k + 1))
    for position in range(k + 1):
        order = axes[:position] + [k] + axes[position:k]
        total += np.transpose(outer, order + [k + 1])
    return total / (k + 1)


class _LevelSetFunctional:
    """Weighted functionals of the level-set jet and their host derivatives.

    With ``S = scale (rho_iso - rho)`` and ``rho = sum P_uv chi_u chi_v``, a
    set of symmetric weights ``W_k(i)`` on the spatial orders ``k = 0..3`` at
    the surface points defines ``Lambda(R, P) = sum_i sum_k W_k(i) : d^k S(r_i)``.
    Everything the host completes a gradient, a Fock matrix or their
    directional tangents from is a derivative of such a functional with the
    weights frozen, so one Leibniz expansion over the AO derivative table
    serves them all. AO derivatives to fourth order are needed and taken once.
    """

    def __init__(self, host, coords):
        self.host = host
        self.mol = host.mol
        self.scale = host.scale
        self.ao = host._ao(coords, 4)
        self.atom_of_ao = np.empty(self.mol.nao, dtype=int)
        for atom, (_, _, lo, hi) in enumerate(host._aoslice):
            self.atom_of_ao[lo:hi] = atom

    def _terms(self, weights):
        """Yield ``(w(ngrid), alpha, beta)`` for every weighted Leibniz split."""
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

    def _ao_at(self, axes):
        from .pyscf import _component_index

        return self.ao[_component_index(axes)]

    @staticmethod
    def _project(left, dm):
        """``left @ dm`` with BLAS's dirty floating-point flags suppressed.

        The same artifact :meth:`PySCFHost.density` documents: for these shapes
        the SIMD tail reads padding lanes, so numpy reports divide-by-zero,
        overflow and invalid from a product that performs none of them. It
        first appears at def2-SVP sizes and not at STO-3G, which is why only
        the density callback had to guard against it so far. The result is
        bit-identical to the BLAS-free einsum, so the flags are cleared rather
        than paying the einsum path's cost.
        """
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            return left @ dm

    def _per_atom(self, per_ao):
        out = np.zeros((3, self.mol.natm))
        for atom, (_, _, lo, hi) in enumerate(self.host._aoslice):
            out[:, atom] = per_ao[:, lo:hi].sum(axis=1)
        return out

    def fock(self, weights):
        """``d Lambda / d P_uv`` at frozen weights, ``(nao, nao)``."""
        fock = np.zeros((self.mol.nao, self.mol.nao))
        for w, alpha, beta in self._terms(weights):
            fock += _contract("i,iu,iv->uv", w, self._ao_at(alpha), self._ao_at(beta))
        return -self.scale * fock

    def gradient(self, weights, dm):
        """``d Lambda / d R_A`` at frozen weights through the AO centres, ``(3, natm)``."""
        per_ao = np.zeros((3, self.mol.nao))
        for w, alpha, beta in self._terms(weights):
            right = self._project(self._ao_at(beta), dm)
            for c in range(3):
                per_ao[c] += np.einsum("i,iu,iu->u", w, self._ao_at(alpha + (c,)), right)
        return 2.0 * self.scale * self._per_atom(per_ao)

    def fock_centre_tangent(self, weights, dirs):
        """Directional derivative of :meth:`fock` along nuclear directions ``dirs (3, natm)``."""
        t = dirs[:, self.atom_of_ao]
        cross = np.zeros((self.mol.nao, self.mol.nao))
        for w, alpha, beta in self._terms(weights):
            left = sum(t[c] * self._ao_at(alpha + (c,)) for c in range(3))
            cross += _contract("i,iu,iv->uv", w, left, self._ao_at(beta))
        return self.scale * (cross + cross.T)

    def gradient_centre_tangent(self, weights, dm, dirs):
        """Directional derivative of :meth:`gradient` along nuclear directions ``dirs (3, natm)``."""
        t = dirs[:, self.atom_of_ao]
        per_ao = np.zeros((3, self.mol.nao))
        for w, alpha, beta in self._terms(weights):
            right = self._project(self._ao_at(beta), dm)
            right_moved = self._project(
                sum(t[d] * self._ao_at(beta + (d,)) for d in range(3)), dm
            )
            for c in range(3):
                left_moved = sum(t[d] * self._ao_at(alpha + (c, d)) for d in range(3))
                per_ao[c] -= np.einsum("i,iu,iu->u", w, left_moved, right)
                per_ao[c] -= np.einsum("i,iu,iu->u", w, self._ao_at(alpha + (c,)), right_moved)
        return 2.0 * self.scale * self._per_atom(per_ao)

    def jet_tangent(self, dm, dm_dir, dirs):
        """Partial tangents of the scaled level-set jet at the fixed points, orders 0..3.

        ``dm_dir`` is the density direction and ``dirs (3, natm)`` the nuclear
        one, along which the AO centres move; returns ``(40, ngrid)`` packed
        by spatial order as full Cartesian tensors in Fortran order.
        """
        ngrid = self.ao.shape[1]
        t = dirs[:, self.atom_of_ao]
        jets = np.zeros((40, ngrid))
        offsets = {0: 0, 1: 1, 2: 4, 3: 13}
        for order in range(4):
            for key, _ in _sorted_multi_indices(order):
                value = np.zeros(ngrid)
                for mask in product((False, True), repeat=order):
                    alpha = tuple(a for a, take in zip(key, mask) if take)
                    beta = tuple(a for a, take in zip(key, mask) if not take)
                    left = self._ao_at(alpha)
                    right = self._ao_at(beta)
                    if dm_dir is not None:
                        value += np.einsum("iu,iu->i", left, self._project(right, dm_dir))
                    if np.any(dirs):
                        moved = sum(t[d] * self._ao_at(alpha + (d,)) for d in range(3))
                        value -= 2.0 * np.einsum(
                            "iu,iu->i", moved, self._project(right, dm)
                        )
                value *= -self.scale
                for axes in set(permutations(key)):
                    index = offsets[order] + sum(a * 3**p for p, a in enumerate(axes))
                    jets[index] = value
        return jets


class _ChargeWeighted:
    """The integrals with the point index summed against a fixed charge set.

    What is left is a tensor over AO pairs alone, which is the whole point: the
    grid is the long index, and every direction of a Hessian-vector call sees
    the same converged charges, so summing it here removes it from the inner
    loop rather than re-walking it per direction.
    """

    def __init__(self, ints, q):
        #: ``sum_i q_i (d_bra uv|i)``, ``(3, nao, nao)``
        self.bra = _contract("i,cuvi->cuv", q, ints.ip1)
        #: ``sum_i q_i (d_bra d_bra uv|i)``, ``(3, 3, nao, nao)``
        self.bra_bra = _contract("i,acuvi->acuv", q, ints.ipip1)
        #: ``sum_i q_i (d_bra uv d_ket|i)``, ``(3, 3, nao, nao)``
        self.bra_ket = _contract("i,acuvi->acuv", q, ints.ipvip1)


class _DensityContracted:
    """The integrals with the AO pair summed against a fixed density.

    The mirror image of :class:`_ChargeWeighted`, and the one that matters most:
    the tensors carrying two AO indices *and* the grid are the large ones, and
    the converged density contracts them down to the grid. Where a derivative
    sits on the bra its AO index is left open, because that is the index a
    nuclear direction is spread over.
    """

    def __init__(self, ints, dm):
        #: ``grad phi_el(r_i)`` up to sign, ``(3, ngrid)``
        self.field = _contract("cuvi,uv->ci", ints.ip2, dm)
        #: ``sum_v (d_bra uv|i) P_uv``, ``(3, nao, ngrid)``
        self.bra = _contract("cuvi,uv->cui", ints.ip1, dm)
        #: ``sum_v (d_bra uv d_point|i) P_uv``, ``(3, 3, nao, ngrid)``
        self.bra_point = _contract("acuvi,uv->acui", ints.ip1ip2, dm)
        # d^2 (uv|i) / dr_a dr_c with d/dr = -(d/dR_u + d/dR_v), contracted with
        # the density; see the class docstring of _PotentialIntegrals for why
        # ipip2 is not used. Assembled from the AO-centre derivatives already
        # contracted rather than from the summed integral, which at def2-SVP
        # sizes is 1.4 GB and is only ever wanted against this density.
        bra_bra = _contract("acuvi,uv->aci", ints.ipip1, dm)
        ket_ket = _contract("acuvi,vu->aci", ints.ipip1, dm)
        bra_ket = _contract("acuvi,uv->aci", ints.ipvip1, dm)
        #: ``sum_uv d^2 (uv|i)/dr_a dr_c P_uv``, ``(3, 3, ngrid)``
        self.point = bra_bra + ket_ket + bra_ket + bra_ket.transpose(1, 0, 2)


class _PotentialIntegrals:
    """Point-charge integrals at the surface points and the electronic potential terms.

    ``(uv|i)`` and its derivatives from PySCF's three-centre integrals over a
    charge basis at the surface points; the electronic potential is
    ``phi_el(r_i) = -sum_uv P_uv (uv|i)``. Derivative conventions, all checked
    against finite differences: the gradient of the potential at the point is
    the plain ``ip2`` contraction, an AO centre moving is minus the bra
    derivative (``ip1``, ``ipip1``, ``ipvip1``), and the point moving in a
    derivative integral is minus the mixed integral (``ip1ip2``).

    The second derivative with respect to the point is **not** taken from
    ``int3c2e_ipip2``. A charge basis is a Gaussian of exponent 1e16 standing
    in for a point charge, and differentiating it twice leaves a contamination
    that survives the limit: measured against finite differences the
    off-diagonal components are exact while every diagonal one is wrong by
    order the value itself. Translational invariance of the three-centre
    integral gives the same tensor from derivatives of the *AO* centres, which
    are ordinary Gaussians, and that form is exact.
    """

    def __init__(self, mol, coords):
        from pyscf import df, gto

        self.mol = mol
        nao = mol.nao
        charges = gto.fakemol_for_charges(coords)

        def ints(name):
            return df.incore.aux_e2(mol, charges, intor=name, aosym="s1")

        self.i0 = ints("int3c2e")
        self.ip1 = ints("int3c2e_ip1")
        self.ip2 = ints("int3c2e_ip2")
        self.ipip1 = ints("int3c2e_ipip1").reshape(3, 3, nao, nao, -1)
        self.ipvip1 = ints("int3c2e_ipvip1").reshape(3, 3, nao, nao, -1)
        self.ip1ip2 = ints("int3c2e_ip1ip2").reshape(3, 3, nao, nao, -1)
        self._reduced = {}
        self.aoslice = mol.aoslice_by_atom()
        self.atom_of_ao = np.empty(nao, dtype=int)
        for atom, (_, _, lo, hi) in enumerate(self.aoslice):
            self.atom_of_ao[lo:hi] = atom

    def _per_atom(self, per_ao):
        out = np.zeros((3, self.mol.natm))
        for atom, (_, _, lo, hi) in enumerate(self.aoslice):
            out[:, atom] = per_ao[:, lo:hi].sum(axis=1)
        return out

    def _reduce(self, slot, key, build):
        """One-slot cache of a reduction over an operand that is fixed here.

        The integrals outlive a single Hessian-vector call -- they are kept in
        the solve's ``shared`` dict -- and both the converged charges and the
        converged density are the same for every call in a solve, so each
        reduction is built once for the whole Hessian and read by the nuclear
        block and every coupled-perturbed iteration after it.

        The key is compared by value rather than identity: the coupling hands
        out a fresh density array per call, so identity would rebuild every
        time, and a reduction left stale against a density that had actually
        changed would return a wrong Hessian rather than fail. One slot is
        enough because these are the operands that do *not* vary; a caller whose
        operand does vary has a separate entry point below.
        """
        cached = self._reduced.get(slot)
        if cached is not None and np.array_equal(cached[0], key):
            return cached[1]
        value = build()
        self._reduced[slot] = (np.array(key, dtype=float, copy=True), value)
        return value

    def charge_weighted(self, q):
        """The point index summed against ``q``; see :class:`_ChargeWeighted`."""
        return self._reduce("q", q, lambda: _ChargeWeighted(self, q))

    def density_contracted(self, dm):
        """The AO pair summed against ``dm``; see :class:`_DensityContracted`."""
        return self._reduce("dm", dm, lambda: _DensityContracted(self, dm))

    def potential(self, dm):
        return -_contract("uvi,uv->i", self.i0, dm)

    def field(self, dm):
        """``grad phi_el(r_i)``, ``(3, ngrid)``."""
        return self.density_contracted(dm).field

    def gradient(self, q, dm):
        """``sum_i q_i d phi_el(r_i)/dR_A``, ``(3, natm)``.

        Factorized with the density as the fixed operand, which is what the
        completion wants: there it is the charges that carry the direction. For
        the other way round use :meth:`gradient_density_tangent`.
        """
        per_ao = 2.0 * _contract("i,cui->cu", q, self.density_contracted(dm).bra)
        return self._per_atom(per_ao)

    def gradient_density_tangent(self, q, dm_dir):
        """:meth:`gradient` with the charges fixed and the density the direction.

        The same bilinear form, factorized the other way round. Splitting it
        from :meth:`gradient` is not cosmetic: the two call sites disagree about
        which operand repeats, and a single entry point would rebuild one of the
        two reductions on every direction.
        """
        per_ao = 2.0 * _contract("cuv,uv->cu", self.charge_weighted(q).bra, dm_dir)
        return self._per_atom(per_ao)

    def gradient_centre_tangent(self, q, dm, dirs):
        """Directional derivative of :meth:`gradient` along ``dirs (3, natm)`` through the AO centres."""
        t = dirs[:, self.atom_of_ao]
        charge = self.charge_weighted(q)
        per_ao = -2.0 * _contract("acuv,cu,uv->au", charge.bra_bra, t, dm)
        per_ao -= 2.0 * _contract("acuv,cv,uv->au", charge.bra_ket, t, dm)
        return self._per_atom(per_ao)

    def gradient_point_tangent(self, q, dm, d_xyz):
        """Directional derivative of :meth:`gradient` through the surface points moving by ``d_xyz``."""
        # q rides on the direction rather than on the integral: the alternative
        # is a charge-weighted copy of ip1ip2, which is the tensor itself again.
        per_ao = -2.0 * _contract("acui,ci->au", self.density_contracted(dm).bra_point, q * d_xyz)
        return self._per_atom(per_ao)

    def fock(self, q):
        """``-sum_i q_i (uv|i)``."""
        return -_contract("i,uvi->uv", q, self.i0)

    def fock_centre_tangent(self, q, dirs):
        t = dirs[:, self.atom_of_ao]
        cross = _contract("cu,cuv->uv", t, self.charge_weighted(q).bra)
        return cross + cross.T

    def fock_point_tangent(self, q, d_xyz):
        return _contract("ci,cuvi->uv", q * d_xyz, self.ip2)

    def potential_tangent(self, dm, dm_dir, dirs, d_xyz):
        """Total tangent of ``phi_el`` at the moving points, ``(ngrid,)``."""
        t = dirs[:, self.atom_of_ao]
        density = self.density_contracted(dm)
        dphi = 2.0 * _contract("cu,cui->i", t, density.bra)
        if dm_dir is not None:
            dphi -= _contract("uvi,uv->i", self.i0, dm_dir)
        dphi += np.einsum("ci,ci->i", density.field, d_xyz)
        return dphi

    def field_tangent(self, dm, dm_dir, dirs, d_xyz):
        """Total tangent of ``grad phi_el`` at the moving points, ``(3, ngrid)``."""
        t = dirs[:, self.atom_of_ao]
        density = self.density_contracted(dm)
        dfield = -2.0 * _contract("cu,caui->ai", t, density.bra_point)
        dfield -= _contract("aci,ci->ai", density.point, d_xyz)
        if dm_dir is not None:
            # The one term of the exchange no reduction reaches: the density is
            # the direction here, so it is threads or nothing.
            dfield += _contract("auvi,uv->ai", self.ip2, dm_dir)
        return dfield


def _shared_host_terms(shared, host, coords, *, level_set, potential):
    """The geometry-only halves of the host exchange, built once per solve.

    ``_LevelSetFunctional`` and ``_PotentialIntegrals`` are functions of the
    host and the surface points alone -- the density enters only through their
    methods -- but a :class:`DirectionalHost` is constructed per
    Hessian-vector call, so without a cache the six ``aux_e2`` passes and the
    fourth-order AO table are rebuilt for the nuclear block and again on every
    coupled-perturbed iteration.

    ``shared`` is a plain dict owned by the caller that spans one solve.
    Passing ``None`` keeps the old behaviour of building them fresh, which is
    what a one-off call wants. The cache is keyed on the host object and the
    point coordinates and rebuilt when either changes: a surface that moved
    under a stale table would return a wrong Hessian rather than fail, so the
    comparison is worth its cost next to the integrals it guards.
    """

    if shared is None:
        shared = {}
    elif shared.get("host") is not host or not np.array_equal(shared.get("coords"), coords):
        shared.clear()
    shared["host"] = host
    shared["coords"] = coords
    if level_set and "lsf" not in shared:
        shared["lsf"] = _LevelSetFunctional(host, coords)
    if potential and "ints" not in shared:
        shared["ints"] = _PotentialIntegrals(host.mol, coords)
    return (shared["lsf"] if level_set else None), (shared["ints"] if potential else None)


class DirectionalHost:
    """PySCF's side of the second-order host exchange, and its completion.

    One object per Hessian-vector call: it holds the direction set -- the
    nuclear parts ``dirs (3, natm, n)`` moist sees and the density parts
    ``dm_dirs (n, nao, nao)`` it does not -- answers moist's two callbacks
    along any block of it, and completes moist's part of the products with
    the terms that run through the AO basis from the response tangent.
    Everything the completion needs is a first or second derivative of the
    frozen-weight functional ``Lambda(R, P) = sum_i q_i phi_el(r_i) + sum_i
    omega_i . jet_i(r_i)`` whose first derivatives are the host's gradient and
    Fock completions.
    """

    def __init__(self, host, dm, coords, dirs, dm_dirs=None, *, level_set=True,
                 potential=True, shared=None):
        self.host = host
        self.dm = np.asarray(dm)
        self.coords = np.asarray(coords)
        self.dirs = np.asarray(dirs, dtype=float)
        n = self.dirs.shape[2]
        self.dm_dirs = None if dm_dirs is None else np.asarray(dm_dirs, dtype=float)
        if self.dm_dirs is not None and self.dm_dirs.shape != (n,) + self.dm.shape:
            raise ValueError("dm_dirs must have shape (ndir, nao, nao)")
        self.lsf, self.ints = _shared_host_terms(
            shared, host, self.coords, level_set=level_set, potential=potential
        )

    def _dm_dir(self, j):
        if self.dm_dirs is None:
            return None
        direction = self.dm_dirs[j]
        return direction if np.any(direction) else None

    # -- the callbacks -------------------------------------------------------

    def level_set_tangent(self, first, dirs, xyz):
        nblk = dirs.shape[2]
        jets = np.empty((40, xyz.shape[1], nblk))

        def one(j):
            jets[:, :, j] = self.lsf.jet_tangent(self.dm, self._dm_dir(first + j), dirs[:, :, j])

        _map_directions(one, nblk)
        return jets

    def field_tangent(self, first, dirs, xyz, d_xyz):
        nblk = dirs.shape[2]
        ngrid = xyz.shape[1]
        # Also where the density-contracted reductions are built: before the
        # loop, so no thread finds them missing and builds a second copy.
        efield = self.ints.field(self.dm)
        dphi = np.empty((ngrid, nblk))
        defield = np.empty((3, ngrid, nblk))

        def one(j):
            dm_dir = self._dm_dir(first + j)
            dphi[:, j] = self.ints.potential_tangent(self.dm, dm_dir, dirs[:, :, j], d_xyz[:, :, j])
            defield[:, :, j] = self.ints.field_tangent(self.dm, dm_dir, dirs[:, :, j], d_xyz[:, :, j])

        _map_directions(one, nblk)
        return efield, dphi, defield

    # -- the completion ------------------------------------------------------

    def complete(self, q, hvp, tangent):
        """Complete moist's products with the host's terms.

        Returns the nuclear columns ``(3, natm, n)`` and the Fock-matrix
        tangents ``(n, nao, nao)`` of every direction.
        """
        n = self.dirs.shape[2]
        nao = self.host.mol.nao
        nuclear = np.array(hvp, dtype=float, copy=True)
        fock = np.zeros((n, nao, nao))
        weights = None
        if self.lsf is not None:
            weights = {0: tangent.w_value, 1: tangent.w_gradient, 2: tangent.w_hessian}
        charged = self.ints is not None and q is not None
        if charged:
            # Same reason as in field_tangent: these are what every direction
            # shares, so they are built here rather than on first use inside
            # the loop.
            self.ints.charge_weighted(q)
            self.ints.density_contracted(self.dm)

        def one(j):
            v_r = self.dirs[:, :, j]
            dm_dir = self._dm_dir(j)
            d_xyz = tangent.xyz[:, :, j]
            if charged:
                dq = tangent.surface_charge[:, j]
                nuclear[:, :, j] += self.ints.gradient(dq, self.dm)
                if dm_dir is not None:
                    nuclear[:, :, j] += self.ints.gradient_density_tangent(q, dm_dir)
                nuclear[:, :, j] += self.ints.gradient_centre_tangent(q, self.dm, v_r)
                nuclear[:, :, j] += self.ints.gradient_point_tangent(q, self.dm, d_xyz)
                fock[j] += self.ints.fock(dq)
                fock[j] += self.ints.fock_centre_tangent(q, v_r)
                fock[j] += self.ints.fock_point_tangent(q, d_xyz)
            if weights is not None:
                dweights = {
                    0: tangent.dw_value[:, j],
                    1: tangent.dw_gradient[:, :, j],
                    2: tangent.dw_hessian[:, :, :, j],
                }
                moving = {
                    1: _symmetric_outer(weights[0], d_xyz),
                    2: _symmetric_outer(weights[1], d_xyz),
                    3: _symmetric_outer(weights[2], d_xyz),
                }
                nuclear[:, :, j] += self.lsf.gradient(dweights, self.dm)
                if dm_dir is not None:
                    nuclear[:, :, j] += self.lsf.gradient(weights, dm_dir)
                nuclear[:, :, j] += self.lsf.gradient_centre_tangent(weights, self.dm, v_r)
                nuclear[:, :, j] += self.lsf.gradient(moving, self.dm)
                fock[j] += self.lsf.fock(dweights)
                fock[j] += self.lsf.fock_centre_tangent(weights, v_r)
                fock[j] += self.lsf.fock(moving)

        _map_directions(one, n)
        return nuclear, fock


def solvent_hvp(evaluation, host, dirs, dm_dirs=None, *, shared=None):
    """Coupled Hessian-vector products of the solvent energy along mixed directions.

    ``dirs (3, natm, n)`` are the nuclear parts and ``dm_dirs (n, nao, nao)``
    the (symmetric) density parts of ``n`` directions. Returns the total
    nuclear columns ``(3, natm, n)`` and the Fock-matrix tangents ``(n, nao,
    nao)`` -- the response of the solvent Fock contribution along each
    direction, with the surface, the charges and the level-set weights
    relaxed and the host density held to its direction.

    ``shared`` is an optional dict, owned by the caller and spanning one solve,
    in which the geometry-only host tables are kept so that repeated calls at
    the same surface do not rebuild them; see :func:`_shared_host_terms`.
    """
    dirs = np.asfortranarray(dirs, dtype=float)
    model = evaluation._model
    density_defined = bool(model.cavity.density_dependent)
    has_charges = evaluation.response.electrostatics is not None
    coupling = evaluation._coupling
    coupling.activate()
    # Points as rows, the layout PySCF's grid evaluators and charge bases take
    coords = model.cavity.snapshot().xyz.T
    dm = coupling.density_matrix
    tangent_host = DirectionalHost(
        host, dm, coords, dirs, dm_dirs, level_set=density_defined,
        potential=has_charges, shared=shared,
    )
    hvp, tangent = evaluation.hvp_coupled(dirs, tangent_host, level_set=density_defined)
    q = evaluation.charges if has_charges else None
    return tangent_host.complete(q, hvp, tangent)


def _solvent_response(mean_field, model, host):
    """Request model-owned second derivatives in host parameter coordinates."""
    _check_restricted(mean_field)
    coupling = host.coupling(mean_field.make_rdm1())
    parameters = coupling.second_order()
    evaluation = model.evaluate(coupling=coupling)
    return evaluation.linearize(parameters, max_memory=mean_field.max_memory), parameters


def _check_restricted(mean_field):
    """Refuse a reference the solvent coupling cannot be attached to.

    The solvent terms themselves know nothing about the electronic method --
    they are second derivatives of the model energy in the nuclear coordinates
    and the AO density matrix. What has to hold is that a single closed-shell
    real density matrix is the whole host parameter, and that PySCF's own
    restricted Hessian is the object the terms are folded into. Restricted KS
    satisfies both, so it is accepted alongside RHF; everything that breaks the
    single-real-symmetric-P assumption is refused by name.
    """
    from pyscf.scf.hf import RHF

    if not mean_field.converged:
        raise ValueError("SCF must converge before evaluating its Hessian")
    if (not isinstance(mean_field, RHF)
            or mean_field.mol.has_ecp() or np.iscomplexobj(mean_field.mo_coeff)
            or np.any((mean_field.mo_occ != 0) & (mean_field.mo_occ != 2))
            or getattr(mean_field, "with_x2c", None) is not None):
        raise NotImplementedError(
            "The analytic coupling requires a conventional real, closed-shell, "
            "all-electron restricted reference (RHF or RKS)"
        )
    fitting = getattr(mean_field, "with_df", None)
    if fitting is not None and not _has_df_hessian(mean_field):
        raise NotImplementedError(
            f"{type(fitting).__name__} has no PySCF second-derivative "
            "implementation, so the electronic half of the Hessian cannot be "
            "formed consistently with the converged reference"
        )


#: Historical name; the check now admits restricted KS as well as RHF.
_check_rhf = _check_restricted


def _has_df_hessian(mean_field) -> bool:
    """Whether PySCF implements second derivatives for this fitted reference.

    Plain density fitting does; the semi-numerical exchange (COSX/SGX) fit does
    not, and differentiating it as if it were plain DF would silently return the
    Hessian of a different Fock operator than the one that was converged.
    """
    from pyscf.df.df_jk import _DFHF

    fitting = getattr(mean_field, "with_df", None)
    if fitting is None:
        return False
    if not isinstance(mean_field, _DFHF):
        return False
    return not type(fitting).__module__.startswith("pyscf.sgx")


def _pyscf_hessian_base(mean_field):
    """PySCF's own restricted Hessian class for this reference.

    Four classes, picked on two independent axes: Kohn-Sham or Hartree-Fock, and
    density-fitted or not.  All four take the same three overrides -- the
    density-fitted pair subclasses the plain one and replaces
    ``partial_hess_elec`` and ``make_h1`` with the same signatures, and every one
    of them inherits ``solve_mo1`` -- and the CPHF kernel fallback goes through
    ``mf.gen_response``, which already carries the XC kernel and the fitted
    Coulomb operator.  So selecting the base class is the whole of the
    generalization; the solvent terms are untouched by it, being derivatives of
    the model energy in the nuclear coordinates and the AO density matrix.
    """
    from pyscf.dft.rks import KohnShamDFT

    kohn_sham = isinstance(mean_field, KohnShamDFT)
    if _has_df_hessian(mean_field):
        if kohn_sham:
            from pyscf.df.hessian import rks as df_rks

            return df_rks.Hessian
        from pyscf.df.hessian import rhf as df_rhf

        return df_rhf.Hessian
    if kohn_sham:
        from pyscf.hessian import rks

        return rks.Hessian
    from pyscf.hessian import rhf

    return rhf.Hessian


def _directional_solvent_terms(mean_field, model, host):
    """The solvent's Hessian terms through the directional protocol.

    Returns the explicit nuclear block ``(nat, nat, 3, 3)``, the Fock tangents
    along the nuclear unit directions ``(nat, 3, nao, nao)`` and a function
    mapping density directions ``(nset, nao, nao)`` to their Fock tangents.
    """
    _check_restricted(mean_field)
    mol = mean_field.mol
    nat = mol.natm
    coupling = host.coupling(mean_field.make_rdm1())
    evaluation = model.evaluate(coupling=coupling)
    if not model.cavity.density_dependent and evaluation.response.electrostatics is None:
        raise NotImplementedError("The model exposes no second-order host channel")

    unit = np.zeros((3, nat, 3 * nat), order="F")
    for atom in range(nat):
        for axis in range(3):
            unit[axis, atom, 3 * atom + axis] = 1.0
    # One cache for the whole solve: the surface is fixed at the converged
    # density, so the nuclear block below and every coupled-perturbed iteration
    # afterwards ask for the host tables at the same points.
    shared = {}
    nuclear, fock = solvent_hvp(evaluation, host, unit, shared=shared)
    hrr = nuclear.reshape(3, nat, nat, 3).transpose(2, 1, 3, 0)
    hrp = fock.reshape(nat, 3, mol.nao, mol.nao)

    def density_fock_tangent(dm_dirs):
        dm_dirs = np.asarray(dm_dirs)
        zero = np.zeros((3, nat, len(dm_dirs)), order="F")
        _, tangent = solvent_hvp(evaluation, host, zero, dm_dirs, shared=shared)
        return tangent

    return hrr, hrp, density_fock_tangent


def rhf_hessian(mean_field, model, host, *, method="dense"):
    """Build a PySCF analytic Hessian using a model-owned solvent response.

    Includes moving AO centers, density-defined surface motion, PCM charge
    response, and solvent terms in both the CPHF right-hand side and kernel.
    Model and host must match the SCF settings; single-branch DROP projections
    and a conventional real all-electron restricted reference are supported.
    A restricted KS reference is accepted as well as RHF: the solvent terms are
    derivatives of the model energy in the nuclear coordinates and the AO
    density matrix, so they are the same either way, and PySCF's own restricted
    KS Hessian supplies the electronic half.

    ``method="dense"`` assembles the solvent second derivatives in the full
    host parameter space through :meth:`Evaluation.linearize`, a reference
    that stores surface partials for every pair of parameters. ``method=
    "directional"`` uses the single protocol shared with the Fortran and C
    layers: moist differentiates along directions and returns the response
    tangent, and the host completes the columns from it, so nothing quadratic
    in the number of density variables is ever stored and the coupled-perturbed
    kernel is applied direction by direction.

    The returned object carries a ``kernel_callback`` attribute, ``None`` by
    default, which a host may set to observe each application of the
    coupled-perturbed kernel; see the attribute's own documentation.
    """
    from pyscf.hessian import rhf

    base_hessian = _pyscf_hessian_base(mean_field)
    if method == "dense":
        solvent, derivatives = _solvent_response(mean_field, model, host)
        nat = mean_field.mol.natm
        hrr = solvent.nuclear.reshape(nat, 3, nat, 3).transpose(0, 2, 1, 3)
        hrp = derivatives.unpack_fock(solvent.mixed).reshape(nat, 3, *derivatives.dm.shape)

        def density_fock_tangent(dm1):
            packed = dm1[:, derivatives.rows, derivatives.cols]
            return derivatives.unpack_fock(solvent.density_response(packed))
    elif method == "directional":
        hrr, hrp, density_fock_tangent = _directional_solvent_terms(mean_field, model, host)
        nat = mean_field.mol.natm
    else:
        raise ValueError(f"unknown Hessian method {method!r}; use 'dense' or 'directional'")
    nr = 3 * nat

    class Hessian(base_hessian):
        #: Optional progress hook, called with the number of density directions
        #: each time the coupled-perturbed kernel is applied -- once per Krylov
        #: iteration per block of atoms. It lives here rather than in the host
        #: because the solvent kernel is what makes an iteration expensive: every
        #: application goes through ``solvent_hvp``, so an iteration count is the
        #: only honest measure of how much work the solve is doing. ``None``
        #: disables it, which is the default and costs one attribute read.
        kernel_callback = None

        def partial_hess_elec(self, mo_energy=None, mo_coeff=None, mo_occ=None,
                              atmlst=None, max_memory=4000, verbose=None):
            atoms = list(range(nat)) if atmlst is None else list(atmlst)
            return super().partial_hess_elec(
                mo_energy, mo_coeff, mo_occ, atoms, max_memory, verbose
            ) + hrr[atoms][:, atoms]

        def make_h1(self, mo_coeff, mo_occ, chkfile=None, atmlst=None, verbose=None):
            h1 = super().make_h1(mo_coeff, mo_occ, chkfile, atmlst, verbose)
            for atom in range(nat) if atmlst is None else atmlst:
                h1[atom] += hrp[atom]
            return h1

        def solve_mo1(self, mo_energy, mo_coeff, mo_occ, h1ao, fx=None,
                      atmlst=None, max_memory=4000, verbose=None):
            if fx is None:
                gas_response = rhf.gen_vind(self.base, mo_coeff, mo_occ)
                occupied = mo_coeff[:, mo_occ > 0]

                def fx(mo1):
                    rotations = mo1.reshape(-1, mo_coeff.shape[1], occupied.shape[1])
                    if self.kernel_callback is not None:
                        self.kernel_callback(len(rotations))
                    dm1 = 2 * _contract("um,pmi,vi->puv", mo_coeff, rotations, occupied)
                    dm1 += dm1.transpose(0, 2, 1).copy()
                    f1 = density_fock_tangent(dm1)
                    return gas_response(mo1) + _contract("um,puv,vi->pmi", mo_coeff, f1, occupied)

            return super().solve_mo1(mo_energy, mo_coeff, mo_occ, h1ao, fx, atmlst, max_memory, verbose)

        def gen_hop(self, mo_energy=None, mo_coeff=None, mo_occ=None, verbose=None):
            # PySCF's generic HVP constructs a gas-phase response kernel itself,
            # bypassing solve_mo1. Use the complete dense Hessian here as well.
            dense = self.kernel(mo_energy, mo_coeff, mo_occ, atmlst=range(nat))
            matrix = dense.transpose(0, 2, 1, 3).reshape(nr, nr)
            return matrix.dot, matrix.diagonal().copy()

    return Hessian(mean_field)


def solvated_rhf_gradient(mean_field, model, host) -> np.ndarray:
    """Return the total analytic RHF gradient, shape ``(natoms, 3)``.

    ``model`` and ``host`` must describe the same geometry and solvent settings
    used to converge ``mean_field``. The MOIST contribution is evaluated at its
    converged density, including the host's explicit nuclear derivatives.

    The gas-phase half comes from the reference's own gradient method, so a
    restricted KS reference contributes its exchange-correlation terms rather
    than silently being differentiated as Hartree-Fock.
    """
    if not mean_field.converged:
        raise ValueError("SCF must converge before evaluating its gradient")
    result = model.evaluate(coupling=host.coupling(mean_field.make_rdm1()))
    return mean_field.nuc_grad_method().kernel() + result.gradient.T


#: Central-difference stencils for a first derivative: offsets in units of the
#: step, their weights, and the common denominator. Order 2 is the textbook
#: pair; order 4 is the five-point rule without its zero-weight centre, whose
#: truncation falls as ``h**4`` rather than ``h**2``. Which one a test wants is
#: set by where the two error walls cross: raising the order lowers truncation
#: at a fixed step, but never moves the round-off floor, which is set by the
#: accuracy of the differenced gradient divided by the step.
_FD_STENCILS = {
    2: ((-1, 1), (-1.0, 1.0), 2.0),
    4: ((-2, -1, 1, 2), (1.0, -8.0, 8.0, -1.0), 12.0),
}


def finite_difference_gradient_tangent(
    gradient: Callable[[np.ndarray], np.ndarray],
    positions: np.ndarray,
    direction: np.ndarray,
    *,
    step: float = 1.0e-3,
    stencil: int = 4,
) -> np.ndarray:
    """Differentiate an analytic nuclear gradient along one direction.

    The directional twin of :func:`finite_difference_hessian`: the Hessian
    contracted with ``direction (natoms, 3)``, for the cost of one stencil
    rather than one per coordinate. ``positions`` is never modified.
    """
    positions = np.asarray(positions, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if direction.shape != positions.shape:
        raise ValueError("direction must have the shape of positions")
    offsets, weights, denominator = _fd_stencil(stencil, step)

    total = None
    for offset, weight in zip(offsets, weights):
        value = np.array(gradient(positions + offset * step * direction), dtype=float, copy=True)
        if value.shape != positions.shape:
            raise ValueError("gradient must return shape (natoms, 3)")
        total = weight * value if total is None else total + weight * value
    return total / (denominator * step)


def _fd_stencil(stencil: int, step: float):
    """Validate a stencil order and step, and return the rule."""
    if stencil not in _FD_STENCILS:
        raise ValueError(f"unknown stencil order {stencil!r}; use {sorted(_FD_STENCILS)}")
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")
    return _FD_STENCILS[stencil]


def finite_difference_hessian(
    gradient: Callable[[np.ndarray], np.ndarray],
    positions: np.ndarray,
    *,
    step: float = 1.0e-3,
    stencil: int = 2,
) -> np.ndarray:
    """Differentiate an analytic nuclear gradient with central differences.

    ``positions`` and the callback's result have shape ``(natoms, 3)`` in
    bohr and Hartree/bohr, respectively. The callback must reconverge the SCF
    at each supplied geometry to include electronic relaxation.

    ``stencil`` selects the central rule: 2 costs two gradients per coordinate
    and truncates at ``h**2``, 4 costs four and truncates at ``h**4``. Against a
    converged analytic gradient 4 is the one worth paying for -- the reference
    is otherwise truncation-limited long before it is noise-limited -- but the
    default stays at 2 so the cost of an existing caller does not change.

    Return ``H[A, B, a, b] = d gradient[A, a] / d positions[B, b]`` in
    Hartree/bohr**2. The result is not symmetrized, so numerical or gradient
    errors remain visible. The input positions are never modified.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must have shape (natoms, 3)")
    offsets, weights, denominator = _fd_stencil(stencil, step)

    natoms = len(positions)
    result = np.empty((natoms, natoms, 3, 3))
    for atom in range(natoms):
        for axis in range(3):
            total = None
            for offset, weight in zip(offsets, weights):
                moved = positions.copy()
                moved[atom, axis] += offset * step
                # Copy before the next call: a scanner may reuse its buffer.
                value = np.array(gradient(moved), dtype=float, copy=True)
                if value.shape != positions.shape:
                    raise ValueError("gradient must return shape (natoms, 3)")
                total = weight * value if total is None else total + weight * value
            result[:, atom, :, axis] = total / (denominator * step)
    return result
