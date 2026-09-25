"""PySCF integration for MOIST.

Import this module to register ``mf.MOIST(cavity=..., components=...)`` on
PySCF mean-field objects.  Everything the integration needs is in this one
file, top to bottom:

1. :class:`PySCFHost` -- the AO-basis integrals: the surface potential and its
   raw derivatives, the density callback an isodensity cavity evaluates, and
   the contractions of MOIST's response items into Fock and nuclear-gradient
   contributions.
2. :class:`GaussianMoments` -- the GOSTSHYP half of the host: the three-centre
   Gaussian moment integrals a ``gaussian_moments`` request asks for, and the
   contraction of the returned amplitudes.
3. :class:`PySCFSolvation` -- the driver.  It owns the model and one coupling
   and runs the same six steps a Fortran host runs: prepare a phase, walk the
   requests with ``for request in coupling``, answer their missing outputs on
   the cavity grid with ``coupling.answer(...)``, read the result and the
   response.
4. ``moist_for_scf`` -- the mean-field hook that plugs the driver into
   PySCF's SCF and nuclear-gradient machinery.

PCM uses Gaussian potentials; density adjoints already include the level-set
scale.  Units are atomic throughout.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Optional

import numpy as np

from .configuration import (
    CFC, DROP, ISwiG, Isodensity, SvdW, CavityConfiguration, _LSF_SETTINGS,
)
from .interface import (
    Coupling,
    DensityResponse,
    GaussianMomentRequest,
    GaussianPotentialRequest,
    GostshypAmplitudeResponse,
    PointPotentialRequest,
    PotentialAdjointResponse,
    Response,
    SolvationModel,
    SolvationModelComponent,
    Structure,
    _immutable_array,
)
from .parameters import ModelParameters, _resolve

__all__ = [
    "DROP", "ISwiG", "SvdW", "CFC", "Isodensity",
    "GPA_TO_AU", "PySCFHost", "GaussianMoments", "PySCFSolvation", "moist_for_scf",
]

#: 1 GPa in Hartree / bohr^3.
GPA_TO_AU = 1.0e9 * 5.29177210903e-11**3 / 4.3597447222071e-18


def _component_index(axes: tuple[int, ...]) -> int:
    """Index of a cartesian derivative in PySCF's ``GTOval_*_deriv`` output.

    ``eval_gto`` returns the derivative orders concatenated, each block ordered
    by descending ``lx`` then descending ``ly``: order 1 is ``x, y, z`` and
    order 2 is ``xx, xy, xz, yy, yz, zz``.  ``axes`` is a tuple of cartesian
    directions, e.g. ``(0, 2)`` for ``d^2/dx dz``; order does not matter.
    """
    order = len(axes)
    base = sum((m + 1) * (m + 2) // 2 for m in range(order))
    lx = sum(1 for a in axes if a == 0)
    ly = sum(1 for a in axes if a == 1)
    position = 0
    for ix in range(order, -1, -1):
        for iy in range(order - ix, -1, -1):
            if (ix, iy) == (lx, ly):
                return base + position
            position += 1
    raise ValueError(f"invalid derivative axes: {axes}")


# -----------------------------------------------------------------------------
# Gaussian-per-grid-point integrals (shared by the potential and the moments)
# -----------------------------------------------------------------------------

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

#: Mirrors ``overlap_floor`` in ``src/moist/model/component/gostshyp.f90``.
#: Nothing on the energy, Fock or gradient path reads this copy; it exists for
#: the pressure-independent diagnostics of :class:`GaussianMoments`, which
#: cannot recover the mask from amplitudes that vanish with the pressure.
_OVERLAP_FLOOR = 1.0e-9


def _fakemol_gaussians(coords: np.ndarray, exponents: np.ndarray, angl: int):
    """One coefficient-1 GTO shell of angular momentum ``angl`` per grid point."""
    from pyscf import gto

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
    """Three-centre one-electron integrals over a Gaussian-per-grid-point fakemol."""
    fakemol = _fakemol_gaussians(centers, omega, angl)
    nbas = mol.nbas
    shls_slice = (0, nbas, 0, nbas, nbas, nbas + fakemol.nbas)
    return (mol + fakemol).intor(intor, shls_slice=shls_slice)


# -----------------------------------------------------------------------------
# The host: AO-basis integrals and the density callback
# -----------------------------------------------------------------------------


class PySCFHost:
    """The AO-basis half of the coupling for one PySCF molecule.

    The density matrix is mutable state (:attr:`dm`) because the isodensity
    surface follows the density: it must be current *before* every cavity
    build, and the level-set callback reads it on every point evaluation.

    :param mol: PySCF molecule.  Coordinates are read in bohr.
    """

    def __init__(self, mol) -> None:
        self.mol = mol
        self.dm: Optional[np.ndarray] = None
        self._gto_prefix = "GTOval_cart_deriv" if mol.cart else "GTOval_sph_deriv"
        self._aoslice = mol.aoslice_by_atom()

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------

    def structure(self) -> Structure:
        """Molecular structure in moist's representation (bohr)."""
        numbers = self.mol.atom_charges()
        if self.mol.has_ecp():
            raise ValueError("effective core potentials are not supported")
        return Structure(numbers, self.mol.atom_coords())

    def _density_matrix(self) -> np.ndarray:
        if self.dm is None:
            raise RuntimeError("set host.dm before using the host")
        return np.asarray(self.dm)

    def _ao(self, coords: np.ndarray, order: int) -> np.ndarray:
        """AO values and derivatives up to ``order``: ``(ncomp, ngrid, nao)``."""
        return self.mol.eval_gto(f"{self._gto_prefix}{order}", np.asarray(coords))

    # ------------------------------------------------------------------
    # level set
    # ------------------------------------------------------------------

    def density(self, point: np.ndarray, order: int):
        """Density callback of an isodensity cavity built with ``source=host``.

        Returns the bare ``rho`` and its spatial derivatives up to ``order``
        only, so the projection's value+gradient phase never pays for the
        density Hessian.  moist owns the level set built from it: it subtracts
        ``rho_iso`` and applies ``scale`` and the DROP sign convention.  Doing
        any of that here as well would move the surface or square the scale.
        """
        dm = self._density_matrix()
        ao = self._ao(np.asarray(point).reshape(1, 3), order)[:, 0, :]
        # t[c1, c2] = sum_uv ao[c1, u] P_uv ao[c2, v]; symmetric in (c1, c2).
        #
        # BLAS leaves dirty floating-point status flags behind for these shapes
        # -- the SIMD tail reads padding lanes -- so numpy reports divide-by-zero
        # and overflow from a product that performs neither.  The results agree
        # with a BLAS-free einsum to rounding, which `test_density_callback_is_finite`
        # pins, so the flags are suppressed here rather than paying the
        # order-of-magnitude cost of the einsum path in the hottest callback.
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            proj = ao @ dm
            t = proj @ np.ascontiguousarray(ao.T)

        rho = t[0, 0]
        i1 = [_component_index((a,)) for a in range(3)]
        drho = 2.0 * t[i1, 0]
        if order < 2:
            return rho, drho

        i2 = [[_component_index((a, b)) for b in range(3)] for a in range(3)]
        d2rho = np.empty((3, 3))
        for a in range(3):
            for b in range(3):
                d2rho[a, b] = 2.0 * (t[i2[a][b], 0] + t[i1[a], i1[b]])
        if order < 3:
            return rho, drho, d2rho

        d3rho = np.empty((3, 3, 3))
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    d3rho[a, b, c] = 2.0 * (
                        t[_component_index((a, b, c)), 0]
                        + t[i2[a][b], i1[c]]
                        + t[i2[a][c], i1[b]]
                        + t[i2[b][c], i1[a]]
                    )
        return rho, drho, d2rho, d3rho

    # ------------------------------------------------------------------
    # electrostatics: what a potential request asks for
    # ------------------------------------------------------------------

    def _gaussian_integrals(self, coords, xi, intor="int3c2e"):
        """Coulomb integrals against normalized Gaussian surface charges."""
        from pyscf import df, gto

        auxiliary = gto.fakemol_for_charges(np.asarray(coords), expnt=np.asarray(xi)**2)
        auxiliary.cart = self.mol.cart
        return df.incore.aux_e2(
            self.mol, auxiliary, intor=self.mol._add_suffix(intor), aosym="s1"
        )

    def surface_potential(self, coords: np.ndarray, xi=None) -> np.ndarray:
        """Total molecular potential; ``xi`` selects the Gaussian charge operator."""
        from scipy.special import erf

        coords = np.asarray(coords)
        dm = self._density_matrix()
        delta = coords[:, None, :] - self.mol.atom_coords()[None, :, :]
        dist = np.linalg.norm(delta, axis=2)
        if xi is None:
            kernel = 1.0 / dist
            electronic = np.einsum("iuv,uv->i", self.mol.intor("int1e_grids", grids=coords), dm)
        else:
            width = np.asarray(xi)[:, None]
            kernel = np.divide(erf(width * dist), dist, out=np.broadcast_to(
                2.0 * width / np.sqrt(np.pi), dist.shape).copy(), where=dist != 0)
            electronic = np.einsum("uvi,uv->i", self._gaussian_integrals(coords, xi), dm)
        return np.einsum("iA,A->i", kernel, self.mol.atom_charges()) - electronic

    def _grad_phi_elec(self, coords: np.ndarray, xi=None) -> np.ndarray:
        """Electronic potential derivative with respect to surface centers, ``(3, ngrid)``."""
        dm = self._density_matrix()
        if xi is not None:
            return np.einsum("kuvi,uv->ki", self._gaussian_integrals(coords, xi, "int3c2e_ip2"), dm)
        tint = self.mol.intor("int1e_grids_ip", grids=coords)
        return -2.0 * np.einsum("kiuv,uv->ki", tint, dm)

    def _grad_phi_nuc(self, coords: np.ndarray, xi=None) -> np.ndarray:
        """Nuclear potential derivative with respect to surface centers, ``(3, ngrid)``."""
        from scipy.special import erf

        delta = np.asarray(coords)[:, None, :] - self.mol.atom_coords()[None, :, :]
        dist = np.linalg.norm(delta, axis=2)
        if xi is None:
            radial = dist**-3
        else:
            width = np.asarray(xi)[:, None]
            x = width * dist
            numerator = erf(x) - 2.0 * x * np.exp(-x*x) / np.sqrt(np.pi)
            radial = np.divide(numerator, dist**3, out=np.zeros_like(dist), where=dist != 0)
            series = 4.0 * width**3 / (3.0*np.sqrt(np.pi)) * (1.0 - 0.6*x*x + 3.0*x**4/14.0)
            radial = np.where(np.abs(x) < 1e-3, series, radial)
        return -np.einsum("A,iAk,iA->ki", self.mol.atom_charges(), delta, radial)

    def surface_potential_gradient(self, coords: np.ndarray, xi=None) -> np.ndarray:
        """``dphi/dr`` at the grid points, ``(ngrid, 3)``: the ``dphi_dr`` output."""
        return np.ascontiguousarray((self._grad_phi_nuc(coords, xi) + self._grad_phi_elec(coords, xi)).T)

    def _dphi_dxi(self, coords: np.ndarray, xi: np.ndarray) -> np.ndarray:
        """Raw derivative of the Gaussian potential with respect to its width."""
        delta = np.asarray(coords)[:, None, :] - self.mol.atom_coords()[None, :, :]
        nuclear = (2.0 / np.sqrt(np.pi)) * (np.exp(
            -np.asarray(xi)[:, None]**2 * np.sum(delta*delta, axis=2)) @ self.mol.atom_charges())
        dm = self._density_matrix()
        if not self.mol.cart:
            c2s = self.mol.cart2sph_coeff(normalized="sp")
            dm = c2s @ dm @ c2s.T
        # The unnormalized s probe carries 1/(2 sqrt(pi)); the kernel
        # derivative carries 2/sqrt(pi), hence the factor four.
        electronic = 4.0 * np.einsum("uvi,uv->i", _int3c1e(
            self.mol, coords, np.asarray(xi)**2, 0), dm)
        return nuclear - electronic

    # ------------------------------------------------------------------
    # contractions of the response items
    # ------------------------------------------------------------------

    def fock_potential(self, coords: np.ndarray, xi: np.ndarray, w_phi: np.ndarray) -> np.ndarray:
        """``sum_i w_phi_i dphi_i/dP``: the potential adjoint contracted with the Gaussian integrals."""
        return -np.einsum("i,uvi->uv", np.asarray(w_phi), self._gaussian_integrals(coords, xi))

    def _fock_lsf(self, coords: np.ndarray, density: DensityResponse) -> np.ndarray:
        """``d rho/dP`` contracted with the density adjoints.

        With ``rho = sum P_uv chi_u chi_v`` every term is a weighted outer
        product of AO derivative blocks summed over the grid.  Nothing here
        carries ``scale``: moist already folded ``dS/drho`` into the weights,
        so these are adjoints of the density, not of the level set.
        """
        ao = self._ao(coords, 2)
        w0 = np.asarray(density.w_rho)
        w1 = np.asarray(density.w_grad_rho).T
        w2 = np.asarray(density.w_hess_rho).T

        def outer(weight, left, right):
            return np.einsum("i,iu,iv->uv", weight, left, right)

        def symmetrised(weight, left, right):
            block = outer(weight, left, right)
            return block + block.T

        a0 = ao[0]
        fock = outer(w0, a0, a0)
        for a in range(3):
            fock += symmetrised(w1[a], ao[_component_index((a,))], a0)
        for a in range(3):
            for b in range(3):
                fock += symmetrised(w2[a, b], ao[_component_index((a, b))], a0)
                fock += symmetrised(
                    w2[a, b], ao[_component_index((a,))], ao[_component_index((b,))]
                )
        return fock

    def _gradient_phi(self, coords: np.ndarray, q: np.ndarray, xi=None) -> np.ndarray:
        """``sum_i q_i d(phi_elec)_i/dR_A`` at fixed grid points and fixed P, ``(natm, 3)``."""
        dm = self._density_matrix()

        # d(r_i|uv)/dR_A = -T[k,i,u,v] delta_{u in A} - T[k,i,v,u] delta_{v in A},
        # and phi_elec carries a further minus sign; P symmetry merges the two
        # halves into a single factor of two.
        if xi is None:
            tint = self.mol.intor("int1e_grids_ip", grids=coords)
            weighted = 2.0 * np.einsum("i,kiuv,uv->ku", q, tint, dm)
        else:
            tint = self._gaussian_integrals(coords, xi, "int3c2e_ip1")
            weighted = 2.0 * np.einsum("i,kuvi,uv->ku", q, tint, dm)
        gradient = np.zeros((3, self.mol.natm))
        for atom in range(self.mol.natm):
            lo, hi = self._aoslice[atom, 2:4]
            gradient[:, atom] = weighted[:, lo:hi].sum(axis=1)
        return np.ascontiguousarray(gradient.T)

    def _gradient_lsf(self, coords: np.ndarray, density: DensityResponse) -> np.ndarray:
        """Density adjoints contracted with ``d rho/dR_A`` at fixed P, ``(natm, 3)``.

        Because ``d/dR_A`` commutes with the spatial derivatives, every order is
        a spatial derivative of the single "displaced density"

        .. math:: u^{Ak}(r) = \\sum_{\\mu \\in A, \\nu} P_{\\mu\\nu}
                  \\partial_k \\chi_\\mu(r) \\chi_\\nu(r)

        with ``drho^{(m)}/dR_Ak = -2 d^m u^{Ak}``.  The factor is exactly that:
        ``scale`` belongs to moist's level set and is already folded into the
        weights.
        """
        dm = self._density_matrix()
        ao = self._ao(coords, 3)
        w0 = np.asarray(density.w_rho)
        w1 = np.asarray(density.w_grad_rho).T
        w2 = np.asarray(density.w_hess_rho).T

        # right[c, i, u] = sum_v ao[c, i, v] P_uv
        right = np.einsum("civ,uv->ciu", ao, dm)

        gradient = np.zeros((3, self.mol.natm))
        for k in range(3):
            # Leibniz expansion of d^m (d_k chi_u . chi_v) for every weighted
            # order m, accumulated per AO so it can be sliced by atom.
            terms = [(w0, (k,), ())]
            for a in range(3):
                terms.append((w1[a], (k, a), ()))
                terms.append((w1[a], (k,), (a,)))
            for a in range(3):
                for b in range(3):
                    terms.append((w2[a, b], (k, a, b), ()))
                    terms.append((w2[a, b], (k, a), (b,)))
                    terms.append((w2[a, b], (k, b), (a,)))
                    terms.append((w2[a, b], (k,), (a, b)))

            per_ao = np.zeros(self.mol.nao)
            for weight, left_axes, right_axes in terms:
                left = ao[_component_index(left_axes)]
                rgt = right[_component_index(right_axes)]
                per_ao += np.einsum("i,iu,iu->u", weight, left, rgt)

            per_ao *= -2.0
            for atom in range(self.mol.natm):
                lo, hi = self._aoslice[atom, 2:4]
                gradient[k, atom] = per_ao[lo:hi].sum()
        return np.ascontiguousarray(gradient.T)


# -----------------------------------------------------------------------------
# GOSTSHYP: the Gaussian moments a pressure component asks for
# -----------------------------------------------------------------------------


class GaussianMoments:
    """The host half of GOSTSHYP: Gaussian moment integrals on the cavity grid.

    GOSTSHYP places an unnormalized Gaussian ``G_j = exp(-omega_j |r-C_j|^2)``
    on every grid point.  The component asks the host for the density moments
    ``gt = <G>``, ``pt = <(r-C) G>``, ``mt = <(r-C)(r-C) G>`` and
    ``rt = <(r-C)|r-C|^2 G>`` (the ``gaussian_moments`` request) and hands back
    the amplitudes ``w_overlap`` and ``w_normal_deriv`` the host contracts as
    ``F += sum_j [w_overlap_j g_j + w_normal_deriv_j f_j]`` with
    ``f_j = n_j . grad_r g_j``.  ``ftilde`` is not exchanged: moist derives it
    as ``-2 omega_j (n_j . pt_j)`` from the same moments it differentiates.

    Every per-grid-point Gaussian carries a normalization that cancels between
    energy and Fock, so it is not formed.  Only the relative s/p/d/f angular
    constants are restored, which is what makes ``f = n . grad g`` exact; they
    are pinned against an independent quadrature in ``test_gostshyp.py``.

    The widths ``omega`` are the component's choice, read from the request
    snapshot; the pressure is the component's too.  Bound to one surface and
    one density by :meth:`bind`; the driver rebinds on every evaluation.
    """

    def __init__(self, host: PySCFHost) -> None:
        self.host = host
        self.mol = host.mol
        self.ngrid = 0
        self._dm: Optional[np.ndarray] = None
        self._traces: Optional[tuple[np.ndarray, np.ndarray]] = None
        self._c2s: Optional[np.ndarray] = None
        self._G: Optional[np.ndarray] = None
        self._F: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # surface and density
    # ------------------------------------------------------------------

    def bind(self, cavity, dm: np.ndarray) -> None:
        """Snapshot the live cavity grid and fix the density the moments are of."""
        self.centers = np.ascontiguousarray(cavity.xyz, dtype=np.float64)
        self.areas = np.ascontiguousarray(cavity.a, dtype=np.float64)
        self.ngrid = int(self.centers.shape[0])

        normals = np.array(cavity.normal0, dtype=np.float64, order="C", copy=True)
        norm = np.linalg.norm(normals, axis=1)
        good = norm > 0.0
        normals[good] /= norm[good, None]
        self.normals = normals

        # omega_j = pi ln2 / a_j, matching `gaussian_width` in the component; a
        # degenerate zero-area grid point is inert.  The request's own widths
        # replace these when a moment request is answered.
        with np.errstate(divide="ignore", invalid="ignore"):
            omega = np.pi * math.log(2.0) / self.areas
        omega[~np.isfinite(omega)] = 0.0
        self.omega = omega

        self._dm = dm
        self._traces = None
        self._G = self._F = None

    def answer(self, coupling: Coupling, request: GaussianMomentRequest) -> None:
        """Answer the missing moments of the current request from the bound surface.

        ``request`` is the snapshot the coupling's loop yielded for its current
        request: it supplies the widths and the moments that are missing.
        """
        self.omega = np.asarray(request.width)
        gt, pt, mt, rt = self._surface_moments(
            self._density_matrix_cart(self._dm), required=request.missing
        )
        coupling.answer(gt=gt, pt=pt, mt=mt, rt=rt)
        if gt is not None and pt is not None:
            self._traces = (
                gt,
                -2.0 * self.omega * np.einsum("ja,ja->j", self.normals, pt, optimize=True),
            )

    # ------------------------------------------------------------------
    # integrals
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

    def _density_matrix_cart(self, dm: np.ndarray) -> np.ndarray:
        """Density matrix in the cartesian AO basis the fakemol blocks use."""
        c2s = self._cart2sph
        # BLAS leaves dirty FP status flags for these shapes; see PySCFHost.density.
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            return c2s @ np.asarray(dm) @ c2s.T

    def f_vector(self) -> np.ndarray:
        """``grad_r g_uv,j`` before the normal projection, ``(nao, nao, ngrid, 3)``.

        Recomputed rather than cached: only the projected ``f`` is needed to
        build a Fock matrix, and this block is three times its size.  Kept as a
        method because it is what pins the p-shell angular constant.
        """
        p_cart = _int3c1e(self.mol, self.centers, self.omega, 1)
        ncart = p_cart.shape[0]
        p_cart = p_cart.reshape(ncart, ncart, self.ngrid, 3)
        # dG/dC_a = 2 omega (r_a - C_a) G, and displacing the field point is the
        # opposite of displacing the center, so grad_r g = -2 omega <(r-C) G>.
        fvec_cart = -_S_OVER_P_NORM * 2.0 * np.einsum(
            "j,pqja->pqja", self.omega, p_cart, optimize=True
        )
        return self._to_spherical(fvec_cart)

    def _build_integrals(self) -> None:
        """The dense ``g`` and ``f`` blocks the amplitudes are contracted with."""
        g_cart = _int3c1e(self.mol, self.centers, self.omega, 0)
        self._G = self._to_spherical(g_cart)
        self._F = np.einsum("uvja,ja->uvj", self.f_vector(), self.normals, optimize=True)

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

    def _surface_moments(self, dm_cart: np.ndarray, centers=None, omega=None, *, required=None):
        """``(gt, pt, mt, rt)`` -- the s/p/d/f Gaussian moments of the density.

        Numpy-ordered: ``(ngrid,)``, ``(ngrid, 3)``, ``(ngrid, 3, 3)`` and
        ``(ngrid, 3)``, the shapes the request takes.  Only the ``required``
        moments are built.  ``centers``/``omega`` default to the bound surface
        and are overridable so a finite difference can rebuild the moments at
        a displaced surface without touching cached state.
        """
        centers = self.centers if centers is None else centers
        omega = self.omega if omega is None else omega
        ngrid = int(np.asarray(omega).size)

        required = {"gt", "pt", "mt", "rt"} if required is None else required
        ncart = dm_cart.shape[0]
        gt = pt = mt = rt = None
        if "gt" in required:
            block = _int3c1e(self.mol, centers, omega, 0)
            gt = np.einsum("pqj,pq->j", block, dm_cart, optimize=True)
        if "pt" in required:
            block = _int3c1e(self.mol, centers, omega, 1)
            pt = self._contract_p_moments(dm_cart, block.reshape(ncart, ncart, ngrid, 3))
        if "mt" in required:
            block = _int3c1e(self.mol, centers, omega, 2)
            mt = self._contract_d_moments(dm_cart, block.reshape(ncart, ncart, ngrid, 6))
        if "rt" in required:
            block = _int3c1e(self.mol, centers, omega, 3)
            rt = self._contract_f_rho2_moments(dm_cart, block.reshape(ncart, ncart, ngrid, 10))
        return gt, pt, mt, rt

    # ------------------------------------------------------------------
    # diagnostics
    # ------------------------------------------------------------------

    def traces(self, dm: np.ndarray, *, centers=None, omega=None, normals=None):
        """``(gtilde, ftilde)`` from this module's own moments.

        The host's copy of the two traces moist works from -- a diagnostic and
        a self-check, never an input to the energy.  The surface parameters
        default to the bound surface and are overridable for finite differences.
        """
        omega = self.omega if omega is None else omega
        normals = self.normals if normals is None else normals
        gt, pt, _mt, _rt = self._surface_moments(
            self._density_matrix_cart(dm), centers, omega, required={"gt", "pt"}
        )
        ftilde = -2.0 * np.asarray(omega) * np.einsum("ja,ja->j", normals, pt, optimize=True)
        return gt, ftilde

    @property
    def live_traces(self) -> tuple[np.ndarray, np.ndarray]:
        """``(gtilde, ftilde)`` on the bound surface, kept from the last answer."""
        if self._dm is None:
            raise RuntimeError("evaluate the model before reading GOSTSHYP traces")
        if self._traces is None:
            self._traces = self.traces(self._dm)
        return self._traces

    @property
    def inactive_count(self) -> int:
        """Grid points the component switched off, out of :attr:`ngrid`.

        Derived from ``ftilde`` rather than from the amplitudes: the amplitudes
        carry the pressure, so at ``p_inp = 0`` they report every point dropped.
        """
        _, ftilde = self.live_traces
        floor = _OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
        return int(np.count_nonzero(np.abs(ftilde) <= floor))

    def effective_volume(self) -> float:
        """``E / p_inp`` (eq 11), evaluated as ``sum_j a_j gtilde_j / ftilde_j``.

        Not the cavity volume, and well defined at ``p_inp = 0`` where the
        energy vanishes with the pressure but the volume does not.
        """
        gt, ftilde = self.live_traces
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = self.areas * gt / ftilde
        floor = _OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
        active = np.abs(ftilde) > floor
        return float(np.sum(np.where(active, ratio, 0.0)))

    # ------------------------------------------------------------------
    # contractions of the amplitudes
    # ------------------------------------------------------------------

    def fock(self, amplitude: GostshypAmplitudeResponse) -> np.ndarray:
        """``sum_j [w_overlap_j g_j + w_normal_deriv_j f_j]`` at a frozen surface."""
        if self._G is None:
            self._build_integrals()
        fock = np.einsum("j,uvj->uv", amplitude.w_overlap, self._G, optimize=True)
        fock += np.einsum("j,uvj->uv", amplitude.w_normal_deriv, self._F, optimize=True)
        return 0.5 * (fock + fock.T)

    def nuclear_gradient(self, dm: np.ndarray, amplitude: GostshypAmplitudeResponse) -> np.ndarray:
        """AO centers move, surface frozen: the ``int3c1e_ip1`` route, ``(natm, 3)``."""
        dm_cart = self._density_matrix_cart(dm)
        ncart = self._cart2sph.shape[0]

        ip1_g = _int3c1e(self.mol, self.centers, self.omega, 0, "int3c1e_ip1_cart")
        ip1_p = _int3c1e(self.mol, self.centers, self.omega, 1, "int3c1e_ip1_cart")
        ip1_g = ip1_g.reshape(3, ncart, ncart, self.ngrid)
        ip1_p = ip1_p.reshape(3, ncart, ncart, self.ngrid, 3)
        ip1_f = -_S_OVER_P_NORM * 2.0 * np.einsum(
            "j,xpqja,ja->xpqj", self.omega, ip1_p, self.normals, optimize=True
        )

        kernel = np.einsum("j,xpqj->xpq", amplitude.w_overlap, ip1_g, optimize=True)
        kernel += np.einsum("j,xpqj->xpq", amplitude.w_normal_deriv, ip1_f, optimize=True)
        row = np.einsum("xpq,pq->xp", kernel, dm_cart, optimize=True)

        # Cartesian AO rows fold onto atoms through the cartesian slices.
        gradient = np.zeros((3, self.mol.natm))
        aoslice = self.mol.aoslice_by_atom(self.mol.ao_loc_nr(cart=True))
        for atom in range(self.mol.natm):
            lo, hi = aoslice[atom, 2:4]
            # dg/dR = -ip1, and both AO legs contribute equally for symmetric P.
            gradient[:, atom] = -2.0 * row[:, lo:hi].sum(axis=1)
        return np.ascontiguousarray(gradient.T)


# -----------------------------------------------------------------------------
# The driver: one model, one coupling, the six-step protocol
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Result:
    """Energy and Fock contribution of one evaluation at a fixed density."""

    energy: float
    fock: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "fock", _immutable_array(np.asarray(self.fock)))


class PySCFSolvation:
    """Drive a MOIST model for a PySCF molecule.

    Owns the model and one coupling and runs the Fortran host loop for every
    phase: ``prepare_*``, then ``for request in coupling`` answering what is
    missing on the cavity grid, then ``get_*``.
    The requests are answered by :class:`PySCFHost` (potentials) and
    :class:`GaussianMoments` (GOSTSHYP moments); the response items are
    contracted by the same two into the Fock matrix and the nuclear gradient.

    :param mol: PySCF molecule.
    :param cavity: a cavity configuration, ``DROP(lsf=...)`` or ``ISwiG(...)``.
        An isodensity level set is bound to this molecule's density.
    :param components: the model components, e.g. ``[ModelComponentCPCM(80.0)]``.
    :param parameters: model logging settings.
    """

    def __init__(self, mol, cavity: CavityConfiguration, components, *,
                 parameters: ModelParameters | None = None) -> None:
        if not isinstance(cavity, CavityConfiguration):
            raise TypeError("cavity must be DROP(lsf=...) or ISwiG(...)")
        items = tuple(components)
        if not items or any(not isinstance(item, SolvationModelComponent) for item in items):
            raise TypeError("components must be a nonempty sequence of MOIST components")
        self.mol = mol
        self.configuration = cavity
        self.components = items
        self.parameters = _resolve(ModelParameters, parameters, {})
        self.host = PySCFHost(mol)
        self.host.structure()  # Validate the molecular representation before use.
        self.model = SolvationModel(
            cavity.build(source=self.host if cavity.density_dependent else None),
            items, parameters=self.parameters,
        )
        self.coupling = None
        #: The GOSTSHYP half, created the first time a moment request appears
        self.moments: Optional[GaussianMoments] = None
        self._dm: Optional[np.ndarray] = None
        self._result: Optional[Result] = None
        self._response: Optional[Response] = None
        self._xi: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # results of the last evaluation
    # ------------------------------------------------------------------

    def _require_result(self) -> Result:
        if self._result is None:
            raise RuntimeError("evaluate the model at a density first")
        return self._result

    @property
    def result(self) -> Optional[Result]:
        """The last :meth:`evaluate` result, or ``None`` before the first."""
        return self._result

    @property
    def energy(self) -> float:
        """Solvation energy of the last evaluation."""
        return self._require_result().energy

    @property
    def fock(self) -> np.ndarray:
        """Fock contribution ``dE/dP`` of the last evaluation, ``(nao, nao)``."""
        return self._require_result().fock

    @property
    def response(self) -> Response:
        """Every item the model handed back in the response phase."""
        self._require_result()
        return self._response

    @property
    def density_matrix(self) -> Optional[np.ndarray]:
        """The density the last evaluation was made at."""
        return self._dm

    # ------------------------------------------------------------------
    # the protocol
    # ------------------------------------------------------------------

    def _answer(self, request) -> None:
        """Answer the missing outputs of the current request on the cavity grid."""
        coupling = self.coupling
        cavity = self.model.cavity
        coords = cavity.xyz
        if isinstance(request, PointPotentialRequest):
            outputs = {}
            if "phi" in request.missing:
                outputs["phi"] = self.host.surface_potential(coords)
            if "dphi_dr" in request.missing:
                outputs["dphi_dr"] = self.host.surface_potential_gradient(coords)
            coupling.answer(**outputs)
        elif isinstance(request, GaussianPotentialRequest):
            xi = self._xi = cavity.xi0
            outputs = {}
            if "phi" in request.missing:
                outputs["phi"] = self.host.surface_potential(coords, xi)
            if "dphi_dr" in request.missing:
                outputs["dphi_dr"] = self.host.surface_potential_gradient(coords, xi)
            if "dphi_dxi" in request.missing:
                outputs["dphi_dxi"] = self.host._dphi_dxi(coords, xi)
            coupling.answer(**outputs)
        elif isinstance(request, GaussianMomentRequest):
            self._bound_moments().answer(coupling, request)
        else:
            raise NotImplementedError(f"PySCF host cannot answer {request.name!r}")

    def _bound_moments(self) -> GaussianMoments:
        """The GOSTSHYP half, bound to the current surface and density.

        Created on first demand -- a model without a moment request never makes
        the host form the three-centre integrals -- and rebound whenever the
        density changes.  A zero-pressure component asks for no moments but
        still hands back (zero) amplitudes, so the response contraction binds
        it as well; its diagnostics then describe the surface, not the pressure.
        """
        if self.moments is None:
            self.moments = GaussianMoments(self.host)
        if self.moments._dm is not self._dm:
            self.moments.bind(self.model.cavity, self._dm)
        return self.moments

    def _answer_requests(self) -> None:
        """One pass over the staged phase, answering every request that misses an output."""
        for request in self.coupling:
            self._answer(request)

    def _density(self, dm) -> np.ndarray:
        density = np.asarray(dm)
        n = self.mol.nao_nr()
        if density.shape == (2, n, n):
            density = density.sum(axis=0)
        if density.shape != (n, n):
            raise ValueError(f"Expected density shape {(n, n)} or {(2, n, n)}")
        if np.iscomplexobj(density):
            raise NotImplementedError("MOIST supports real density matrices only")
        if not np.all(np.isfinite(density)):
            raise ValueError("Density must contain finite values")
        return _immutable_array(np.array(density, dtype=np.float64))

    def evaluate(self, dm) -> Result:
        """Energy and Fock contribution at the density ``dm``.

        Rebuilds a density-dependent cavity, then runs the energy and the
        response phase.  The result is cached until the density changes.
        """
        density = self._density(dm)
        if self._result is not None and np.array_equal(density, self._dm):
            return self._result
        self._result = self._response = None
        self._dm = density
        self.host.dm = density
        try:
            if self.coupling is None or self.configuration.density_dependent:
                self.model.update(self.host.structure())
            if self.coupling is None:
                self.coupling = self.model.new_coupling()
            cpl = self.coupling

            self.model.prepare_energy(cpl)
            self._answer_requests()
            energy = np.array(0.0)
            self.model.get_energy(cpl, energy)

            self.model.prepare_response(cpl)
            self._answer_requests()
            response = self.model.get_response(cpl)
        except Exception:
            self._dm = None
            raise
        self._response = response
        self._result = Result(float(energy), self._fock(response))
        return self._result

    def _fock(self, response: Response, *, with_density: bool = True) -> np.ndarray:
        """Contract the response items into the Fock matrix."""
        nao = self.mol.nao
        coords = self.model.cavity.xyz
        fock = np.zeros((nao, nao))
        for item in response:
            if isinstance(item, PotentialAdjointResponse):
                fock += self.host.fock_potential(coords, self._xi, item.w_phi)
            elif isinstance(item, GostshypAmplitudeResponse):
                fock += self._bound_moments().fock(item)
            elif isinstance(item, DensityResponse):
                # Present exactly when the surface follows the density; the
                # frozen-surface Fock matrix leaves this route out on purpose.
                if with_density:
                    fock += self.host._fock_lsf(coords, item)
            else:
                raise NotImplementedError(f"Unsupported response item '{item.name}'")
        return fock

    def frozen_fock(self) -> np.ndarray:
        """The Fock contribution with the surface held fixed.

        Omits the density-weight contraction, so on a density-dependent cavity
        this is *not* ``dE/dP``: the surface moves with the density.
        """
        return self._fock(self.response, with_density=False)

    def gradient_channels(self, dm) -> dict[str, np.ndarray]:
        """The nuclear gradient at fixed ``dm``, split by route, each ``(natm, 3)``.

        ``model``
            Cavity motion and every term moist owns, from the gradient phase.
        ``potential``
            The basis-centre derivative of the potential, contracted with the
            potential adjoint.
        ``moments``
            The basis-centre derivative of the Gaussian moments, contracted
            with the GOSTSHYP amplitudes (present only with a moment request).
        ``density``
            The basis-centre derivative of the level set, contracted with the
            density weights (present only on a density-dependent cavity).
        """
        self.evaluate(dm)
        self.host.dm = self._dm
        cpl = self.coupling
        self.model.prepare_gradient(cpl)
        self._answer_requests()
        model = np.zeros((self.mol.natm, 3))
        self.model.get_gradient(cpl, model)

        # The response phase's items, not the ones get_gradient just returned.
        # Both carry the same potential adjoint and amplitudes, but only the
        # response phase emits `density`, and its weights (dE/drho at fixed
        # nuclei) are what the basis-centre route below contracts.  evaluate()
        # above ran that phase, so self.response holds them.
        response = self.response
        coords = self.model.cavity.xyz
        channels = {"model": model}
        for item in response:
            if isinstance(item, PotentialAdjointResponse):
                channels["potential"] = self.host._gradient_phi(
                    coords, np.asarray(item.w_phi), xi=self._xi)
            elif isinstance(item, GostshypAmplitudeResponse):
                channels["moments"] = self._bound_moments().nuclear_gradient(self._dm, item)
            elif isinstance(item, DensityResponse):
                channels["density"] = self.host._gradient_lsf(coords, item)
            else:
                raise NotImplementedError(f"Unsupported response item '{item.name}'")
        return channels

    def gradient(self, dm) -> np.ndarray:
        """Total nuclear gradient of the solvation energy at fixed ``dm``, ``(natm, 3)``."""
        return sum(self.gradient_channels(dm).values())


# -----------------------------------------------------------------------------
# The mean-field hook
# -----------------------------------------------------------------------------


class _MoistState:
    """One calculation's configuration and its driver.

    Change configuration with ``set``; this discards native state and results.
    ``result`` is the latest :class:`Result`, or ``None`` before use.
    """

    def __init__(self, mol, *, cavity, components, parameters=None):
        self.mol = mol
        self._cavity = None
        self._components = ()
        self._parameters = _resolve(ModelParameters, parameters, {})
        self.set(cavity=cavity, components=components)

    @property
    def cavity(self):
        return self._cavity

    @property
    def components(self):
        return self._components

    @property
    def parameters(self):
        return self._parameters

    @property
    def cavity_options(self):
        return self._cavity.options

    @property
    def solvation(self) -> Optional[PySCFSolvation]:
        """The driver of the current molecule, or ``None`` before the first use."""
        return self._solvation

    @property
    def result(self) -> Optional[Result]:
        return None if self._solvation is None else self._solvation.result

    @property
    def e(self):
        result = self.result
        return None if result is None else result.energy

    @property
    def v(self):
        result = self.result
        return None if result is None else result.fock

    def set(self, *, cavity=None, components=None, parameters=None, **options):
        """Replace model settings or cavity configuration, clearing cached results.

        Plain keywords are DROP/ISwiG parameters; ``lsf=`` replaces the level
        set.  Prefer ``set(cavity=replace(config, parameters=...))``.
        """
        config = self._cavity if cavity is None else cavity
        if not isinstance(config, CavityConfiguration):
            raise TypeError("cavity must be DROP(lsf=...) or ISwiG(...)")
        items = self._components if components is None else tuple(components)
        if not items or any(not isinstance(item, SolvationModelComponent) for item in items):
            raise TypeError("components must be a nonempty sequence of MOIST components")
        changes = {}
        if "lsf" in options:
            if not isinstance(config, DROP):
                raise TypeError("ISwiG does not accept an LSF")
            changes["lsf"] = options.pop("lsf")
        if options:
            if options.keys() & _LSF_SETTINGS:
                raise TypeError("Configure surface settings on the LSF")
            changes["parameters"] = replace(config.parameters, **options)
        config = replace(config, **changes) if changes else config
        model_parameters = self.parameters if parameters is None else _resolve(ModelParameters, parameters, {})
        self._cavity, self._components, self._parameters = config, items, model_parameters
        return self.reset()

    def reset(self, mol=None):
        if mol is not None:
            self.mol = mol
        self._fingerprint = None
        self._solvation = None
        return self

    def copy(self):
        return type(self)(self.mol, cavity=self.cavity, components=self.components,
                          parameters=self.parameters)

    def _molecule_key(self):
        mol = self.mol
        return (mol._atm.tobytes(), mol._bas.tobytes(), mol._env.tobytes(),
                mol.cart, mol.charge, mol.spin)

    def _driver(self) -> PySCFSolvation:
        key = self._molecule_key()
        if self._fingerprint != key:
            self.reset()
        if self._solvation is None:
            if self.mol.has_ecp():
                raise ValueError("effective core potentials are not supported")
            self._solvation = PySCFSolvation(
                self.mol, self.cavity, self.components, parameters=self.parameters)
            self._fingerprint = key
        return self._solvation

    def evaluate(self, dm) -> Result:
        try:
            return self._driver().evaluate(dm)
        except Exception:
            self.reset()
            raise

    def kernel(self, dm):
        result = self.evaluate(dm)
        return result.energy, result.fock

    def gradient(self, dm):
        try:
            return self._driver().gradient(dm)
        except Exception:
            self.reset()
            raise


class _MoistSCF:
    _keys = {"with_moist"}

    def undo_moist(self):
        from pyscf import lib

        obj = lib.view(self, lib.drop_class(self.__class__, _MoistSCF, "MOIST"))
        del obj.with_moist
        return obj

    def copy(self):
        obj = super().copy()
        obj.with_moist = self.with_moist.copy()
        return obj

    def as_scanner(self):
        from pyscf import lib, scf

        if isinstance(self, lib.SinglePointScanner):
            return self
        return scf.hf.as_scanner(self.copy())

    def reset(self, mol=None):
        self.with_moist.reset(mol)
        return super().reset(mol)

    def dump_flags(self, verbose=None):
        from pyscf.lib import logger

        super().dump_flags(verbose)
        logger.info(self, "MOIST cavity: %s", self.with_moist.cavity)
        return self

    def get_veff(self, mol=None, dm=None, *args, **kwargs):
        from pyscf import lib

        if dm is None:
            dm = self.make_rdm1()
        mol = self.mol if mol is None else mol
        if mol is not self.with_moist.mol:
            self.with_moist.reset(mol)
        veff = super().get_veff(mol, dm, *args, **kwargs)
        energy, potential = self.with_moist.kernel(dm)
        return lib.tag_array(veff, e_moist=energy, v_moist=potential)

    def get_fock(self, h1e=None, s1e=None, vhf=None, dm=None, *args, **kwargs):
        if dm is None:
            dm = self.make_rdm1()
        if getattr(vhf, "v_moist", None) is None:
            vhf = self.get_veff(self.mol, dm)
        return super().get_fock(h1e, s1e, vhf + vhf.v_moist, dm, *args, **kwargs)

    def energy_elec(self, dm=None, h1e=None, vhf=None):
        if dm is None:
            dm = self.make_rdm1()
        if getattr(vhf, "e_moist", None) is None:
            vhf = self.get_veff(self.mol, dm)
        energy, coulomb = super().energy_elec(dm, h1e, vhf)
        self.scf_summary["e_moist"] = float(vhf.e_moist)
        return energy + vhf.e_moist, coulomb

    def get_grad(self, mo_coeff, mo_occ, fock=None):
        if fock is None:
            fock = self.get_fock(dm=self.make_rdm1(mo_coeff, mo_occ))
        return super().get_grad(mo_coeff, mo_occ, fock)

    def nuc_grad_method(self):
        from pyscf import lib

        grad = self.undo_moist().nuc_grad_method()
        grad.base = self
        return lib.set_class(grad, (_MoistGrad, grad.__class__), "MOIST" + grad.__class__.__name__)

    Gradients = nuc_grad_method

    def density_fit(self, *args, **kwargs):
        base = self.undo_moist().density_fit(*args, **kwargs)
        return moist_for_scf(base, cavity=self.with_moist.cavity,
                             components=self.with_moist.components, parameters=self.with_moist.parameters)

    def _unsupported(self, *args, **kwargs):
        raise NotImplementedError("MOIST currently supports ground-state SCF and nuclear gradients only")

    newton = stability = gen_response = Hessian = _unsupported
    TDA = TDHF = TDDFT = CasidaTDDFT = _unsupported
    MP2 = CISD = CCSD = CASCI = CASSCF = _unsupported
    to_gpu = to_rhf = to_uhf = to_rks = to_uks = to_ghf = to_gks = _unsupported
    PCM = ddPCM = ddCOSMO = SMD = _unsupported


class _MoistGrad:
    def grad_elec(self, mo_energy=None, mo_coeff=None, mo_occ=None, atmlst=None):
        electronic = super().grad_elec(mo_energy, mo_coeff, mo_occ, atmlst)
        if mo_coeff is None:
            mo_coeff = self.base.mo_coeff
        if mo_occ is None:
            mo_occ = self.base.mo_occ
        dm = self.base.make_rdm1(mo_coeff, mo_occ)
        if self.base.with_moist.mol is not self.mol:
            self.base.with_moist.reset(self.mol)
        correction = self.base.with_moist.gradient(dm)
        if atmlst is not None:
            correction = correction[np.asarray(atmlst, dtype=int)]
        return electronic + correction

    def to_gpu(self, *args, **kwargs):
        raise NotImplementedError("MOIST GPU gradients are not supported")


def moist_for_scf(mf, *, cavity, components, parameters=None):
    """Attach MOIST to RHF/RKS/UHF/UKS without running SCF.

    Use ``DROP(lsf=SvdW(...))``, ``DROP(lsf=CFC(...))``,
    ``DROP(lsf=Isodensity(...))`` or ``ISwiG(...)``.
    Importing :mod:`moist.pyscf` also registers this function as ``mf.MOIST``.
    """
    from pyscf import lib, scf
    from pyscf.soscf.newton_ah import _CIAH_SOSCF

    if isinstance(mf, _MoistSCF):
        raise ValueError("MOIST is already attached; combine terms in components")
    if isinstance(mf, (scf.rohf.ROHF, _CIAH_SOSCF)) or not isinstance(mf, (scf.hf.RHF, scf.uhf.UHF)):
        raise NotImplementedError("MOIST supports RHF, RKS, UHF and UKS calculations")
    if getattr(mf, "with_solvent", None) is not None:
        raise ValueError("Cannot attach MOIST to a calculation with an existing solvent")
    if mf.mol.has_ecp():
        raise ValueError("effective core potentials are not supported")
    state = _MoistState(mf.mol, cavity=cavity, components=components, parameters=parameters)
    obj = mf.copy()
    obj.with_moist = state
    # Results from an earlier gas-phase calculation are not MOIST results.
    obj.converged = False
    obj.scf_summary = dict(mf.scf_summary)
    return lib.set_class(obj, (_MoistSCF, mf.__class__), "MOIST" + mf.__class__.__name__)


def _register():
    try:
        from pyscf import scf
    except ImportError:
        return
    scf.hf.SCF.MOIST = moist_for_scf


_register()
