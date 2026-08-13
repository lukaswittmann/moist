"""GOSTSHYP hydrostatic pressure on a moist isodensity cavity.

GOSTSHYP simulates hydrostatic pressure by placing an unnormalized Gaussian
potential on every cavity grid-point and fixing its amplitude so the force the
wall exerts on the electron density matches ``p_inp`` times the grid-point area:

.. math::

    \\omega_j &= \\pi \\ln 2 / a_j                                        \\\\
    g_{\\mu\\nu,j} &= \\langle \\mu | e^{-\\omega_j |r-C_j|^2} | \\nu \\rangle  \\\\
    f_{\\mu\\nu,j} &= n_j \\cdot \\nabla_r g_{\\mu\\nu,j}                    \\\\
    p_j &= p_\\mathrm{inp} a_j / \\tilde f_j                              \\\\
    E^\\mathrm{GOST} &= \\sum_j p_j \\tilde g_j

with ``gtilde``/``ftilde`` the density contractions of ``g``/``f``.

This module is the host part and forms the AO-basis three-center integrals

The energy, the amplitudes, and the derivatives with respect to the cavity
parameters are in moist's ``gostshyp`` solvation-model component in
``src/moist/model/component/gostshyp.f90``

This requires one exchange per cavity update::
    host -> moist    gt, Pt, Mt, Rt        Gaussian moments of the density
    moist -> host    w_overlap, w_normal_deriv  amplitudes for the host's Fock

where the moments are ``<G>``, ``<(r-C) G>``, ``<(r-C)(r-C) G>`` and
``<(r-C) |r-C|^2 G>``, and the host computes its Fock contribution as
``F += sum_j [w_overlap_j g_j + w_normal_deriv_j f_j]``

Note ``ftilde`` is not exchanged: moist derives it as ``-2 omega_j (n_j . Pt_j)``
from the same moments it differentiates (normal convention)

Every per-grid-point Gaussian carries a normalization constant ``N_j`` that
cancels exactly between energy and Fock, so it is not explicitly formed.
Only the relative s/p/d/f angular constants are restored, which is what makes
``f = n . grad g`` hold exactly. Those constants are the one convention that
remains host-side, and they are pinned against an independent quadrature in
``test_gostshyp.py``.

The nuclear gradient has three routes:

* integral: the AO centers move at a frozen surface.  Host-side, from
  ``int3c1e_ip1``.
* field: the density, and hence the level set, moves with the nuclei.
  Host-side, through ``PySCFHost._gradient_lsf`` contracted with the
  level-set adjoints moist returns.
* surface: the grid points are dragged by their anchor atoms and the
  surface itself responds.  Entirely moist's, through the reverse-mode path;
  this module never sees a surface weight.

Unit of ``pressure`` in Hartree/bohr^3.
"""

from __future__ import annotations

import math
from typing import Optional
import warnings

import numpy as np
from pyscf import gto

from .interface import (
    CavityDROPIsodensity,
    CavitySnapshot,
    CouplingTransaction,
    Evaluation,
    Response,
    GostshypMoments,
    ModelComponentGOSTSHYP,
    SolvationModel,
)
from .pyscf import PySCFHost

__all__ = ["GostshypModel", "GostshypWall", "GPA_TO_AU"]

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

#: Mirrors ``overlap_floor`` in ``src/moist/model/component/gostshyp.f90``.
#:
#: The floor that decides which grid points still carry a usable density overlap
#: belongs beside the amplitudes it masks, and that is where it acts: nothing on
#: the energy, Fock or gradient path reads this copy.  It exists only for
#: :meth:`GostshypModel.effective_volume`, which has to report a
#: pressure-*independent* quantity and so cannot recover the mask from
#: amplitudes that vanish with the pressure.  The linearity test pins the two
#: copies against each other wherever the pressure is nonzero.
_OVERLAP_FLOOR = 1.0e-9


def _fakemol_gaussians(coords: np.ndarray, exponents: np.ndarray, angl: int) -> gto.Mole:
    """One coefficient-1 GTO shell of angular momentum ``angl`` per grid-point."""

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
    """Three-center one-electron integrals over a Gaussian-per-grid-point fakemol."""

    fakemol = _fakemol_gaussians(centers, omega, angl)
    nbas = mol.nbas
    shls_slice = (0, nbas, 0, nbas, nbas, nbas + fakemol.nbas)
    return (mol + fakemol).intor(intor, shls_slice=shls_slice)


class _PySCFGostshyp:
    """PySCF-side Gaussian integrals and response for a GOSTSHYP component."""

    def __init__(self, host: PySCFHost) -> None:
        self._initialize_host(host)

    def _initialize_host(self, host: PySCFHost) -> None:
        self.host = host
        self.mol = host.mol
        self.energy = 0.0
        self.ngrid = 0
        self._response: Optional[Response] = None
        self._evaluation: Optional[Evaluation] = None
        self._traces: Optional[tuple[np.ndarray, np.ndarray]] = None
        self._c2s: Optional[np.ndarray] = None
        self._G: Optional[np.ndarray] = None
        self._F: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # cavity and grid point geometry
    # ------------------------------------------------------------------

    def _set_grid_points(self, result: CavitySnapshot) -> None:
        """Snapshot the live cavity into the arrays the integrals need."""

        self.centers = np.ascontiguousarray(result.xyz.T, dtype=np.float64)
        self.areas = np.ascontiguousarray(result.a, dtype=np.float64)
        self.ngrid = int(self.centers.shape[0])
        self.nsph = int(result.nsph)

        normals = np.array(result.normal0.T, dtype=np.float64, order="C", copy=True)
        norm = np.linalg.norm(normals, axis=1)
        good = norm > 0.0
        normals[good] /= norm[good, None]
        self.normals = normals

        # omega_j = pi ln2 / a_j, matching `gaussian_width` in the component; a
        # degenerate zero-area grid point is inert.
        with np.errstate(divide="ignore", invalid="ignore"):
            omega = np.pi * math.log(2.0) / self.areas
        omega[~np.isfinite(omega)] = 0.0
        self.omega = omega

    def _prepare(self, transaction: CouplingTransaction, dm: np.ndarray) -> None:
        """Build and supply Gaussian moments for one coupling transaction."""
        self._set_grid_points(transaction.cavity)
        self._build_integrals()
        gt, pt, mt, rt = self._surface_moments(self._density_matrix_cart(dm))
        transaction.supply_gostshyp(
            GostshypMoments(
                gt,
                np.asfortranarray(pt.T),
                np.asfortranarray(mt.transpose(1, 2, 0)),
                np.asfortranarray(rt.T),
            )
        )
        self._traces = (
            gt,
            -2.0
            * self.omega
            * np.einsum("ja,ja->j", self.normals, pt, optimize=True),
        )

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

    def f_vector(self) -> np.ndarray:
        """``grad_r g_uv,j`` before the normal projection, ``(nao, nao, ngrid, 3)``.

        Recomputed rather than cached: only the projected ``f`` is needed to
        build a Fock matrix, and this block is three times its size.  Kept as a
        method because it is what pins the p-shell angular constant -- the
        projection that produces ``f`` cannot distinguish a wrong constant from
        a wrong normal.
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
        self._F = np.einsum(
            "uvja,ja->uvj", self.f_vector(), self.normals, optimize=True
        )

    def traces(
        self, dm: np.ndarray, *, centers=None, omega=None, normals=None
    ) -> tuple[np.ndarray, np.ndarray]:
        """``(gtilde, ftilde)`` from this module's own moments.

        The host's copy of the two traces moist works from.  It is a diagnostic
        and a self-check, never an input to the energy: ``ftilde`` is
        deliberately not part of the exchange, so this reproduces moist's
        ``-2 omega (n . Pt)`` rather than reading it back.  The linearity test
        pins the two against each other.

        The surface parameters default to the live cavity and are overridable
        for finite differences; see :meth:`_surface_moments`.

        Always rebuilt from scratch, which costs four dense integral blocks.
        The live-surface values are cached by :meth:`update`, so anything that
        just wants *those* -- :meth:`effective_volume`, say -- should read
        :attr:`live_traces` instead of calling this with default arguments.
        """

        omega = self.omega if omega is None else omega
        normals = self.normals if normals is None else normals
        gt, pt, _mt, _rt = self._surface_moments(
            self._density_matrix_cart(dm), centers, omega
        )
        ftilde = -2.0 * np.asarray(omega) * np.einsum(
            "ja,ja->j", normals, pt, optimize=True
        )
        return gt, ftilde

    def _density_matrix_cart(self, dm: np.ndarray) -> np.ndarray:
        """Density matrix in the cartesian AO basis the fakemol blocks use.

        BLAS leaves dirty floating-point status flags behind for these shapes --
        the SIMD tail reads padding lanes -- so numpy reports divide-by-zero and
        overflow from a product that performs neither.  Suppressed here for the
        same reason :meth:`~moist.pyscf.PySCFHost.lsf` suppresses them.
        """

        c2s = self._cart2sph
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            return c2s @ np.asarray(dm) @ c2s.T

    # ------------------------------------------------------------------
    # the Gaussian moments -- this module's product
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

    def _surface_moments(self, dm_cart: np.ndarray, centers=None, omega=None):
        """``(gt, Pt, Mt, Rt)`` — the s/p/d/f Gaussian moments of the density.

        Numpy-ordered: ``(ngrid,)``, ``(ngrid, 3)``, ``(ngrid, 3, 3)`` and
        ``(ngrid, 3)``.  :meth:`update` transposes them into the Fortran layout
        the C API expects.

        ``centers``/``omega`` default to the live surface.  They are overridable
        so a finite difference can rebuild the moments at a displaced surface
        without touching cached state -- the Gaussians live *on* the grid
        points, so a reference that held them fixed would differentiate the
        wrong function.
        """

        centers = self.centers if centers is None else centers
        omega = self.omega if omega is None else omega
        ngrid = int(np.asarray(omega).size)

        g_cart = _int3c1e(self.mol, centers, omega, 0)
        p_cart = _int3c1e(self.mol, centers, omega, 1)
        d_cart = _int3c1e(self.mol, centers, omega, 2)
        f_cart = _int3c1e(self.mol, centers, omega, 3)
        ncart = g_cart.shape[0]

        gt = np.einsum("pqj,pq->j", g_cart, dm_cart, optimize=True)
        pt = self._contract_p_moments(dm_cart, p_cart.reshape(ncart, ncart, ngrid, 3))
        mt = self._contract_d_moments(dm_cart, d_cart.reshape(ncart, ncart, ngrid, 6))
        rt = self._contract_f_rho2_moments(dm_cart, f_cart.reshape(ncart, ncart, ngrid, 10))
        return gt, pt, mt, rt

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------

    def _require_response(self) -> Response:
        if self._response is None:
            raise RuntimeError("call update(dm) successfully before reading GOSTSHYP results")
        return self._response

    @property
    def response(self) -> Response:
        """Every channel group moist returned for the last :meth:`update`.

        Carries both the amplitudes the host contracts into its Fock matrix and
        the level-set adjoints in ``lsf`` that
        :meth:`~moist.pyscf.PySCFHost._fock_lsf` consumes.
        """

        return self._require_response()

    @property
    def amplitudes(self) -> Optional[np.ndarray]:
        """The per-grid-point amplitudes ``p_j`` of the last :meth:`update`."""

        return None if self._response is None else self._response.gostshyp.w_overlap

    @property
    def inactive_count(self) -> int:
        """Grid points the component switched off, out of :attr:`ngrid`.

        The component drops points whose density overlap is a relative
        round-off away from zero, and says nothing about it.  A sizeable count
        is normal -- 7 of 71 for water at 50 GPa, 24 of 159 for fluoroacetate --
        so this is a number to watch across a trajectory rather than a pass/fail
        one: a systematically wrong moment supply shrinks every ``ftilde`` and
        shows up here as points quietly leaving the surface.

        Derived from ``ftilde`` rather than from the amplitudes, for the reason
        :meth:`effective_volume` gives: the amplitudes carry the pressure, so at
        ``p_inp = 0`` they report every point as dropped.
        """

        self._require_response()
        _, ftilde = self.live_traces
        floor = _OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
        return int(np.count_nonzero(np.abs(ftilde) <= floor))

    @property
    def live_traces(self) -> tuple[np.ndarray, np.ndarray]:
        """``(gtilde, ftilde)`` on the surface of the last :meth:`update`.

        Kept from the moment build that update already performed, so reading
        them costs nothing.  Use :meth:`traces` only for a *different* surface
        or density.
        """

        if self._traces is None:
            raise RuntimeError("call update(dm) successfully before reading GOSTSHYP results")
        return self._traces

    def effective_volume(self) -> float:
        """``E / p_inp`` (eq 11).  Not the cavity volume.

        Evaluated as ``sum_j a_j gtilde_j / ftilde_j`` rather than by dividing
        the energy: that ratio is what the quantity *is*, and it stays well
        defined at ``p_inp = 0``, where the energy vanishes with the pressure
        but the volume it reports does not.

        Reads the traces cached by :meth:`update` rather than rebuilding them:
        this is a diagnostic, and it should not cost four integral blocks.  It
        is the one place the host reproduces moist's ``ftilde``, never an input
        to the energy -- and the linearity test pins it against
        ``energy / pressure``, which fails if the two conventions ever drift.
        """

        self._require_response()
        gt, ftilde = self.live_traces
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = self.areas * gt / ftilde
        # The mask cannot be recovered from the amplitudes: they are proportional
        # to the pressure, so at p = 0 they say every grid point is inactive.
        floor = _OVERLAP_FLOOR * float(np.max(np.abs(ftilde), initial=0.0))
        active = np.abs(ftilde) > floor
        return float(np.sum(np.where(active, ratio, 0.0)))

    def fock(self, dm: np.ndarray, *, include_cavity_response: bool = True) -> np.ndarray:
        """Solvation contribution to the Fock matrix, ``(nao, nao)``.

        :param include_cavity_response: contract the level-set adjoints.  With
            ``False`` this is the eq-16 matrix at a frozen surface, which is
            *not* the derivative of the energy an isodensity cavity reports —
            the surface moves with the density.
        """

        return self._fock_from(
            self._require_response(),
            include_cavity_response=include_cavity_response,
        )

    def _fock_from(
        self,
        response: Response,
        *,
        include_cavity_response: bool,
    ) -> np.ndarray:
        """Contract a captured response with the current integral blocks."""
        fock = np.einsum(
            "j,uvj->uv", response.gostshyp.w_overlap, self._G, optimize=True
        )
        fock += np.einsum(
            "j,uvj->uv", response.gostshyp.w_normal_deriv, self._F, optimize=True
        )
        fock = 0.5 * (fock + fock.T)
        if include_cavity_response:
            fock = fock + self.host._fock_lsf(self.centers, response)
        return fock

    # ------------------------------------------------------------------
    # nuclear gradient
    # ------------------------------------------------------------------

    def _integral_nuclear_gradient(
        self,
        dm: np.ndarray,
        response: Optional[Response] = None,
    ) -> np.ndarray:
        """AO centers move, surface frozen.  Fortran ``(3, natm)``."""

        if response is None:
            response = self._require_response()
        dm_cart = self._density_matrix_cart(dm)
        ncart = self._cart2sph.shape[0]

        ip1_g = _int3c1e(self.mol, self.centers, self.omega, 0, "int3c1e_ip1_cart")
        ip1_p = _int3c1e(self.mol, self.centers, self.omega, 1, "int3c1e_ip1_cart")
        ip1_g = ip1_g.reshape(3, ncart, ncart, self.ngrid)
        ip1_p = ip1_p.reshape(3, ncart, ncart, self.ngrid, 3)
        ip1_f = -_S_OVER_P_NORM * 2.0 * np.einsum(
            "j,xpqja,ja->xpqj", self.omega, ip1_p, self.normals, optimize=True
        )

        kernel = np.einsum(
            "j,xpqj->xpq", response.gostshyp.w_overlap, ip1_g, optimize=True
        )
        kernel += np.einsum(
            "j,xpqj->xpq", response.gostshyp.w_normal_deriv, ip1_f, optimize=True
        )
        row = np.einsum("xpq,pq->xp", kernel, dm_cart, optimize=True)

        # Cartesian AO rows fold onto atoms through the cartesian slices.
        gradient = np.zeros((3, self.mol.natm))
        aoslice = self.mol.aoslice_by_atom(self.mol.ao_loc_nr(cart=True))
        for atom in range(self.mol.natm):
            lo, hi = aoslice[atom, 2:4]
            # dg/dR = -ip1, and both AO legs contribute equally for symmetric P.
            gradient[:, atom] = -2.0 * row[:, lo:hi].sum(axis=1)
        return np.asfortranarray(gradient)

    def _field_nuclear_gradient(self, dm: np.ndarray) -> np.ndarray:
        """The level set's own dependence on the nuclei, through the AO basis."""

        self.host.dm = dm
        return self.host._gradient_lsf(self.centers, self._require_response())

    def _nuclear_gradient_from(
        self,
        dm: np.ndarray,
        response: Response,
        model_gradient,
    ) -> np.ndarray:
        """Assemble every gradient route for one captured evaluation."""
        self.host.dm = dm
        return (
            self._integral_nuclear_gradient(dm, response)
            + self.host._gradient_lsf(self.centers, response)
            + model_gradient()
        )

class GostshypModel(_PySCFGostshyp):
    """Deprecated all-in-one wrapper for an isodensity GOSTSHYP model.

    Use :class:`~moist.interface.ModelComponentGOSTSHYP` in a regular
    :class:`~moist.interface.SolvationModel` and evaluate it with
    :meth:`moist.pyscf.PySCFHost.coupling` instead.
    """

    def __init__(
        self,
        host: PySCFHost,
        pressure: float,
        **cavity_kwargs,
    ) -> None:
        warnings.warn(
            "GostshypModel is deprecated; use ModelComponentGOSTSHYP as a "
            "SolvationModel component",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(host)
        self.pressure = float(pressure)
        self.component = ModelComponentGOSTSHYP(self.pressure)
        self.model = SolvationModel(
            CavityDROPIsodensity(host, **cavity_kwargs), [self.component]
        )

    def evaluate(self, dm: np.ndarray) -> Evaluation:
        """Return one coherent energy, Fock, gradient, and cavity evaluation.

        The cached results are dropped *before* the rebuild.  The host density
        and the cavity are already mutated by the time anything downstream can
        fail, so keeping the previous state would leave every accessor
        reporting coherent-looking numbers for a surface that no longer exists.
        """

        self._response = None
        self._traces = None
        self._evaluation = None
        self.energy = 0.0

        coupling = self.host.coupling(dm)
        coupling._set_gostshyp_response(self)
        try:
            result = self.model.evaluate(coupling=coupling)
        except Exception:
            # ``prepare`` necessarily builds transient host data before the
            # native model can finish.  The compatibility wrapper publishes
            # none of it when the enclosing evaluation fails.
            self._response = None
            self._traces = None
            self._evaluation = None
            self.energy = 0.0
            raise
        assert self._traces is not None
        self.energy = result.energy
        self._response = result.response
        self._evaluation = result
        return result

    def update(self, dm: np.ndarray) -> float:
        """Compatibility method returning only :meth:`evaluate`'s energy."""
        return self.evaluate(dm).energy

    def _surface_nuclear_gradient(self) -> np.ndarray:
        """The cavity's own response, contracted by moist in reverse mode.

        Covers both the rigid drag of each grid point by its anchor atom and the
        surface's response at a frozen field.  The component states ``w_a``,
        ``w_xyz`` and ``w_n``; the cavity contracts them against its own nuclear
        derivatives without ever building a Jacobian, so nothing here folds
        normals or Hessians by hand.
        """

        self._require_response()
        return self.model.gradient()

    def nuclear_gradient(self, dm: np.ndarray) -> np.ndarray:
        """Total ``dE^GOST/dR`` at fixed ``dm``, Fortran ``(3, natm)``.

        Three routes a nuclear displacement takes: the AO centers move at a
        frozen surface (integral), the density and hence the level set moves
        (field), and the cavity responds (surface).  At a converged variational
        SCF this is also the total pressure-wall gradient, with no orbital
        response.
        """

        return self._nuclear_gradient_from(
            dm,
            self._require_response(),
            self.model.gradient,
        )


class GostshypWall(GostshypModel):
    """Deprecated compatibility name for :class:`GostshypModel`."""
