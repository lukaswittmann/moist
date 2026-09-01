"""PySCF host bindings for moist.

moist never sees an AO basis: it hands back *adjoint weights* and expects the
host to finish the chain rule with its own density derivatives.  This module is
the reference implementation of that host side for PySCF, covering both the
plain solute-vdW cavity and the isodensity cavity whose surface follows the
electron density.

The level set moist is given is

.. math::  S(r) = \\mathrm{scale} \\cdot (\\rho_\\mathrm{iso} - \\rho(r))

so that the interior is negative and the exterior positive, matching the DROP
sign convention.

Two chain rules have to be completed by the host.  For the Fock matrix,

.. math::

    F_{\\mu\\nu} = \\sum_i q_i \\frac{\\partial \\phi_i}{\\partial P_{\\mu\\nu}}
    + \\sum_i \\Big[ w^{(0)}_i \\frac{\\partial S_i}{\\partial P_{\\mu\\nu}}
    + w^{(1)}_i \\cdot \\frac{\\partial \\nabla S_i}{\\partial P_{\\mu\\nu}}
    + w^{(2)}_i : \\frac{\\partial \\nabla^2 S_i}{\\partial P_{\\mu\\nu}} \\Big]

and for the nuclear gradient the same weights are contracted with
``dS/dR_A`` instead, alongside the AO-derivative part of ``dphi/dR_A``.  The
remaining routes -- the nuclear field, the surface motion, the A-matrix and the
anchor/switching geometry -- belong to moist.  A
:class:`~moist.interface.Evaluation` composes both halves into its ``fock`` and
``gradient`` results.

Conventions this module depends on, each verified against finite differences by
``test_pyscf.py``:

* ``phi`` is the **bare** point potential.  moist builds the nuclear half of the
  surface-motion adjoint from an unblurred ``Z_A (r_i - R_A)/r^3``; a
  Gaussian-blurred ``phi`` would be inconsistent with it.
* ``qefield_i = q_i * grad_r phi_elec(r_i)`` -- a gradient of the electronic
  potential.  Despite the C header calling it "the electronic field" it is not
  negated.
* The ``lsf`` weights are adjoints of the **scaled** level set, while the
  callback returns the density, so the host chain rule carries ``scale``
  (and the sign flip) that moist applied when it built ``S`` from ``rho``.
* Nuclear charges come from the structure's atomic numbers, so ECPs are not
  supported here.
"""

from __future__ import annotations

from typing import Optional
import warnings

import numpy as np

from .interface import (
    CavityDROPIsodensity,
    CavitySnapshot,
    CouplingChannel,
    CouplingTransaction,
    Electrostatics,
    Response,
    ModelComponentCPCM,
    SolvationCoupling,
    SolvationModel,
    Structure,
    _immutable_array,
)

__all__ = ["PySCFCoupling", "PySCFHost", "PySCFIsodensityHost", "solvated_rhf"]

#: Default isodensity contour in electrons/bohr^3.
DEFAULT_RHO_ISO = 4.0e-4


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


class PySCFHost:
    """Adapt a PySCF molecule to MOIST host and isodensity operations.

    The density matrix is mutable state (:attr:`dm`) because the isodensity
    surface follows the density: it must be current *before* every cavity
    build, and the level-set callback reads it on every point evaluation.

    :param mol: PySCF molecule.  Coordinates are read in bohr.
    :param rho_iso: Density contour defining the surface.
    :param scale: Constant multiplier moist applies to the level set; must match
        the value passed to :class:`~moist.interface.CavityDROPIsodensity`.
    """

    def __init__(
        self,
        mol,
        rho_iso: float = DEFAULT_RHO_ISO,
        scale: float = 1000.0,
    ) -> None:
        self.mol = mol
        self.rho_iso = float(rho_iso)
        self.scale = float(scale)
        self.dm: Optional[np.ndarray] = None
        self._gto_prefix = "GTOval_cart_deriv" if mol.cart else "GTOval_sph_deriv"
        self._aoslice = mol.aoslice_by_atom()

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------

    def structure(self) -> Structure:
        """Molecular structure in moist's representation (bohr)."""
        charges = self.mol.atom_charges()
        numbers = np.asarray(self.mol.atom_charges(), dtype=np.int64)
        if not np.array_equal(charges, numbers):
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
        """Density callback for :class:`~moist.interface.CavityDROPIsodensity`.

        Returns the bare ``rho`` and its spatial derivatives up to ``order``
        only, so the projection's value+gradient phase never pays for the
        density Hessian.  moist owns the level set built from it: it subtracts
        :attr:`rho_iso` and applies :attr:`scale` and the DROP sign convention.
        Doing any of that here as well would move the surface or square the
        scale -- and a squared scale leaves the zero level set, and hence the
        surface, unchanged while every adjoint comes back wrong by a factor of
        ``scale``.
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

    def make_cavity(self, **kwargs) -> CavityDROPIsodensity:
        """Deprecated compatibility factory for ``CavityDROPIsodensity(self)``."""
        warnings.warn(
            "host.make_cavity() is deprecated; use CavityDROPIsodensity(host, ...)",
            DeprecationWarning,
            stacklevel=2,
        )
        return CavityDROPIsodensity(self, **kwargs)

    def coupling(self, density_matrix: np.ndarray) -> "PySCFCoupling":
        """Bind one density matrix to this host for a coherent evaluation."""
        return PySCFCoupling(self, density_matrix)

    # ------------------------------------------------------------------
    # electrostatics
    # ------------------------------------------------------------------

    def surface_potential(self, coords: np.ndarray) -> np.ndarray:
        """Bare molecular electrostatic potential at the  grid points, ``(ngrid,)``."""
        coords = np.asarray(coords)
        dm = self._density_matrix()
        centers = self.mol.atom_coords()
        charges = self.mol.atom_charges().astype(float)
        delta = coords[:, None, :] - centers[None, :, :]
        dist = np.linalg.norm(delta, axis=2)
        phi_nuc = (charges[None, :] / dist).sum(axis=1)
        vgrids = self.mol.intor("int1e_grids", grids=coords)
        return phi_nuc - np.einsum("iuv,uv->i", vgrids, dm)

    def _grad_phi_elec(self, coords: np.ndarray) -> np.ndarray:
        """``grad_r phi_elec(r)`` at each grid point, ``(3, ngrid)``.

        The derivative with respect to the grid origin follows from
        translational invariance: shifting the operator center by ``d`` is the
        same as shifting both AO centers by ``-d``, giving
        ``d/dC (r_i|uv) = T_uv + T_vu`` with ``T`` the bra-derivative integral.
        """
        dm = self._density_matrix()
        tint = self.mol.intor("int1e_grids_ip", grids=coords)
        return -2.0 * np.einsum("kiuv,uv->ki", tint, dm)

    def _grad_phi_nuc(self, coords: np.ndarray) -> np.ndarray:
        """``grad_r phi_nuc(r)`` at each grid point, ``(3, ngrid)``."""
        centers = self.mol.atom_coords()
        charges = self.mol.atom_charges().astype(float)
        delta = coords[:, None, :] - centers[None, :, :]
        dist = np.linalg.norm(delta, axis=2)
        return -np.einsum("A,iAk,iA->ki", charges, delta, dist**-3)

    def qefield(self, coords: np.ndarray, q: np.ndarray) -> np.ndarray:
        """``q_i * grad_r phi_elec(r_i)``, Fortran ``(3, ngrid)``.

        For the **gradient** path.  Only the electronic half: moist rebuilds the
        nuclear half itself from the atomic numbers and adds it, forming
        ``w_xyz = qefield - q_i E_nuc``.  Passing the total here would count the
        nuclear field twice.
        """
        coords = np.asarray(coords)
        return np.asfortranarray(np.asarray(q)[None, :] * self._grad_phi_elec(coords))

    def surface_position_weights(self, coords: np.ndarray, q: np.ndarray) -> np.ndarray:
        """``q_i * grad_r phi_total(r_i)``, Fortran ``(3, ngrid)``.

        For the **potential** path, where it must be supplied as ``w_xyz``.
        When the density changes, the isodensity surface moves and ``phi_i =
        phi(r_i)`` changes with it; that route is the dominant part of the
        cavity response and moist cannot see it, because ``phi`` is the host's
        function.  Omitting it does not fail -- it silently returns ``lsf``
        weights that are missing their largest contribution.

        Unlike :meth:`qefield` this is the **total** potential gradient.  moist
        has no nuclear-field reconstruction on the potential path, so the
        nuclear part has to be included here.
        """
        coords = np.asarray(coords)
        gradient = self._grad_phi_nuc(coords) + self._grad_phi_elec(coords)
        return np.asfortranarray(np.asarray(q)[None, :] * gradient)

    def solve(self, model, coords: np.ndarray):
        """Supply electrostatics in the order the two derivative paths need.

        The charges are needed to build the response weights, but the response
        weights have to be in place before the potential is read, so the
        potential is supplied twice: once bare to obtain ``q``, then again with
        ``w_xyz`` (consumed by the potential path) and ``qefield`` (consumed by
        the gradient path).

        Returns ``(energy, response)``.
        """
        coords = np.asarray(coords)
        phi = self.surface_potential(coords)
        model.supply_electrostatics(phi)
        q = model.trace_response().electrostatics.surface_charge
        model.supply_electrostatics(
            phi,
            w_xyz=self.surface_position_weights(coords, q),
            qefield=self.qefield(coords, q),
        )
        return model.get_energy(), model.response()

    # ------------------------------------------------------------------
    # analytic derivatives
    # ------------------------------------------------------------------

    def fock(self, coords: np.ndarray, response, *, include_lsf: bool = True) -> np.ndarray:
        """Solvation contribution to the Fock matrix, ``(nao, nao)``.

        :param response: the :class:`~moist.interface.Response` from
            :meth:`SolvationModel.response`.
        :param include_lsf: contract the level-set adjoints.  Must be ``False``
            for a density-independent cavity, whose ``lsf`` weights describe a
            level set that has nothing to do with the electron density.
        """
        coords = np.asarray(coords)
        nao = self.mol.nao
        if response.electrostatics is None:
            # No electrostatic component, so there is no charge term at all --
            # distinct from a component that returned zero charges.
            fock = np.zeros((nao, nao))
        else:
            q = np.asarray(response.electrostatics.surface_charge)
            vgrids = self.mol.intor("int1e_grids", grids=coords)
            fock = -np.einsum("i,iuv->uv", q, vgrids)
        if include_lsf:
            fock += self._fock_lsf(coords, response)
        return fock

    def _fock_lsf(self, coords: np.ndarray, response) -> np.ndarray:
        """``dS/dP`` contracted with the level-set adjoints.

        With ``S = scale (rho_iso - rho)`` and ``rho = sum P_uv chi_u chi_v``,
        ``dS/dP_uv = -scale chi_u chi_v``, so every term is a weighted outer
        product of AO derivative blocks summed over the grid.
        """
        ao = self._ao(coords, 2)
        w0 = np.asarray(response.lsf.w_value)
        w1 = np.asarray(response.lsf.w_gradient)
        w2 = np.asarray(response.lsf.w_hessian)

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
        return -self.scale * fock

    def gradient(
        self,
        coords: np.ndarray,
        response,
        *,
        include_lsf: bool = True,
    ) -> np.ndarray:
        """Host-side nuclear gradient terms, Fortran ``(3, natm)``.

        These are exactly the routes that run through the AO basis and which
        moist therefore cannot see: the basis-center derivative of ``phi``, and
        -- for an isodensity cavity, whose level set reports zero nuclear
        partials by construction -- the basis-center derivative of the level
        set.  Add the result to
        :meth:`SolvationModel.get_gradient`.
        """
        coords = np.asarray(coords)
        if response.electrostatics is None:
            gradient = np.zeros((3, self.mol.natm))
        else:
            gradient = self._gradient_phi(
                coords, np.asarray(response.electrostatics.surface_charge)
            )
        if include_lsf:
            gradient += self._gradient_lsf(coords, response)
        return gradient

    def _gradient_phi(self, coords: np.ndarray, q: np.ndarray) -> np.ndarray:
        """``sum_i q_i d(phi_elec)_i/dR_A`` at fixed  grid points and fixed P."""
        dm = self._density_matrix()
        tint = self.mol.intor("int1e_grids_ip", grids=coords)
        # d(r_i|uv)/dR_A = -T[k,i,u,v] delta_{u in A} - T[k,i,v,u] delta_{v in A},
        # and phi_elec carries a further minus sign; P symmetry merges the two
        # halves into a single factor of two.
        weighted = 2.0 * np.einsum("i,kiuv,uv->ku", q, tint, dm)
        gradient = np.zeros((3, self.mol.natm))
        for atom in range(self.mol.natm):
            lo, hi = self._aoslice[atom, 2:4]
            gradient[:, atom] = weighted[:, lo:hi].sum(axis=1)
        return np.asfortranarray(gradient)

    def _gradient_lsf(self, coords: np.ndarray, response) -> np.ndarray:
        """Level-set adjoints contracted with ``dS/dR_A`` at fixed P.

        Because ``d/dR_A`` commutes with the spatial derivatives, every order is
        a spatial derivative of the single "displaced density"

        .. math:: u^{Ak}(r) = \\sum_{\\mu \\in A, \\nu} P_{\\mu\\nu}
                  \\partial_k \\chi_\\mu(r) \\chi_\\nu(r)

        with ``drho^{(m)}/dR_Ak = -2 d^m u^{Ak}`` and hence
        ``dS^{(m)}/dR_Ak = +2 scale d^m u^{Ak}``.
        """
        dm = self._density_matrix()
        ao = self._ao(coords, 3)
        w0 = np.asarray(response.lsf.w_value)
        w1 = np.asarray(response.lsf.w_gradient)
        w2 = np.asarray(response.lsf.w_hessian)

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

            per_ao *= 2.0 * self.scale
            for atom in range(self.mol.natm):
                lo, hi = self._aoslice[atom, 2:4]
                gradient[k, atom] = per_ao[lo:hi].sum()
        return np.asfortranarray(gradient)


class PySCFIsodensityHost(PySCFHost):
    """Deprecated compatibility name for :class:`PySCFHost`."""

    def __init__(self, *args, **kwargs) -> None:
        warnings.warn(
            "PySCFIsodensityHost is deprecated; use PySCFHost",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)


class PySCFCoupling(SolvationCoupling):
    """Adapter for one moist evaluation at a fixed PySCF density matrix.

    The adapter supplies whichever host channels the model requests.  For
    electrostatics it hides the two-pass surface-charge exchange; for GOSTSHYP
    it builds the Gaussian density moments.  It completes all returned adjoints
    into host Fock and nuclear-gradient contributions.

    Parameters
    ----------
    host
        A :class:`PySCFHost` providing the PySCF host operations and, for an
        isodensity cavity, its level-set source.
    density_matrix
        Density matrix held fixed for this evaluation.
    """

    channels = frozenset(
        {CouplingChannel.ELECTROSTATICS, CouplingChannel.GOSTSHYP}
    )

    def __init__(
        self,
        host: PySCFHost,
        density_matrix: np.ndarray,
    ) -> None:
        self.host = host
        density = np.ascontiguousarray(density_matrix)
        expected = (host.mol.nao_nr(), host.mol.nao_nr())
        if density.shape != expected:
            raise ValueError(f"density_matrix must have shape {expected}")
        self._density_matrix = _immutable_array(density)
        # The snapshot handed to fock()/gradient() deliberately carries values
        # only, so density dependence cannot be asked of it there; prepare()
        # runs first and leaves the answer here.
        self._include_lsf = True
        self._gostshyp = None

    @property
    def density_matrix(self) -> np.ndarray:
        """The fixed, immutable density associated with this evaluation."""
        return self._density_matrix

    @property
    def structure(self) -> Structure:
        return self.host.structure()

    def activate(self) -> None:
        # The isodensity callback runs during model.update(), before prepare().
        self.host.dm = self.density_matrix

    def _set_gostshyp_response(self, response) -> None:
        """Select host response state for a deprecated all-in-one wrapper."""
        self._gostshyp = response

    def prepare(self, transaction: CouplingTransaction) -> None:
        self._include_lsf = transaction.density_dependent
        if transaction.requires(CouplingChannel.ELECTROSTATICS):
            cavity = transaction.cavity
            coords = cavity.xyz.T
            phi = self.host.surface_potential(coords)

            def electrostatics(
                _cavity: CavitySnapshot,
                trace: Optional[Response],
            ) -> Electrostatics:
                if trace is None:
                    return Electrostatics(phi)
                return Electrostatics(
                    phi,
                    w_xyz=self.host.surface_position_weights(
                        coords, trace.electrostatics.surface_charge
                    ),
                    qefield=self.host.qefield(
                        coords, trace.electrostatics.surface_charge
                    ),
                )

            transaction.exchange_electrostatics(electrostatics)

        if transaction.requires(CouplingChannel.GOSTSHYP):
            if self._gostshyp is None:
                # Local import keeps the public PySCF adapter independent of
                # the optional host-side integral implementation at import time.
                from .gostshyp import _PySCFGostshyp

                self._gostshyp = _PySCFGostshyp(self.host)
            self._gostshyp._prepare(transaction, self.density_matrix)
        else:
            self._gostshyp = None

    def fock(
        self,
        cavity: CavitySnapshot,
        response: Response,
    ) -> np.ndarray:
        self.activate()
        fock = self.host.fock(
            cavity.xyz.T,
            response,
            include_lsf=self._include_lsf,
        )
        if self._gostshyp is not None:
            fock = fock + self._gostshyp._fock_from(
                response,
                include_cavity_response=False,
            )
        return fock

    def gradient(
        self,
        cavity: CavitySnapshot,
        response: Response,
        model_gradient,
    ) -> np.ndarray:
        self.activate()
        gradient = model_gradient() + self.host.gradient(
            cavity.xyz.T,
            response,
            include_lsf=self._include_lsf,
        )
        if self._gostshyp is not None:
            gradient = gradient + self._gostshyp._integral_nuclear_gradient(
                self.density_matrix,
                response,
            )
        return gradient

def solvated_rhf(
    mol,
    epsilon: float,
    *,
    rho_iso: float = DEFAULT_RHO_ISO,
    scale: float = 1000.0,
    conv_tol: float = 1e-13,
    conv_tol_grad: float = 1e-9,
    **cavity_kwargs,
):
    """Restricted Hartree-Fock in a CPCM isodensity cavity, solved to self-consistency.

    The surface follows the density, so the cavity is rebuilt from scratch on
    every SCF iteration.  Because :meth:`PySCFHost.fock` is the exact
    derivative of the solvation energy, the SCF remains a stationary-point
    search for ``E_HF + E_solv`` and the converged density is variational.

    Returns the converged PySCF mean-field object.
    """
    from pyscf import lib, scf

    host = PySCFHost(mol, rho_iso=rho_iso, scale=scale)

    class _SolvatedRHF(scf.hf.RHF):
        """RHF carrying the solvation response as a tagged extra potential."""

        def _solvent(self, dm):
            model = SolvationModel(
                CavityDROPIsodensity(host, **cavity_kwargs), [ModelComponentCPCM(epsilon)]
            )
            result = model.evaluate(coupling=host.coupling(dm))
            return result.energy, result.fock

        def get_veff(self, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
            veff = super().get_veff(mol, dm, dm_last, vhf_last, hermi)
            energy, potential = self._solvent(dm)
            return lib.tag_array(veff + potential, e_solv=energy, v_solv=potential)

        def energy_elec(self, dm=None, h1e=None, vhf=None):
            if dm is None:
                dm = self.make_rdm1()
            if vhf is None or getattr(vhf, "e_solv", None) is None:
                vhf = self.get_veff(self.mol, dm)
            energy, coulomb = super().energy_elec(dm, h1e, vhf)
            # The base class folded 0.5 Tr[P v_solv] into the Coulomb term, which
            # is the double-counting correction for a *linear* response.  The
            # isodensity solvation energy is not that -- the cavity itself moves
            # with P -- so undo it and add the energy moist actually reported.
            double_counted = 0.5 * float(np.einsum("uv,uv->", dm, vhf.v_solv))
            return energy - double_counted + vhf.e_solv, coulomb - double_counted

    mean_field = _SolvatedRHF(mol)
    # A finite difference of the converged energy divides by ~12h, so residual
    # SCF error is amplified by ~10^3.  The orbital-gradient threshold is the
    # one that matters and has to be set explicitly: PySCF defaults it to
    # sqrt(conv_tol), which would leave it at ~1e-7 and dominate the residual.
    # conv_tol itself is kept at 1e-13 -- for a -75 Ha energy 1e-14 is already
    # at machine precision, and the converged energy is unchanged either way.
    mean_field.conv_tol = conv_tol
    mean_field.conv_tol_grad = conv_tol_grad
    mean_field.kernel()
    return mean_field
