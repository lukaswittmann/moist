"""Tests for two-stage PySCF Hessian-vector products and analytic SCF Hessians."""

from __future__ import annotations

import os

import numpy as np
import pytest

try:
    from pyscf import gto
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from moist import (
    CFC,
    Context,
    DROP,
    DROPParameters,
    ISwiG,
    ISwiGParameters,
    Isodensity,
    IsodensityParameters,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentDPCM,
    ModelComponentDPCMOriginal,
    ModelComponentGOSTSHYP,
    ModelComponentIEFPCM,
    ModelComponentLPB,
    ModelComponentNonlocal,
    ModelComponentPV,
    ModelComponentSLIC,
    ModelComponentSSVPE,
    PCMParameters,
    SvdW,
)
from moist._pyscf_hessian import _DensityTables, _ElectrostaticTables
from moist.interface import DensityTangentResponse, PotentialAdjointTangentResponse
from moist.pyscf import PySCFHost, PySCFSolvation  # registers .MOIST()

CONTEXT = Context()
DROP_26 = DROPParameters(nleb=26)
DENSE = PCMParameters(solver="dense")
DENSE_POINT = PCMParameters(solver="dense", point_potential=True)

_SYM2 = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
_SYM3 = (
    (0, 0, 0), (0, 0, 1), (0, 0, 2),
    (0, 1, 1), (0, 1, 2), (0, 2, 2),
    (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2),
)


def finite_difference_gradient_tangent(gradient_fn, positions, direction, *, step=1.0e-3):
    """Four-point central difference of ``gradient_fn(pos)`` along ``direction``."""
    total = None
    for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
        val = np.asarray(gradient_fn(positions + off * step * direction), dtype=np.float64)
        total = wt * val if total is None else total + wt * val
    return total / (12.0 * step)


def finite_difference_hessian(gradient_fn, positions, *, step=1.0e-3):
    """Assemble ``(natm, natm, 3, 3)`` Hessian from four-point finite differences of gradients."""
    natm = positions.shape[0]
    hess = np.empty((natm, natm, 3, 3), dtype=np.float64)
    for atom in range(natm):
        for axis in range(3):
            direction = np.zeros_like(positions)
            direction[atom, axis] = 1.0
            hess[atom, :, axis, :] = finite_difference_gradient_tangent(
                gradient_fn, positions, direction, step=step
            )
    return hess


def _make_mol(atom: str, *, basis: str = "sto-3g", charge: int = 0, spin: int = 0, cart: bool = False):
    return gto.M(atom=atom, unit="bohr", basis=basis, charge=charge, spin=spin, cart=cart, verbose=0)


_ASYM_WATER = "O 0.0 0.0 0.0; H 1.8 0.2 0.1; H -0.3 1.7 0.4"
_WATER_NE = "O 0.0 0.0 0.0; H 0.0 -1.4 1.1; H 0.0 1.4 1.1; Ne 0.0 0.0 4.2"


# ---------------------------------------------------------------------------
# Host integral and density table derivative checks (def2-SVP)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gaussian", [True, False], ids=["gaussian_xi", "point"])
@pytest.mark.parametrize("cart", [False, True], ids=["sph", "cart"])
def test_electrostatic_tables_match_fd_def2svp(gaussian: bool, cart: bool) -> None:
    """Verify host electrostatic tangent and completion tables against 4-pt FD on def2-SVP."""
    mol = _make_mol("O 0 0 0; H 0 1.4 -1.1; H 0 -1.4 -1.1", basis="def2-svp", cart=cart)
    mf = mol.RHF().set(conv_tol=1e-12).run()
    dm = mf.make_rdm1()
    host = PySCFHost(mol)
    host.dm = dm

    coords = np.array([[1.2, -0.7, 0.9], [-0.8, 1.5, 0.4], [0.2, 0.1, -1.6]], dtype=np.float64)
    xi = np.array([1.2, 0.85, 1.45], dtype=np.float64) if gaussian else None
    tables = _ElectrostaticTables(host, coords, xi)

    rng = np.random.default_rng(13)
    natm, ngrid = mol.natm, len(coords)
    dirs = rng.standard_normal((1, natm, 3))
    dm_dir = rng.standard_normal(dm.shape)
    dm_dir = 0.5 * (dm_dir + dm_dir.T)
    dm_dirs = dm_dir[None, :, :]
    d_xyz = rng.standard_normal((1, ngrid, 3))
    d_xi = rng.standard_normal((1, ngrid)) if gaussian else None

    missing = frozenset(
        ["d_phi", "d_dphi_dr", "d_dphi_dxi"] if gaussian else ["d_phi", "d_dphi_dr"]
    )
    tangents = tables.potential_tangent(dm, dirs, dm_dirs, d_xyz, d_xi, missing)

    def eval_pot(t: float):
        mol_t = mol.set_geom_(mol.atom_coords() + t * dirs[0], unit="bohr", inplace=False)
        h_t = PySCFHost(mol_t)
        h_t.dm = dm + t * dm_dir
        c_t = coords + t * d_xyz[0]
        xi_t = None if xi is None else xi + t * d_xi[0]
        phi = h_t.surface_potential(c_t, xi_t)
        dphi_dr = h_t.surface_potential_gradient(c_t, xi_t)
        if xi_t is None:
            return phi, dphi_dr
        return phi, dphi_dr, h_t._dphi_dxi(c_t, xi_t)

    step = 2.0e-4
    fd_phi = fd_grad = fd_dxi = 0.0
    for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
        vals = eval_pot(off * step)
        fd_phi = fd_phi + wt * vals[0] / (12.0 * step)
        fd_grad = fd_grad + wt * vals[1] / (12.0 * step)
        if gaussian:
            fd_dxi = fd_dxi + wt * vals[2] / (12.0 * step)

    np.testing.assert_allclose(tangents["d_phi"][0], fd_phi, atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(tangents["d_dphi_dr"][0], fd_grad, atol=1e-9, rtol=1e-8)
    if gaussian:
        np.testing.assert_allclose(tangents["d_dphi_dxi"][0], fd_dxi, atol=1e-9, rtol=1e-8)

    # Completion check: differentiate w_phi * phi_el(R, P, r(t), xi(t)) wrt R and P
    w_phi = rng.standard_normal(ngrid)
    dw_phi = rng.standard_normal((1, ngrid))
    item = PotentialAdjointTangentResponse(
        w_phi=w_phi,
        dw_phi=dw_phi,
        d_xyz=d_xyz,
        d_xi=d_xi,
    )
    hvp_nuc = np.zeros((1, natm, 3), dtype=np.float64)
    fock_tan = np.zeros((1, mol.nao, mol.nao), dtype=np.float64)
    tables.complete_potential(dm, dirs, dm_dirs, item, hvp_nuc, fock_tan)

    def eval_first_derivs(t: float):
        mol_t = mol.set_geom_(mol.atom_coords() + t * dirs[0], unit="bohr", inplace=False)
        h_t = PySCFHost(mol_t)
        h_t.dm = dm + t * dm_dir
        c_t = coords + t * d_xyz[0]
        xi_t = None if xi is None else xi + t * d_xi[0]
        w_t = w_phi + t * dw_phi[0]
        f = h_t.fock_potential(c_t, xi_t, w_t)
        g = h_t._gradient_phi(c_t, w_t, xi=xi_t)
        return g, f

    fd_hvp = fd_fock = 0.0
    for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
        g_t, f_t = eval_first_derivs(off * step)
        fd_hvp = fd_hvp + wt * g_t / (12.0 * step)
        fd_fock = fd_fock + wt * f_t / (12.0 * step)

    np.testing.assert_allclose(hvp_nuc[0], fd_hvp, atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(fock_tan[0], fd_fock, atol=1e-9, rtol=1e-8)


def test_density_tables_match_fd_def2svp() -> None:
    """Verify bare _DensityTables tangent and completion terms against 4-pt FD on def2-SVP."""
    from moist.interface import DensityResponse

    mol = _make_mol(_ASYM_WATER, basis="def2-svp")
    mf = mol.RHF().set(conv_tol=1e-12).run()
    dm = mf.make_rdm1()
    host = PySCFHost(mol)
    host.dm = dm

    coords = np.array([[1.1, -0.5, 0.7], [-0.9, 1.3, -0.3]], dtype=np.float64)
    tables = _DensityTables(host, coords)
    rng = np.random.default_rng(19)
    natm, ngrid = mol.natm, len(coords)
    dirs = rng.standard_normal((1, natm, 3))
    dm_dir = rng.standard_normal(dm.shape)
    dm_dir = 0.5 * (dm_dir + dm_dir.T)
    dm_dirs = dm_dir[None, :, :]
    d_xyz = rng.standard_normal((1, ngrid, 3))

    tan = tables.density_tangent(
        dm, dirs, dm_dirs, frozenset(["d_rho", "d_grad_rho", "d_hess_rho", "d_d3_rho"])
    )

    def bare_rho_jet(mol_cur, dm_cur, pts):
        ao = mol_cur.eval_gto("GTOval_sph_deriv3", pts)
        from moist.pyscf import _component_index

        i1 = [_component_index((a,)) for a in range(3)]
        i2 = [[_component_index((a, b)) for b in range(3)] for a in range(3)]
        # Ao eval_gto for multiple grid points: (ncomp, ngrid, nao)
        rho = np.einsum("iu,uv,iv->i", ao[0], dm_cur, ao[0])
        grad = np.zeros((len(pts), 3))
        hess = np.zeros((len(pts), 6))
        d3 = np.zeros((len(pts), 10))
        for a in range(3):
            grad[:, a] = 2.0 * np.einsum("iu,uv,iv->i", ao[i1[a]], dm_cur, ao[0])
        for k, (a, b) in enumerate(_SYM2):
            hess[:, k] = 2.0 * (
                np.einsum("iu,uv,iv->i", ao[i2[a][b]], dm_cur, ao[0])
                + np.einsum("iu,uv,iv->i", ao[i1[a]], dm_cur, ao[i1[b]])
            )
        for k, (a, b, c) in enumerate(_SYM3):
            d3[:, k] = 2.0 * (
                np.einsum("iu,uv,iv->i", ao[_component_index((a, b, c))], dm_cur, ao[0])
                + np.einsum("iu,uv,iv->i", ao[i2[a][b]], dm_cur, ao[i1[c]])
                + np.einsum("iu,uv,iv->i", ao[i2[a][c]], dm_cur, ao[i1[b]])
                + np.einsum("iu,uv,iv->i", ao[i2[b][c]], dm_cur, ao[i1[a]])
            )
        return rho, grad, hess, d3

    step = 2.0e-4
    fd_j = [0.0, 0.0, 0.0, 0.0]
    for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
        mol_t = mol.set_geom_(mol.atom_coords() + off * step * dirs[0], unit="bohr", inplace=False)
        vals = bare_rho_jet(mol_t, dm + off * step * dm_dir, coords)
        for q in range(4):
            fd_j[q] = fd_j[q] + wt * vals[q] / (12.0 * step)

    np.testing.assert_allclose(tan["d_rho"][0], fd_j[0], atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(tan["d_grad_rho"][0], fd_j[1], atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(tan["d_hess_rho"][0], fd_j[2], atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(tan["d_d3_rho"][0], fd_j[3], atol=1e-9, rtol=1e-8)

    # Completion test against FD of PySCFHost._gradient_lsf / _fock_lsf
    w_rho = rng.standard_normal(ngrid)
    w_grad = rng.standard_normal((ngrid, 3))
    w_hess_raw = rng.standard_normal((ngrid, 3, 3))
    w_hess = 0.5 * (w_hess_raw + w_hess_raw.transpose(0, 2, 1))
    dw_rho = rng.standard_normal((1, ngrid))
    dw_grad = rng.standard_normal((1, ngrid, 3))
    dw_hess = rng.standard_normal((1, ngrid, 6))
    dw_hess_mat = np.empty((ngrid, 3, 3))
    for k, (a, b) in enumerate(_SYM2):
        dw_hess_mat[:, a, b] = dw_hess_mat[:, b, a] = dw_hess[0, :, k]

    item = DensityTangentResponse(
        w_rho=w_rho,
        w_grad_rho=w_grad,
        w_hess_rho=w_hess,
        dw_rho=dw_rho,
        dw_grad_rho=dw_grad,
        dw_hess_rho=dw_hess,
        d_xyz=d_xyz,
    )
    hvp_nuc = np.zeros((1, natm, 3), dtype=np.float64)
    fock_tan = np.zeros((1, mol.nao, mol.nao), dtype=np.float64)
    tables.complete_density(dm, dirs, dm_dirs, item, hvp_nuc, fock_tan)

    fd_hvp = fd_fock = 0.0
    for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
        t = off * step
        mol_t = mol.set_geom_(mol.atom_coords() + t * dirs[0], unit="bohr", inplace=False)
        h_t = PySCFHost(mol_t)
        h_t.dm = dm + t * dm_dir
        c_t = coords + t * d_xyz[0]
        d_resp = DensityResponse(
            w_rho=w_rho + t * dw_rho[0],
            w_grad_rho=w_grad + t * dw_grad[0],
            w_hess_rho=w_hess + t * dw_hess_mat,
        )
        fd_hvp = fd_hvp + wt * h_t._gradient_lsf(c_t, d_resp) / (12.0 * step)
        fd_fock = fd_fock + wt * h_t._fock_lsf(c_t, d_resp) / (12.0 * step)

    np.testing.assert_allclose(hvp_nuc[0], fd_hvp, atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(fock_tan[0], fd_fock, atol=1e-9, rtol=1e-8)


# ---------------------------------------------------------------------------
# Fixed-density PySCFSolvation.hvp & PySCFSolvation.hessian across models
# ---------------------------------------------------------------------------


def _build_solv_case(mol, dm, cavity, components):
    return PySCFSolvation(mol, cavity, components, context=CONTEXT)


@pytest.mark.parametrize(
    "component_factory",
    [
        lambda: [ModelComponentPV(1.2e-3)],
        lambda: [ModelComponentCPCM(32.0, parameters=DENSE)],
        lambda: [ModelComponentCOSMO(32.0, parameters=DENSE)],
        lambda: [ModelComponentIEFPCM(32.0, parameters=DENSE)],
        lambda: [ModelComponentSSVPE(32.0, parameters=DENSE)],
        lambda: [ModelComponentLPB(32.0, 0.1, parameters=DENSE)],
        lambda: [ModelComponentCPCM(32.0, parameters=DENSE_POINT)],
    ],
    ids=["PV", "CPCM", "COSMO", "IEFPCM", "SSVPE", "LPB", "CPCM_point"],
)
@pytest.mark.parametrize(
    "lsf_factory",
    [
        lambda: DROP(lsf=SvdW(), parameters=DROP_26),
        lambda: DROP(lsf=CFC(), parameters=DROP_26),
        lambda: DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26),
    ],
    ids=["SvdW", "CFC", "Isodensity"],
)
def test_fixed_density_hvp_matches_fd_all_scalar_pcms_and_lsfs(component_factory, lsf_factory) -> None:
    """Verify PySCFSolvation.hvp against 4-pt FD across PV, all 5 scalar PCMs, and all 3 DROP LSFs."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g")
    mf = mol.RHF().set(conv_tol=1e-12).run()
    dm = mf.make_rdm1()
    rng = np.random.default_rng(29)
    dirs = rng.standard_normal((2, mol.natm, 3))
    dm_dir = rng.standard_normal(dm.shape)
    dm_dir = 0.5 * (dm_dir + dm_dir.T)
    dm_dirs = np.stack([dm_dir, -0.5 * dm_dir], axis=0)

    cavity = lsf_factory()
    comps = component_factory()
    solv = PySCFSolvation(mol, cavity, comps, context=CONTEXT)
    hvp_nuc, fock_tan = solv.hvp(dm, dirs, dm_dirs)

    step = 2.0e-4
    fd_nuc = np.zeros_like(hvp_nuc)
    fd_fock = np.zeros_like(fock_tan)
    for j in range(2):
        for off, wt in ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)):
            t = off * step
            mol_t = mol.set_geom_(mol.atom_coords() + t * dirs[j], unit="bohr", inplace=False)
            dm_t = dm + t * dm_dirs[j]
            s_t = PySCFSolvation(mol_t, cavity, comps, context=CONTEXT)
            res_t = s_t.evaluate(dm_t)
            grad_t = s_t.gradient(dm_t)
            fd_nuc[j] += wt * grad_t / (12.0 * step)
            fd_fock[j] += wt * res_t.fock / (12.0 * step)

    np.testing.assert_allclose(hvp_nuc, fd_nuc, atol=1e-8, rtol=1e-7)
    np.testing.assert_allclose(fock_tan, fd_fock, atol=1e-8, rtol=1e-7)


def test_fixed_density_hessian_multibranch_and_symmetry() -> None:
    """Verify PySCFSolvation.hessian on multi-branch H2O+Ne (CPCM + PV + Isodensity)."""
    mol = _make_mol(_WATER_NE, basis="sto-3g")
    mf = mol.RHF().set(conv_tol=1e-12).run()
    dm = mf.make_rdm1()
    cavity = DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26)
    comps = [ModelComponentCPCM(32.0, parameters=DENSE), ModelComponentPV(1.0e-3)]
    solv = PySCFSolvation(mol, cavity, comps, context=CONTEXT)

    hess = solv.hessian(dm)
    assert hess.shape == (mol.natm, mol.natm, 3, 3)
    np.testing.assert_allclose(hess, hess.transpose(1, 0, 3, 2), atol=1e-10, rtol=1e-10)

    def fixed_dm_grad(pos):
        mol_t = mol.set_geom_(pos, unit="bohr", inplace=False)
        return PySCFSolvation(mol_t, cavity, comps, context=CONTEXT).gradient(dm)

    fd_hess = finite_difference_hessian(fixed_dm_grad, mol.atom_coords(), step=2.0e-4)
    np.testing.assert_allclose(hess, fd_hess, atol=1e-8, rtol=1e-7)


def test_reciprocity_and_block_symmetry_iefpcm_pv_isodensity() -> None:
    """Verify H_RR = H_RR^T, H_PP = H_PP^T, and H_RP = H_PR^T via unit nuclear and density HVPs."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g")
    mf = mol.RHF().set(conv_tol=1e-12).run()
    dm = mf.make_rdm1()
    nao = mol.nao_nr()
    nat = mol.natm

    cavity = DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26)
    comps = [ModelComponentIEFPCM(10.0, parameters=DENSE), ModelComponentPV(1.2e-3)]
    solv = PySCFSolvation(mol, cavity, comps, context=CONTEXT)

    unit_r = np.eye(3 * nat, dtype=np.float64).reshape(3 * nat, nat, 3)
    h_rr, f_r = solv.hvp(dm, unit_r)
    h_rr_mat = h_rr.reshape(3 * nat, 3 * nat)

    sym_basis = []
    for u in range(nao):
        for v in range(u, nao):
            mat = np.zeros((nao, nao), dtype=np.float64)
            mat[u, v] = mat[v, u] = 1.0
            sym_basis.append(mat)
    dm_dirs = np.stack(sym_basis, axis=0)
    h_pr, f_pp = solv.hvp(dm, dm_dirs=dm_dirs)

    n_p = len(sym_basis)
    h_rp_mat = np.empty((3 * nat, n_p), dtype=np.float64)
    h_pr_mat = h_pr.reshape(n_p, 3 * nat)
    h_pp_mat = np.empty((n_p, n_p), dtype=np.float64)
    for k, ek in enumerate(sym_basis):
        h_rp_mat[:, k] = np.einsum("juv,uv->j", f_r, ek)
        h_pp_mat[:, k] = np.einsum("kuv,uv->k", f_pp, ek)

    np.testing.assert_allclose(h_rr_mat, h_rr_mat.T, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(h_pp_mat, h_pp_mat.T, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(h_rp_mat, h_pr_mat.T, atol=1e-10, rtol=1e-10)


# ---------------------------------------------------------------------------
# Relaxed SCF mf.Hessian() vs 4-pt FD of mf.Gradients()
# ---------------------------------------------------------------------------


def _build_base_scf(mol, method: str, *, df: bool = False):
    base = getattr(mol, method)()
    if method.endswith("KS"):
        base.xc = "pbe" if method == "RKS" else "lda,vwn"
        base.grids.level = 4
        base.grids.prune = None
    if df:
        base = base.density_fit()
    base.conv_tol = 1e-13
    base.conv_tol_grad = 1e-10
    base.max_cycle = 200
    return base


def _run_scf(mol, method: str, cavity, comps, *, df: bool = False):
    base = _build_base_scf(mol, method, df=df)
    mf = base.MOIST(cavity=cavity, components=comps, context=CONTEXT)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-10
    mf.conv_tol_cpscf = 1e-11
    mf.max_cycle = 200
    mf.run()
    assert mf.converged
    return mf


def _run_gas_scf(mol, method: str, *, df: bool = False):
    mf = _build_base_scf(mol, method, df=df)
    mf.conv_tol_cpscf = 1e-11
    mf.run()
    assert mf.converged
    return mf


@pytest.mark.parametrize(
    "method,charge,spin,df,cavity,comps",
    [
        (
            "RHF",
            0,
            0,
            False,
            DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26),
            [ModelComponentCPCM(32.0, parameters=DENSE), ModelComponentPV(1.0e-3)],
        ),
        (
            "RHF",
            0,
            0,
            True,
            DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26),
            [ModelComponentIEFPCM(16.0, parameters=DENSE)],
        ),
        (
            "RKS",
            0,
            0,
            False,
            DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26),
            [ModelComponentCPCM(32.0, parameters=DENSE)],
        ),
        (
            "UHF",
            1,
            1,
            False,
            DROP(lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)), parameters=DROP_26),
            [ModelComponentCPCM(32.0, parameters=DENSE), ModelComponentPV(1.0e-3)],
        ),
        (
            "UKS",
            1,
            1,
            False,
            DROP(lsf=SvdW(), parameters=DROP_26),
            [ModelComponentCPCM(32.0, parameters=DENSE)],
        ),
    ],
    ids=["RHF_iso_cpcm_pv", "DF_RHF_iso_iefpcm", "RKS_pbe_iso_cpcm", "UHF_iso_cpcm_pv", "UKS_svdw_cpcm"],
)
def test_relaxed_scf_hessian_matches_fd(method, charge, spin, df, cavity, comps) -> None:
    """Verify full relaxed mf.Hessian().kernel() against 4-pt central FD of mf.Gradients()."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g", charge=charge, spin=spin)
    mf = _run_scf(mol, method, cavity, comps, df=df)
    hess_obj = mf.Hessian()
    calls = []
    hess_obj.kernel_callback = calls.append
    hess = hess_obj.kernel()
    assert hess.shape == (mol.natm, mol.natm, 3, 3)
    assert len(calls) > 0 and all(n > 0 for n in calls)
    is_dft = method in ("RKS", "UKS")
    sym_atol, sym_rtol = (5e-8, 1e-7) if is_dft else (1e-9, 1e-9)
    np.testing.assert_allclose(hess, hess.transpose(1, 0, 3, 2), atol=sym_atol, rtol=sym_rtol)

    def grad_at(pos):
        mol_t = mol.set_geom_(pos, unit="bohr", inplace=False)
        mf_t = _run_scf(mol_t, method, cavity, comps, df=df)
        return mf_t.Gradients().kernel()

    fd_hess = finite_difference_hessian(grad_at, mol.atom_coords(), step=1.0e-3)
    if is_dft:
        np.testing.assert_allclose(hess, fd_hess, atol=5e-5, rtol=1e-4)
        gas_hess = _run_gas_scf(mol, method, df=df).Hessian().kernel()

        def gas_grad_at(pos):
            mol_t = mol.set_geom_(pos, unit="bohr", inplace=False)
            mf_t = _run_gas_scf(mol_t, method, df=df)
            return mf_t.Gradients().kernel()

        solv_hess = hess - gas_hess
        assert np.max(np.abs(solv_hess)) > 1e-3
        solv_fd_hess = fd_hess - finite_difference_hessian(
            gas_grad_at, mol.atom_coords(), step=1.0e-3
        )
        np.testing.assert_allclose(solv_hess, solv_fd_hess, atol=2e-6, rtol=1e-5)
    else:
        np.testing.assert_allclose(hess, fd_hess, atol=5e-7, rtol=1e-6)


def test_hessian_atmlst_gen_hop_and_state_wrapper() -> None:
    """Verify atmlst subsetting, gen_hop HVP, and _MoistState.{hvp, hessian}."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g")
    cavity = DROP(lsf=SvdW(), parameters=DROP_26)
    comps = [ModelComponentCPCM(32.0, parameters=DENSE), ModelComponentPV(1.0e-3)]
    mf = _run_scf(mol, "RHF", cavity, comps)

    hess_obj = mf.Hessian()
    full = hess_obj.kernel()
    sub = hess_obj.kernel(atmlst=[0, 2])
    np.testing.assert_allclose(sub, full[np.ix_([0, 2], [0, 2])], atol=1e-11, rtol=1e-11)

    hop, hdiag = hess_obj.gen_hop()
    rng = np.random.default_rng(42)
    vec = rng.standard_normal(3 * mol.natm)
    np.testing.assert_allclose(
        hop(vec),
        full.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm) @ vec,
        atol=1e-12,
    )
    np.testing.assert_allclose(hdiag, np.einsum("aaxx->ax", full).ravel(), atol=1e-12)

    dm = mf.make_rdm1()
    dirs = rng.standard_normal((2, mol.natm, 3))
    hvp_state, f_state = mf.with_moist.hvp(dm, dirs)
    hess_state = mf.with_moist.hessian(dm)
    assert hvp_state.shape == (2, mol.natm, 3) and f_state.shape == (2, mol.nao, mol.nao)
    np.testing.assert_allclose(hess_state, mf.with_moist.solvation.hessian(dm), atol=1e-12)


# ---------------------------------------------------------------------------
# Guard checks for unsupported models and SCF references
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "component",
    [
        ModelComponentDPCM(32.0),
        ModelComponentDPCMOriginal(32.0),
        ModelComponentSLIC(32.0, 0.2, 0.5, 0.8, 0.1, 0.0),
        ModelComponentNonlocal(32.0, 1.8, 3.0, 1.0),
        ModelComponentGOSTSHYP(1.0e-5),
    ],
    ids=["DPCM", "DPCMOriginal", "SLIC", "Nonlocal", "GOSTSHYP"],
)
def test_unsupported_components_raise_in_pyscf_hessian(component) -> None:
    """Normal-trace PCMs and GOSTSHYP raise NotImplementedError in PySCF Hessian coupling."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g")
    mf = _run_scf(mol, "RHF", DROP(lsf=SvdW(), parameters=DROP_26), [component])
    with pytest.raises(NotImplementedError):
        mf.Hessian().kernel()


def test_unsupported_cavity_and_reference_guards(monkeypatch) -> None:
    """ISwiG cavities, unconverged SCF, ECP/X2C/SGX raise clean errors."""
    mol = _make_mol(_ASYM_WATER, basis="sto-3g")
    unconv = mol.RHF().MOIST(
        cavity=DROP(lsf=SvdW(), parameters=DROP_26),
        components=[ModelComponentCPCM(32.0)],
        context=CONTEXT,
    )
    with pytest.raises(ValueError, match="converge"):
        unconv.Hessian()

    mf_iswig = _run_scf(mol, "RHF", ISwiG(parameters=ISwiGParameters(nleb=26)), [ModelComponentCPCM(32.0)])
    with pytest.raises(NotImplementedError):
        mf_iswig.Hessian().kernel()

    mf_ok = _run_scf(mol, "RHF", DROP(lsf=SvdW(), parameters=DROP_26), [ModelComponentCPCM(32.0)])
    monkeypatch.setattr(mf_ok.mol, "has_ecp", lambda: True)
    with pytest.raises(NotImplementedError):
        mf_ok.Hessian()
    monkeypatch.undo()

    monkeypatch.setattr(mf_ok, "with_x2c", object(), raising=False)
    with pytest.raises(NotImplementedError):
        mf_ok.Hessian()
    monkeypatch.undo()

    monkeypatch.setattr(mf_ok, "with_df", object(), raising=False)
    with pytest.raises(NotImplementedError):
        mf_ok.Hessian()
