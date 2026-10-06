"""Public PySCF wrapper: variational SCF, gradients and state isolation."""

import os
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

try:
    from pyscf import gto
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from moist import ModelComponentCPCM, ModelComponentPV, ModelComponentGOSTSHYP
from moist.interface import SolvationModel
from moist.pyscf import CFC, DROP, ISwiG, Isodensity, SvdW  # registers .MOIST()
from moist import (
    Context, DROPParameters, ISwiGParameters, IsodensityParameters, PCMParameters,
    CustomRadii,
)
from moist.pyscf import PySCFHost, PySCFSolvation

#: Run context shared by every cavity and model in this module
CONTEXT = Context()

CAVITIES = {"svdw-drop": DROP(lsf=SvdW(), parameters=DROPParameters(nleb=26)),
            "cfc-drop": DROP(lsf=CFC(), parameters=DROPParameters(nleb=26)),
            "iswig": ISwiG(parameters=ISwiGParameters(nleb=26)),
            "rho-drop": DROP(lsf=Isodensity(), parameters=DROPParameters(nleb=26))}


FD_ATOL = 1e-10
FD_RTOL = 1e-9
FD_SCF_ATOL = 5e-10  # Total gradients through displaced SCF solutions
FD_SCF_RTOL = 5e-9
FD_DENSITY_STEP = 1e-3  # Density-matrix direction amplitude
FD_NUCLEAR_STEP = 1e-3  # Bohr
FD_SCF_STEP = 1e-2  # Bohr
FD_SCF_ENERGY_TOL = 1e-14
FD_SCF_GRAD_TOL = 1e-10


def assert_fd(analytic, reference, atol=FD_ATOL, rtol=FD_RTOL):
    assert np.isfinite(analytic) and np.isfinite(reference)
    assert abs(analytic - reference) <= max(atol, rtol * abs(reference))


def fd4(energy, step):
    minus2, minus1, plus1, plus2 = [energy(k * step) for k in (-2, -1, 1, 2)]
    return ((plus1 - minus1) * 8 - (plus2 - minus2)) / (12 * step)


def fd6(energy, step):
    minus3, minus2, minus1, plus1, plus2, plus3 = [
        energy(k * step) for k in (-3, -2, -1, 1, 2, 3)]
    return (45 * (plus1 - minus1) - 9 * (plus2 - minus2) +
            (plus3 - minus3)) / (60 * step)


@pytest.fixture
def mol():
    return gto.M(atom="O 0 0 0; H 0 -1.4 1.1; H 0 1.4 1.1",
                 unit="bohr", basis="sto-3g", verbose=0)


def attach(mf, cavity="svdw-drop", components=None):
    if components is None:
        components = [ModelComponentCPCM(32.0)]
    return mf.MOIST(cavity=CAVITIES[cavity], components=components, context=CONTEXT)


@pytest.mark.parametrize("cavity", ["svdw-drop", "cfc-drop", "iswig", "rho-drop"])
def test_names_and_deferred_construction(mol, cavity):
    base = mol.RHF().set(conv_tol=2e-10, max_cycle=71)
    mf = attach(base, cavity)
    assert not hasattr(base, "with_moist")
    assert mf.with_moist.result is None
    assert mf.mo_coeff is None
    assert mf.conv_tol == base.conv_tol and mf.max_cycle == 71
    assert mf.with_moist.cavity.name == CAVITIES[cavity].name
    dm = base.get_init_guess()
    result = mf.with_moist.evaluate(dm)
    assert np.isfinite(result.energy)
    assert result.fock.shape == dm.shape


@pytest.mark.parametrize("method", ["RHF", "RKS", "UHF", "UKS"])
@pytest.mark.parametrize("cavity", ["iswig", "rho-drop"])
def test_scf_energy_fock_and_spin_contract(mol, method, cavity):
    if method.startswith("U"):
        mol = mol.copy()
        mol.charge, mol.spin = 1, 1
        mol.build()
    base = getattr(mol, method)()
    if method.endswith("KS"):
        base.xc = "lda,vwn"
        base.grids.level = 0
    mf = attach(base, cavity, [ModelComponentCPCM(32), ModelComponentPV(1e-5)])
    dm = base.get_init_guess()
    base_v = base.get_veff(dm=dm)
    v = mf.get_veff(dm=dm)
    np.testing.assert_allclose(v, base_v, atol=1e-12)
    result = mf.with_moist.evaluate(dm)
    np.testing.assert_allclose(mf.get_fock(dm=dm, vhf=v), base.get_fock(dm=dm, vhf=base_v) + result.fock,
                               atol=1e-12)
    assert mf.energy_elec(dm, vhf=v)[0] == pytest.approx(base.energy_elec(dm, vhf=base_v)[0] + result.energy)
    if method.endswith("KS"):
        assert v.exc == base_v.exc  # Checks that DFT exc metadata is preserved.
    assert mf.scf_summary["e_moist"] == result.energy


@pytest.mark.parametrize("cavity,builds", [("svdw-drop", 1), ("rho-drop", 2)])
def test_cache_by_density_value_and_geometry(mol, monkeypatch, cavity, builds):
    calls = []
    update = SolvationModel.update
    def counted(self, structure):
        calls.append(structure)
        return update(self, structure)
    monkeypatch.setattr(SolvationModel, "update", counted)
    mf = attach(mol.RHF(), cavity)
    dm = mf.get_init_guess()
    first = mf.with_moist.evaluate(dm)
    assert mf.with_moist.evaluate(dm.copy()) is first
    dm[0, 0] += 1e-3
    second = mf.with_moist.evaluate(dm)
    assert second is not first
    assert len(calls) == builds
    assert np.isfinite(mf.with_moist.gradient(dm)).all()
    mol.set_geom_(mol.atom_coords() + [0, 0, .1], unit="Bohr")
    assert mf.with_moist.evaluate(dm) is not second
    assert len(calls) == builds + 1
    current = mf.with_moist.cavity
    mf.with_moist.set(cavity=replace(
        current, parameters=replace(current.parameters, nleb=50)))
    assert mf.with_moist.result is None
    mf.with_moist.evaluate(dm)
    assert len(calls) == builds + 2


def test_copy_reset_and_configuration(mol):
    mf = attach(mol.RHF(), "rho-drop")
    dm = mf.get_init_guess()
    result = mf.with_moist.evaluate(dm)
    other = mf.copy()
    assert other.with_moist is not mf.with_moist
    other_current = other.with_moist.cavity
    other.with_moist.set(cavity=replace(
        other_current, lsf=Isodensity(parameters=IsodensityParameters(rho_iso=8e-4))))
    assert mf.with_moist.result is result
    assert mf.with_moist.cavity.lsf.parameters == Isodensity().parameters
    other.with_moist.evaluate(dm)
    assert other.with_moist.solvation.model is not mf.with_moist.solvation.model
    other.reset(mol.copy())
    assert other.with_moist.result is None
    assert mf.with_moist.result is result
    with pytest.raises(FrozenInstanceError):
        mf.with_moist.cavity.parameters.nleb = 50


@pytest.mark.parametrize("cavity", ["svdw-drop", "rho-drop"])
def test_nonlinear_fock_is_energy_derivative(mol, cavity):
    mf = attach(mol.RHF(), cavity, [ModelComponentCPCM(32), ModelComponentPV(1e-4)])
    dm = mf.get_init_guess()
    rng = np.random.default_rng(45)
    direction = rng.normal(size=dm.shape)
    direction += direction.T
    direction /= np.linalg.norm(direction)
    analytic = np.einsum("ij,ji", mf.with_moist.evaluate(dm).fock, direction)
    numerical = fd4(lambda offset: mf.with_moist.evaluate(dm + offset*direction).energy,
                    FD_DENSITY_STEP)
    assert_fd(analytic, numerical)


@pytest.mark.parametrize("method", ["RHF", "RKS", "UHF", "UKS"])
@pytest.mark.parametrize("cavity", ["svdw-drop", "rho-drop"])
def test_total_gradient_and_subset(mol, method, cavity):
    if method.startswith("U"):
        mol.charge, mol.spin = 1, 1
        mol.build()
    base = getattr(mol, method)().set(conv_tol=FD_SCF_ENERGY_TOL,
                                      conv_tol_grad=FD_SCF_GRAD_TOL, max_cycle=150)
    if method.endswith("KS"):
        base.xc = "lda,vwn"
        base.grids.level = 0
    mf = attach(base, cavity, [ModelComponentCPCM(32), ModelComponentPV(1e-5)]).run()
    assert mf.converged
    grad = mf.nuc_grad_method()
    if method.endswith("KS"):
        grad.grid_response = True
    analytic = grad.kernel()
    subset = grad.kernel(atmlst=[2, 0])
    np.testing.assert_allclose(subset, analytic[[2, 0]], atol=1e-10)
    def energy(offset):
        xyz = mol.atom_coords()
        xyz[1, 2] += offset
        moved = mol.set_geom_(xyz, unit="Bohr", inplace=False)
        displaced = mf.copy().reset(moved).run(dm0=mf.make_rdm1())
        assert displaced.converged
        return displaced.e_tot
    assert_fd(analytic[1, 2], fd6(energy, FD_SCF_STEP), FD_SCF_ATOL, FD_SCF_RTOL)


@pytest.mark.parametrize("cavity", ["svdw-drop", "rho-drop"])
def test_direct_scf_and_density_fitting(mol, cavity):
    energies = []
    for direct in [True, False]:
        mf = attach(mol.RHF().set(direct_scf=direct, conv_tol=1e-11), cavity).run()
        assert mf.converged
        energies.append(mf.e_tot)
    assert energies[0] == pytest.approx(energies[1], abs=1e-10)
    before = attach(mol.RHF().density_fit(), cavity).run(conv_tol=1e-11)
    after = attach(mol.RHF(), cavity).density_fit().run(conv_tol=1e-11)
    assert before.converged and after.converged
    assert before.e_tot == pytest.approx(after.e_tot, abs=1e-10)
    np.testing.assert_allclose(before.Gradients().kernel(), after.Gradients().kernel(), atol=1e-9)


@pytest.mark.parametrize("cavity", ["iswig", "rho-drop"])
def test_gradient_scanner(mol, cavity):
    mf = attach(mol.RHF().set(conv_tol=1e-11), cavity).run()
    scanner = mf.Gradients().as_scanner()
    assert scanner.base.with_moist is not mf.with_moist
    xyz = mol.atom_coords()
    xyz[1, 2] += .02
    moved = mol.set_geom_(xyz, unit="Bohr", inplace=False)
    energy, gradient = scanner(moved)
    reference = mf.copy().reset(moved).run()
    assert energy == pytest.approx(reference.e_tot, abs=1e-10)
    np.testing.assert_allclose(gradient, reference.Gradients().kernel(), atol=1e-7)


@pytest.mark.parametrize("cavity", ["svdw-drop", "cfc-drop", "iswig", "rho-drop"])
def test_gostshyp_composition(mol, cavity):
    mf = attach(mol.RHF(), cavity, [ModelComponentCPCM(32), ModelComponentGOSTSHYP(1e-4)])
    dm = mf.get_init_guess()
    result = mf.with_moist.evaluate(dm)
    assert np.isfinite(result.energy)
    assert np.isfinite(result.fock).all()
    analytic_gradient = mf.with_moist.gradient(dm)
    assert np.isfinite(analytic_gradient).all()
    direction = np.eye(dm.shape[0])
    density_fd = fd4(lambda offset: mf.with_moist.evaluate(dm + offset*direction).energy,
                     FD_DENSITY_STEP)
    assert_fd(np.einsum("ij,ji", result.fock, direction), density_fd)
    def energy(offset):
        xyz = mol.atom_coords()
        xyz[1, 2] += offset
        moved = mol.set_geom_(xyz, unit="Bohr", inplace=False)
        return mf.copy().reset(moved).with_moist.evaluate(dm).energy
    assert_fd(analytic_gradient[1, 2], fd4(energy, FD_NUCLEAR_STEP))


def test_unsupported_operations_are_explicit(mol, monkeypatch):
    with pytest.raises(TypeError, match="cavity must"):
        mol.RHF().MOIST(cavity="invalid", components=[ModelComponentCPCM(32)], context=CONTEXT)
    with pytest.raises(TypeError, match="keyword"):
        ISwiG(rho_iso=4e-4)
    with pytest.raises(TypeError, match="keyword"):
        Isodensity(source=object())
    mf = attach(mol.RHF())
    with pytest.raises(ValueError, match="already attached"):
        attach(mf)
    with pytest.raises(ValueError, match="existing solvent"):
        attach(mol.RHF().PCM())
    with pytest.raises(NotImplementedError):
        attach(mol.ROHF())
    with pytest.raises(NotImplementedError):
        attach(mol.RHF().newton())
    monkeypatch.setattr(mol, "has_ecp", lambda: True)
    with pytest.raises(ValueError, match="effective core potentials"):
        attach(mol.RHF())
    for name in ["newton", "stability", "gen_response", "Hessian", "TDA", "TDHF", "TDDFT", "MP2", "CCSD", "to_gpu"]:
        with pytest.raises(NotImplementedError):
            getattr(mf, name)()


def test_explicit_molecule_does_not_poison_default_evaluation(mol):
    mf = attach(mol.RHF())
    dm = mf.get_init_guess()
    expected = mf.get_veff(dm=dm).e_moist
    xyz = mol.atom_coords()
    xyz[1, 2] += .1
    moved = mol.set_geom_(xyz, unit="Bohr", inplace=False)
    mf.get_veff(mol=moved, dm=dm)
    assert mf.get_veff(dm=dm).e_moist == pytest.approx(expected, abs=1e-12)


def test_basis_change_clears_cache(mol):
    mf = attach(mol.RHF())
    dm = mf.get_init_guess()
    result = mf.with_moist.evaluate(dm)
    mol.basis = "6-31g"
    mol.build()
    changed = mf.with_moist.evaluate(mf.get_init_guess())
    assert changed is not result
    assert changed.fock.shape == (mol.nao_nr(), mol.nao_nr())


def test_no_interaction_matches_gas_phase(mol):
    base = mol.RHF().run(conv_tol=1e-11)
    wrapped = attach(mol.RHF(), "rho-drop", [ModelComponentCPCM(1)]).run(conv_tol=1e-11)
    assert wrapped.e_tot == pytest.approx(base.e_tot, abs=1e-10)
    np.testing.assert_allclose(wrapped.Gradients().kernel(), base.Gradients().kernel(), atol=1e-8)


def test_explicit_molecule_does_not_poison_gradient(mol):
    mf = attach(mol.RHF()).run(conv_tol=1e-11)
    expected = mf.Gradients().kernel()
    xyz = mol.atom_coords()
    xyz[1, 2] += .1
    moved = mol.set_geom_(xyz, unit="Bohr", inplace=False)
    mf.get_veff(mol=moved, dm=mf.make_rdm1())
    np.testing.assert_allclose(mf.Gradients().kernel(), expected, atol=1e-12)


def test_reusable_cavity_and_lsf_configuration(mol):
    surface = Isodensity()
    config = DROP(lsf=surface, parameters=DROPParameters(nleb=26))
    first = mol.RHF().MOIST(cavity=config, components=[ModelComponentCPCM(32)], context=CONTEXT)
    second = mol.RHF().MOIST(cavity=config, components=[ModelComponentCPCM(32)], context=CONTEXT)
    first.with_moist.set(cavity=replace(
        config, parameters=replace(config.parameters, nleb=50)))
    assert config.parameters.nleb == 26
    assert second.with_moist.cavity.parameters.nleb == 26
    with pytest.raises(FrozenInstanceError):
        config.parameters.nleb = 50


@pytest.mark.parametrize("make", [
    lambda: DROP(), lambda: DROP(SvdW()), lambda: DROP(lsf=None),
    lambda: DROP(lsf="SvdW"), lambda: ISwiG(lsf=SvdW()),
])
def test_cavity_and_lsf_settings_are_separate(make):
    with pytest.raises(TypeError):
        make()


def test_lsf_change_invalidates_results_and_preserves_grid_options(mol):
    mf = attach(mol.RHF(), "rho-drop")
    dm = mf.get_init_guess()
    old = mf.with_moist.evaluate(dm)
    current = mf.with_moist.cavity
    mf.with_moist.set(cavity=replace(
        current, lsf=Isodensity(parameters=IsodensityParameters(rho_iso=8e-4))))
    assert mf.with_moist.result is None
    assert mf.with_moist.cavity.parameters.nleb == 26
    new = mf.with_moist.evaluate(dm)
    assert new.energy != pytest.approx(old.energy, abs=1e-8)
    with pytest.raises(FrozenInstanceError):
        mf.with_moist.cavity.lsf.parameters.rho_iso = 4e-4
    with pytest.raises(TypeError, match="keyword"):
        mf.with_moist.set(rho_iso=4e-4)
    assert mf.with_moist.result is new


@pytest.mark.parametrize("density_dependent", [False, True])
def test_core_configuration_is_shared_with_pyscf(mol, density_dependent):
    import moist

    assert DROP is moist.DROP and Isodensity is moist.Isodensity
    parameters = DROPParameters(nleb=26, proj_level=2, tolerance=2e-10)
    surface = Isodensity(parameters=IsodensityParameters(rho_iso=8e-4)) if density_dependent else SvdW()
    radii = CustomRadii([3., 2.5, 2.5])
    config = DROP(lsf=surface, parameters=parameters, radii=radii)
    terms = [ModelComponentCPCM(32, parameters=PCMParameters(solver="lu"))]
    mf = mol.RHF().MOIST(cavity=config, components=terms, context=CONTEXT)
    dm = mf.get_init_guess()
    result = mf.with_moist.evaluate(dm)
    reference = PySCFSolvation(mol, config, terms, context=CONTEXT)
    expected = reference.evaluate(dm)
    assert result.energy == pytest.approx(expected.energy, abs=1e-12)
    np.testing.assert_allclose(result.fock, expected.fock, atol=1e-12)
    np.testing.assert_allclose(mf.with_moist.gradient(dm), reference.gradient(dm), atol=1e-11)
    assert mf.with_moist.solvation.model.cavity.configuration == config
    updated = replace(config, parameters=replace(parameters, nleb=50))
    mf.with_moist.set(cavity=updated)
    assert mf.with_moist.result is None
    assert config.parameters.nleb == 26
    assert mf.copy().with_moist.cavity == updated


def test_density_fit_preserves_context(mol):
    context = Context(verbosity=0, debug=True)
    mf = mol.RHF().MOIST(cavity=CAVITIES["iswig"],
                         components=[ModelComponentPV(1e-4)], context=context)
    fitted = mf.density_fit()
    assert fitted.with_moist.context is context
    assert fitted.with_moist.cavity == mf.with_moist.cavity
    assert fitted.with_moist.components == mf.with_moist.components


def test_host_rejects_fractional_atomic_numbers(mol, monkeypatch):
    host = PySCFHost(mol)
    monkeypatch.setattr(mol, "atom_charges", lambda: np.array([8.5, 1., 1.]))
    with pytest.raises(ValueError, match="integer vector"):
        host.structure()
