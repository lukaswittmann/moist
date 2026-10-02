"""SCF Hessians against finite differences of analytic nuclear gradients

Solvated comparison includes
- electronic density response
- motion of the isodensity DROP surface
- PCM charge response

Both sides use fully converged SCF states; the reference differentiates analytic
gradients only
"""

import os

import numpy as np
import pytest

try:
    from pyscf import gto, lib, scf
except ImportError as exc:
    if os.environ.get("MOIST_REQUIRE_PYSCF"):
        raise
    pytest.skip(f"pyscf is unavailable: {exc}", allow_module_level=True)

from .hessian import finite_difference_hessian, solvated_rhf_gradient
from .interface import CavityDROPIsodensityCallback, ModelComponentCPCM, ModelComponentCOSMO, ModelComponentPV, SolvationModel
from .pyscf import PySCFHost, solvated_rhf


pytestmark = pytest.mark.scf

#: Four-point central difference at a deliberately large step
#: - truncation falls as ``h**4``, round-off as ``eps max|g| / h``
#: - both sit far below the floor these references actually hit
FD_STEP = 1.0e-2
FD_STENCIL = 4

#: Agreement bound, ~7x the worst measured deviation for this system
#: - 7.3e-9 solvated, 7.1e-9 for PySCF's gas phase alone
#: - floor is PySCF's own gradient/Hessian consistency, not moist's
#: - see :mod:`moist.test_hessian_directional` docstring for the measurements and
#:   the knobs the floor does not respond to
FD_ATOL = 5.0e-8


@pytest.fixture
def water():
    # Asymmetric geometry: exercises off-diagonal atom and Cartesian blocks
    return gto.M(
        atom="O 0 0 -0.3893; H 0.7629 0 0.1947; H -0.7991 0.0953 0.2223",
        basis="sto-3g",
        unit="Angstrom",
        verbose=0,
    )


def rhf(mol):
    mean_field = scf.RHF(mol)
    mean_field.conv_tol = 1.0e-13
    mean_field.conv_tol_grad = 1.0e-10
    mean_field.kernel()
    assert mean_field.converged, "SCF must converge before differentiating"
    return mean_field


def test_gas_phase_scf_hessian_matches_fd(water):
    """Every relaxed RHF Hessian entry matches FD of analytic gradients"""
    with lib.with_omp_threads(1):
        mean_field = rhf(water)
        analytic = mean_field.Hessian().kernel()

        def gradient(positions):
            moved = water.set_geom_(positions, unit="Bohr", inplace=False)
            return rhf(moved).nuc_grad_method().kernel()

        numerical = finite_difference_hessian(
            gradient, water.atom_coords(), step=FD_STEP, stencil=FD_STENCIL
        )

    assert analytic.shape == (water.natm, water.natm, 3, 3)
    np.testing.assert_allclose(analytic, numerical, rtol=0, atol=FD_ATOL)
    np.testing.assert_allclose(
        analytic, analytic.transpose(1, 0, 3, 2), rtol=0, atol=1.0e-10
    )


@pytest.mark.parametrize("component_specs", [
    [(ModelComponentCPCM, 80.0)],
    [(ModelComponentCOSMO, 12.0)],
    [(ModelComponentCPCM, 4.0), (ModelComponentCOSMO, 12.0)],
])
def test_solvated_scf_hessian_matches_fd(water, component_specs):
    """Every relaxed DROP/CPCM RHF Hessian entry matches gradient FD"""
    cavity_options = dict(nleb=50, tolerance=1.0e-13)

    def model_factory(host):
        return SolvationModel(
            CavityDROPIsodensityCallback(host, **cavity_options),
            [component(epsilon) for component, epsilon in component_specs],
        )

    def gradient(positions):
        moved = water.set_geom_(positions, unit="Bohr", inplace=False)
        mean_field = solvated_rhf(
            moved, model_factory=model_factory, conv_tol_grad=1.0e-10
        )
        host = PySCFHost(moved)
        model = model_factory(host)
        return solvated_rhf_gradient(mean_field, model, host)

    with lib.with_omp_threads(1):
        mean_field = solvated_rhf(water, model_factory=model_factory, conv_tol_grad=1.0e-10)
        assert mean_field.converged
        hessian = mean_field.Hessian()
        analytic = hessian.kernel()
        numerical = finite_difference_hessian(
            gradient, water.atom_coords(), step=FD_STEP, stencil=FD_STENCIL
        )
        gas = rhf(water).Hessian().kernel()
        # PySCF's generic HVP has its own response path; must include solvent
        hop, diagonal = hessian.gen_hop()
        matrix = analytic.transpose(0, 2, 1, 3).reshape(3 * water.natm, -1)
        direction = np.arange(3 * water.natm, dtype=float)
        np.testing.assert_allclose(hop(direction), matrix.dot(direction), rtol=0, atol=1.0e-10)
        np.testing.assert_allclose(diagonal, matrix.diagonal(), rtol=0, atol=1.0e-10)

    np.testing.assert_allclose(analytic, numerical, rtol=0, atol=FD_ATOL)
    assert np.isfinite(numerical).all()
    np.testing.assert_allclose(
        numerical, numerical.transpose(1, 0, 3, 2), rtol=0, atol=FD_ATOL
    )
    assert np.max(np.abs(numerical - gas)) > 1.0e-4


def test_solvated_hessian_vacuum_limit():
    """Zero-solvent branch has no charge response and no division by zero"""
    mol = gto.M(atom="H 0 0 0; H 0.2 0.1 1.3", unit="Bohr", basis="sto-3g", verbose=0)
    with lib.with_omp_threads(1):
        vacuum = solvated_rhf(mol, 1.0, nleb=26, tolerance=1.0e-13)
        np.testing.assert_allclose(
            vacuum.Hessian().kernel(), rhf(mol).Hessian().kernel(), rtol=0, atol=1.0e-8
        )


@pytest.mark.parametrize("components", [
    [ModelComponentCPCM(80.0)],
    [ModelComponentCOSMO(12.0)],
    [ModelComponentCPCM(4.0), ModelComponentCOSMO(12.0)],
])
def test_solvent_density_response_matches_fd(water, components):
    """Mixed nuclear/density and pure density blocks match FD, pre-CPHF"""
    with lib.with_omp_threads(1):
        mean_field = rhf(water)
        # Nondefault level-set scale, arbitrary DM response
        host = PySCFHost(water, scale=2000.0)
        model = SolvationModel(
            CavityDROPIsodensityCallback(host, nleb=50, tolerance=1.0e-13),
            components,
        )
        coupling = host.coupling(mean_field.make_rdm1())
        parameters = coupling.second_order()
        analytic = model.evaluate(coupling=coupling).linearize()
        np.testing.assert_allclose(analytic.nuclear, analytic.nuclear.T, rtol=0, atol=1.0e-9)
        pp = analytic.density_response(np.eye(len(parameters.rows)))
        np.testing.assert_allclose(pp, pp.T, rtol=0, atol=1.0e-9)
        direction = np.random.default_rng(17).normal(size=len(parameters.rows))
        direction /= np.linalg.norm(direction)
        dm_direction = np.einsum("p,puv->uv", direction, parameters.dp[parameters.nr:])

        def gradient(offset):
            result = model.evaluate(coupling=host.coupling(mean_field.make_rdm1() + offset * dm_direction))
            fock = result.fock[parameters.rows, parameters.cols]
            fock = fock * np.where(parameters.rows == parameters.cols, 1, 2)
            return np.concatenate((result.gradient.T.ravel(), fock))

        step = 1.0e-5
        numerical = (gradient(step) - gradient(-step)) / (2 * step)
        expected = np.concatenate((
            np.einsum("pq,q->p", analytic.mixed, direction),
            analytic.density_response(direction),
        ))
    np.testing.assert_allclose(expected, numerical, rtol=1.0e-5, atol=2.0e-7)
    assert np.max(np.abs(expected[parameters.nr:])) > 1.0e-4


def test_model_linearization_capabilities_and_epoch(water):
    """Public model path rejects incomplete responses and stale evaluations"""
    with lib.with_omp_threads(1):
        mean_field = rhf(water)
        host = PySCFHost(water)
        model = SolvationModel(
            CavityDROPIsodensityCallback(host, nleb=50), [ModelComponentPV(1.0e-5)]
        )
        coupling = host.coupling(mean_field.make_rdm1())
        evaluation = model.evaluate(coupling=coupling)
        with pytest.raises(NotImplementedError, match="ModelComponentPV"):
            evaluation.linearize()
        model.evaluate(coupling=coupling)
        with pytest.raises(RuntimeError, match="superseded"):
            evaluation.linearize()


def test_composed_model_first_derivatives_are_additive(water):
    """Charge-weighted host fields are counted once across components"""
    with lib.with_omp_threads(1):
        mean_field = rhf(water)
        host = PySCFHost(water)
        coupling = host.coupling(mean_field.make_rdm1())
        components = [ModelComponentCPCM(4.0), ModelComponentCOSMO(12.0)]

        def evaluate(items):
            model = SolvationModel(CavityDROPIsodensityCallback(host, nleb=50), items)
            result = model.evaluate(coupling=coupling)
            return result.energy, result.fock, result.gradient

        combined = evaluate(components)
        separate = [evaluate([component]) for component in components]
        for actual, first, second in zip(combined, *separate):
            np.testing.assert_allclose(actual, first + second, rtol=0, atol=1.0e-10)


# ---------------------------------------------------------------------------
# PySCF Hessian class the solvent terms are folded into
#
# Solvent terms are method-agnostic
# - second derivatives of the model energy in nuclear coordinates and AO density
# - what must be right is the object they are added to
# - pure dispatch over two axes: KS or HF, density-fitted or not
# These pin dispatch and refusals; composed numbers: FD tests above
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("build,expected", [
    (lambda mol: scf.RHF(mol), "pyscf.hessian.rhf"),
    (lambda mol: scf.RKS(mol), "pyscf.hessian.rks"),
    (lambda mol: scf.RHF(mol).density_fit(), "pyscf.df.hessian.rhf"),
    (lambda mol: scf.RKS(mol).density_fit(), "pyscf.df.hessian.rks"),
])
def test_hessian_base_class_follows_the_reference(water, build, expected):
    """Each restricted reference maps to PySCF's Hessian class for it"""
    from .hessian import _pyscf_hessian_base

    assert _pyscf_hessian_base(build(water)).__module__ == expected


def test_restricted_ks_and_density_fitting_are_accepted(water):
    """Restricted KS and plain density fitting both carry a PySCF Hessian"""
    from .hessian import _check_restricted

    for mean_field in (scf.RKS(water), scf.RHF(water).density_fit()):
        mean_field.conv_tol = 1.0e-12
        mean_field.kernel()
        assert mean_field.converged
        _check_restricted(mean_field)


def test_semi_numerical_exchange_is_refused_by_name(water):
    """COSX is refused: no PySCF second derivative, cannot be treated as DF

    Treating it as plain DF returns the Hessian of a different Fock operator than
    the converged one, the failure of a permissive ``with_df`` check
    """
    from pyscf import sgx

    from .hessian import _check_restricted

    mean_field = sgx.sgx_fit(scf.RHF(water))
    mean_field.conv_tol = 1.0e-12
    mean_field.kernel()
    assert mean_field.converged
    with pytest.raises(NotImplementedError, match="no PySCF second-derivative"):
        _check_restricted(mean_field)


def test_unrestricted_reference_is_refused(water):
    """Unrestricted reference is refused: coupling needs one real symmetric density"""
    from .hessian import _check_restricted

    mean_field = scf.UHF(water)
    mean_field.conv_tol = 1.0e-12
    mean_field.kernel()
    with pytest.raises(NotImplementedError, match="closed-shell"):
        _check_restricted(mean_field)
