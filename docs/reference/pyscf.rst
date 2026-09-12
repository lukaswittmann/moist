MOIST-QM Coupling (PySCF)
=========================

:mod:`moist.pyscf` is the reference implementation for using MOIST together
with PySCF.

Interface
---------

The PySCF interface separates four roles:

``PySCFHost``
   Adapts a PySCF molecule and density to the density, electrostatic, Fock, and nuclear-gradient operations MOIST needs.

``CavityDROPIsodensity(host, ...)``
   Constructs the cavity from the host's density.
   The cavity owns its native state and configuration while retaining the host that supplies ``density(point, order)`` and the matching ``rho_iso`` and ``scale``.

``host.coupling(dm)``
   Captures one density matrix and returns a :class:`~moist.pyscf.PySCFCoupling`.
   The adapter layer supplies the host data requested by the model components: electrostatics for PCM components such as   CPCM and COSMO, and Gaussian density moments for GOSTSHYP
   It also completes the Fock and gradient contributions.

``model.evaluate(...)``
   Rebuilds the cavity and completes the host exchange.
   It returns an immutable :class:`~moist.interface.Evaluation` containing the energy, potential, Fock contribution, and cavity.

Using CPCM and COSMO
--------------------

The :class:`~moist.interface.ModelComponentCPCM` model can be used as follows.

.. code-block:: python

   from moist import CavityDROPIsodensity, ModelComponentCPCM, SolvationModel
   from moist.pyscf import PySCFHost

   host = PySCFHost(mol)
   model = SolvationModel(
      cavity=CavityDROPIsodensity(host, nleb=194),
      components=[ModelComponentCPCM(80.0)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

With a density-independent cavity, a component such as :class:`~moist.interface.ModelComponentPV` requires no host data and
``model.evaluate(structure)`` is sufficient.
An isodensity cavity still needs ``host.coupling(dm)`` to get the density used by its level set callback, even when none of its model components needs an additional host channel.
Non-PySCF hosts that already own native-order arrays can use :class:`~moist.interface.ArrayCoupling` instead.

:class:`~moist.interface.ModelComponentCOSMO` uses a dfferent constructor but with the same electrostatic host coupling as :class:`~moist.interface.ModelComponentCPCM`.
CPCM applies :math:`(\epsilon-1)/\epsilon`, whereas COSMO applies :math:`(\epsilon-1)/(\epsilon+1/2)`.

Both CPCM and COSMO accept the optional ``solver`` argument.
It may be ``"inversion"``, ``"lu"``, ``"cholesky"`` (the default), ``"iterative"``, or the corresponding :class:`~moist.interface.PCMSolver` value.


Choosing SvdW-DROP, CFC-DROP, or iSwiG
--------------------------------------

The radii-based cavities are density-independent and use the same model and PySCF coupling.
Only the cavity constructor changes.
For example, the following are three alternative CPCM models; keep the constructor for the surface you want:

.. code-block:: python

   from moist import (
       CavityDROPCFC,
       CavityDROPSvdW,
       CavityISwiG,
       ModelComponentCPCM,
       SolvationModel,
   )
   from moist.pyscf import PySCFHost

   host = PySCFHost(mol)

   # Smooth van der Waals (SvdW) DROP cavity
   cavity = CavityDROPSvdW(
       nleb=194,
       blend_k=5.5,
       blend_1b=1.0,
       blend_2b=0.0,
       blend_3b=3.0,
       debug=False,
       verbosity=0,
       do_fine=False,
       tolerance=1.0e-10,
       proj_maxiter=150,
       proj_level=3,
       branch_weight_s=0.0025,
       rho_grid_h=1.0,
       wleb_prune_level=0,
   )

   # COSMO Fine Cavity (CFC) DROP cavity
   cavity = CavityDROPCFC(
       nleb=194,
       a1=-15.0,
       a2=-9.0,
       c=5.0,
       m=4,
       screen_k=3.0,
       debug=False,
       verbosity=0,
       do_fine=False,
       tolerance=1.0e-10,
       proj_maxiter=150,
       proj_level=3,
       branch_weight_s=0.0025,
       rho_grid_h=1.0,
       wleb_prune_level=0,
   )

   # Improved Switching-Gaussian cavity
   cavity = CavityISwiG(
       nleb=194,
       cut_a=0.0,
       cut_f=1.0e-10,
       debug=False,
       verbosity=0,
   )

   # Setup solvation model using the cavity
   model = SolvationModel(cavity=cavity, components=[ModelComponentCPCM(80.0)])
   result = model.evaluate(coupling=host.coupling(dm))

All constructor arguments are optional.
The shared ``nleb`` argument controls the Lebedev grid.
``debug`` enables native diagnostic output and ``verbosity`` sets its level.
The two DROP constructors also accept ``do_fine`` to request all optional surface properties and ``tolerance`` as their master numerical tolerance.
Their remaining shared controls are ``proj_maxiter`` and ``proj_level`` for surface projection, ``branch_weight_s`` for weighting competing projection branches, ``rho_grid_h`` for the grid-density kernel, and ``wleb_prune_level`` (0--6) for smooth pruning of negligible quadrature weights.
See :doc:`/cavities/drop` for details.

For :class:`~moist.interface.CavityDROPSvdW` (:doc:`/cavities/svdw`), ``blend_k`` and ``blend_1b``/``blend_2b``/``blend_3b`` control the smooth one-, two-, and three-body blend.
For :class:`~moist.interface.CavityDROPCFC` (:doc:`/cavities/cfc`), ``a1``, ``a2``, ``c``, and ``m`` are the CFC pseudo-density parameters; ``screen_k`` controls only neighbour screening.
For :class:`~moist.interface.CavityISwiG` (:doc:`/cavities/iswig`), a positive ``cut_a`` selects area-based pruning; otherwise ``cut_f`` is the switching-function cutoff.

Isodensity ρ-DROP + CPCM
------------------------

CPCM needs the molecular electrostatic potential on the current surface.
The PySCF coupling computes it, obtains the induced surface charges from MOIST, and then supplies the charge-dependent surface and nuclear response terms.
The caller only constructs the objects and asks for one evaluation:

.. code-block:: python

   from pyscf import gto, scf
   from moist import CavityDROPIsodensity, ModelComponentCPCM, SolvationModel
   from moist.pyscf import PySCFHost

   mol = gto.M(atom="O 0 0 -0.7357; H 1.4418 0 0.3679; H -1.4418 0 0.3679",
               basis="def2-svp", unit="bohr")

   dm = scf.RHF(mol).run().make_rdm1()

   host = PySCFHost(mol)
   model = SolvationModel(
      cavity=CavityDROPIsodensity(host, nleb=194),
      components=[ModelComponentCPCM(80.0)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

   energy = result.energy
   fock = result.fock
   gradient = result.gradient
   surface = result.cavity
   grid_points = surface.xyz          # (3, ngrid) in bohr

Components compose, and the host side is unchanged by which ones are present.
Adding a pressure term makes the cavity-shape response dominant rather than a small correction, because the PV energy depends on the density *only* through the volume of the surface:

.. code-block:: python

   from moist import ModelComponentPV

   pressure = 1.0e10 / 2.9421015697e13        # 10 GPa in E_h / a_0^3
   model = SolvationModel(
      cavity=CavityDROPIsodensity(host, nleb=194),
      components=[ModelComponentCPCM(80.0), ModelComponentPV(pressure)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

With ``ModelComponentPV`` alone the surface charges vanish and the entire Fock contribution is the ``lsf`` contraction;
useful when checking a host implementation, since nothing else can mask an error in it.

For a self-consistent calculation, :func:`~moist.pyscf.solvated_rhf` wraps the whole cycle, rebuilding the cavity on every SCF iteration because the surface follows the density:

.. code-block:: python

   from moist.pyscf import solvated_rhf

   mf = solvated_rhf(mol, epsilon=80.0, nleb=194)
   print(mf.e_tot)

Analytic SCF Hessians
---------------------

The converged solvated RHF object provides the total analytic nuclear Hessian:

.. code-block:: python

   mol = gto.M(atom="O 0 0 -0.3893; H 0.7629 0 0.1947; H -0.7991 0.0953 0.2223",
               basis="sto-3g", unit="Angstrom")
   mf = solvated_rhf(mol, epsilon=80.0, nleb=50, tolerance=1.0e-13)
   hessian = mf.Hessian().kernel()   # (natoms, natoms, 3, 3), E_h / a_0^2

This includes the analytic AO derivatives, motion of the density-defined DROP
surface, PCM charge response, and solvent contributions to both the
coupled-perturbed SCF right-hand side and response kernel. The native DROP and
PCM second-derivative kernels are reused through host-parameter derivative
interfaces. No finite differences enter the analytic result.

There are two backends, selected by ``mf.Hessian(method)``. Both give the same
Hessian; they differ in what the host and moist exchange to get it.

``"directional"`` is the protocol the Fortran, C and Python layers share. moist
differentiates along a batch of directions and the host answers for its own
data along the same directions, so nothing quadratic in the number of density
variables is ever formed and the coupled-perturbed kernel is applied direction
by direction. A direction has a nuclear part moist sees and a host-private part
-- here a density-matrix direction -- that it never does, so one call serves
the RR, RP, PR and PP blocks alike.

Per block of directions moist asks the host for two things, through
:class:`~moist.interface.CouplingTangent`:

``level_set_tangent(first, dirs, xyz)``
   partial tangents of the scaled level set and its spatial derivatives to
   third order at the fixed surface points, ``(40, ngrid, nblk)``. Requested
   only when the level set is the host's, which for an isodensity cavity it
   is; moist's own level set reports no nuclear partials there, so a Hessian
   asked for without this is refused rather than served an unphysical zero.
``field_tangent(first, dirs, xyz, d_xyz)``
   the electronic potential and field at the *moving* points, in the
   ``qefield`` convention. Requested only by an electrostatic component whose
   potential the host supplies; moist forms the nuclear halves itself.

and hands back a :class:`~moist.interface.ResponseTangent` per direction: the
tangents of the surface charges and of the surface points, and the
gradient-path level-set adjoint weights with their tangents. The host completes
the Hessian columns and the Fock-matrix tangents from those, exactly as it
completes a gradient and a Fock matrix from a
:class:`~moist.interface.Response`. ``evaluation.hvp_coupled(dirs, tangent)``
is the entry point, ``moist.hessian.solvent_hvp`` the PySCF completion around
it, and ``moist_general_model_get_hvp_coupled()`` the same protocol in C with
the two callbacks as function pointers.

``"dense"`` is the reference backend. It requests a model-owned response from
``evaluation.linearize()``, containing ``nuclear`` (explicit RR derivatives),
``mixed`` (RP derivatives) and ``density_response(direction)`` (the PP action
in independent host density coordinates). PCM charge relaxation is included in
these quantities; electronic relaxation is handled by PySCF's CPHF solver.
Responses own their data and remain valid after the model is evaluated again.
Linearizing a superseded evaluation raises an error.

Components implement ``second_order(transaction)`` and cavities provide their
own surface derivatives. The PySCF Hessian adapter does not inspect component
types, dielectric constants, or PCM matrices. CPCM, COSMO, and sums of those
components currently implement the coupled second-order capability. Components
without that capability, including PV and GOSTSHYP, raise an explicit error.

Use a model factory to select components consistently for SCF and Hessians:

.. code-block:: python

   def model_factory(host):
       return SolvationModel(
           CavityDROPIsodensity(host, nleb=50, tolerance=1.0e-13),
           [ModelComponentCOSMO(80.0)],
       )

   mf = solvated_rhf(mol, model_factory=model_factory)
   hessian = mf.Hessian("directional").kernel()

The dense backend in :mod:`moist.second_order` stores dense surface second
partials for single-branch isodensity DROP. Its storage grows quadratically with
the number of independent AO density-matrix elements, so it is intended for small
systems; the directional backend has no such term. The charge-response part is
applied in factored form without assembling its PP block. Both backends support
conventional real, closed-shell, all-electron RHF and single-branch DROP
projections. PV has a second-order channel on the directional path only.
The native ``model.hessian()`` remains a fixed-host model derivative; the total
relaxed SCF Hessian is obtained from ``mf.Hessian()``.

``test_hessian_directional.py`` pins the directional protocol: the completed
nuclear columns and Fock tangents against central differences of the analytic
gradient and Fock matrix along a mixed nuclear/density direction, every channel
of the response tangent against differences of the quantity it claims to move,
the total SCF Hessian against the dense backend, and one configuration against
finite differences of reconverged SCF gradients. Run it with ``meson test -C
build pyscf_hessian_directional``.

``test_hessian.py`` compares all Hessian entries with central differences of
analytic gradients, reconverging SCF at every displaced geometry. It also checks
the gas-phase and vacuum limits and the mixed nuclear/density and pure density
response blocks. Run it with ``meson test -C build pyscf_hessian``.

Isodensity ρ-DROP + GOSTSHYP
----------------------------

GOSTSHYP inverts the usual direction of the coupling.
The host computes AO-basis three-center integrals *for* MOIST—the Gaussian moments of the density about every grid point—and receives amplitudes to fold into its Fock matrix.
It is nevertheless an ordinary :class:`~moist.interface.SolvationModelComponent`:
the cavity owns ``nleb``, the component owns the pressure, and the PySCF coupling performs the host-specific moment exchange.

.. code-block:: python

   from pyscf import gto, scf
   from moist import (
       CavityDROPIsodensity,
       ModelComponentCPCM,
       ModelComponentGOSTSHYP,
       SolvationModel,
   )
   from moist.gostshyp import GPA_TO_AU
   from moist.pyscf import PySCFHost

   mol = gto.M(atom="O 0 0 -0.7357; H 1.4418 0 0.3679; H -1.4418 0 0.3679",
               basis="def2-svp", unit="bohr")
   dm = scf.RHF(mol).run().make_rdm1()

   host = PySCFHost(mol)
   model = SolvationModel(
       cavity=CavityDROPIsodensity(host, nleb=194),
       components=[ModelComponentGOSTSHYP(50.0 * GPA_TO_AU)],
   )

   result = model.evaluate(coupling=host.coupling(dm))
   energy = result.energy
   fock = result.fock
   gradient = result.gradient
   surface = result.cavity

GOSTSHYP composes with other components without changing the host interface.
For example, an isodensity CPCM solvent under GOSTSHYP pressure uses the same
coupling:

.. code-block:: python

   model = SolvationModel(
       cavity=CavityDROPIsodensity(host, nleb=194),
       components=[ModelComponentCPCM(80.0), ModelComponentGOSTSHYP(50.0 * GPA_TO_AU)],
   )
   result = model.evaluate(coupling=host.coupling(dm))

The Gaussians sit *on* the grid points, so their moments are valid only for the surface that produced them.
:meth:`~moist.interface.SolvationModel.evaluate` therefore rebuilds the cavity and completes every requested host exchange in one transaction.

The returned :attr:`~moist.interface.Evaluation.gradient` is complete even though GOSTSHYP has no forward-mode path.  The coupling assembles the three routes a displacement takes:
the AO centers moving at a frozen surface, the level set following the density, and the surface response, which is MOIST's reverse-mode path.
As with CPCM, the density captured for the evaluation is immutable and restored before a lazy gradient.
The former all-in-one :class:`~moist.gostshyp.GostshypModel` and :class:`~moist.gostshyp.GostshypWall` interfaces remain as deprecated compatibility wrappers.

Conventions
-----------

These are easy to get wrong and are verified against finite differences in ``moist/test_pyscf.py`` and ``moist/test_gostshyp.py``:

``phi`` is the bare point potential
   moist builds the nuclear half of the surface-motion adjoint from an    unblurred ``Z_A (r_i - R_A)/r^3``.  A Gaussian-blurred potential would be inconsistent with it.
   The widths from :meth:`~moist.interface.CavityDROP.get_gaussian` enter only the A-matrix.

``qefield`` and ``w_xyz`` are different quantities
   Both are charge-weighted potential gradients at the  grid points, but    ``qefield`` (gradient path) carries the *electronic* part only -- moist adds the nuclear part itself -- while ``w_xyz`` (potential path) must be the    *total*.
   Passing one where the other is expected is silent.

``w_xyz`` is required for a correct Fock matrix
   When the density changes the  grid points move and ``phi(r_i)`` moves with them.
   That route dominates the cavity response and moist cannot see it, so it must be supplied before the potential is read.  :class:`~moist.pyscf.PySCFCoupling` handles the ordering inside :meth:`~moist.interface.SolvationModel.evaluate`.

The callback returns the bare density; moist forms ``S = scale * (rho_iso - rho)`` itself.
Subtracting the isovalue in the callback as well moves the surface, and scaling there as well leaves the surface unchanged -- the zero level set is scale-invariant -- while making every adjoint wrong by a factor of ``scale``.

API
---

.. automodule:: moist.pyscf
   :members:
   :undoc-members:
   :show-inheritance:
