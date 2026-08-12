MOIST-QM Coupling (PySCF)
=========================

:mod:`moist.pyscf` is the reference implementation for using MOIST together
with PySCF.

Interface
---------

The PySCF interface separates four roles:

``PySCFHost``
   Adapts a PySCF molecule and density to the level set, electrostatic, Fock, and nuclear-gradient operations MOIST needs.

``IsodensityDROPCavity(host, ...)``
   Constructs the cavity from the host's level set.
   The cavity owns its native state and configuration while retaining the host that supplies ``lsf(point, order)`` and the matching scale.

``host.coupling(dm)``
   Captures one density matrix and returns a :class:`~moist.pyscf.PySCFCoupling`.
   The adapter layer supplies the host data requested by the model components: electrostatics for PCM components such as   CPCM and COSMO, and Gaussian density moments for GOSTSHYP
   It also completes the Fock and gradient contributions.

``model.evaluate(...)``
   Rebuilds the cavity and completes the host exchange.
   It returns an immutable :class:`~moist.interface.Evaluation` containing the energy, potential, Fock contribution, and cavity.

Using CPCM and COSMO
--------------------

The :class:`~moist.interface.CPCM` model can be used as follows.

.. code-block:: python

   from moist import CPCM, IsodensityDROPCavity, SolvationModel
   from moist.pyscf import PySCFHost

   host = PySCFHost(mol)
   model = SolvationModel(
      cavity=IsodensityDROPCavity(host, nleb=194),
      components=[CPCM(80.0)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

With a density-independent cavity, a component such as :class:`~moist.interface.PV` requires no host data and
``model.evaluate(structure)`` is sufficient.
An isodensity cavity still needs ``host.coupling(dm)`` to get the density used by its level set callback, even when none of its model components needs an additional host channel.
Non-PySCF hosts that already own native-order arrays can use :class:`~moist.interface.ArrayCoupling` instead.

:class:`~moist.interface.COSMO` uses a dfferent constructor but with the same electrostatic host coupling as :class:`~moist.interface.CPCM`.
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
       CPCM,
       CFCDROPCavity,
       ISwiGCavity,
       SvdWDROPCavity,
       SolvationModel,
   )
   from moist.pyscf import PySCFHost

   host = PySCFHost(mol)

   # Smooth van der Waals (SvdW) DROP cavity
   cavity = SvdWDROPCavity(
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
       branch_weight_s=0.05,
       rho_grid_h=1.0,
       wleb_prune_level=0,
   )

   # COSMO Fine Cavity (CFC) DROP cavity
   cavity = CFCDROPCavity(
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
       branch_weight_s=0.05,
       rho_grid_h=1.0,
       wleb_prune_level=0,
   )

   # Improved Switching-Gaussian cavity
   cavity = ISwiGCavity(
       nleb=194,
       cut_a=0.0,
       cut_f=1.0e-10,
       debug=False,
       verbosity=0,
   )

   # Setup solvation model using the cavity
   model = SolvationModel(cavity=cavity, components=[CPCM(80.0)])
   result = model.evaluate(coupling=host.coupling(dm))

All constructor arguments are optional.
The shared ``nleb`` argument controls the Lebedev grid.
``debug`` enables native diagnostic output and ``verbosity`` sets its level.
The two DROP constructors also accept ``do_fine`` to request all optional surface properties and ``tolerance`` as their master numerical tolerance.
Their remaining shared controls are ``proj_maxiter`` and ``proj_level`` for surface projection, ``branch_weight_s`` for weighting competing projection branches, ``rho_grid_h`` for the grid-density kernel, and ``wleb_prune_level`` (0--6) for smooth pruning of negligible quadrature weights.
See :doc:`/cavities/drop` for details.

For :class:`~moist.interface.SvdWDROPCavity` (:doc:`/cavities/svdw`), ``blend_k`` and ``blend_1b``/``blend_2b``/``blend_3b`` control the smooth one-, two-, and three-body blend.
For :class:`~moist.interface.CFCDROPCavity` (:doc:`/cavities/cfc`), ``a1``, ``a2``, ``c``, and ``m`` are the CFC pseudo-density parameters; ``screen_k`` controls only neighbour screening.
For :class:`~moist.interface.ISwiGCavity` (:doc:`/cavities/iswig`), a positive ``cut_a`` selects area-based pruning; otherwise ``cut_f`` is the switching-function cutoff.

Isodensity ρ-DROP + CPCM
------------------------

CPCM needs the molecular electrostatic potential on the current surface.
The PySCF coupling computes it, obtains the induced surface charges from MOIST, and then supplies the charge-dependent surface and nuclear response terms.
The caller only constructs the objects and asks for one evaluation:

.. code-block:: python

   from pyscf import gto, scf
   from moist import CPCM, IsodensityDROPCavity, SolvationModel
   from moist.pyscf import PySCFHost

   mol = gto.M(atom="O 0 0 -0.7357; H 1.4418 0 0.3679; H -1.4418 0 0.3679",
               basis="def2-svp", unit="bohr")

   dm = scf.RHF(mol).run().make_rdm1()

   host = PySCFHost(mol)
   model = SolvationModel(
      cavity=IsodensityDROPCavity(host, nleb=194),
      components=[CPCM(80.0)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

   energy = result.energy
   fock = result.fock
   gradient = result.gradient
   surface = result.cavity
   grid_points = surface.grid_points

Components compose, and the host side is unchanged by which ones are present.
Adding a pressure term makes the cavity-shape response dominant rather than a small correction, because the PV energy depends on the density *only* through the volume of the surface:

.. code-block:: python

   from moist import PV

   pressure = 1.0e10 / 2.9421015697e13        # 10 GPa in E_h / a_0^3
   model = SolvationModel(
      cavity=IsodensityDROPCavity(host, nleb=194),
      components=[CPCM(80.0), PV(pressure)]
      )
   result = model.evaluate(coupling=host.coupling(dm))

With ``PV`` alone the surface charges vanish and the entire Fock contribution is the ``w_lsf`` contraction;
useful when checking a host implementation, since nothing else can mask an error in it.

For a self-consistent calculation, :func:`~moist.pyscf.solvated_rhf` wraps the whole cycle, rebuilding the cavity on every SCF iteration because the surface follows the density:

.. code-block:: python

   from moist.pyscf import solvated_rhf

   mf = solvated_rhf(mol, epsilon=80.0, nleb=194)
   print(mf.e_tot)

Isodensity ρ-DROP + GOSTSHYP
----------------------------

GOSTSHYP inverts the usual direction of the coupling.
The host computes AO-basis three-center integrals *for* MOIST—the Gaussian moments of the density about every grid point—and receives amplitudes to fold into its Fock matrix.
It is nevertheless an ordinary :class:`~moist.interface.SolvationComponent`:
the cavity owns ``nleb``, the component owns the pressure, and the PySCF coupling performs the host-specific moment exchange.

.. code-block:: python

   from pyscf import gto, scf
   from moist import CPCM, Gostshyp, IsodensityDROPCavity, SolvationModel
   from moist.gostshyp import GPA_TO_AU
   from moist.pyscf import PySCFHost

   mol = gto.M(atom="O 0 0 -0.7357; H 1.4418 0 0.3679; H -1.4418 0 0.3679",
               basis="def2-svp", unit="bohr")
   dm = scf.RHF(mol).run().make_rdm1()

   host = PySCFHost(mol)
   model = SolvationModel(
       cavity=IsodensityDROPCavity(host, nleb=194),
       components=[Gostshyp(50.0 * GPA_TO_AU)],
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
       cavity=IsodensityDROPCavity(host, nleb=194),
       components=[CPCM(80.0), Gostshyp(50.0 * GPA_TO_AU)],
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
   The widths from :meth:`~moist.interface.DROPCavity.get_gaussian` enter only the A-matrix.

``qefield`` and ``w_xyz`` are different quantities
   Both are charge-weighted potential gradients at the  grid points, but    ``qefield`` (gradient path) carries the *electronic* part only -- moist adds the nuclear part itself -- while ``w_xyz`` (potential path) must be the    *total*.
   Passing one where the other is expected is silent.

``w_xyz`` is required for a correct Fock matrix
   When the density changes the  grid points move and ``phi(r_i)`` moves with them.
   That route dominates the cavity response and moist cannot see it, so it must be supplied before the potential is read.  :class:`~moist.pyscf.PySCFCoupling` handles the ordering inside :meth:`~moist.interface.SolvationModel.evaluate`.

The callback returns the *unscaled* level set moist applies ``scale`` itself.
Scaling in the callback as well leaves the surface unchanged -- the zero level set is scale-invariant -- while making every adjoint wrong by a factor of ``scale``.

API
---

.. automodule:: moist.pyscf
   :members:
   :undoc-members:
   :show-inheritance:
