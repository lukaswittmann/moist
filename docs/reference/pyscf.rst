PySCF integration
=================

Import ``moist.pyscf`` and attach MOIST to an existing calculation. Attachment
does not run SCF; use PySCF's usual ``kernel`` or ``run`` methods.

.. code-block:: python

   from pyscf import gto
   from moist import ModelComponentCPCM, DROPParameters
   from moist.pyscf import DROP, SvdW, Isodensity

   mol = gto.M(atom="O 0 0 0; H 0 -1.4 1.1; H 0 1.4 1.1",
               basis="def2-svp", unit="bohr")
   mf = mol.RHF().MOIST(
       cavity=DROP(lsf=SvdW(), parameters=DROPParameters(nleb=194)),
       components=[ModelComponentCPCM(80.0)],
   )
   mf.kernel()
   gradient = mf.nuc_grad_method().kernel()  # (natoms, 3)

The whole integration is one module, ``moist.pyscf``, which reads top to
bottom: the AO-basis integrals (:class:`~moist.pyscf.PySCFHost`), the GOSTSHYP
moment integrals (:class:`~moist.pyscf.GaussianMoments`), the driver that runs
the :doc:`coupling protocol <coupling>` (:class:`~moist.pyscf.PySCFSolvation`),
and the mean-field hook (:func:`~moist.pyscf.moist_for_scf`).

Cavities and components
-----------------------

Choose ``DROP(lsf=..., ...)`` or ``ISwiG(...)``; ``moist.pyscf`` re-exports
both, and importing it registers the ``.MOIST`` method on PySCF mean-field
objects.
``DROP`` requires a keyword-only ``lsf``: ``SvdW()``, ``CFC()`` or
``Isodensity()``. Put typed surface parameters on the LSF and discretization
parameters on ``DROP``. For example:

.. code-block:: python

   from moist import CFC, ISwiG, SvdWParameters, CFCParameters, ISwiGParameters

   grid = DROPParameters(nleb=194)
   cavity = DROP(lsf=SvdW(parameters=SvdWParameters(blend_k=5.5)), parameters=grid)
   cavity = DROP(lsf=CFC(parameters=CFCParameters(a1=-15.0)), parameters=grid)
   cavity = ISwiG(parameters=ISwiGParameters(nleb=194))

Configurations are immutable and reusable. See :doc:`/cavities/index` for
surface definitions and numerical controls.

For an isodensity surface, change the cavity arguments:

.. code-block:: python

   from moist import IsodensityParameters

   mf = mol.RHF().MOIST(
       cavity=DROP(
           lsf=Isodensity(parameters=IsodensityParameters(rho_iso=4e-4)),
           parameters=DROPParameters(nleb=194),
       ),
       components=[ModelComponentCPCM(80.0)],
   )
   mf.kernel()

The wrapper binds the current density before rebuilding the isodensity surface.
Static surfaces are reused across density changes. Energy and Fock include all
required cavity-response terms. The default ``rho_iso`` is ``1e-3``
electrons/bohr**3, as in C and Fortran. Radius models can be selected with
``radii=...`` on either cavity configuration.

CPCM, COSMO, PV and :doc:`GOSTSHYP </models/components/gostshyp>` components
share the selected cavity. For example:

.. code-block:: python

   from moist import ModelComponentGOSTSHYP, ModelComponentPV
   from moist.pyscf import GPA_TO_AU

   mf = mol.RHF().MOIST(
       cavity=DROP(lsf=Isodensity(), parameters=DROPParameters(nleb=194)),
       components=[ModelComponentCPCM(80.0),
                   ModelComponentGOSTSHYP(50.0 * GPA_TO_AU)],
   )
   # For a pV term, use ModelComponentPV(10.0 * GPA_TO_AU).

Configuration and results
-------------------------

``mf.with_moist`` holds the configuration and the latest MOIST ``result``
(energy and Fock contribution); ``e`` and ``v`` expose the two, and all three
are ``None`` before the first evaluation. Replace the cavity configuration to
change settings; this clears cached results:

.. code-block:: python

   from dataclasses import replace

   config = mf.with_moist.cavity
   mf.with_moist.set(cavity=replace(
       config, parameters=replace(config.parameters, nleb=302),
   ))

``mf.MOIST(..., parameters=ModelParameters(...))`` and
``mf.with_moist.set(parameters=...)`` configure model logging. ``set`` takes
only ``cavity``, ``components`` and ``parameters``; cavity settings live on the
cavity configuration, so ``set(nleb=302)`` raises ``TypeError``.

RHF, RKS, UHF and UKS are supported, including density fitting, nuclear gradients,
atom subsets, and energy/gradient scanners. Existing solver settings are
preserved. ``copy()`` creates independent MOIST state; ``reset(mol)`` clears it.
Geometry, basis and density changes invalidate the appropriate caches.

Attach MOIST once and combine terms in ``components``. ECPs, ROHF, existing
solvent attachments, Newton SCF, stability, TDDFT, Hessians, post-HF methods and
GPU execution are not supported by this wrapper.

Driving a model at a fixed density
----------------------------------

Outside an SCF, :class:`~moist.pyscf.PySCFSolvation` evaluates a model at a
given density matrix. It owns the model and one coupling and runs the six
steps of the :doc:`coupling` protocol itself, answering every request with
PySCF integrals on the cavity grid:

.. code-block:: python

   from moist.pyscf import PySCFSolvation

   solvation = PySCFSolvation(
       mol,
       DROP(lsf=Isodensity(), parameters=DROPParameters(nleb=194)),
       [ModelComponentCPCM(80.0)],
   )
   result = solvation.evaluate(dm)      # energy and Fock contribution
   gradient = solvation.gradient(dm)    # (natoms, 3), at fixed dm
   response = solvation.response        # the items MOIST handed back

``gradient_channels(dm)`` returns the same gradient split by route: the terms
MOIST owns (``model``), the basis-center derivative of the potential
(``potential``), of the level set (``density``, isodensity cavities only) and
of the Gaussian moments (``moments``, GOSTSHYP only). The host's own
contractions walk the response and dispatch on each item:
``potential_adjoint``, ``density`` and ``gaussian_amplitude``, each present
only when the model produces it; an item the driver does not know raises.

API
---

.. autofunction:: moist.pyscf.moist_for_scf

.. autoclass:: moist.pyscf.PySCFSolvation
   :members:

.. autoclass:: moist.pyscf.PySCFHost
   :members:

.. autoclass:: moist.pyscf.GaussianMoments
   :members:

``from moist.pyscf import DROP`` and ``from moist import DROP`` name the same
class. The shared configurations are documented once in :doc:`python`:
:py:class:`moist.DROP`, :py:class:`moist.ISwiG`, :py:class:`moist.SvdW`,
:py:class:`moist.CFC` and :py:class:`moist.Isodensity`.
