Continuum Components
====================

Components are reusable energy components that can be combined in a ``model_continuum_type`` and evaluated on a shared cavity.

Supported features
------------------

.. table::
   :class: matrix

   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | Component                           | Energies                    | Nuclear derivatives | Cavity response [#iso]_ | Complex variant [#complex]_ |
   |                                     +--------+--------------------+----------+----------+                         |                             |
   |                                     | Energy | Fock               | Gradient | Hessian  |                         |                             |
   +=====================================+========+====================+==========+==========+=========================+=============================+
   | :doc:`CPCM <pcm/cpcm>`              | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`COSMO <pcm/cosmo>`            | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`IEF-PCM <pcm/iefpcm>`         | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`SS(V)PE <pcm/ssvpe>`          | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`DPCM <pcm/dpcm>`              | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`LPB <pcm/lpb>`                | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |light|                     |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`SLIC <pcm/slic>`              | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |no|                        |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`NLE <pcm/nonlocal>`           | |yes|  | |yes|              | |yes|    | |yes|    | |yes|                   | |no|                        |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`PV <pressure/pv>`             | |yes|  | |light| [#pvfock]_ | |yes|    | |yes|    | |yes|                   | |no|                        |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+
   | :doc:`GOSTSHYP <pressure/gostshyp>` | |yes|  | |yes|              | |yes|    | |no|     | |yes|                   | |no|                        |
   +-------------------------------------+--------+--------------------+----------+----------+-------------------------+-----------------------------+

|yes| supported, |light| limited (see the footnotes), |no| not supported, |soon| in development.

.. [#iso] Needed for density-dependent cavities like :doc:`ρ-DROP </cavities/isodensity>`.
.. [#pvfock] PV itself does not depend on the density; only via the cavity response in case of an density-dependent cavity.
.. [#complex] Available in the generated kernels, not yet exposed through the APIs.

Energies
   The energy and its Fock-matrix contribution in the SCF.

Nuclear derivatives
   Analytic nuclear gradients and Hessians (on DROP cavities; all nine PCM
   components and PV expose directional HVPs and second-order surface-weight
   tangents in the channel API).

Cavity response
   Surface weights, the derivatives of the energy with respect to the surface quantities.
   They make the energy variational on a :doc:`ρ-DROP </cavities/isodensity>` cavity, whose surface follows the density.

Complex variant
   The same equations with a complex potential and permittivity :math:`\varepsilon(\omega)`,
   as needed for, e.g., frequency-dependent response, where the solvent polarization lags behind the field.
   

Components
----------

.. toctree::
   :hidden:

   pcm/index
   pressure/index

- :doc:`pcm/index`

  - :doc:`pcm/cpcm`
  - :doc:`pcm/cosmo`
  - :doc:`pcm/iefpcm`
  - :doc:`pcm/ssvpe`
  - :doc:`pcm/dpcm`
  - :doc:`pcm/lpb`
  - :doc:`pcm/slic`
  - :doc:`pcm/nonlocal`

- :doc:`pressure/index`

  - :doc:`pressure/pv`
  - :doc:`pressure/gostshyp`
