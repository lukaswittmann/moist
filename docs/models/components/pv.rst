Pressure-Volume Component
=========================

The ``solvation_model_component_pv`` component adds a pressure-volume contribution :cite:p:`spooner2014compressed,zeller2025pressure`.

.. math::

   E_\mathrm{PV} = pV,

where ``pressure`` is supplied to ``new_component_pv`` and :math:`V` is the total cavity volume.
The component requires an updated cavity that provides ``total_volume``.

The PV surface weights are available and can be combined with the :doc:`ρ-DROP isodensity cavity </cavities/isodensity>` for a fully variational, self-consistent isodensity PV energy.
Analytic nuclear gradients are available with both static and isodensity cavities.
