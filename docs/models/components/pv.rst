Pressure-Volume Component
=========================

The ``pv`` component adds a pressure-volume contribution :cite:p:`spooner2014compressed,zeller2025pressure`.

.. math::

   E_\mathrm{PV} = pV,

where ``pressure`` is supplied to ``new_pv`` and :math:`V` is the total cavity volume.
The component requires an updated cavity that provides ``total_volume``.

The PV surface weights are available and can be combiend with the :doc:`ρ-DROP isodensity cavity </cavities/isodensity>` for a fully self-consistent isodensity-PV.
Analytic nuclear gradients (with static and isodensity cavities) are also available.
