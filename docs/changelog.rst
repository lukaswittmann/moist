Release Notes
=============

1.2.1-alpha
-----------

Breaking changes
   - Solvent molecular data: IDs 181 to 187 (argon, bromine, hydrogen chloride, helium, neon,
     krypton, xenon), which have no descriptors, moved to 182 to 188 to make room for ammonia.
   - NLE component: arguments ``eps_solvent``, ``eps_inf`` and ``eps_solute`` renamed to
     ``epsilon``, ``epsilon_inf`` and ``epsilon_solute`` (Fortran, C, Python).

New features
   - Ammonia solvent (ID 181).
   - Solvent aliases ``propanonitrile`` and ``propanenitrile`` (propionitrile) and ``butanonitrile``
     (butyronitrile).

Fixes
   - Python: version lookup no longer fails (API reference build).
   - Phenol (ID 180): ``solvent`` charges and multipoles recomputed with its own permittivity
     (12.4).

Documentation
   - Restructured cavities, models and components with support tables and recommendations.
   - New Continuum Models, Pressure Models and Solvents sections.
   - Expanded PCM overview (model choice, physical inputs, numerical options) and GOSTSHYP derivation.
   - Updated contributing guidelines.
