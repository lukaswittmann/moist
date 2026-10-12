Release Notes
=============

Unreleased
----------

* Changed the PCM C options layout and solver IDs (dense = 1, iterative = 2).
  Callers must rebuild against the new header; the old, smaller options struct
  is rejected. The memory budget is a positive integer in GiB.

* Added generated-kernel CPCM/COSMO, IEF-PCM, SS(V)PE, DPCM (including the
  original discretization), LPB, SLIC and nonlocal PCM components. All use a
  shared dense or restarted-GMRES solver configuration.
* Extended the host potential protocol with spatial Hessians, the Gaussian
  mixed spatial/width derivative and ``w_dphi_dr`` adjoint weights. PySCF
  supplies and contracts these quantities for PCM energy, Fock and gradients.
* Added cavity normal and curvature requirements for PCM boundary models and
  preserved CPCM/COSMO regression values from the native A-matrix
  implementation on iSwiG and DROP cavities.
* Removed the public PCM A-matrix C helpers and the legacy conjugate-gradient
  solver choice.
