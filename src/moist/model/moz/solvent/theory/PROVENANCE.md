# DRISM correction

The Perkyns-Pettitt (1992) dielectric correction is a mathematical port of
RISMiCal's `drism.f` and `mathlib.f:diprot`, following the Fortran pysofft port
in `moist_dev`, revision `c15bd85a`,
`src/moist/model/rism/1d/equation.f90:compute_chi_drism`.

The RISMiCal attribution in that implementation is retained in `drism.f90`.
The dielectric amplitude, centre of absolute charge, dipole projection and
Gaussian/Bessel form factors follow that implementation. MOIST adds a complete
molecular frame, stable small-argument Bessel series and an analytical
`zeta/k^2` limit. Selection of ordinary RISM is explicit; `dielectric = 1`
continues to mean DRISM with that requested dielectric.

The native VV solve and the reconstructed solvent susceptibility use one
fixed theory configuration. All quantities are in atomic units; `smear` is
the Gaussian width in bohr in `exp(-smear^2*k^2/2)`.
