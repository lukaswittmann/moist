MOZ Closures
============

``moist_model_moz_closure_type`` defines the local closure relation, its first
and second derivatives and its free-energy integrand. The abstract interface, concrete
closures, construction and input checks live in
``src/moist/model/moz/closure/closure.f90``. Numerical formulas live in
``kernel.f90``, including the private PSE inverse used by its energy integrand.

Closures use physical quantities. The caller handles potential splitting
and grid layout; the MOZ thermodynamics module integrates the local energy.

Supported closures
------------------

Define the physical indirect correlation ``t = h - c`` and the closure
exponent ``x = -beta*lambda*u + t``. The factory accepts the names below
in uppercase or lowercase, with surrounding whitespace removed.

.. list-table::
   :header-rows: 1
   :widths: 15 55 30

   * - Name
     - Direct correlation
     - Local energy integrand
   * - HNC
     - ``c = exp(x) - 1 - t``
     - Available
   * - KH :cite:p:`kovalenko1999kh`
     - HNC for ``x <= 0``; ``c = x - t`` for ``x > 0``
     - Available
   * - PSE2, PSE3, PSE4 :cite:p:`kast2008pse`
     - HNC for ``x <= 0``; ``c = sum(x**k/k!, k=0..n) - 1 - t``
       for ``x > 0``
     - Available

API and conventions
-------------------

The inputs are the physical, unsplit pair potential ``u`` and physical
``t = h - c``. Arrays have shape ``(npts, npair)``; the caller flattens
spatial points. ``beta > 0`` is the inverse thermal energy in reciprocal
units of ``u``; ``lambda`` scales the potential and lies in ``[0, 1]``.

``new_closure(closure, name, error, nthreads)``
   Allocate the requested concrete closure. Unknown names leave the output
   unallocated and report an error. The PSE type exposes ``order``, which
   must be 2, 3 or 4. ``nthreads`` sizes the closure's OpenMP teams; models
   should pass ``ctx%get_num_threads()``. Absent or zero takes
   ``omp_get_max_threads()`` at each call, a negative count is an error, and
   sweeps below 65536 values run on one thread. ``closure%team_size(n)``
   returns the team size used for ``n`` values.

``closure%apply(u, t, beta, lambda, values, error)``
   Evaluate the physical direct correlation ``c`` with the closure kernels.
   All closures use block sweeps, with OpenMP parallelism for large arrays.

``closure%derivative(u, t, beta, lambda, values, error)``
   Evaluate ``dc/dt`` at fixed potential, temperature and coupling.

``closure%derivative2(u, t, beta, lambda, values, error)``
   Evaluate ``d2c/dt2`` in parallel. KH jumps at ``t* = 0``; ``c`` and
   ``dc/dt`` are continuous there, and PSE-n is continuous. Custom closures
   report it as unsupported unless they override it.

``closure%potential_derivative(u, t, beta, lambda, values, error)``
   Evaluate ``dc/du`` at fixed physical ``t``, temperature and coupling.

``closure%free_energy_derivatives(h, c, dh, dc, error)``
   Evaluate partial derivatives of the local integrand with respect to
   independent physical ``h`` and ``c``. Density and grid weights are excluded.

``closure%free_energy_derivative_fields(h, c, dh, dc, error)``
   The same derivatives over ``(npts, npair)`` arrays, with the input checks
   of ``free_energy_integrands``. Built-in closures use parallel sweeps;
   custom closures default to their checked scalar derivative operation.
   The UV response contracts the derivatives with its own weights.

``closure%local_energy_derivatives(h, c, dh, dc)``
   Pure subroutine returning both derivatives without input checks, the
   counterpart of ``local_integrand``. Optional for custom closures; the
   default returns NaNs, and the checked derivative operations report an error.
   Overriding this hook supplies the default checked scalar derivatives too.

``closure%apply_and_derivative(u, t, beta, lambda, values, derivs, error)``
   Evaluate ``c`` and ``dc/dt`` in one parallel sweep. HNC, KH and PSE share
   one exponential per point; ``c`` equals ``apply`` exactly. Custom closures
   default to ``apply`` followed by ``derivative``.

``closure%free_energy_integrands(h, c, phi, error)``
   Evaluate local integrands over ``(npts, npair)`` arrays in parallel, with
   the same input checks as the scalar method.

``closure%local_integrand(h, c)``
   Pure function returning the local integrand without input checks, for use
   inside parallel or pure reductions whose caller already validated ``h``
   and ``c``.

``closure%free_energy_integrand(h, c, phi, error)``
   Evaluate the dimensionless local integrand for the selected closure.
   Correlations must be finite and ``1 + h >= -1e-12``; the floor admits
   hard-core roundoff. The caller supplies density
   factors and integration weights.

The ``closure_*_c(u, t, beta, lambda, c)`` kernel subroutines evaluate a
one-dimensional block. They keep exponential arithmetic in the loop for
vectorization. Arrays must have matching sizes and separate output storage.
Pointwise elemental functions give derivatives and integrands; the
``closure_dc_sweep``, ``closure_d2c_sweep`` and ``closure_c_dc_sweep`` block
kernels (kind 0 HNC, 1 KH, 2-4 PSE order) back the parallel object methods.
Custom closures without
potential or integrand derivatives can still solve; response requests report
those unsupported operations explicitly. Energy-only extensions implement
``apply``, ``derivative``, ``free_energy_integrand`` and ``local_integrand``.
Integrand derivative support can be added through ``free_energy_derivatives``
or the pure ``local_energy_derivatives`` hook.
The exponential range reduction requires IEEE binary64,
round-to-nearest and compilation without unsafe math transformations such as
``-ffast-math``.

For HNC, ``Phi = h*h/2 - c - h*c/2``. KH includes the ``h*h/2``
term only for ``h <= 0``. PSE-n uses the HNC expression for
``h <= 0``; for positive ``h``, it subtracts
``z**(n+1)/(n+1)!``, where ``sum(z**k/k!, k=0..n) = 1 + h``.
The kernel recovers ``z`` by Newton iteration.

Input arrays and output storage must have matching shapes. Non-finite
inputs, invalid parameter domains and non-finite kernel results are reported
through ``error``. Outputs are defined only on success. The
exponential approximation clamps its argument to ``[-708, 710]``; overflow
is reported as a non-finite result.

The MOZ potentials return the Ng-split tables ``u_sr``, ``ur_lr`` and
``uk_lr`` (:doc:`potential`). HNC, KH and PSE take the short-range inputs
directly: ``apply(u_sr, t_s, ...)`` returns the short-range direct
correlation ``c_s``, because the split leaves the closure exponent
``-beta*lambda*u + t`` unchanged (``t_s = t - beta*lambda*u_lr``). The local
energy integrand needs the physical ``h`` and ``c = c_s - beta*lambda*u_lr``.

Tests
-----

``test/unit/test_moz_closure.f90`` checks factory construction, exponential
and polynomial reference values, coupling limits, branch continuity,
analytic first derivatives against finite differences, local energy
references, input errors and strided array sections. Block sweeps are compared
with direct kernel calls across incomplete blocks, single and multiple pair
columns, and strided storage with one and four OpenMP threads.

``test/unit/test_moz_closure_kernel.f90`` checks the kernel derivatives,
compares sweeps with independent mathematical references, and checks
branch boundaries, reversed/empty sections and exponential accuracy.
