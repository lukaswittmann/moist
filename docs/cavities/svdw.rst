SvdW-DROP Cavity
================

The smooth van der Waals surface :cite:p:`wittmann2026drop` is the default DROP level set.
It recovers solvent-excluded-surface-like features while avoiding the geometric singularities and crevices of a plain sphere union.
Neighbouring atomic contributions are combined through a smooth one-/two-/three-body blend.
The default parameters reproduce a probe radius of around 1.4 Å.

Optional settings:

``blend_k`` (real, default ``5.5``)
   Smoothing ``k`` in the ``exp(-k * d)`` kernel.
   Larger values give sharper features; smaller values smooth crevices more aggressively.

``blend_1b`` (real, default ``1.0``)
   One-body smoothing.

``blend_2b`` (real, default ``0.0``)
   Two-body smoothing.

``blend_3b`` (real, default ``3.0``)
   three-body smoothing.
