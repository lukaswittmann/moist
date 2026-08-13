CFC-DROP Cavity
===============

The COSMO Fine Cavity (CFC) is a radii-based pseudo-density surface following :cite:t:`klamt2018cfc`, originally discretized via a marching tetrahedron algorithm.
A pseudo-density ``PD(r)`` is assembled from atomic and pairwise terms; the level set is ``-log PD(r)`` so the interior remains negative.

Optional settings:

``a1`` (real, default ``-15.0``)
   Atomic-term exponent.

``a2`` (real, default ``-9.0``)
   Pair-term exponent.

``c`` (real, default ``5.0``)
   Pair-term coupling constant.

``m`` (integer, default ``4``)
   Pair-term polynomial power.
   The kernel is generated for ``m = 4``; other values are currently ignored.

``screen_k`` (real, default ``3.0``)
   Sharpness of the conservative screening using the SvdW SSD.
   This affects computational cost, but not the CFC pseudo-density itself.
