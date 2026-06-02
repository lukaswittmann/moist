!> Proxy module re-exporting the small-matrix and tensor linear-algebra utilities.
!>
!> Implementations live in the linalg/ sub-modules (decomp, outer, symmetrize,
!> geometry); this module aggregates their public interfaces so that consumers
!> can simply `use moist_math_linalg`.
module moist_math_linalg
   use moist_math_linalg_decomp, only: mat3x3_inv, eig_2x2_symmetric
   use moist_math_linalg_outer, only: outer_matrix, outer3, outer3_linear, outer4
   use moist_math_linalg_symmetrize, only: sym3_21, sym4_31, sym4_22, sym4_211
   use moist_math_linalg_geometry, only: cross_product, setup_tangent_frame, logaddexp
   implicit none
   public
end module moist_math_linalg
