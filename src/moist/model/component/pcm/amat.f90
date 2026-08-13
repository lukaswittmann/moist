!> Gaussian PCM interaction-matrix interface
!>
!> Matrix assembly and derivative contractions are implemented in
!> [[moist_model_component_pcm_amat_assembly]] and
!> [[moist_model_component_pcm_amat_adjoint]] submodules
!> Generated pair mathematics in [[moist_model_component_pcm_amat_kernel]]
module moist_model_component_pcm_amat
   use mctc_env, only: wp, error_type
   implicit none(type, external)
   private

   public :: assemble_pcm_amat
   public :: assemble_pcm_amat_with_gradient
   public :: pcm_amat_surface_weights
   public :: pcm_amat_nuclear_gradient

   !> Floor applied to squared separations
   real(wp), parameter :: r2_floor = 1.0e-200_wp

   interface

      !> Validate the common Gaussian surface arrays
      !>
      !> @param[in]  xi     Gaussian widths
      !> @param[in]  f      Gaussian switching factors
      !> @param[in]  xyz    Surface positions
      !> @param[out] error  Error handling
      module subroutine validate_pcm_surface(xi, f, xyz, error)
         implicit none(type, external)
         real(wp), intent(in) :: xi(:)
         real(wp), intent(in) :: f(:)
         real(wp), intent(in) :: xyz(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine validate_pcm_surface

      !> Compute per-point saturation bounds for the near/far test
      !>
      !> @param[in]  xi     Gaussian widths
      !> @param[out] bound  Per-point saturation bounds
      module pure subroutine saturation_bounds(xi, bound)
         implicit none(type, external)
         real(wp), intent(in) :: xi(:)
         real(wp), intent(out) :: bound(:)
      end subroutine saturation_bounds

      !> Assemble the Gaussian PCM interaction matrix
      !>
      !> @param[in]  xi     Gaussian widths
      !> @param[in]  f      Gaussian switching factors
      !> @param[in]  xyz    Surface positions
      !> @param[out] amat   Interaction matrix
      !> @param[out] error  Error handling
      module subroutine assemble_pcm_amat(xi, f, xyz, amat, error)
         implicit none(type, external)
         real(wp), intent(in) :: xi(:)
         real(wp), intent(in) :: f(:)
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(out) :: amat(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_pcm_amat

      !> Assemble the Gaussian PCM matrix and its nuclear derivative tensor
      !>
      !> @param[in]  xi       Gaussian widths
      !> @param[in]  f        Gaussian switching factors
      !> @param[in]  xyz      Surface positions
      !> @param[in]  xi1_rA   Nuclear derivatives of widths
      !> @param[in]  f1_rA    Nuclear derivatives of switching factors
      !> @param[in]  xyz1_rA  Nuclear derivatives of surface positions
      !> @param[out] amat     Interaction matrix
      !> @param[out] amat1_rA Nuclear derivative tensor
      !> @param[out] error    Error handling
      module subroutine assemble_pcm_amat_with_gradient(xi, f, xyz, xi1_rA, f1_rA, &
                                                        xyz1_rA, amat, amat1_rA, error)
         implicit none(type, external)
         real(wp), intent(in) :: xi(:)
         real(wp), intent(in) :: f(:)
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: xi1_rA(:, :, :)
         real(wp), intent(in) :: f1_rA(:, :, :)
         real(wp), intent(in) :: xyz1_rA(:, :, :, :)
         real(wp), intent(out) :: amat(:, :)
         real(wp), intent(out) :: amat1_rA(:, :, :, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_pcm_amat_with_gradient

      !> Contract a Gaussian PCM matrix derivative to surface-variable weights
      !>
      !> @param[in]  xi     Gaussian widths
      !> @param[in]  f      Gaussian switching factors
      !> @param[in]  xyz    Surface positions
      !> @param[in]  q1     Left contraction vector
      !> @param[in]  q2     Right contraction vector
      !> @param[out] w_xi   Width weights
      !> @param[out] w_f    Switching-factor weights
      !> @param[out] w_xyz  Position weights
      !> @param[out] error  Error handling
      module subroutine pcm_amat_surface_weights(xi, f, xyz, q1, q2, w_xi, w_f, &
                                                 w_xyz, error)
         implicit none(type, external)
         real(wp), intent(in) :: xi(:)
         real(wp), intent(in) :: f(:)
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: q1(:), q2(:)
         real(wp), intent(out) :: w_xi(:), w_f(:)
         real(wp), intent(out) :: w_xyz(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine pcm_amat_surface_weights

      !> Contract surface-variable weights with nuclear derivative arrays
      !>
      !> @param[in]  xi1_rA      Width derivatives
      !> @param[in]  f1_rA       Switching-factor derivatives
      !> @param[in]  xyz1_rA     Surface-position derivatives
      !> @param[in]  w_xi        Width weights
      !> @param[in]  w_f         Switching-factor weights
      !> @param[in]  w_xyz       Position weights
      !> @param[out] grad_rA     Nuclear gradient
      !> @param[out] error       Error handling
      module subroutine pcm_amat_nuclear_gradient(xi1_rA, f1_rA, xyz1_rA, w_xi, &
                                                   w_f, w_xyz, grad_rA, error)
         implicit none(type, external)
         real(wp), intent(in) :: xi1_rA(:, :, :), f1_rA(:, :, :), xyz1_rA(:, :, :, :)
         real(wp), intent(in) :: w_xi(:), w_f(:), w_xyz(:, :)
         real(wp), intent(out) :: grad_rA(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine pcm_amat_nuclear_gradient

   end interface

end module moist_model_component_pcm_amat
