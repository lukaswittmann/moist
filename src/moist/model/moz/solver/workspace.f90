!> Owned scratch storage for sequential primal and transpose Krylov solves
module moist_model_moz_solver_workspace
   use mctc_env, only: wp
   implicit none(type, external)
   private
   public :: moist_model_moz_solver_workspace_type, exchange_solver_workspace

   !> Allocatable arrays keep their layout when ownership moves between solvers.
   !> Solver objects retain convergence settings, counters and validity flags.
   !> Recycled vectors/images are stored here but invalidated by the solver on transfer.
   type :: moist_model_moz_solver_workspace_type
      !> Shape of the real-space fields; zero before allocation
      integer :: ws_nflat = 0, ws_npts = 0, ws_npair = 0
      !> Shared scratch retains capacity for both requested restart dimensions
      logical :: shared = .false.
      !> Newton trial, correction and residual fields
      real(wp), allocatable :: ws_gamma_trial(:, :), ws_delta(:, :)
      real(wp), allocatable :: ws_f_curr(:, :), ws_f_trial(:, :)
      !> Allocated only when a finite-difference product is actually needed
      real(wp), allocatable :: ws_gamma_pert(:, :), ws_f_pert(:, :)
      !> Flat right-hand side, correction and Arnoldi work
      real(wp), allocatable :: ws_rhs_flat(:), ws_delta_flat(:), ws_r(:), ws_w(:)
      !> Optional deflation work
      real(wp), allocatable :: ws_z(:)
      !> Arnoldi basis (nflat, capacity + 1) and Hessenberg/Givens work
      !> Capacity follows the requested restart, not the number of steps actually taken
      real(wp), allocatable :: ws_v(:, :), ws_h(:, :)
      real(wp), allocatable :: ws_cs(:), ws_sn(:), ws_g(:), ws_y(:)
      real(wp), allocatable :: ws_proj(:), ws_proj2(:)
      !> Recycled basis/image storage; validity remains in the solver state
      real(wp), allocatable :: ws_u(:, :), ws_c(:, :)
      real(wp), allocatable :: ws_hbar(:, :), ws_bproj(:, :)
   end type moist_model_moz_solver_workspace_type

contains

   !> Swap complete allocations while preserving array layout and shared capacity
   !>
   !> @param[inout] first Owned scratch of one completed solver
   !> @param[inout] second Owned scratch of a distinct solver or solve context
   subroutine exchange_solver_workspace(first, second)
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: first, second
      type(moist_model_moz_solver_workspace_type), allocatable :: held

      call move_alloc(first, held)
      call move_alloc(second, first)
      call move_alloc(held, second)
      if (allocated(first)) first%shared = .true.
      if (allocated(second)) second%shared = .true.
   end subroutine exchange_solver_workspace
end module moist_model_moz_solver_workspace
