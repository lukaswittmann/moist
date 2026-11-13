!> Abstract seeder of the batched DROP projection
!>
!> A seeder proposes extra Newton starting points for an anchor whose local
!> solve failed or whose level searches every anchor
!>
!> - input: the anchor, its owning atom, the level set jet at the anchor, the
!>   admissible window, the radius from a verified local foot and the fallback
!>   search radius ([[seed_request_type]])
!> - output: `points(:, 1:n_points)`, which the caller appends to the open run
!>   of its batch with zero multipliers
!> - the count is known before anything is appended, so the caller can flush
!>   its batch first and never split an anchor's run
module moist_cavity_drop_projection_seed
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   implicit none(type, external)
   private

   public :: seed_request_type, projection_seeder_type

   !> Per-anchor seeding request

   ! TODO: good idea to introduce a new type here with drop parameters again?
   type :: seed_request_type
      !> Anchor slot, for traces
      integer :: index = 0
      !> Anchor position
      real(wp) :: anchor(3) = 0.0_wp
      !> Owning atom center
      real(wp) :: owner_center(3) = 0.0_wp
      !> Level set value at the anchor
      real(wp) :: lsf0_anchor = 0.0_wp
      !> Level set gradient at the anchor
      real(wp) :: grad_anchor(3) = 0.0_wp
      !> Admissible radius from a verified local foot, `sqrt(rho_1^2 + slack)`;
      !> non-positive when no such foot exists
      real(wp) :: rho_start = -1.0_wp
      !> Squared-distance window over the closest branch (Bohr^2)
      real(wp) :: rho2_slack = 0.0_wp
      !> Search radius around the anchor when neither `rho_start` nor a surface
      !> trace bounds the search (Bohr)
      real(wp) :: rho_fallback = 15.0_wp
      !> Seeds ordered by distance to the anchor, nearest first; seeders without an
      !> order of their own ignore it
      logical :: nearest_first = .false.
   end type seed_request_type

   !> Abstract seeder
   type, abstract :: projection_seeder_type
      !> Seeds of the last request (3, capacity)
      real(wp), allocatable :: points(:, :)
      !> Seeds produced by the last request
      integer :: n_points = 0
      !> Distance beyond the searched ball at which the last request's seeds still
      !> stand for surface inside it (Bohr); 0 when the seeder searches no ball
      real(wp) :: reach = 0.0_wp
   contains
      !> Propose seeds for one anchor into `points`
      procedure(seeder_seed_iface), deferred :: seed
      !> Append one seed to `points`, growing it geometrically
      procedure :: push => seeder_push
   end type projection_seeder_type

   abstract interface
      !> Propose seeds for one anchor
      !>
      !> @param[inout] self      Seeder
      !> @param[inout] lsf       Thread-local level set
      !> @param[in]    request   Anchor request
      !> @param[out]   error     Build-fatal failure (LSF, certificate budget)
      subroutine seeder_seed_iface(self, lsf, request, error)
         import :: projection_seeder_type, moist_cavity_drop_lsf_type, seed_request_type, error_type
         implicit none(type, external)
         class(projection_seeder_type), intent(inout), target :: self
         class(moist_cavity_drop_lsf_type), intent(inout), target :: lsf
         type(seed_request_type), intent(in) :: request
         type(error_type), allocatable, intent(out) :: error
      end subroutine seeder_seed_iface
   end interface

contains

   !> Append one seed to `points`
   !>
   !> @param[inout] self Seeder
   !> @param[in]    x    Seed
   subroutine seeder_push(self, x)
      class(projection_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: x(3)

      real(wp), allocatable :: tmp(:, :)
      integer :: cap

      if (.not. allocated(self%points)) allocate (self%points(3, 32))
      if (self%n_points + 1 > size(self%points, 2)) then
         cap = 2*size(self%points, 2)
         allocate (tmp(3, cap))
         tmp(:, 1:self%n_points) = self%points(:, 1:self%n_points)
         call move_alloc(tmp, self%points)
      end if
      self%n_points = self%n_points + 1
      self%points(:, self%n_points) = x
   end subroutine seeder_push

end module moist_cavity_drop_projection_seed
