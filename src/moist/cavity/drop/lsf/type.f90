!> Abstract level-set function (LSF) base type for the DROP scheme
!>
!> The DROP cavity discretization is independent of *which* LSF defines the cavity surface.
!> This module provides the abstract base any concrete LSF must extend.
!>
!> A concrete LSF supplies:
!>   - its own constructor (per-LSF parameters; see e.g. SvdW's [[moist_cavity_drop_lsf_svdw_type]]%new)
!>   - implementations of the deferred procedures below
!>   - any LSF-specific public methods (e.g. SvdW's pou_*, normalized_* and f4_* derivatives used only by its FD tests)
!>
!> The cavity holds a `class(moist_cavity_drop_lsf_type), allocatable` model from which thread-local clones are made.
!> The [[lsf_thread_slot]] wrapper enables arrays of polymorphic LSF clones, since Fortran has no class(...), allocatable :: arr(:)
module moist_cavity_drop_lsf_base
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   implicit none
   private

   public :: moist_cavity_drop_lsf_type
   public :: lsf_thread_slot

   !> Abstract LSF base
   type, abstract :: moist_cavity_drop_lsf_type
      !> Number of atomic centers.
      integer :: ncenters = 0
      !> Molecular structure data.
      type(structure_type) :: mol
      !> Atomic radii (per-atom).
      real(wp), allocatable :: radii(:)
      !> Screening threshold below which the LSF contribution is treated
      !> as zero. Owned by the cavity (couples to the projection tolerance);
      !> the LSF concrete reads this value to size its internal screening
      !> caches whenever `update` runs. Direct users (tests calling the
      !> concrete without a cavity) may set it before calling `update`.
      real(wp) :: screening_threshold = 0.0_wp
   contains
      !> Bind molecular geometry. Concrete LSFs may override to refresh
      !> additional caches (e.g. SSD); the override should call this
      !> base implementation first via
      !> `call self%moist_cavity_drop_lsf_type%update(mol, radii)`.
      procedure :: update => lsf_base_update
      !> Cache per-point screening / state. Called once per evaluation
      !> point before any derivative method. No-op-safe.
      procedure(lsf_prepare_iface), deferred :: prepare
      !> Cache per-point state with a caller-provided candidate list.
      procedure(lsf_prepare_subset_iface), deferred :: prepare_subset
      !> Configure the highest spatial derivative order required.
      procedure(lsf_set_max_deriv_iface), deferred :: set_max_deriv
      !> Number of atoms currently active after `prepare`/`prepare_subset`.
      procedure(lsf_active_count_iface), deferred :: active_count
      !> User-space atom id of the i-th currently active atom.
      procedure(lsf_active_atom_iface), deferred :: active_atom
      !> LSF value only (lowest-cost path; used by marching cubes).
      procedure(lsf_f0_iface), deferred :: f0_screened
      !> Combined value/gradient/Hessian (any subset via optional args).
      procedure(lsf_f012_r_iface), deferred :: f012_r_screened
      !> Third spatial derivative (plus optionally lower-order outputs).
      procedure(lsf_f3_rrr_iface), deferred :: f3_rrr_screened
      !> Mixed third derivative: spatial Hessian w.r.t. nuclear positions.
      procedure(lsf_f3_rr_rA_iface), deferred :: f3_rr_rA_screened
      !> Per-atom radial offset beyond which this LSF's contribution
      !> falls below its internal screening threshold. Used by the cavity
      !> to size the per-atom cell-grid reach independently of the
      !> concrete LSF.
      procedure(lsf_neighbor_cutoff_iface), deferred :: neighbor_cutoff
   end type moist_cavity_drop_lsf_type

   !> Wrapper struct for arrays of polymorphic-allocatable LSF clones.
   !>
   !> Fortran does (afaik) not allow `class(...), allocatable :: arr(:)`, so per-thread LSF copies are
   !> stored as `type(lsf_thread_slot) :: arr(:)` with the polymorphic clone living inside the wrapper
   type :: lsf_thread_slot
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
   end type lsf_thread_slot

   abstract interface

      !> @param[inout] self   LSF instance
      !> @param[in]    point  Evaluation point (3,)
      subroutine lsf_prepare_iface(self, point)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
      end subroutine lsf_prepare_iface

      !> @param[inout] self              LSF instance
      !> @param[in]    point             Evaluation point (3,)
      !> @param[in]    candidate_indices Atom ids to consider (user-space)
      subroutine lsf_prepare_subset_iface(self, point, candidate_indices)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
         integer, intent(in) :: candidate_indices(:)
      end subroutine lsf_prepare_subset_iface

      !> @param[inout] self LSF instance
      !> @param[in]    n    Requested max derivative order
      subroutine lsf_set_max_deriv_iface(self, n)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         integer, intent(in) :: n
      end subroutine lsf_set_max_deriv_iface

      !> @param[in]  self LSF instance
      !> @returns         Number of currently active atoms
      pure function lsf_active_count_iface(self) result(n)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer :: n
      end function lsf_active_count_iface

      !> @param[in]  self LSF instance
      !> @param[in]  i    Active-list index (1 <= i <= active_count())
      !> @returns         User-space atom id of the i-th active atom
      pure function lsf_active_atom_iface(self, i) result(idx)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer, intent(in) :: i
         integer :: idx
      end function lsf_active_atom_iface

      !> @param[in]  self LSF instance
      !> @param[out] val  LSF value at the current evaluation point
      subroutine lsf_f0_iface(self, val)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out) :: val
      end subroutine lsf_f0_iface

      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Gradient w.r.t. spatial coords (optional)
      !> @param[out] lsf2_rr  Hessian w.r.t. spatial coords (optional)
      subroutine lsf_f012_r_iface(self, lsf0, lsf1_r, lsf2_rr)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:, :)
      end subroutine lsf_f012_r_iface

      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Spatial gradient (optional)
      !> @param[out] lsf2_rr  Spatial Hessian (optional)
      !> @param[out] lsf3_rrr Third spatial derivative tensor
      subroutine lsf_f3_rrr_iface(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:, :)
         real(wp), allocatable, intent(out) :: lsf3_rrr(:, :, :)
      end subroutine lsf_f3_rrr_iface

      !> @param[in]  self       LSF instance
      !> @param[out] lsf1_rA    LSF gradient w.r.t. nuclear positions (optional)
      !> @param[out] lsf2_r_rA  Mixed second derivative (optional)
      !> @param[out] lsf3_rr_rA Mixed third derivative (allocatable)
      subroutine lsf_f3_rr_rA_iface(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf1_rA(:, :)
         real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
         real(wp), allocatable, intent(out) :: lsf3_rr_rA(:, :, :, :)
      end subroutine lsf_f3_rr_rA_iface

      !> Conservative radial offset (measured from the atom surface)
      !> beyond which this LSF's contribution falls below its own
      !> screening threshold. Concrete LSFs supply the math; the cavity
      !> uses the result to size per-atom cell-grid reaches without
      !> knowing the LSF concrete.
      !>
      !> @param[in] self    LSF instance
      !> @param[in] radius  Atom radius (Bohr)
      !> @returns           Radial offset from atom surface (Bohr)
      pure function lsf_neighbor_cutoff_iface(self, radius) result(d)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(in) :: radius
         real(wp) :: d
      end function lsf_neighbor_cutoff_iface

   end interface

contains

   !> Default `update` implementation: copy geometry and radii into the common LSF state
   !>
   !> Concrete LSFs that need extra work (e.g. refresh screening caches) override this and call back
   !> to it via `call self%moist_cavity_drop_lsf_type%update(mol, radii)`
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Per-atom radii (size mol%nat)
   subroutine lsf_base_update(self, mol, radii)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      self%mol = mol
      self%ncenters = mol%nat
      if (allocated(self%radii)) deallocate (self%radii)
      self%radii = radii
   end subroutine lsf_base_update

end module moist_cavity_drop_lsf_base
