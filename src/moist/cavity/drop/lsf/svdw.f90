!> Provides the Smooth van der Waals LSF function and its derivatives
!>
!> All implementations use screening-accelerated versions backed by SSD systems
!> to iterate over active nodes only (n_active << ncenters)
module moist_cavity_drop_lsf_svdw
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_svdw_ssd, only: moist_cavity_drop_lsf_svdw_ssd_type, ssd0
   use moist_math_linalg, only: sym3_21, outer3, outer3_linear, outer_matrix, &
                                outer4, sym4_31, sym4_22, sym4_211
   implicit none

   integer, parameter :: ndim = 3

   !> Smooth van der Waals LSF
   !>
   !> Concrete level-set function: takes its blending parameters directly
   !> through [[new]], owns the screening SSD system, refreshes geometry caches
   !> via [[update]], and primes per-point screening via [[prepare]]. Inherits
   !> common atom-LSF state (ncenters, mol, radii) from
   !> [[moist_cavity_drop_lsf_type]].
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_svdw_type

      !> SvdW parameters (blending k/1b/2b/3b + SSD screening threshold)
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      !> Screened SSD system, owned by the LSF
      !>
      !> Carries per-point active atom list and pre-computed derivative tensors.
      !> Lifecycle: allocated by [[new]], geometry-refreshed by [[update]], and
      !> recomputed at each evaluation point via [[prepare]]
      type(moist_cavity_drop_lsf_svdw_ssd_type) :: ssd_system

      !> Cached `[1..ncenters]` for prepare() to feed to ssd_system%compute.
      !> Built in [[update]] alongside the geometry refresh.
      integer, allocatable :: all_indices(:)
   contains
      !> Constructor: configure blending parameters and allocate the SSD system
      procedure, public :: new => lsf_new
      !> Bind molecular geometry and refresh SSD geometry caches
      procedure, public :: update => lsf_update
      !> Relabel cell-grid candidate ids into the SSD's spatially-sorted index
      !> space (so the screen loop drops its per-candidate orig->sorted gather)
      procedure, public :: remap_candidate_grid => lsf_svdw_remap_candidate_grid
      !> Point preparation: runs SSD screening at the evaluation point
      procedure, public :: prepare => lsf_prepare
      !> Point preparation variant using a caller-provided candidate atom list
      procedure, public :: prepare_subset => lsf_prepare_subset
      !> Configure the maximum SSD derivative order (re/de-allocates SoA arrays)
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      !> Number of atoms currently active after the latest prepare/prepare_subset
      procedure, public :: active_count => lsf_active_count
      !> User-space atom index of the i-th active atom (1 <= i <= active_count())
      procedure, public :: active_atom => lsf_active_atom
      !> Screened value (uses only base-type node fields, works with max_deriv=0)
      procedure, public :: f0_screened => lsf_f0_screened
      !> Screened combined value, gradient, and Hessian
      procedure, public :: f012_r_screened => lsf_f012_r_screened
      !> Screened third derivative w.r.t. spatial coordinates
      procedure, public :: f3_rrr_screened => lsf_f3_rrr_screened
      !> Screened mixed third derivative: spatial Hessian w.r.t. nuclear positions
      procedure, public :: f3_rr_rA_screened => lsf_f3_rr_rA_screened
      !> Screened partition of unity weights (spatial derivatives)
      procedure, public :: pou_f012_r_screened => lsf_pou_f012_r_screened
      !> Screened nuclear derivative of POU spatial gradient
      procedure, public :: pou_f2_r_rA_screened => lsf_pou_f2_r_rA_screened
      !> Screened normalized LSF: f0 / ||f1_r|| and its nuclear derivatives
      procedure, public :: normalized_f01_rA_screened => lsf_normalized_f01_rA_screened
      !> Screened pure nuclear Hessian: d^2S / dR_A dR_B
      procedure, public :: f2_rArB_screened => lsf_f2_rArB_screened
      !> Screened mixed third derivative: d^3S / dr dR_A dR_B
      procedure, public :: f3_r_rArB_screened => lsf_f3_r_rArB_screened
      !> Screened pure spatial fourth derivative: d^4S / dr^4
      procedure, public :: f4_rrrr_screened => lsf_f4_rrrr_screened
      !> Screened mixed fourth derivative: d^4S / dr^3 dR_A
      procedure, public :: f4_rrr_rA_screened => lsf_f4_rrr_rA_screened
      !> Screened mixed fourth derivative: d^4S / dr^2 dR_A dR_B
      procedure, public :: f4_rr_rArB_screened => lsf_f4_rr_rArB_screened
      !> Per-atom screening cutoff (radial offset from atom surface where
      !> the SvdW contribution falls below the inherited `screening_threshold`)
      procedure, public :: neighbor_cutoff => lsf_neighbor_cutoff
      !> Finalizer
      final :: finalize_lsf_svdw
   end type moist_cavity_drop_lsf_svdw_type

   public :: moist_cavity_drop_lsf_svdw_type

contains

!* ================================================================================= *!
!*                              LSF lifecycle methods                                *!
!* ================================================================================= *!

!> Configure LSF blending parameters
!>
!> @param[inout] self     LSF instance
!> @param[in]    blend_k  Blending sharpness k (optional)
!> @param[in]    blend_1b One-body weight (optional)
!> @param[in]    blend_2b Two-body weight (optional)
!> @param[in]    blend_3b Three-body weight (optional)
   subroutine lsf_new(self, blend_k, blend_1b, blend_2b, blend_3b)
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Blending sharpness k (optional override)
      real(wp), intent(in), optional :: blend_k
      !> One-body weight (optional override)
      real(wp), intent(in), optional :: blend_1b
      !> Two-body weight (optional override)
      real(wp), intent(in), optional :: blend_2b
      !> Three-body weight (optional override)
      real(wp), intent(in), optional :: blend_3b

      call self%param%new(blend_k=blend_k, blend_1b=blend_1b, &
                          blend_2b=blend_2b, blend_3b=blend_3b)
   end subroutine lsf_new

!> Bind molecular geometry and refresh SSD geometry caches
!>
!> @param[inout] self   LSF instance
!> @param[in]    mol    Molecular structure
!> @param[in]    radii  Per-atom radii (size mol%nat)
   subroutine lsf_update(self, mol, radii)
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      integer :: i

      integer :: prior_max_deriv

      self%mol = mol
      self%radii = radii
      self%ncenters = mol%nat

      !> (Re)initialise the screening SSD system from the LSF's current
      !> parameters and the inherited `screening_threshold` and preserve any max derivative order
      prior_max_deriv = self%ssd_system%max_deriv
      if (prior_max_deriv < 2) prior_max_deriv = 2
      call self%ssd_system%new( &
         k=self%param%blend_k, &
         threshold=self%screening_threshold, &
         max_deriv=prior_max_deriv)
      call self%ssd_system%update(mol%xyz, radii)

      ! Full-scan candidate list, in the SSD's spatially-sorted index space so
      ! the screen loop can index sorted_centers directly
      if (allocated(self%all_indices)) deallocate (self%all_indices)
      allocate (self%all_indices(mol%nat))
      do i = 1, mol%nat
         self%all_indices(i) = self%ssd_system%orig_to_sorted(i)
      end do
   end subroutine lsf_update

!> Relabel cell-grid candidate atom ids into the SSD's spatially-sorted index space
!>
!> @param[in]    self      SvdW LSF instance (must have been `update`d)
!> @param[inout] cell_nlat Flat cell-grid candidate atom-id list
   subroutine lsf_svdw_remap_candidate_grid(self, cell_nlat)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      integer, intent(inout) :: cell_nlat(:)
      integer :: i

      if (.not. allocated(self%ssd_system%orig_to_sorted)) return
      do i = 1, size(cell_nlat)
         cell_nlat(i) = self%ssd_system%orig_to_sorted(cell_nlat(i))
      end do
   end subroutine lsf_svdw_remap_candidate_grid

!> Run SSD screening at the evaluation point
!>
!> Populates the active-atom SoA arrays for the subsequent derivative call
!>
!> @param[inout] self  LSF instance
!> @param[in]    point Evaluation point (3,)
   subroutine lsf_prepare(self, point)
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

      call self%ssd_system%compute(point, self%all_indices)
   end subroutine lsf_prepare

!> Run SSD screening at the evaluation point using the caller-provided candidate atom list
!>
!> @param[inout] self              LSF instance
!> @param[in]    point             Evaluation point (3,)
!> @param[in]    candidate_indices Atom ids to screen, in the SSD's sorted-index
!>                                 space (see remap_candidate_grid)
   subroutine lsf_prepare_subset(self, point, candidate_indices)
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)

      call self%ssd_system%compute(point, candidate_indices)
   end subroutine lsf_prepare_subset

!> Configure the maximum spatial derivative order for SSD precomputation
!>
!> Allocates/deallocates `f3_rrr_arr` and `f4_rrrr_arr` storage so that the
!> next prepare() call fills only the tensors actually needed
!>
!> @param[inout] self LSF instance
!> @param[in]    n    Requested max derivative order (0..4)
   subroutine lsf_set_max_deriv(self, n)
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      integer, intent(in) :: n

      integer :: n_alloc

      self%ssd_system%max_deriv = n

      if (.not. allocated(self%ssd_system%k3f0_arr)) return
      n_alloc = size(self%ssd_system%k3f0_arr)

      ! f2_rr is saved only at max_deriv >= 3
      if (n >= 3) then
         if (.not. allocated(self%ssd_system%f2_rr_arr)) &
            allocate (self%ssd_system%f2_rr_arr(3, 3, n_alloc))
         if (.not. allocated(self%ssd_system%f3_rrr_arr)) &
            allocate (self%ssd_system%f3_rrr_arr(3, 3, 3, n_alloc))
      else
         if (allocated(self%ssd_system%f2_rr_arr)) deallocate (self%ssd_system%f2_rr_arr)
         if (allocated(self%ssd_system%f3_rrr_arr)) deallocate (self%ssd_system%f3_rrr_arr)
      end if

      if (n >= 4) then
         if (.not. allocated(self%ssd_system%f4_rrrr_arr)) &
            allocate (self%ssd_system%f4_rrrr_arr(3, 3, 3, 3, n_alloc))
      else
         if (allocated(self%ssd_system%f4_rrrr_arr)) deallocate (self%ssd_system%f4_rrrr_arr)
      end if
   end subroutine lsf_set_max_deriv

!> Number of atoms currently active after the latest prepare/prepare_subset.
   pure integer function lsf_active_count(self) result(n)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      n = self%ssd_system%n_active
   end function lsf_active_count

!> User-space atom index of the i-th currently active atom.
   pure integer function lsf_active_atom(self, i) result(idx)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      integer, intent(in) :: i
      idx = self%ssd_system%atom_indices(i)
   end function lsf_active_atom

!* ================================================================================= *!
!*                          Screened derivative methods                              *!
!* ================================================================================= *!

!> Screened value-only LSF evaluation
   subroutine lsf_f0_screened(self, val)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(out) :: val

      real(wp) :: k, s_1, s_2, s_3, Z
      integer :: n_active

      val = 0.0_wp
      n_active = self%ssd_system%n_active
      if (n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      call compute_lsf_z0_screened(self%ssd_system, k, s_1, s_2, s_3, Z)

      val = -log(Z)/k
   end subroutine lsf_f0_screened

!> Compute value, gradient, and Hessian of LSF function (screened)
   subroutine lsf_f012_r_screened(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: invZ, invZ2

      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%ssd_system%n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      if (present(lsf2_rr)) then
         call compute_lsf_z012_rr_screened(self%ssd_system, k, s_1, s_2, s_3, Z, gradZ, hessZ)
      else
         call compute_lsf_z012_rr_screened(self%ssd_system, k, s_1, s_2, s_3, Z, gradZ)
      end if

      if (present(lsf0)) lsf0 = -1.0_wp/k*log(Z)
      invZ = 1.0_wp/max(Z, 1.0e-100_wp)
      if (present(lsf1_r)) lsf1_r = -gradZ*invZ/k
      if (present(lsf2_rr)) then
         invZ2 = invZ*invZ
         lsf2_rr = (outer_matrix(gradZ, gradZ)*invZ2 - hessZ*invZ)/k
      end if
   end subroutine lsf_f012_r_screened

!> Compute third derivative of LSF function w.r.t. spatial coordinates (screened)
   subroutine lsf_f3_rrr_screened(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), allocatable, intent(out) :: lsf3_rrr(:, :, :)

      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim), thirdZ(ndim, ndim, ndim)
      real(wp) :: invZ, invZ2, invZ3

      allocate (lsf3_rrr(ndim, ndim, ndim))
      lsf3_rrr = 0.0_wp
      if (self%ssd_system%n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      call compute_lsf_z0123_rrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                          Z, gradZ, hessZ, thirdZ)

      if (present(lsf0)) lsf0 = -1.0_wp/k*log(Z)
      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ
      if (present(lsf1_r)) lsf1_r = -gradZ*invZ/k
      if (present(lsf2_rr)) then
         lsf2_rr = (outer_matrix(gradZ, gradZ)*invZ2 - hessZ*invZ)/k
      end if
      lsf3_rrr = (-thirdZ*invZ + sym3_21(hessZ, gradZ)*invZ2 &
                  - 2.0_wp*outer3(gradZ)*invZ3)/k
   end subroutine lsf_f3_rrr_screened

!> Compute mixed third derivative: spatial Hessian w.r.t. nuclear positions (screened)
   subroutine lsf_f3_rr_rA_screened(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: lsf3_rr_rA(:, :, :, :)

      real(wp), allocatable :: dS(:, :)
      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp), allocatable :: dH_Z(:, :, :, :)
      real(wp), allocatable :: grad_s_nuc(:, :, :)
      real(wp) :: invZ, invZ2, invZ3, dZ, dS_axis
      real(wp) :: dgradZ(ndim), dhess(ndim, ndim)
      integer :: atom, axis, i

      allocate (lsf3_rr_rA(ndim, ndim, ndim, self%ncenters))
      lsf3_rr_rA = 0.0_wp
      if (self%ncenters == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      call compute_lsf_z012_rr_screened(self%ssd_system, k, s_1, s_2, s_3, Z, gradZ, hessZ)
      call compute_lsf_z3_rr_rA_screened(self%ssd_system, k, s_1, s_2, s_3, dH_Z)

      call compute_lsf_f1_rA_screened(self%ssd_system, k, s_1, s_2, s_3, Z, dS)
      if (present(lsf1_rA)) then
         lsf1_rA = 0.0_wp
         do i = 1, self%ssd_system%n_active
            lsf1_rA(:, self%ssd_system%atom_indices(i)) = dS(:, i)
         end do
      end if

      call compute_lsf_f2_r_rA_screened(self%ssd_system, k, s_1, s_2, s_3, Z, grad_s_nuc)
      if (present(lsf2_r_rA)) then
         lsf2_r_rA = 0.0_wp
         do i = 1, self%ssd_system%n_active
            lsf2_r_rA(:, :, self%ssd_system%atom_indices(i)) = grad_s_nuc(:, :, i)
         end do
      end if

      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ

      do i = 1, self%ssd_system%n_active
         atom = self%ssd_system%atom_indices(i)
         do axis = 1, ndim
            dS_axis = dS(axis, i)
            dZ = -k*Z*dS_axis
            dgradZ = -k*Z*grad_s_nuc(:, axis, i) + gradZ*(dZ*invZ)
            dhess = dH_Z(:, :, axis, i)
            lsf3_rr_rA(:, :, axis, atom) = (outer_matrix(dgradZ, gradZ) &
                                            + outer_matrix(gradZ, dgradZ))*invZ2
            lsf3_rr_rA(:, :, axis, atom) = lsf3_rr_rA(:, :, axis, atom) &
                                           + outer_matrix(gradZ, gradZ)*(-2.0_wp*dZ*invZ3)
            lsf3_rr_rA(:, :, axis, atom) = lsf3_rr_rA(:, :, axis, atom) - dhess*invZ
            lsf3_rr_rA(:, :, axis, atom) = lsf3_rr_rA(:, :, axis, atom) + hessZ*(dZ*invZ2)
            lsf3_rr_rA(:, :, axis, atom) = lsf3_rr_rA(:, :, axis, atom)/k
         end do
      end do
   end subroutine lsf_f3_rr_rA_screened

!* ================================================================================= *!
!*                     Fourth-order primitives (f4_*_screened)                       *!
!* ================================================================================= *!

!> Compute pure spatial fourth derivative of LSF function (screened)
!> via Z-level quotient rule extension of the third-derivative formula
   subroutine lsf_f4_rrrr_screened(self, lsf4_rrrr)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: lsf4_rrrr(:, :, :, :)

      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: thirdZ(ndim, ndim, ndim), fourthZ(ndim, ndim, ndim, ndim)
      real(wp) :: invZ, invZ2, invZ3, invZ4

      allocate (lsf4_rrrr(ndim, ndim, ndim, ndim))
      lsf4_rrrr = 0.0_wp
      if (self%ssd_system%n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      call compute_lsf_z01234_rrrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                            Z, gradZ, hessZ, thirdZ, fourthZ)

      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ
      invZ4 = invZ2*invZ2

      lsf4_rrrr = (-fourthZ*invZ &
                   + sym4_31(gradZ, thirdZ)*invZ2 &
                   + sym4_22(hessZ, hessZ)*invZ2 &
                   - 2.0_wp*sym4_211(gradZ, hessZ)*invZ3 &
                   + 6.0_wp*outer4(gradZ)*invZ4)/k
   end subroutine lsf_f4_rrrr_screened

!> Compute mixed fourth derivative d^4S / (dr_j dr_k dr_l dR_{A,s_1}) (screened)
   subroutine lsf_f4_rrr_rA_screened(self, lsf4_rrr_rA)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: lsf4_rrr_rA(:, :, :, :, :)

      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: thirdZ(ndim, ndim, ndim), fourthZ(ndim, ndim, ndim, ndim)
      real(wp) :: invZ, invZ2, invZ3, invZ4

      ! Per-active-atom Z-tensor R-derivatives (s_1 is R-axis):
      ! dZ_rA(s_1) = dZ/ dR_{A,s_1}
      ! dgradZ_rA(j, s_1) = d(Z_j)/ dR_{A,s_1}
      ! dhessZ_rA(j, k, s_1) = d(Z_jk)/ dR_{A,s_1}
      ! dthirdZ_rA(j, k, l, s_1) = d(Z_jkl)/ dR_{A,s_1}
      real(wp), allocatable :: dZ_rA(:, :)
      real(wp), allocatable :: dgradZ_rA(:, :, :)
      real(wp), allocatable :: dhessZ_rA(:, :, :, :)
      real(wp), allocatable :: dthirdZ_rA(:, :, :, :, :)

      real(wp) :: g(ndim), h(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: dZ_val, dg(ndim), dh(ndim, ndim), dt3(ndim, ndim, ndim)
      real(wp) :: blk(ndim, ndim, ndim)
      integer :: n_active, iA, atom, axis

      allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, self%ncenters))
      lsf4_rrr_rA = 0.0_wp
      if (self%ncenters == 0) return
      n_active = self%ssd_system%n_active
      if (n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b

      call compute_lsf_z01234_rrrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                            Z, gradZ, hessZ, thirdZ, fourthZ)

      call compute_lsf_zR_derivs_rrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                              dZ_rA, dgradZ_rA, dhessZ_rA, dthirdZ_rA)

      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ
      invZ4 = invZ2*invZ2

      g = gradZ
      h = hessZ
      t3 = thirdZ

      do iA = 1, n_active
         atom = self%ssd_system%atom_indices(iA)
         do axis = 1, ndim
            dZ_val = dZ_rA(axis, iA)
            dg = dgradZ_rA(:, axis, iA)
            dh = dhessZ_rA(:, :, axis, iA)
            dt3 = dthirdZ_rA(:, :, :, axis, iA)

            blk = -dt3*invZ + t3*(dZ_val*invZ2)
            blk = blk + sym3_21(dh, g)*invZ2
            blk = blk + sym3_21(h, dg)*invZ2
            blk = blk - sym3_21(h, g)*(2.0_wp*dZ_val*invZ3)
            blk = blk - 2.0_wp*outer3_linear(g, dg)*invZ3
            blk = blk + 6.0_wp*outer3(g)*(dZ_val*invZ4)

            lsf4_rrr_rA(:, :, :, axis, atom) = blk/k
         end do
      end do

      deallocate (dZ_rA, dgradZ_rA, dhessZ_rA, dthirdZ_rA)
   end subroutine lsf_f4_rrr_rA_screened

!> Compute mixed fourth derivative d^4S / (dr_j dr_k dR_{A,s_1} dR_{B,s_2}) (screened)
   subroutine lsf_f4_rr_rArB_screened(self, lsf1_rA, lsf2_r_rA, lsf4_rr_rArB)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(in) :: lsf1_rA(:, :)
      real(wp), intent(in) :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: lsf4_rr_rArB(:, :, :, :, :, :)

      real(wp) :: k, s_1, s_2, s_3, inv_k
      real(wp) :: Z, invZ, invZ2, invZ3, invZ4
      real(wp) :: gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: thirdZ(ndim, ndim, ndim), fourthZ(ndim, ndim, ndim, ndim)

      ! Z R-derivatives (3D output): atoms x axes
      real(wp), allocatable :: dZ_rA(:, :)               ! (s_1, iA)
      real(wp), allocatable :: dgradZ_rA(:, :, :)        ! (j, s_1, iA)
      real(wp), allocatable :: dhessZ_rA(:, :, :, :)     ! (j, k, s_1, iA)
      real(wp), allocatable :: dthirdZ_rA(:, :, :, :, :) ! (j, k, l, s_1, iA)

      ! Z RR-derivatives (pair): d^2 Z / (dR_A dR_B), etc.
      real(wp), allocatable :: d2Z_rArB(:, :, :, :)       ! (s_1, iA, s_2, iB)
      real(wp), allocatable :: d2gradZ_rArB(:, :, :, :, :) ! (j, s_1, iA, s_2, iB)
      real(wp), allocatable :: d2hessZ_rArB(:, :, :, :, :, :) ! (j, k, s_1, iA, s_2, iB)

      real(wp) :: ZA, ZB, ZAB
      real(wp) :: gA(ndim), gB(ndim)
      real(wp) :: hjkA(ndim, ndim), hjkB(ndim, ndim)
      real(wp) :: kappa(ndim, ndim)
      real(wp) :: g(ndim), h(ndim, ndim)
      integer :: n_active, iA, iB, axa, axb, j, kk

      ! Note: lsf1_rA and lsf2_r_rA are inputs kept for API symmetry with
      ! lsf3_r_rArB. They are not used internally because we re-derive everything
      ! at the Z level; declared inputs satisfy the calling convention.
      if (size(lsf1_rA) < 0) return
      if (size(lsf2_r_rA) < 0) return

      n_active = self%ssd_system%n_active
      allocate (lsf4_rr_rArB(ndim, ndim, ndim, n_active, ndim, n_active))
      lsf4_rr_rArB = 0.0_wp
      if (n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b
      inv_k = 1.0_wp/k

      call compute_lsf_z01234_rrrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                            Z, gradZ, hessZ, thirdZ, fourthZ)
      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ
      invZ4 = invZ2*invZ2

      call compute_lsf_zR_derivs_rrr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                              dZ_rA, dgradZ_rA, dhessZ_rA, dthirdZ_rA)
      call compute_lsf_zRR_derivs_rr_screened(self%ssd_system, k, s_1, s_2, s_3, &
                                              d2Z_rArB, d2gradZ_rArB, d2hessZ_rArB)

      g = gradZ; h = hessZ

      do iB = 1, n_active
         do axb = 1, ndim
            ZB = dZ_rA(axb, iB)
            gB = dgradZ_rA(:, axb, iB)
            hjkB = dhessZ_rA(:, :, axb, iB)

            do iA = 1, n_active
               do axa = 1, ndim
                  ZA = dZ_rA(axa, iA)
                  gA = dgradZ_rA(:, axa, iA)
                  hjkA = dhessZ_rA(:, :, axa, iA)
                  ZAB = d2Z_rArB(axa, iA, axb, iB)

                  do kk = 1, ndim
                     do j = 1, ndim
                        kappa(j, kk) = d2hessZ_rArB(j, kk, axa, iA, axb, iB)*invZ
                        kappa(j, kk) = kappa(j, kk) &
                                       - (hjkA(j, kk)*ZB + hjkB(j, kk)*ZA &
                                          + d2gradZ_rArB(j, axa, iA, axb, iB)*g(kk) &
                                          + d2gradZ_rArB(kk, axa, iA, axb, iB)*g(j))*invZ2
                        kappa(j, kk) = kappa(j, kk) &
                                       - (h(j, kk)*ZAB &
                                          + gA(j)*gB(kk) &
                                          + gB(j)*gA(kk))*invZ2
                        kappa(j, kk) = kappa(j, kk) &
                                       + 2.0_wp*( &
                                       g(j)*g(kk)*ZAB &
                                       + ZA*ZB*h(j, kk) &
                                       + g(j)*ZA*gB(kk) &
                                       + g(j)*ZB*gA(kk) &
                                       + g(kk)*ZA*gB(j) &
                                       + g(kk)*ZB*gA(j) &
                                       )*invZ3
                        kappa(j, kk) = kappa(j, kk) &
                                       - 6.0_wp*g(j)*g(kk)*ZA*ZB*invZ4
                     end do
                  end do

                  lsf4_rr_rArB(:, :, axa, iA, axb, iB) = -kappa*inv_k
               end do
            end do
         end do
      end do

      deallocate (dZ_rA, dgradZ_rA, dhessZ_rA, dthirdZ_rA)
      deallocate (d2Z_rArB, d2gradZ_rArB, d2hessZ_rArB)
   end subroutine lsf_f4_rr_rArB_screened

!* ================================================================================= *!
!*                     Unified Z routines (power-sum based)                          *!
!* ================================================================================= *!

!> Compute Z = s_1*e1 + s_2*e2 + s_3*e3 value only, using O(N) power sums
!> The shared SSD cache stores b_i = exp(-(k/3) d_i), from which: e1_i = b_i^3,  e2_i = b_i * sqrt(b_i),  e3_i = b_i
   subroutine compute_lsf_z0_screened(ssd_system, k, s_1, s_2, s_3, Z)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), intent(out) :: Z

      integer :: i, n_active
      real(wp) :: base, base_sq, e1, e2
      real(wp) :: sum_e1, sum_e2, sum_e3, sum_e3_sq

      Z = 0.0_wp
      n_active = ssd_system%n_active
      if (n_active == 0) return

      sum_e1 = 0.0_wp
      sum_e2 = 0.0_wp
      sum_e3 = 0.0_wp
      sum_e3_sq = 0.0_wp
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         sum_e1 = sum_e1 + e1
         sum_e2 = sum_e2 + e2
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq
      end do

      Z = s_1*sum_e1
      if (n_active >= 2) Z = Z + s_2*(sum_e2*sum_e2 - sum_e1)*0.5_wp
      if (n_active >= 3) then
         Z = Z + s_3*(sum_e3**3 - 3.0_wp*sum_e3*sum_e3_sq + 2.0_wp*sum_e1)/6.0_wp
      end if
   end subroutine compute_lsf_z0_screened

!> Compute Z, gradZ, and (optionally) hessZ using O(N) power sums
   subroutine compute_lsf_z012_rr_screened(ssd_system, k, s_1, s_2, s_3, Z, gradZ, hessZ)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), intent(out) :: Z, gradZ(ndim)
      real(wp), intent(out), optional :: hessZ(ndim, ndim)

      integer :: i, j, ii, n_active
      real(wp) :: lambda1, lambda2, lambda3, lambda3_sq
      real(wp) :: lambda1_sq, lambda2_sq, lambda3_sq_sq
      real(wp) :: base, base_sq, e1, e2
      real(wp) :: sum_e1, sum_e2, sum_e3, sum_e3_sq
      real(wp) :: dsum_e1(ndim), dsum_e2(ndim), dsum_e3(ndim), dsum_e3_sq(ndim)
      real(wp) :: grad_d(ndim), hess_d(ndim, ndim), inv_x
      !> Cached outer product of grad_d with itself
      real(wp) :: nn(ndim, ndim)
      !> Per-kind blending coefficients combining s_1, s_2, s_3 with scalar
      !> sums. a_k is the coefficient of d2sum_e_k (and dsum_e_k) in the final
      !> assembled hessZ (and gradZ), derived from the elementary symmetric
      !> polynomial expansion.
      real(wp) :: a1, a2, a3, a4
      !> Fused outer-product weight coefficients (q_k = a_k * lambda_k^2)
      real(wp) :: q1, q2, q3, q4
      !> Fused hessian weight coefficients (p_k = a_k * lambda_k)
      real(wp) :: p1, p2, p3, p4
      !> Per-atom fused weights for outer-product and hessian contributions
      real(wp) :: w_nn, w_hd
      !> Whether to run the Hessian pass (pass 2) and Hessian cross-terms
      logical :: do_hess

      do_hess = present(hessZ)

      Z = 0.0_wp; gradZ = 0.0_wp
      if (do_hess) hessZ = 0.0_wp
      n_active = ssd_system%n_active
      if (n_active == 0) return

      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      lambda3_sq = 2.0_wp*k/3.0_wp
      lambda1_sq = lambda1*lambda1
      lambda2_sq = lambda2*lambda2
      lambda3_sq_sq = lambda3_sq*lambda3_sq

      ! === Pass 1: Scalar sums and gradient sums ===
      ! 4 scalar + 4x3 gradient accumulators = 16 values (fits in registers)
      sum_e1 = 0.0_wp; sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      dsum_e1 = 0.0_wp; dsum_e2 = 0.0_wp; dsum_e3 = 0.0_wp; dsum_e3_sq = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)

         sum_e1 = sum_e1 + e1
         sum_e2 = sum_e2 + e2
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq

         dsum_e1 = dsum_e1 - (lambda1*e1)*grad_d
         dsum_e2 = dsum_e2 - (lambda2*e2)*grad_d
         dsum_e3 = dsum_e3 - (lambda3*base)*grad_d
         dsum_e3_sq = dsum_e3_sq - (lambda3_sq*base_sq)*grad_d
      end do

      ! === Blending coefficients ===
      ! a_k = coefficient of d2sum_e_k in final hessZ assembly, derived by
      ! collecting like terms across the s_1, s_2, s_3 polynomial blocks.
      a1 = s_1
      a2 = 0.0_wp; a3 = 0.0_wp; a4 = 0.0_wp
      if (n_active >= 2) then
         a1 = a1 - 0.5_wp*s_2
         a2 = s_2*sum_e2
      end if
      if (n_active >= 3) then
         a1 = a1 + (1.0_wp/3.0_wp)*s_3
         a3 = 0.5_wp*s_3*(sum_e3*sum_e3 - sum_e3_sq)
         a4 = -0.5_wp*s_3*sum_e3
      end if

      ! Fused-weight constants: q_k = a_k * lambda_k^2, p_k = a_k * lambda_k
      ! These pre-blend the 4 exponential kinds so the hessian loop accumulates
      ! into a single (ndim, ndim) array instead of four.
      if (do_hess) then
         q1 = a1*lambda1_sq; p1 = a1*lambda1
         q2 = a2*lambda2_sq; p2 = a2*lambda2
         q3 = a3*lambda3*lambda3; p3 = a3*lambda3
         q4 = a4*lambda3_sq_sq; p4 = a4*lambda3_sq

         ! === Pass 2: Fused hessian accumulation ===
         ! 9 hessian accumulators (vs 36 in the unfused single-pass approach)
         do i = 1, n_active
            base = ssd_system%k3f0_arr(i)
            base_sq = base*base
            e1 = base_sq*base
            e2 = base*ssd_system%sqrt_k3f0_arr(i)

            ! Per-atom fused weights blend all 4 exponential kinds into one scalar
            w_nn = q1*e1 + q2*e2 + q3*base + q4*base_sq
            w_hd = p1*e1 + p2*e2 + p3*base + p4*base_sq

            grad_d = ssd_system%f1_r_arr(:, i)

            ! Reconstruct hess_d = (I - n n^T)/x from f1_r (= n) and inv_x instead of reading a stored f2_rr_arr
            inv_x = ssd_system%inv_x_arr(i)
            do j = 1, ndim
               do ii = 1, ndim
                  hess_d(ii, j) = -grad_d(ii)*grad_d(j)*inv_x
               end do
               hess_d(j, j) = hess_d(j, j) + inv_x
            end do

            nn(:, 1) = grad_d*grad_d(1)
            nn(:, 2) = grad_d*grad_d(2)
            nn(:, 3) = grad_d*grad_d(3)

            hessZ = hessZ + w_nn*nn - w_hd*hess_d
         end do
      end if

      ! === Assemble scalar Z from power sums ===
      Z = s_1*sum_e1
      if (n_active >= 2) then
         Z = Z + s_2*(sum_e2*sum_e2 - sum_e1)*0.5_wp
      end if
      if (n_active >= 3) then
         Z = Z + s_3*(sum_e3**3 - 3.0_wp*sum_e3*sum_e3_sq &
                      + 2.0_wp*sum_e1)/6.0_wp
      end if

      ! === Assemble gradient via blending coefficients ===
      gradZ = a1*dsum_e1 + a2*dsum_e2 + a3*dsum_e3 + a4*dsum_e3_sq

      ! === Add hessian cross-terms (rank-1 updates from gradient sums) ===
      if (do_hess) then
         if (n_active >= 2) then
            hessZ = hessZ + s_2*outer_matrix(dsum_e2, dsum_e2)
         end if
         if (n_active >= 3) then
            hessZ = hessZ + s_3*( &
                    sum_e3*outer_matrix(dsum_e3, dsum_e3) &
                    - 0.5_wp*(outer_matrix(dsum_e3, dsum_e3_sq) &
                              + outer_matrix(dsum_e3_sq, dsum_e3)))
         end if
      end if
   end subroutine compute_lsf_z012_rr_screened

!> Compute Z and all spatial derivatives up to third order using O(N) power sums
!> Three-pass fused-weight approach
!>  - scalar+gradient sums
!>  - fused hessian plus individual d2sum for third-order cross-terms
!>  - fused third-order accumulation (eliminates heap allocations and BLAS calls)
   subroutine compute_lsf_z0123_rrr_screened(ssd_system, k, s_1, s_2, s_3, &
                                             Z, gradZ, hessZ, thirdZ)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), intent(out) :: Z, gradZ(ndim), hessZ(ndim, ndim), thirdZ(ndim, ndim, ndim)

      integer :: i, j, n_active
      real(wp) :: lambda1, lambda2, lambda3, lambda3_sq
      real(wp) :: lambda1_sq, lambda2_sq, lambda3_sq_sq
      real(wp) :: lambda1_cub, lambda2_cub, lambda3_cub, lambda3_sq_cub
      real(wp) :: base, base_sq, e1, e2
      real(wp) :: sum_e1, sum_e2, sum_e3, sum_e3_sq
      real(wp) :: dsum_e1(ndim), dsum_e2(ndim), dsum_e3(ndim), dsum_e3_sq(ndim)
      !> Individual hessian sums for e2, e3, e3_sq kinds (needed by thirdZ
      !> cross-terms). The e1 contribution is absorbed into the fused hessZ.
      real(wp) :: d2sum_e2(ndim, ndim), d2sum_e3(ndim, ndim), d2sum_e3_sq(ndim, ndim)
      real(wp) :: grad_d(ndim), hess_d(ndim, ndim)
      !> Cached outer product of grad_d with itself
      real(wp) :: nn(ndim, ndim)
      !> Per-atom weighted coefficients for individual hessian contributions
      real(wp) :: c_nn, c_hd
      !> Local rank-3 tensors for sym3_21 and outer3 in third-order pass
      real(wp) :: shg_local(ndim, ndim, ndim), o3_local(ndim, ndim, ndim)
      !> Per-kind blending coefficients (a_k = coefficient of d2sum_e_k / d3sum_e_k
      !> in the final assembly, derived from the elementary symmetric polynomials)
      real(wp) :: a1, a2, a3, a4
      !> Fused outer-product weight coefficients (q_k = a_k * lambda_k^2)
      real(wp) :: q1, q2, q3, q4
      !> Fused hessian weight coefficients (p_k = a_k * lambda_k)
      real(wp) :: p1, p2, p3, p4
      !> Fused third-order outer3 weight coefficients (r_k = a_k * lambda_k^3)
      real(wp) :: r1, r2, r3, r4
      !> Per-atom fused weights for hessian and third-order contributions
      real(wp) :: w_nn, w_hd, w_t, w_s, w_o

      Z = 0.0_wp; gradZ = 0.0_wp; hessZ = 0.0_wp; thirdZ = 0.0_wp
      n_active = ssd_system%n_active
      if (n_active == 0) return

      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      lambda3_sq = 2.0_wp*k/3.0_wp
      lambda1_sq = lambda1*lambda1
      lambda2_sq = lambda2*lambda2
      lambda3_sq_sq = lambda3_sq*lambda3_sq
      lambda1_cub = lambda1_sq*lambda1
      lambda2_cub = lambda2_sq*lambda2
      lambda3_cub = lambda3*lambda3*lambda3
      lambda3_sq_cub = lambda3_sq_sq*lambda3_sq

      ! === Pass 1: Scalar sums and gradient sums ===
      ! 4 scalar + 4x3 gradient accumulators = 16 values (fits in registers)
      sum_e1 = 0.0_wp; sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      dsum_e1 = 0.0_wp; dsum_e2 = 0.0_wp; dsum_e3 = 0.0_wp; dsum_e3_sq = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)

         sum_e1 = sum_e1 + e1
         sum_e2 = sum_e2 + e2
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq

         dsum_e1 = dsum_e1 - (lambda1*e1)*grad_d
         dsum_e2 = dsum_e2 - (lambda2*e2)*grad_d
         dsum_e3 = dsum_e3 - (lambda3*base)*grad_d
         dsum_e3_sq = dsum_e3_sq - (lambda3_sq*base_sq)*grad_d
      end do

      ! === Blending coefficients ===
      a1 = s_1
      a2 = 0.0_wp; a3 = 0.0_wp; a4 = 0.0_wp
      if (n_active >= 2) then
         a1 = a1 - 0.5_wp*s_2
         a2 = s_2*sum_e2
      end if
      if (n_active >= 3) then
         a1 = a1 + (1.0_wp/3.0_wp)*s_3
         a3 = 0.5_wp*s_3*(sum_e3*sum_e3 - sum_e3_sq)
         a4 = -0.5_wp*s_3*sum_e3
      end if

      ! Fused-weight constants
      q1 = a1*lambda1_sq; p1 = a1*lambda1
      q2 = a2*lambda2_sq; p2 = a2*lambda2
      q3 = a3*lambda3*lambda3; p3 = a3*lambda3
      q4 = a4*lambda3_sq_sq; p4 = a4*lambda3_sq
      r1 = a1*lambda1_cub; r2 = a2*lambda2_cub
      r3 = a3*lambda3_cub; r4 = a4*lambda3_sq_cub

      ! === Pass 2: Fused hessian + individual d2sum for cross-terms ===
      ! hessZ(9) fused + d2sum_e2/e3/e3_sq(27) individual = 36 accumulators
      d2sum_e2 = 0.0_wp; d2sum_e3 = 0.0_wp; d2sum_e3_sq = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)
         hess_d = ssd_system%f2_rr_arr(:, :, i)

         nn(:, 1) = grad_d*grad_d(1)
         nn(:, 2) = grad_d*grad_d(2)
         nn(:, 3) = grad_d*grad_d(3)

         ! Fused hessian (blends all 4 kinds into one accumulator)
         w_nn = q1*e1 + q2*e2 + q3*base + q4*base_sq
         w_hd = p1*e1 + p2*e2 + p3*base + p4*base_sq
         hessZ = hessZ + w_nn*nn - w_hd*hess_d

         ! Individual d2sum for e2, e3, e3_sq (needed by thirdZ cross-terms)
         c_nn = lambda2_sq*e2; c_hd = lambda2*e2
         d2sum_e2 = d2sum_e2 + c_nn*nn - c_hd*hess_d
         c_nn = lambda3*lambda3*base; c_hd = lambda3*base
         d2sum_e3 = d2sum_e3 + c_nn*nn - c_hd*hess_d
         c_nn = lambda3_sq_sq*base_sq; c_hd = lambda3_sq*base_sq
         d2sum_e3_sq = d2sum_e3_sq + c_nn*nn - c_hd*hess_d
      end do

      ! === Pass 3: Fused third-order accumulation ===
      ! 27 thirdZ accumulators - replaces 6 heap allocations + 3 DGEMMs.
      ! Per-atom weights blend all 4 exponential kinds into 3 scalars:
      !   w_t (f3_rrr), w_s (sym3_21), w_o (outer3)
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)

         w_t = -(p1*e1 + p2*e2 + p3*base + p4*base_sq)
         w_s = q1*e1 + q2*e2 + q3*base + q4*base_sq
         w_o = -(r1*e1 + r2*e2 + r3*base + r4*base_sq)

         grad_d = ssd_system%f1_r_arr(:, i)
         hess_d = ssd_system%f2_rr_arr(:, :, i)

         nn(:, 1) = grad_d*grad_d(1)
         nn(:, 2) = grad_d*grad_d(2)
         nn(:, 3) = grad_d*grad_d(3)

         shg_local = sym3_21(hess_d, grad_d)
         do j = 1, ndim
            o3_local(:, :, j) = nn*grad_d(j)
         end do

         thirdZ = thirdZ + w_t*ssd_system%f3_rrr_arr(:, :, :, i) &
                  + w_s*shg_local + w_o*o3_local
      end do

      ! === Assemble scalar Z from power sums ===
      Z = s_1*sum_e1
      if (n_active >= 2) then
         Z = Z + s_2*(sum_e2*sum_e2 - sum_e1)*0.5_wp
      end if
      if (n_active >= 3) then
         Z = Z + s_3*(sum_e3**3 - 3.0_wp*sum_e3*sum_e3_sq &
                      + 2.0_wp*sum_e1)/6.0_wp
      end if

      ! === Assemble gradient via blending coefficients ===
      gradZ = a1*dsum_e1 + a2*dsum_e2 + a3*dsum_e3 + a4*dsum_e3_sq

      ! === Add hessian cross-terms (rank-1 updates from gradient sums) ===
      if (n_active >= 2) then
         hessZ = hessZ + s_2*outer_matrix(dsum_e2, dsum_e2)
      end if
      if (n_active >= 3) then
         hessZ = hessZ + s_3*( &
                 sum_e3*outer_matrix(dsum_e3, dsum_e3) &
                 - 0.5_wp*(outer_matrix(dsum_e3, dsum_e3_sq) &
                           + outer_matrix(dsum_e3_sq, dsum_e3)))
      end if

      ! === Add third-order cross-terms (d2sum x dsum interactions) ===
      if (n_active >= 2) then
         thirdZ = thirdZ + s_2*sym3_21(d2sum_e2, dsum_e2)
      end if
      if (n_active >= 3) then
         thirdZ = thirdZ + s_3*( &
                  outer3(dsum_e3) &
                  + sum_e3*sym3_21(d2sum_e3, dsum_e3) &
                  - 0.5_wp*(sym3_21(d2sum_e3, dsum_e3_sq) &
                            + sym3_21(d2sum_e3_sq, dsum_e3)))
      end if
   end subroutine compute_lsf_z0123_rrr_screened

!> Compute nuclear derivatives of Z Hessian (screened version)
   subroutine compute_lsf_z3_rr_rA_screened(ssd_system, k, s_1, s_2, s_3, deriv)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), allocatable, intent(out) :: deriv(:, :, :, :)

      integer :: axis, i, j, n_active
      real(wp) :: lambda1, lambda2, lambda3, lambda3_sq
      real(wp) :: lambda1_sq, lambda2_sq, lambda3_sq_sq
      real(wp) :: lambda1_cub, lambda2_cub, lambda3_cub, lambda3_sq_cub
      real(wp) :: base, base_sq, e1, e2
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: dsum_e2(ndim), dsum_e3(ndim), dsum_e3_sq(ndim)
      real(wp) :: d2sum_e2(ndim, ndim), d2sum_e3(ndim, ndim), d2sum_e3_sq(ndim, ndim)
      real(wp) :: grad_d(ndim), hess_d(ndim, ndim), third_d(ndim, ndim, ndim)
      !> Cached outer product of grad_d with itself
      real(wp) :: nn(ndim, ndim)
      !> Per-kind blending coefficients combining s_1, s_2, s_3 with scalar
      !> sums. a_k is the coefficient of d2sum_e_k in the final assembled hessZ,
      !> derived from the elementary symmetric polynomial expansion.
      real(wp) :: a1, a2, a3, a4
      !> Fused-weight constants: p_k = a_k*lambda_k, q_k = a_k*lambda_k^2,
      !> r_k = a_k*lambda_k^3. Pre-blend the 4 exponential kinds so pass 2
      !> computes 3 scalar weights instead of 4 separate hc matrices.
      real(wp) :: p1, p2, p3, p4, q1, q2, q3, q4, r1, r2, r3, r4
      !> Per-atom fused scalar weights:
      !>  w_A = sum_k r_k * e_k  (for nn in hc_combined)
      !>  w_B = sum_k q_k * e_k  (for hess_d in hc_combined, and nn_mixed)
      !>  w_C = sum_k p_k * e_k  (for third_d contribution)
      real(wp) :: w_A, w_B, w_C
      !> Per-atom cross-term scalars: P_k = -lambda_k^2 * e_k (grad_d(axis)
      !> coefficient in dgrad_e_k), Q_k = -lambda_k * e_k (h coefficient)
      real(wp) :: P_e2, Q_e2, P_e3, Q_e3, P_e3sq, Q_e3sq
      !> Precomputed constant vectors (from gradient sums, constant across atoms)
      real(wp) :: s_2_dsum_e2(ndim), s_3_v3(ndim), s_3_v3sq(ndim)
      !> Precomputed constant matrix (from hessian sums, constant across atoms)
      real(wp) :: M3_const(ndim, ndim)
      !> Per-atom fused outer-product vector U = sum_k P_k * C_k => sym_outer(grad_d, U) added to T
      real(wp) :: U(ndim)
      !> Per-atom fused outer-product vector V = w_B * grad_d + sum_k Q_k * C_k => sym_outer(h, V) per axis
      real(wp) :: V(ndim)
      !> Combined per-atom matrix: hc_combined + M2 + M3 + sym_outer(grad_d, U)
      real(wp) :: T(ndim, ndim)
      !> Per-axis column of -hess_d
      real(wp) :: h(ndim)

      n_active = ssd_system%n_active
      allocate (deriv(ndim, ndim, ndim, n_active))
      deriv = 0.0_wp
      if (n_active == 0) return

      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      lambda3_sq = 2.0_wp*k/3.0_wp
      lambda1_sq = lambda1*lambda1
      lambda2_sq = lambda2*lambda2
      lambda3_sq_sq = lambda3_sq*lambda3_sq
      lambda1_cub = lambda1_sq*lambda1
      lambda2_cub = lambda2_sq*lambda2
      lambda3_cub = lambda3*lambda3*lambda3
      lambda3_sq_cub = lambda3_sq_sq*lambda3_sq

      ! === Pass 1a: Scalar sums and gradient sums ===
      ! 3 scalar + 3x3 gradient accumulators = 12 values (fits in registers).
      ! sum_e1 / dsum_e1 are not needed: this routine computes the nuclear
      ! derivative of hessZ, and the e1 contribution is fully fused in pass 2.
      sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      dsum_e2 = 0.0_wp; dsum_e3 = 0.0_wp; dsum_e3_sq = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)

         sum_e2 = sum_e2 + e2
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq

         dsum_e2 = dsum_e2 - (lambda2*e2)*grad_d
         dsum_e3 = dsum_e3 - (lambda3*base)*grad_d
         dsum_e3_sq = dsum_e3_sq - (lambda3_sq*base_sq)*grad_d
      end do

      ! === Blending coefficients ===
      a1 = s_1
      a2 = 0.0_wp; a3 = 0.0_wp; a4 = 0.0_wp
      s_2_dsum_e2 = 0.0_wp
      s_3_v3 = 0.0_wp; s_3_v3sq = 0.0_wp
      if (n_active >= 2) then
         a1 = a1 - 0.5_wp*s_2
         a2 = s_2*sum_e2
         s_2_dsum_e2 = s_2*dsum_e2
      end if
      if (n_active >= 3) then
         a1 = a1 + (1.0_wp/3.0_wp)*s_3
         a3 = 0.5_wp*s_3*(sum_e3*sum_e3 - sum_e3_sq)
         a4 = -0.5_wp*s_3*sum_e3
         s_3_v3 = s_3*(sum_e3*dsum_e3 - 0.5_wp*dsum_e3_sq)
         s_3_v3sq = s_3*(-0.5_wp*dsum_e3)
      end if

      ! Fused-weight constants for pass 2 (same convention as z012/z0123)
      p1 = a1*lambda1; q1 = a1*lambda1_sq
      p2 = a2*lambda2; q2 = a2*lambda2_sq
      p3 = a3*lambda3; q3 = a3*lambda3*lambda3
      p4 = a4*lambda3_sq; q4 = a4*lambda3_sq_sq
      r1 = a1*lambda1_cub; r2 = a2*lambda2_cub
      r3 = a3*lambda3_cub; r4 = a4*lambda3_sq_cub

      ! === Pass 1b: Hessian sums (d2sum_e2, d2sum_e3, d2sum_e3_sq) ===
      ! 27 accumulators. d2sum_e1 is not needed: its contribution is captured
      ! entirely by the fused scalars w_A, w_B, w_C in pass 2.
      d2sum_e2 = 0.0_wp; d2sum_e3 = 0.0_wp; d2sum_e3_sq = 0.0_wp
      if (n_active >= 2) then
         do i = 1, n_active
            base = ssd_system%k3f0_arr(i)
            base_sq = base*base
            e2 = base*ssd_system%sqrt_k3f0_arr(i)
            grad_d = ssd_system%f1_r_arr(:, i)
            hess_d = ssd_system%f2_rr_arr(:, :, i)

            nn(:, 1) = grad_d*grad_d(1)
            nn(:, 2) = grad_d*grad_d(2)
            nn(:, 3) = grad_d*grad_d(3)

            d2sum_e2 = d2sum_e2 + e2*(lambda2_sq*nn - lambda2*hess_d)
            if (n_active >= 3) then
               d2sum_e3 = d2sum_e3 &
                          + base*(lambda3*lambda3*nn - lambda3*hess_d)
               d2sum_e3_sq = d2sum_e3_sq &
                             + base_sq*(lambda3_sq_sq*nn - lambda3_sq*hess_d)
            end if
         end do
      end if

      ! Precompute M3_const from hessian sums (constant across atoms)
      M3_const = 0.0_wp
      if (n_active >= 3) then
         M3_const = sum_e3*d2sum_e3 - 0.5_wp*d2sum_e3_sq
         M3_const(:, 1) = M3_const(:, 1) + dsum_e3*dsum_e3(1)
         M3_const(:, 2) = M3_const(:, 2) + dsum_e3*dsum_e3(2)
         M3_const(:, 3) = M3_const(:, 3) + dsum_e3*dsum_e3(3)
      end if

      ! === Pass 2: Per-atom nuclear derivatives (fused weights) ===
      ! Per-axis work is reduced from 8 outer_matrix calls + nn_mixed to a single
      ! sym_outer(h, V), by decomposing dgrad_e_k = P_k*grad_d(axis)*grad_d + Q_k*h
      ! and collecting like terms across all exponential kinds.
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)
         hess_d = ssd_system%f2_rr_arr(:, :, i)
         third_d = ssd_system%f3_rrr_arr(:, :, :, i)

         ! Fused scalar weights: 3 scalars replace 4 separate hc matrices
         w_A = r1*e1 + r2*e2 + r3*base + r4*base_sq
         w_B = q1*e1 + q2*e2 + q3*base + q4*base_sq
         w_C = p1*e1 + p2*e2 + p3*base + p4*base_sq

         nn(:, 1) = grad_d*grad_d(1)
         nn(:, 2) = grad_d*grad_d(2)
         nn(:, 3) = grad_d*grad_d(3)

         ! T = fused hc_combined: replaces w_hc1*hc1 + ... + w_hc3sq*hc3sq
         T = w_A*nn - w_B*hess_d

         ! Add 2-body and 3-body scalar*matrix terms to T
         if (n_active >= 2) then
            T = T + s_2*lambda2*e2*d2sum_e2
         end if
         if (n_active >= 3) then
            T = T + s_3*(lambda3*base*M3_const &
                         - 0.5_wp*lambda3_sq*base_sq*d2sum_e3)
         end if

         ! Fuse all cross-term outer products into T (via U) and V:
         ! Each dgrad_e_k = P_k * grad_d(axis) * grad_d + Q_k * h, so
         ! sym_outer(dgrad_e_k, C_k) splits into an axis-independent part
         ! (P_k * grad_d(axis) * sym_outer(grad_d, C_k), absorbed into T)
         ! and an axis-dependent part (Q_k * sym_outer(h, C_k), fused into V).
         U = 0.0_wp
         V = w_B*grad_d
         if (n_active >= 2) then
            P_e2 = -lambda2_sq*e2
            Q_e2 = -lambda2*e2
            U = P_e2*s_2_dsum_e2
            V = V + Q_e2*s_2_dsum_e2
         end if
         if (n_active >= 3) then
            P_e3 = -lambda3*lambda3*base
            Q_e3 = -lambda3*base
            P_e3sq = -lambda3_sq_sq*base_sq
            Q_e3sq = -lambda3_sq*base_sq
            U = U + P_e3*s_3_v3 + P_e3sq*s_3_v3sq
            V = V + Q_e3*s_3_v3 + Q_e3sq*s_3_v3sq
         end if

         ! Add SG = sym_outer(grad_d, U) to T
         do j = 1, ndim
            T(:, j) = T(:, j) + grad_d*U(j) + U*grad_d(j)
         end do

         ! Per-axis: deriv = grad_d(axis)*T + w_C*third_d + sym_outer(h, V)
         do axis = 1, ndim
            h = -hess_d(:, axis)
            deriv(:, :, axis, i) = grad_d(axis)*T + w_C*third_d(:, :, axis)
            do j = 1, ndim
               deriv(:, j, axis, i) = deriv(:, j, axis, i) + h*V(j) + V*h(j)
            end do
         end do
      end do
   end subroutine compute_lsf_z3_rr_rA_screened

!> Compute Z and all spatial r-derivatives up to fourth order using per-kind
!> Faa di Bruno on the four "kinds" of exponentials.
!> TODO: This is not yet very aggressively fused; so here is still optimization potential
   subroutine compute_lsf_z01234_rrrr_screened(ssd_system, k, s_1, s_2, s_3, Z, gradZ, hessZ, thirdZ, fourthZ)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), intent(out) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp), intent(out) :: thirdZ(ndim, ndim, ndim)
      real(wp), intent(out) :: fourthZ(ndim, ndim, ndim, ndim)

      integer :: i, n_active
      real(wp) :: lambda(4)
      real(wp) :: base, base_sq, sqrt_base, e_val(4)
      real(wp) :: S(4)                                    ! S1, S2, S3, T3
      real(wp) :: dS(ndim, 4)
      real(wp) :: d2S(ndim, ndim, 4)
      real(wp) :: d3S(ndim, ndim, ndim, 4)
      real(wp) :: d4S(ndim, ndim, ndim, ndim, 4)
      real(wp) :: f1(ndim), f2(ndim, ndim)
      real(wp) :: f3(ndim, ndim, ndim), f4(ndim, ndim, ndim, ndim)
      real(wp) :: v0, v1(ndim), v2(ndim, ndim)
      real(wp) :: v3(ndim, ndim, ndim), v4(ndim, ndim, ndim, ndim)
      real(wp) :: C1, half_s_2, s_36, neg_half_s_3
      integer :: kk

      Z = 0.0_wp; gradZ = 0.0_wp; hessZ = 0.0_wp
      thirdZ = 0.0_wp; fourthZ = 0.0_wp
      n_active = ssd_system%n_active
      if (n_active == 0) return

      ! lambda(1)=k (for e1), lambda(2)=k/2 (e2), lambda(3)=k/3 (e3),
      ! lambda(4)=2k/3 (e3^2; treat T3 = sum e3^2 = sum exp(-2k/3 d) as the 4th "kind")
      lambda(1) = k
      lambda(2) = 0.5_wp*k
      lambda(3) = k/3.0_wp
      lambda(4) = 2.0_wp*k/3.0_wp

      S = 0.0_wp; dS = 0.0_wp; d2S = 0.0_wp; d3S = 0.0_wp; d4S = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = ssd_system%sqrt_k3f0_arr(i)
         e_val(1) = base_sq*base       ! e1
         e_val(2) = base*sqrt_base     ! e2
         e_val(3) = base                  ! e3
         e_val(4) = base_sq               ! e3^2

         f1 = ssd_system%f1_r_arr(:, i)
         f2 = ssd_system%f2_rr_arr(:, :, i)
         f3 = ssd_system%f3_rrr_arr(:, :, :, i)
         f4 = ssd_system%f4_rrrr_arr(:, :, :, :, i)

         do kk = 1, 4
            call atom_exp_derivs(e_val(kk), lambda(kk), f1, f2, f3, f4, v0, v1, v2, v3, v4)
            S(kk) = S(kk) + v0
            dS(:, kk) = dS(:, kk) + v1
            d2S(:, :, kk) = d2S(:, :, kk) + v2
            d3S(:, :, :, kk) = d3S(:, :, :, kk) + v3
            d4S(:, :, :, :, kk) = d4S(:, :, :, :, kk) + v4
         end do
      end do

      ! Coefficients of the polynomial
      !   Z = C1 * S1 + (s_2/2) * S2^2 + (s_3/6) * S3^3 - (s_3/2) * S3 * T3
      C1 = s_1
      if (n_active >= 2) C1 = C1 - 0.5_wp*s_2
      if (n_active >= 3) C1 = C1 + s_3/3.0_wp
      half_s_2 = 0.5_wp*s_2
      s_36 = s_3/6.0_wp
      neg_half_s_3 = -0.5_wp*s_3

      ! Linear term: C1 * S1
      Z = Z + C1*S(1)
      gradZ = gradZ + C1*dS(:, 1)
      hessZ = hessZ + C1*d2S(:, :, 1)
      thirdZ = thirdZ + C1*d3S(:, :, :, 1)
      fourthZ = fourthZ + C1*d4S(:, :, :, :, 1)

      ! Quadratic 2-body term: (s_2/2) * S2^2  (only when n_active >= 2)
      if (n_active >= 2) then
         call accum_uv_0(Z, half_s_2, S(2), S(2))
         call accum_uv_1(gradZ, half_s_2, S(2), dS(:, 2), S(2), dS(:, 2))
         call accum_uv_2(hessZ, half_s_2, S(2), dS(:, 2), d2S(:, :, 2), &
                         S(2), dS(:, 2), d2S(:, :, 2))
         call accum_uv_3(thirdZ, half_s_2, S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2), &
                         S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2))
         call accum_uv_4(fourthZ, half_s_2, S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2), d4S(:, :, :, :, 2), &
                         S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2), d4S(:, :, :, :, 2))
      end if

      ! Cubic 3-body S3^3 term: (s_3/6) * S3^3.
      ! Use double product (S3 * S3) * S3 with two-step Leibniz.
      ! Step 1: build P = S3 * S3 and its derivatives (call them PS).
      ! Step 2: assemble S3 * P contributions.
      if (n_active >= 3) then
         call accum_S3_cubed(fourthZ, thirdZ, hessZ, gradZ, Z, s_36, &
                             S(3), dS(:, 3), d2S(:, :, 3), d3S(:, :, :, 3), d4S(:, :, :, :, 3))

         ! Mixed term: -(s_3/2) * S3 * T3
         call accum_uv_0(Z, neg_half_s_3, S(3), S(4))
         call accum_uv_1(gradZ, neg_half_s_3, S(3), dS(:, 3), S(4), dS(:, 4))
         call accum_uv_2(hessZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), &
                         S(4), dS(:, 4), d2S(:, :, 4))
         call accum_uv_3(thirdZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), d3S(:, :, :, 3), &
                         S(4), dS(:, 4), d2S(:, :, 4), d3S(:, :, :, 4))
         call accum_uv_4(fourthZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), d3S(:, :, :, 3), d4S(:, :, :, :, 3), &
                         S(4), dS(:, 4), d2S(:, :, 4), d3S(:, :, :, 4), d4S(:, :, :, :, 4))
      end if
   end subroutine compute_lsf_z01234_rrrr_screened

!> Compute Z r-derivative tensors and their first R-derivatives
!>
!> Outputs (allocated by this routine):
!>   dZ_rA(s_1, iA) = dZ/ dR_{A,s_1} (n_active=iA)
!>   dgradZ_rA(j, s_1, iA) = d(Z_j)/ dR_{A,s_1}
!>   dhessZ_rA(j, k, s_1, iA) = d(Z_{jk})/ dR_{A,s_1}
!>   dthirdZ_rA(j, k, l, s_1, iA) = d(Z_{jkl})/ dR_{A,s_1}
!>
!> Strategy: for each active atom A we (re-)assemble Z r-derivatives but with
!> one factor of S_k replaced by its "atom-A R-derivative" version (which is
!> the negative of the per-atom-A contribution to the next r-derivative of
!> S_k).  By Leibniz this gives the full R_A-derivative of the Z tensors.
   subroutine compute_lsf_zR_derivs_rrr_screened(ssd_system, k, s_1, s_2, s_3, &
                                                 dZ_rA, dgradZ_rA, dhessZ_rA, dthirdZ_rA)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), allocatable, intent(out) :: dZ_rA(:, :)
      real(wp), allocatable, intent(out) :: dgradZ_rA(:, :, :)
      real(wp), allocatable, intent(out) :: dhessZ_rA(:, :, :, :)
      real(wp), allocatable, intent(out) :: dthirdZ_rA(:, :, :, :, :)

      integer :: i, n_active, iA, kk, axis
      real(wp) :: lambda(4)
      real(wp) :: base, base_sq, sqrt_base, e_val(4)
      real(wp) :: S(4)
      real(wp) :: dS(ndim, 4)
      real(wp) :: d2S(ndim, ndim, 4)
      real(wp) :: d3S(ndim, ndim, ndim, 4)
      ! Per-atom Faa di Bruno tensors of e^{-lambda_k d_A}:
      !   sA(kind, iA), dsA(:, kind, iA), d2sA(:,:, kind, iA), d3sA(:,:,:, kind, iA),
      !   d4sA(:,:,:,:, kind, iA).  Needed: up to order 4, since the R-derivative
      !   of (d^3 S_k) requires the per-atom-A piece of (d^4 S_k).
      real(wp), allocatable :: sA(:, :)
      real(wp), allocatable :: dsA(:, :, :)
      real(wp), allocatable :: d2sA(:, :, :, :)
      real(wp), allocatable :: d3sA(:, :, :, :, :)
      real(wp), allocatable :: d4sA(:, :, :, :, :, :)
      real(wp) :: f1(ndim), f2(ndim, ndim)
      real(wp) :: f3(ndim, ndim, ndim), f4(ndim, ndim, ndim, ndim)
      real(wp) :: v0, v1(ndim), v2(ndim, ndim)
      real(wp) :: v3(ndim, ndim, ndim), v4(ndim, ndim, ndim, ndim)

      ! "R-derivative" of S_k at atom A, axis s_1 (denoted RS_k below) is a
      ! scalar; "R-derivative of dS_k_a" is a vector; etc.  We assemble these
      ! per (iA, axis) and feed them via the Leibniz accumulators.
      real(wp) :: RS(4)
      real(wp) :: RdS(ndim, 4)
      real(wp) :: Rd2S(ndim, ndim, 4)
      real(wp) :: Rd3S(ndim, ndim, ndim, 4)
      real(wp) :: C1, half_s_2, s_36, neg_half_s_3
      real(wp) :: dZ_val, dgradZ(ndim), dhessZ(ndim, ndim), dthirdZ(ndim, ndim, ndim)

      n_active = ssd_system%n_active
      allocate (dZ_rA(ndim, n_active))
      allocate (dgradZ_rA(ndim, ndim, n_active))
      allocate (dhessZ_rA(ndim, ndim, ndim, n_active))
      allocate (dthirdZ_rA(ndim, ndim, ndim, ndim, n_active))
      dZ_rA = 0.0_wp; dgradZ_rA = 0.0_wp; dhessZ_rA = 0.0_wp; dthirdZ_rA = 0.0_wp
      if (n_active == 0) return

      lambda(1) = k
      lambda(2) = 0.5_wp*k
      lambda(3) = k/3.0_wp
      lambda(4) = 2.0_wp*k/3.0_wp

      ! === Pass 1: build per-atom Faa di Bruno tensors AND global sums ===
      allocate (sA(4, n_active))
      allocate (dsA(ndim, 4, n_active))
      allocate (d2sA(ndim, ndim, 4, n_active))
      allocate (d3sA(ndim, ndim, ndim, 4, n_active))
      allocate (d4sA(ndim, ndim, ndim, ndim, 4, n_active))
      S = 0.0_wp; dS = 0.0_wp; d2S = 0.0_wp; d3S = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = ssd_system%sqrt_k3f0_arr(i)
         e_val(1) = base_sq*base
         e_val(2) = base*sqrt_base
         e_val(3) = base
         e_val(4) = base_sq

         f1 = ssd_system%f1_r_arr(:, i)
         f2 = ssd_system%f2_rr_arr(:, :, i)
         f3 = ssd_system%f3_rrr_arr(:, :, :, i)
         f4 = ssd_system%f4_rrrr_arr(:, :, :, :, i)

         do kk = 1, 4
            call atom_exp_derivs(e_val(kk), lambda(kk), f1, f2, f3, f4, v0, v1, v2, v3, v4)
            sA(kk, i) = v0
            dsA(:, kk, i) = v1
            d2sA(:, :, kk, i) = v2
            d3sA(:, :, :, kk, i) = v3
            d4sA(:, :, :, :, kk, i) = v4
            S(kk) = S(kk) + v0
            dS(:, kk) = dS(:, kk) + v1
            d2S(:, :, kk) = d2S(:, :, kk) + v2
            d3S(:, :, :, kk) = d3S(:, :, :, kk) + v3
         end do
      end do

      C1 = s_1
      if (n_active >= 2) C1 = C1 - 0.5_wp*s_2
      if (n_active >= 3) C1 = C1 + s_3/3.0_wp
      half_s_2 = 0.5_wp*s_2
      s_36 = s_3/6.0_wp
      neg_half_s_3 = -0.5_wp*s_3

      ! === Pass 2: for each (iA, axis), build the R-derivative tensors ===
      ! R-derivative S^R_k(axis) := dS_k/ dR_{A,axis} = -(per-atom-A piece of dS_k_axis)
      !                            = -dsA(axis, kk, iA)
      ! R-derivative of d^n S_k at remaining n indices: -d^{n+1}sA[axis, n indices, kk, iA]
      ! We compute -d^{(n+1)}sA[axis, ..., kk, iA] for each axis.
      do iA = 1, n_active
         do axis = 1, ndim
            ! RS_k = -dsA(axis, kk, iA)
            do kk = 1, 4
               RS(kk) = -dsA(axis, kk, iA)
               RdS(:, kk) = -d2sA(:, axis, kk, iA)
               Rd2S(:, :, kk) = -d3sA(:, :, axis, kk, iA)
               Rd3S(:, :, :, kk) = -d4sA(:, :, :, axis, kk, iA)
            end do

            dZ_val = 0.0_wp
            dgradZ = 0.0_wp
            dhessZ = 0.0_wp
            dthirdZ = 0.0_wp

            ! Linear term: C1 * S1.  R-derivative: C1 * RS1.
            dZ_val = dZ_val + C1*RS(1)
            dgradZ = dgradZ + C1*RdS(:, 1)
            dhessZ = dhessZ + C1*Rd2S(:, :, 1)
            dthirdZ = dthirdZ + C1*Rd3S(:, :, :, 1)

            ! Quadratic (s_2/2) S2^2:  R-derivative is 2 * (s_2/2) * S2 * RS2 etc.
            ! Use accum_uv_n with (u=S2_R-derivative, v=S2) plus (u=S2, v=S2_R-derivative).
            ! Equivalently: 2 * accum_uv_n with (u=RS2_*, v=S2_*).
            ! Here RS2's "value" at the chosen axis is a scalar RS(2); its "first
            ! r-derivative" tensor along (j) is RdS(j,2); etc.  We can use uv_n
            ! with u being the R-derivative tensor of S2 (which has its own
            ! r-derivatives RS, RdS, Rd2S, Rd3S) and v=S2 (with S, dS, d2S, d3S).
            if (n_active >= 2) then
               call accum_uv_0(dZ_val, half_s_2, RS(2), S(2))
               call accum_uv_0(dZ_val, half_s_2, S(2), RS(2))
               call accum_uv_1(dgradZ, half_s_2, RS(2), RdS(:, 2), S(2), dS(:, 2))
               call accum_uv_1(dgradZ, half_s_2, S(2), dS(:, 2), RS(2), RdS(:, 2))
               call accum_uv_2(dhessZ, half_s_2, RS(2), RdS(:, 2), Rd2S(:, :, 2), &
                               S(2), dS(:, 2), d2S(:, :, 2))
               call accum_uv_2(dhessZ, half_s_2, S(2), dS(:, 2), d2S(:, :, 2), &
                               RS(2), RdS(:, 2), Rd2S(:, :, 2))
               call accum_uv_3(dthirdZ, half_s_2, RS(2), RdS(:, 2), Rd2S(:, :, 2), Rd3S(:, :, :, 2), &
                               S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2))
               call accum_uv_3(dthirdZ, half_s_2, S(2), dS(:, 2), d2S(:, :, 2), d3S(:, :, :, 2), &
                               RS(2), RdS(:, 2), Rd2S(:, :, 2), Rd3S(:, :, :, 2))
            end if

            if (n_active >= 3) then
               ! Cubic S3^3:  R-derivative is 3 * S3^2 * RS3.  Implement via two-step
               ! product (S * S) * S, using Leibniz with R-derivative on each factor.
               ! Easier: compute R-derivative of S3^3 as 3 * S3^2 * RS3 directly,
               ! and its r-derivatives via Leibniz on the (S3^2)(RS3) and similar.
               call accum_R_S3_cubed(dthirdZ, dhessZ, dgradZ, dZ_val, s_36, &
                                     S(3), dS(:, 3), d2S(:, :, 3), d3S(:, :, :, 3), &
                                     RS(3), RdS(:, 3), Rd2S(:, :, 3), Rd3S(:, :, :, 3))

               ! Mixed -(s_3/2) S3 T3.  R-derivative = -(s_3/2) (RS3 T3 + S3 RT3).
               call accum_uv_0(dZ_val, neg_half_s_3, RS(3), S(4))
               call accum_uv_0(dZ_val, neg_half_s_3, S(3), RS(4))
               call accum_uv_1(dgradZ, neg_half_s_3, RS(3), RdS(:, 3), S(4), dS(:, 4))
               call accum_uv_1(dgradZ, neg_half_s_3, S(3), dS(:, 3), RS(4), RdS(:, 4))
               call accum_uv_2(dhessZ, neg_half_s_3, RS(3), RdS(:, 3), Rd2S(:, :, 3), &
                               S(4), dS(:, 4), d2S(:, :, 4))
               call accum_uv_2(dhessZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), &
                               RS(4), RdS(:, 4), Rd2S(:, :, 4))
               call accum_uv_3(dthirdZ, neg_half_s_3, RS(3), RdS(:, 3), Rd2S(:, :, 3), Rd3S(:, :, :, 3), &
                               S(4), dS(:, 4), d2S(:, :, 4), d3S(:, :, :, 4))
               call accum_uv_3(dthirdZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), d3S(:, :, :, 3), &
                               RS(4), RdS(:, 4), Rd2S(:, :, 4), Rd3S(:, :, :, 4))
            end if

            dZ_rA(axis, iA) = dZ_val
            dgradZ_rA(:, axis, iA) = dgradZ
            dhessZ_rA(:, :, axis, iA) = dhessZ
            dthirdZ_rA(:, :, :, axis, iA) = dthirdZ
         end do
      end do

      deallocate (sA, dsA, d2sA, d3sA, d4sA)
   end subroutine compute_lsf_zR_derivs_rrr_screened

!> Compute Z r-derivatives AND their second R-derivatives (a pair of R atoms).
!>
!> Outputs:
!>   d2Z_rArB(s_1A, iA, s_2B, iB) = d^2 Z/ dR_A dR_B
!>   d2gradZ_rArB(j, s_1A, iA, s_2B, iB) = d^2 (Z_j)/ dR_A dR_B
!>   d2hessZ_rArB(j, k, s_1A, iA, s_2B, iB) = d^2 (Z_{jk})/ dR_A dR_B
!>
!> Strategy: P_A acts on atoms equal to A; P_B acts on atoms equal to B.  When
!> iA == iB they coincide (single atom, second r-derivative); otherwise the two
!> restrictions are independent and the cross term is a product of single-atom
!> restrictions in each factor of the polynomial.
   subroutine compute_lsf_zRR_derivs_rr_screened(ssd_system, k, s_1, s_2, s_3, &
                                                 d2Z_rArB, d2gradZ_rArB, d2hessZ_rArB)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3
      real(wp), allocatable, intent(out) :: d2Z_rArB(:, :, :, :)
      real(wp), allocatable, intent(out) :: d2gradZ_rArB(:, :, :, :, :)
      real(wp), allocatable, intent(out) :: d2hessZ_rArB(:, :, :, :, :, :)

      integer :: i, n_active, iA, iB, axa, axb, kk
      real(wp) :: lambda(4)
      real(wp) :: base, base_sq, sqrt_base, e_val(4)
      real(wp) :: S(4), dS(ndim, 4), d2S(ndim, ndim, 4)
      real(wp), allocatable :: sA(:, :), dsA(:, :, :), d2sA(:, :, :, :)
      real(wp), allocatable :: d3sA(:, :, :, :, :), d4sA(:, :, :, :, :, :)
      real(wp) :: f1(ndim), f2(ndim, ndim)
      real(wp) :: f3(ndim, ndim, ndim), f4(ndim, ndim, ndim, ndim)
      real(wp) :: v0, v1(ndim), v2(ndim, ndim)
      real(wp) :: v3(ndim, ndim, ndim), v4(ndim, ndim, ndim, ndim)
      real(wp) :: C1, half_s_2, s_36, neg_half_s_3
      ! Per-axis R-derivative tensors (atom A, axis s_1A) for each kind kk:
      real(wp) :: RSa(4), RdSa(ndim, 4), Rd2Sa(ndim, ndim, 4)
      real(wp) :: RSb(4), RdSb(ndim, 4), Rd2Sb(ndim, ndim, 4)
      ! Second R-derivatives (s_1A, s_2B) for each kind kk:
      real(wp) :: RRS(4), RRdS(ndim, 4), RRd2S(ndim, ndim, 4)
      real(wp) :: rrZ_val, rrGZ(ndim), rrHZ(ndim, ndim)

      n_active = ssd_system%n_active
      allocate (d2Z_rArB(ndim, n_active, ndim, n_active))
      allocate (d2gradZ_rArB(ndim, ndim, n_active, ndim, n_active))
      allocate (d2hessZ_rArB(ndim, ndim, ndim, n_active, ndim, n_active))
      d2Z_rArB = 0.0_wp
      d2gradZ_rArB = 0.0_wp
      d2hessZ_rArB = 0.0_wp
      if (n_active == 0) return

      lambda(1) = k
      lambda(2) = 0.5_wp*k
      lambda(3) = k/3.0_wp
      lambda(4) = 2.0_wp*k/3.0_wp

      allocate (sA(4, n_active))
      allocate (dsA(ndim, 4, n_active))
      allocate (d2sA(ndim, ndim, 4, n_active))
      allocate (d3sA(ndim, ndim, ndim, 4, n_active))
      allocate (d4sA(ndim, ndim, ndim, ndim, 4, n_active))
      S = 0.0_wp; dS = 0.0_wp; d2S = 0.0_wp

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = ssd_system%sqrt_k3f0_arr(i)
         e_val(1) = base_sq*base
         e_val(2) = base*sqrt_base
         e_val(3) = base
         e_val(4) = base_sq

         f1 = ssd_system%f1_r_arr(:, i)
         f2 = ssd_system%f2_rr_arr(:, :, i)
         f3 = ssd_system%f3_rrr_arr(:, :, :, i)
         f4 = ssd_system%f4_rrrr_arr(:, :, :, :, i)

         do kk = 1, 4
            call atom_exp_derivs(e_val(kk), lambda(kk), f1, f2, f3, f4, v0, v1, v2, v3, v4)
            sA(kk, i) = v0
            dsA(:, kk, i) = v1
            d2sA(:, :, kk, i) = v2
            d3sA(:, :, :, kk, i) = v3
            d4sA(:, :, :, :, kk, i) = v4
            S(kk) = S(kk) + v0
            dS(:, kk) = dS(:, kk) + v1
            d2S(:, :, kk) = d2S(:, :, kk) + v2
         end do
      end do

      C1 = s_1
      if (n_active >= 2) C1 = C1 - 0.5_wp*s_2
      if (n_active >= 3) C1 = C1 + s_3/3.0_wp
      half_s_2 = 0.5_wp*s_2
      s_36 = s_3/6.0_wp
      neg_half_s_3 = -0.5_wp*s_3

      do iB = 1, n_active
         do axb = 1, ndim
            do kk = 1, 4
               RSb(kk) = -dsA(axb, kk, iB)
               RdSb(:, kk) = -d2sA(:, axb, kk, iB)
               Rd2Sb(:, :, kk) = -d3sA(:, :, axb, kk, iB)
            end do
            do iA = 1, n_active
               do axa = 1, ndim
                  do kk = 1, 4
                     RSa(kk) = -dsA(axa, kk, iA)
                     RdSa(:, kk) = -d2sA(:, axa, kk, iA)
                     Rd2Sa(:, :, kk) = -d3sA(:, :, axa, kk, iA)
                     ! Second R-derivative ( d^2 S_k/ dR_A dR_B) is nonzero only when iA==iB,
                     ! and equals +d^2 sA at (axa, axb) restricted to that atom:
                     !    d^2 S_k/ dR_A dR_B = d/ dR_B [-(dsA(axa,kk,iA))]
                     !                  = +(d^2 sA(axb, axa, kk, iA)) when iA==iB
                     !                  = 0 otherwise
                     if (iA == iB) then
                        RRS(kk) = d2sA(axb, axa, kk, iA)
                        RRdS(:, kk) = d3sA(:, axb, axa, kk, iA)
                        RRd2S(:, :, kk) = d4sA(:, :, axb, axa, kk, iA)
                     else
                        RRS(kk) = 0.0_wp
                        RRdS(:, kk) = 0.0_wp
                        RRd2S(:, :, kk) = 0.0_wp
                     end if
                  end do

                  ! Assemble d^2 Z/ dR_A dR_B and r-derivatives via Leibniz on the polynomial:
                  !   Z = C1*S1 + (s_2/2)*S2^2 + (s_3/6)*S3^3 - (s_3/2)*S3*T3
                  rrZ_val = 0.0_wp
                  rrGZ = 0.0_wp
                  rrHZ = 0.0_wp

                  ! Linear: C1 * S1. d^2 (C1*S1) = C1 * RRS1.
                  rrZ_val = rrZ_val + C1*RRS(1)
                  rrGZ = rrGZ + C1*RRdS(:, 1)
                  rrHZ = rrHZ + C1*RRd2S(:, :, 1)

                  ! Quadratic (s_2/2) S2^2:
                  !    d^2 (S2^2) = 2 (RRS2)*S2 + 2 (RS2_A)*(RS2_B)
                  ! For the r-derivative tensors, apply Leibniz with each pairing of
                  ! R_A/R_B replacements onto a single factor.
                  if (n_active >= 2) then
                     ! Term: 2 * (S2 * RRS2_AB) ... i.e., both R-derivatives onto one factor
                     ! There are 2 ways to pick which of the two S2-factors carries
                     ! the second-order RR-derivative (rest is bare S2):
                     call accum_uv_0(rrZ_val, 2.0_wp*half_s_2, RRS(2), S(2))
                     call accum_uv_1(rrGZ, 2.0_wp*half_s_2, RRS(2), RRdS(:, 2), &
                                     S(2), dS(:, 2))
                     call accum_uv_2(rrHZ, 2.0_wp*half_s_2, RRS(2), RRdS(:, 2), RRd2S(:, :, 2), &
                                     S(2), dS(:, 2), d2S(:, :, 2))
                     ! Term: 2 * (RS2_A * RS2_B): one R-derivative on each factor (2 orderings,
                     ! both give same product by symmetry).
                     call accum_uv_0(rrZ_val, 2.0_wp*half_s_2, RSa(2), RSb(2))
                     call accum_uv_1(rrGZ, 2.0_wp*half_s_2, RSa(2), RdSa(:, 2), &
                                     RSb(2), RdSb(:, 2))
                     call accum_uv_2(rrHZ, 2.0_wp*half_s_2, RSa(2), RdSa(:, 2), Rd2Sa(:, :, 2), &
                                     RSb(2), RdSb(:, 2), Rd2Sb(:, :, 2))
                  end if

                  if (n_active >= 3) then
                     ! Cubic S3^3 -> d^2 (S3^3) = 3*S3^2*RRS3 + 6*S3*RS3_A*RS3_B
                     call accum_RR_S3_cubed(rrHZ, rrGZ, rrZ_val, s_36, &
                                            S(3), dS(:, 3), d2S(:, :, 3), &
                                            RSa(3), RdSa(:, 3), Rd2Sa(:, :, 3), &
                                            RSb(3), RdSb(:, 3), Rd2Sb(:, :, 3), &
                                            RRS(3), RRdS(:, 3), RRd2S(:, :, 3))

                     ! Mixed -(s_3/2) S3*T3:
                     !  d^2 (S3*T3) = RRS3 * T3 + S3 * RRT3 + RS3_A * RT3_B + RS3_B * RT3_A
                     call accum_uv_0(rrZ_val, neg_half_s_3, RRS(3), S(4))
                     call accum_uv_0(rrZ_val, neg_half_s_3, S(3), RRS(4))
                     call accum_uv_0(rrZ_val, neg_half_s_3, RSa(3), RSb(4))
                     call accum_uv_0(rrZ_val, neg_half_s_3, RSb(3), RSa(4))

                     call accum_uv_1(rrGZ, neg_half_s_3, RRS(3), RRdS(:, 3), S(4), dS(:, 4))
                     call accum_uv_1(rrGZ, neg_half_s_3, S(3), dS(:, 3), RRS(4), RRdS(:, 4))
                     call accum_uv_1(rrGZ, neg_half_s_3, RSa(3), RdSa(:, 3), RSb(4), RdSb(:, 4))
                     call accum_uv_1(rrGZ, neg_half_s_3, RSb(3), RdSb(:, 3), RSa(4), RdSa(:, 4))

                     call accum_uv_2(rrHZ, neg_half_s_3, RRS(3), RRdS(:, 3), RRd2S(:, :, 3), &
                                     S(4), dS(:, 4), d2S(:, :, 4))
                     call accum_uv_2(rrHZ, neg_half_s_3, S(3), dS(:, 3), d2S(:, :, 3), &
                                     RRS(4), RRdS(:, 4), RRd2S(:, :, 4))
                     call accum_uv_2(rrHZ, neg_half_s_3, RSa(3), RdSa(:, 3), Rd2Sa(:, :, 3), &
                                     RSb(4), RdSb(:, 4), Rd2Sb(:, :, 4))
                     call accum_uv_2(rrHZ, neg_half_s_3, RSb(3), RdSb(:, 3), Rd2Sb(:, :, 3), &
                                     RSa(4), RdSa(:, 4), Rd2Sa(:, :, 4))
                  end if

                  d2Z_rArB(axa, iA, axb, iB) = rrZ_val
                  d2gradZ_rArB(:, axa, iA, axb, iB) = rrGZ
                  d2hessZ_rArB(:, :, axa, iA, axb, iB) = rrHZ
               end do
            end do
         end do
      end do

      deallocate (sA, dsA, d2sA, d3sA, d4sA)
   end subroutine compute_lsf_zRR_derivs_rr_screened

!* ================================================================================= *!
!*                     Nuclear derivative routines                                    *!
!* ================================================================================= *!

!> Compute LSF nuclear gradient (screened version)
   subroutine compute_lsf_f1_rA_screened(ssd_system, k, s_1, s_2, s_3, Z, result)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3, Z
      real(wp), allocatable, intent(out) :: result(:, :)

      integer :: i, n_active
      real(wp) :: base, base_sq, sqrt_base, e1, e2, e3, e3_sq
      real(wp) :: invZ, sum_e2, sum_e3, sum_e3_sq
      real(wp) :: sum_e3_excl, sum_e3_sq_excl, pair_sum_excl, w_i

      n_active = ssd_system%n_active
      allocate (result(ndim, n_active))
      result = 0.0_wp
      if (n_active == 0) return

      invZ = 1.0_wp/Z

      sum_e2 = 0.0_wp
      sum_e3 = 0.0_wp
      sum_e3_sq = 0.0_wp
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         sum_e2 = sum_e2 + base*ssd_system%sqrt_k3f0_arr(i)
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq
      end do

      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = ssd_system%sqrt_k3f0_arr(i)
         e1 = base_sq*base
         e2 = base*sqrt_base
         e3 = base
         e3_sq = base_sq
         w_i = s_1*e1*invZ
         if (n_active >= 2) w_i = w_i + 0.5_wp*s_2*e2*(sum_e2 - e2)*invZ
         if (n_active >= 3) then
            sum_e3_excl = sum_e3 - e3
            sum_e3_sq_excl = sum_e3_sq - e3_sq
            pair_sum_excl = 0.5_wp*(sum_e3_excl*sum_e3_excl - sum_e3_sq_excl)
            w_i = w_i + (s_3/3.0_wp)*e3*pair_sum_excl*invZ
         end if
         result(:, i) = -w_i*ssd_system%f1_r_arr(:, i)
      end do
   end subroutine compute_lsf_f1_rA_screened

!> Compute LSF mixed second derivative d^2f/(dr dR_A) (screened version)
!>
!> Uses O(n_active) algorithm by decomposing the pairwise coupling matrix C_{w,p} = d(M_w)/d(R_p) / f1_r(:,p)
!> into separable terms, so the inner sum sigma_p = sum_w f1_r(:,w) * C_{w,p} is expressible via precomputed
!> global vector sums F_e2, F_e3, F_e3sq with per-atom self-exclusion
   subroutine compute_lsf_f2_r_rA_screened(ssd_system, k, s_1, s_2, s_3, Z, result)
      type(moist_cavity_drop_lsf_svdw_ssd_type), intent(in) :: ssd_system
      real(wp), intent(in) :: k, s_1, s_2, s_3, Z
      real(wp), allocatable, intent(out) :: result(:, :, :)

      integer :: i, j, ii, n_active
      real(wp) :: lambda1, lambda2, lambda3
      real(wp) :: base, base_sq, e1, e2
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: sum_e3_excl, sum_e3_sq_excl, pair_sum_excl
      real(wp) :: invZ, M_i, inv_x, hess_col(ndim)
      !> Global weighted gradient sums: F_ek(m) = sum_i e_k_i * f1_r(m, i)
      real(wp) :: F_e2(ndim), F_e3(ndim), F_e3sq(ndim)
      !> Per-atom partition weights and their spatial moment
      real(wp), allocatable :: W_vals(:)
      real(wp) :: gradW(ndim)
      !> Per-atom assembly scalars for the separable decomposition:
      !>  C_diag = C_{p,p} (diagonal coupling coefficient)
      !>  g1, g2, g3 = off-diagonal coupling coefficients for each F vector
      !>  C_local = C_diag - g1*e2_p - g2*e3_p + g3*e3_sq_p (self-exclusion)
      real(wp) :: C_diag, C_local, g1, g2, g3
      real(wp) :: sigma(ndim), S_p(ndim)
      real(wp) :: grad_d(ndim)

      n_active = ssd_system%n_active
      allocate (result(ndim, ndim, n_active))
      result = 0.0_wp
      if (n_active == 0) return

      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      invZ = 1.0_wp/Z

      ! === Pass 1: Scalar sums and weighted gradient sums ===
      ! 3 scalar + 3x3 vector accumulators = 12 values (fits in registers)
      sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      F_e2 = 0.0_wp; F_e3 = 0.0_wp; F_e3sq = 0.0_wp
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)
         sum_e2 = sum_e2 + e2
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq
         F_e2 = F_e2 + e2*grad_d
         F_e3 = F_e3 + base*grad_d
         F_e3sq = F_e3sq + base_sq*grad_d
      end do

      ! === Pass 2: Compute W_vals and gradW ===
      allocate (W_vals(n_active))
      gradW = 0.0_wp
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         M_i = s_1*e1
         if (n_active >= 2) M_i = M_i + 0.5_wp*s_2*e2*(sum_e2 - e2)
         if (n_active >= 3) then
            sum_e3_excl = sum_e3 - base
            sum_e3_sq_excl = sum_e3_sq - base_sq
            pair_sum_excl = 0.5_wp*(sum_e3_excl*sum_e3_excl - sum_e3_sq_excl)
            M_i = M_i + (s_3/3.0_wp)*base*pair_sum_excl
         end if
         W_vals(i) = M_i*invZ
         gradW = gradW + W_vals(i)*ssd_system%f1_r_arr(:, i)
      end do

      ! === Pass 3: Per-atom result in O(1) per atom ===
      ! For each perturbed atom p, the sum sigma_p(m) = sum_w f1_r(m,w) * C_{w,p}
      ! decomposes via separable C_{w,p} into precomputed F vectors with self-
      ! exclusion. Result: outer(S_p, f1_r(:,p)) - W_p * f2_rr(:,:,p).
      do i = 1, n_active
         base = ssd_system%k3f0_arr(i)
         base_sq = base*base
         e1 = base_sq*base
         e2 = base*ssd_system%sqrt_k3f0_arr(i)
         grad_d = ssd_system%f1_r_arr(:, i)

         ! Diagonal coupling coefficient C_{p,p}
         C_diag = s_1*lambda1*e1

         ! Off-diagonal coupling coefficients for the separable decomposition:
         ! C_{w,p} = g1*e2_w + g2*e3_w - g3*e3_sq_w  (for w != p)
         g1 = 0.0_wp; g2 = 0.0_wp; g3 = 0.0_wp
         if (n_active >= 2) then
            C_diag = C_diag + 0.5_wp*s_2*lambda2*e2*(sum_e2 - e2)
            g1 = 0.5_wp*s_2*lambda2*e2
         end if
         if (n_active >= 3) then
            sum_e3_excl = sum_e3 - base
            sum_e3_sq_excl = sum_e3_sq - base_sq
            pair_sum_excl = 0.5_wp*(sum_e3_excl*sum_e3_excl - sum_e3_sq_excl)
            C_diag = C_diag + (s_3/3.0_wp)*lambda3*base*pair_sum_excl
            g2 = (s_3/3.0_wp)*lambda3*(sum_e3*base - base_sq)
            g3 = (s_3/3.0_wp)*lambda3*base
         end if

         ! sigma = C_local * f1_r(:,p) + g1*F_e2 + g2*F_e3 - g3*F_e3sq
         ! where C_local absorbs the diagonal and self-exclusion correction
         C_local = C_diag - g1*e2 - g2*base + g3*base_sq
         sigma = C_local*grad_d + g1*F_e2 + g2*F_e3 - g3*F_e3sq

         ! S_p = sigma/Z - k * W_p * gradW
         S_p = sigma*invZ - k*W_vals(i)*gradW

         ! result = outer(S_p, f1_r(:,p)) - W_p * f2_rr(:,:,p), reconstructing
         ! each f2_rr column = (I - n n^T)/x from f1_r (= grad_d) and inv_x
         inv_x = ssd_system%inv_x_arr(i)
         do j = 1, ndim
            do ii = 1, ndim
               hess_col(ii) = -grad_d(ii)*grad_d(j)*inv_x
            end do
            hess_col(j) = hess_col(j) + inv_x
            result(:, j, i) = S_p*grad_d(j) - W_vals(i)*hess_col
         end do
      end do
   end subroutine compute_lsf_f2_r_rA_screened

!> Compute pure nuclear Hessian of LSF: d^2 S / dR_{A,s_1} dR_{B,s_2} (screened).
!>
!> Formula:  lsf2_rArB(s_1,A,s_2,B) = D_{AB} * f1_r(s_1,A) * f1_r(s_2,B)
!>                                        + delta_{AB} * W_A * f2_rr(s_1,s_2,A)
!> where D_{AB} = -(1/k) * (d^2 Z / dd_A dd_B) / Z + k * W_A * W_B
!> and W_A = M_A / Z are Boltzmann partition weights.
!>
!> @param[in]  self        LSF primitive with cached parameters (uses self%ssd_system)
!> @param[out] result Nuclear Hessian (3, nsph, 3, nsph) - sparse, only active blocks filled
   subroutine lsf_f2_rArB_screened(self, result)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), allocatable, intent(out) :: result(:, :, :, :)

      real(wp) :: k_val, s_1_c, s_2_c, s_3_c
      real(wp) :: Z, invZ
      integer :: n_active, i, j, iA, iB
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: base_A, base_B, e2_A, e2_B, e1_A
      real(wp) :: base_sq_A
      real(wp) :: sum_e3_excl_A, sum_e3_sq_excl_A, pair_excl_A
      real(wp) :: d2Z_dd, D_AB
      real(wp), allocatable :: W_arr(:), M_arr(:)

      n_active = self%ssd_system%n_active
      allocate (result(ndim, self%ncenters, ndim, self%ncenters))
      result = 0.0_wp
      if (n_active == 0) return

      k_val = self%param%blend_k
      s_1_c = self%param%blend_1b
      s_2_c = self%param%blend_2b
      s_3_c = self%param%blend_3b

      call compute_lsf_z0_screened(self%ssd_system, k_val, s_1_c, s_2_c, s_3_c, Z)
      invZ = 1.0_wp/Z

      ! Compute partition weights W_A = M_A / Z
      allocate (W_arr(n_active), M_arr(n_active))
      sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      do i = 1, n_active
         base_A = self%ssd_system%k3f0_arr(i)
         sum_e2 = sum_e2 + base_A*self%ssd_system%sqrt_k3f0_arr(i)
         sum_e3 = sum_e3 + base_A
         sum_e3_sq = sum_e3_sq + base_A*base_A
      end do

      do i = 1, n_active
         base_A = self%ssd_system%k3f0_arr(i)
         base_sq_A = base_A*base_A
         e1_A = base_sq_A*base_A
         e2_A = base_A*self%ssd_system%sqrt_k3f0_arr(i)
         M_arr(i) = s_1_c*e1_A
         if (n_active >= 2) M_arr(i) = M_arr(i) + 0.5_wp*s_2_c*e2_A*(sum_e2 - e2_A)
         if (n_active >= 3) then
            sum_e3_excl_A = sum_e3 - base_A
            sum_e3_sq_excl_A = sum_e3_sq - base_sq_A
            pair_excl_A = 0.5_wp*(sum_e3_excl_A**2 - sum_e3_sq_excl_A)
            M_arr(i) = M_arr(i) + (s_3_c/3.0_wp)*base_A*pair_excl_A
         end if
         W_arr(i) = M_arr(i)*invZ
      end do

      ! Assemble result for each (A, B) pair of active atoms
      do iA = 1, n_active
         i = self%ssd_system%atom_indices(iA)
         base_A = self%ssd_system%k3f0_arr(iA)
         base_sq_A = base_A*base_A
         e1_A = base_sq_A*base_A
         e2_A = base_A*self%ssd_system%sqrt_k3f0_arr(iA)

         do iB = 1, n_active
            j = self%ssd_system%atom_indices(iB)
            base_B = self%ssd_system%k3f0_arr(iB)
            e2_B = base_B*self%ssd_system%sqrt_k3f0_arr(iB)

            ! Compute d^2 Z / dd_A dd_B
            if (iA == iB) then
               sum_e3_excl_A = sum_e3 - base_A
               sum_e3_sq_excl_A = sum_e3_sq - base_sq_A
               pair_excl_A = 0.5_wp*(sum_e3_excl_A**2 - sum_e3_sq_excl_A)
               d2Z_dd = s_1_c*k_val**2*e1_A &
                        + s_2_c*(0.5_wp*k_val)**2*e2_A*(sum_e2 - e2_A) &
                        + s_3_c*(k_val/3.0_wp)**2*base_A*pair_excl_A
            else
               ! Off-diagonal: d^2Z/dd_A dd_B
               d2Z_dd = s_2_c*(0.5_wp*k_val)**2*e2_A*e2_B
               if (n_active >= 3) then
                  d2Z_dd = d2Z_dd + s_3_c*(k_val/3.0_wp)**2 &
                           *base_A*base_B*(sum_e3 - base_A - base_B)
               end if
            end if

            ! D_{AB} = -(1/k) * (d^2Z/dd_A dd_B) / Z + k * W_A * W_B
            D_AB = -(1.0_wp/k_val)*d2Z_dd*invZ + k_val*W_arr(iA)*W_arr(iB)

            ! result(s_1, A, s_2, B) = D_AB * f1_r(s_1,A) * f1_r(s_2,B)
            !                            + delta_{AB} * W_A * f2_rr(s_1,s_2,A)
            call assemble_f2_rArB_block(result(:, i, :, j), &
                                        D_AB, self%ssd_system%f1_r_arr(:, iA), self%ssd_system%f1_r_arr(:, iB), &
                                        W_arr(iA), self%ssd_system%f2_rr_arr(:, :, iA), iA == iB)
         end do
      end do

      deallocate (W_arr, M_arr)

   contains

      pure subroutine assemble_f2_rArB_block(block, D, gA, gB, WA, hA, is_diag)
         real(wp), intent(inout) :: block(ndim, ndim)
         real(wp), intent(in) :: D, gA(ndim), gB(ndim), WA, hA(ndim, ndim)
         logical, intent(in) :: is_diag
         integer :: a, b
         do b = 1, ndim
            do a = 1, ndim
               block(a, b) = D*gA(a)*gB(b)
            end do
         end do
         if (is_diag) then
            do b = 1, ndim
               do a = 1, ndim
                  block(a, b) = block(a, b) + WA*hA(a, b)
               end do
            end do
         end if
      end subroutine assemble_f2_rArB_block

   end subroutine lsf_f2_rArB_screened

!> Compute mixed third derivative d^3 S / dr_j dR_{A,s_1} dR_{B,s_2} (screened).
!>
!> Uses the Z-level quotient rule decomposition. Since S = -(1/k)*ln(Z), the
!> third derivative is:
!>   lsf3 = -(1/k)*[d2gZ/Z - dgZ_A*dZ_B/Z^2 - dgZ_B*dZ_A/Z^2
!>                   - gZ*d2Z_AB/Z^2 + 2*gZ*dZ_A*dZ_B/Z^3]
!> where gZ = dZ/dr, dgZ_A = d^2Z/(dr dR_A), dZ_A = dZ/dR_A, etc.
!> All Z-level derivatives are computed from the SSD exponential sums.
!>
!> @param[in]  self        LSF primitive (uses self%ssd_system, needs max_deriv >= 3)
!> @param[in]  lsf1_rA     First nuclear derivatives (3, n_active), unused but kept for API
!> @param[in]  lsf2_r_rA   Mixed second derivatives (3, 3, n_active), unused but kept for API
!> @param[out] result      (j, s_1, iA, s_2, iB) with active indices iA, iB = 1..n_active
   subroutine lsf_f3_r_rArB_screened(self, lsf1_rA, lsf2_r_rA, result)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(in) :: lsf1_rA(:, :)
      real(wp), intent(in) :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: result(:, :, :, :, :)

      real(wp) :: k_val, s_1_c, s_2_c, s_3_c, k2, k3
      real(wp) :: Z, invZ, invZ2, invZ3, inv_k
      integer :: n_active, iA, iB, ax_a, ax_b
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: base_A, base_B, bsq_A, bsq_B, e1_A, e2_A, e2_B
      real(wp) :: s_excl, ssq_excl, pair_excl
      real(wp) :: M_val, C_local, g1_v, g2_v, g3_v

      !> Precomputed sums
      real(wp) :: F_e2(ndim), F_e3(ndim), F_e3sq(ndim), gradZ_vec(ndim)

      !> Per-atom arrays
      real(wp), allocatable :: W_arr(:), Z1_arr(:), sigma_Z(:, :)

      !> Per-pair arrays
      real(wp), allocatable :: Z2_mat(:, :), tau_vec(:, :, :)

      !> Intermediates for assembly
      real(wp) :: f1_A(ndim), f1_B(ndim), f2_A(ndim, ndim), f2_B(ndim, ndim)
      real(wp) :: dgZ_rA(ndim), dgZ_rB(ndim), d2gZ(ndim)
      real(wp) :: dZ_rA_val, dZ_rB_val, d2Z_rArB_val
      real(wp) :: Z_AAA, Z_AAB, Fexcl_e2(ndim), Fexcl_e3(ndim), Fexcl_e3sq(ndim)

      n_active = self%ssd_system%n_active
      allocate (result(ndim, ndim, n_active, ndim, n_active))
      result = 0.0_wp
      if (n_active == 0) return

      k_val = self%param%blend_k
      s_1_c = self%param%blend_1b
      s_2_c = self%param%blend_2b
      s_3_c = self%param%blend_3b
      k2 = 0.5_wp*k_val
      k3 = k_val/3.0_wp
      inv_k = 1.0_wp/k_val

      call compute_lsf_z0_screened(self%ssd_system, k_val, s_1_c, s_2_c, s_3_c, Z)
      invZ = 1.0_wp/Z
      invZ2 = invZ*invZ
      invZ3 = invZ2*invZ

      ! === Pass 0: Exponential sums and weighted gradient sums ===
      sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      F_e2 = 0.0_wp; F_e3 = 0.0_wp; F_e3sq = 0.0_wp

      do iA = 1, n_active
         base_A = self%ssd_system%k3f0_arr(iA)
         bsq_A = base_A*base_A
         e2_A = base_A*self%ssd_system%sqrt_k3f0_arr(iA)
         f1_A = self%ssd_system%f1_r_arr(:, iA)
         sum_e2 = sum_e2 + e2_A
         sum_e3 = sum_e3 + base_A
         sum_e3_sq = sum_e3_sq + bsq_A
         F_e2 = F_e2 + e2_A*f1_A
         F_e3 = F_e3 + base_A*f1_A
         F_e3sq = F_e3sq + bsq_A*f1_A
      end do

      ! === Pass 1: Per-atom W, Z_A, sigma_Z, gradZ ===
      allocate (W_arr(n_active), Z1_arr(n_active), sigma_Z(ndim, n_active))
      gradZ_vec = 0.0_wp

      do iA = 1, n_active
         base_A = self%ssd_system%k3f0_arr(iA)
         bsq_A = base_A*base_A
         e1_A = bsq_A*base_A
         e2_A = base_A*self%ssd_system%sqrt_k3f0_arr(iA)
         f1_A = self%ssd_system%f1_r_arr(:, iA)

         M_val = s_1_c*e1_A
         if (n_active >= 2) M_val = M_val + 0.5_wp*s_2_c*e2_A*(sum_e2 - e2_A)
         if (n_active >= 3) then
            s_excl = sum_e3 - base_A
            ssq_excl = sum_e3_sq - bsq_A
            pair_excl = 0.5_wp*(s_excl**2 - ssq_excl)
            M_val = M_val + (s_3_c/3.0_wp)*base_A*pair_excl
         end if
         W_arr(iA) = M_val*invZ
         Z1_arr(iA) = -k_val*M_val

         gradZ_vec = gradZ_vec + Z1_arr(iA)*f1_A

         ! sigma_Z_A(j) = sum_J Z_{JA}*f1_r(j,J) = k * sigma_code_A(j)
         C_local = s_1_c*k_val*e1_A
         g1_v = 0.0_wp; g2_v = 0.0_wp; g3_v = 0.0_wp
         if (n_active >= 2) then
            C_local = C_local + 0.5_wp*s_2_c*k2*e2_A*(sum_e2 - e2_A)
            g1_v = 0.5_wp*s_2_c*k2*e2_A
            C_local = C_local - g1_v*e2_A
         end if
         if (n_active >= 3) then
            s_excl = sum_e3 - base_A
            ssq_excl = sum_e3_sq - bsq_A
            pair_excl = 0.5_wp*(s_excl**2 - ssq_excl)
            C_local = C_local + (s_3_c/3.0_wp)*k3*base_A*pair_excl
            g2_v = (s_3_c/3.0_wp)*k3*(sum_e3*base_A - bsq_A)
            g3_v = (s_3_c/3.0_wp)*k3*base_A
            C_local = C_local - g2_v*base_A + g3_v*bsq_A
         end if
         sigma_Z(:, iA) = k_val*(C_local*f1_A + g1_v*F_e2 + g2_v*F_e3 - g3_v*F_e3sq)
      end do

      ! === Pass 2: Per-pair Z_{AB} and tau_{AB} ===
      allocate (Z2_mat(n_active, n_active), tau_vec(ndim, n_active, n_active))

      do iB = 1, n_active
         base_B = self%ssd_system%k3f0_arr(iB)
         bsq_B = base_B*base_B
         e2_B = base_B*self%ssd_system%sqrt_k3f0_arr(iB)
         f1_B = self%ssd_system%f1_r_arr(:, iB)

         do iA = 1, n_active
            base_A = self%ssd_system%k3f0_arr(iA)
            bsq_A = base_A*base_A
            e1_A = bsq_A*base_A
            e2_A = base_A*self%ssd_system%sqrt_k3f0_arr(iA)
            f1_A = self%ssd_system%f1_r_arr(:, iA)

            if (iA == iB) then
               ! --- Diagonal Z_{AA} ---
               s_excl = sum_e3 - base_A
               ssq_excl = sum_e3_sq - bsq_A
               pair_excl = 0.5_wp*(s_excl**2 - ssq_excl)
               Z2_mat(iA, iB) = s_1_c*k_val**2*e1_A &
                                + s_2_c*k2**2*e2_A*(sum_e2 - e2_A) &
                                + s_3_c*k3**2*base_A*pair_excl

               ! --- Diagonal Z_{AAA} and tau_{AA} ---
               Z_AAA = -s_1_c*k_val**3*e1_A &
                       - s_2_c*k2**3*e2_A*(sum_e2 - e2_A) &
                       - s_3_c*k3**3*base_A*pair_excl

               Fexcl_e2 = F_e2 - e2_A*f1_A
               Fexcl_e3 = F_e3 - base_A*f1_A
               Fexcl_e3sq = F_e3sq - bsq_A*f1_A

               tau_vec(:, iA, iB) = Z_AAA*f1_A &
                                    - s_2_c*k2**3*e2_A*Fexcl_e2
               if (n_active >= 3) then
                  tau_vec(:, iA, iB) = tau_vec(:, iA, iB) &
                                       - s_3_c*k3**3*base_A*(s_excl*Fexcl_e3 - Fexcl_e3sq)
               end if
            else
               ! --- Off-diagonal Z_{AB} ---
               Z2_mat(iA, iB) = s_2_c*k2**2*e2_A*e2_B
               if (n_active >= 3) then
                  Z2_mat(iA, iB) = Z2_mat(iA, iB) &
                                   + s_3_c*k3**2*base_A*base_B*(sum_e3 - base_A - base_B)
               end if

               ! --- Off-diagonal tau_{AB} ---
               ! Z_{AAB} = Z_{ABB} for A != B (symmetric product structure)
               Z_AAB = -s_2_c*k2**3*e2_A*e2_B
               if (n_active >= 3) then
                  Z_AAB = Z_AAB - s_3_c*k3**3*base_A*base_B &
                          *(sum_e3 - base_A - base_B)
               end if
               tau_vec(:, iA, iB) = Z_AAB*(f1_A + f1_B)
               if (n_active >= 3) then
                  tau_vec(:, iA, iB) = tau_vec(:, iA, iB) &
                                       - s_3_c*k3**3*base_A*base_B &
                                       *(F_e3 - base_A*f1_A - base_B*f1_B)
               end if
            end if
         end do
      end do

      ! === Pass 3: Assembly via quotient rule ===
      ! lsf3(j) = -(1/k)*[d2gZ(j)/Z - dgZ_A(j)*dZ_B/Z^2 - dgZ_B(j)*dZ_A/Z^2
      !                    - gZ(j)*d2Z_AB/Z^2 + 2*gZ(j)*dZ_A*dZ_B/Z^3]
      do iB = 1, n_active
         f1_B = self%ssd_system%f1_r_arr(:, iB)
         f2_B = self%ssd_system%f2_rr_arr(:, :, iB)

         do iA = 1, n_active
            f1_A = self%ssd_system%f1_r_arr(:, iA)
            f2_A = self%ssd_system%f2_rr_arr(:, :, iA)

            do ax_b = 1, ndim
               ! dZ/dR_{B,s_2} = Z_B * (-f1_r(s_2,B))
               dZ_rB_val = -Z1_arr(iB)*f1_B(ax_b)

               ! dgradZ/dR_{B,s_2}(j) = -f1_r(s_2,B)*sigma_Z_B(j) - Z_B*f2_rr(j,s_2,B)
               dgZ_rB = -f1_B(ax_b)*sigma_Z(:, iB) - Z1_arr(iB)*f2_B(:, ax_b)

               do ax_a = 1, ndim
                  ! dZ/dR_{A,s_1}
                  dZ_rA_val = -Z1_arr(iA)*f1_A(ax_a)

                  ! dgradZ/dR_{A,s_1}(j)
                  dgZ_rA = -f1_A(ax_a)*sigma_Z(:, iA) - Z1_arr(iA)*f2_A(:, ax_a)

                  ! d^2Z / (dR_A dR_B)
                  d2Z_rArB_val = Z2_mat(iA, iB)*f1_A(ax_a)*f1_B(ax_b)
                  if (iA == iB) then
                     d2Z_rArB_val = d2Z_rArB_val + Z1_arr(iA)*f2_A(ax_a, ax_b)
                  end if

                  ! d^2(gradZ) / (dR_A dR_B)(j)
                  d2gZ = f1_A(ax_a)*f1_B(ax_b)*tau_vec(:, iA, iB) &
                         + f1_A(ax_a)*Z2_mat(iA, iB)*f2_B(:, ax_b) &
                         + Z2_mat(iA, iB)*f1_B(ax_b)*f2_A(:, ax_a)
                  if (iA == iB) then
                     d2gZ = d2gZ + f2_A(ax_a, ax_b)*sigma_Z(:, iA) &
                            + Z1_arr(iA)*self%ssd_system%f3_rrr_arr(:, ax_a, ax_b, iA)
                  end if

                  ! Assemble via quotient rule
                  result(:, ax_a, iA, ax_b, iB) = -inv_k*( &
                                                  d2gZ*invZ &
                                                  - dgZ_rA*dZ_rB_val*invZ2 &
                                                  - dgZ_rB*dZ_rA_val*invZ2 &
                                                  - gradZ_vec*d2Z_rArB_val*invZ2 &
                                                  + 2.0_wp*gradZ_vec*dZ_rA_val*dZ_rB_val*invZ3)
               end do
            end do
         end do
      end do

      deallocate (W_arr, Z1_arr, sigma_Z, Z2_mat, tau_vec)
   end subroutine lsf_f3_r_rArB_screened

!* ================================================================================= *!
!*                     POU (Partition of Unity) routines                             *!
!* ================================================================================= *!

!> These are not needed anymore currently, but may be useful again in the future..

!> Compute partition-of-unity weight and spatial derivatives for a single owner (screened)
   subroutine lsf_pou_f012_r_screened(self, owner_id, weight, dweight_r, d2weight_rr)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      integer, intent(in) :: owner_id
      real(wp), intent(out) :: weight
      real(wp), intent(out), optional :: dweight_r(:)
      real(wp), intent(out), optional :: d2weight_rr(:, :)

      real(wp) :: k, s_1, s_2, s_3
      real(wp) :: lambda1, lambda2, lambda3, lambda3_sq
      real(wp) :: lambda1_sq, lambda2_sq, lambda3_sq_sq
      real(wp) :: Z, invZ, invZ2, invZ3
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: e1_owner, e2_owner, e3_owner, e3_sq_owner
      real(wp) :: sum_e3_excl, sum_e3_sq_excl, pair_sum_excl
      real(wp) :: M_i, grad_M_i(ndim), hess_M_i(ndim, ndim)
      real(wp) :: grad_d(ndim), hess_d(ndim, ndim)
      real(wp) :: grad_e1(ndim), grad_e2(ndim), grad_e3(ndim), grad_e3_sq(ndim)
      real(wp) :: hess_e1(ndim, ndim), hess_e2(ndim, ndim)
      real(wp) :: hess_e3(ndim, ndim), hess_e3_sq(ndim, ndim)
      real(wp) :: dsum_e2(ndim), dsum_e3(ndim), dsum_e3_sq(ndim)
      real(wp) :: d2sum_e2(ndim, ndim), d2sum_e3(ndim, ndim), d2sum_e3_sq(ndim, ndim)
      real(wp) :: dsum_e3_excl(ndim), dsum_e3_sq_excl(ndim)
      real(wp) :: d2sum_e3_excl(ndim, ndim), d2sum_e3_sq_excl(ndim, ndim)
      real(wp) :: dpair_sum_excl(ndim), d2pair_sum_excl(ndim, ndim)
      real(wp) :: base, base_sq, sqrt_base, inv_x
      logical :: need_grad, need_hess, compute_grad
      integer :: i, j, ii, atom, n_active

      need_grad = present(dweight_r)
      need_hess = present(d2weight_rr)
      compute_grad = need_grad .or. need_hess

      weight = 0.0_wp
      if (need_grad) dweight_r = 0.0_wp
      if (need_hess) d2weight_rr = 0.0_wp
      if (self%ncenters == 0) return
      if (owner_id < 1 .or. owner_id > self%ncenters) return
      n_active = self%ssd_system%n_active
      if (n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b
      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      lambda3_sq = 2.0_wp*k/3.0_wp
      lambda1_sq = lambda1*lambda1
      lambda2_sq = lambda2*lambda2
      lambda3_sq_sq = lambda3_sq*lambda3_sq

      sum_e2 = 0.0_wp
      sum_e3 = 0.0_wp
      sum_e3_sq = 0.0_wp
      dsum_e2 = 0.0_wp
      dsum_e3 = 0.0_wp
      dsum_e3_sq = 0.0_wp
      d2sum_e2 = 0.0_wp
      d2sum_e3 = 0.0_wp
      d2sum_e3_sq = 0.0_wp
      e1_owner = 0.0_wp
      e2_owner = 0.0_wp
      e3_owner = 0.0_wp
      e3_sq_owner = 0.0_wp
      grad_e1 = 0.0_wp
      grad_e2 = 0.0_wp
      grad_e3 = 0.0_wp
      grad_e3_sq = 0.0_wp
      hess_e1 = 0.0_wp
      hess_e2 = 0.0_wp
      hess_e3 = 0.0_wp
      hess_e3_sq = 0.0_wp

      do i = 1, n_active
         atom = self%ssd_system%atom_indices(i)
         base = self%ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = self%ssd_system%sqrt_k3f0_arr(i)
         sum_e2 = sum_e2 + base*sqrt_base
         sum_e3 = sum_e3 + base
         sum_e3_sq = sum_e3_sq + base_sq
         if (compute_grad) then
            grad_d = self%ssd_system%f1_r_arr(:, i)
            dsum_e2 = dsum_e2 - lambda2*(base*sqrt_base)*grad_d
            dsum_e3 = dsum_e3 - lambda3*base*grad_d
            dsum_e3_sq = dsum_e3_sq - lambda3_sq*base_sq*grad_d
         end if
         if (need_hess) then
            ! Reconstruct hess_d = (I - n n^T)/x from f1_r (= grad_d) and inv_x
            inv_x = self%ssd_system%inv_x_arr(i)
            do j = 1, ndim
               do ii = 1, ndim
                  hess_d(ii, j) = -grad_d(ii)*grad_d(j)*inv_x
               end do
               hess_d(j, j) = hess_d(j, j) + inv_x
            end do
            d2sum_e2 = d2sum_e2 + (base*sqrt_base)*(lambda2_sq*outer_matrix(grad_d, grad_d) - lambda2*hess_d)
            d2sum_e3 = d2sum_e3 + base*(lambda3*lambda3*outer_matrix(grad_d, grad_d) - lambda3*hess_d)
            d2sum_e3_sq = d2sum_e3_sq + base_sq*(lambda3_sq_sq*outer_matrix(grad_d, grad_d) - lambda3_sq*hess_d)
         end if
         if (atom == owner_id) then
            e1_owner = base_sq*base
            e2_owner = base*sqrt_base
            e3_owner = base
            e3_sq_owner = base_sq
            if (compute_grad) then
               grad_e1 = -lambda1*e1_owner*grad_d
               grad_e2 = -lambda2*e2_owner*grad_d
               grad_e3 = -lambda3*e3_owner*grad_d
               grad_e3_sq = -lambda3_sq*e3_sq_owner*grad_d
            end if
            if (need_hess) then
               hess_e1 = e1_owner*(lambda1_sq*outer_matrix(grad_d, grad_d) - lambda1*hess_d)
               hess_e2 = e2_owner*(lambda2_sq*outer_matrix(grad_d, grad_d) - lambda2*hess_d)
               hess_e3 = e3_owner*(lambda3*lambda3*outer_matrix(grad_d, grad_d) - lambda3*hess_d)
               hess_e3_sq = e3_sq_owner*(lambda3_sq_sq*outer_matrix(grad_d, grad_d) - lambda3_sq*hess_d)
            end if
         end if
      end do

      call compute_lsf_z012_rr_screened(self%ssd_system, k, s_1, s_2, s_3, Z, gradZ, hessZ)
      invZ = 1.0_wp/max(Z, 1.0e-100_wp)
      invZ2 = invZ*invZ; invZ3 = invZ2*invZ

      M_i = s_1*e1_owner
      if (compute_grad) grad_M_i = s_1*grad_e1
      if (need_hess) hess_M_i = s_1*hess_e1

      if (n_active >= 2) then
         M_i = M_i + 0.5_wp*s_2*e2_owner*(sum_e2 - e2_owner)
         if (compute_grad) then
            grad_M_i = grad_M_i + 0.5_wp*s_2*(grad_e2*(sum_e2 - 2.0_wp*e2_owner) + e2_owner*dsum_e2)
         end if
         if (need_hess) then
            hess_M_i = hess_M_i + 0.5_wp*s_2*( &
                       hess_e2*(sum_e2 - 2.0_wp*e2_owner) &
                       + outer_matrix(grad_e2, dsum_e2) &
                       + outer_matrix(dsum_e2 - 2.0_wp*grad_e2, grad_e2) &
                       + e2_owner*d2sum_e2)
         end if
      end if

      if (n_active >= 3) then
         sum_e3_excl = sum_e3 - e3_owner
         sum_e3_sq_excl = sum_e3_sq - e3_sq_owner
         pair_sum_excl = 0.5_wp*(sum_e3_excl*sum_e3_excl - sum_e3_sq_excl)
         M_i = M_i + (s_3/3.0_wp)*e3_owner*pair_sum_excl
         if (compute_grad) then
            dsum_e3_excl = dsum_e3 - grad_e3
            dsum_e3_sq_excl = dsum_e3_sq - grad_e3_sq
            dpair_sum_excl = sum_e3_excl*dsum_e3_excl - 0.5_wp*dsum_e3_sq_excl
            grad_M_i = grad_M_i + (s_3/3.0_wp)*(grad_e3*pair_sum_excl + e3_owner*dpair_sum_excl)
         end if
         if (need_hess) then
            d2sum_e3_excl = d2sum_e3 - hess_e3
            d2sum_e3_sq_excl = d2sum_e3_sq - hess_e3_sq
            d2pair_sum_excl = outer_matrix(dsum_e3_excl, dsum_e3_excl) + sum_e3_excl*d2sum_e3_excl &
                              - 0.5_wp*d2sum_e3_sq_excl
            hess_M_i = hess_M_i + (s_3/3.0_wp)*( &
                       pair_sum_excl*hess_e3 &
                       + outer_matrix(grad_e3, dpair_sum_excl) &
                       + outer_matrix(dpair_sum_excl, grad_e3) &
                       + e3_owner*d2pair_sum_excl)
         end if
      end if

      weight = M_i*invZ
      if (need_grad) dweight_r = grad_M_i*invZ - M_i*gradZ*invZ2
      if (need_hess) then
         d2weight_rr = hess_M_i*invZ &
                       - (outer_matrix(grad_M_i, gradZ) + outer_matrix(gradZ, grad_M_i))*invZ2 &
                       - M_i*hessZ*invZ2 &
                       + 2.0_wp*M_i*outer_matrix(gradZ, gradZ)*invZ3
      end if
   end subroutine lsf_pou_f012_r_screened

!> Compute nuclear derivative of POU spatial gradient for one owner (screened)
   subroutine lsf_pou_f2_r_rA_screened(self, owner_id, d2weight_r_rA)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      integer, intent(in) :: owner_id
      real(wp), intent(out) :: d2weight_r_rA(:, :, :)

      real(wp), allocatable :: e1_val(:), e2_val(:), e3_val(:), e3_sq_val(:)
      real(wp), allocatable :: grad_e1_arr(:, :), grad_e2_arr(:, :), grad_e3_arr(:, :), grad_e3_sq_arr(:, :)
      real(wp), allocatable :: hess_e1_arr(:, :, :), hess_e2_arr(:, :, :), hess_e3_arr(:, :, :), hess_e3_sq_arr(:, :, :)
      real(wp) :: Z, gradZ(ndim), hessZ(ndim, ndim)
      real(wp) :: sum_e2, sum_e3, sum_e3_sq
      real(wp) :: dsum_e2(ndim), dsum_e3(ndim), dsum_e3_sq(ndim)
      real(wp) :: M_i, grad_M_i(ndim)
      real(wp) :: dZ_dRA(ndim), d_gradZ_dRA(ndim, ndim)
      real(wp) :: dM_i_dRA(ndim), d_gradM_i_dRA(ndim, ndim)
      real(wp) :: de1_dRA(ndim), de2_dRA(ndim), de3_dRA(ndim), de3_sq_dRA(ndim)
      real(wp) :: dgrad_e1_dRA(ndim, ndim), dgrad_e2_dRA(ndim, ndim)
      real(wp) :: dgrad_e3_dRA(ndim, ndim), dgrad_e3_sq_dRA(ndim, ndim)
      real(wp) :: de1_owner_dRA(ndim), de2_owner_dRA(ndim), de3_owner_dRA(ndim), de3_sq_owner_dRA(ndim)
      real(wp) :: dgrad_e1_owner_dRA(ndim, ndim), dgrad_e2_owner_dRA(ndim, ndim)
      real(wp) :: dgrad_e3_owner_dRA(ndim, ndim), dgrad_e3_sq_owner_dRA(ndim, ndim)
      real(wp) :: sum_e3_excl, sum_e3_sq_excl, pair_sum_excl
      real(wp) :: dsum_e3_excl(ndim), dsum_e3_sq_excl(ndim), dpair_sum_excl(ndim)
      real(wp) :: dsum_e3_excl_dRA(ndim), dsum_e3_sq_excl_dRA(ndim)
      real(wp) :: d_dsum_e3_excl_dRA(ndim, ndim), d_dsum_e3_sq_excl_dRA(ndim, ndim)
      real(wp) :: dpair_sum_excl_dRA(ndim), d_dpair_sum_excl_dRA(ndim, ndim)
      real(wp) :: k, s_1, s_2_p, s_3
      real(wp) :: lambda1, lambda2, lambda3
      real(wp) :: invZ, invZ2, invZ3
      real(wp) :: lambda1_sq, lambda2_sq, lambda3_sq, lambda3_sq_sq
      real(wp) :: grad_d(ndim), hess_d(ndim, ndim), base, base_sq, sqrt_base, inv_x
      integer :: i, j, ii, atom, atomA, n, n_active

      d2weight_r_rA = 0.0_wp
      n = self%ncenters
      if (n == 0) return
      if (owner_id < 1 .or. owner_id > n) return
      n_active = self%ssd_system%n_active
      if (n_active == 0) return

      k = self%param%blend_k
      s_1 = self%param%blend_1b
      s_2_p = self%param%blend_2b
      s_3 = self%param%blend_3b
      lambda1 = k
      lambda2 = 0.5_wp*k
      lambda3 = k/3.0_wp
      lambda1_sq = lambda1*lambda1
      lambda2_sq = lambda2*lambda2
      lambda3_sq = 2.0_wp*k/3.0_wp
      lambda3_sq_sq = lambda3_sq*lambda3_sq

      allocate (e1_val(n), e2_val(n), e3_val(n), e3_sq_val(n))
      allocate (grad_e1_arr(ndim, n), grad_e2_arr(ndim, n), grad_e3_arr(ndim, n), grad_e3_sq_arr(ndim, n))
      allocate (hess_e1_arr(ndim, ndim, n), hess_e2_arr(ndim, ndim, n), hess_e3_arr(ndim, ndim, n), hess_e3_sq_arr(ndim, ndim, n))
      e1_val = 0.0_wp; e2_val = 0.0_wp; e3_val = 0.0_wp; e3_sq_val = 0.0_wp
      grad_e1_arr = 0.0_wp; grad_e2_arr = 0.0_wp; grad_e3_arr = 0.0_wp; grad_e3_sq_arr = 0.0_wp
      hess_e1_arr = 0.0_wp; hess_e2_arr = 0.0_wp; hess_e3_arr = 0.0_wp; hess_e3_sq_arr = 0.0_wp
      sum_e2 = 0.0_wp; sum_e3 = 0.0_wp; sum_e3_sq = 0.0_wp
      dsum_e2 = 0.0_wp; dsum_e3 = 0.0_wp; dsum_e3_sq = 0.0_wp

      do i = 1, n_active
         atom = self%ssd_system%atom_indices(i)
         base = self%ssd_system%k3f0_arr(i)
         base_sq = base*base
         sqrt_base = self%ssd_system%sqrt_k3f0_arr(i)
         grad_d = self%ssd_system%f1_r_arr(:, i)

         ! Reconstruct hess_d = (I - n n^T)/x from f1_r (= grad_d) and inv_x
         inv_x = self%ssd_system%inv_x_arr(i)
         do j = 1, ndim
            do ii = 1, ndim
               hess_d(ii, j) = -grad_d(ii)*grad_d(j)*inv_x
            end do
            hess_d(j, j) = hess_d(j, j) + inv_x
         end do
         e1_val(atom) = base_sq*base
         e2_val(atom) = base*sqrt_base
         e3_val(atom) = base
         e3_sq_val(atom) = base_sq
         grad_e1_arr(:, atom) = -lambda1*e1_val(atom)*grad_d
         grad_e2_arr(:, atom) = -lambda2*e2_val(atom)*grad_d
         grad_e3_arr(:, atom) = -lambda3*e3_val(atom)*grad_d
         grad_e3_sq_arr(:, atom) = -lambda3_sq*e3_sq_val(atom)*grad_d
         hess_e1_arr(:, :, atom) = e1_val(atom)*(lambda1_sq*outer_matrix(grad_d, grad_d) - lambda1*hess_d)
         hess_e2_arr(:, :, atom) = e2_val(atom)*(lambda2_sq*outer_matrix(grad_d, grad_d) - lambda2*hess_d)
         hess_e3_arr(:, :, atom) = e3_val(atom)*(lambda3*lambda3*outer_matrix(grad_d, grad_d) - lambda3*hess_d)
         hess_e3_sq_arr(:, :, atom) = e3_sq_val(atom)*(lambda3_sq_sq*outer_matrix(grad_d, grad_d) - lambda3_sq*hess_d)
         sum_e2 = sum_e2 + e2_val(atom)
         sum_e3 = sum_e3 + e3_val(atom)
         sum_e3_sq = sum_e3_sq + e3_sq_val(atom)
         dsum_e2 = dsum_e2 + grad_e2_arr(:, atom)
         dsum_e3 = dsum_e3 + grad_e3_arr(:, atom)
         dsum_e3_sq = dsum_e3_sq + grad_e3_sq_arr(:, atom)
      end do

      call compute_lsf_z012_rr_screened(self%ssd_system, k, s_1, s_2_p, s_3, Z, gradZ, hessZ)
      invZ = 1.0_wp/max(Z, 1.0e-100_wp)
      invZ2 = invZ*invZ; invZ3 = invZ2*invZ

      M_i = s_1*e1_val(owner_id)
      grad_M_i = s_1*grad_e1_arr(:, owner_id)
      if (n_active >= 2) then
         M_i = M_i + 0.5_wp*s_2_p*e2_val(owner_id)*(sum_e2 - e2_val(owner_id))
         grad_M_i = grad_M_i + 0.5_wp*s_2_p*( &
                    grad_e2_arr(:, owner_id)*(sum_e2 - 2.0_wp*e2_val(owner_id)) + e2_val(owner_id)*dsum_e2)
      end if
      if (n_active >= 3) then
         sum_e3_excl = sum_e3 - e3_val(owner_id)
         sum_e3_sq_excl = sum_e3_sq - e3_sq_val(owner_id)
         pair_sum_excl = 0.5_wp*(sum_e3_excl*sum_e3_excl - sum_e3_sq_excl)
         dsum_e3_excl = dsum_e3 - grad_e3_arr(:, owner_id)
         dsum_e3_sq_excl = dsum_e3_sq - grad_e3_sq_arr(:, owner_id)
         dpair_sum_excl = sum_e3_excl*dsum_e3_excl - 0.5_wp*dsum_e3_sq_excl
         M_i = M_i + (s_3/3.0_wp)*e3_val(owner_id)*pair_sum_excl
         grad_M_i = grad_M_i + (s_3/3.0_wp)*( &
                    grad_e3_arr(:, owner_id)*pair_sum_excl + e3_val(owner_id)*dpair_sum_excl)
      end if

      do atomA = 1, n
         de1_dRA = -grad_e1_arr(:, atomA)
         de2_dRA = -grad_e2_arr(:, atomA)
         de3_dRA = -grad_e3_arr(:, atomA)
         de3_sq_dRA = -grad_e3_sq_arr(:, atomA)
         dgrad_e1_dRA = -hess_e1_arr(:, :, atomA)
         dgrad_e2_dRA = -hess_e2_arr(:, :, atomA)
         dgrad_e3_dRA = -hess_e3_arr(:, :, atomA)
         dgrad_e3_sq_dRA = -hess_e3_sq_arr(:, :, atomA)

         dZ_dRA = s_1*de1_dRA
         d_gradZ_dRA = s_1*dgrad_e1_dRA
         if (n_active >= 2) then
            dZ_dRA = dZ_dRA + s_2_p*(sum_e2*de2_dRA - 0.5_wp*de1_dRA)
            d_gradZ_dRA = d_gradZ_dRA + s_2_p*( &
                          outer_matrix(dsum_e2, de2_dRA) + sum_e2*dgrad_e2_dRA - 0.5_wp*dgrad_e1_dRA)
         end if
         if (n_active >= 3) then
            dZ_dRA = dZ_dRA + s_3*( &
                     0.5_wp*sum_e3*sum_e3*de3_dRA &
                     - 0.5_wp*(de3_dRA*sum_e3_sq + sum_e3*de3_sq_dRA) &
                     + (1.0_wp/3.0_wp)*de1_dRA)
            d_gradZ_dRA = d_gradZ_dRA + s_3*( &
                          sum_e3*outer_matrix(dsum_e3, de3_dRA) &
                          + 0.5_wp*sum_e3*sum_e3*dgrad_e3_dRA &
                          - 0.5_wp*( &
                          outer_matrix(dsum_e3, de3_sq_dRA) + sum_e3_sq*dgrad_e3_dRA &
                          + outer_matrix(dsum_e3_sq, de3_dRA) + sum_e3*dgrad_e3_sq_dRA) &
                          + (1.0_wp/3.0_wp)*dgrad_e1_dRA)
         end if

         if (atomA == owner_id) then
            de1_owner_dRA = de1_dRA
            de2_owner_dRA = de2_dRA
            de3_owner_dRA = de3_dRA
            de3_sq_owner_dRA = de3_sq_dRA
            dgrad_e1_owner_dRA = dgrad_e1_dRA
            dgrad_e2_owner_dRA = dgrad_e2_dRA
            dgrad_e3_owner_dRA = dgrad_e3_dRA
            dgrad_e3_sq_owner_dRA = dgrad_e3_sq_dRA
         else
            de1_owner_dRA = 0.0_wp
            de2_owner_dRA = 0.0_wp
            de3_owner_dRA = 0.0_wp
            de3_sq_owner_dRA = 0.0_wp
            dgrad_e1_owner_dRA = 0.0_wp
            dgrad_e2_owner_dRA = 0.0_wp
            dgrad_e3_owner_dRA = 0.0_wp
            dgrad_e3_sq_owner_dRA = 0.0_wp
         end if

         dM_i_dRA = s_1*de1_owner_dRA
         d_gradM_i_dRA = s_1*dgrad_e1_owner_dRA

         if (n_active >= 2) then
            dM_i_dRA = dM_i_dRA + 0.5_wp*s_2_p*( &
                       de2_owner_dRA*(sum_e2 - 2.0_wp*e2_val(owner_id)) + e2_val(owner_id)*de2_dRA)
            d_gradM_i_dRA = d_gradM_i_dRA + 0.5_wp*s_2_p*( &
                            dgrad_e2_owner_dRA*(sum_e2 - 2.0_wp*e2_val(owner_id)) &
                            + outer_matrix(grad_e2_arr(:, owner_id), de2_dRA - 2.0_wp*de2_owner_dRA) &
                            + outer_matrix(dsum_e2, de2_owner_dRA) &
                            + e2_val(owner_id)*dgrad_e2_dRA)
         end if

         if (n_active >= 3) then
            dsum_e3_excl_dRA = de3_dRA - de3_owner_dRA
            dsum_e3_sq_excl_dRA = de3_sq_dRA - de3_sq_owner_dRA
            dpair_sum_excl_dRA = sum_e3_excl*dsum_e3_excl_dRA - 0.5_wp*dsum_e3_sq_excl_dRA
            dM_i_dRA = dM_i_dRA + (s_3/3.0_wp)*(de3_owner_dRA*pair_sum_excl + e3_val(owner_id)*dpair_sum_excl_dRA)
            d_dsum_e3_excl_dRA = dgrad_e3_dRA - dgrad_e3_owner_dRA
            d_dsum_e3_sq_excl_dRA = dgrad_e3_sq_dRA - dgrad_e3_sq_owner_dRA
            d_dpair_sum_excl_dRA = outer_matrix(dsum_e3_excl, dsum_e3_excl_dRA) &
                                   + sum_e3_excl*d_dsum_e3_excl_dRA - 0.5_wp*d_dsum_e3_sq_excl_dRA
            d_gradM_i_dRA = d_gradM_i_dRA + (s_3/3.0_wp)*( &
                            dgrad_e3_owner_dRA*pair_sum_excl &
                            + outer_matrix(grad_e3_arr(:, owner_id), dpair_sum_excl_dRA) &
                            + outer_matrix(dpair_sum_excl, de3_owner_dRA) &
                            + e3_val(owner_id)*d_dpair_sum_excl_dRA)
         end if

         d2weight_r_rA(:, :, atomA) = d_gradM_i_dRA*invZ &
                                      - outer_matrix(grad_M_i, dZ_dRA)*invZ2 &
                                      - (outer_matrix(gradZ, dM_i_dRA) + M_i*d_gradZ_dRA)*invZ2 &
                                      + 2.0_wp*M_i*outer_matrix(gradZ, dZ_dRA)*invZ3
      end do
   end subroutine lsf_pou_f2_r_rA_screened

!> Compute normalized POU first nuclear derivative (screened)
   subroutine lsf_normalized_f01_rA_screened(self, val, deriv_rA)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      real(wp), intent(out) :: val
      real(wp), allocatable, intent(out), optional :: deriv_rA(:, :)

      real(wp) :: lsf_val, grad_norm, inv_grad_norm, inv_grad_norm3
      real(wp) :: lsf_grad(ndim), gd(ndim)
      real(wp), allocatable :: lsf_deriv_rA(:, :), grad_hess_rA(:, :, :)
      real(wp) :: k, s_1, s_2, s_3, Z
      logical :: need_deriv
      integer :: atom, i, j
      real(wp), parameter :: eps = 1.0e-14_wp

      need_deriv = present(deriv_rA)
      call self%f012_r_screened(lsf_val, lsf_grad)
      grad_norm = norm2(lsf_grad)

      if (grad_norm < eps) then
         val = 0.0_wp
         if (need_deriv) allocate (deriv_rA(ndim, self%ncenters), source=0.0_wp)
         return
      end if

      inv_grad_norm = 1.0_wp/grad_norm
      val = lsf_val*inv_grad_norm

      if (need_deriv) then
         allocate (deriv_rA(ndim, self%ncenters))
         k = self%param%blend_k
         s_1 = self%param%blend_1b
         s_2 = self%param%blend_2b
         s_3 = self%param%blend_3b

         call compute_lsf_z0_screened(self%ssd_system, k, s_1, s_2, s_3, Z)
         call compute_lsf_f1_rA_screened(self%ssd_system, k, s_1, s_2, s_3, Z, lsf_deriv_rA)
         call compute_lsf_f2_r_rA_screened(self%ssd_system, k, s_1, s_2, s_3, Z, grad_hess_rA)

         inv_grad_norm3 = inv_grad_norm*inv_grad_norm*inv_grad_norm
         deriv_rA = 0.0_wp
         do i = 1, self%ssd_system%n_active
            atom = self%ssd_system%atom_indices(i)
            gd = 0.0_wp
            do j = 1, ndim
               gd = gd + lsf_grad(j)*grad_hess_rA(j, :, i)
            end do
            deriv_rA(:, atom) = lsf_deriv_rA(:, i)*inv_grad_norm &
                                - lsf_val*gd*inv_grad_norm3
         end do
      end if
   end subroutine lsf_normalized_f01_rA_screened

!* ================================================================================= *!
!*         Private helpers: Faa di Bruno, tensor Leibniz, S^3 accumulators           *!
!* ================================================================================= *!
!
! Strategy: compute global per-kind r-derivative tensors of the four power-sum
! quantities (S1, S2, S3, T3 = sum e3^2), where each sum is built from
! e_{k,i} = exp(-lambda_k * d_i) over active atoms i.  All four lambda values
! are stored locally.  Each kind k yields tensors S_k, dS_k, d2S_k, d3S_k, d4S_k
! by Faa di Bruno on the per-atom exponentials. Z and its r-derivatives are
! assembled by tensorial Leibniz from the polynomial
!     Z = C1*S1 + (s_2/2)*S2^2 + (s_3/6)*S3^3 - (s_3/2)*S3*T3
! where C1 = s_1 - (s_2/2)*[n>=2] + (s_3/3)*[n>=3].
!
! R-derivative helpers reuse the same per-atom Faa di Bruno tensors, but instead
! of summing over all active atoms they "restrict to atom A only" and apply one
! additional r-derivative (sign-flipped). This expresses d/ dR_{A,s_1} as
! -P_A applied to d/ dr_s_1, where P_A restricts S_k accumulators to atom A.

!> Compute exp-derivative tensors for a single atom with given lambda.
!> Returns the n-th r-derivative tensors of e^{-lambda * d} via Faa di Bruno
!> applied to chain f1, f2, f3, f4 of d.
!>
!> Outputs:
!>   v0 = e (scalar value)
!>   v1(:) = de/ dr = -lambda * e * f1
!>   v2(:,:) = d^2 e/ dr^2 = e * (lambda^2 f1 f1 - lambda f2)
!>   v3(:,:,:) = d^3 e/ dr^3 = e * (-lambda^3 f1^3 + lambda^2 [sym f1 f2] - lambda f3)
!>   v4(:,:,:,:) = d^4 e/ dr^4 (full Faa di Bruno fourth)
   pure subroutine atom_exp_derivs(e, lam, f1, f2, f3, f4, v0, v1, v2, v3, v4)
      real(wp), intent(in) :: e, lam
      real(wp), intent(in) :: f1(ndim), f2(ndim, ndim)
      real(wp), intent(in) :: f3(ndim, ndim, ndim), f4(ndim, ndim, ndim, ndim)
      real(wp), intent(out) :: v0
      real(wp), intent(out) :: v1(ndim)
      real(wp), intent(out) :: v2(ndim, ndim)
      real(wp), intent(out) :: v3(ndim, ndim, ndim)
      real(wp), intent(out) :: v4(ndim, ndim, ndim, ndim)
      real(wp) :: lam2, lam3, lam4
      integer :: i, j, k_, l

      lam2 = lam*lam
      lam3 = lam2*lam
      lam4 = lam2*lam2

      v0 = e
      v1 = -lam*e*f1
      v2 = e*(lam2*outer_matrix(f1, f1) - lam*f2)

      ! v3_{ijk} = e * (-lam^3 f1_i f1_j f1_k
      !                 + lam^2 (f2_{ij} f1_k + f2_{ik} f1_j + f2_{jk} f1_i)
      !                 - lam f3_{ijk})
      do k_ = 1, ndim
         do j = 1, ndim
            do i = 1, ndim
               v3(i, j, k_) = e*( &
                              -lam3*f1(i)*f1(j)*f1(k_) &
                              + lam2*(f2(i, j)*f1(k_) + f2(i, k_)*f1(j) + f2(j, k_)*f1(i)) &
                              - lam*f3(i, j, k_))
            end do
         end do
      end do

      ! v4_{ijkl} (Faa di Bruno fourth derivative of e^{-lam d}):
      !   = e * [ lam^4 f1_i f1_j f1_k f1_l
      !           - lam^3 sym6{f2_{ij} f1_k f1_l}     (6 unordered pairs of 4)
      !           + lam^2 ( sym4{f3_{ijk} f1_l}        (4 splits of 3+1)
      !                   + sym3{f2_{ij} f2_{kl}} )    (3 unordered pair partitions)
      !           - lam f4_{ijkl} ]
      do l = 1, ndim
         do k_ = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  v4(i, j, k_, l) = e*( &
                                    lam4*f1(i)*f1(j)*f1(k_)*f1(l) &
                                    - lam3*( &
                                    f2(i, j)*f1(k_)*f1(l) &
                                    + f2(i, k_)*f1(j)*f1(l) &
                                    + f2(i, l)*f1(j)*f1(k_) &
                                    + f2(j, k_)*f1(i)*f1(l) &
                                    + f2(j, l)*f1(i)*f1(k_) &
                                    + f2(k_, l)*f1(i)*f1(j)) &
                                    + lam2*( &
                                    f3(i, j, k_)*f1(l) &
                                    + f3(i, j, l)*f1(k_) &
                                    + f3(i, k_, l)*f1(j) &
                                    + f3(j, k_, l)*f1(i) &
                                    + f2(i, j)*f2(k_, l) &
                                    + f2(i, k_)*f2(j, l) &
                                    + f2(i, l)*f2(j, k_)) &
                                    - lam*f4(i, j, k_, l))
               end do
            end do
         end do
      end do
   end subroutine atom_exp_derivs

!> Tensor Leibniz: d^n_{ab...}(u v) given all derivatives up to order n of u, v.
!> Here implemented up to n=4. The result tensors are accumulators (intent inout
!> += operation), suitable for assembling Z r-derivatives.

!> Add C * u * v to a scalar accumulator.
   pure subroutine accum_uv_0(acc, c, u, v)
      real(wp), intent(inout) :: acc
      real(wp), intent(in) :: c, u, v
      acc = acc + c*u*v
   end subroutine accum_uv_0

!> Add C * d(uv) = C * (u_a v + u v_a) to a rank-1 accumulator.
   pure subroutine accum_uv_1(acc, c, u, u1, v, v1)
      real(wp), intent(inout) :: acc(ndim)
      real(wp), intent(in) :: c, u, v
      real(wp), intent(in) :: u1(ndim), v1(ndim)
      acc = acc + c*(u1*v + u*v1)
   end subroutine accum_uv_1

!> Add C * d^2 (uv) = C * (u_ab v + u_a v_b + u_b v_a + u v_ab) to rank-2 acc.
   pure subroutine accum_uv_2(acc, c, u, u1, u2, v, v1, v2)
      real(wp), intent(inout) :: acc(ndim, ndim)
      real(wp), intent(in) :: c, u, v
      real(wp), intent(in) :: u1(ndim), v1(ndim)
      real(wp), intent(in) :: u2(ndim, ndim), v2(ndim, ndim)
      acc = acc + c*(u2*v + outer_matrix(u1, v1) + outer_matrix(v1, u1) + u*v2)
   end subroutine accum_uv_2

!> Add C * d^3 (uv) = C * (u_abc v + 3sym{u_ab v_c} + 3sym{u_a v_bc} + u v_abc)
   pure subroutine accum_uv_3(acc, c, u, u1, u2, u3, v, v1, v2, v3)
      real(wp), intent(inout) :: acc(ndim, ndim, ndim)
      real(wp), intent(in) :: c, u, v
      real(wp), intent(in) :: u1(ndim), v1(ndim)
      real(wp), intent(in) :: u2(ndim, ndim), v2(ndim, ndim)
      real(wp), intent(in) :: u3(ndim, ndim, ndim), v3(ndim, ndim, ndim)
      acc = acc + c*( &
            u3*v &
            + sym3_21(u2, v1) &
            + sym3_21(v2, u1) &
            + u*v3)
   end subroutine accum_uv_3

!> Add C * d^4 (uv) = C * (u_abcd v + 4sym{u_abc v_d} + 6sym{u_ab v_cd}
!>                        + 4sym{u_a v_bcd} + u v_abcd)
   pure subroutine accum_uv_4(acc, c, u, u1, u2, u3, u4, v, v1, v2, v3, v4)
      real(wp), intent(inout) :: acc(ndim, ndim, ndim, ndim)
      real(wp), intent(in) :: c, u, v
      real(wp), intent(in) :: u1(ndim), v1(ndim)
      real(wp), intent(in) :: u2(ndim, ndim), v2(ndim, ndim)
      real(wp), intent(in) :: u3(ndim, ndim, ndim), v3(ndim, ndim, ndim)
      real(wp), intent(in) :: u4(ndim, ndim, ndim, ndim), v4(ndim, ndim, ndim, ndim)
      ! Tensor Leibniz: 1+4+6+4+1 = 16 ordered ways to distribute 4 derivatives.
      ! The 2-2 partition contributes 6 ordered terms; sym4_22 only covers
      ! 3 of them (those with index i on the "u" side), so we add it twice with
      ! swapped arguments to cover all 6 ordered partitions.
      acc = acc + c*( &
            u4*v &
            + sym4_31(v1, u3) &
            + sym4_22(u2, v2) &
            + sym4_22(v2, u2) &
            + sym4_31(u1, v3) &
            + u*v4)
   end subroutine accum_uv_4

!> Accumulate (coef * S^3) and all r-derivatives via two-step Leibniz:
!>   P := S * S; then Z += coef * P * S.
!> We compute derivatives of P up to order 4 first, then apply accum_uv_n
!> with (P, S).
   pure subroutine accum_S3_cubed(d4_acc, d3_acc, d2_acc, d1_acc, d0_acc, &
                                  coef, s0, s1, s2, s3, s4)
      real(wp), intent(inout) :: d4_acc(ndim, ndim, ndim, ndim)
      real(wp), intent(inout) :: d3_acc(ndim, ndim, ndim)
      real(wp), intent(inout) :: d2_acc(ndim, ndim)
      real(wp), intent(inout) :: d1_acc(ndim)
      real(wp), intent(inout) :: d0_acc
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: s0
      real(wp), intent(in) :: s1(ndim)
      real(wp), intent(in) :: s2(ndim, ndim)
      real(wp), intent(in) :: s3(ndim, ndim, ndim)
      real(wp), intent(in) :: s4(ndim, ndim, ndim, ndim)

      real(wp) :: p0, p1(ndim), p2(ndim, ndim)
      real(wp) :: p3(ndim, ndim, ndim), p4(ndim, ndim, ndim, ndim)

      ! P = s * s
      p0 = 0.0_wp; p1 = 0.0_wp; p2 = 0.0_wp; p3 = 0.0_wp; p4 = 0.0_wp
      call accum_uv_0(p0, 1.0_wp, s0, s0)
      call accum_uv_1(p1, 1.0_wp, s0, s1, s0, s1)
      call accum_uv_2(p2, 1.0_wp, s0, s1, s2, s0, s1, s2)
      call accum_uv_3(p3, 1.0_wp, s0, s1, s2, s3, s0, s1, s2, s3)
      call accum_uv_4(p4, 1.0_wp, s0, s1, s2, s3, s4, s0, s1, s2, s3, s4)

      ! Z += coef * P * S
      call accum_uv_0(d0_acc, coef, p0, s0)
      call accum_uv_1(d1_acc, coef, p0, p1, s0, s1)
      call accum_uv_2(d2_acc, coef, p0, p1, p2, s0, s1, s2)
      call accum_uv_3(d3_acc, coef, p0, p1, p2, p3, s0, s1, s2, s3)
      call accum_uv_4(d4_acc, coef, p0, p1, p2, p3, p4, s0, s1, s2, s3, s4)
   end subroutine accum_S3_cubed

!> Accumulate R-derivative of (coef * S^3) and r-derivatives up to order 3 via:
!>   d/dR(S^3) = 3 S^2 (dS/dR), so we add coef * d/dR[S*S*S] which by Leibniz
!>   is coef * (RS * S * S + S * RS * S + S * S * RS) = 3 * coef * S^2 * RS.
!> For r-derivatives we expand the product P3 = S * S * S and take d/dR via
!> single-factor replacement.
   pure subroutine accum_R_S3_cubed(d3_acc, d2_acc, d1_acc, d0_acc, coef, &
                                    s0, s1, s2, s3, rs0, rs1, rs2, rs3)
      real(wp), intent(inout) :: d3_acc(ndim, ndim, ndim)
      real(wp), intent(inout) :: d2_acc(ndim, ndim)
      real(wp), intent(inout) :: d1_acc(ndim)
      real(wp), intent(inout) :: d0_acc
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: s0, rs0
      real(wp), intent(in) :: s1(ndim), rs1(ndim)
      real(wp), intent(in) :: s2(ndim, ndim), rs2(ndim, ndim)
      real(wp), intent(in) :: s3(ndim, ndim, ndim), rs3(ndim, ndim, ndim)

      real(wp) :: p0, p1(ndim), p2(ndim, ndim), p3(ndim, ndim, ndim)
      ! P = S * S; build derivatives of P
      p0 = 0.0_wp; p1 = 0.0_wp; p2 = 0.0_wp; p3 = 0.0_wp
      call accum_uv_0(p0, 1.0_wp, s0, s0)
      call accum_uv_1(p1, 1.0_wp, s0, s1, s0, s1)
      call accum_uv_2(p2, 1.0_wp, s0, s1, s2, s0, s1, s2)
      call accum_uv_3(p3, 1.0_wp, s0, s1, s2, s3, s0, s1, s2, s3)

      ! d/dR (S*S*S) replaces one factor with RS. Since all three factors are the
      ! same function S, the three results are the identical product RS*S*S:
      !   (RS)*S*S + S*(RS)*S + S*S*(RS) = 3 * RS * S^2.
      ! Their full r-derivative tensors are identical too (the n-th derivative of a
      ! product does not depend on how the factors are grouped), so the factor of 3
      ! is pure multiplicity, not three distinct values. We evaluate the same
      ! quantity via two equivalent groupings and sum to 3 * (RS * S * S):
      !   term1:       (RS) * P     with P = S*S
      !   term2+term3: S * Q  (x2)  with Q = RS*S

      ! Term 1: RS * S * S = RS * P
      call accum_uv_0(d0_acc, coef, rs0, p0)
      call accum_uv_1(d1_acc, coef, rs0, rs1, p0, p1)
      call accum_uv_2(d2_acc, coef, rs0, rs1, rs2, p0, p1, p2)
      call accum_uv_3(d3_acc, coef, rs0, rs1, rs2, rs3, p0, p1, p2, p3)

      ! Term 2: S * RS * S = S * Q where Q = RS * S
      call accum_R_S3_term_with_Q_RS_S(d3_acc, d2_acc, d1_acc, d0_acc, coef, &
                                       s0, s1, s2, s3, rs0, rs1, rs2, rs3)
      ! Term 3: S * S * RS = S * Q where Q = S * RS (identical to Q above).  Add same.
      call accum_R_S3_term_with_Q_RS_S(d3_acc, d2_acc, d1_acc, d0_acc, coef, &
                                       s0, s1, s2, s3, rs0, rs1, rs2, rs3)
   end subroutine accum_R_S3_cubed

!> Helper: add coef * S * (RS * S) and its r-derivatives up to order 3.
   pure subroutine accum_R_S3_term_with_Q_RS_S(d3_acc, d2_acc, d1_acc, d0_acc, coef, &
                                               s0, s1, s2, s3, rs0, rs1, rs2, rs3)
      real(wp), intent(inout) :: d3_acc(ndim, ndim, ndim)
      real(wp), intent(inout) :: d2_acc(ndim, ndim)
      real(wp), intent(inout) :: d1_acc(ndim)
      real(wp), intent(inout) :: d0_acc
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: s0, rs0
      real(wp), intent(in) :: s1(ndim), rs1(ndim)
      real(wp), intent(in) :: s2(ndim, ndim), rs2(ndim, ndim)
      real(wp), intent(in) :: s3(ndim, ndim, ndim), rs3(ndim, ndim, ndim)

      real(wp) :: q0, q1(ndim), q2(ndim, ndim), q3(ndim, ndim, ndim)
      q0 = 0.0_wp; q1 = 0.0_wp; q2 = 0.0_wp; q3 = 0.0_wp
      ! Q = RS * S
      call accum_uv_0(q0, 1.0_wp, rs0, s0)
      call accum_uv_1(q1, 1.0_wp, rs0, rs1, s0, s1)
      call accum_uv_2(q2, 1.0_wp, rs0, rs1, rs2, s0, s1, s2)
      call accum_uv_3(q3, 1.0_wp, rs0, rs1, rs2, rs3, s0, s1, s2, s3)
      ! Add coef * S * Q
      call accum_uv_0(d0_acc, coef, s0, q0)
      call accum_uv_1(d1_acc, coef, s0, s1, q0, q1)
      call accum_uv_2(d2_acc, coef, s0, s1, s2, q0, q1, q2)
      call accum_uv_3(d3_acc, coef, s0, s1, s2, s3, q0, q1, q2, q3)
   end subroutine accum_R_S3_term_with_Q_RS_S

!> Add d^2 (coef * S^3)/ dR_A dR_B and r-derivatives up to order 2 via:
!>    d^2 (S^3) = 3 S^2 RRS + 6 S (RS_A)(RS_B)
!> applied at each r-derivative order with Leibniz on all factors.
   pure subroutine accum_RR_S3_cubed(d2_acc, d1_acc, d0_acc, coef, &
                                     s0, s1, s2, &
                                     rsa0, rsa1, rsa2, &
                                     rsb0, rsb1, rsb2, &
                                     rrs0, rrs1, rrs2)
      real(wp), intent(inout) :: d2_acc(ndim, ndim)
      real(wp), intent(inout) :: d1_acc(ndim)
      real(wp), intent(inout) :: d0_acc
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: s0
      real(wp), intent(in) :: s1(ndim)
      real(wp), intent(in) :: s2(ndim, ndim)
      real(wp), intent(in) :: rsa0, rsb0, rrs0
      real(wp), intent(in) :: rsa1(ndim), rsb1(ndim), rrs1(ndim)
      real(wp), intent(in) :: rsa2(ndim, ndim), rsb2(ndim, ndim), rrs2(ndim, ndim)

      real(wp) :: p0, p1(ndim), p2(ndim, ndim)
      real(wp) :: ssa0, ssa1(ndim), ssa2(ndim, ndim)  ! RSa * RSb product and derivatives
      !  d^2 (S^3) via expanding (a)(b)(c) with two R-derivatives:
      ! 1) both on same factor: 3 ways x (RRS * S * S) -> total 3 (RRS)(S^2)
      ! 2) one on each of two different factors: 3 pairs x 2 orderings = 6 (S)(RSa)(RSb)
      !
      ! Group by structure:
      !   3 * (RRS * S^2) + 6 * (S * RSa * RSb)
      ! Each of these is a triple product whose r-derivatives we evaluate via two
      ! Leibniz steps.

      ! P := S * S
      p0 = 0.0_wp; p1 = 0.0_wp; p2 = 0.0_wp
      call accum_uv_0(p0, 1.0_wp, s0, s0)
      call accum_uv_1(p1, 1.0_wp, s0, s1, s0, s1)
      call accum_uv_2(p2, 1.0_wp, s0, s1, s2, s0, s1, s2)

      ! 3 * coef * (RRS * P)
      call accum_uv_0(d0_acc, 3.0_wp*coef, rrs0, p0)
      call accum_uv_1(d1_acc, 3.0_wp*coef, rrs0, rrs1, p0, p1)
      call accum_uv_2(d2_acc, 3.0_wp*coef, rrs0, rrs1, rrs2, p0, p1, p2)

      ! 6 * coef * (S * RSa * RSb) via Q = RSa * RSb, then S * Q
      ssa0 = 0.0_wp; ssa1 = 0.0_wp; ssa2 = 0.0_wp
      call accum_uv_0(ssa0, 1.0_wp, rsa0, rsb0)
      call accum_uv_1(ssa1, 1.0_wp, rsa0, rsa1, rsb0, rsb1)
      call accum_uv_2(ssa2, 1.0_wp, rsa0, rsa1, rsa2, rsb0, rsb1, rsb2)
      call accum_uv_0(d0_acc, 6.0_wp*coef, s0, ssa0)
      call accum_uv_1(d1_acc, 6.0_wp*coef, s0, s1, ssa0, ssa1)
      call accum_uv_2(d2_acc, 6.0_wp*coef, s0, s1, s2, ssa0, ssa1, ssa2)
   end subroutine accum_RR_S3_cubed

!> Finalizer for LSF primitive LSF: free owned allocatable components.
   subroutine finalize_lsf_svdw(self)
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self

      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%all_indices)) deallocate (self%all_indices)

      ! Deallocate structure_type allocatable components
      if (allocated(self%mol%id)) deallocate (self%mol%id)
      if (allocated(self%mol%num)) deallocate (self%mol%num)
      if (allocated(self%mol%sym)) deallocate (self%mol%sym)
      if (allocated(self%mol%xyz)) deallocate (self%mol%xyz)
      if (allocated(self%mol%lattice)) deallocate (self%mol%lattice)
      if (allocated(self%mol%periodic)) deallocate (self%mol%periodic)
      if (allocated(self%mol%bond)) deallocate (self%mol%bond)
      if (allocated(self%mol%comment)) deallocate (self%mol%comment)
      if (allocated(self%mol%sdf)) deallocate (self%mol%sdf)
      if (allocated(self%mol%pdb)) deallocate (self%mol%pdb)

   end subroutine finalize_lsf_svdw

!> Conservative shell thickness from the SvdW screening threshold.
!>
!> Finds the radial offset `delta >= 0` for which
!> `exp(-k/3 * ssd0(r=radius+delta)) = 0.1 * threshold`. The returned
!> distance is measured from the sphere surface. Used by the DROP cavity
!> via the abstract LSF `neighbor_cutoff` to size per-atom cell-grid
!> reaches without knowing the LSF concrete.
!>
!> @param[in] self    SvdW LSF instance (reads param%blend_k, screening_threshold)
!> @param[in] radius  Atom radius (Bohr)
!> @returns           Radial offset from atom surface (Bohr)
   pure function lsf_neighbor_cutoff(self, radius) result(cutoff_distance)
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from atom surface (Bohr)
      real(wp) :: cutoff_distance

      real(wp) :: center(3), point(3)
      real(wp) :: lower, upper, middle
      real(wp) :: screening_target
      real(wp) :: weight_upper
      real(wp) :: k_local, threshold
      integer :: iter

      k_local = self%param%blend_k
      threshold = self%screening_threshold

      center = 0.0_wp
      screening_target = max(0.1_wp*threshold, tiny(1.0_wp))

      lower = 0.0_wp
      upper = 1.0_wp
      do
         point = [radius + upper, 0.0_wp, 0.0_wp]
         weight_upper = exp(-(k_local/3.0_wp)*ssd0(point, center, radius))
         if (weight_upper <= screening_target) exit
         upper = 2.0_wp*upper
         if (upper > 1.0e4_wp) exit
      end do

      do iter = 1, 64
         middle = 0.5_wp*(lower + upper)
         point = [radius + middle, 0.0_wp, 0.0_wp]
         if (exp(-(k_local/3.0_wp)*ssd0(point, center, radius)) > screening_target) then
            lower = middle
         else
            upper = middle
         end if
      end do

      cutoff_distance = upper
   end function lsf_neighbor_cutoff

end module moist_cavity_drop_lsf_svdw
