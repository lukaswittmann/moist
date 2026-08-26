!> Smooth van der Waals (SvdW) level set function
!>
!> The mathematics lives in the code-generated module
!> [[moist_cavity_drop_lsf_svdw_kernel]]; this module is the orchestration layer
!> between it and the LSF contract of [[moist_cavity_drop_lsf_base]].
!>
!> Division of labour
!> ------------------
!> The kernel is symbolic in the *power sums* of the per-atom screening factors
!>
!>    u_A = exp(-(k/3)(||r - R_A|| - R_a)) ,   P_m = sum_A u_A**m ,
!>
!> for m = 1, 2, 3 and 3/2, from which `Z` (and hence `S = -(1/k) log Z`) and
!> every derivative follow in O(n_active) work. So:
!>
!>   - `prepare` / `prepare_subset` run the base screening gate and cache the
!>     *minimum*: per active atom its user-space id, its displacement `r - R_A`,
!>     its radius and the distance `||r - R_A||`, plus the accumulated power-sum
!>     tensors. Every one of those buffers is sized once per `update` and never
!>     reallocated per point.
!>   - every accessor rebuilds whatever tensors it needs on demand by calling the
!>     kernel. Nothing derivative-shaped is stored between points, which is why
!>     this module no longer owns an SSD derivative cache.
!>
!> Buffers, never allocations
!> --------------------------
!> No accessor returns an `allocatable, intent(out)` result. Every result is a
!> caller-provided buffer whose nuclear extent the caller sizes from
!> [[active_count]]; the accessor writes its first `active_count()` slots and
!> leaves the rest alone. This is what keeps the OpenMP projection path free of
!> per-call heap traffic.
!>
!> Index space of the nuclear outputs
!> ----------------------------------
!> Every nuclear index is an *active-list* index: slot `i` belongs to the atom
!> `active_atom(i)`. Screened-away atoms have no slot at all, which is also why
!> the outputs no longer have to be zeroed over the whole molecule on entry.
!>
!> Derivative-order contract
!> -------------------------
!> `prepare` is the only thing that can be under-provisioned, and the only thing
!> it provisions is the power sums, up to `max_deriv`. Each accessor therefore
!> asks [[require_deriv]] for exactly the power-sum order the kernel routine it
!> calls reads -- not for the order of the tensor it returns, which it rebuilds
!> itself. `f4_rrr_rA`, for instance, needs `ps3` and builds the order-4 per-atom
!> tensors on the spot, so it runs correctly at `set_max_deriv(3)`.
!>
!> Two things the kernel deliberately does not handle, and are handled here
!> ------------------------------------------------------------------------
!>   1. **A point exactly on a nucleus.** `svdw_atom_eval` divides by
!>      `||r - R_A||`. [[atom_tensors]] and [[atom_tangent_tensors]] intercept
!>      that: at `x == 0` the value tensor is still the kernel's own `u_A`, and
!>      every derivative order is zero -- the same convention the retired SSD
!>      fill used (`n = (r - R_A)/x` taken as zero).
!>   2. **The low-`n_active` guards.** The historical Z assembly switched the
!>      two- and three-body blocks off below two and three active atoms. Those
!>      elementary symmetric polynomials vanish *identically* there, so the
!>      guards never changed a value in exact arithmetic -- but they did suppress
!>      the floating-point residue of `q**2 - p3` and `p1**3 - 3 p1 p2 + 2 p3`
!>      cancelling to zero. [[svdw_weights]] reproduces them exactly and for free
!>      by handing the kernel `s_2 = 0` below two atoms and `s_3 = 0` below three:
!>      multiplying the vanishing block by an exact zero is the same arithmetic as
!>      skipping it.
module moist_cavity_drop_lsf_svdw
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_sorted
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_svdw_kernel, only: nkind, &
                                                svdw_atom_eval, svdw_atom_tangent_eval, &
                                                svdw_spatial_eval, svdw_nuclear_eval, &
                                                svdw_radius_eval, svdw_radius_hvp_eval, &
                                                svdw_atom_radius_tangent_eval, &
                                                svdw_pair_eval, svdw_pair_diag_eval, &
                                                svdw_radpair_eval, svdw_radpair_diag_eval, &
                                                svdw_nucrad_eval, svdw_nucrad_diag_eval, &
                                                svdw_tangent_eval, svdw_hvp_eval, &
                                                svdw_vjp_eval, svdw_radius_vjp_eval, &
                                                svdw_normalized_eval, svdw_powersums
   implicit none (type, external)
   private

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Smooth van der Waals LSF
   !>
   !> Concrete level set function: takes its blending parameters directly through
   !> [[new]], caches the screened per-atom geometry and the power sums via
   !> [[prepare]], and rebuilds every derivative tensor on demand from the
   !> generated kernel. Inherits the common atom-LSF state (ncenters, mol, radii,
   !> screening caches) from [[moist_cavity_drop_lsf_type]].
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_svdw_type

      !> SvdW blending parameters (k, one-, two- and three-body weights)
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      !> Highest spatial-derivative order the power sums are accumulated to
      integer :: max_deriv = 2

      !* --------------------------- Per-point active list -------------------------- *!

      !? These should be in the base lsf type (we need to add these to the isodensity ones)

      !> Number of atoms that survived screening at the cached point
      integer :: n_active = 0
      !> Candidate-space ids of the survivors, first `n_active` valid [ncenters]
      integer, allocatable :: active_cand(:)
      !> User-space atom id of active slot i [ncenters]
      integer, allocatable :: act_atom(:)
      !> Displacement `r - R_A` of active slot i [ndim, ncenters]
      real(wp), allocatable :: act_d(:, :)
      !> Radius of active slot i [ncenters]
      real(wp), allocatable :: act_radius(:)
      !> Distance `||r - R_A||` of active slot i [ncenters]
      !>
      !> The reciprocal distance, the unit vector `(r - R_A)/x` and the screening
      !> factor `u_A` used to be cached next to it. They are deliberately not:
      !> the generated kernel takes `(d, R_a, k)` and recomputes `x` and `u_A`
      !> from them, so caching those three bought nothing and cost an extra
      !> `exp`, a division and three multiplies per atom on the single hottest
      !> loop in the projection. `act_x` stays because the on-nucleus guard
      !> reads it.
      real(wp), allocatable :: act_x(:)

      !* ------------------------------- Power sums ---------------------------------- *!
      !> Orders above `max_deriv` are stale after `prepare`; the accessors guard
      !> on [[require_deriv]] rather than on the contents.

      !? Why save these always?

      !> Order-0 power sums; last index selects the kind (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp) :: ps0(nkind) = 0.0_wp
      !> Order-1 power-sum tensors
      real(wp) :: ps1(ndim, nkind) = 0.0_wp
      !> Order-2 power-sum tensors
      real(wp) :: ps2(ndim, ndim, nkind) = 0.0_wp
      !> Order-3 power-sum tensors
      real(wp) :: ps3(ndim, ndim, ndim, nkind) = 0.0_wp
      !> Order-4 power-sum tensors
      real(wp) :: ps4(ndim, ndim, ndim, ndim, nkind) = 0.0_wp
   contains
      !> Constructor: configure blending parameters and declare the candidate space
      procedure, public :: new => lsf_new
      !> Bind molecular geometry and resize the per-atom caches
      procedure, public :: update => lsf_update
      !> Point preparation: screening plus power-sum accumulation
      procedure, public :: prepare => lsf_prepare
      !> Point preparation with a caller-provided candidate atom list
      procedure, public :: prepare_subset => lsf_prepare_subset
      !> Configure the highest power-sum order `prepare` accumulates
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      !> Number of atoms active after the latest prepare/prepare_subset
      procedure, public :: active_count => lsf_active_count
      !> User-space atom index of the i-th active atom
      procedure, public :: active_atom => lsf_active_atom
      !> Value only (lowest-cost path)
      procedure, public :: f0 => lsf_f0
      !> Combined value/gradient/Hessian
      procedure, public :: f012_r => lsf_f012_r
      !> Third spatial derivative
      procedure, public :: f3_rrr => lsf_f3_rrr
      !> Mixed third derivative: spatial Hessian w.r.t. nuclear positions
      procedure, public :: f3_rr_rA => lsf_f3_rr_rA
      !> Pure spatial fourth derivative
      procedure, public :: f4_rrrr => lsf_f4_rrrr
      !> Mixed fourth derivative d^4S / dr^3 dR_A
      procedure, public :: f4_rrr_rA => lsf_f4_rrr_rA
      !> Radius derivative d^3S / dr^2 dR_a
      procedure, public :: f3_rr_rad => lsf_f3_rr_rad
      !> Pure nuclear Hessian d^2S / dR_A dR_B
      procedure, public :: f2_rArB => lsf_f2_rArB
      !> Mixed third derivative d^3S / dr dR_A dR_B
      procedure, public :: f3_r_rArB => lsf_f3_r_rArB
      !> Mixed fourth derivative d^4S / dr^2 dR_A dR_B
      procedure, public :: f4_rr_rArB => lsf_f4_rr_rArB
      !> Pure radius Hessian d^2S / dR_a dR_b
      procedure, public :: f2_radrad => lsf_f2_radrad
      !> Mixed third derivative d^3S / dr dR_a dR_b
      procedure, public :: f3_r_radrad => lsf_f3_r_radrad
      !> Mixed fourth derivative d^4S / dr^2 dR_a dR_b
      procedure, public :: f4_rr_radrad => lsf_f4_rr_radrad
      !> Mixed nuclear-radius Hessian d^2S / dR_A dR_b
      procedure, public :: f2_rA_rad => lsf_f2_rA_rad
      !> Mixed third derivative d^3S / dr dR_A dR_b
      procedure, public :: f3_r_rA_rad => lsf_f3_r_rA_rad
      !> Mixed fourth derivative d^4S / dr^2 dR_A dR_b
      procedure, public :: f4_rr_rA_rad => lsf_f4_rr_rA_rad
      !> Normalized level set f0 / ||f1_r|| and its nuclear derivative
      procedure, public :: normalized_f01_rA => lsf_normalized_f01_rA
      !> Directional nuclear derivative of the value
      procedure, public :: tangent_f0 => lsf_tangent_f0
      !> Directional nuclear derivative of the spatial gradient
      procedure, public :: tangent_f1_r => lsf_tangent_f1_r
      !> Directional nuclear derivative of the spatial Hessian
      procedure, public :: tangent_f2_rr => lsf_tangent_f2_rr
      !> Directional nuclear derivative of the third spatial derivative
      procedure, public :: tangent_f3_rrr => lsf_tangent_f3_rrr
      !> Nuclear Hessian-vector product
      procedure, public :: hvp_f1_rA => lsf_hvp_f1_rA
      !> Directional nuclear derivative of `f2_r_rA`
      procedure, public :: hvp_f2_r_rA => lsf_hvp_f2_r_rA
      !> Directional nuclear derivative of `f3_rr_rA`
      procedure, public :: hvp_f3_rr_rA => lsf_hvp_f3_rr_rA
      !> Radius row of the joint Hessian-vector product
      procedure, public :: hvp_f1_rad => lsf_hvp_f1_rad
      !> Joint directional derivative of `f2_r_rad`
      procedure, public :: hvp_f2_r_rad => lsf_hvp_f2_r_rad
      !> Joint directional derivative of `f3_rr_rad`
      procedure, public :: hvp_f3_rr_rad => lsf_hvp_f3_rr_rad
      !> Jet-contracted nuclear gradient (reverse mode)
      procedure, public :: vjp_f1_rA => lsf_vjp_f1_rA
      !> Jet-contracted radius gradient (reverse mode)
      procedure, public :: vjp_f1_rad => lsf_vjp_f1_rad
      !> Exact radial offset where the SvdW weight equals `screening_threshold`
      procedure, public :: screening_offset => lsf_screening_offset
      !> Exact surface-free radius from the 1-Lipschitz property
      procedure, public :: exclusion_radius => lsf_svdw_exclusion_radius
      !> Finalizer
      final :: finalize_lsf_svdw
   end type moist_cavity_drop_lsf_svdw_type

   public :: moist_cavity_drop_lsf_svdw_type

contains

   !* ================================================================================= *!
   !*                              LSF lifecycle methods                                *!
   !* ================================================================================= *!

   !> Configure LSF blending parameters and declare the candidate index space
   !>
   !> @param[inout] self     LSF instance
   !> @param[in]    blend_k  Blending sharpness k (optional)
   !> @param[in]    blend_1b One-body weight (optional)
   !> @param[in]    blend_2b Two-body weight (optional)
   !> @param[in]    blend_3b Three-body weight (optional)
   subroutine lsf_new(self, blend_k, blend_1b, blend_2b, blend_3b)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Blending sharpness k (optional override)
      real(wp), intent(in), optional :: blend_k
      !> One-body weight (optional override)
      real(wp), intent(in), optional :: blend_1b
      !> Two-body weight (optional override)
      real(wp), intent(in), optional :: blend_2b
      !> Three-body weight (optional override)
      real(wp), intent(in), optional :: blend_3b

      ! The screen loop indexes the base's candidate-space geometry mirror
      ! directly, so candidate ids must arrive spatially sorted.
      self%candidate_space = lsf_candidate_space_sorted

      call self%param%new(blend_k=blend_k, blend_1b=blend_1b, &
                          blend_2b=blend_2b, blend_3b=blend_3b)
   end subroutine lsf_new

   !> Bind molecular geometry and resize the per-atom caches
   !>
   !> The base handles the spatial sort, the candidate-space geometry mirror and
   !> the screening bounds; this override only sizes the per-point buffers, once,
   !> to the molecule's atom count.
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    mol    Molecular structure
   !> @param[in]    radii  Per-atom radii (size mol%nat)
   subroutine lsf_update(self, mol, radii)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)

      !> Atom capacity of the per-point buffers
      integer :: n_alloc

      call lsf_base_update(self, mol, radii)

      n_alloc = mol%nat
      if (allocated(self%active_cand)) deallocate (self%active_cand)
      if (allocated(self%act_atom)) deallocate (self%act_atom)
      if (allocated(self%act_d)) deallocate (self%act_d)
      if (allocated(self%act_radius)) deallocate (self%act_radius)
      if (allocated(self%act_x)) deallocate (self%act_x)

      allocate (self%active_cand(n_alloc))
      allocate (self%act_atom(n_alloc))
      allocate (self%act_d(ndim, n_alloc))
      allocate (self%act_radius(n_alloc))
      allocate (self%act_x(n_alloc))

      self%n_active = 0
      self%prepared_deriv = -1
   end subroutine lsf_update

   !> Configure the highest power-sum order `prepare` accumulates
   !>
   !> Nothing is allocated here: the power sums are fixed-size components and the
   !> per-atom buffers are sized by [[update]].
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Requested max derivative order (0..4)
   subroutine lsf_set_max_deriv(self, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Requested max derivative order
      integer, intent(in) :: n

      self%max_deriv = min(4, max(0, n))
   end subroutine lsf_set_max_deriv

   !> Screen every atom at the evaluation point and refresh the caches
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point (3,)
   !> @param[out]   error Evaluation failure (never set by SvdW)
   subroutine lsf_prepare(self, point, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      call lsf_svdw_screen(self, point, self%full_scan_cand)
   end subroutine lsf_prepare

   !> Screen a caller-provided candidate list and refresh the caches
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point (3,)
   !> @param[in]    candidate_indices Atom ids in the base's sorted candidate space
   !> @param[out]   error             Evaluation failure (never set by SvdW)
   subroutine lsf_prepare_subset(self, point, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      call lsf_svdw_screen(self, point, candidate_indices)
   end subroutine lsf_prepare_subset

   !> Run the base screening gate and rebuild the per-point caches
   !>
   !> Pass one is the base's allocation-free reject test; pass two touches only
   !> the survivors and fills the per-atom geometry; pass three accumulates the
   !> power-sum tensors up to `max_deriv`.
   !>
   !> The accumulation has two spellings of the same sum. The fast one is the
   !> kernel's own [[svdw_powersums]], a *fused* accumulator: it never
   !> materialises a per-atom tensor, keeps one scalar per independent tensor
   !> component for the whole loop, and skips the `q` kind when the two-body
   !> weight is zero. That matters because this is the single hottest thing the
   !> projection does -- the unfused spelling below writes and re-reads 52
   !> doubles per atom at `max_deriv = 2` alone, and `svdw_atom_eval` is far too
   !> large for gfortran to inline, so those buffers cannot be optimised away.
   !> The fused form cannot handle a point sitting exactly on a nucleus, so the
   !> fill records whether any survivor is at zero distance and only then falls
   !> back to the guarded loop below. That branch is taken per *point*, not per
   !> atom, and in production it is never taken at all.
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point [ndim]
   !> @param[in]    candidate_indices Candidate ids in the sorted candidate space
   subroutine lsf_svdw_screen(self, point, candidate_indices)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)

      !> Distance from the evaluation point to the active atom
      real(wp) :: x
      !> Loop counter, candidate id, survivor count and derivative order
      integer :: i, c, n, md
      !> Does any survivor sit exactly on the evaluation point?
      logical :: on_nucleus

      self%n_active = 0
      self%prepared_deriv = self%max_deriv
      if (.not. allocated(self%cand_screen)) return

      call self%screen_candidates(point, candidate_indices, self%active_cand, n)
      self%n_active = n

      md = self%max_deriv
      self%ps0 = 0.0_wp
      if (md >= 1) self%ps1 = 0.0_wp
      if (md >= 2) self%ps2 = 0.0_wp
      if (md >= 3) self%ps3 = 0.0_wp
      if (md >= 4) self%ps4 = 0.0_wp
      if (n == 0) return

      on_nucleus = .false.
      do i = 1, n
         c = self%active_cand(i)
         self%act_atom(i) = self%cand_to_user(c)
         self%act_radius(i) = self%cand_radii(c)
         self%act_d(1, i) = point(1) - self%cand_screen(1, c)
         self%act_d(2, i) = point(2) - self%cand_screen(2, c)
         self%act_d(3, i) = point(3) - self%cand_screen(3, c)
         x = sqrt(self%act_d(1, i)*self%act_d(1, i) &
                  + self%act_d(2, i)*self%act_d(2, i) &
                  + self%act_d(3, i)*self%act_d(3, i))
         self%act_x(i) = x
         if (x <= 0.0_wp) on_nucleus = .true.
      end do

      if (on_nucleus) then
         call powersums_guarded(self, md)
      else
         call svdw_powersums(self%act_d(:, 1:n), self%act_radius(1:n), &
                             self%param%blend_k, self%param%blend_2b /= 0.0_wp, md, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4)
      end if
   end subroutine lsf_svdw_screen

   !> Accumulate the power sums with the on-nucleus guard in the loop
   !>
   !> Only reached when some active atom sits exactly at the evaluation point;
   !> see [[atom_tensors]] for the convention applied there.
   !>
   !> Unlike the fused path this always accumulates the `q` kind, even when the
   !> two-body weight is zero and the fused path would leave those slots at zero.
   !> Both are correct: every consumer reaches that kind only through the
   !> two-body weight, so the slots are multiplied by zero either way. Computing
   !> them here keeps the cold path free of one more special case.
   !>
   !> @param[inout] self      LSF instance with the per-atom geometry filled
   !> @param[in]    max_deriv Highest spatial-derivative order to accumulate
   subroutine powersums_guarded(self, max_deriv)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: max_deriv

      !> Per-atom kind tensors of one active atom
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Active-list index
      integer :: ia

      do ia = 1, self%n_active
         call atom_tensors(self, ia, max_deriv, at0, at1, at2, at3, at4)
         self%ps0 = self%ps0 + at0
         if (max_deriv >= 1) self%ps1 = self%ps1 + at1
         if (max_deriv >= 2) self%ps2 = self%ps2 + at2
         if (max_deriv >= 3) self%ps3 = self%ps3 + at3
         if (max_deriv >= 4) self%ps4 = self%ps4 + at4
      end do
   end subroutine powersums_guarded

   !> Number of atoms currently active after the latest prepare/prepare_subset
   !>
   !> @param[in] self LSF instance
   !> @returns   n    Active-atom count
   pure integer function lsf_active_count(self) result(n)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      n = self%n_active
   end function lsf_active_count

   !> User-space atom index of the i-th currently active atom
   !>
   !> @param[in] self LSF instance
   !> @param[in] i    Active-list index (1 <= i <= active_count())
   !> @returns   idx  User-space atom id
   pure integer function lsf_active_atom(self, i) result(idx)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: i
      idx = self%act_atom(i)
   end function lsf_active_atom

   !* ================================================================================= *!
   !*                          Kernel-facing private helpers                            *!
   !* ================================================================================= *!

   !> Blending weights the kernel is called with, low-`n_active` guards folded in
   !>
   !> See the module header: below two (three) active atoms the two-body
   !> (three-body) elementary symmetric polynomial vanishes identically, so
   !> switching the weight to an exact zero is the same value as skipping the
   !> block -- and, unlike evaluating it, carries none of its cancellation noise.
   !>
   !> @param[in]  self LSF instance
   !> @param[out] s_1  One-body weight
   !> @param[out] s_2  Two-body weight, zero below two active atoms
   !> @param[out] s_3  Three-body weight, zero below three active atoms
   pure subroutine svdw_weights(self, s_1, s_2, s_3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> One-body weight
      real(wp), intent(out) :: s_1
      !> Two-body weight
      real(wp), intent(out) :: s_2
      !> Three-body weight
      real(wp), intent(out) :: s_3

      s_1 = self%param%blend_1b
      s_2 = 0.0_wp
      s_3 = 0.0_wp
      if (self%n_active >= 2) s_2 = self%param%blend_2b
      if (self%n_active >= 3) s_3 = self%param%blend_3b
   end subroutine svdw_weights

   !> Per-atom kind tensors of active slot `ia`, guarded against `x == 0`
   !>
   !> The kernel divides by the distance from the first derivative onwards. On a
   !> nucleus that distance is zero and the direction `(r - R_A)/x` is undefined;
   !> the convention here -- the same one the retired SSD fill used -- is to take
   !> it as zero, so the value tensor still carries `u_A = exp(k R_a/3)` while
   !> every derivative order of this atom drops out of the power sums.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..4)
   !> @param[out] at0       Order-0 kind tensor
   !> @param[out] at1       Order-1 kind tensor
   !> @param[out] at2       Order-2 kind tensor
   !> @param[out] at3       Order-3 kind tensor
   !> @param[out] at4       Order-4 kind tensor
   pure subroutine atom_tensors(self, ia, max_deriv, at0, at1, at2, at3, at4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 kind tensor
      real(wp), intent(out) :: at0(nkind)
      !> Order-1 kind tensor
      real(wp), intent(out) :: at1(ndim, nkind)
      !> Order-2 kind tensor
      real(wp), intent(out) :: at2(ndim, ndim, nkind)
      !> Order-3 kind tensor
      real(wp), intent(out) :: at3(ndim, ndim, ndim, nkind)
      !> Order-4 kind tensor
      real(wp), intent(out) :: at4(ndim, ndim, ndim, ndim, nkind)

      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                             max_deriv, at0, at1, at2, at3, at4)
         return
      end if

      ! On the nucleus: order 0 is still the kernel's own u_A, the rest is zero
      call svdw_atom_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                          0, at0, at1, at2, at3, at4)
      if (max_deriv >= 1) at1 = 0.0_wp
      if (max_deriv >= 2) at2 = 0.0_wp
      if (max_deriv >= 3) at3 = 0.0_wp
      if (max_deriv >= 4) at4 = 0.0_wp
   end subroutine atom_tensors

   !> Direction-contracted per-atom kind tensors of active slot `ia`
   !>
   !> Same `x == 0` convention as [[atom_tensors]]: every order of the contracted
   !> tensor carries a nuclear derivative, so all of them vanish there.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  v_a       Nuclear displacement direction of this atom [ndim]
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..3)
   !> @param[out] aw0       Order-0 contracted kind tensor
   !> @param[out] aw1       Order-1 contracted kind tensor
   !> @param[out] aw2       Order-2 contracted kind tensor
   !> @param[out] aw3       Order-3 contracted kind tensor
   pure subroutine atom_tangent_tensors(self, ia, v_a, max_deriv, aw0, aw1, aw2, aw3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Nuclear displacement direction of this atom
      real(wp), intent(in) :: v_a(ndim)
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 contracted kind tensor
      real(wp), intent(out) :: aw0(nkind)
      !> Order-1 contracted kind tensor
      real(wp), intent(out) :: aw1(ndim, nkind)
      !> Order-2 contracted kind tensor
      real(wp), intent(out) :: aw2(ndim, ndim, nkind)
      !> Order-3 contracted kind tensor
      real(wp), intent(out) :: aw3(ndim, ndim, ndim, nkind)

      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%blend_k, v_a, max_deriv, aw0, aw1, aw2, aw3)
         return
      end if

      aw0 = 0.0_wp
      if (max_deriv >= 1) aw1 = 0.0_wp
      if (max_deriv >= 2) aw2 = 0.0_wp
      if (max_deriv >= 3) aw3 = 0.0_wp
   end subroutine atom_tangent_tensors

   !> Radius-contracted per-atom kind tensors of active slot `ia`
   !>
   !> The radius half of the joint `(v, vr)` contraction:
   !> `awr{n} = (blend_k m_kind / 3) vr_a at{n}`. It lands in the same family as
   !> [[atom_tangent_tensors]], so the caller accumulates both into one `ws*` and
   !> every aggregate kernel routine picks the radius term up unchanged.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  vr_a      Radius direction of this atom
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..3)
   !> @param[out] awr0      Order-0 radius-contracted kind tensor
   !> @param[out] awr1      Order-1 radius-contracted kind tensor
   !> @param[out] awr2      Order-2 radius-contracted kind tensor
   !> @param[out] awr3      Order-3 radius-contracted kind tensor
   pure subroutine atom_radius_tangent_tensors(self, ia, vr_a, max_deriv, &
                                               awr0, awr1, awr2, awr3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Radius direction of this atom
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 radius-contracted kind tensor
      real(wp), intent(out) :: awr0(nkind)
      !> Order-1 radius-contracted kind tensor
      real(wp), intent(out) :: awr1(ndim, nkind)
      !> Order-2 radius-contracted kind tensor
      real(wp), intent(out) :: awr2(ndim, ndim, nkind)
      !> Order-3 radius-contracted kind tensor
      real(wp), intent(out) :: awr3(ndim, ndim, ndim, nkind)

      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_radius_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                            self%param%blend_k, vr_a, max_deriv, &
                                            awr0, awr1, awr2, awr3)
         return
      end if

      ! On the nucleus: order 0 is the kernel's own scaled u_A, the rest is zero
      call svdw_atom_radius_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                         self%param%blend_k, vr_a, 0, &
                                         awr0, awr1, awr2, awr3)
      if (max_deriv >= 1) awr1 = 0.0_wp
      if (max_deriv >= 2) awr2 = 0.0_wp
      if (max_deriv >= 3) awr3 = 0.0_wp
   end subroutine atom_radius_tangent_tensors

   !> Accumulate the direction-contracted power sums for a nuclear direction field
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [ndim, ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order to accumulate (0..3)
   !> @param[out] ws0       Order-0 contracted power sums
   !> @param[out] ws1       Order-1 contracted power sums
   !> @param[out] ws2       Order-2 contracted power sums
   !> @param[out] ws3       Order-3 contracted power sums
   !> @param[in]  vrad      Radius directions [ncenters] (optional)
   pure subroutine tangent_powersums(self, v, max_deriv, ws0, ws1, ws2, ws3, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: max_deriv
      !> Order-0 contracted power sums
      real(wp), intent(out) :: ws0(nkind)
      !> Order-1 contracted power sums
      real(wp), intent(out) :: ws1(ndim, nkind)
      !> Order-2 contracted power sums
      real(wp), intent(out) :: ws2(ndim, ndim, nkind)
      !> Order-3 contracted power sums
      real(wp), intent(out) :: ws3(ndim, ndim, ndim, nkind)
      !> Radius half of the joint direction. Absent leaves every operation below
      !> untouched, so a geometry-independent radius model stays bit-for-bit
      real(wp), intent(in), optional :: vrad(:)

      !> Per-atom contracted kind tensors
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      !> Per-atom radius-contracted kind tensors
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Active-list index
      integer :: ia

      ws0 = 0.0_wp
      if (max_deriv >= 1) ws1 = 0.0_wp
      if (max_deriv >= 2) ws2 = 0.0_wp
      if (max_deriv >= 3) ws3 = 0.0_wp

      do ia = 1, self%n_active
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), max_deriv, &
                                   aw0, aw1, aw2, aw3)
         ws0 = ws0 + aw0
         if (max_deriv >= 1) ws1 = ws1 + aw1
         if (max_deriv >= 2) ws2 = ws2 + aw2
         if (max_deriv >= 3) ws3 = ws3 + aw3
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), &
                                             max_deriv, awr0, awr1, awr2, awr3)
            ws0 = ws0 + awr0
            if (max_deriv >= 1) ws1 = ws1 + awr1
            if (max_deriv >= 2) ws2 = ws2 + awr2
            if (max_deriv >= 3) ws3 = ws3 + awr3
         end if
      end do
   end subroutine tangent_powersums

   !> Build the per-atom kind-tensor cache the pair loops consume
   !>
   !> The pair kernels read atom A's tensors once per partner, so the cache turns
   !> an O(n_active**2) exponential count back into O(n_active). It is a local
   !> allocation rather than a persistent buffer because only the three pair
   !> accessors need it, none of them sits on the projection hot path, and a
   !> persistent order-4 cache would cost 3.9 kB per atom per thread.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Highest spatial-derivative order to fill (2..4)
   !> @param[out] ca0       Order-0 kind tensors of every active atom
   !> @param[out] ca1       Order-1 kind tensors of every active atom
   !> @param[out] ca2       Order-2 kind tensors of every active atom
   !> @param[out] ca3       Order-3 kind tensors of every active atom
   !> @param[out] ca4       Order-4 kind tensors of every active atom
   subroutine pair_atom_cache(self, max_deriv, ca0, ca1, ca2, ca3, ca4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Highest spatial-derivative order to fill
      integer, intent(in) :: max_deriv
      !> Order-0 kind tensors
      real(wp), allocatable, intent(out) :: ca0(:, :)
      !> Order-1 kind tensors
      real(wp), allocatable, intent(out) :: ca1(:, :, :)
      !> Order-2 kind tensors
      real(wp), allocatable, intent(out) :: ca2(:, :, :, :)
      !> Order-3 kind tensors
      real(wp), allocatable, intent(out) :: ca3(:, :, :, :, :)
      !> Order-4 kind tensors
      real(wp), allocatable, intent(out) :: ca4(:, :, :, :, :, :)

      !> Active-list index and active-atom count
      integer :: ia, n

      n = self%n_active
      allocate (ca0(nkind, n))
      allocate (ca1(ndim, nkind, n))
      allocate (ca2(ndim, ndim, nkind, n))
      allocate (ca3(ndim, ndim, ndim, nkind, n))
      allocate (ca4(ndim, ndim, ndim, ndim, nkind, n))

      do ia = 1, n
         call atom_tensors(self, ia, max_deriv, ca0(:, ia), ca1(:, :, ia), &
                           ca2(:, :, :, ia), ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia))
      end do
   end subroutine pair_atom_cache

   !* ================================================================================= *!
   !*                             Pure spatial derivatives                              *!
   !* ================================================================================= *!

   !> Level-set value only
   !>
   !> @param[in]  self LSF instance
   !> @param[out] val  LSF value at the prepared point
   subroutine lsf_f0(self, val)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out) :: val

      !> Blending weights and unused higher-order kernel outputs
      real(wp) :: s_1, s_2, s_3
      real(wp) :: d1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)

      val = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(0, "f0")

      call svdw_weights(self, s_1, s_2, s_3)
      call svdw_spatial_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4, 0, &
                             val, d1, d2, d3, d4)
   end subroutine lsf_f0

   !> Level-set value, spatial gradient and spatial Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    LSF value (optional)
   !> @param[out] lsf1_r  Spatial gradient (optional)
   !> @param[out] lsf2_rr Spatial Hessian (optional)
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      !> Blending weights and kernel outputs
      real(wp) :: s_1, s_2, s_3, f0
      real(wp) :: f1_r(ndim), f2_rr(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Highest spatial-derivative order requested
      integer :: md

      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%n_active == 0) return

      md = 0
      if (present(lsf1_r)) md = 1
      if (present(lsf2_rr)) md = 2
      call self%require_deriv(md, "f012_r")

      call svdw_weights(self, s_1, s_2, s_3)
      call svdw_spatial_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4, md, &
                             f0, f1_r, f2_rr, d3, d4)

      if (present(lsf0)) lsf0 = f0
      if (present(lsf1_r)) lsf1_r = f1_r
      if (present(lsf2_rr)) lsf2_rr = f2_rr
   end subroutine lsf_f012_r

   !> Third spatial derivative, plus optionally the lower orders
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     LSF value (optional)
   !> @param[out] lsf1_r   Spatial gradient (optional)
   !> @param[out] lsf2_rr  Spatial Hessian (optional)
   !> @param[out] lsf3_rrr Third spatial derivative [3, 3, 3]
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      !> Third spatial derivative
      real(wp), intent(out) :: lsf3_rrr(:, :, :)

      !> Blending weights and kernel outputs
      real(wp) :: s_1, s_2, s_3, f0
      real(wp) :: f1_r(ndim), f2_rr(ndim, ndim), f3_rrr(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)

      lsf3_rrr = 0.0_wp
      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(3, "f3_rrr")

      call svdw_weights(self, s_1, s_2, s_3)
      call svdw_spatial_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4, 3, &
                             f0, f1_r, f2_rr, f3_rrr, d4)

      if (present(lsf0)) lsf0 = f0
      if (present(lsf1_r)) lsf1_r = f1_r
      if (present(lsf2_rr)) lsf2_rr = f2_rr
      lsf3_rrr = f3_rrr
   end subroutine lsf_f3_rrr

   !> Fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative [3, 3, 3, 3]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Fourth spatial derivative
      real(wp), intent(out) :: lsf4_rrrr(:, :, :, :)

      !> Blending weights and kernel outputs
      real(wp) :: s_1, s_2, s_3, f0
      real(wp) :: f1_r(ndim), f2_rr(ndim, ndim), f3_rrr(ndim, ndim, ndim)
      real(wp) :: f4_rrrr(ndim, ndim, ndim, ndim)

      lsf4_rrrr = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(4, "f4_rrrr")

      call svdw_weights(self, s_1, s_2, s_3)
      call svdw_spatial_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4, 4, &
                             f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
      lsf4_rrrr = f4_rrrr
   end subroutine lsf_f4_rrrr

   !* ================================================================================= *!
   !*                          One-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Mixed third derivative d^3S / (dr^2 dR_A) and its lower orders
   !>
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`.
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    dS/dR_A [3, >= active_count()] (optional)
   !> @param[out] lsf2_r_rA  d^2S/(dr dR_A) [3, 3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rA d^3S/(dr^2 dR_A) [3, 3, 3, >= active_count()]
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear gradient
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(2, "f3_rr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 2, &
                                f1_rA, f2_r_rA, f3_rr_rA, d4)
         if (present(lsf1_rA)) lsf1_rA(:, ia) = f1_rA
         if (present(lsf2_r_rA)) lsf2_r_rA(:, :, ia) = f2_r_rA
         lsf3_rr_rA(:, :, :, ia) = f3_rr_rA
      end do
   end subroutine lsf_f3_rr_rA

   !* ================================================================================= *!
   !*                              Radius derivatives                                   *!
   !* ================================================================================= *!

   !> Radius derivative d^3S / (dr^2 dR_a) and its lower orders
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf1_rad     dS/dR_a [>= active_count()] (optional)
   !> @param[out] lsf2_r_rad   d^2S/(dr dR_a) [3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rad  d^3S/(dr^2 dR_a) [3, 3, >= active_count()]
   subroutine lsf_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Radius gradient
      real(wp), intent(out), optional :: lsf1_rad(:)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_rr_rad(:, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rad, f2_r_rad(ndim), f3_rr_rad(ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         ! Nothing is active, so no slot is owned and "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather
         ! than return bare: the results are `intent(out)`, so a bare return
         ! hands the caller a buffer it is not allowed to read.
         if (present(lsf1_rad)) lsf1_rad = 0.0_wp
         if (present(lsf2_r_rad)) lsf2_r_rad = 0.0_wp
         lsf3_rr_rad = 0.0_wp
         return
      end if
      call self%require_deriv(2, "f3_rr_rad")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         ! `d/dR_a` consumes no derivative order, so order 2 is all that is
         ! needed here -- one order less than `f3_rr_rA`. `atom_tensors` declares
         ! `at3` `intent(out)` and then leaves it unset at `max_deriv = 2`;
         ! `svdw_radius_eval` never reads it, but it is still passed by
         ! reference, so define it rather than hand an `intent(in)` dummy
         ! undefined memory.
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         at3 = 0.0_wp
         call svdw_radius_eval(self%param%blend_k, s_1, s_2, s_3, &
                               self%ps0, self%ps1, self%ps2, self%ps3, &
                               at0, at1, at2, at3, 2, &
                               f1_rad, f2_r_rad, f3_rr_rad, d4)
         if (present(lsf1_rad)) lsf1_rad(ia) = f1_rad
         if (present(lsf2_r_rad)) lsf2_r_rad(:, ia) = f2_r_rad
         lsf3_rr_rad(:, :, ia) = f3_rr_rad
      end do
   end subroutine lsf_f3_rr_rad

   !> Mixed fourth derivative d^4S / (dr^3 dR_A), active-indexed
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA d^4S/(dr^3 dR_A) [3, 3, 3, 3, >= active_count()]
   subroutine lsf_f4_rrr_rA(self, lsf4_rrr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: f4_rrr_rA(ndim, ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(3, "f4_rrr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 4, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 3, &
                                f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
         lsf4_rrr_rA(:, :, :, :, ia) = f4_rrr_rA
      end do
   end subroutine lsf_f4_rrr_rA

   !> Normalized level set S/||grad S|| and its nuclear gradient
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] val      S/||grad S||
   !> @param[out] deriv_rA d/dR_A of S/||grad S|| [3, >= active_count()] (optional)
   subroutine lsf_normalized_f01_rA(self, val, deriv_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Normalized level-set value
      real(wp), intent(out) :: val
      !> Nuclear gradient of the normalized value
      real(wp), intent(out), optional :: deriv_rA(:, :)

      !> Blending weights and pure spatial outputs
      real(wp) :: s_1, s_2, s_3, f0
      real(wp) :: f1_r(ndim), f2_rr(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> One-nuclear-index outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), d3n(ndim, ndim, ndim)
      real(wp) :: d4n(ndim, ndim, ndim, ndim)
      !> Normalized outputs of one atom
      real(wp) :: norm0, norm1_rA(ndim)
      !> Active-list index
      integer :: ia

      val = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(1, "normalized_f01_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call svdw_spatial_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, self%ps4, 1, &
                             f0, f1_r, f2_rr, d3, d4)

      ! The value half needs no nuclear input; the zero dummies keep the single
      ! shared gradient-norm branch of the kernel in charge of the degenerate case
      f1_rA = 0.0_wp
      f2_r_rA = 0.0_wp
      call svdw_normalized_eval(f0, f1_r, f1_rA, f2_r_rA, norm0, norm1_rA)
      val = norm0
      if (.not. present(deriv_rA)) return

      do ia = 1, self%n_active
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 1, &
                                f1_rA, f2_r_rA, d3n, d4n)
         call svdw_normalized_eval(f0, f1_r, f1_rA, f2_r_rA, norm0, norm1_rA)
         deriv_rA(:, ia) = norm1_rA
      end do
   end subroutine lsf_normalized_f01_rA

   !* ================================================================================= *!
   !*                          Two-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Pure nuclear Hessian d^2S / (dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf2_rArB d^2S/(dR_A dR_B) [3, >= n_act, 3, >= n_act]
   subroutine lsf_f2_rArB(self, lsf2_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear Hessian
      real(wp), intent(out) :: lsf2_rArB(:, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices and active-atom count
      integer :: ia, ib, n

      if (self%n_active == 0) return
      call self%require_deriv(0, "f2_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 2, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        0, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 0, b2, b3, b4)
            end if
            lsf2_rArB(:, ia, :, ib) = b2
         end do
      end do
   end subroutine lsf_f2_rArB

   !> Mixed third derivative d^3S / (dr dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf3_r_rArB d^3S/(dr dR_A dR_B) [3, 3, >= n_act, 3, >= n_act]
   subroutine lsf_f3_r_rArB(self, lsf3_r_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rArB(:, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices, active-atom count and nuclear components
      integer :: ia, ib, n, s, t

      if (self%n_active == 0) return
      call self%require_deriv(1, "f3_r_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 3, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        1, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 1, b2, b3, b4)
            end if
            do t = 1, ndim
               do s = 1, ndim
                  lsf3_r_rArB(:, s, ia, t, ib) = b3(:, s, t)
               end do
            end do
         end do
      end do
   end subroutine lsf_f3_r_rArB

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf4_rr_rArB d^4S/(dr^2 dR_A dR_B) [3, 3, 3, >= n_act, 3, >= n_act]
   subroutine lsf_f4_rr_rArB(self, lsf4_rr_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices, active-atom count and tensor components
      integer :: ia, ib, n, s, t, j, k

      if (self%n_active == 0) return
      call self%require_deriv(2, "f4_rr_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 4, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        2, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 2, b2, b3, b4)
            end if
            do t = 1, ndim
               do s = 1, ndim
                  do k = 1, ndim
                     do j = 1, ndim
                        lsf4_rr_rArB(j, k, s, ia, t, ib) = b4(j, k, s, t)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine lsf_f4_rr_rArB

   !* ================================================================================= *!
   !*                     Two-radius and nuclear-radius derivatives                     *!
   !* ================================================================================= *!
   
   !> Shared pair loop of the two-radius block
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Spatial-derivative order of the requested output
   !> @param[in]  caller    Name reported by `require_deriv` on an under-request
   !> @param[out] res2      d^2S/(dR_a dR_b) (optional)
   !> @param[out] res3      d^3S/(dr dR_a dR_b) (optional)
   !> @param[out] res4      d^4S/(dr^2 dR_a dR_b) (optional)
   subroutine svdw_radpair_block(self, max_deriv, caller, res2, res3, res4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Spatial-derivative order of the requested output
      integer, intent(in) :: max_deriv
      !> Name reported on an under-request
      character(len=*), intent(in) :: caller
      !> Radius Hessian
      real(wp), intent(out), optional :: res2(:, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: res3(:, :, :)
      !> Mixed fourth derivative
      real(wp), intent(out), optional :: res4(:, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair; a radius slot adds no Cartesian index
      real(wp) :: b2, b3(ndim), b4(ndim, ndim)
      !> Active-list indices and active-atom count
      integer :: ia, ib, n

      if (self%n_active == 0) then
         ! Nothing is active, so no slot is owned and "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather
         ! than return bare: the results are `intent(out)`, so a bare return
         ! hands the caller a buffer it is not allowed to read.
         if (present(res2)) res2 = 0.0_wp
         if (present(res3)) res3 = 0.0_wp
         if (present(res4)) res4 = 0.0_wp
         return
      end if
      call self%require_deriv(max_deriv, caller)

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, max_deriv, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_radpair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                           self%ps0, self%ps1, self%ps2, &
                                           ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                           max_deriv, b2, b3, b4)
            else
               call svdw_radpair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                      self%ps0, self%ps1, self%ps2, &
                                      ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                      ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                      max_deriv, b2, b3, b4)
            end if
            if (present(res2)) res2(ia, ib) = b2
            if (present(res3)) res3(:, ia, ib) = b3
            if (present(res4)) res4(:, :, ia, ib) = b4
         end do
      end do
   end subroutine svdw_radpair_block

   !> Shared pair loop of the nuclear-radius block
   !>
   !> The first atom index carries the *position* derivative and the second the
   !> *radius* one, so unlike [[svdw_radpair_block]] the result is not symmetric
   !> under exchanging them; the loop therefore visits every ordered pair.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Spatial-derivative order of the requested output
   !> @param[in]  caller    Name reported by `require_deriv` on an under-request
   !> @param[out] res2      d^2S/(dR_A dR_b) (optional)
   !> @param[out] res3      d^3S/(dr dR_A dR_b) (optional)
   !> @param[out] res4      d^4S/(dr^2 dR_A dR_b) (optional)
   subroutine svdw_nucrad_block(self, max_deriv, caller, res2, res3, res4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Spatial-derivative order of the requested output
      integer, intent(in) :: max_deriv
      !> Name reported on an under-request
      character(len=*), intent(in) :: caller
      !> Nuclear-radius Hessian
      real(wp), intent(out), optional :: res2(:, :, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: res3(:, :, :, :)
      !> Mixed fourth derivative
      real(wp), intent(out), optional :: res4(:, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair; only the position slot adds an index
      real(wp) :: b2(ndim), b3(ndim, ndim), b4(ndim, ndim, ndim)
      !> Active-list indices and active-atom count
      integer :: ia, ib, n

      if (self%n_active == 0) then
         ! Nothing is active, so no slot is owned and "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather
         ! than return bare: the results are `intent(out)`, so a bare return
         ! hands the caller a buffer it is not allowed to read.
         if (present(res2)) res2 = 0.0_wp
         if (present(res3)) res3 = 0.0_wp
         if (present(res4)) res4 = 0.0_wp
         return
      end if
      call self%require_deriv(max_deriv, caller)

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, max_deriv + 1, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_nucrad_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                          self%ps0, self%ps1, self%ps2, &
                                          ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                          ca3(:, :, :, :, ia), max_deriv, b2, b3, b4)
            else
               call svdw_nucrad_eval(self%param%blend_k, s_1, s_2, s_3, &
                                     self%ps0, self%ps1, self%ps2, &
                                     ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                     ca3(:, :, :, :, ia), &
                                     ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                     max_deriv, b2, b3, b4)
            end if
            if (present(res2)) res2(:, ia, ib) = b2
            if (present(res3)) res3(:, :, ia, ib) = b3
            if (present(res4)) res4(:, :, :, ia, ib) = b4
         end do
      end do
   end subroutine svdw_nucrad_block

   !> Pure radius Hessian d^2S / (dR_a dR_b), active-indexed in both atoms
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_radrad d^2S/(dR_a dR_b) [>= n_act, >= n_act]
   subroutine lsf_f2_radrad(self, lsf2_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Radius Hessian
      real(wp), intent(out) :: lsf2_radrad(:, :)

      call svdw_radpair_block(self, 0, "f2_radrad", res2=lsf2_radrad)
   end subroutine lsf_f2_radrad

   !> Mixed third derivative d^3S / (dr dR_a dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_radrad d^3S/(dr dR_a dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f3_r_radrad(self, lsf3_r_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_radrad(:, :, :)

      call svdw_radpair_block(self, 1, "f3_r_radrad", res3=lsf3_r_radrad)
   end subroutine lsf_f3_r_radrad

   !> Mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_f4_rr_radrad(self, lsf4_rr_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :, :)

      call svdw_radpair_block(self, 2, "f4_rr_radrad", res4=lsf4_rr_radrad)
   end subroutine lsf_f4_rr_radrad

   !> Mixed nuclear-radius Hessian d^2S / (dR_A dR_b)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_rA_rad d^2S/(dR_A dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f2_rA_rad(self, lsf2_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear-radius Hessian
      real(wp), intent(out) :: lsf2_rA_rad(:, :, :)

      call svdw_nucrad_block(self, 0, "f2_rA_rad", res2=lsf2_rA_rad)
   end subroutine lsf_f2_rA_rad

   !> Mixed third derivative d^3S / (dr dR_A dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_rA_rad d^3S/(dr dR_A dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_f3_r_rA_rad(self, lsf3_r_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rA_rad(:, :, :, :)

      call svdw_nucrad_block(self, 1, "f3_r_rA_rad", res3=lsf3_r_rA_rad)
   end subroutine lsf_f3_r_rA_rad

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) [3, 3, 3, >= n_act, >= n_act]
   subroutine lsf_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :, :)

      call svdw_nucrad_block(self, 2, "f4_rr_rA_rad", res4=lsf4_rr_rA_rad)
   end subroutine lsf_f4_rr_rA_rad

   !* ================================================================================= *!
   !*                     Direction-contracted nuclear derivatives                      *!
   !* ================================================================================= *!

   !> Directional nuclear derivative of the level-set value
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . dS/dR_B
   subroutine lsf_tangent_f0(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted value derivative
      real(wp), intent(out) :: res

      !> Blending weights, contracted power sums and unused kernel outputs
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: d1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(0, "tangent_f0")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 0, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 0, res, d1, d2, d3)
   end subroutine lsf_tangent_f0

   !> Directional nuclear derivative of the spatial gradient
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of dS/dr [3]
   subroutine lsf_tangent_f1_r(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted gradient derivative
      real(wp), intent(out) :: res(:)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(1, "tangent_f1_r")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 1, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 1, t0, t1, d2, d3)
      res = t1
   end subroutine lsf_tangent_f1_r

   !> Directional nuclear derivative of the spatial Hessian
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^2S/dr^2 [3, 3]
   subroutine lsf_tangent_f2_rr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted Hessian derivative
      real(wp), intent(out) :: res(:, :)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), t2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(2, "tangent_f2_rr")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 2, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 2, t0, t1, t2, d3)
      res = t2
   end subroutine lsf_tangent_f2_rr

   !> Directional nuclear derivative of the third spatial derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^3S/dr^3 [3, 3, 3]
   subroutine lsf_tangent_f3_rrr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted third-derivative derivative
      real(wp), intent(out) :: res(:, :, :)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(3, "tangent_f3_rrr")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 3, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 3, t0, t1, t2, t3)
      res = t3
   end subroutine lsf_tangent_f3_rrr

   !> Nuclear Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^2S/(dR_A dR_B) [3, >= active_count()]
   !> @param[in]  vrad Radius directions [ncenters] (optional). Supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, turning `res` from the nuclear-nuclear
   !>                  block into the nuclear row of the joint Hessian
   subroutine lsf_hvp_f1_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted nuclear Hessian
      real(wp), intent(out) :: res(:, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(0, "hvp_f1_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 0, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 1, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 1, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 1, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 0, h1, d2, d3)
         res(:, ia) = h1
      end do
   end subroutine lsf_hvp_f1_rA

   !> Directional nuclear derivative of `f2_r_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^3S/(dr dR_A dR_B) [3, 3, >= active_count()]
   !> @param[in]  vrad Radius directions [ncenters] (optional). Supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, exactly as for [[lsf_hvp_f1_rA]]
   subroutine lsf_hvp_f2_r_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(1, "hvp_f2_r_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 1, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 2, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 2, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
            aw2 = aw2 + awr2
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 1, h1, h2, d3)
         res(:, :, ia) = h2
      end do
   end subroutine lsf_hvp_f2_r_rA

   !> Directional nuclear derivative of `f3_rr_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^4S/(dr^2 dR_A dR_B) [3, 3, 3, >= active_count()]
   !> @param[in]  vrad Radius directions [ncenters] (optional). Supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, exactly as for [[lsf_hvp_f1_rA]]
   subroutine lsf_hvp_f3_rr_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: res(:, :, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(2, "hvp_f3_rr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 2, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 3, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 3, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
            aw2 = aw2 + awr2
            aw3 = aw3 + awr3
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 2, h1, h2, h3)
         res(:, :, :, ia) = h3
      end do
   end subroutine lsf_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                    Radius row of the joint Hessian-vector product                 *!
   !* ================================================================================= *!

   !> Radius row of the joint Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  sum_B (v_B . d/dR_B + vr_B d/dR_b) dS/dR_a
   !>                  [>= active_count()]
   subroutine lsf_hvp_f1_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted radius Hessian row
      real(wp), intent(out) :: res(:)

      call svdw_radius_hvp(self, v, vrad, 0, "hvp_f1_rad", res1=res)
   end subroutine lsf_hvp_f1_rad

   !> Joint directional derivative of `f2_r_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed third derivative [3, >= active_count()]
   subroutine lsf_hvp_f2_r_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :)

      call svdw_radius_hvp(self, v, vrad, 1, "hvp_f2_r_rad", res2=res)
   end subroutine lsf_hvp_f2_r_rad

   !> Joint directional derivative of `f3_rr_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed fourth derivative [3, 3, >= active_count()]
   subroutine lsf_hvp_f3_rr_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: res(:, :, :)

      call svdw_radius_hvp(self, v, vrad, 2, "hvp_f3_rr_rad", res3=res)
   end subroutine lsf_hvp_f3_rr_rad

   !> Shared driver of the three radius-row Hessian-vector products
   !>
   !> The three public entry points differ only in `max_deriv` and in which
   !> output they keep, and unlike the `_rA` ladder there is no reason to spell
   !> the loop out three times: the radius channel has one code path. Exactly one
   !> of `res1`/`res2`/`res3` must be present, matching `max_deriv`.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad      Radius directions [ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order (0..2)
   !> @param[in]  caller    Name used by `require_deriv` on failure
   !> @param[out] res1      Order-0 result [>= active_count()] (optional)
   !> @param[out] res2      Order-1 result [3, >= active_count()] (optional)
   !> @param[out] res3      Order-2 result [3, 3, >= active_count()] (optional)
   subroutine svdw_radius_hvp(self, v, vrad, max_deriv, caller, res1, res2, res3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Caller name for the derivative-order check
      character(len=*), intent(in) :: caller
      !> Order-0 result
      real(wp), intent(out), optional :: res1(:)
      !> Order-1 result
      real(wp), intent(out), optional :: res2(:, :)
      !> Order-2 result
      real(wp), intent(out), optional :: res3(:, :, :)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1, h2(ndim), h3(ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         ! Nothing is active, so no slot is owned and "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather
         ! than return bare: the results are `intent(out)`, so a bare return
         ! hands the caller a buffer it is not allowed to read.
         if (present(res1)) res1 = 0.0_wp
         if (present(res2)) res2 = 0.0_wp
         if (present(res3)) res3 = 0.0_wp
         return
      end if
      call self%require_deriv(max_deriv, caller)

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, max_deriv, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         ! `d/dR_a` adds no index, so unlike the `_rA` ladder the per-atom
         ! tensors are needed only to `max_deriv`, not `max_deriv + 1`.
         call atom_tensors(self, ia, max_deriv, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), max_deriv, &
                                   aw0, aw1, aw2, aw3)
         call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), &
                                          max_deriv, awr0, awr1, awr2, awr3)
         aw0 = aw0 + awr0
         if (max_deriv >= 1) aw1 = aw1 + awr1
         if (max_deriv >= 2) aw2 = aw2 + awr2
         ! `at3`/`at4`/`aw3`/`awr3` exist only because the three tensor
         ! builders take them as non-optional `intent(out)` dummies; the eval
         ! below never sees them, so they are left undefined on purpose.
         call svdw_radius_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                                   at0, at1, at2, aw0, aw1, aw2, max_deriv, &
                                   h1, h2, h3)
         if (present(res1)) res1(ia) = h1
         if (present(res2)) res2(:, ia) = h2
         if (present(res3)) res3(:, :, ia) = h3
      end do
   end subroutine svdw_radius_hvp

   !* ================================================================================= *!
   !*                        Jet-contracted nuclear derivative                          *!
   !* ================================================================================= *!

   !> Adjoint-weighted nuclear gradient, active-indexed
   !>
   !> Contracts one evaluation point's adjoint weights against the nuclear
   !> ladder and keeps only the nuclear index:
   !>
   !>    res(s, i) = w0*lsf1_rA(s, i) + sum_a w1(a)*lsf2_r_rA(a, s, i)
   !>                + sum_a sum_b w2(a, b)*lsf3_rr_rA(a, b, s, i) .
   !>
   !> This is the `f1_rA` rung with the jet indices contracted away, exactly as
   !> [[lsf_hvp_f1_rA]] is that same rung contracted with a nuclear direction --
   !> the reverse-mode mirror of the `tangent_*` family, which contracts the
   !> nuclear index and keeps the spatial ones.
   !>
   !> `w2` is a general 3x3: all nine entries are contracted, no symmetry is
   !> assumed and no factor of two is folded into the off-diagonals.
   !>
   !> The contraction happens inside the kernel, on the weighted jet rather than
   !> on its 3 + 9 + 27 components, so an adjoint pass never has to materialize
   !> the `(3, 3, 3, >= active_count())` tensor [[lsf_f3_rr_rA]] would hand it.
   !> The prepared order is the same as for [[lsf_f3_rr_rA]], since the same
   !> `lsf3_rr_rA` rung enters the contraction.
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the level-set value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Adjoint weights of the spatial Hessian [3, 3]
   !> @param[out] res  Jet-contracted nuclear gradient [3, >= active_count()]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)
      !> Jet-contracted nuclear gradient
      real(wp), intent(out) :: res(:, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel output of one atom
      real(wp) :: vjp_f1_rA(ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(2, "vjp_f1_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call svdw_vjp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, &
                            at0, at1, at2, at3, w0, w1, w2, vjp_f1_rA)
         res(:, ia) = vjp_f1_rA
      end do
   end subroutine lsf_vjp_f1_rA

   !> Adjoint-weighted radius gradient, active-indexed
   !>
   !> The radius twin of [[lsf_vjp_f1_rA]]: the same evaluation point's adjoint
   !> weights, contracted against the radius ladder instead of the nuclear one,
   !>
   !>    res(i) = w0*lsf1_rad(i) + sum_a w1(a)*lsf2_r_rad(a, i)
   !>             + sum_a sum_b w2(a, b)*lsf3_rr_rad(a, b, i) .
   !>
   !> A radius is a scalar, so -- exactly as [[lsf_f3_rr_rad]] carries one rank
   !> less than [[lsf_f3_rr_rA]] -- there is no index left once the jet indices
   !> are contracted away: one number per active atom, against the three of
   !> [[lsf_vjp_f1_rA]] and the 1 + 3 + 9 of [[lsf_f3_rr_rad]].
   !>
   !> `w2` is a general 3x3: all nine entries are contracted, no symmetry is
   !> assumed and no factor of two is folded into the off-diagonals.
   !>
   !> The contraction happens inside the kernel, on the weighted jet rather than
   !> on its components. The prepared order is the same as for [[lsf_f3_rr_rad]]:
   !> `d/dR_a` consumes no derivative order, so the kernel reads the power sums
   !> and per-atom tensors only up to order 2.
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the level-set value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Adjoint weights of the spatial Hessian [3, 3]
   !> @param[out] res  Jet-contracted radius gradient [>= active_count()]
   subroutine lsf_vjp_f1_rad(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)
      !> Jet-contracted radius gradient
      real(wp), intent(out) :: res(:)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel output of one atom
      real(wp) :: vjp_f1_rad
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(2, "vjp_f1_rad")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         ! Order 2 is enough: the radius derivative bumps no index, so the kernel
         ! never reaches `at3`/`at4` and does not take them as arguments -- unlike
         ! `svdw_radius_eval`, which does and therefore needs `at3` defined.
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         call svdw_radius_vjp_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   at0, at1, at2, w0, w1, w2, vjp_f1_rad)
         res(ia) = vjp_f1_rad
      end do
   end subroutine lsf_vjp_f1_rad

   !* ================================================================================= *!
   !*                                Screening offset                                   *!
   !* ================================================================================= *!

   !> Exact radial offset where the SvdW screening weight equals the threshold
   !>
   !> An atom at distance `x` from the evaluation point enters the blend with weight
   !> `w = exp(-(k/3) * (x - R))`, so `w = threshold` at
   !>
   !>    `x - R = -3*ln(threshold)/k`
   !>
   !> which is independent of the radius.
   !>
   !> This is the *screening* bound query. The companion `exclusion_radius` below is
   !> the *surface* bound query; see the module header of
   !> [[moist_cavity_drop_lsf_base]] for why the two cannot be expressed in terms of
   !> each other.
   !>
   !> @param[in] self    SvdW LSF instance (reads param%blend_k, screening_threshold)
   !> @param[in] radius  Atom radius (Bohr, unused: the SvdW offset is radius-independent)
   !> @returns   offset  Radial offset from the atom surface (Bohr)
   pure function lsf_screening_offset(self, radius) result(offset)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from atom surface (Bohr)
      real(wp) :: offset

      !> Blending sharpness and screening threshold
      real(wp) :: k_local, threshold

      k_local = self%param%blend_k
      threshold = self%screening_threshold

      ! With threshold <= 0, k <= 0 or a nonsensical radius the criterion is
      ! ill-defined; report an unbounded reach so screening is effectively disabled.
      if (threshold <= 0.0_wp .or. k_local <= 0.0_wp .or. radius < 0.0_wp) then
         offset = huge(0.0_wp)
         return
      end if

      offset = max(0.0_wp, -3.0_wp*log(threshold)/k_local)
   end function lsf_screening_offset

   !* ================================================================================= *!
   !*                        Surface-free exclusion certificate                         *!
   !* ================================================================================= *!

   !> Exact surface-free radius from the 1-Lipschitz property
   !>
   !> `S = -(1/k) ln Z` is a log-sum-exp soft minimum of atomic signed
   !> distances. Its gradient is the blend-weighted mean of the gradients of
   !> the per-term exponents, and each of those is itself an average of unit
   !> vectors `(r - r_I)/||r - r_I||`, so every term contributes a vector of
   !> norm at most one. The blend weights sum to one, hence `||grad S|| <= 1`
   !> and `S` cannot reach zero from `S(x)` in less than `|S(x)|`: the ball
   !> `B(x, |S(x)|)` contains no point of the zero level set. This is the
   !> tightest radius a Lipschitz argument can give, and it costs nothing
   !> beyond the value the caller already holds.
   !>
   !> The bound needs the blend weights to be non-negative, which holds exactly
   !> when the many-body coefficients are. A negative coefficient turns the
   !> weighted mean into an extrapolation and the bound is lost, so such a
   !> parameterization certifies nothing and gets the safe answer of zero.
   !>
   !> @param[in] self LSF instance
   !> @param[in] lsf0 LSF value at the evaluation point
   !> @returns   r    Radius of a ball around the point free of surface
   pure function lsf_svdw_exclusion_radius(self, lsf0) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value at the evaluation point
      real(wp), intent(in) :: lsf0
      !> Surface-free radius
      real(wp) :: r

      r = 0.0_wp
      if (self%param%blend_1b < 0.0_wp) return
      if (self%param%blend_2b < 0.0_wp) return
      if (self%param%blend_3b < 0.0_wp) return

      r = abs(lsf0)
   end function lsf_svdw_exclusion_radius

   !* ================================================================================= *!
   !*                                    Finalizer                                      *!
   !* ================================================================================= *!

   !> Release the allocatable components of an SvdW LSF instance
   !>
   !> @param[inout] self LSF instance being finalized
   subroutine finalize_lsf_svdw(self)
      !> LSF instance
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self

      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%active_cand)) deallocate (self%active_cand)
      if (allocated(self%act_atom)) deallocate (self%act_atom)
      if (allocated(self%act_d)) deallocate (self%act_d)
      if (allocated(self%act_radius)) deallocate (self%act_radius)
      if (allocated(self%act_x)) deallocate (self%act_x)

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

end module moist_cavity_drop_lsf_svdw
