!> COSMO Fine Cavity (CFC) level-set function
!>
!> Concrete LSF implementing the Diedenhofen & Klamt 2018 pseudo-density
!>   PD(r) = sum_a   exp{ a1 (s_a - 1) }                                     (atomic)
!>         + sum_{a<b} c (1 - vec(s_a).vec(s_b))^m exp{ a2 (s_a + s_b - 2) } (pair)
!>   with s_a = ||r - r_a|| / R_a and vec(s_a) = (r - r_a) / R_a.
!>
!> The level set returned by this LSF is
!>   PD*(r) = -log PD(r)
!>
!> Thus matching the SvdW sign convention (interior negative)
!>
!> The math kernel itself lives in the code-generated module [[moist_cavity_drop_lsf_cfc_kernel]]
!>
!> TODO: (screening): the SSD uses `screen_k = 3`, which matches the pair-term decay rate (|a2|/3 = 3)
!> but is *over*-conservative for the atomic-term decay (|a1|/3 = 5)
module moist_cavity_drop_lsf_cfc
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw_ssd, only: moist_cavity_drop_lsf_svdw_ssd_type, ssd0
   use moist_cavity_drop_lsf_cfc_kernel, only: cfc_atomic_term_eval, &
                                               cfc_pair_term_eval, cfc_log_lift
   implicit none
   private

   integer, parameter :: ndim = 3

   !> Diedenhofen-Klamt 2018 defaults for the four CFC shape parameters.
   real(wp), parameter :: a1_default = -15.0_wp
   real(wp), parameter :: a2_default = -9.0_wp
   real(wp), parameter :: c_default = 5.0_wp
   integer, parameter :: m_default = 4
   !> Default SSD screening sharpness (sized so the SSD's exp(-k/3 * d)
   !> filter drops faster than the loosest CFC contribution). The CFC
   !> neighbour cutoff is the analytic source of truth; this k_screen
   !> only controls the SSD's coarse pre-cull and need not be tuned.
   real(wp), parameter :: screen_k_default = 3.0_wp

   !> COSMO Fine Cavity LSF: radii-based pseudo-density level set.
   !>
   !> Holds the four CFC shape parameters, a screening SSD system, and
   !> per-evaluation-point caches for PD and its spatial derivatives plus
   !> per-atom accumulators used by the nuclear-derivative routine.
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_cfc_type
      !> Atomic-term exponent (Diedenhofen-Klamt 2018: -15)
      real(wp) :: a1 = a1_default
      !> Pair-term exponent (Diedenhofen-Klamt 2018: -9)
      real(wp) :: a2 = a2_default
      !> Pair-term coupling constant (Diedenhofen-Klamt 2018: 5)
      real(wp) :: c = c_default
      !> Pair-term polynomial power (Diedenhofen-Klamt 2018: 4)
      integer  :: m = m_default

      !> SSD screening sharpness (k in `exp(-k/3 * d_I)`)
      !> NOT a CFC math parameter; only steers the SSD's coarse atom cull
      real(wp) :: screen_k = screen_k_default

      !> Highest spatial-derivative order the caches are sized for
      integer  :: max_deriv = 0

      !> Screening SSD: reused unchanged from the SvdW LSF
      type(moist_cavity_drop_lsf_svdw_ssd_type) :: ssd_system

      !> Cached `[1..ncenters]` for `prepare()` to feed to `ssd_system%compute`
      integer, allocatable :: all_indices(:)

      ! ------- Per-evaluation-point caches (filled by prepare) ---------- *!

      !> Pseudo-density value at the cached evaluation point
      real(wp) :: PD0 = 0.0_wp
      !> Pseudo-density spatial gradient at the cached point [ndim]
      real(wp), allocatable :: PD1(:)
      !> Pseudo-density spatial Hessian [ndim, ndim]
      real(wp), allocatable :: PD2(:, :)
      !> Pseudo-density third spatial derivative [ndim, ndim, ndim]
      real(wp), allocatable :: PD3(:, :, :)

      ! --- Per-atom accumulators for nuclear-derivative assembly ------- *!

      !> Q1^A(alpha) = dPD/d(d_A,alpha) for atom A (collected from atom-A
      !> self-term and every pair involving A) [ndim, n_alloc]
      real(wp), allocatable :: Q1(:, :)
      !> Q2^A(alpha, j) = d^2 PD/(d(d_A,alpha) d r_j) per atom; [ndim, ndim, n_alloc]
      real(wp), allocatable :: Q2(:, :, :)
      !> Q3^A(alpha, j, k) = d^3 PD/(d(d_A,alpha) d r_j d r_k) per atom [ndim, ndim, ndim, n_alloc]
      real(wp), allocatable :: Q3(:, :, :, :)
   contains
      procedure, public :: new => lsf_new
      procedure, public :: update => lsf_update
      procedure, public :: prepare => lsf_prepare
      procedure, public :: prepare_subset => lsf_prepare_subset
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      procedure, public :: active_count => lsf_active_count
      procedure, public :: active_atom => lsf_active_atom
      procedure, public :: f0_screened => lsf_f0_screened
      procedure, public :: f012_r_screened => lsf_f012_r_screened
      procedure, public :: f3_rrr_screened => lsf_f3_rrr_screened
      procedure, public :: f3_rr_rA_screened => lsf_f3_rr_rA_screened
      procedure, public :: neighbor_cutoff => lsf_neighbor_cutoff
      final :: finalize_lsf_cfc
   end type moist_cavity_drop_lsf_cfc_type

   public :: moist_cavity_drop_lsf_cfc_type

contains

   !* ================================================================================= *!
   !*                              LSF lifecycle methods                                *!
   !* ================================================================================= *!

   !> Configure CFC shape parameters
   !>
   !> @param[inout] self      LSF instance
   !> @param[in]    a1        Atomic-term exponent (optional, default -15)
   !> @param[in]    a2        Pair-term exponent (optional, default -9)
   !> @param[in]    c         Pair-term coupling (optional, default 5)
   !> @param[in]    m         Pair-term power (optional, default 4; the
   !>                         kernel hardcodes m=4 so passing m /= 4 is
   !>                         currently ignored)
   !> @param[in]    screen_k  SSD screening sharpness (optional, default 3)
   subroutine lsf_new(self, a1, a2, c, m, screen_k)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Atomic-term exponent override (optional)
      real(wp), intent(in), optional :: a1
      !> Pair-term exponent override (optional)
      real(wp), intent(in), optional :: a2
      !> Pair-term coupling override (optional)
      real(wp), intent(in), optional :: c
      !> Pair-term power override (optional; currently advisory only)
      integer, intent(in), optional :: m
      !> SSD screening sharpness override (optional)
      real(wp), intent(in), optional :: screen_k

      if (present(a1)) self%a1 = a1
      if (present(a2)) self%a2 = a2
      if (present(c)) self%c = c
      if (present(m)) self%m = m
      if (present(screen_k)) self%screen_k = screen_k
   end subroutine lsf_new

   !> Bind molecular geometry and refresh the screening SSD cache
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    mol    Molecular structure
   !> @param[in]    radii  Per-atom radii (size mol%nat)
   subroutine lsf_update(self, mol, radii)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      integer :: i, prior_max_deriv

      self%mol = mol
      self%radii = radii
      self%ncenters = mol%nat

      prior_max_deriv = self%max_deriv
      if (prior_max_deriv < 2) prior_max_deriv = 2

      call self%ssd_system%new( &
         k=self%screen_k, &
         threshold=self%screening_threshold, &
         max_deriv=0)
      call self%ssd_system%update(mol%xyz, radii)

      if (allocated(self%all_indices)) deallocate (self%all_indices)
      allocate (self%all_indices(mol%nat))
      do i = 1, mol%nat
         self%all_indices(i) = i
      end do

      ! Per-atom Q caches are sized by atom count
      if (allocated(self%Q1)) deallocate (self%Q1)
      if (allocated(self%Q2)) deallocate (self%Q2)
      if (allocated(self%Q3)) deallocate (self%Q3)

      call lsf_set_max_deriv(self, prior_max_deriv)
   end subroutine lsf_update

   !> Configure the highest spatial-derivative order to cache
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    n     Requested max derivative order (0..3)
   subroutine lsf_set_max_deriv(self, n)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      integer, intent(in) :: n

      integer :: n_alloc, want
      logical :: have_ssd

      want = n
      if (want < 0) want = 0
      if (want > 3) want = 3
      self%max_deriv = want

      have_ssd = allocated(self%ssd_system%k3f0_arr)
      if (have_ssd) then
         n_alloc = size(self%ssd_system%k3f0_arr)
      else
         n_alloc = 0
      end if

      ! Global PD spatial caches
      if (want >= 1) then
         if (.not. allocated(self%PD1)) allocate (self%PD1(ndim))
      else
         if (allocated(self%PD1)) deallocate (self%PD1)
      end if

      if (want >= 2) then
         if (.not. allocated(self%PD2)) allocate (self%PD2(ndim, ndim))
      else
         if (allocated(self%PD2)) deallocate (self%PD2)
      end if

      if (want >= 3) then
         if (.not. allocated(self%PD3)) allocate (self%PD3(ndim, ndim, ndim))
      else
         if (allocated(self%PD3)) deallocate (self%PD3)
      end if

      ! Per-atom Q caches (gated on the SSD having been allocated)
      if (have_ssd .and. want >= 1) then
         if (.not. allocated(self%Q1)) allocate (self%Q1(ndim, n_alloc))
      else
         if (allocated(self%Q1)) deallocate (self%Q1)
      end if

      if (have_ssd .and. want >= 2) then
         if (.not. allocated(self%Q2)) allocate (self%Q2(ndim, ndim, n_alloc))
      else
         if (allocated(self%Q2)) deallocate (self%Q2)
      end if

      if (have_ssd .and. want >= 3) then
         if (.not. allocated(self%Q3)) allocate (self%Q3(ndim, ndim, ndim, n_alloc))
      else
         if (allocated(self%Q3)) deallocate (self%Q3)
      end if
   end subroutine lsf_set_max_deriv

   !> Run screening at the evaluation point and refresh CFC caches
   subroutine lsf_prepare(self, point)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

      call self%ssd_system%compute(point, self%all_indices)
      call cfc_compute_caches(self, point)
   end subroutine lsf_prepare

   !> Run screening + cache refresh for a caller-provided candidate list
   subroutine lsf_prepare_subset(self, point, candidate_indices)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)

      call self%ssd_system%compute(point, candidate_indices)
      call cfc_compute_caches(self, point)
   end subroutine lsf_prepare_subset

   !> Number of atoms currently active
   pure integer function lsf_active_count(self) result(n)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      n = self%ssd_system%n_active
   end function lsf_active_count

   !> User-space atom id of the i-th active atom
   pure integer function lsf_active_atom(self, i) result(idx)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      integer, intent(in) :: i
      idx = self%ssd_system%atom_indices(i)
   end function lsf_active_atom

   !* ================================================================================= *!
   !*                       Internal cache filler (atom + pair sweep)                   *!
   !* ================================================================================= *!

   !> Populate `PD0..PD3` and `Q1..Q3` from the SSD-supplied active list
   !>
   !> Strategy:
   !>   - Zero the requested-order caches
   !>   - Atomic-term sweep (one pass over active atoms)
   !>   - Pair-term sweep (i<j over the active list)
   !>
   !> @param[inout] self   LSF instance with SSD active list populated
   !> @param[in]    point  Evaluation point (3,)
   subroutine cfc_compute_caches(self, point)
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

      integer  :: n_active, i, j, atom_i, atom_j, max_deriv
      real(wp) :: d_a(ndim), d_b(ndim)
      real(wp) :: R_a, R_b
      real(wp) :: pd0_a, pd1_a_atomic(ndim)
      real(wp) :: pd2_a_atomic(ndim, ndim)
      real(wp) :: pd3_a_atomic(ndim, ndim, ndim)
      real(wp) :: pd0_p
      real(wp) :: p1_a(ndim), p1_b(ndim)
      real(wp) :: p2_aa(ndim, ndim), p2_ab(ndim, ndim), p2_bb(ndim, ndim)
      real(wp) :: p3_aaa(ndim, ndim, ndim), p3_aab(ndim, ndim, ndim)
      real(wp) :: p3_abb(ndim, ndim, ndim), p3_bbb(ndim, ndim, ndim)
      integer  :: ii, jj, kk

      max_deriv = self%max_deriv

      self%PD0 = 0.0_wp
      if (allocated(self%PD1)) self%PD1 = 0.0_wp
      if (allocated(self%PD2)) self%PD2 = 0.0_wp
      if (allocated(self%PD3)) self%PD3 = 0.0_wp
      if (allocated(self%Q1)) self%Q1 = 0.0_wp
      if (allocated(self%Q2)) self%Q2 = 0.0_wp
      if (allocated(self%Q3)) self%Q3 = 0.0_wp

      n_active = self%ssd_system%n_active
      if (n_active == 0) return

      ! ------------------------ Atomic-term sweep ----------------------- *!
      do i = 1, n_active
         atom_i = self%ssd_system%atom_indices(i)
         d_a = point - self%ssd_system%centers(:, atom_i)
         R_a = self%radii(atom_i)

         pd0_a = 0.0_wp
         pd1_a_atomic = 0.0_wp
         pd2_a_atomic = 0.0_wp
         pd3_a_atomic = 0.0_wp

         call cfc_atomic_term_eval(d_a, R_a, self%a1, max_deriv, &
                                   pd0_a, pd1_a_atomic, pd2_a_atomic, pd3_a_atomic)

         self%PD0 = self%PD0 + pd0_a
         if (max_deriv >= 1) then
            self%PD1 = self%PD1 + pd1_a_atomic
            self%Q1(:, i) = self%Q1(:, i) + pd1_a_atomic
         end if
         if (max_deriv >= 2) then
            self%PD2 = self%PD2 + pd2_a_atomic
            self%Q2(:, :, i) = self%Q2(:, :, i) + pd2_a_atomic
         end if
         if (max_deriv >= 3) then
            self%PD3 = self%PD3 + pd3_a_atomic
            self%Q3(:, :, :, i) = self%Q3(:, :, :, i) + pd3_a_atomic
         end if
      end do

      ! ------------------------ Pair-term sweep ------------------------- *!
      do i = 1, n_active
         atom_i = self%ssd_system%atom_indices(i)
         d_a = point - self%ssd_system%centers(:, atom_i)
         R_a = self%radii(atom_i)

         do j = i + 1, n_active
            atom_j = self%ssd_system%atom_indices(j)
            d_b = point - self%ssd_system%centers(:, atom_j)
            R_b = self%radii(atom_j)

            pd0_p = 0.0_wp
            p1_a = 0.0_wp; p1_b = 0.0_wp
            p2_aa = 0.0_wp; p2_ab = 0.0_wp; p2_bb = 0.0_wp
            p3_aaa = 0.0_wp; p3_aab = 0.0_wp
            p3_abb = 0.0_wp; p3_bbb = 0.0_wp

            call cfc_pair_term_eval(d_a, d_b, R_a, R_b, self%a2, self%c, &
                                    max_deriv, pd0_p, p1_a, p1_b, p2_aa, p2_ab, p2_bb, &
                                    p3_aaa, p3_aab, p3_abb, p3_bbb)

            self%PD0 = self%PD0 + pd0_p

            if (max_deriv >= 1) then
               self%PD1 = self%PD1 + p1_a + p1_b
               self%Q1(:, i) = self%Q1(:, i) + p1_a
               self%Q1(:, j) = self%Q1(:, j) + p1_b
            end if

            if (max_deriv >= 2) then
               do jj = 1, ndim
                  do ii = 1, ndim
                     self%PD2(ii, jj) = self%PD2(ii, jj) &
                                        + p2_aa(ii, jj) + p2_bb(ii, jj) &
                                        + p2_ab(ii, jj) + p2_ab(jj, ii)
                     self%Q2(ii, jj, i) = self%Q2(ii, jj, i) &
                                          + p2_aa(ii, jj) + p2_ab(ii, jj)
                     self%Q2(ii, jj, j) = self%Q2(ii, jj, j) &
                                          + p2_bb(ii, jj) + p2_ab(jj, ii)
                  end do
               end do
            end if

            if (max_deriv >= 3) then
               do kk = 1, ndim
                  do jj = 1, ndim
                     do ii = 1, ndim
                        self%PD3(ii, jj, kk) = self%PD3(ii, jj, kk) &
                                               + p3_aaa(ii, jj, kk) + p3_bbb(ii, jj, kk) &
                                               + p3_aab(ii, jj, kk) + p3_aab(ii, kk, jj) &
                                               + p3_aab(jj, kk, ii) &
                                               + p3_abb(ii, jj, kk) + p3_abb(jj, ii, kk) &
                                               + p3_abb(kk, ii, jj)
                        self%Q3(ii, jj, kk, i) = self%Q3(ii, jj, kk, i) &
                                                 + p3_aaa(ii, jj, kk) &
                                                 + p3_aab(ii, jj, kk) + p3_aab(ii, kk, jj) &
                                                 + p3_abb(ii, jj, kk)
                        self%Q3(ii, jj, kk, j) = self%Q3(ii, jj, kk, j) &
                                                 + p3_aab(jj, kk, ii) &
                                                 + p3_abb(jj, ii, kk) + p3_abb(kk, ii, jj) &
                                                 + p3_bbb(ii, jj, kk)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end subroutine cfc_compute_caches

   !* ================================================================================= *!
   !*                          Screened derivative methods                              *!
   !* ================================================================================= *!

   !> Screened value-only LSF: returns `-log PD` at the cached point
   subroutine lsf_f0_screened(self, val)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      real(wp), intent(out) :: val

      val = 0.0_wp
      if (self%ssd_system%n_active == 0) return
      if (self%PD0 <= 0.0_wp) return

      val = -log(self%PD0)
   end subroutine lsf_f0_screened

   !> Screened value/gradient/Hessian w.r.t. spatial coordinates
   subroutine lsf_f012_r_screened(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      real(wp) :: pd1_in(ndim), pd2_in(ndim, ndim), pd3_dummy(ndim, ndim, ndim)
      real(wp) :: lpd0, lpd1(ndim), lpd2(ndim, ndim), lpd3_dummy(ndim, ndim, ndim)
      integer  :: max_lift

      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%ssd_system%n_active == 0) return
      if (self%PD0 <= 0.0_wp) return

      max_lift = 0
      pd1_in = 0.0_wp
      pd2_in = 0.0_wp
      pd3_dummy = 0.0_wp
      if (allocated(self%PD1)) then
         pd1_in = self%PD1
         max_lift = 1
      end if
      if (allocated(self%PD2)) then
         pd2_in = self%PD2
         max_lift = 2
      end if

      call cfc_log_lift(self%PD0, pd1_in, pd2_in, pd3_dummy, max_lift, &
                        lpd0, lpd1, lpd2, lpd3_dummy)

      if (present(lsf0)) lsf0 = -lpd0
      if (present(lsf1_r)) lsf1_r = -lpd1
      if (present(lsf2_rr)) lsf2_rr = -lpd2
   end subroutine lsf_f012_r_screened

   !> Screened third spatial derivative (plus optional lower orders)
   subroutine lsf_f3_rrr_screened(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), allocatable, intent(out) :: lsf3_rrr(:, :, :)

      real(wp) :: pd1_in(ndim), pd2_in(ndim, ndim), pd3_in(ndim, ndim, ndim)
      real(wp) :: lpd0, lpd1(ndim), lpd2(ndim, ndim), lpd3(ndim, ndim, ndim)
      integer  :: max_lift

      allocate (lsf3_rrr(ndim, ndim, ndim))
      lsf3_rrr = 0.0_wp
      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%ssd_system%n_active == 0) return
      if (self%PD0 <= 0.0_wp) return

      max_lift = 0
      pd1_in = 0.0_wp
      pd2_in = 0.0_wp
      pd3_in = 0.0_wp
      if (allocated(self%PD1)) then
         pd1_in = self%PD1
         max_lift = 1
      end if
      if (allocated(self%PD2)) then
         pd2_in = self%PD2
         max_lift = 2
      end if
      if (allocated(self%PD3)) then
         pd3_in = self%PD3
         max_lift = 3
      end if

      call cfc_log_lift(self%PD0, pd1_in, pd2_in, pd3_in, max_lift, &
                        lpd0, lpd1, lpd2, lpd3)

      if (present(lsf0)) lsf0 = -lpd0
      if (present(lsf1_r)) lsf1_r = -lpd1
      if (present(lsf2_rr)) lsf2_rr = -lpd2
      lsf3_rrr = -lpd3
   end subroutine lsf_f3_rrr_screened

   !> Screened mixed third derivative: spatial Hessian w.r.t. nuclei
   subroutine lsf_f3_rr_rA_screened(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: lsf3_rr_rA(:, :, :, :)

      integer  :: i, atom, alpha, j, k
      real(wp) :: invPD, invPD2, invPD3
      real(wp) :: q1a, q2aj, q2ak, q3ajk
      real(wp) :: gj, gk, hjk

      allocate (lsf3_rr_rA(ndim, ndim, ndim, self%ncenters))
      lsf3_rr_rA = 0.0_wp
      if (present(lsf1_rA)) lsf1_rA = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA = 0.0_wp
      if (self%ssd_system%n_active == 0) return
      if (self%PD0 <= 0.0_wp) return

      invPD = 1.0_wp/self%PD0
      invPD2 = invPD*invPD
      invPD3 = invPD2*invPD

      ! The chain rule `d/d(R_A,alpha) = -d/d(d_A,alpha)` already contributes one minus to each term;
      ! the outer `LSF = -log PD` contributes a second minus; thus the assembly is
      ! +Q3/PD - (...) /PD^2 + 2(...)/PD^3 with no extra outer sign

      ! lsf1_rA: requires Q1 only
      if (present(lsf1_rA) .and. allocated(self%Q1)) then
         do i = 1, self%ssd_system%n_active
            atom = self%ssd_system%atom_indices(i)
            do alpha = 1, ndim
               lsf1_rA(alpha, atom) = invPD*self%Q1(alpha, i)
            end do
         end do
      end if

      ! lsf2_r_rA: requires Q1 and Q2
      if (present(lsf2_r_rA) .and. allocated(self%Q1) .and. allocated(self%Q2)) then
         do i = 1, self%ssd_system%n_active
            atom = self%ssd_system%atom_indices(i)
            do alpha = 1, ndim
               q1a = self%Q1(alpha, i)
               do j = 1, ndim
                  q2aj = self%Q2(alpha, j, i)
                  gj = self%PD1(j)
                  lsf2_r_rA(j, alpha, atom) = &
                     invPD*q2aj - invPD2*gj*q1a
               end do
            end do
         end do
      end if

      ! lsf3_rr_rA: requires Q1, Q2, Q3
      if (.not. (allocated(self%Q1) .and. allocated(self%Q2) &
                 .and. allocated(self%Q3))) return

      do i = 1, self%ssd_system%n_active
         atom = self%ssd_system%atom_indices(i)
         do alpha = 1, ndim
            q1a = self%Q1(alpha, i)
            do k = 1, ndim
               gk = self%PD1(k)
               q2ak = self%Q2(alpha, k, i)
               do j = 1, ndim
                  gj = self%PD1(j)
                  hjk = self%PD2(j, k)
                  q2aj = self%Q2(alpha, j, i)
                  q3ajk = self%Q3(alpha, j, k, i)
                  lsf3_rr_rA(j, k, alpha, atom) = &
                     invPD*q3ajk &
                     - invPD2*(hjk*q1a + q2aj*gk + q2ak*gj) &
                     + 2.0_wp*invPD3*gj*gk*q1a
               end do
            end do
         end do
      end do
   end subroutine lsf_f3_rr_rA_screened

   !* ================================================================================= *!
   !*                                 Neighbour cutoff                                  *!
   !* ================================================================================= *!

   !> Conservative shell thickness from the CFC screening threshold
   !>
   !> Atomic term:  exp(a1 (s-1)) <= thr  with s = 1 + delta/R, a1 < 0
   !>     => delta_atomic = -log(thr) / |a1| * R
   !> Pair term:    upper bound c * 2^m * exp(a2 (s-1)) <= thr (taking
   !>     the partner at s_b = 1 and (1 - vec(s_a).vec(s_b))^m <= 2^m)
   !>     => delta_pair  = (log(thr) - log(c * 2^m)) / a2 * R
   !>                    = R * log(thr / (c * 2^m)) / a2
   !>
   !> If `threshold <= 0` or any exponent is non-negative the routine
   !> returns `huge(0.0_wp)` (mirroring the SSD guard in svdw_ssd.f90).
   pure function lsf_neighbor_cutoff(self, radius) result(cutoff_distance)
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      real(wp), intent(in) :: radius
      real(wp) :: cutoff_distance

      real(wp) :: thr, log_thr, delta_atomic, delta_pair, cap
      real(wp) :: two_pow_m

      thr = self%screening_threshold
      if (thr <= 0.0_wp .or. self%a1 >= 0.0_wp .or. self%a2 >= 0.0_wp) then
         cutoff_distance = huge(0.0_wp)
         return
      end if

      log_thr = log(thr)
      delta_atomic = -log_thr/abs(self%a1)*radius

      two_pow_m = 2.0_wp**self%m
      cap = self%c*two_pow_m
      delta_pair = radius*(log_thr - log(cap))/self%a2

      cutoff_distance = max(delta_atomic, delta_pair)
   end function lsf_neighbor_cutoff

   ! ================================================================================= *!
   !                                    Finalizer                                      *!
   ! ================================================================================= *!

   !> Release allocatable components of a CFC LSF instance
   subroutine finalize_lsf_cfc(self)
      type(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self

      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%all_indices)) deallocate (self%all_indices)
      if (allocated(self%PD1)) deallocate (self%PD1)
      if (allocated(self%PD2)) deallocate (self%PD2)
      if (allocated(self%PD3)) deallocate (self%PD3)
      if (allocated(self%Q1)) deallocate (self%Q1)
      if (allocated(self%Q2)) deallocate (self%Q2)
      if (allocated(self%Q3)) deallocate (self%Q3)

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
   end subroutine finalize_lsf_cfc

end module moist_cavity_drop_lsf_cfc
