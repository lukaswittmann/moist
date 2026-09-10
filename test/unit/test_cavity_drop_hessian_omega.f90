!> The DROP surface Hessian with a moving surface adjoint
!>
!> `get_hessian` and `get_surface_hessian` differentiate `J^T omega` at frozen
!> raw adjoints. A model's adjoints are not frozen: they are functions of the
!> surface observables, and the model hands their response along the nuclear
!> directions back to the cavity through a [[surface_adjoint_response_type]].
!> This suite drives that channel with the simplest functional that has a
!> nonzero response, the enclosed volume
!>
!>     V = sum_i a_i (r_i . n_i) / 3,
!>
!> whose adjoint `omega = (dV/da, dV/dr, dV/dn)` moves with every one of the
!> three observables it is built from. The reference is the shipped
!> `get_surface_gradient` with the adjoint **rebuilt at every displaced
!> geometry**, central-differenced, so it is by definition the derivative of
!> the gradient a model would actually return.
!>
!> ## What the assertions separate
!>
!>   * `volume_hessian_fd_*`: the analytic block with the response object
!>     against that reference, at two steps, absolute *or* relative.
!>   * `frozen_misses_response_*`: the same block **without** the response
!>     object against the same reference. It must fail by a wide margin -- the
!>     frozen-adjoint Hessian is a different object -- so a response channel
!>     that silently contributes nothing is caught here, not by the check
!>     above passing for a trivial reason.
!>   * `hvp_matches_dense_*`: the per-direction entry with the response object
!>     reproduces the dense block, to the bit on unit directions and to a few
!>     ulp on general ones.
!>   * `curvature_response_refused`: a response carrying curvature weights on
!>     top of adjoints that carry none is refused rather than served the zero
!>     the prologue would give it.
!>
!> The frozen block is symmetric on its own -- it is the Hessian of a linear
!> functional of the observables -- and so is the response term, so symmetry
!> cannot tell the two apart; it is checked as a sanity bound only.
!>
!> ## Step and tolerance
!>
!> The steps and the bound are those of the fixed-adjoint end-to-end suite,
!> whose reference is the same differenced gradient on the same fixture. The
!> worst measured deviations are on `OMEGA_TOL`.
module test_cavity_drop_hessian_omega
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_type, only: cavity_type, surface_adjoint_response_type
   use moist_cavity_drop, only: cavity_type_drop, drop_hvp_per_dir_max
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type
   use moist_context, only: moist_context_type
   use test_helpers, only: drop_fixture_geometry, build_drop_test_cavity, &
                           LSF_SVDW, LSF_CFC, FIX_PLAIN

   implicit none(type, external)
   private

   public :: collect_cavity_drop_hessian_omega

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Central-difference steps of the reference; two, so that a value that
   !> agrees at one step only still fails
   real(wp), parameter :: FD_STEPS(*) = [3.0E-4_wp, 2.5E-4_wp]

   !> Agreement bound, absolute *and* relative: a component fails only when
   !> it misses both
   !>
   !> Measured with the response object, worst over the block and both steps:
   !> `9.4e-11` (SvdW) and `1.3e-10` (CFC), the reference's own round-off
   !> floor, as on the fixed-adjoint suite.
   real(wp), parameter :: OMEGA_TOL = 3.0E-10_wp

   !> Symmetry bound of the analytic block, relative to its largest entry
   real(wp), parameter :: SYM_TOL = 1.0E-11_wp

   !> Below this the reference carries no information
   real(wp), parameter :: VACUITY_THR = 1.0E-4_wp

   !> How far the frozen block must miss the reference, in units of `OMEGA_TOL`
   !>
   !> Measured: the frozen block misses by `11.3` (SvdW) and `10.0` (CFC)
   !> absolute on a block whose largest entry is `35` / `45`, eleven decades
   !> over the bound.
   real(wp), parameter :: TEETH_FACTOR = 1.0E2_wp

   !> Agreement of the per-direction entry with the contracted dense block
   real(wp), parameter :: HVP_UNIT_TOL = 0.0_wp
   real(wp), parameter :: HVP_GEN_TOL = 1.0E-13_wp

   !> Adjoint response of the enclosed volume
   !>
   !> `omega = p (r.n/3, a n/3, a r/3)` on `(a, r, n)`, so along a surface
   !> tangent `(d_a, d_r, d_n)` its response is
   !>
   !>     d(w_a)   = p (d_r . n + r . d_n) / 3
   !>     d(w_xyz) = p (d_a n + a d_n) / 3
   !>     d(w_n)   = p (d_a r + a d_r) / 3
   type, extends(surface_adjoint_response_type) :: volume_response_type
      !> Scale of the functional
      real(wp) :: pressure = 1.0_wp
      !> Whether to add a curvature weight the primal adjoint does not carry
      logical :: inject_curvature = .false.
   contains
      procedure :: apply => volume_response_apply
   end type volume_response_type

contains

   !> Collect the suite
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_cavity_drop_hessian_omega(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("volume_hessian_fd_svdw", test_volume_svdw), &
                  new_unittest("volume_hessian_fd_cfc", test_volume_cfc), &
                  new_unittest("frozen_misses_response_svdw", test_teeth_svdw), &
                  new_unittest("frozen_misses_response_cfc", test_teeth_cfc), &
                  new_unittest("hvp_matches_dense_svdw", test_hvp_svdw), &
                  new_unittest("hvp_matches_dense_cfc", test_hvp_cfc), &
                  new_unittest("curvature_response_refused", test_curvature_refused) &
                  ]
   end subroutine collect_cavity_drop_hessian_omega

   !* ================================================================================= *!
   !*                            The response of the volume                            *!
   !* ================================================================================= *!

   !> Adjoint response of the volume functional along a block of directions
   subroutine volume_response_apply(self, cavity, dirs, tangent, dacc, hvp_direct, error)
      !> Response object
      class(volume_response_type), intent(inout) :: self
      !> Cavity the tangent was taken on
      class(cavity_type), intent(in) :: cavity
      !> Nuclear directions of the block
      real(wp), intent(in) :: dirs(:, :, :)
      !> Surface tangent of the block
      type(cavity_surface_tangent_type), intent(in) :: tangent
      !> Surface-adjoint response per direction
      type(cavity_surface_adjoint_type), intent(inout) :: dacc(:)
      !> Non-surface Hessian columns of the block
      real(wp), intent(inout) :: hvp_direct(:, :, :)
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: error

      real(wp), allocatable :: dw_a(:), dw_xyz(:, :), dw_n(:, :), dw_k1(:)
      real(wp) :: third
      integer :: ngrid, ndir, idir, igrid

      ngrid = cavity%ngrid
      ndir = size(dirs, 3)
      third = self%pressure/3.0_wp
      allocate (dw_a(ngrid), dw_xyz(ndim, ngrid), dw_n(ndim, ngrid))

      do idir = 1, ndir
         do igrid = 1, ngrid
            dw_a(igrid) = third*(dot_product(tangent%d_xyz(:, igrid, idir), &
                                             cavity%normal0(:, igrid)) &
                                 + dot_product(cavity%xyz(:, igrid), tangent%d_n(:, igrid, idir)))
            dw_xyz(:, igrid) = third*(tangent%d_a(igrid, idir)*cavity%normal0(:, igrid) &
                                      + cavity%a(igrid)*tangent%d_n(:, igrid, idir))
            dw_n(:, igrid) = third*(tangent%d_a(igrid, idir)*cavity%xyz(:, igrid) &
                                    + cavity%a(igrid)*tangent%d_xyz(:, igrid, idir))
         end do
         call dacc(idir)%add_surface_weights(error, w_a=dw_a, w_xyz=dw_xyz, w_n=dw_n)
         if (allocated(error)) return
         if (self%inject_curvature) then
            allocate (dw_k1(ngrid), source=1.0_wp)
            call dacc(idir)%add_surface_weights(error, w_k1=dw_k1)
            if (allocated(error)) return
            deallocate (dw_k1)
         end if
      end do

   end subroutine volume_response_apply

   !> Surface adjoint of the volume on the current geometry
   !>
   !> @param[in]  cavity   Cavity supplying the surface
   !> @param[in]  pressure Scale of the functional
   !> @param[out] acc      Surface-adjoint accumulator
   !> @param[out] error    Error handle
   subroutine volume_adjoint(cavity, pressure, acc, error)
      !> Cavity supplying the surface
      type(cavity_type_drop), intent(in) :: cavity
      !> Scale of the functional
      real(wp), intent(in) :: pressure
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: add_error
      real(wp), allocatable :: w_a(:), w_xyz(:, :), w_n(:, :)
      integer :: ngrid, igrid

      ngrid = cavity%ngrid
      allocate (w_a(ngrid), w_xyz(ndim, ngrid), w_n(ndim, ngrid))
      do igrid = 1, ngrid
         w_a(igrid) = pressure/3.0_wp*dot_product(cavity%xyz(:, igrid), cavity%normal0(:, igrid))
         w_xyz(:, igrid) = pressure*cavity%a(igrid)/3.0_wp*cavity%normal0(:, igrid)
         w_n(:, igrid) = pressure*cavity%a(igrid)/3.0_wp*cavity%xyz(:, igrid)
      end do

      call acc%init(ngrid)
      call acc%add_surface_weights(add_error, w_a=w_a, w_xyz=w_xyz, w_n=w_n)
      if (allocated(add_error)) then
         call test_failed(error, "failed to seed the volume adjoint: "//add_error%message)
         return
      end if

   end subroutine volume_adjoint

   !* ================================================================================= *!
   !*                          Finite-difference reference                             *!
   !* ================================================================================= *!

   !> SvdW: analytic block with the response against the differenced gradient
   subroutine test_volume_svdw(error)
      type(error_type), allocatable, intent(out) :: error

      call run_volume_fd(LSF_SVDW, .true., "svdw", error)
   end subroutine test_volume_svdw

   !> CFC: analytic block with the response against the differenced gradient
   subroutine test_volume_cfc(error)
      type(error_type), allocatable, intent(out) :: error

      call run_volume_fd(LSF_CFC, .true., "cfc", error)
   end subroutine test_volume_cfc

   !> SvdW: the frozen block must miss the differenced gradient
   subroutine test_teeth_svdw(error)
      type(error_type), allocatable, intent(out) :: error

      call run_volume_fd(LSF_SVDW, .false., "svdw", error)
   end subroutine test_teeth_svdw

   !> CFC: the frozen block must miss the differenced gradient
   subroutine test_teeth_cfc(error)
      type(error_type), allocatable, intent(out) :: error

      call run_volume_fd(LSF_CFC, .false., "cfc", error)
   end subroutine test_teeth_cfc

   !> Compare the analytic volume Hessian, with or without the response, to
   !> the differenced volume gradient
   !>
   !> @param[in]  lsf_kind      Level-set model
   !> @param[in]  with_response Whether the response object is handed in
   !> @param[in]  label         Human-readable case description
   !> @param[out] error         Error handle
   subroutine run_volume_fd(lsf_kind, with_response, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Whether the response object is handed in
      logical, intent(in) :: with_response
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_surface_adjoint_type) :: acc
      type(volume_response_type) :: omega
      type(mctc_error), allocatable :: cav_error

      real(wp), allocatable :: hess(:, :, :, :), ref(:, :, :, :)
      real(wp) :: worst, worst_ref, scale, sym
      integer :: istep, nsph, i1, i2, i3, i4
      logical :: grid_ok

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, lsf_kind, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      nsph = cavity%nsph

      call volume_adjoint(cavity, 1.0_wp, acc, error)
      if (allocated(error)) return

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      if (with_response) then
         call cavity%get_hessian(acc, hess, cav_error, omega_v=omega)
      else
         call cavity%get_hessian(acc, hess, cav_error)
      end if
      if (allocated(cav_error)) then
         call test_failed(error, "get_hessian failed ("//label//"): "//cav_error%message)
         return
      end if

      scale = maxval(abs(hess))
      if (scale <= VACUITY_THR) then
         call test_failed(error, "analytic volume Hessian is vacuous ("//label//")")
         return
      end if

      ! Sanity: symmetric with or without the response, both pieces are on
      ! their own. Not translationally invariant: `sum_i a_i n_i` is the vector
      ! area of the discretised surface, exactly zero for a closed surface and
      ! only small for a switched, filtered grid, and its nuclear gradient is
      ! what the column sum of this block measures; the reference carries the
      ! same sum, and the comparison below covers it.
      sym = 0.0_wp
      do i4 = 1, nsph
         do i3 = 1, ndim
            do i2 = 1, nsph
               do i1 = 1, ndim
                  sym = max(sym, abs(hess(i1, i2, i3, i4) - hess(i3, i4, i1, i2)))
               end do
            end do
         end do
      end do
      if (sym > SYM_TOL*scale) then
         call test_failed(error, "analytic volume Hessian is asymmetric ("//label// &
                          "): defect "//to_string(sym)//" against "//to_string(scale))
         return
      end if

      do istep = 1, size(FD_STEPS)
         call numerical_volume_hessian(mol, lsf_kind, FD_STEPS(istep), ref, grid_ok, error)
         if (allocated(error)) return
         if (.not. grid_ok) then
            call test_failed(error, "grid changed under the step "//to_string(FD_STEPS(istep))// &
                             " ("//label//")")
            return
         end if
         if (maxval(abs(ref)) <= VACUITY_THR) then
            call test_failed(error, "reference volume Hessian is vacuous ("//label//")")
            return
         end if

         call worst_deviation(hess, ref, worst, worst_ref)
         if (with_response) then
            if (worst > OMEGA_TOL .and. worst > OMEGA_TOL*worst_ref) then
               call test_failed(error, "volume Hessian mismatch ("//label//", h = "// &
                                to_string(FD_STEPS(istep))//"): worst deviation "// &
                                to_string(worst)//" against reference "//to_string(worst_ref))
               return
            end if
         else
            if (worst <= TEETH_FACTOR*OMEGA_TOL) then
               call test_failed(error, "the frozen-adjoint Hessian reproduces the moving-"// &
                                "adjoint reference ("//label//", h = "// &
                                to_string(FD_STEPS(istep))//"): worst deviation "// &
                                to_string(worst)//"; the response channel has no teeth")
               return
            end if
         end if
      end do

   end subroutine run_volume_fd

   !> Five-point central difference of the reverse-mode volume gradient
   !>
   !> The adjoint is rebuilt on every displaced cavity, so the difference is
   !> that of the gradient a model would return, moving adjoints included.
   !>
   !> @param[in]  mol      Base structure
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  step     Central-difference step
   !> @param[out] hess     Differenced gradient `(3, nsph, 3, nsph)`
   !> @param[out] grid_ok  Whether every displaced grid matched the base one
   !> @param[out] error    Error handle
   subroutine numerical_volume_hessian(mol, lsf_kind, step, hess, grid_ok, error)
      !> Base structure
      type(structure_type), intent(in) :: mol
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Central-difference step
      real(wp), intent(in) :: step
      !> Differenced gradient
      real(wp), allocatable, intent(out) :: hess(:, :, :, :)
      !> Whether every displaced grid matched the base one
      logical, intent(out) :: grid_ok
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(cavity_type_drop), allocatable :: ref_cav, cavity
      type(moist_context_type), target :: ref_ctx, ctx
      type(structure_type) :: mol_disp
      type(cavity_surface_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: grad(:, :)
      integer :: nsph, batom, baxis, ioff

      grid_ok = .true.
      call build_drop_test_cavity(ref_cav, ref_ctx, mol, FIX_PLAIN, lsf_kind, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      nsph = ref_cav%nsph
      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      allocate (grad(ndim, nsph))

      do batom = 1, nsph
         do baxis = 1, ndim
            do ioff = 1, size(OFFSET)
               mol_disp = mol
               mol_disp%xyz(baxis, batom) = mol%xyz(baxis, batom) + real(OFFSET(ioff), wp)*step
               call build_drop_test_cavity(cavity, ctx, mol_disp, FIX_PLAIN, lsf_kind, error, &
                                           want_fine=.true.)
               if (allocated(error)) return
               if (.not. same_grid(ref_cav, cavity)) then
                  grid_ok = .false.
                  return
               end if

               call volume_adjoint(cavity, 1.0_wp, acc, error)
               if (allocated(error)) return
               grad = 0.0_wp
               call cavity%get_surface_gradient(acc, grad, cav_error)
               if (allocated(cav_error)) then
                  call test_failed(error, "surface gradient failed: "//cav_error%message)
                  return
               end if
               hess(:, :, baxis, batom) = hess(:, :, baxis, batom) + COEFF(ioff)*grad/step
               deallocate (cavity)
            end do
         end do
      end do

   end subroutine numerical_volume_hessian

   !> Worst deviation between two blocks, and the reference entry it sits on
   subroutine worst_deviation(hess, ref, worst, worst_ref)
      real(wp), intent(in) :: hess(:, :, :, :), ref(:, :, :, :)
      real(wp), intent(out) :: worst, worst_ref

      real(wp) :: diff
      integer :: i1, i2, i3, i4

      worst = 0.0_wp
      worst_ref = 0.0_wp
      do i4 = 1, size(hess, 4)
         do i3 = 1, size(hess, 3)
            do i2 = 1, size(hess, 2)
               do i1 = 1, size(hess, 1)
                  diff = abs(hess(i1, i2, i3, i4) - ref(i1, i2, i3, i4))
                  ! A component fails only when it misses the absolute and the
                  ! relative bound, so the worst one is the one that misses
                  ! the relative bound by most among those over the absolute
                  if (diff > OMEGA_TOL .and. diff > OMEGA_TOL*abs(ref(i1, i2, i3, i4))) then
                     if (diff > worst) then
                        worst = diff
                        worst_ref = abs(ref(i1, i2, i3, i4))
                     end if
                  end if
               end do
            end do
         end do
      end do
   end subroutine worst_deviation

   !> Do two cavities carry the very same grid?
   function same_grid(ref, cav) result(same)
      type(cavity_type_drop), intent(in) :: ref, cav
      logical :: same

      integer :: igrid

      same = .false.
      if (cav%ngrid /= ref%ngrid) return
      do igrid = 1, ref%ngrid
         if (cav%numbering(igrid) /= ref%numbering(igrid)) return
         if (cav%owner(igrid) /= ref%owner(igrid)) return
      end do
      if (allocated(ref%branch_count) .and. allocated(cav%branch_count)) then
         do igrid = 1, ref%ngrid
            if (cav%branch_count(igrid) /= ref%branch_count(igrid)) return
         end do
      end if
      same = .true.
   end function same_grid

   !* ================================================================================= *!
   !*                     Hessian-vector products with the response                    *!
   !* ================================================================================= *!

   !> SvdW: per-direction products with the response reproduce the dense block
   subroutine test_hvp_svdw(error)
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(LSF_SVDW, "svdw", error)
   end subroutine test_hvp_svdw

   !> CFC: per-direction products with the response reproduce the dense block
   subroutine test_hvp_cfc(error)
      type(error_type), allocatable, intent(out) :: error

      call run_hvp(LSF_CFC, "cfc", error)
   end subroutine test_hvp_cfc

   !> Unit directions to the bit, general directions to a few ulp
   !>
   !> @param[in]  lsf_kind Level-set model
   !> @param[in]  label    Human-readable case description
   !> @param[out] error    Error handle
   subroutine run_hvp(lsf_kind, label, error)
      !> Level-set model
      integer, intent(in) :: lsf_kind
      !> Case description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_surface_adjoint_type) :: acc
      type(volume_response_type) :: omega
      type(mctc_error), allocatable :: cav_error

      real(wp), allocatable :: hess(:, :, :, :), dirs(:, :, :), hvp(:, :, :), expect(:, :, :)
      integer :: nsph, ndir, idir, iatom, iaxis
      real(wp) :: scale

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, lsf_kind, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      nsph = cavity%nsph

      call volume_adjoint(cavity, 1.0_wp, acc, error)
      if (allocated(error)) return

      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      call cavity%get_hessian(acc, hess, cav_error, omega_v=omega)
      if (allocated(cav_error)) then
         call test_failed(error, "get_hessian failed ("//label//"): "//cav_error%message)
         return
      end if
      scale = maxval(abs(hess))

      !* --------------------- Unit directions: the dense columns --------------------- *!
      ndir = ndim*nsph
      allocate (dirs(ndim, nsph, ndir), source=0.0_wp)
      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do
      allocate (hvp(ndim, nsph, ndir), source=0.0_wp)
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error, omega_v=omega)
      if (allocated(cav_error)) then
         call test_failed(error, "get_surface_hessian failed ("//label//"): "// &
                          cav_error%message)
         return
      end if
      do iatom = 1, nsph
         do iaxis = 1, ndim
            idir = ndim*(iatom - 1) + iaxis
            if (maxval(abs(hvp(:, :, idir) - hess(:, :, iaxis, iatom))) > HVP_UNIT_TOL) then
               call test_failed(error, "unit-direction product differs from the dense"// &
                                " column ("//label//"): "// &
                                to_string(maxval(abs(hvp(:, :, idir) - hess(:, :, iaxis, iatom)))))
               return
            end if
         end do
      end do
      deallocate (dirs, hvp)

      !* ----------------- General directions: the contracted block ------------------- *!
      ndir = min(drop_hvp_per_dir_max, ndim*nsph - 1)
      call build_directions(nsph, ndir, dirs)
      allocate (hvp(ndim, nsph, ndir), source=0.0_wp)
      allocate (expect(ndim, nsph, ndir), source=0.0_wp)
      call cavity%get_surface_hessian(acc, dirs, hvp, cav_error, omega_v=omega)
      if (allocated(cav_error)) then
         call test_failed(error, "get_surface_hessian failed on general directions ("// &
                          label//"): "//cav_error%message)
         return
      end if
      do idir = 1, ndir
         do iatom = 1, nsph
            do iaxis = 1, ndim
               expect(:, :, idir) = expect(:, :, idir) &
                                    + hess(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
            end do
         end do
      end do
      if (maxval(abs(hvp - expect)) > HVP_GEN_TOL*scale) then
         call test_failed(error, "general-direction product differs from the contracted"// &
                          " block ("//label//"): "//to_string(maxval(abs(hvp - expect)))// &
                          " against "//to_string(scale))
         return
      end if

   end subroutine run_hvp

   !> Deterministic, dense, non-unit directions
   subroutine build_directions(nsph, ndir, dirs)
      integer, intent(in) :: nsph, ndir
      real(wp), allocatable, intent(out) :: dirs(:, :, :)

      integer :: idir, iatom, iaxis

      allocate (dirs(ndim, nsph, ndir))
      do idir = 1, ndir
         do iatom = 1, nsph
            do iaxis = 1, ndim
               dirs(iaxis, iatom, idir) = 0.5_wp + 0.4_wp*sin(0.37_wp*real(iaxis, wp) &
                                                                + 0.61_wp*real(iatom, wp) &
                                                                + 1.1_wp*real(idir, wp))
            end do
         end do
      end do
   end subroutine build_directions

   !* ================================================================================= *!
   !*                                     Guards                                       *!
   !* ================================================================================= *!

   !> A curvature response on a curvature-free adjoint is refused
   subroutine test_curvature_refused(error)
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(cavity_surface_adjoint_type) :: acc
      type(volume_response_type) :: omega
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: hess(:, :, :, :)

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      call volume_adjoint(cavity, 1.0_wp, acc, error)
      if (allocated(error)) return

      omega%inject_curvature = .true.
      allocate (hess(ndim, cavity%nsph, ndim, cavity%nsph), source=0.0_wp)
      call cavity%get_hessian(acc, hess, cav_error, omega_v=omega)
      if (.not. allocated(cav_error)) then
         call test_failed(error, "a curvature response on curvature-free adjoints was accepted")
         return
      end if
      if (index(cav_error%message, "curvature") == 0) then
         call test_failed(error, "refusal does not name the curvature channel: "// &
                          cav_error%message)
         return
      end if
      if (maxval(abs(hess)) /= 0.0_wp) then
         call test_failed(error, "the accumulator was touched on failure")
         return
      end if

   end subroutine test_curvature_refused

end module test_cavity_drop_hessian_omega
