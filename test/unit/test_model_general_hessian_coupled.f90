!> The second-order host exchange of the general model: the response tangent
!>
!> `get_hvp` and `get_hessian` accept a [[response_tangent_type]] that
!> receives, per direction, what a host completes its own Hessian columns and
!> Fock-matrix tangents from: the tangent of the surface charges, the motion
!> of the surface points, and the gradient-path level-set weights with their
!> tangents. This suite pins that container against references that exist
!> independently of it, on the DROP fixture with PV and CPCM (point-charge
!> source) components:
!>
!>   * `*_columns_unchanged`: asking for the response tangent does not move
!>     the Hessian-vector products. On general directions both runs take the
!>     per-direction fixed channel and must agree to the bit; on the Cartesian
!>     basis the forced per-direction form is compared with the rank-4 dense
!>     block to a few ulp;
!>   * `xyz_matches_tangent`: the position tangent is the cavity's own surface
!>     tangent;
!>   * `cpcm_charge_tangent_fd`: the charge tangent against the differenced
!>     surface charges over rebuilt models along general directions;
!>   * `*_weight_tangent_fd`: the level-set weight tangents against the
!>     differenced base weights of the same container; and, for PV, the base
!>     weights against the potential-path level-set response;
!>   * `guards`: an uninitialised or mis-sized container is refused, and a
!>     container without level-set channels leaves the rank-4 form in charge.
module test_model_general_hessian_coupled
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: new_unittest, unittest_type, error_type, to_string, test_failed
   use moist_channels, only: coupling_type, response_type, response_tangent_type
   use moist_cavity_drop, only: cavity_type_drop, drop_hvp_per_dir_max
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type
   use moist_model_component_pcm_type, only: solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, new_component_cpcm
   use moist_model_components, only: solvation_model_component_pv, new_component_pv
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_context, only: moist_context_type
   use test_helpers, only: drop_fixture_geometry, build_drop_test_cavity, &
                           make_charge_coupling, LSF_SVDW, FIX_PLAIN

   implicit none(type, external)
   private

   public :: collect_model_general_hessian_coupled

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Component selectors
   integer, parameter :: MODEL_PV = 1, MODEL_CPCM = 2

   !> Pressure, dielectric constant and point charges of the components
   real(wp), parameter :: PRESSURE = 2.5_wp
   real(wp), parameter :: EPSILON = 32.0_wp
   real(wp), parameter :: QAT(3) = [1.2_wp, -1.8_wp, 0.6_wp]

   !> Central-difference steps of the references; two, so that a value that
   !> agrees at one step only still fails
   real(wp), parameter :: FD_STEPS(*) = [3.0E-4_wp, 2.5E-4_wp]

   !> Agreement bounds of the differenced references, absolute *and* relative
   real(wp), parameter :: FD_TOL = 1.0E-8_wp

   !> Agreement of the forced per-direction form with the rank-4 dense block,
   !> and of the position tangent with the materialised surface tangent,
   !> relative to the largest entry
   real(wp), parameter :: ULP_TOL = 1.0E-13_wp

   !> Below this a reference carries no information
   real(wp), parameter :: VACUITY_THR = 1.0E-6_wp

contains

   !> Collect the suite
   subroutine collect_model_general_hessian_coupled(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pv_columns_unchanged", test_pv_columns), &
                  new_unittest("cpcm_columns_unchanged", test_cpcm_columns), &
                  new_unittest("xyz_matches_tangent", test_xyz), &
                  new_unittest("cpcm_charge_tangent_fd", test_charge_fd), &
                  new_unittest("pv_weight_tangent_fd", test_pv_weights), &
                  new_unittest("cpcm_weight_tangent_fd", test_cpcm_weights), &
                  new_unittest("guards", test_guards) &
                  ]
   end subroutine collect_model_general_hessian_coupled

   !* ================================================================================= *!
   !*                                    Fixture                                       *!
   !* ================================================================================= *!

   !> The model of one selector on the DROP fixture, with its coupling
   subroutine build_model(kind, mol, model, ctx, coupling, error)
      integer, intent(in) :: kind
      type(structure_type), intent(in) :: mol
      type(solvation_model_general), intent(out) :: model
      type(moist_context_type), intent(inout), target :: ctx
      type(coupling_type), intent(out) :: coupling
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(solvation_model_component_pv) :: pv
      type(solvation_model_component_cpcm) :: cpcm
      type(mctc_error), allocatable :: err

      call build_drop_test_cavity(cavity, ctx, mol, FIX_PLAIN, LSF_SVDW, error, &
                                  want_fine=.true.)
      if (allocated(error)) return
      call new_model_general(model, cavity, ctx, err)
      if (.not. allocated(err)) then
         select case (kind)
         case (MODEL_PV)
            call new_component_pv(pv, PRESSURE)
            call model%add_component(pv, err)
         case (MODEL_CPCM)
            call make_charge_coupling(QAT, coupling)
            call new_component_cpcm(cpcm, ctx, EPSILON, solver=solver_type%cholesky, &
                                    phi_source=potential_source%charges, error=err)
            if (.not. allocated(err)) call model%add_component(cpcm, err)
         case default
            call test_failed(error, "unknown component selector")
            return
         end select
      end if
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "model setup failed: "//err%message)
         return
      end if
   end subroutine build_model

   !> A set of general directions, few enough for the per-direction channel
   subroutine general_directions(nsph, dirs)
      integer, intent(in) :: nsph
      real(wp), allocatable, intent(out) :: dirs(:, :, :)
      integer :: ndir, idir, iatom, iaxis

      ndir = min(drop_hvp_per_dir_max, ndim*nsph - 1)
      allocate (dirs(ndim, nsph, ndir))
      do idir = 1, ndir
         do iatom = 1, nsph
            do iaxis = 1, ndim
               dirs(iaxis, iatom, idir) = 0.5_wp + 0.4_wp*sin(0.37_wp*iaxis + 0.61_wp*iatom &
                                                                + 1.1_wp*idir)
            end do
         end do
      end do
   end subroutine general_directions

   !> The Cartesian unit directions
   subroutine unit_directions(nsph, dirs)
      integer, intent(in) :: nsph
      real(wp), allocatable, intent(out) :: dirs(:, :, :)
      integer :: iatom, iaxis

      allocate (dirs(ndim, nsph, ndim*nsph), source=0.0_wp)
      do iatom = 1, nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do
   end subroutine unit_directions

   !> The persistent grid identity of the model's DROP cavity
   subroutine grid_identity(model, numbering, owner, error)
      type(solvation_model_general), intent(in) :: model
      integer, allocatable, intent(out) :: numbering(:), owner(:)
      type(error_type), allocatable, intent(out) :: error

      select type (cav => model%cavity)
      type is (cavity_type_drop)
         numbering = cav%numbering(1:cav%ngrid)
         owner = cav%owner(1:cav%ngrid)
      class default
         call test_failed(error, "the model cavity is not a DROP cavity")
      end select
   end subroutine grid_identity

   !> The symmetric part of every `(3, 3)` block of a `(3, 3, n)` array
   pure function symmetrised(w) result(s)
      real(wp), intent(in) :: w(:, :, :)
      real(wp) :: s(size(w, 1), size(w, 2), size(w, 3))
      integer :: i

      do i = 1, size(w, 3)
         s(:, :, i) = 0.5_wp*(w(:, :, i) + transpose(w(:, :, i)))
      end do
   end function symmetrised

   !> Worst deviation among the entries that miss both an absolute and a
   !> relative bound, and the reference entry it was measured against
   subroutine worst_deviation(val, ref, tol, worst, worst_ref)
      real(wp), intent(in) :: val(:), ref(:)
      real(wp), intent(in) :: tol
      real(wp), intent(out) :: worst, worst_ref
      real(wp) :: diff
      integer :: i

      worst = 0.0_wp
      worst_ref = 0.0_wp
      do i = 1, size(val)
         diff = abs(val(i) - ref(i))
         if (diff > tol .and. diff > tol*abs(ref(i))) then
            if (diff > worst) then
               worst = diff
               worst_ref = abs(ref(i))
            end if
         end if
      end do
   end subroutine worst_deviation

   !* ================================================================================= *!
   !*                      The columns do not move, the positions match                *!
   !* ================================================================================= *!

   subroutine test_pv_columns(error)
      type(error_type), allocatable, intent(out) :: error
      call run_columns(MODEL_PV, error)
   end subroutine test_pv_columns

   subroutine test_cpcm_columns(error)
      type(error_type), allocatable, intent(out) :: error
      call run_columns(MODEL_CPCM, error)
   end subroutine test_cpcm_columns

   !> Products with and without the response tangent
   subroutine run_columns(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(response_tangent_type) :: rt
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), plain(:, :, :)
      real(wp), allocatable :: hess(:, :, :, :), dense(:, :, :, :)
      integer :: nsph, ngrid

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(kind, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      ngrid = model%cavity%ngrid

      ! General directions: both runs are per direction, to the bit
      call general_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      allocate (plain(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      call model%get_hvp(coupling, dirs, plain, err)
      if (.not. allocated(err)) then
         call rt%init(ngrid, size(dirs, 3), want_lsf=.true.)
         call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      end if
      if (allocated(err)) then
         call test_failed(error, "Hessian-vector products failed: "//err%message)
         return
      end if
      if (maxval(abs(plain)) <= VACUITY_THR) then
         call test_failed(error, "Hessian-vector products are vacuous")
         return
      end if
      if (any(hvp /= plain)) then
         call test_failed(error, "the response tangent moved the per-direction columns by "// &
                          to_string(maxval(abs(hvp - plain))))
         return
      end if
      if (maxval(abs(rt%xyz)) <= VACUITY_THR .or. maxval(abs(rt%dw_value)) <= VACUITY_THR) then
         call test_failed(error, "the response tangent is vacuous")
         return
      end if
      if (kind == MODEL_CPCM .and. maxval(abs(rt%surface_charge)) <= VACUITY_THR) then
         call test_failed(error, "the charge tangent is vacuous")
         return
      end if
      if (kind == MODEL_PV .and. maxval(abs(rt%surface_charge)) /= 0.0_wp) then
         call test_failed(error, "a PV model deposited a charge tangent")
         return
      end if

      ! The Cartesian basis: the forced per-direction dense block against the
      ! rank-4 one
      allocate (hess(ndim, nsph, ndim, nsph), source=0.0_wp)
      allocate (dense(ndim, nsph, ndim, nsph), source=0.0_wp)
      call model%get_hessian(coupling, dense, err)
      if (.not. allocated(err)) then
         call rt%init(ngrid, ndim*nsph, want_lsf=.true.)
         call model%get_hessian(coupling, hess, err, rt=rt)
      end if
      if (allocated(err)) then
         call test_failed(error, "dense Hessian failed: "//err%message)
         return
      end if
      if (maxval(abs(hess - dense)) > ULP_TOL*maxval(abs(dense))) then
         call test_failed(error, "the forced per-direction dense block differs from the"// &
                          " rank-4 one by "//to_string(maxval(abs(hess - dense)))// &
                          " against "//to_string(maxval(abs(dense))))
         return
      end if
   end subroutine run_columns

   !> The position tangent is the cavity's own surface tangent
   subroutine test_xyz(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(response_tangent_type) :: rt
      type(cavity_surface_tangent_type) :: tangent
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :)
      integer :: nsph, ngrid

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(MODEL_PV, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      ngrid = model%cavity%ngrid
      call general_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      call rt%init(ngrid, size(dirs, 3), want_lsf=.false.)
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (.not. allocated(err)) then
         call tangent%init(ngrid, size(dirs, 3), .false.)
         call model%cavity%get_surface_tangent(dirs, tangent, err)
      end if
      if (allocated(err)) then
         call test_failed(error, "tangent failed: "//err%message)
         return
      end if
      if (maxval(abs(tangent%d_xyz)) <= VACUITY_THR) then
         call test_failed(error, "the surface tangent is vacuous")
         return
      end if
      if (maxval(abs(rt%xyz - tangent%d_xyz)) > ULP_TOL*maxval(abs(tangent%d_xyz))) then
         call test_failed(error, "the position tangent differs from the surface tangent by "// &
                          to_string(maxval(abs(rt%xyz - tangent%d_xyz))))
         return
      end if
   end subroutine test_xyz

   !* ================================================================================= *!
   !*                           Finite-difference references                           *!
   !* ================================================================================= *!

   !> The charge tangent against the differenced surface charges
   subroutine test_charge_fd(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(response_tangent_type) :: rt
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), ref(:, :), q(:)
      real(wp) :: worst, worst_ref
      integer :: nsph, ngrid, istep, idir

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(MODEL_CPCM, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      ngrid = model%cavity%ngrid
      call general_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      call rt%init(ngrid, size(dirs, 3), want_lsf=.false.)
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (allocated(err)) then
         call test_failed(error, "Hessian-vector products failed: "//err%message)
         return
      end if

      allocate (ref(ngrid, size(dirs, 3)))
      do istep = 1, size(FD_STEPS)
         do idir = 1, size(dirs, 3)
            call differenced(MODEL_CPCM, mol, dirs(:, :, idir), FD_STEPS(istep), &
                             charges_at, ref(:, idir), error)
            if (allocated(error)) return
         end do
         if (maxval(abs(ref)) <= VACUITY_THR) then
            call test_failed(error, "the differenced charges are vacuous")
            return
         end if
         call worst_deviation(reshape(rt%surface_charge, [size(ref)]), &
                              reshape(ref, [size(ref)]), FD_TOL, worst, worst_ref)
         if (worst > 0.0_wp) then
            call test_failed(error, "charge tangent mismatch (h = "//to_string(FD_STEPS(istep))// &
                             "): worst deviation "//to_string(worst)//" against reference "// &
                             to_string(worst_ref))
            return
         end if
      end do
      if (allocated(q)) deallocate (q)
   end subroutine test_charge_fd

   subroutine test_pv_weights(error)
      type(error_type), allocatable, intent(out) :: error
      call run_weight_fd(MODEL_PV, error)
   end subroutine test_pv_weights

   subroutine test_cpcm_weights(error)
      type(error_type), allocatable, intent(out) :: error
      call run_weight_fd(MODEL_CPCM, error)
   end subroutine test_cpcm_weights

   !> The level-set weight tangents against the differenced base weights
   subroutine run_weight_fd(kind, error)
      integer, intent(in) :: kind
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(response_tangent_type) :: rt
      type(response_type) :: response
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), ref(:, :), val(:, :), base(:)
      real(wp) :: worst, worst_ref
      integer :: nsph, ngrid, istep, idir, nw

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(kind, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      ngrid = model%cavity%ngrid
      call general_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      call rt%init(ngrid, size(dirs, 3), want_lsf=.true.)
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (allocated(err)) then
         call test_failed(error, "Hessian-vector products failed: "//err%message)
         return
      end if

      ! PV: the gradient-path base weights are the potential-path ones
      if (kind == MODEL_PV) then
         call model%get_response(coupling, response, err)
         if (allocated(err)) then
            call test_failed(error, "model response failed: "//err%message)
            return
         end if
         if (.not. allocated(response%lsf%w_value)) then
            call test_failed(error, "the potential path returned no level-set weights")
            return
         end if
         ! The Hessian weight is compared symmetrised: the response's entries
         ! are meaningful only as transpose pairs, the tangent's are handed out
         ! already symmetrised
         if (maxval(abs(rt%w_value - response%lsf%w_value)) > ULP_TOL*maxval(abs(rt%w_value)) &
             .or. maxval(abs(rt%w_gradient - response%lsf%w_gradient)) &
             > ULP_TOL*maxval(abs(rt%w_gradient)) &
             .or. maxval(abs(rt%w_hessian - symmetrised(response%lsf%w_hessian))) &
             > ULP_TOL*maxval(abs(rt%w_hessian))) then
            call test_failed(error, "the base level-set weights differ from the potential-path"// &
                             " response")
            return
         end if
      end if

      nw = 13*ngrid
      allocate (ref(nw, size(dirs, 3)), val(nw, size(dirs, 3)))
      do idir = 1, size(dirs, 3)
         val(1:ngrid, idir) = rt%dw_value(:, idir)
         val(ngrid + 1:4*ngrid, idir) = reshape(rt%dw_gradient(:, :, idir), [3*ngrid])
         val(4*ngrid + 1:, idir) = reshape(rt%dw_hessian(:, :, :, idir), [9*ngrid])
      end do
      do istep = 1, size(FD_STEPS)
         do idir = 1, size(dirs, 3)
            call differenced(kind, mol, dirs(:, :, idir), FD_STEPS(istep), weights_at, &
                             ref(:, idir), error)
            if (allocated(error)) return
         end do
         if (maxval(abs(ref)) <= VACUITY_THR) then
            call test_failed(error, "the differenced weights are vacuous")
            return
         end if
         call worst_deviation(reshape(val, [size(val)]), reshape(ref, [size(ref)]), FD_TOL, &
                              worst, worst_ref)
         if (worst > 0.0_wp) then
            call test_failed(error, "weight tangent mismatch (h = "//to_string(FD_STEPS(istep))// &
                             "): worst deviation "//to_string(worst)//" against reference "// &
                             to_string(worst_ref))
            return
         end if
      end do
      if (allocated(base)) deallocate (base)
   end subroutine run_weight_fd

   !> Five-point central difference of a per-point quantity along one direction
   !>
   !> The model is rebuilt at `R + t v` for the four offsets of the stencil,
   !> its grid identity checked against the base geometry, and the quantity
   !> read through `probe`.
   subroutine differenced(kind, mol, dir, step, probe, deriv, error)
      integer, intent(in) :: kind
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: dir(:, :)
      real(wp), intent(in) :: step
      interface
         subroutine probe(model, coupling, values, error)
            import :: solvation_model_general, coupling_type, wp, error_type
            type(solvation_model_general), intent(inout) :: model
            type(coupling_type), intent(in) :: coupling
            real(wp), allocatable, intent(out) :: values(:)
            type(error_type), allocatable, intent(out) :: error
         end subroutine probe
      end interface
      real(wp), intent(out) :: deriv(:)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: OFFSET(4) = [-2, -1, 1, 2]
      real(wp), parameter :: COEFF(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]/12.0_wp

      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_disp
      type(coupling_type) :: coupling
      type(mctc_error), allocatable :: err
      integer, allocatable :: ref_numbering(:), ref_owner(:), numbering(:), owner(:)
      real(wp), allocatable :: values(:)
      integer :: ioff

      call build_model(kind, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      call grid_identity(model, ref_numbering, ref_owner, error)
      if (allocated(error)) return

      deriv = 0.0_wp
      do ioff = 1, size(OFFSET)
         mol_disp = mol
         mol_disp%xyz = mol%xyz + real(OFFSET(ioff), wp)*step*dir
         call model%update(mol_disp, err)
         if (allocated(err)) then
            call test_failed(error, "model update on the displaced geometry failed: "// &
                             err%message)
            return
         end if
         call grid_identity(model, numbering, owner, error)
         if (allocated(error)) return
         if (size(numbering) /= size(ref_numbering)) then
            call test_failed(error, "grid changed under the step "//to_string(step))
            return
         end if
         if (any(numbering /= ref_numbering) .or. any(owner /= ref_owner)) then
            call test_failed(error, "grid changed under the step "//to_string(step))
            return
         end if
         call probe(model, coupling, values, error)
         if (allocated(error)) return
         if (size(values) /= size(deriv)) then
            call test_failed(error, "the probed quantity changed size under the step")
            return
         end if
         deriv = deriv + COEFF(ioff)*values/step
      end do
   end subroutine differenced

   !> The surface charges of the model at its current geometry
   subroutine charges_at(model, coupling, values, error)
      type(solvation_model_general), intent(inout) :: model
      type(coupling_type), intent(in) :: coupling
      real(wp), allocatable, intent(out) :: values(:)
      type(error_type), allocatable, intent(out) :: error

      type(response_type) :: response
      type(mctc_error), allocatable :: err

      call model%get_trace_response(coupling, response, err)
      if (allocated(err)) then
         call test_failed(error, "trace response failed: "//err%message)
         return
      end if
      if (.not. allocated(response%electrostatics%surface_charge)) then
         call test_failed(error, "no surface charges were returned")
         return
      end if
      values = response%electrostatics%surface_charge
   end subroutine charges_at

   !> The gradient-path level-set weights of the model at its current geometry
   !>
   !> Read off a response tangent of a single direction, so that the base
   !> weights differenced here are the very ones the tangents claim to move.
   subroutine weights_at(model, coupling, values, error)
      type(solvation_model_general), intent(inout) :: model
      type(coupling_type), intent(in) :: coupling
      real(wp), allocatable, intent(out) :: values(:)
      type(error_type), allocatable, intent(out) :: error

      type(response_tangent_type) :: rt
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :)
      integer :: ngrid

      ngrid = model%cavity%ngrid
      allocate (dirs(ndim, model%cavity%nsph, 1), source=0.0_wp)
      dirs(1, 1, 1) = 1.0_wp
      allocate (hvp(ndim, model%cavity%nsph, 1), source=0.0_wp)
      call rt%init(ngrid, 1, want_lsf=.true.)
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (allocated(err)) then
         call test_failed(error, "Hessian-vector product on the displaced geometry failed: "// &
                          err%message)
         return
      end if
      allocate (values(13*ngrid))
      values(1:ngrid) = rt%w_value
      values(ngrid + 1:4*ngrid) = reshape(rt%w_gradient, [3*ngrid])
      values(4*ngrid + 1:) = reshape(rt%w_hessian, [9*ngrid])
   end subroutine weights_at

   !* ================================================================================= *!
   !*                                     Guards                                       *!
   !* ================================================================================= *!

   !> Mis-sized containers are refused; one without level-set channels keeps
   !> the rank-4 form and its bit-exact columns
   subroutine test_guards(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(solvation_model_general) :: model
      type(moist_context_type), target :: ctx
      type(coupling_type) :: coupling
      type(response_tangent_type) :: rt
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: dirs(:, :, :), hvp(:, :, :), plain(:, :, :)
      integer :: nsph, ngrid

      call drop_fixture_geometry(FIX_PLAIN, mol)
      call build_model(MODEL_CPCM, mol, model, ctx, coupling, error)
      if (allocated(error)) return
      nsph = model%cavity%nsph
      ngrid = model%cavity%ngrid
      call unit_directions(nsph, dirs)
      allocate (hvp(ndim, nsph, size(dirs, 3)), source=0.0_wp)
      allocate (plain(ndim, nsph, size(dirs, 3)), source=0.0_wp)

      ! Uninitialised
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (.not. allocated(err)) then
         call test_failed(error, "an uninitialised response tangent was accepted")
         return
      end if
      deallocate (err)
      if (any(hvp /= 0.0_wp)) then
         call test_failed(error, "a refused call touched the accumulator")
         return
      end if

      ! Wrong direction count
      call rt%init(ngrid, size(dirs, 3) - 1, want_lsf=.false.)
      call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      if (.not. allocated(err)) then
         call test_failed(error, "a mis-sized response tangent was accepted")
         return
      end if
      deallocate (err)

      ! Charge and position channels alone: rank-4 stays, columns to the bit
      call model%get_hvp(coupling, dirs, plain, err)
      if (.not. allocated(err)) then
         call rt%init(ngrid, size(dirs, 3), want_lsf=.false.)
         call model%get_hvp(coupling, dirs, hvp, err, rt=rt)
      end if
      if (allocated(err)) then
         call test_failed(error, "Hessian-vector products failed: "//err%message)
         return
      end if
      if (any(hvp /= plain)) then
         call test_failed(error, "a response tangent without level-set channels moved the"// &
                          " rank-4 columns")
         return
      end if
      if (rt%want_lsf()) then
         call test_failed(error, "level-set channels appeared unasked")
         return
      end if
      if (maxval(abs(rt%surface_charge)) <= VACUITY_THR) then
         call test_failed(error, "the charge tangent is vacuous")
         return
      end if
   end subroutine test_guards

end module test_model_general_hessian_coupled
