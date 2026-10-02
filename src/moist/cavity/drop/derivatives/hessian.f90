!> Public DROP Hessian accessors
!>
!> - `d/dv [J^T omega] = (dJ^T/dv) omega + J^T (d omega/dv)`: fixed plus
!>   response channel, both run by [[drop_hessian_traverse]]
!> - Owns only the form of the fixed channel and the two entry points,
!>   [[get_surface_hessian_drop]] (products) and [[get_hessian_drop]] (dense)
!> - Fixed channel forms: direction-free rank-4 block or one column per
!>   supplied direction, chosen by [[hvp_fixed_mode]]
!> - Per direction up to `drop_hvp_per_dir_max` directions, rank-4 beyond;
!>   measured crossover at 12 to 18 directions
!> - Constant bound errs towards rank-4, whose cost is bounded by the dense path
!> - Rank-4 memory: dense `(3, nsph, 3, nsph)` staging block plus the
!>   traversal's per-thread sparse accumulators
!> - `host` or level-set weight tangents in `rt` force the per-direction form
!> - Dense block: all `3 nsph` unit directions in one response batch,
!>   bit-for-bit the product path along the same unit directions
!> - Surface adjoints folded once, in [[surface_hessian_halves]]: fixed half
!>   takes `eff` alone, response half takes `eff` and the raw `acc`
!> - Multi-branch grids: response half moves `branch_phi_adj`, fixed half
!>   differentiates at fixed adjoint with the branch weight moving
!> - Neither entry point restricts the grid
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !* ================================================================================= *!
   !*                            Public Hessian accessors                               *!
   !* ================================================================================= *!

   !> Hessian-vector products of the DROP surface contribution
   !>
   !> - One gradient column `d/dv [J^T omega]` per supplied direction
   !> - Added to `hvp`; `hvp` untouched on failure
   !> - With `omega_v`: columns of the full model Hessian; without it:
   !>   frozen-adjoint columns
   !> - `host` or level-set weights in `rt` force the per-direction fixed half
   !>   whatever the count, see [[host_exchange_per_dir]]
   !> - `host` and `rt` need `omega_v`, see [[check_host_exchange]]
   !>
   !> @param[in]     self     DROP cavity instance holding a projected grid
   !> @param[in]     acc      accumulated surface-observable adjoints
   !> @param[in]     dirs     nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] hvp      Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]    error    error object, allocated on failure
   !> @param[in,out] omega_v  surface-adjoint response of the model, optional
   !> @param[in,out] host     second-order host exchange, optional
   !> @param[in,out] rt       response tangent of the direction set, optional
   module subroutine get_surface_hessian_drop(self, acc, dirs, hvp, error, omega_v, host, rt)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v
      !> Second-order host exchange
      class(coupling_tangent_type), intent(inout), optional :: host
      !> Response tangent of the direction set
      type(response_tangent_type), intent(inout), optional :: rt

      !> Direction-free fixed half (rank-4 form only) and staged columns
      real(wp), allocatable :: hess_fixed(:, :, :, :), total(:, :, :)
      !> Form of the fixed channel
      integer :: fixed_mode
      !> Extents and loop indices
      integer :: ndir, idir, iatom, iaxis

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      call check_direction_set(self, dirs, hvp, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      call check_host_exchange(self, dirs, "get_surface_hessian_drop", error, omega_v, host, rt)
      if (allocated(error)) return
      if (self%ngrid <= 0) return
      ndir = size(dirs, 3)

      !* --------------------------------- Both halves -------------------------------- *!
      fixed_mode = hvp_fixed_mode(ndir, self%nsph)
      if (host_exchange_per_dir(host, rt)) fixed_mode = drop_fixed_per_dir
      allocate (total(ndim, self%nsph, ndir), source=0.0_wp)

      if (fixed_mode == drop_fixed_per_dir) then
         ! Both channels land in `total`
         call surface_hessian_halves(self, acc, dirs, fixed_mode, "get_surface_hessian_drop", &
                                     total, error, omega_v=omega_v, host=host, rt=rt)
         if (allocated(error)) return
      else
         allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
         call surface_hessian_halves(self, acc, dirs, fixed_mode, "get_surface_hessian_drop", &
                                     total, error, hess_fixed=hess_fixed, omega_v=omega_v, &
                                     rt=rt)
         if (allocated(error)) return

         !* -------------------------- Contract the fixed half ------------------------ *!
         ! On top of the response half already in `total`
         do idir = 1, ndir
            do iatom = 1, self%nsph
               do iaxis = 1, ndim
                  total(:, :, idir) = total(:, :, idir) &
                                      + hess_fixed(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
               end do
            end do
         end do
      end if

      hvp = hvp + total
   end subroutine get_surface_hessian_drop

   !> Dense nuclear Hessian of the DROP surface contribution
   !>
   !> - Added to `hessian`; `hessian` untouched on failure
   !> - Column `(beta, B)`: Hessian-vector product along unit direction
   !>   `e_(beta, B)`
   !> - All `3 nsph` unit directions go to the response half in one call
   !> - Rank-4 fixed half added column by column, no contraction
   !> - `host` or level-set weights in `rt` run the fixed half per direction
   !>
   !> @param[in]     self     DROP cavity instance holding a projected grid
   !> @param[in]     acc      accumulated surface-observable adjoints
   !> @param[in,out] hessian  nuclear-Hessian accumulator `(3, nsph, 3, nsph)`
   !> @param[out]    error    error object, allocated on failure
   !> @param[in,out] omega_v  surface-adjoint response of the model, optional
   !> @param[in,out] host     second-order host exchange, optional
   !> @param[in,out] rt       response tangent of the Cartesian basis, optional
   module subroutine get_hessian_drop(self, acc, hessian, error, omega_v, host, rt)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v
      !> Second-order host exchange
      class(coupling_tangent_type), intent(inout), optional :: host
      !> Response tangent of the Cartesian basis
      type(response_tangent_type), intent(inout), optional :: rt

      !> Cartesian unit directions, one per nuclear degree of freedom
      real(wp), allocatable :: dirs(:, :, :)
      !> Direction-free fixed half and per-column response half
      real(wp), allocatable :: hess_fixed(:, :, :, :), resp(:, :, :)
      !> Extents and loop indices
      integer :: ndir, idir, iatom, iaxis

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_hessian_drop", error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [ndim, self%nsph, ndim, self%nsph])) then
         call fatal_error(error, "get_hessian_drop: hessian shape mismatch")
         return
      end if
      if (self%ngrid <= 0 .or. self%nsph <= 0) return

      !* -------------------------- Cartesian unit directions ------------------------- *!
      ndir = ndim*self%nsph
      allocate (dirs(ndim, self%nsph, ndir), source=0.0_wp)
      do iatom = 1, self%nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do
      call check_host_exchange(self, dirs, "get_hessian_drop", error, omega_v, host, rt)
      if (allocated(error)) return

      !* --------------------------------- Both halves -------------------------------- *!
      allocate (resp(ndim, self%nsph, ndir), source=0.0_wp)

      if (host_exchange_per_dir(host, rt)) then
         ! Fixed half per direction: both halves land in `resp`, one block
         ! column per unit direction
         call surface_hessian_halves(self, acc, dirs, drop_fixed_per_dir, "get_hessian_drop", &
                                     resp, error, omega_v=omega_v, host=host, rt=rt)
         if (allocated(error)) return
         do iatom = 1, self%nsph
            do iaxis = 1, ndim
               idir = ndim*(iatom - 1) + iaxis
               hessian(:, :, iaxis, iatom) = hessian(:, :, iaxis, iatom) + resp(:, :, idir)
            end do
         end do
         return
      end if

      allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
      call surface_hessian_halves(self, acc, dirs, drop_fixed_rank4, "get_hessian_drop", &
                                  resp, error, hess_fixed=hess_fixed, omega_v=omega_v, rt=rt)
      if (allocated(error)) return

      do iatom = 1, self%nsph
         do iaxis = 1, ndim
            idir = ndim*(iatom - 1) + iaxis
            hessian(:, :, iaxis, iatom) = hessian(:, :, iaxis, iatom) &
                                          + hess_fixed(:, :, iaxis, iatom) + resp(:, :, idir)
         end do
      end do
   end subroutine get_hessian_drop

   !> Whether the second-order host exchange forces the per-direction fixed half
   !>
   !> - Host jet tangents are not contractions of the point tensors with the
   !>   direction, as the rank-4 chain basis assumes
   !> - Level-set weight tangents of `rt` exist only in the per-direction chain
   !>
   !> @param[in] host  second-order host exchange, optional
   !> @param[in] rt    response tangent, optional
   logical function host_exchange_per_dir(host, rt) result(per_dir)
      !> Second-order host exchange
      class(coupling_tangent_type), intent(in), optional :: host
      !> Response tangent
      type(response_tangent_type), intent(in), optional :: rt

      per_dir = present(host)
      if (present(rt)) per_dir = per_dir .or. rt%want_lsf()
   end function host_exchange_per_dir

   !> Check the second-order host exchange against the direction set
   !>
   !> - `host` and `rt` both need `omega_v`: it carries the host field tangents
   !>   to the components and their charge tangents back
   !> - `rt` must be initialised for this grid and direction set
   !> - No check without `host` and `rt`
   !>
   !> @param[in]  self     DROP cavity instance
   !> @param[in]  dirs     nuclear directions `(3, nsph, ndir)`
   !> @param[in]  context  calling routine, prefix of the diagnostics
   !> @param[out] error    error object, allocated on a mismatch
   !> @param[in]  omega_v  surface-adjoint response of the model, optional
   !> @param[in]  host     second-order host exchange, optional
   !> @param[in]  rt       response tangent, optional
   subroutine check_host_exchange(self, dirs, context, error, omega_v, host, rt)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(in), optional :: omega_v
      !> Second-order host exchange
      class(coupling_tangent_type), intent(in), optional :: host
      !> Response tangent
      type(response_tangent_type), intent(in), optional :: rt

      if (.not. (present(host) .or. present(rt))) return
      if (.not. present(omega_v)) then
         call fatal_error(error, context//": the second-order host exchange needs the"// &
                          " model's adjoint response")
         return
      end if
      if (present(rt)) then
         if (.not. rt%is_initialized(self%ngrid, size(dirs, 3))) then
            call fatal_error(error, context//": the response tangent is not initialised"// &
                             " for this grid and direction set")
            return
         end if
      end if
   end subroutine check_host_exchange

   !* ================================================================================= *!
   !*                              Composition of the halves                            *!
   !* ================================================================================= *!

   !> Form of the fixed channel for a Hessian-vector product of `ndir` directions
   !>
   !> - `drop_fixed_per_dir` for `ndir <= drop_hvp_per_dir_max` and
   !>   `ndir < 3 nsph`, `drop_fixed_rank4` otherwise
   !> - Bound explained in the module header
   !>
   !> @param[in] ndir  directions asked for
   !> @param[in] nsph  spheres of the cavity
   pure function hvp_fixed_mode(ndir, nsph) result(mode)
      !> Directions asked for
      integer, intent(in) :: ndir
      !> Spheres of the cavity
      integer, intent(in) :: nsph
      !> Form of the fixed channel, `drop_fixed_per_dir` or `drop_fixed_rank4`
      integer :: mode

      if (ndir <= drop_hvp_per_dir_max .and. ndir < ndim*nsph) then
         mode = drop_fixed_per_dir
      else
         mode = drop_fixed_rank4
      end if
   end function hvp_fixed_mode

   !> Evaluate both halves of the surface Hessian into caller-owned buffers
   !>
   !> - Only place the two halves meet and the surface adjoints are folded
   !> - Buffers are added to: zeroed on entry, the caller's staging buffers and
   !>   not its public accumulators
   !> - `columns`: response half, plus the fixed half in per-direction mode
   !> - `hess_fixed`: rank-4 fixed half, present exactly in rank-4 mode
   !> - `context`: public entry point, named in the error of a singular
   !>   bordered system
   !>
   !> @param[in]     self        DROP cavity instance
   !> @param[in]     acc         accumulated surface-observable adjoints
   !> @param[in]     dirs        nuclear directions `(3, nsph, ndir)`
   !> @param[in]     fixed_mode  form of the fixed channel
   !> @param[in]     context     calling routine, prefix of the diagnostics
   !> @param[in,out] columns     per-direction half or halves `(3, nsph, ndir)`
   !> @param[out]    error       error object, allocated on failure
   !> @param[in,out] hess_fixed  rank-4 fixed half `(3, nsph, 3, nsph)`, optional
   !> @param[in,out] omega_v     surface-adjoint response of the model, optional
   !> @param[in,out] host        second-order host exchange, optional
   !> @param[in,out] rt          response tangent of the direction set, optional
   subroutine surface_hessian_halves(self, acc, dirs, fixed_mode, context, columns, error, &
                                     hess_fixed, omega_v, host, rt)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Form of the fixed channel
      integer, intent(in) :: fixed_mode
      !> Calling routine, the public entry point named on failure
      character(len=*), intent(in) :: context
      !> Per-direction columns
      real(wp), intent(inout) :: columns(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Rank-4 fixed half
      real(wp), intent(inout), optional :: hess_fixed(:, :, :, :)
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v
      !> Second-order host exchange
      class(coupling_tangent_type), intent(inout), optional :: host
      !> Response tangent of the direction set
      type(response_tangent_type), intent(inout), optional :: rt

      !> Folded surface adjoints of the base geometry, read by both channels
      type(drop_surface_weights_type) :: eff

      ! `fold_switching = .true.` as on the nuclear path: nuclear motion moves
      ! `f`, so the area channel's `da/df` enters the switching adjoint
      call prepare_surface_weights(self, acc, .true., eff)

      ! Fixed channel `(dPhi/dv) . eff` reads `eff` only; response channel
      ! `Phi . (d eff/dv)` also needs the raw `acc`
      ! Absent `hess_fixed` passes through as absent
      call drop_hessian_traverse(self, eff, fixed_mode, .true., context, &
                                 acc=acc, dirs=dirs, hess_fixed=hess_fixed, hvp=columns, &
                                 error=error, omega_v=omega_v, host=host, rt=rt)
   end subroutine surface_hessian_halves

end submodule moist_cavity_drop_derivatives_hessian
