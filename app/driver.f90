
!> Entry point for running single point calculations with moist
module moist_driver
   use, intrinsic :: iso_fortran_env, only: output_unit, input_unit
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type, read_structure, filetype
   use mctc_io_utils, only: to_lower
   use moist_cli, only: run_config
   use moist_output_ascii, only: moist_header, moist_build_header, cavity_header
   use moist_data_solvents, only: get_solvent_id
   use moist_data_solvents, only: solvation_system_parameters, new_solvation_system_parameters
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_marchingcubes, only: cavity_type_marchingcubes, new_cavity_marchingcubes
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_radii, only: radius_type, new_radii
   use moist_type, only: coupling_type, cavity_type, solvation_model
#ifdef WITH_RISM
   ! use moist_model_rism1d, only: rism1d_model, new_rism1d_model
   ! use moist_model_rism3d, only: rism3d_model, new_rism3d_model
#endif
   use moist_context, only: moist_context_type, new_context
!$ use omp_lib
! #ifdef WITH_MKL
! !$ use mkl_service
! #endif

   implicit none(type, external)
   private

   public :: main

contains

!> Main entry point for the driver
   subroutine main(config, error)

      !> Configuration for this driver
      type(run_config), intent(in) :: config

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call run_main(config, error)
   end subroutine main

!> Entry point for the single point driver
   subroutine run_main(config, error)

      !> Configuration for this driver
      type(run_config), intent(in) :: config

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      !> Solvation model
      class(solvation_model), allocatable :: sm

      !> Shared run context owned for the whole run borrowed by everything
      type(moist_context_type), target :: ctx

      !> Solvation system type
      type(solvation_system_parameters), allocatable :: system

      character(len=:), allocatable :: filename
      character(len=:), allocatable :: solvent

      real(wp) :: energy
      type(coupling_type) :: coupling

      integer :: solvent_id

      !> File reading and parsing
      integer :: unit, stat
      character(len=256) :: tmp
      real(wp) :: tmp_wp

      !> Open the run context; every cavity/model borrows it, so all sub-timers
      !> started later nest under the "total" node opened here.
      call new_context(ctx, verbosity=config%verbosity, debug=config%debug)
      call ctx%timer%start("total")

      !* ---------------------------- Thread configuration --------------------------- *!
      !> Routed through the context so it stays the single place the thread budget
      !> is changed -- `get_num_threads` and `print_settings` then report what the
      !> run is actually using, and the pin is released again on `ctx%delete`.
      if (config%num_threads > 0) then
!$       if (.false.) then
            write (ctx%unit, '(a)') &
               "[Warn] Program compiled without OpenMP support, ignoring --threads"
!$       else
!$          call ctx%set_num_threads(config%num_threads)
! #ifdef WITH_MKL
! !$       call mkl_set_num_threads(config%num_threads)
! #endif
!$          if (config%verbosity > 0) then
!$             write (ctx%unit, '(a,i0,a)') &
!$                "[Info] OpenMP threads set to ", config%num_threads, ""
!$          end if
!$       end if
      end if

      !* ---------------------------- Solvent inspection ----------------------------- *!
      if (to_lower(config%mode) == "solvent") then
         if (len_trim(config%solvent) == 0) then
            call fatal_error(error, "No solvent specified")
            return
         end if

         call get_solvent_id(trim(config%solvent), solvent_id, error)
         if (allocated(error)) return
         allocate (system)
         call new_solvation_system_parameters(system, solvent_id, &
                                              temperature=config%temperature, pressure_si=config%pressure_si, &
                                              error=error)
         if (allocated(error)) return
         call system%print()
         return
      end if

      !* -------------------------------- Read input file -------------------------------- *!

      ! Read input file
      if (config%input == "-") then
         if (.not. allocated(config%input_format)) then
            call read_structure(mol, input_unit, filetype%xyz, error)
         else
            call read_structure(mol, input_unit, config%input_format, error)
         end if
      else
         call read_structure(mol, config%input, error, config%input_format)
      end if
      if (allocated(error)) return

      !* ---------------------------------- Read charge ---------------------------------- *!

      ! Get charge
      if (allocated(config%charge)) then
         mol%charge = config%charge
      else
         filename = join(dirname(config%input), ".CHRG")
         if (exists(filename)) then
            open (file=filename, newunit=unit)
            read (unit, *, iostat=stat) tmp_wp
            if (stat == 0) then
               mol%charge = tmp_wp
               if (config%verbosity > 0) write (output_unit, '(a)') &
                  "[Info] Molecular charge read from '"//filename//"'"
            else
               if (config%verbosity > 0) write (output_unit, '(a)') &
                  "[Warn] Could not read molecular charge read from '"//filename//"'"
            end if
            close (unit)
         end if
      end if

      !* ================================================================================= *!
      !*                                      Cavities                                     *!
      !* ================================================================================= *!

      ! Cavity-only modes: build cavity grid and optionally write output
      if (to_lower(config%mode) == "numsa" .or. &
          to_lower(config%mode) == "drop" .or. &
          to_lower(config%mode) == "mc" .or. &
          to_lower(config%mode) == "iswig") then

         block
            class(cavity_type), allocatable :: cavity
            class(radius_type), allocatable :: radius_model

            ! Print header
            call cavity_header(output_unit, trim(config%mode))

            ! Instantiate the appropriate cavity type
            if (to_lower(config%mode) == "numsa") then
               block
                  type(cavity_type_numsa), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  call new_cavity_numsa(tmp_cavity, ctx, nleb=config%nleb, &
                                        radii=radius_model, error=error)
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "iswig") then
               block
                  type(cavity_type_iswig), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  call new_cavity_iswig(tmp_cavity, ctx, nleb=config%nleb, &
                                        radius_model=radius_model, error=error)
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "drop") then
               block
                  type(cavity_type_drop), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  if (to_lower(config%drop_variant) == "svdw") then
                     block
                        type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
                        !> The cavity couples the LSF screening threshold to
                        !> its own tolerance (passed below as `tolerance`),
                        !> so we only forward the shape parameters here.
                        call svdw_template%new( &
                           blend_k=config%drop_blend_k, &
                           blend_1b=config%drop_blend_1b, &
                           blend_2b=config%drop_blend_2b, &
                           blend_3b=config%drop_blend_3b)
                        call new_cavity_drop(tmp_cavity, ctx, &
                                             nleb=config%nleb, &
                                             tolerance=config%drop_tol, proj_level=config%drop_proj_level, &
                                             wleb_prune_level=config%drop_wleb_prune_level, &
                                             radius_model=radius_model, &
                                             lsf_model=svdw_template, error=error)
                     end block
                  else if (to_lower(config%drop_variant) == "cfc") then
                     block
                        type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
                        call cfc_template%new(a1=config%cfc_a1, a2=config%cfc_a2, &
                                              c=config%cfc_c, m=config%cfc_m, screen_k=config%cfc_screen_k)
                        call new_cavity_drop(tmp_cavity, ctx, &
                                             nleb=config%nleb, &
                                             tolerance=config%drop_tol, proj_level=config%drop_proj_level, &
                                             wleb_prune_level=config%drop_wleb_prune_level, &
                                             radius_model=radius_model, &
                                             lsf_model=cfc_template, error=error)
                     end block
                  else
                     call fatal_error(error, "Unknown DROP variant: "//trim(config%drop_variant))
                  end if
                  if (allocated(error)) return
                  call tmp_cavity%properties(do_fine=config%cavity_fine)
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "mc") then
               block
                  type(cavity_type_marchingcubes), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  !> Marching cubes writes its mesh while integrating, so the
                  !> export paths have to be known before update runs.
                  if (to_lower(config%drop_variant) == "svdw") then
                     block
                        type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
                        call svdw_template%new( &
                           blend_k=config%drop_blend_k, &
                           blend_1b=config%drop_blend_1b, &
                           blend_2b=config%drop_blend_2b, &
                           blend_3b=config%drop_blend_3b)
                        call new_mc_cavity(tmp_cavity, ctx, radius_model, svdw_template, &
                                           config%cavity_mc_spacing, config%dump, error)
                     end block
                  else if (to_lower(config%drop_variant) == "cfc") then
                     block
                        type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
                        call cfc_template%new(a1=config%cfc_a1, a2=config%cfc_a2, &
                                              c=config%cfc_c, m=config%cfc_m, screen_k=config%cfc_screen_k)
                        call new_mc_cavity(tmp_cavity, ctx, radius_model, cfc_template, &
                                           config%cavity_mc_spacing, config%dump, error)
                     end block
                  else
                     call fatal_error(error, "Unknown level set function: "//trim(config%drop_variant))
                  end if
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            end if

            ! Use polymorphic cavity methods
            call cavity%update(mol, error=error)
            if (allocated(error)) return

            ! Gradient if asked for
            if (config%grad) then
               call cavity%get_gradient()
               if (allocated(cavity%error)) then
                  allocate (error, source=cavity%error)
                  return
               end if
            end if

            ! Print results; with no unit the cavity follows its own context
            call cavity%print()

            ! Write cavity files (xyz, csv, pqr) only when --dump is given.
            ! Marching cubes has no grid points to dump; it already wrote its
            ! triangle mesh (cavity.obj/cavity.pqr) during update.
            if (config%dump .and. to_lower(config%mode) /= "mc") then
               call cavity%write_xyz_debug('cavity.xyz', error=error)
               if (allocated(error)) return
               call cavity%write_csv_debug('cavity.csv', error=error)
               if (allocated(error)) return
               call cavity%write_pqr_debug('cavity.pqr', error=error)
               if (allocated(error)) return
            end if

         end block
         call report_run_timings()
         call ctx%delete()
         return
      end if

      !* ---------------------------------- Get Solvent ---------------------------------- *!

      ! Get solvent
      if (len_trim(config%solvent) > 0) then
         solvent = trim(config%solvent)
      else
         filename = join(dirname(config%input), ".SOLVENT")
         if (exists(filename)) then
            open (file=filename, newunit=unit)
            read (unit, '(A)', iostat=stat) tmp
            close (unit)
            if (stat == 0) then
               ! allocate exactly to the trimmed length
               allocate (character(len=len_trim(tmp)) :: solvent)
               solvent = trim(tmp)
            end if
         end if
      end if
      if (.not. allocated(solvent)) then
         call fatal_error(error, "No solvent specified")
         return
      end if

      call get_solvent_id(solvent, solvent_id, error)
      if (allocated(error)) return
      allocate (system)
      call new_solvation_system_parameters(system, solvent_id, &
                                           temperature=config%temperature, pressure_si=config%pressure_si, &
                                           error=error)
      if (allocated(error)) return

      !* ================================================================================= *!
      !*                                  Solvation models                                 *!
      !* ================================================================================= *!

      ! Exit, no models implemented in the current preview version
      if ((to_lower(config%mode) == "gems") .or. to_lower(config%mode) == "rism1d" .or. &
          to_lower(config%mode) == "rism3d" .or. to_lower(config%mode) == "alpb") then
         call fatal_error(error, "No solvation models implemented in the current preview version")
         return
      end if

      call report_run_timings()
      call ctx%delete()

   contains

      !> Stop the top-level timer and print the accumulated timing tree. Every
      !> cavity/model sub-timer nests under "total" (opened at run start), so this
      !> prints the full hierarchical breakdown for the run.
      subroutine report_run_timings()
         call ctx%timer%stop("total")
         if (config%verbosity > 0) &
            call ctx%timer%write(output_unit, "moist", max_depth=ctx%report_depth())
      end subroutine report_run_timings

   end subroutine run_main

   !> Construct a marching-cubes cavity, wiring up mesh export when --dump is set
   !>
   !> The integrator writes its triangle mesh during `update`, so the export paths
   !> have to be decided at construction; this wrapper keeps that branch out of the
   !> per-level-set blocks above.
   !>
   !> @param[inout] cavity        Cavity instance to initialize
   !> @param[in]    ctx           Shared run context
   !> @param[in]    radius_model  Atomic radius model
   !> @param[in]    lsf_model     Level set function template
   !> @param[in]    spacing       Finest grid spacing in bohr
   !> @param[in]    dump          Write the triangle mesh to cavity.obj/cavity.pqr
   !> @param[out]   error         Error handling
   subroutine new_mc_cavity(cavity, ctx, radius_model, lsf_model, spacing, dump, error)
      type(cavity_type_marchingcubes), intent(inout) :: cavity
      type(moist_context_type), intent(in), target :: ctx
      class(radius_type), intent(in) :: radius_model
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      real(wp), intent(in) :: spacing
      logical, intent(in) :: dump
      type(error_type), allocatable, intent(out) :: error

      if (dump) then
         call new_cavity_marchingcubes(cavity, ctx, radius_model=radius_model, &
                                       lsf_model=lsf_model, spacing=spacing, &
                                       obj_file='cavity.obj', pqr_file='cavity.pqr', &
                                       error=error)
      else
         call new_cavity_marchingcubes(cavity, ctx, radius_model=radius_model, &
                                       lsf_model=lsf_model, spacing=spacing, &
                                       error=error)
      end if

   end subroutine new_mc_cavity

!> Construct path by joining strings with os file separator
   function join(a1, a2) result(path)
      use mctc_env_system, only: is_windows
      character(len=*), intent(in) :: a1, a2
      character(len=:), allocatable :: path
      character :: filesep

      if (is_windows()) then
         filesep = '\'
      else
         filesep = '/'
      end if

      path = a1//filesep//a2
   end function join

!> test if pathname already exists
   function exists(filename)
      character(len=*), intent(in) :: filename
      logical :: exists
      inquire (file=filename, exist=exists)
   end function exists

!> Extract dirname from path
   function dirname(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: dirname

      dirname = filename(1:scan(filename, "/\", back=.true.))
      if (len_trim(dirname) == 0) dirname = "."
   end function dirname

end module moist_driver
