
!> Entry point for running single point calculations with moist
module moist_driver
   use, intrinsic :: iso_fortran_env, only: output_unit, input_unit
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type, read_structure, filetype
   use mctc_io_utils, only: to_lower
   use moist_cli, only: run_config
   use moist_output_ascii, only: moist_header, moist_build_header, gems_header, cavity_header
   use moist_data_solvents, only: get_solvent_id
   use moist_data_solvents, only: solvation_system_parameters, new_solvation_system_parameters
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_radii, only: radius_type, new_radii
   use moist_type, only: wavefunction_type, cavity_type, solvation_model
   ! use moist_model_gems, only: gems_model, new_gems_model
#ifdef WITH_RISM
   ! use moist_model_rism1d, only: rism1d_model, new_rism1d_model
   ! use moist_model_rism3d, only: rism3d_model, new_rism3d_model
#endif
   use moist_utils_timer, only: timer_type
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

      !> Timer for performance measurement
      type(timer_type) :: timer

      !> Solvation system type
      type(solvation_system_parameters), allocatable :: system

      character(len=:), allocatable :: filename
      character(len=:), allocatable :: solvent

      real(wp) :: energy
      type(wavefunction_type) :: wfn

      integer :: solvent_id

      !> File reading and parsing
      integer :: unit, stat
      character(len=256) :: tmp
      real(wp) :: tmp_wp

      call timer%new(2, .true.)
      call timer%measure(1, 'total')

      !* ---------------------------- Thread configuration --------------------------- *!
      if (config%num_threads > 0) then
!$       if (.false.) then
            write (output_unit, '(a)') &
               "[Warn] Program compiled without OpenMP support, ignoring --threads"
!$       else
!$          call omp_set_num_threads(config%num_threads)
! #ifdef WITH_MKL
! !$       call mkl_set_num_threads(config%num_threads)
! #endif
!$          if (config%verbosity > 0) then
!$             write (output_unit, '(a,i0,a)') &
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
                  call new_cavity_numsa(tmp_cavity, nleb=config%nleb, &
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
                  call new_cavity_iswig(tmp_cavity, nleb=config%nleb, &
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
                           blend_k =config%drop_blend_k, &
                           blend_1b=config%drop_blend_1b, &
                           blend_2b=config%drop_blend_2b, &
                           blend_3b=config%drop_blend_3b)
                        call new_cavity_drop(tmp_cavity, &
                                            debug=config%debug, verbose=config%verbosity, &
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
                        call new_cavity_drop(tmp_cavity, &
                                            debug=config%debug, verbose=config%verbosity, &
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
                  call tmp_cavity%properties(do_fine=config%cavity_fine, &
                                             do_mc=config%cavity_mc)
                  call move_alloc(tmp_cavity, cavity)
               end block
            end if

            ! Use polymorphic cavity methods
            call cavity%update(mol, error=error)
            if (allocated(error)) return

            ! Gradient if asked for
            if (config%grad) then
               call cavity%get_gradient()
            end if

            ! Print results
            call cavity%print(output_unit)
            if (config%verbosity > 2) call cavity%print_fine(output_unit)

            ! Write cavity files (xyz, csv, pqr) only when --dump is given
            if (config%dump) then
               call cavity%write_xyz_debug('cavity.xyz', error=error)
               if (allocated(error)) return
               call cavity%write_csv_debug('cavity.csv', error=error)
               if (allocated(error)) return
               call cavity%write_pqr_debug('cavity.pqr', error=error)
               if (allocated(error)) return
            end if

         end block
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

!       ! gems model
!       if (to_lower(config%mode) == "gems") then
!          call gems_header(output_unit)
!          block
!             type(gems_model), allocatable :: tmp
!             allocate (tmp)

!             ! todo: here one needs to pass also the info about the system, i.e. solvent, temperature, etc.
!             call new_gems_model(error, tmp, system, &
!                                  debug=config%debug, verbosity=config%verbosity, &
!                                  read_parameters=config%read_parameters, parameter_file=config%parameters_path &
!                                  )

!             call move_alloc(tmp, sm)
!          end block

!          ! rism1d model
!       else if (to_lower(config%mode) == "rism1d") then
! #ifdef WITH_RISM
!          block
!             type(rism1d_model), allocatable :: tmp
!             allocate (tmp)
!             call new_rism1d_model(error, tmp, &
!                & theory=trim(config%theory), closure=trim(config%closure), &
!                & solver=trim(config%solver), &
!                & verbosity=config%verbosity, system_parameters=system)
!             if (allocated(error)) return
!             call move_alloc(tmp, sm)
!          end block
! #else
!          call fatal_error(error, "RISM support not enabled at build time (-Drism=true)")
!          return
! #endif

!          ! rism3d model
!       else if (to_lower(config%mode) == "rism3d") then
! #ifdef WITH_RISM
!          block
!             type(rism3d_model), allocatable :: tmp
!             allocate (tmp)
!             call new_rism3d_model(error, tmp, &
!                & theory=trim(config%theory), closure=trim(config%closure), &
!                & solver=trim(config%solver), &
!                & verbosity=config%verbosity, system_parameters=system)
!             if (allocated(error)) return
!             call move_alloc(tmp, sm)
!          end block
! #else
!          call fatal_error(error, "RISM support not enabled at build time (-Drism=true)")
!          return
! #endif

!          ! alpb model
!       else if (to_lower(config%mode) == "alpb") then
!          call fatal_error(error, "ALPB not implemented yet")
!          return

!          ! else if(to_lower(config%mode) == "smd") then
!          !    block
!          !       type(smd_model), allocatable :: tmp
!          !       allocate(tmp)
!          !       call new_smd_model(error, tmp, system, debug=config%debug, verbosity=config%verbosity)
!          !       call move_alloc(tmp, sm)
!          !    end block

!       else
!          call fatal_error(error, "Unknown model selected")
!          return
!       end if

!       !> Initialize the solvation model
!       call sm%update(mol, error)
!       if (allocated(error)) return

!       !> Get the solvation energy
!       call sm%get_energy(wfn, energy, error)
!       if (allocated(error)) return

!       call timer%measure(1)

!       !> Print the solvation energy
!       if (config%verbosity > 0) write (output_unit, '(a,f20.10,a,f20.6,a)') &
!          "Final solvation free energy ", energy, &
!          " Eh ", energy*627.509, " kcal/mol", ""

!       if (config%writeenergy) then
!          !> Write the energy to a file
!          open (file=".GSOLV", newunit=unit, status='replace', action='write', iostat=stat)
!          if (stat /= 0) then
!             call fatal_error(error, "Could not open file for writing: "//trim(filename))
!             return
!          end if
!          write (unit, '(f20.12)') energy
!          close (unit)
!       end if

      call timer%write_timing(output_unit, 1, 'total', .true.)

   end subroutine run_main

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
