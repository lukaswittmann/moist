
!> Entry point for running single point calculations with moist
module moist_driver
   use moist_cavity_mc, only: moist_cavity_mc_parameters_type
   use, intrinsic :: iso_fortran_env, only: output_unit, input_unit
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type, read_structure, filetype
   use mctc_io_utils, only: to_lower
   use moist_cli, only: run_config
   use moist_app_rism_input, only: rism_input_type, new_rism_input
   use moist_app_rism_setup, only: new_rism_model
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   use moist_model_type, only: solvation_model_type
   use mctc_io_convert, only: autokcal
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_output_ascii, only: cavity_header
   use moist_data_solvents, only: get_solvent_id
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_mc, only: cavity_type_mc, new_cavity_mc
   use moist_cavity_mt, only: cavity_type_mt, new_cavity_mt
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_radii, only: radius_type, new_radii
   use moist_cavity_type, only: cavity_type
   use moist_context, only: moist_context_type, new_context
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
!>
!> @param[in] config  configuration for this driver
!> @param[out] error  error handling
   subroutine run_main(config, error)

      !> Configuration for this driver
      type(run_config), intent(in) :: config

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      !> Shared run context owned for the whole run borrowed by everything
      type(moist_context_type), target :: ctx

      !> Solvent system
      type(solvent_system_type), allocatable :: system
      !> TOML run settings, used only by the RISM commands
      type(rism_input_type) :: rism_input
      !> Whether the command is a RISM model (rism1d or rism3d)
      logical :: rism_mode

      character(len=:), allocatable :: filename
      character(len=:), allocatable :: solvent

      integer :: solvent_id

      !> File reading and parsing
      integer :: unit, stat
      character(len=256) :: tmp
      real(wp) :: tmp_wp

      !> Open the run context; every cavity/model borrows it, so all sub-timers
      !> started later nest under the "total" node opened here. A positive
      !> --threads is applied to the OpenMP runtime by the context.
      call new_context(ctx, nthreads=config%num_threads, verbosity=config%verbosity, debug=config%debug, &
         & do_profile=config%profile)
      call ctx%timer%start("total")

      !* ---------------------------- Thread configuration --------------------------- *!
      if (config%num_threads > 0) then
!$       if (.false.) then
            call ctx%warning("Program compiled without OpenMP support, ignoring --threads")
!$       else
! #ifdef WITH_MKL
! !$       call mkl_set_num_threads(config%num_threads)
! #endif
!$          if (ctx%verbosity(1)) then
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
         call new_solvent_system(system, solvent_id, &
                                              temperature=config%temperature, pressure_si=config%pressure_si, &
                                              error=error)
         if (allocated(error)) return
         if (ctx%verbosity(1)) call system%print()
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
               if (ctx%verbosity(1)) write (output_unit, "(a)") &
                  "[Info] Molecular charge read from '"//filename//"'"
            else
               call ctx%warning("Could not read molecular charge from '"//filename//"'")
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
          to_lower(config%mode) == "mt" .or. &
          to_lower(config%mode) == "iswig") then

         block
            class(cavity_type), allocatable :: cavity
            class(radius_type), allocatable :: radius_model

            ! Print header
            if (ctx%verbosity(1)) call cavity_header(output_unit, trim(config%mode))

            ! Instantiate the appropriate cavity type
            if (to_lower(config%mode) == "numsa") then
               block
                  type(cavity_type_numsa), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  call new_cavity_numsa(tmp_cavity, radii=radius_model, error=error, &
                     param=config%numsa, ctx=ctx)
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "iswig") then
               block
                  type(cavity_type_iswig), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  call new_cavity_iswig(tmp_cavity, radius_model=radius_model, error=error, &
                     param=config%iswig, ctx=ctx)
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
                        !> its own tolerance, so only the shape parameters go here
                        call svdw_template%new(param=config%svdw)
                        call new_cavity_drop(tmp_cavity, radius_model=radius_model, &
                           lsf_model=svdw_template, error=error, &
                           param=config%drop, ctx=ctx)
                     end block
                  else if (to_lower(config%drop_variant) == "cfc") then
                     block
                        type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
                        call cfc_template%new(param=config%cfc)
                        call new_cavity_drop(tmp_cavity, radius_model=radius_model, &
                           lsf_model=cfc_template, error=error, &
                           param=config%drop, ctx=ctx)
                     end block
                  else
                     call fatal_error(error, "Unknown DROP variant: "//trim(config%drop_variant))
                  end if
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "mc") then
               block
                  type(cavity_type_mc), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  !> Marching cubes writes its mesh while integrating, so the
                  !> export paths have to be known before update runs.
                  if (to_lower(config%drop_variant) == "svdw") then
                     block
                        type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
                        call svdw_template%new(param=config%svdw)
                        call new_mc_cavity(tmp_cavity, ctx, radius_model, svdw_template, &
                                           config%mc, config%dump, error)
                     end block
                  else if (to_lower(config%drop_variant) == "cfc") then
                     block
                        type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
                        call cfc_template%new(param=config%cfc)
                        call new_mc_cavity(tmp_cavity, ctx, radius_model, cfc_template, &
                                           config%mc, config%dump, error)
                     end block
                  else
                     call fatal_error(error, "Unknown level set function: "//trim(config%drop_variant))
                  end if
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            else if (to_lower(config%mode) == "mt") then
               block
                  type(cavity_type_mt), allocatable :: tmp_cavity
                  allocate (tmp_cavity)
                  call new_radii(config%radii, radius_model, error)
                  if (allocated(error)) return
                  if (to_lower(config%drop_variant) == "svdw") then
                     block
                        type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
                        call svdw_template%new(param=config%svdw)
                        call new_cavity_mt(tmp_cavity, ctx=ctx, radius_model=radius_model, &
                           lsf_model=svdw_template, param=config%mt, error=error)
                     end block
                  else if (to_lower(config%drop_variant) == "cfc") then
                     block
                        type(moist_cavity_drop_lsf_cfc_type) :: cfc_template
                        call cfc_template%new(param=config%cfc)
                        call new_cavity_mt(tmp_cavity, ctx=ctx, radius_model=radius_model, &
                           lsf_model=cfc_template, param=config%mt, error=error)
                     end block
                  else
                     call fatal_error(error, "Unknown level set function: "//trim(config%drop_variant))
                  end if
                  if (allocated(error)) return
                  call move_alloc(tmp_cavity, cavity)
               end block
            end if

            ! Settings of the constructed cavity
            if (ctx%verbosity(1)) call cavity%print_parameters()

            ! Use polymorphic cavity methods
            call cavity%update(mol, error=error)
            if (allocated(error)) return

            ! Radii of the structure; the update itself prints none
            if (ctx%verbosity(2)) call cavity%print_radii()

            ! Gradient if asked for
            if (config%grad) then
               call cavity%get_gradient(error)
               if (allocated(error)) return
            end if

            ! Print results; with no unit the cavity follows its own context
            if (ctx%verbosity(1)) call cavity%print()

            ! Write cavity files (xyz, csv, pqr) only when --dump is given.
            ! Marching cubes has no grid points to dump; it already wrote its
            ! triangle mesh (cavity.obj/cavity.pqr) during update.
            if (config%dump .and. to_lower(config%mode) /= "mc") then
               call cavity%write_xyz_debug("cavity.xyz", error=error)
               if (allocated(error)) return
               call cavity%write_csv_debug("cavity.csv", error=error)
               if (allocated(error)) return
               call cavity%write_pqr_debug("cavity.pqr", error=error)
               if (allocated(error)) return
            end if

         end block
         call ctx%timer%stop("total")
         call ctx%timer%print("moist")
         call ctx%print_citations()
         call ctx%delete()
         return
      end if

      ! Read native UV/VV settings plus app setup; CLI values override only when supplied
      rism_mode = to_lower(config%mode) == "rism1d" .or. to_lower(config%mode) == "rism3d"
      if (rism_mode) then
         call new_rism_input(rism_input, trim(config%mode), error)
         if (.not. allocated(error) .and. config%read_parameters) then
            call rism_input%read_file(trim(config%parameters_path), error, strict=.true.)
         end if
         if (allocated(error)) return
         if (allocated(config%temperature)) rism_input%system%temperature = config%temperature
         if (allocated(config%pressure_si)) rism_input%system%pressure = config%pressure_si
         if (allocated(rism_input%system%charge)) mol%charge = rism_input%system%charge
         if (allocated(config%charge)) mol%charge = config%charge
         rism_input%system%charge = mol%charge
      end if

      !* ---------------------------------- Get Solvent ---------------------------------- *!

      ! Get solvent
      if (len_trim(config%solvent) > 0) then
         solvent = trim(config%solvent)
      else if (rism_mode .and. len_trim(rism_input%system%solvent) > 0) then
         solvent = trim(rism_input%system%solvent)
      else
         filename = join(dirname(config%input), ".SOLVENT")
         if (exists(filename)) then
            open (file=filename, newunit=unit)
            read (unit, "(A)", iostat=stat) tmp
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

      if (rism_mode) then
         rism_input%system%solvent = solvent
         call run_rism(config, ctx, mol, rism_input, error)
         if (allocated(error)) return
         call ctx%timer%stop("total")
         call ctx%timer%print("moist")
         call ctx%delete()
         return
      end if

      call get_solvent_id(solvent, solvent_id, error)
      if (allocated(error)) return
      allocate (system)
      call new_solvent_system(system, solvent_id, &
                                           temperature=config%temperature, pressure_si=config%pressure_si, &
                                           error=error)
      if (allocated(error)) return

      !* ================================================================================= *!
      !*                                  Solvation models                                 *!
      !* ================================================================================= *!

      ! Reserved model commands without an implementation
      call fatal_error(error, "The selected model command is not implemented: "//trim(config%mode))
      return

   end subroutine run_main

   !> Run a standalone UV model through the domain-independent model interface
   !>
   !> @param[in] config    CLI output and gradient options
   !> @param[in] ctx       borrowed run context
   !> @param[in] mol       solute with the resolved total charge
   !> @param[inout] input  native model and application settings
   !> @param[out] error    setup, solve or output error
   subroutine run_rism(config, ctx, mol, input, error)
      !> CLI output and gradient options
      type(run_config), intent(in) :: config
      !> Borrowed run context
      type(moist_context_type), intent(in), target :: ctx
      !> Solute with the resolved total charge
      type(structure_type), intent(in) :: mol
      !> Native model and application settings
      type(rism_input_type), intent(inout) :: input
      !> Setup, solve or output error
      type(error_type), allocatable, intent(out) :: error
      class(solvation_model_type), allocatable, target :: model
      type(prettyprinter) :: pp
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: energy
      integer :: iat, unit, stat, close_stat

      ! Main settings at verbosity 2, every setting of the selected solvers from 3
      if (ctx%verbosity(2)) then
         pp = new_prettyprinter(unit=ctx%unit)
         call pp%section("Resolved RISM input (lengths in bohr, energies in Hartree)")
         call input%print_table(pp, ctx%verbosity())
         call pp%blank()
      end if
      call new_rism_model(input, ctx, mol, model, error)
      if (allocated(error)) return
      if (config%grad) allocate(gradient(3, mol%nat))
      if (config%grad .and. .not. config%numgrad) then
         call evaluate_rism_model(model, energy, error, gradient)
      else
         call evaluate_rism_model(model, energy, error)
         if (.not. allocated(error) .and. config%numgrad) then
            call numerical_rism_gradient(model, mol, gradient, error)
         end if
      end if
      if (allocated(error)) return
      if (ctx%verbosity(1)) then
         write(ctx%unit, "(a,f20.10,a,f16.6,a)") &
            & "Final solvation free energy ", energy, " Eh ", energy*autokcal, " kcal/mol"
         if (allocated(gradient)) then
            write(ctx%unit, "(a)") "Solvation gradient (Hartree/bohr):"
            do iat = 1, mol%nat
               write(ctx%unit, "(i6,3es22.12)") iat, gradient(:, iat)
            end do
         end if
      end if
      if (config%writeenergy) then
         open(newunit=unit, file=".GSOLV", status="replace", action="write", iostat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot open .GSOLV for writing")
            return
         end if
         write(unit, "(es24.16)", iostat=stat) energy
         close(unit, iostat=close_stat)
         if (stat /= 0 .or. close_stat /= 0) call fatal_error(error, "Cannot write solvation energy to .GSOLV")
      end if
   end subroutine run_rism

   !> Evaluate a parameter-only model through the shared coupling lifecycle
   !>
   !> @param[inout] model   updated solvation model
   !> @param[out] energy    solvation free energy, Hartree
   !> @param[out] error     coupling, energy or gradient failure
   !> @param[out] gradient  optional nuclear gradient (3, nat), Hartree/bohr
   subroutine evaluate_rism_model(model, energy, error, gradient)
      !> Updated solvation model
      class(solvation_model_type), intent(inout), target :: model
      !> Solvation free energy, Hartree
      real(wp), intent(out) :: energy
      !> Coupling, energy or gradient failure
      type(error_type), allocatable, intent(out) :: error
      !> Optional nuclear gradient (3, nat), Hartree/bohr
      real(wp), intent(out), optional :: gradient(:, :)
      type(coupling_type), pointer :: coupling
      type(response_type) :: response

      energy = 0.0_wp
      call model%new_coupling(coupling, error)
      if (allocated(error)) return
      call model%prepare_energy(coupling, error)
      if (.not. allocated(error)) call model%get_energy(coupling, energy, error)
      if (.not. allocated(error) .and. present(gradient)) then
         gradient = 0.0_wp
         call model%prepare_gradient(coupling, error)
         if (.not. allocated(error)) call model%get_gradient(coupling, response, gradient, error)
      end if
      call model%release_coupling(coupling)
   end subroutine evaluate_rism_model

   !> Central-difference CLI gradient with fixed charges, bulk and grid settings
   !>
   !> @param[inout] model   constructed model, restored to mol before returning
   !> @param[in] mol        original solute
   !> @param[out] gradient  nuclear gradient (3, nat), Hartree/bohr
   !> @param[out] error     failed displaced or restored solve
   subroutine numerical_rism_gradient(model, mol, gradient, error)
      !> Constructed model, restored to mol before returning
      class(solvation_model_type), intent(inout), target :: model
      !> Original solute
      type(structure_type), intent(in) :: mol
      !> Nuclear gradient (3, nat), Hartree/bohr
      real(wp), intent(out) :: gradient(:, :)
      !> Failed displaced or restored solve
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: trial
      type(error_type), allocatable :: restore_error
      real(wp), parameter :: step = 1.0e-4_wp
      real(wp) :: plus, minus, restored
      integer :: iat, axis

      if (any(shape(gradient) /= [3, mol%nat])) then
         call fatal_error(error, "Numerical RISM gradient shape mismatch")
         return
      end if
      gradient = 0.0_wp
      atoms: do iat = 1, mol%nat
         do axis = 1, 3
            trial = mol
            trial%xyz(axis, iat) = mol%xyz(axis, iat) + step
            call model%update(trial, error)
            if (.not. allocated(error)) call evaluate_rism_model(model, plus, error)
            if (allocated(error)) exit atoms
            trial%xyz(axis, iat) = mol%xyz(axis, iat) - step
            call model%update(trial, error)
            if (.not. allocated(error)) call evaluate_rism_model(model, minus, error)
            if (allocated(error)) exit atoms
            gradient(axis, iat) = (plus - minus)/(2.0_wp*step)
         end do
      end do atoms
      call model%update(mol, restore_error)
      if (.not. allocated(error) .and. .not. allocated(restore_error)) then
         call evaluate_rism_model(model, restored, restore_error)
      end if
      if (.not. allocated(error)) call move_alloc(restore_error, error)
   end subroutine numerical_rism_gradient

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
   !> @param[in]    param         Marching-cubes settings
   !> @param[in]    dump          Write the triangle mesh to cavity.obj/cavity.pqr
   !> @param[out]   error         Error handling
   subroutine new_mc_cavity(cavity, ctx, radius_model, lsf_model, param, dump, error)
      type(cavity_type_mc), intent(inout) :: cavity
      type(moist_context_type), intent(in), target :: ctx
      class(radius_type), intent(in) :: radius_model
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      type(moist_cavity_mc_parameters_type), intent(in) :: param
      logical, intent(in) :: dump
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_mc_parameters_type) :: settings

      settings = param
      if (dump) then
         settings%obj_file = "cavity.obj"
         settings%pqr_file = "cavity.pqr"
      end if
      call new_cavity_mc(cavity, radius_model=radius_model, lsf_model=lsf_model, &
         error=error, param=settings, ctx=ctx)

   end subroutine new_mc_cavity

!> Construct path by joining strings with os file separator
   function join(a1, a2) result(path)
      use mctc_env_system, only: is_windows
      character(len=*), intent(in) :: a1, a2
      character(len=:), allocatable :: path
      character :: filesep

      if (is_windows()) then
         filesep = "\"
      else
         filesep = "/"
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
