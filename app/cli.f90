

!> fclap-based command line interface for moist
module moist_cli
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_convert, only: aatoau
   use mctc_io, only : get_filetype
   use moist, only : get_moist_version
   use moist_output_ascii, only : moist_header, moist_build_header
   use moist_output_citations, only : print_citations
   use moist_output_license, only : print_license
   use fclap, only : ArgumentParser, Namespace, not_less_than
   implicit none (type, external)
   private

   public :: run_config, get_arguments

   abstract interface
      !> Callback signature for adding arguments to a subsubparser
      subroutine subsubparser_args_callback(p)
         import :: ArgumentParser
         type(ArgumentParser), intent(inout) :: p
      end subroutine subsubparser_args_callback
   end interface

   !> Declarative configuration for a subsubparser entry
   type :: subsubparser_spec
      character(len=16) :: name = ''
      character(len=64) :: help_text = ''
      character(len=64) :: prog = ''
      character(len=128) :: description = ''
      procedure(subsubparser_args_callback), pointer, nopass :: add_specific_args => null()
   contains
      procedure :: init => subsubparser_spec_init
   end type subsubparser_spec


   !> Configuration data for running stand-alone calculations
   type :: run_config

      character(256) :: input = ''
      integer, allocatable :: input_format

      character(32) :: mode = ''

      ! System info
      real(wp), allocatable :: charge
      real(wp) :: temperature = 298.15_wp
      real(wp) :: pressure_si = 101325.0_wp
      character(64) :: solvent = ''

      ! RISM model selection (rism1d/rism3d subcommands)
      character(32) :: closure = 'KH'
      character(32) :: theory = 'DRISM'
      character(32) :: solver = 'gmres'

      logical :: json = .false.

      logical :: grad = .false.
      logical :: numgrad = .false.

      logical :: read_parameters = .false.
      character(256) :: parameters_path = ''

      logical :: writeenergy = .false.

      integer :: verbosity = 2
      logical :: debug = .false.

      !> Number of OpenMP threads (0 = use default)
      integer :: num_threads = 0

      !> Cavity discretization (Lebedev points per sphere)
      integer :: nleb = 194

      !> Radii set selection (cpcm, smd, d3, cosmo, bondi)
      character(32) :: radii = 'cpcm'

      !> DROP cavity settings
      character(32) :: drop_variant = ''
      real(wp) :: drop_tol = 1.0E-10_wp
      integer :: drop_proj_level = 3
      integer :: drop_wleb_prune_level = 0
      real(wp) :: drop_blend_k = 5.5_wp
      real(wp) :: drop_blend_1b = 1.0_wp
      real(wp) :: drop_blend_2b = 0.0_wp
      real(wp) :: drop_blend_3b = 3.0_wp
      real(wp) :: cfc_a1 = -15.0_wp
      real(wp) :: cfc_a2 = -9.0_wp
      real(wp) :: cfc_c = 5.0_wp
      integer :: cfc_m = 4
      real(wp) :: cfc_screen_k = 3.0_wp

      !> Enable all optional cavity properties (curvature, normals, etc.)
      logical :: cavity_fine = .true.
      !> Enable marching-cubes reference area/volume computation
      logical :: cavity_mc = .false.
      !> Write cavity files (xyz, pqr, csv) to disk
      logical :: dump = .false.

   end type run_config

contains


!> Parse command line arguments and build run configuration
!> @param[out] config Parsed runtime configuration for the driver
!> @param[out] error  Allocated on command line validation failures
subroutine get_arguments(config, error)

   !> Configuation data
   type(run_config), intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(ArgumentParser), save :: parser
   type(ArgumentParser), save :: general_parent
   type(ArgumentParser), save :: model_parser
   type(ArgumentParser), save :: model_subparsers(4)
   type(ArgumentParser), save :: cavity_parser
   type(ArgumentParser), save :: cavity_subparsers(4)
   type(ArgumentParser), save :: drop_parser
   type(ArgumentParser), save :: drop_subparsers(2)
   type(ArgumentParser), save :: solvent_parser
   type(Namespace) :: args
   type(subsubparser_spec) :: model_specs(4)
   type(subsubparser_spec) :: cavity_specs(3)

   character(len=:), allocatable :: version_string
   character(len=32) :: command
   character(len=256) :: str_tmp
   integer :: verbose_count, quiet_count
   logical :: show_citation, show_license

   config = run_config()

   call moist_header(output_unit)
   call moist_build_header(output_unit)

   call get_moist_version(string=version_string)

   call parser%init( &
      prog='moist', &
      description='Modular and Open-source Implicit Solvation Toolkit', &
      epilog='For more information, see the documentation at https://github.com/lukaswittmann/moist', &
      version='v'//version_string)

   call parser%add_argument('--citation', action='store_true', dest='citation', &
      help='Print citation information and exit')
   call parser%add_argument('--license', action='store_true', dest='license', &
      help='Print full license information and exit')

   call parser%add_subparsers(title='subcommands', dest='command')

   call model_specs(1)%init('gems', 'GEMS model', 'moist model gems', &
      'Run the GEMS solvation model')
   call model_specs(2)%init('alpb', 'ALPB model', 'moist model alpb', &
      'Run the ALPB solvation model')
   call model_specs(3)%init('rism1d', 'RISM1D model', 'moist model rism1d', &
      'Run the RISM1D solvation model')
   call model_specs(4)%init('rism3d', 'RISM3D model', 'moist model rism3d', &
      'Run the RISM3D solvation model')
   model_specs(3)%add_specific_args => add_model_rism_arguments
   model_specs(4)%add_specific_args => add_model_rism_arguments

   call cavity_specs(1)%init('numsa', 'NUMSA cavity', 'moist cavity numsa', &
      'Construct NUMSA cavities')
   call cavity_specs(2)%init('iswig', 'ISWIG cavity', 'moist cavity iswig', &
      'Construct ISWIG cavities')
   call cavity_specs(3)%init('drop', 'DROP cavity', 'moist cavity drop', &
      'Construct DROP/DROP cavities')

   call setup_general_parent(general_parent)
   call setup_model_parser(model_parser, general_parent, model_specs, model_subparsers)
   call setup_cavity_parser(cavity_parser, general_parent, cavity_specs, cavity_subparsers, &
      drop_parser, drop_subparsers)
   call setup_solvent_parser(solvent_parser, general_parent)

   call parser%add_parser('model', model_parser, &
      help_text='Run a full solvation model (gems, alpb, rism1d, rism3d)')
   call parser%add_parser('cavity', cavity_parser, &
      help_text='Construct a cavity (numsa, iswig, drop)')
   call parser%add_parser('solvent', solvent_parser, &
      help_text='Inspect solvent properties by name or alias')

   args = parser%parse_args()

   show_citation = .false.
   show_license = .false.
   if (args%has_key('citation')) call args%get('citation', show_citation)
   if (args%has_key('license')) call args%get('license', show_license)

   if (show_citation) then
      call print_citations(output_unit)
      stop
   end if

   if (show_license) then
      call print_license(output_unit)
      stop
   end if

   if (.not. args%has_key('command')) then
      call fatal_error(error, 'Select one of: model, cavity, solvent subcommands')
      return
   end if

   call args%get('command', command)

   verbose_count = 0
   quiet_count = 0
   if (args%has_key('verbose')) call args%get('verbose', verbose_count)
   if (args%has_key('quiet')) call args%get('quiet', quiet_count)
   if (args%has_key('debug')) call args%get('debug', config%debug)
   if (args%has_key('threads')) call args%get('threads', config%num_threads)
   config%verbosity = config%verbosity + verbose_count - quiet_count

   select case(trim(command))
   case('model')
      if (.not. args%has_key('model_mode')) then
         call fatal_error(error, 'Select one of: gems, alpb, rism1d, rism3d models')
         return
      end if

      call args%get('model_mode', config%mode)
      call args%get('input', config%input)

      if (args%has_key('closure')) then
         call args%get('closure', config%closure)
         config%closure = trim(adjustl(config%closure))
      end if

      if (args%has_key('theory')) then
         call args%get('theory', config%theory)
         config%theory = trim(adjustl(config%theory))
      end if

      if (args%has_key('solver')) then
         call args%get('solver', config%solver)
         config%solver = trim(adjustl(config%solver))
      end if

      if (args%has_key('charge')) then
         allocate(config%charge)
         call args%get('charge', config%charge)
      end if

      if (args%has_key('solvent')) then
         call args%get('solvent', config%solvent)
         config%solvent = trim(adjustl(config%solvent))
      end if

      if (args%has_key('temperature')) then
         call args%get('temperature', config%temperature)
      end if

      if (args%has_key('pressure')) then
         call args%get('pressure', config%pressure_si)
      end if

      if (args%has_key('numgrad')) then
         call args%get('numgrad', config%numgrad)
      end if
      if (config%numgrad) config%grad = .true.

      if (args%has_key('writeenergy')) then
         call args%get('writeenergy', config%writeenergy)
      end if

      if (args%has_key('parameters_path')) then
         call args%get('parameters_path', str_tmp)
         config%parameters_path = trim(adjustl(str_tmp))
         if (len_trim(config%parameters_path) > 0) config%read_parameters = .true.
      end if

   case('cavity')
      if (.not. args%has_key('cavity_mode')) then
         call fatal_error(error, 'Select one of: numsa, iswig, drop cavity types')
         return
      end if

      call args%get('cavity_mode', config%mode)
      if (trim(config%mode) == 'drop') then
         if (.not. args%has_key('drop_variant')) then
            call fatal_error(error, 'Select one of: svdw, cfc DROP variants')
            return
         end if
         call args%get('drop_variant', config%drop_variant)
      end if
      call args%get('input', config%input)

      if (args%has_key('nleb')) then
         call args%get('nleb', config%nleb)
      end if

      if (args%has_key('radii')) then
         call args%get('radii', config%radii)
      end if
      config%radii = trim(adjustl(config%radii))

      if (args%has_key('drop_tolerance')) then
         call args%get('drop_tolerance', config%drop_tol)
      end if

      if (args%has_key('drop_proj_level')) then
         call args%get('drop_proj_level', config%drop_proj_level)
      end if

      if (args%has_key('drop_wleb_prune_level')) then
         call args%get('drop_wleb_prune_level', config%drop_wleb_prune_level)
      end if

      if (args%has_key('drop_blendk')) then
         call args%get('drop_blendk', config%drop_blend_k)
      end if

      if (args%has_key('drop_blend1b')) then
         call args%get('drop_blend1b', config%drop_blend_1b)
      end if

      if (args%has_key('drop_blend2b')) then
         call args%get('drop_blend2b', config%drop_blend_2b)
      end if

      if (args%has_key('drop_blend3b')) then
         call args%get('drop_blend3b', config%drop_blend_3b)
      end if

      if (args%has_key('cfc_a1')) then
         call args%get('cfc_a1', config%cfc_a1)
      end if

      if (args%has_key('cfc_a2')) then
         call args%get('cfc_a2', config%cfc_a2)
      end if

      if (args%has_key('cfc_c')) then
         call args%get('cfc_c', config%cfc_c)
      end if

      if (args%has_key('cfc_m')) then
         call args%get('cfc_m', config%cfc_m)
      end if

      if (args%has_key('cfc_screen_k')) then
         call args%get('cfc_screen_k', config%cfc_screen_k)
      end if

      if (args%has_key('grad')) then
         call args%get('grad', config%grad)
      end if

      if (args%has_key('nofine')) then
         call args%get('nofine', config%cavity_fine)
      end if
      if (args%has_key('mc')) then
         call args%get('mc', config%cavity_mc)
      end if
      if (args%has_key('dump')) then
         call args%get('dump', config%dump)
      end if

   case('solvent')
      call args%get('solvent', config%solvent)
      config%solvent = trim(adjustl(config%solvent))
      config%mode = 'solvent'
      config%input = ''

   case default
      call fatal_error(error, 'Unknown subcommand selected')
      return
   end select

   config%mode = trim(config%mode)
   config%drop_variant = trim(config%drop_variant)
   config%input = trim(config%input)

   if (trim(command) /= 'solvent') then
      config%input_format = get_filetype(trim(config%input))
   end if

contains

   !> Configure shared parent parser with global runtime options
   !> @param[inout] parent Parent parser receiving shared options
   subroutine setup_general_parent(parent)
      type(ArgumentParser), intent(inout) :: parent
      integer :: grp_general

      call parent%init(add_help=.false.)

      grp_general = parent%add_argument_group('General settings', &
         'Shared runtime configuration')

      call parent%add_argument('-v', '--verbose', action='count', dest='verbose', &
         help='Increase output verbosity (repeatable)', &
         group_idx=grp_general)

      call parent%add_argument('-q', '--quiet', action='count', dest='quiet', &
         help='Decrease output verbosity (repeatable)', &
         group_idx=grp_general)

      call parent%add_argument('-d', '--debug', action='store_true', dest='debug', &
         help='Enable debug output', &
         group_idx=grp_general)

      call parent%add_argument('-t', '--threads', data_type='integer', &
         action=not_less_than(0), &
         dest='threads', metavar='INT', &
         help='Number of OpenMP threads (default: system setting)', &
         group_idx=grp_general)
   end subroutine setup_general_parent


   !> Build and register a set of subsubparsers from declarative metadata
   !> @param[inout] p               Parent parser receiving subsubparsers
   !> @param[in]    parent          Parser providing shared parent arguments
   !> @param[in]    subparser_title Title for the subsubparser help section
   !> @param[in]    subparser_dest  Namespace key receiving selected subsubparser name
   !> @param[in]    specs           Metadata table describing subsubparsers
   !> @param[inout] subparsers      Storage for subsubparser parser instances
   !> @param[in]    add_shared_args Callback adding arguments common to all entries
   subroutine setup_subsubparsers(p, parent, subparser_title, subparser_dest, specs, subparsers, add_shared_args)
      type(ArgumentParser), intent(inout) :: p
      type(ArgumentParser), intent(in) :: parent
      character(len=*), intent(in) :: subparser_title, subparser_dest
      type(subsubparser_spec), intent(in) :: specs(:)
      type(ArgumentParser), intent(inout) :: subparsers(:)
      procedure(subsubparser_args_callback) :: add_shared_args
      integer :: i, n

      call p%add_subparsers(title=trim(subparser_title), dest=trim(subparser_dest))

      n = min(size(specs), size(subparsers))
      do i = 1, n
         call subparsers(i)%init_with_parents([parent], &
            prog=trim(specs(i)%prog), &
            description=trim(specs(i)%description))

         call add_shared_args(subparsers(i))
         if (associated(specs(i)%add_specific_args)) then
            call specs(i)%add_specific_args(subparsers(i))
         end if

         call p%add_parser(trim(specs(i)%name), subparsers(i), help_text=trim(specs(i)%help_text))
      end do
   end subroutine setup_subsubparsers


   !> Configure the model subcommand parser
   !> @param[inout] p          Model subcommand parser
   !> @param[in]    parent     Parent parser providing shared options
   !> @param[in]    specs      Model subsubparser metadata
   !> @param[inout] subparsers Model subsubparser instances
   subroutine setup_model_parser(p, parent, specs, subparsers)
      type(ArgumentParser), intent(inout) :: p
      type(ArgumentParser), intent(in) :: parent
      type(subsubparser_spec), intent(in) :: specs(:)
      type(ArgumentParser), intent(inout) :: subparsers(:)

      call p%init_with_parents([parent], &
         prog='moist model', &
         description='Run solvation models with model-specific subcommands')

      call setup_subsubparsers(p, parent, 'models', 'model_mode', specs, subparsers, add_model_shared_arguments)
   end subroutine setup_model_parser


   !> Add model arguments shared by all model type subcommands
   !> @param[inout] p Model type parser receiving shared options
   subroutine add_model_shared_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_system, grp_io, grp_advanced

      grp_system = p%add_argument_group('System config', &
         'Molecule and solvent thermodynamic settings')
      grp_io = p%add_argument_group('Input/Output', &
         'Input model and output behavior')
      grp_advanced = p%add_argument_group('Advanced', &
         'Advanced and developer options')
      call p%add_argument('input', &
         help='Input structure file', metavar='INPUT', &
         group_idx=grp_io)

      call p%add_argument('-c', '--charge', &
         data_type='real', &
         metavar='REAL', &
         help='Molecular charge, overwrites .CHRG file', &
         group_idx=grp_system)
      call p%add_argument('-s', '--solvent', &
         default_val='', &
         print_default=.false., &
         metavar='SOLVENT', &
         help='Solvent name', &
         group_idx=grp_system)
      call p%add_argument('--temperature', &
         default_val=298.15_wp, &
         action=not_less_than(0.0_wp), &
         data_type='real', &
         metavar='REAL', &
         help='Temperature in Kelvin', &
         group_idx=grp_system)
      call p%add_argument('--pressure', &
         default_val=101325.0_wp, &
         action=not_less_than(0.0_wp), &
         data_type='real', dest='pressure', &
         metavar='REAL', &
         help='Pressure in Pascal', &
         group_idx=grp_system)

      call p%add_argument('-w', '--writeenergy', action='store_true', &
         help='Write energy to .GSOLV file', &
         group_idx=grp_io)
      call p%add_argument('-p', '--parameters', &
         default_val='', &
         print_default=.false., &
         dest='parameters_path', &
         metavar='PATH', &
         help='Path to custom model parameters file', &
         group_idx=grp_io)

      call p%add_argument('--numgrad', action='store_true', &
         help='Calculate numerical gradients', &
         group_idx=grp_advanced)
   end subroutine add_model_shared_arguments


   !> Add RISM-only arguments (theory, closure, solver) to a model subsubparser
   !> @param[inout] p RISM model parser (rism1d or rism3d)
   subroutine add_model_rism_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_rism

      grp_rism = p%add_argument_group('RISM settings', &
         'RISM theory, closure relation, and iterative solver')

      call p%add_argument('--closure', default_val='KH', &
         print_choices=.true., &
         choices=[character(len=4) :: 'HNC', 'KH', 'PY', &
                  'PSE1', 'PSE2', 'PSE3', 'PSE4'], &
         metavar='CLOSURE', &
         help='RISM closure relation', &
         group_idx=grp_rism)
      call p%add_argument('--theory', default_val='DRISM', &
         print_choices=.true., &
         choices=[character(len=5) :: 'DRISM', 'XRISM'], &
         metavar='THEORY', &
         help='RISM theory variant', &
         group_idx=grp_rism)
      call p%add_argument('--solver', default_val='gmres', &
         print_choices=.true., &
         choices=[character(len=6) :: 'picard', 'mdiis', 'gmres', 'hybrid', 'lbfgs'], &
         metavar='SOLVER', &
         help='RISM iterative solver', &
         group_idx=grp_rism)
   end subroutine add_model_rism_arguments


   !> Configure the cavity subcommand parser
   !> @param[inout] p               Cavity subcommand parser
   !> @param[in]    parent          Parent parser providing shared options
   !> @param[in]    specs           Cavity subsubparser metadata
   !> @param[inout] subparsers      Cavity subsubparser instances
   !> @param[inout] drop_parser     Nested DROP parser
   !> @param[inout] drop_subparsers Nested DROP variant parsers
   subroutine setup_cavity_parser(p, parent, specs, subparsers, drop_parser, drop_subparsers)
      type(ArgumentParser), intent(inout) :: p
      type(ArgumentParser), intent(in) :: parent
      type(subsubparser_spec), intent(in) :: specs(:)
      type(ArgumentParser), intent(inout) :: subparsers(:)
      type(ArgumentParser), intent(inout) :: drop_parser
      type(ArgumentParser), intent(inout) :: drop_subparsers(:)
      integer :: i, n

      call p%init_with_parents([parent], &
         prog='moist cavity', &
         description='Run cavity-only workflows with cavity-type subcommands')

      call p%add_subparsers(title='cavity types', dest='cavity_mode')

      n = min(size(specs), size(subparsers))
      do i = 1, n
         if (trim(specs(i)%name) == 'drop') then
            call setup_drop_parser(drop_parser, parent, drop_subparsers)
            call p%add_parser('drop', drop_parser, help_text=trim(specs(i)%help_text))
         else
            call subparsers(i)%init_with_parents([parent], &
               prog=trim(specs(i)%prog), &
               description=trim(specs(i)%description))

            call add_cavity_shared_arguments(subparsers(i))
            if (associated(specs(i)%add_specific_args)) then
               call specs(i)%add_specific_args(subparsers(i))
            end if

            call p%add_parser(trim(specs(i)%name), subparsers(i), help_text=trim(specs(i)%help_text))
         end if
      end do
   end subroutine setup_cavity_parser


   !> Configure the nested DROP cavity parser
   !> @param[inout] p          DROP parser receiving nested variants
   !> @param[in]    parent     Parent parser providing shared options
   !> @param[inout] subparsers DROP variant parser instances
   subroutine setup_drop_parser(p, parent, subparsers)
      type(ArgumentParser), intent(inout) :: p
      type(ArgumentParser), intent(in) :: parent
      type(ArgumentParser), intent(inout) :: subparsers(:)

      call p%init_with_parents([parent], &
         prog='moist cavity drop', &
         description='Construct DROP/DROP cavities with variant-specific level sets')

      call p%add_subparsers(title='DROP variants', dest='drop_variant')

      call subparsers(1)%init_with_parents([parent], &
         prog='moist cavity drop svdw', &
         description='Construct DROP/DROP cavities with the SvdW level set')
      call add_cavity_shared_arguments(subparsers(1))
      call add_cavity_drop_arguments(subparsers(1))
      call add_cavity_drop_svdw_arguments(subparsers(1))
      call p%add_parser('svdw', subparsers(1), help_text='SvdW DROP cavity')

      call subparsers(2)%init_with_parents([parent], &
         prog='moist cavity drop cfc', &
         description='Construct DROP/DROP cavities with the CFC level set')
      call add_cavity_shared_arguments(subparsers(2))
      call add_cavity_drop_arguments(subparsers(2))
      call add_cavity_drop_cfc_arguments(subparsers(2))
      call p%add_parser('cfc', subparsers(2), help_text='CFC DROP cavity')
   end subroutine setup_drop_parser


   !> Add cavity arguments shared by all cavity type subcommands
   !> @param[inout] p Cavity type parser receiving shared options
   subroutine add_cavity_shared_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_input, grp_technical

      grp_input = p%add_argument_group('Input/Output', &
         'Input cavity model and structure')
      grp_technical = p%add_argument_group('Technical settings', &
         'Cavity discretization and numerical control')

      call p%add_argument('input', &
         help='Input coordinate file', metavar='COORD', &
         group_idx=grp_input)
      call p%add_argument('--nleb', data_type='integer', &
         action=not_less_than(1), &
         default_val=194, &
         metavar='INT', &
         help='Lebedev grid points per atom', &
         group_idx=grp_technical)
      call p%add_argument('--radii', default_val='cpcm', &
         print_choices=.true., &
         choices=[character(len=5) :: 'cpcm', 'smd', 'd3', 'cosmo', 'bondi'], &
         metavar='RADII', &
         help='Atomic radii', &
         group_idx=grp_technical)
      call p%add_argument('-g', '--grad', action='store_true', &
         help='Calculate cavity gradients', &
         group_idx=grp_technical)
      call p%add_argument('--dump', action='store_true', &
         dest='dump', &
         help='Write cavity files (xyz, pqr, csv) to disk', &
         group_idx=grp_input)
   end subroutine add_cavity_shared_arguments


   !> Add DROP/DROP arguments shared by all DROP variants
   !> @param[inout] p DROP variant parser
   subroutine add_cavity_drop_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_technical

      grp_technical = p%add_argument_group('DROP settings', &
         'DROP projection and property controls')

      call p%add_argument('--tol', data_type='real', &
         default_val=1.0E-10_wp, &
         action=not_less_than(0.0_wp), &
         dest='drop_tolerance', metavar='REAL', &
         help='Numerical tolerance', &
         group_idx=grp_technical)
      call p%add_argument('--proj-level', data_type='integer', &
         default_val=2, &
         action=not_less_than(1), &
         dest='drop_proj_level', metavar='INT', &
         help='Projection solver level', &
         choices=[character(len=26) :: "1=SLSQP", "2=SLSQP+Newton", &
            "3=Cond. multi-tangent", "4=Cond. SLSQP-deflation", &
            "5=SLSQP-deflation", "6=Newton-deflation", &
            "7=Multistart", "8=Fine multistart"], &
         print_choices=.true., &
         group_idx=grp_technical)
      call p%add_argument('--wleb-switch', data_type='integer', &
         default_val=0, &
         action=not_less_than(0), &
         dest='drop_wleb_prune_level', metavar='INT', &
         help='Smooth weight switching level (0=off, 1-4=increasing)', &
         choices=[character(len=26) :: "0=off", "1=1E-12/1E-10", &
            "2=1E-10/1E-8", "3=1E-8/1E-6", "4=1E-6/1E-4"], &
         print_choices=.true., &
         group_idx=grp_technical)
      call p%add_argument('--nofine', action='store_false', &
         dest='nofine', &
         help='Skip optional cavity properties (curvature, normals, etc.)', &
         group_idx=grp_technical)
      call p%add_argument('--mc', action='store_true', &
         dest='mc', &
         help='Compute marching-cubes reference area and volume', &
         group_idx=grp_technical)
   end subroutine add_cavity_drop_arguments


   !> Add SvdW DROP arguments to a DROP variant parser
   !> @param[inout] p SvdW DROP parser
   subroutine add_cavity_drop_svdw_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_technical

      grp_technical = p%add_argument_group('SvdW settings', &
         'SvdW level-set blending controls')

      call p%add_argument('--blendk', data_type='real', &
         default_val=5.5_wp, &
         action=not_less_than(0.0_wp), &
         dest='drop_blendk', metavar='REAL', &
         help='DROP blending sharpness k', &
         group_idx=grp_technical)
      call p%add_argument('--blend1b', data_type='real', &
         default_val=1.0_wp, &
         dest='drop_blend1b', metavar='REAL', &
         help='DROP one-body contribution', &
         group_idx=grp_technical)
      call p%add_argument('--blend2b', data_type='real', &
         default_val=0.0_wp, &
         dest='drop_blend2b', metavar='REAL', &
         help='DROP two-body contribution', &
         group_idx=grp_technical)
      call p%add_argument('--blend3b', data_type='real', &
         default_val=3.0_wp, &
         dest='drop_blend3b', metavar='REAL', &
         help='DROP three-body contribution', &
         group_idx=grp_technical)
   end subroutine add_cavity_drop_svdw_arguments


   !> Add CFC DROP arguments to a DROP variant parser
   !> @param[inout] p CFC DROP parser
   subroutine add_cavity_drop_cfc_arguments(p)
      type(ArgumentParser), intent(inout) :: p
      integer :: grp_technical

      grp_technical = p%add_argument_group('CFC settings', &
         'CFC level-set controls')

      call p%add_argument('--a1', data_type='real', &
         default_val=-15.0_wp, &
         dest='cfc_a1', metavar='REAL', &
         help='CFC atomic-term exponent', &
         group_idx=grp_technical)
      call p%add_argument('--a2', data_type='real', &
         default_val=-9.0_wp, &
         dest='cfc_a2', metavar='REAL', &
         help='CFC pair-term exponent', &
         group_idx=grp_technical)
      call p%add_argument('--c', data_type='real', &
         default_val=5.0_wp, &
         dest='cfc_c', metavar='REAL', &
         help='CFC pair-term coupling', &
         group_idx=grp_technical)
      call p%add_argument('--m', data_type='integer', &
         default_val=4, &
         action=not_less_than(1), &
         dest='cfc_m', metavar='INT', &
         help='CFC pair-term power', &
         group_idx=grp_technical)
      call p%add_argument('--screen-k', data_type='real', &
         default_val=3.0_wp, &
         action=not_less_than(0.0_wp), &
         dest='cfc_screen_k', metavar='REAL', &
         help='CFC SSD screening sharpness k', &
         group_idx=grp_technical)
   end subroutine add_cavity_drop_cfc_arguments


   !> Configure the solvent subcommand parser
   !> @param[inout] p      Solvent subcommand parser
   !> @param[in]    parent Parent parser providing shared options
   subroutine setup_solvent_parser(p, parent)
      type(ArgumentParser), intent(inout) :: p
      type(ArgumentParser), intent(in) :: parent
      integer :: grp_system

      call p%init_with_parents([parent], &
         prog='moist solvent', &
         description='Inspect solvent parameters by name/alias')

      grp_system = p%add_argument_group('System config', &
         'Solvent query settings')

      call p%add_argument('solvent', help='Solvent name or alias', metavar='SOLVENT', &
         group_idx=grp_system)
   end subroutine setup_solvent_parser


end subroutine get_arguments


!> Initialize one subsubparser specification entry
!> @param[inout] self        Subsubparser metadata entry
!> @param[in]    name        Subparser name
!> @param[in]    help_text   Short help text for command listing
!> @param[in]    prog        Full command name used in usage output
!> @param[in]    description Description text for the subparser help page
subroutine subsubparser_spec_init(self, name, help_text, prog, description)
   class(subsubparser_spec), intent(inout) :: self
   character(len=*), intent(in) :: name, help_text, prog, description

   self%name = ''
   self%help_text = ''
   self%prog = ''
   self%description = ''
   nullify(self%add_specific_args)

   self%name = trim(name)
   self%help_text = trim(help_text)
   self%prog = trim(prog)
   self%description = trim(description)
end subroutine subsubparser_spec_init


end module moist_cli
