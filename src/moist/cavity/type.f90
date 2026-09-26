!> Abstract cavity: the surface grid every discretization shares
module moist_cavity_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   use moist_context, only: moist_context_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_response, only: response_type
   use moist_channels_coupling, only: coupling_type
   use moist_cavity_fields, only: cavity_field_query_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none(type, external)
   private

   public :: cavity_type, list_cavity_fields_base
   public :: write_cavity_xyz_debug, write_cavity_csv_debug, write_cavity_pqr_debug

   !> Abstract base type containing minimal cavity/surface information
   !>
   !> Cavities within moist are per default discretized using Gaussians
   type, abstract :: cavity_type
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller, never allocated or freed by the cavity
      type(moist_context_type), pointer :: ctx => null()

      !> Sphere radii, bohr (nat)
      real(wp), allocatable :: radii(:)
      !> Radii model used to update cached radii
      class(radius_type), allocatable :: radius_model

      !> Number of atomic spheres
      integer :: nsph = 0

      !> Area per sphere (nsph)
      real(wp), allocatable :: asph(:)

      !> Number of cavity points
      integer :: ngrid = 0

      !> Cartesian coordinates of atomic sphere centers (3, nsph)
      real(wp), allocatable :: sphxyz(:, :)

      !> Owner of each grid point (ngrid)
      integer, allocatable :: owner(:)

      !> Cartesian coords of points (3,ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Point area, bohr**2 (ngrid)
      real(wp), allocatable :: a(:)
      !> Gaussian switching factor of each surface point (ngrid)
      real(wp), allocatable :: f(:)
      !> Gaussian width of each surface point (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Outward unit normal of each surface point (3, ngrid)
      real(wp), allocatable :: normal0(:, :)
      !> Point volume element, bohr**3 (ngrid)
      !>
      !> Divergence-theorem partition of the enclosed volume,
      !> v_i = a_i (r_i . n_i)/3, so that `total_volume` is `sum(v)`
      real(wp), allocatable :: v(:)

      !> Total surface area, bohr**2
      real(wp), allocatable :: total_area
      !> Total cavity volume, bohr**3
      real(wp), allocatable :: total_volume

      !> Nuclear derivatives of surface Gaussian widths (3, nsph, ngrid)
      real(wp), allocatable :: xi1_rA(:, :, :)
      !> Nuclear derivatives of Gaussian switching factors (3, nsph, ngrid)
      real(wp), allocatable :: f1_rA(:, :, :)
      !> Nuclear derivatives of surface positions (3, 3, nsph, ngrid)
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
      !> Nuclear derivatives of surface point volumes (3, nsph, ngrid)
      !>
      !> - contracting over the grid gives the total-volume nuclear gradient
      real(wp), allocatable :: v1_rA(:, :, :)
   contains
      procedure(update_cavity), deferred :: update
      procedure(get_cavity_gradient), deferred :: get_gradient
      !> Declare the readable results this cavity currently holds
      procedure :: list_fields => list_cavity_fields_base
      !> Map accumulated surface-observable adjoints to host response channels
      procedure :: get_surface_response => get_cavity_surface_response_default
      !> Contract accumulated surface-observable adjoints into the nuclear gradient
      procedure :: get_surface_gradient => get_cavity_surface_gradient_default
      !> Declare the cavity's own host requests (none by default)
      procedure :: declare_coupling => declare_cavity_coupling_default
      !> Whether the surface geometry responds to the host field (DROP), so that
      !> host surface weights are consumed in the response phase
      procedure :: has_field_dependent_geometry => cavity_field_dependent_default
      !> Write grid to XYZ file for visualization
      procedure :: write_xyz_debug => write_cavity_xyz_debug
      !> Write grid to CSV file for visualization
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> Write grid to PQR file for visualization
      procedure :: write_pqr_debug => write_cavity_pqr_debug
      !> Print basic cavity information
      procedure :: print => print_cavity_info
   end type cavity_type

   ! Abstract interfaces for deferred procedures
   abstract interface

      subroutine update_cavity(self, mol, error)
         import :: cavity_type, structure_type, wp, error_type
         implicit none(type, external)
         class(cavity_type), intent(inout) :: self
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_cavity

      subroutine get_cavity_gradient(self, error)
         import :: cavity_type, error_type
         implicit none(type, external)
         class(cavity_type), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_cavity_gradient

   end interface

contains

   !> Output unit for a cavity: the borrowed run context's unit when one is
   !> attached, otherwise the standard output unit
   !>
   !> - the base procedures below are reachable on a cavity that never went
   !>   through a constructor, so the association is checked, not assumed
   pure function cavity_unit(self) result(iunit)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Unit to write to
      integer :: iunit

      iunit = output_unit
      if (associated(self%ctx)) iunit = self%ctx%unit

   end function cavity_unit

   !* ================================================================================= *!
   !*                                 Readable results                                *!
   !* ================================================================================= *!

   !> Declare the results every discretized cavity carries
   !>
   !> @param[in]    self   Cavity instance
   !> @param[inout] query  Walker collecting or fetching the declarations
   subroutine list_cavity_fields_base(self, query)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(cavity_field_query_type), intent(inout) :: query

      call query%add_int_value("ngrid", "Number of surface grid points", self%ngrid)
      call query%add_int_value("nsph", "Number of atomic spheres", self%nsph)
      call query%add_real_scalar("area", "Total surface area, bohr**2", self%total_area)
      call query%add_real_scalar("volume", "Total enclosed volume, bohr**3", self%total_volume)

      call query%add_real2("xyz", "Grid point coordinates, bohr (3, ngrid)", self%xyz)
      call query%add_real("a", "Grid point area, bohr**2 (ngrid)", self%a)
      call query%add_real("f", "Gaussian switching factor per grid point (ngrid)", self%f)
      call query%add_real("xi0", "Gaussian width per grid point (ngrid)", self%xi0)
      call query%add_real2("normal0", "Outward unit normal per grid point (3, ngrid)", self%normal0)
      call query%add_real("v", "Grid point volume element, bohr**3 (ngrid)", self%v)
      call query%add_int("owner", "Sphere each grid point belongs to, 0-based (ngrid)", &
         & self%owner, zero_based=.true.)

      call query%add_real("radii", "Sphere radii, bohr (nsph)", self%radii)
      call query%add_real("asph", "Surface area per sphere, bohr**2 (nsph)", self%asph)
      call query%add_real2("sphxyz", "Sphere centre coordinates, bohr (3, nsph)", self%sphxyz)

   end subroutine list_cavity_fields_base

   !* ================================================================================= *!
   !*                              Surface adjoint hooks                              *!
   !* ================================================================================= *!

   !> Default surface-response hook for cavities without field-dependent geometry
   !>
   !> @param[inout] self      Cavity instance, unchanged
   !> @param[in]    acc       Surface-observable adjoints, unused
   !> @param[inout] response  Response accumulator, unchanged
   !> @param[out]   error     Error handling
   subroutine get_cavity_surface_response_default(self, acc, response, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_cavity_surface_response_default

   !> Default reverse-mode nuclear-gradient hook
   !>
   !> Cavities that do not implement the surface-adjoint contraction must be
   !> reached through the forward path instead; returning silently here would
   !> hand back a zero gradient, so it errors
   !>
   !> @param[in]    self     Cavity instance
   !> @param[in]    acc      Surface-observable adjoints, unused
   !> @param[inout] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error handling
   subroutine get_cavity_surface_gradient_default(self, acc, gradient, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "This cavity does not provide a reverse-mode surface gradient")

   end subroutine get_cavity_surface_gradient_default

   !> Default geometry query: the surface does not respond to the host field
   !>
   !> A cavity that overrides `get_surface_response` (DROP) overrides this to
   !> `.true.`; it decides whether host surface weights are declared for the
   !> response phase or only for the gradient
   !>
   !> @param[in] self Cavity instance
   function cavity_field_dependent_default(self) result(field_dependent)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Whether the surface geometry depends on the host field
      logical :: field_dependent

      field_dependent = .false.

   end function cavity_field_dependent_default

   !* ================================================================================= *!
   !*                                  Host coupling                                  *!
   !* ================================================================================= *!

   !> Default cavity-side declaration: a cavity requests nothing of its own
   !>
   !> @param[in]    self     Cavity instance
   !> @param[inout] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine declare_cavity_coupling_default(self, coupling, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Coupling being declared
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine declare_cavity_coupling_default

   !* ================================================================================= *!
   !*                                   Diagnostics                                   *!
   !* ================================================================================= *!

   !> Write grid points to an XYZ file as helium atoms (debug visualization)
   subroutine write_cavity_xyz_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_xyz_debug: cavity grid not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_xyz_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open XYZ file for writing: "//trim(filename))
         return
      end if

      write (unit, "(i0)") self%ngrid
      write (unit, "(a)") "drop cavity grid points as He (Angstrom)"
      do i = 1, self%ngrid
         write (unit, "(a2,1x,3f16.8)") "He", &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_xyz_debug

   !> Write grid points to a CSV file (debug visualization)
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: stat, i, unit

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_csv_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_csv_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_csv_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_csv_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open CSV file for writing: "//trim(filename))
         return
      end if

      write (unit, "(a)") "ngrid,x,y,z,owner,area"

      do i = 1, self%ngrid
         write (unit, "(i0,7(',',g0))") i, &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), self%a(i)
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_csv_debug

   !> Write grid points to a PQR file (debug visualization)
   !>
   !> Grid points are written as HETATM records with:
   !> - positions converted from bohr to Angstrom
   !> - charge set to 0.0
   !> - radius set to the final adapted integration weight `a(i)` (area element
   !>   with switching function applied), also converted to Angstrom
   !>
   !> @param[in]  self      Cavity instance
   !> @param[in]  filename  Output PQR file path
   subroutine write_cavity_pqr_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_pqr_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_pqr_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_pqr_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_pqr_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open PQR file for writing: "//trim(filename))
         return
      end if

      do i = 1, self%ngrid
         write (unit, "(a6,i5,1x,a4,a1,a3,1x,a1,i4,4x,3f8.3,f8.4,f7.4)") &
            "HETATM", i, "GP  ", " ", "GRD", "A", self%owner(i), &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa, &
            0.0_wp, &
            (sqrt(self%a(i)/(2.0_wp*pi))*autoaa + 0.0001_wp)
      end do
      write (unit, "(a)") "END"
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity PQR to", trim(filename)

   end subroutine write_cavity_pqr_debug

   !> Print basic cavity information (grid points, total area, total volume)
   subroutine print_cavity_info(self, unit)
      class(cavity_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer :: iunit
      type(prettyprinter) :: pp

      iunit = cavity_unit(self)
      if (present(unit)) iunit = unit

      if (.not. allocated(self%total_area) .or. .not. allocated(self%total_volume)) then
         write (iunit, "(a)") "[Warning] Cavity not fully initialized"
         return
      end if

      pp = new_prettyprinter(unit=iunit, fmt_len=20)

      call pp%blank()
      call pp%push("Results:")
      call pp%kv("Cavity points", self%ngrid)
      call pp%kv("Total area", self%total_area, "bohr^2")
      call pp%kv("Total volume", self%total_volume, "bohr^3")
      call pp%pop()
      call pp%blank()

   end subroutine print_cavity_info

end module moist_cavity_type
