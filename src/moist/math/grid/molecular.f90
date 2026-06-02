!> Atom-centered molecular integration grid.
!>
!> Combines Chebyshev-2 radial quadrature (`moist_math_grid_radial`),
!> Lebedev angular grids (`moist_math_grid_lebedev`), and Becke fuzzy-cell
!> partitioning (`moist_math_grid_becke`) into a single product grid that
!> approximates molecular volume integrals:
!>
!>    integral f(r) dV  approx  sum_i weights(i) * f(xyz(:,i))
!>
!> Two public constructors are provided:
!>
!>   `new_molecular_grid(self, mol, error)`
!>        Per-element default sizes (literature Pople-style defaults).
!>
!>   `new_molecular_grid_uniform(self, mol, nrad, nang, error, rmin, rmax)`
!>        Uniform `(nrad, nang)` per atom; `nang` is the RAW Lebedev point
!>        count (6, 14, 26, ..., 5810). Optional `rmin`/`rmax` clamp the
!>        radial shells in bohr.
!>
!> Units: all lengths are in bohr. Covalent radii come from
!> `moist_data_atomicrad` which stores them already in bohr.
module moist_math_grid_molecular
   use iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use moist_data_atomicrad, only: covalent_rad
   use moist_math_grid_lebedev, only: get_angular_grid, lebedev_order_from_num
   use moist_math_grid_radial, only: chebyshev2_radii
   use moist_math_grid_becke, only: becke_weights
   implicit none
   private

   public :: molecular_grid_type
   public :: new_molecular_grid
   public :: new_molecular_grid_uniform
   public :: default_grid_sizes
   public :: integrand_3d

   !> Default weight-pruning threshold (bohr^3); points with |w| below this
   !> are dropped from the final grid.
   real(wp), parameter :: default_wthr = 1.0e-14_wp

   !> Abstract interface for scalar-valued 3D integrands used by
   !> `molecular_grid_type%integrate`.
   abstract interface
      pure function integrand_3d(r) result(val)
         import :: wp
         !> Point in bohr
         real(wp), intent(in) :: r(3)
         !> Function value at r
         real(wp) :: val
      end function integrand_3d
   end interface

   !> Atom-centered molecular integration grid
   !> (Chebyshev-2 radial x Lebedev angular x Becke fuzzy-cell partitioning).
   type :: molecular_grid_type
      !> Total number of retained grid points (after Becke pruning)
      integer :: npts = 0
      !> Grid point coordinates in bohr, shape (3, npts)
      real(wp), allocatable :: xyz(:, :)
      !> Integration weights, shape (npts); include r^2 dr Jacobian,
      !> 4*pi solid-angle factor, and Becke atomic weight w_A.
      real(wp), allocatable :: weights(:)
      !> CSR-style atom ownership: points owned by atom i are stored
      !> contiguously at indices [atom_offset(i), atom_offset(i+1)-1];
      !> shape (nat+1). atom_offset(1) = 1, atom_offset(nat+1) = npts+1.
      integer, allocatable :: atom_offset(:)
      !> Per-atom radial size actually used, shape (nat)
      integer, allocatable :: nrad_per_atom(:)
      !> Per-atom Lebedev point count actually used, shape (nat)
      integer, allocatable :: nang_per_atom(:)
   contains
      !> Free all array components (idempotent)
      procedure :: destroy => molecular_grid_destroy
      !> Print a short summary to the given unit (default output_unit)
      procedure :: info => molecular_grid_info
      !> Evaluate sum_i w_i * f(xyz(:,i))
      procedure :: integrate => molecular_grid_integrate
   end type molecular_grid_type

contains

   !> Literature-standard per-element default grid sizes (Pople "fine" /
   !> "ultrafine" style):
   !>   H, He            -> (50, 302)
   !>   Li..Ne (row 2)   -> (75, 302)
   !>   Na..Ar (row 3)   -> (75, 434)
   !>   K  and beyond    -> (99, 590)
   !>
   !> `nang` is the raw Lebedev point count; it is always one of the
   !> supported sizes.
   pure subroutine default_grid_sizes(iz, nrad, nang)
      !> Atomic number
      integer, intent(in)  :: iz
      !> Number of radial points
      integer, intent(out) :: nrad
      !> Raw Lebedev point count
      integer, intent(out) :: nang

      if (iz <= 2) then
         nrad = 50; nang = 302
      else if (iz <= 10) then
         nrad = 75; nang = 302
      else if (iz <= 18) then
         nrad = 75; nang = 434
      else
         nrad = 99; nang = 590
      end if
   end subroutine default_grid_sizes

   !> Build a molecular grid using element-dependent default sizes.
   !>
   !> @param[out] self   Initialised molecular grid.
   !> @param[in]  mol    Molecular structure (mctc_io `structure_type`).
   !> @param[out] error  Propagated error (invalid grid sizes, alloc, ...).
   subroutine new_molecular_grid(self, mol, error)
      !> Grid to initialise
      type(molecular_grid_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: nrad_atom(:), nang_atom(:)
      integer :: iat, iz

      allocate (nrad_atom(mol%nat), nang_atom(mol%nat))
      do iat = 1, mol%nat
         iz = mol%num(mol%id(iat))
         call default_grid_sizes(iz, nrad_atom(iat), nang_atom(iat))
      end do

      call build_molecular_grid(self, mol, nrad_atom, nang_atom, &
         & default_wthr, error)
   end subroutine new_molecular_grid

   !> Build a molecular grid with uniform `(nrad, nang)` per atom.
   !>
   !> @param[out] self   Initialised molecular grid.
   !> @param[in]  mol    Molecular structure.
   !> @param[in]  nrad   Number of radial points per atom (>= 1).
   !> @param[in]  nang   Raw Lebedev point count per atom (must be one
   !>                    of the supported sizes: 6, 14, 26, ..., 5810).
   !> @param[out] error  Propagated error.
   !> @param[in]  rmin   Optional minimum radial shell radius (bohr).
   !> @param[in]  rmax   Optional maximum radial shell radius (bohr).
   subroutine new_molecular_grid_uniform(self, mol, nrad, nang, error, rmin, rmax)
      !> Grid to initialise
      type(molecular_grid_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Number of radial points per atom
      integer, intent(in)  :: nrad
      !> Raw Lebedev point count per atom
      integer, intent(in)  :: nang
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Minimum radius (bohr)
      real(wp), optional, intent(in)  :: rmin
      !> Maximum radius (bohr)
      real(wp), optional, intent(in)  :: rmax

      integer, allocatable :: nrad_atom(:), nang_atom(:)

      if (nrad < 1) then
         call fatal_error(error, "molecular grid: nrad must be >= 1")
         return
      end if

      allocate (nrad_atom(mol%nat), nang_atom(mol%nat))
      nrad_atom = nrad
      nang_atom = nang

      call build_molecular_grid(self, mol, nrad_atom, nang_atom, &
         & default_wthr, error, rmin=rmin, rmax=rmax)
   end subroutine new_molecular_grid_uniform

   !> Internal worker: build atom-centered molecular grid with per-atom
   !> (nrad, nang) arrays and optional radial clamp.
   subroutine build_molecular_grid(self, mol, nrad_atom, nang_atom, wthr, &
         & error, rmin, rmax)
      !> Grid to initialise
      type(molecular_grid_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Per-atom radial sizes, shape (nat)
      integer, intent(in)  :: nrad_atom(:)
      !> Per-atom raw Lebedev point counts, shape (nat)
      integer, intent(in)  :: nang_atom(:)
      !> Pruning threshold (bohr^3)
      real(wp), intent(in)  :: wthr
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional radial clamps (bohr)
      real(wp), optional, intent(in)  :: rmin, rmax

      integer  :: nat, iat, iz, ir, il, ig, nr, nl, max_pts, order
      integer, allocatable :: numbers(:)
      real(wp) :: p, r, full_weight
      real(wp), allocatable :: radii(:), rad_w(:)
      real(wp), allocatable :: ang_xyz(:, :), ang_w(:)
      real(wp), allocatable :: tmp_xyz(:, :), tmp_w(:)
      integer, allocatable :: tmp_atom(:)
      real(wp) :: bweights_buf(mol%nat)
      real(wp), parameter :: four_pi = 4.0_wp*pi

      nat = mol%nat

      ! Resolve atomic numbers once (mctc_io indexes through mol%id).
      allocate (numbers(nat))
      do iat = 1, nat
         numbers(iat) = mol%num(mol%id(iat))
      end do

      ! Upper bound on grid size before pruning.
      max_pts = 0
      do iat = 1, nat
         max_pts = max_pts + nrad_atom(iat)*nang_atom(iat)
      end do

      allocate (tmp_xyz(3, max_pts), tmp_w(max_pts), tmp_atom(max_pts))

      ig = 0
      do iat = 1, nat
         iz = numbers(iat)
         nr = nrad_atom(iat)
         nl = nang_atom(iat)

         ! Chebyshev-2 radial scale p (Becke 1988 convention).
         if (iz == 1) then
            p = covalent_rad(1)
         else
            p = 0.5_wp*covalent_rad(iz)
         end if

         allocate (radii(nr), rad_w(nr))
         call chebyshev2_radii(nr, p, radii, rad_w)

         ! Lebedev angular grid (use existing umbrella getter).
         call lebedev_order_from_num(nl, order, error)
         if (allocated(error)) then
            deallocate (radii, rad_w)
            return
         end if
         allocate (ang_xyz(3, nl), ang_w(nl))
         call get_angular_grid(order, ang_xyz, ang_w, error)
         if (allocated(error)) then
            deallocate (radii, rad_w, ang_xyz, ang_w)
            return
         end if

         ! Combine radial x angular, shift to atom center, apply Becke weight.
         do ir = 1, nr
            r = radii(ir)
            if (present(rmin)) then
               if (r < rmin) cycle
            end if
            if (present(rmax)) then
               if (r > rmax) cycle
            end if
            do il = 1, nl
               ig = ig + 1
               tmp_xyz(1, ig) = r*ang_xyz(1, il) + mol%xyz(1, iat)
               tmp_xyz(2, ig) = r*ang_xyz(2, il) + mol%xyz(2, iat)
               tmp_xyz(3, ig) = r*ang_xyz(3, il) + mol%xyz(3, iat)
               full_weight = rad_w(ir)*ang_w(il)*four_pi
               call becke_weights(tmp_xyz(:, ig), nat, mol%xyz, numbers, &
                  & bweights_buf)
               tmp_w(ig) = full_weight*bweights_buf(iat)
               tmp_atom(ig) = iat
            end do
         end do

         deallocate (radii, rad_w, ang_xyz, ang_w)
      end do

      call finalise_grid(self, nat, nrad_atom, nang_atom, &
         & ig, tmp_xyz, tmp_w, tmp_atom, wthr)

      deallocate (tmp_xyz, tmp_w, tmp_atom, numbers)
   end subroutine build_molecular_grid

   !> Compact the raw grid buffer into the final `molecular_grid_type`,
   !> dropping points with |w| < wthr and rebuilding `atom_offset`.
   pure subroutine finalise_grid(self, nat, nrad_atom, nang_atom, &
         & nraw, raw_xyz, raw_w, raw_atom, wthr)
      type(molecular_grid_type), intent(out) :: self
      integer, intent(in) :: nat
      integer, intent(in) :: nrad_atom(:)
      integer, intent(in) :: nang_atom(:)
      integer, intent(in) :: nraw
      real(wp), intent(in) :: raw_xyz(:, :)
      real(wp), intent(in) :: raw_w(:)
      integer, intent(in) :: raw_atom(:)
      real(wp), intent(in) :: wthr

      integer :: i, iat, npts

      ! First pass: count retained points.
      npts = 0
      do i = 1, nraw
         if (abs(raw_w(i)) >= wthr) npts = npts + 1
      end do

      self%npts = npts
      allocate (self%xyz(3, npts), self%weights(npts))
      allocate (self%atom_offset(nat + 1))
      allocate (self%nrad_per_atom(nat), self%nang_per_atom(nat))

      self%nrad_per_atom = nrad_atom
      self%nang_per_atom = nang_atom

      ! Second pass: copy retained points and build atom_offset.
      ! raw_atom is monotonically non-decreasing because we emit atom-by-atom.
      ! Default every entry to the end sentinel so atoms with zero points
      ! get a well-defined (empty) range. Overwritten as points arrive.
      self%atom_offset = npts + 1
      self%atom_offset(1) = 1
      npts = 0
      iat = 1
      do i = 1, nraw
         ! Advance atom_offset for any atoms that start at this position.
         do while (iat < raw_atom(i))
            iat = iat + 1
            self%atom_offset(iat) = npts + 1
         end do
         if (abs(raw_w(i)) >= wthr) then
            npts = npts + 1
            self%xyz(:, npts) = raw_xyz(:, i)
            self%weights(npts) = raw_w(i)
         end if
      end do
      ! Fill remaining (possibly empty) tail atoms.
      do while (iat < nat)
         iat = iat + 1
         self%atom_offset(iat) = npts + 1
      end do
      self%atom_offset(nat + 1) = npts + 1
   end subroutine finalise_grid

   !> Free all component arrays. Idempotent; safe to call repeatedly.
   pure subroutine molecular_grid_destroy(self)
      !> Grid instance
      class(molecular_grid_type), intent(inout) :: self

      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%weights)) deallocate (self%weights)
      if (allocated(self%atom_offset)) deallocate (self%atom_offset)
      if (allocated(self%nrad_per_atom)) deallocate (self%nrad_per_atom)
      if (allocated(self%nang_per_atom)) deallocate (self%nang_per_atom)
      self%npts = 0
   end subroutine molecular_grid_destroy

   !> Write a short grid summary to the given unit.
   subroutine molecular_grid_info(self, unit)
      !> Grid instance
      class(molecular_grid_type), intent(in) :: self
      !> Output unit (default `output_unit`)
      integer, optional, intent(in) :: unit

      integer :: iunit, iat, nat

      iunit = output_unit
      if (present(unit)) iunit = unit

      write (iunit, '(a)') "moist molecular_grid_type"
      write (iunit, '(a,i0)') "  total points   : ", self%npts
      if (.not. allocated(self%xyz)) then
         write (iunit, '(a)') "  (grid is uninitialised)"
         return
      end if
      write (iunit, '(a,es12.4,a,es12.4)') &
         & "  weight range   : min = ", minval(self%weights), &
         & ", max = ", maxval(self%weights)
      write (iunit, '(a,es14.6)') "  sum(weights)   : ", sum(self%weights)
      nat = size(self%nrad_per_atom)
      write (iunit, '(a,i0)') "  atoms          : ", nat
      write (iunit, '(a)') "  per-atom counts (atom, npts, nrad, nang):"
      do iat = 1, nat
         write (iunit, '(4x,i6,3x,i8,3x,i5,3x,i5)') iat, &
            & self%atom_offset(iat + 1) - self%atom_offset(iat), &
            & self%nrad_per_atom(iat), self%nang_per_atom(iat)
      end do
   end subroutine molecular_grid_info

   !> Evaluate sum_i weights(i) * f(xyz(:,i)).
   !>
   !> @param[in]  f       Pure function with signature `f(r) -> real(wp)`.
   !> @param[out] result  Quadrature result.
   subroutine molecular_grid_integrate(self, f, result)
      !> Grid instance
      class(molecular_grid_type), intent(in)  :: self
      !> Integrand
      procedure(integrand_3d)                 :: f
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      ! Reduce with the sum() intrinsic over the per-point products so that for
      ! f == 1 (where w_i*1 == w_i exactly in IEEE) the result is bit-identical
      ! to sum(weights), regardless of how the compiler orders the reduction. A
      ! scalar accumulation loop is reassociated differently from sum() and
      ! breaks that identity.
      result = sum([(self%weights(i)*f(self%xyz(:, i)), i=1, self%npts)])
   end subroutine molecular_grid_integrate

end module moist_math_grid_molecular
