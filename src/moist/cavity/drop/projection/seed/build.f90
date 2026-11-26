!> Lattice seeder build: lattice placement, stages in order
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_build
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                       Build                                       *!
   !* ================================================================================= *!

   !> Build the seed lattice for one geometry
   !>
   !> - the anchors' largest search balls bound the region, the atom centroid
   !>   places the lattice, `geometry%axes` turns it; without an anchor the
   !>   lattice is empty
   !>
   !> @param[inout] self     Seeder (initialized)
   !> @param[in]    lsf      Level set, updated for the geometry; cloned per thread
   !> @param[in]    geometry Geometry and the search ahead
   !> @param[out]   error    Level set failure, lattice out of range, invalid input,
   !>                         axes not a rotation
   module subroutine seeder_build(self, lsf, geometry, error)
      class(projection_seeder_type), intent(inout) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(seed_build_type), intent(in) :: geometry
      type(error_type), allocatable, intent(out) :: error

      type(build_type) :: work
      integer :: nat

      call clear_lattice(self)
      if (.not. self%initialized) then
         call fatal_error(error, "Lattice seeder: build before init")
         return
      end if
      nat = -1
      if (allocated(geometry%centers) .and. allocated(geometry%radii)) nat = size(geometry%centers, 2)
      if (nat >= 0) then
         if (size(geometry%radii) /= nat) nat = -1
      end if
      if (nat < 0) then
         call fatal_error(error, "Lattice seeder: one radius per atom center required")
         return
      end if
      if (.not. (geometry%reach > 0.0_wp)) then
         call fatal_error(error, "Lattice seeder: reach must be positive")
         return
      end if
      if (.not. is_rotation(geometry%axes)) then
         call fatal_error(error, "Lattice seeder: lattice axes must be a proper rotation")
         return
      end if

      if (.not. allocated(geometry%anchors)) then
         self%built = .true.
         return
      end if
      if (size(geometry%anchors, 2) == 0) then
         self%built = .true.
         return
      end if
      ! No centroid to place the lattice on
      if (nat == 0) then
         call fatal_error(error, "Lattice seeder: anchors without an atom center")
         return
      end if

      work%probe%nthreads = max(1, geometry%nthreads)
      allocate (work%probe%slot(work%probe%nthreads))

      call setup_lattice(self, lsf, geometry, work)
      self%frame = work%frame

      call discover_volume(self, lsf, geometry, work, error)
      if (allocated(error)) return

      call refine_cells(self, lsf, work, error)
      if (allocated(error)) return
      call retract_samples(self, lsf, work, error)
      if (allocated(error)) return
      call build_bins(self, max(discovery_size, self%sample_size))
      call count_bytes(self)
      self%built = .true.
   end subroutine seeder_build

   !> Drop the built lattice
   !>
   !> @param[inout] self Seeder
   module subroutine clear_lattice(self)
      type(projection_seeder_type), intent(inout) :: self

      self%built = .false.
      self%frame = lattice_frame_type()
      self%nsample = 0
      self%cover = 0.0_wp
      self%spacing_max = 0.0_wp
      if (allocated(self%y)) deallocate (self%y)
      if (allocated(self%normal)) deallocate (self%normal)
      if (allocated(self%form)) deallocate (self%form)
      if (allocated(self%kmax)) deallocate (self%kmax)
      if (allocated(self%spacing)) deallocate (self%spacing)
      if (allocated(self%edge_mid)) deallocate (self%edge_mid)
      if (allocated(self%edge_axis)) deallocate (self%edge_axis)
      if (allocated(self%bin_first)) deallocate (self%bin_first)
      call table_init(self%bins, 0)
      self%n_probes = 0
      self%n_cross = 0
      self%n_order1 = 0
      self%n_order2 = 0
      self%n_unresolved = 0
      self%lattice_bytes = 0_int64
   end subroutine clear_lattice

   !> Lattice origin, spacings and levels
   !>
   !> - raises the coarse spacing cap up to [[discovery_size_certified]] when the
   !>   level set certifies surface-free balls large enough to clear coarse cells
   !>   (its radius at `S = 10`, or around an unreached point)
   !>
   !> @param[in]    self     Seeder
   !> @param[in]    lsf      Level set (exclusion certificate and its strength)
   !> @param[in]    geometry Build input
   !> @param[inout] work     Build scratch
   subroutine setup_lattice(self, lsf, geometry, work)
      type(projection_seeder_type), intent(in) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(seed_build_type), intent(in) :: geometry
      type(build_type), intent(inout) :: work

      integer :: nat
      real(wp) :: h_cap, r_cert

      nat = size(geometry%centers, 2)
      call work%frame%init(sum(geometry%centers, dim=2)/real(nat, wp), geometry%axes)
      h_cap = discovery_size
      ! Zero for a level set without a surface-free radius: the cap stays
      r_cert = max(lsf%far_exclusion_radius(), lsf%exclusion_radius(10.0_wp))
      if (r_cert >= sqrt(3.0_wp)*discovery_size_certified) then
         h_cap = discovery_size_certified
      else if (r_cert >= sqrt(3.0_wp)*2.0_wp*discovery_size) then
         h_cap = 2.0_wp*discovery_size
      end if
      work%m = 0
      do while (self%sample_size*2.0_wp**(work%m + 1) <= h_cap*(1.0_wp + epsilon(1.0_wp)))
         work%m = work%m + 1
      end do
      work%depth = work%m + self%refine_levels
      work%h0 = self%sample_size*2.0_wp**work%m
      work%unit = self%sample_size/2.0_wp**self%refine_levels
   end subroutine setup_lattice

   !> Bytes held by the built lattice
   !>
   !> @param[inout] self Seeder
   subroutine count_bytes(self)
      type(projection_seeder_type), intent(inout) :: self

      integer(int64) :: b

      b = int(self%nsample, int64)*int(14*storage_size(1.0_wp)/8 + storage_size(1)/8, int64)
      if (allocated(self%bins%key)) b = b + int(size(self%bins%key), int64)*12_int64
      if (allocated(self%bin_first)) b = b + int(size(self%bin_first), int64)*4_int64
      self%lattice_bytes = b
   end subroutine count_bytes

end submodule moist_cavity_drop_projection_seed_build
