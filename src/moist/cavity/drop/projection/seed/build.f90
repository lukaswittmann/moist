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
   !> - the level set clones live for the whole build ([[drop_worker_slots_type]]);
   !>   a clone may borrow from `lsf` until the build returns
   !> - on failure the seeder is left as after `init`, without a lattice
   !>
   !> @param[inout] self     Seeder (initialized)
   !> @param[in]    lsf      Level set, updated for the geometry; cloned per thread
   !> @param[in]    geometry Geometry and the search ahead
   !> @param[out]   error    Level set failure, lattice out of range, invalid input,
   !>                         axes not a rotation
   module subroutine seeder_build(self, lsf, geometry, error)
      !> Seeder (initialized)
      class(projection_seeder_type), intent(inout) :: self
      !> Level set, updated for the geometry; cloned per thread
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
      !> Geometry and the search ahead
      type(seed_build_type), intent(in) :: geometry
      !> Level set failure, lattice out of range, invalid input,
      !> axes not a rotation
      type(error_type), allocatable, intent(out) :: error

      type(build_type) :: work
      integer :: nat, nanc, islot

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
      if (.not. all(ieee_is_finite(geometry%centers))) then
         call fatal_error(error, "Lattice seeder: atom centers must be finite")
         return
      end if
      if (.not. (geometry%reach > 0.0_wp)) then
         call fatal_error(error, "Lattice seeder: reach must be positive")
         return
      end if
      if (.not. ieee_is_finite(geometry%reach)) then
         call fatal_error(error, "Lattice seeder: reach must be finite")
         return
      end if
      if (.not. is_rotation(geometry%axes)) then
         call fatal_error(error, "Lattice seeder: lattice axes must be a proper rotation")
         return
      end if
      nanc = 0
      if (allocated(geometry%anchors)) then
         nanc = size(geometry%anchors, 2)
         if (.not. all(ieee_is_finite(geometry%anchors))) then
            call fatal_error(error, "Lattice seeder: anchors must be finite")
            return
         end if
      end if
      if (allocated(geometry%radius_max)) then
         if (size(geometry%radius_max) /= nanc) then
            call fatal_error(error, "Lattice seeder: one largest search radius per anchor required")
            return
         end if
         if (.not. all(ieee_is_finite(geometry%radius_max))) then
            call fatal_error(error, "Lattice seeder: largest search radii must be finite")
            return
         end if
      end if

      if (nanc == 0) then
         self%built = .true.
         return
      end if
      ! No centroid to place the lattice on
      if (nat == 0) then
         call fatal_error(error, "Lattice seeder: anchors without an atom center")
         return
      end if

      call setup_lattice(self, lsf, geometry, work, error)
      if (.not. allocated(error)) then
         call work%slots%init(geometry%nthreads)
         call build_stages(self, lsf, geometry, work, error)
         do islot = 1, work%slots%nthreads
            call work%slots%release(islot)
         end do
      end if
      if (allocated(error)) then
         call clear_lattice(self)
         return
      end if
      self%built = .true.
   end subroutine seeder_build

   !> Stages of a build on a placed lattice: discovery, refinement, samples, bins
   !>
   !> @param[inout] self     Seeder, lattice written
   !> @param[in]    lsf      Level set, the source of the clones
   !> @param[in]    geometry Build input
   !> @param[inout] work     Build scratch, lattice placed, slots sized
   !> @param[out]   error    Level set failure or lattice out of range
   subroutine build_stages(self, lsf, geometry, work, error)
      !> Seeder, lattice written
      type(projection_seeder_type), intent(inout) :: self
      !> Level set, the source of the clones
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
      !> Build input
      type(seed_build_type), intent(in) :: geometry
      !> Build scratch, lattice placed, slots sized
      type(build_type), intent(inout) :: work
      !> Level set failure or lattice out of range
      type(error_type), allocatable, intent(out) :: error

      self%frame = work%frame
      call discover_volume(self, lsf, geometry, work, error)
      if (allocated(error)) return
      call refine_cells(self, lsf, work, error)
      if (allocated(error)) return
      call retract_samples(self, lsf, work, error)
      if (allocated(error)) return
      call build_bins(self, max(discovery_size, self%sample_size))
      call count_bytes(self)
   end subroutine build_stages

   !> Drop the built lattice
   !>
   !> @param[inout] self Seeder
   module subroutine clear_lattice(self)
      !> Seeder
      type(projection_seeder_type), intent(inout) :: self
      !> Default frame; ifx/ifort reject the empty constructor (error #6595)
      type(lattice_frame_type) :: frame

      self%built = .false.
      self%frame = frame
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
      self%n_order1 = 0
      self%n_order2 = 0
      self%n_unresolved = 0
      self%lattice_bytes = 0_int64
   end subroutine clear_lattice

   !> Lattice origin, spacings and levels
   !>
   !> - raises the coarse spacing cap up to [[discovery_size_certified]] when the
   !>   level set certifies surface-free balls large enough to clear coarse cells
   !>   (its radius at the unscaled level 10, `S = 10 scale`, or around an unreached point)
   !> - refuses a depth whose coarse step leaves no coordinate in the key range
   !>   (`depth > key_bits - 3`)
   !>
   !> @param[in]    self     Seeder
   !> @param[in]    lsf      Level set (exclusion certificate and its strength)
   !> @param[in]    geometry Build input
   !> @param[inout] work     Build scratch
   !> @param[out]   error    Lattice out of range
   subroutine setup_lattice(self, lsf, geometry, work, error)
      !> Seeder
      type(projection_seeder_type), intent(in) :: self
      !> Level set (exclusion certificate and its strength)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Build input
      type(seed_build_type), intent(in) :: geometry
      !> Build scratch
      type(build_type), intent(inout) :: work
      !> Lattice out of range
      type(error_type), allocatable, intent(out) :: error

      integer :: nat
      real(wp) :: h_cap, r_cert

      nat = size(geometry%centers, 2)
      call work%frame%init(sum(geometry%centers, dim=2)/real(nat, wp), geometry%axes)
      h_cap = discovery_size
      ! Zero for a level set without a surface-free radius: the cap stays
      r_cert = max(lsf%far_exclusion_radius(), lsf%exclusion_radius(10.0_wp*lsf%scale))
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
      ! Coarse step `2^depth` within the half key range of the doubled edge keys
      if (work%depth > key_bits - 3) then
         call range_error(error)
         return
      end if
      work%h0 = self%sample_size*2.0_wp**work%m
      work%unit = self%sample_size/2.0_wp**self%refine_levels
   end subroutine setup_lattice

   !> Bytes held by the built lattice
   !>
   !> @param[inout] self Seeder
   subroutine count_bytes(self)
      !> Seeder
      type(projection_seeder_type), intent(inout) :: self

      integer(int64) :: b

      b = int(self%nsample, int64)*int(14*storage_size(1.0_wp)/8 + storage_size(1)/8, int64)
      if (allocated(self%bins%key)) b = b + int(size(self%bins%key), int64)*12_int64
      if (allocated(self%bin_first)) b = b + int(size(self%bin_first), int64)*4_int64
      self%lattice_bytes = b
   end subroutine count_bytes

end submodule moist_cavity_drop_projection_seed_build
