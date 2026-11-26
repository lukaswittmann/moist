!> Lattice seeder build: lattice node positions, values, certificates
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_nodes
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                   Lattice nodes                                   *!
   !* ================================================================================= *!

   !> Lattice step of level `s` in finest units
   !>
   !> @param[in] work Build scratch
   !> @param[in] s    Level
   pure module function level_step(work, s) result(step)
      type(build_type), intent(in) :: work
      integer, intent(in) :: s
      integer :: step

      step = 2**(work%depth - s)
   end function level_step

   !> Lattice-frame position of a lattice node
   !>
   !> - the lab position when the lattice follows the lab axes
   !>
   !> @param[in] work Build scratch
   !> @param[in] ijk  Finest-unit coordinates
   pure module function node_position(work, ijk) result(x)
      type(build_type), intent(in) :: work
      integer, intent(in) :: ijk(3)
      real(wp) :: x(3)

      x = work%frame%origin + work%unit*real(ijk, wp)
   end function node_position

   !> Values and certificates of the nodes waiting for them
   !>
   !> - in slices of [[build_chunk]] points
   !>
   !> @param[inout] self  Seeder (tallies)
   !> @param[in]    lsf   Level set (certificates)
   !> @param[inout] work  Build scratch
   !> @param[out]   error Level set failure
   module subroutine evaluate_nodes(self, lsf, work, error)
      type(projection_seeder_type), intent(inout) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(build_type), intent(inout) :: work
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: x(:, :)
      integer :: k, n0, n1

      n0 = work%nodes%nval + 1
      n1 = work%nodes%n
      if (n1 < n0) return
      allocate (x(3, n1 - n0 + 1))
      do k = n0, n1
         x(:, k - n0 + 1) = work%frame%to_lab(node_position(work, work%nodes%ijk(:, k)))
      end do
      call evaluate_slices(work%probe, lsf, x, 0, build_chunk, work%nodes%val(n0:n1), work%nodes%reached(n0:n1), &
                           error)
      if (allocated(error)) return
      do k = n0, n1
         if (work%nodes%reached(k)) then
            work%nodes%free(k) = lsf%exclusion_radius(work%nodes%val(k))
         else
            work%nodes%free(k) = lsf%far_exclusion_radius()
         end if
      end do
      self%n_probes = self%n_probes + (n1 - n0 + 1)
      work%nodes%nval = n1
   end subroutine evaluate_nodes

   !> Add a lattice node unless known; its value waits for [[evaluate_nodes]]
   !>
   !> - coordinates within half the key range, so that the doubled ones of the
   !>   edges fit too: the coarse nodes are checked by [[discover_volume]], every
   !>   later node lies between them
   !>
   !> @param[inout] nodes Node set
   !> @param[in]    ijk   Finest-unit coordinates
   module subroutine add_node(nodes, ijk)
      type(node_set_type), intent(inout) :: nodes
      integer, intent(in) :: ijk(3)

      integer(int64) :: key
      integer :: cap

      key = lattice_key(ijk, key_bias)
      if (table_find(nodes%table, key) > 0) return
      if (.not. allocated(nodes%ijk)) then
         allocate (nodes%ijk(3, 1024), nodes%val(1024), nodes%free(1024), nodes%reached(1024))
      end if
      if (nodes%n + 1 > size(nodes%val)) then
         cap = 2*size(nodes%val)
         call grow_array(nodes%ijk, 3, cap)
         call grow_array(nodes%val, cap)
         call grow_array(nodes%free, cap)
         call grow_array(nodes%reached, cap)
      end if
      nodes%n = nodes%n + 1
      nodes%ijk(:, nodes%n) = ijk
      call table_insert(nodes%table, key, nodes%n)
   end subroutine add_node

   !> Lattice coordinates out of the key range
   !>
   !> @param[out] error Refusal
   module subroutine range_error(error)
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Lattice seeder: the region does not fit the lattice keys; "// &
                       "raise the sample size or split the geometry")
   end subroutine range_error

end submodule moist_cavity_drop_projection_seed_nodes
