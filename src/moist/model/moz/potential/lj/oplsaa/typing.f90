!> Native OPLS-AA LJ typing from an explicit molecular graph
module moist_model_moz_potential_lj_oplsaa_typing
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_symbols, only: to_number
   use moist_model_moz_potential_lj_oplsaa_parameters, only: oplsaa_key_len
   use moist_model_moz_potential_lj_oplsaa_rules, only: oplsaa_rules, oplsaa_type_len
   implicit none(type, external)
   private

   public :: generate_oplsaa_atomtypes

   !> Maximum atoms in the bundled environment patterns
   integer, parameter :: max_pattern_atoms = 32

   !> Compiled connected atom environment
   type :: pattern_type
      !> Number of atoms
      integer :: nat = 0
      !> Atom predicates in traversal order
      character(len=96) :: expression(max_pattern_atoms) = ""
      !> Required bonds between pattern atoms
      logical :: edge(max_pattern_atoms, max_pattern_atoms) = .false.
   end type pattern_type

   !> Molecular graph in structure order
   type :: graph_type
      !> Atomic number per atom
      integer, allocatable :: num(:)
      !> Formal charge per atom
      integer, allocatable :: charge(:)
      !> Explicit neighbour count per atom
      integer, allocatable :: degree(:)
      !> Neighbour indices, (maximum degree, nat)
      integer, allocatable :: neighbor(:, :)
      !> Bond orders aligned with the neighbours
      integer, allocatable :: order(:, :)
      !> Simple cycles of length 3-6 containing each atom
      integer, allocatable :: rings(:)
      !> Ring-size membership, (nat, 6)
      logical, allocatable :: cyclic(:, :)
   end type graph_type

contains

   !> Assign OPLS-AA rule labels and existing LJ keys
   !>
   !> @param[in] mol Explicit hydrogens and bonds (3, nbd), orders 1-4
   !> @param[out] keys LJ keys in structure order, allocated only on success
   !> @param[out] error Invalid graph, unresolved environment or missing LJ
   !> @param[out] atomtypes Optional original rule labels in structure order
   subroutine generate_oplsaa_atomtypes(mol, keys, error, atomtypes)
      !> Molecular structure
      class(structure_type), intent(in) :: mol
      !> Existing LJ lookup keys
      character(len=oplsaa_key_len), allocatable, intent(out) :: keys(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Original environment labels
      character(len=oplsaa_type_len), allocatable, intent(out), optional :: atomtypes(:)
      type(graph_type) :: graph
      type(pattern_type), allocatable :: patterns(:)
      logical, allocatable :: matched(:, :), updated(:, :), excluded(:)
      character(len=oplsaa_key_len), allocatable :: typed_keys(:)
      character(len=oplsaa_type_len), allocatable :: typed_names(:)
      character(len=:), allocatable :: alternatives
      character(len=32) :: label
      integer :: irule, iat, iteration, chosen, first, last, override
      logical :: valid

      call build_graph(mol, graph, error)
      if (allocated(error)) return
      allocate (patterns(size(oplsaa_rules)))
      do irule = 1, size(patterns)
         call compile_pattern(trim(oplsaa_rules(irule)%smarts), patterns(irule), valid)
         if (.not. valid) then
            call fatal_error(error, "Invalid bundled OPLS-AA rule '"//trim(oplsaa_rules(irule)%name)//"'")
            return
         end if
      end do

      allocate (matched(size(patterns), mol%nat), updated(size(patterns), mol%nat))
      matched = .false.
      do iteration = 1, size(patterns) + 2
         do iat = 1, mol%nat
            do irule = 1, size(patterns)
               updated(irule, iat) = matches_pattern(patterns(irule), graph, iat, matched)
            end do
         end do
         if (all(updated .eqv. matched)) exit
         matched = updated
      end do
      if (iteration > size(patterns) + 2) then
         call fatal_error(error, "OPLS-AA type-reference rules did not converge")
         return
      end if

      allocate (excluded(size(patterns)), typed_keys(mol%nat), typed_names(mol%nat))
      do iat = 1, mol%nat
         excluded = .false.
         do irule = 1, size(patterns)
            if (.not. matched(irule, iat)) cycle
            first = 1
            do while (first <= len_trim(oplsaa_rules(irule)%overrides))
               last = index(oplsaa_rules(irule)%overrides(first:), ",")
               if (last == 0) then
                  last = len_trim(oplsaa_rules(irule)%overrides)
               else
                  last = first + last - 2
               end if
               override = rule_index(oplsaa_rules(irule)%overrides(first:last))
               if (override > 0) excluded(override) = .true.
               first = last + 2
            end do
         end do
         write (label, "(i0)") iat
         chosen = 0
         alternatives = ""
         do irule = 1, size(patterns)
            if (.not. matched(irule, iat) .or. excluded(irule)) cycle
            if (chosen > 0) alternatives = alternatives//", "
            alternatives = alternatives//trim(oplsaa_rules(irule)%name)
            chosen = irule
         end do
         typed_keys(iat) = ""
         typed_names(iat) = ""
         ! No rule, an unsupported charge state or a rule without imported
         ! parameters leaves the atom untyped for a later term of the potential
         if (chosen == 0) cycle
         if (count(matched(:, iat) .and. .not. excluded) /= 1) then
            call fatal_error(error, "Ambiguous OPLS-AA environment for atom "//trim(label)//": "//alternatives)
            return
         end if
         if (.not. valid_rule_charge(oplsaa_rules(chosen)%name, graph%charge(iat))) cycle
         if (len_trim(oplsaa_rules(chosen)%key) == 0) cycle
         typed_keys(iat) = oplsaa_rules(chosen)%key
         typed_names(iat) = oplsaa_rules(chosen)%name
      end do
      call move_alloc(typed_keys, keys)
      if (present(atomtypes)) call move_alloc(typed_names, atomtypes)
   end subroutine generate_oplsaa_atomtypes

   !> Validate explicit connectivity and build neighbour and ring tables
   !>
   !> @param[in] mol Explicit molecular graph
   !> @param[out] graph Validated graph
   !> @param[out] error Missing topology, inconsistent annotations or valence
   subroutine build_graph(mol, graph, error)
      !> Molecular structure
      class(structure_type), intent(in) :: mol
      !> Graph
      type(graph_type), intent(out) :: graph
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: charge(:)
      character(len=32) :: label
      integer :: ibond, iat, i, j, kind, depth, path(6)

      if (mol%nat < 1) then
         call fatal_error(error, "OPLS-AA typing needs at least one atom")
         return
      end if
      if (.not. allocated(mol%id) .or. .not. allocated(mol%num)) then
         call fatal_error(error, "OPLS-AA typing needs initialized atomic numbers")
         return
      end if
      if (size(mol%id) /= mol%nat) then
         call fatal_error(error, "OPLS-AA atom identity count does not match the structure")
         return
      end if
      if (any(mol%id < 1) .or. any(mol%id > size(mol%num))) then
         call fatal_error(error, "OPLS-AA atom identities are outside the species table")
         return
      end if
      if (.not. allocated(mol%bond)) then
         call fatal_error(error, "OPLS-AA typing needs explicit bonds and bond orders; coordinates alone are insufficient")
         return
      end if
      if (size(mol%bond, 1) /= 3 .or. size(mol%bond, 2) /= mol%nbd) then
         call fatal_error(error, "OPLS-AA bonds need shape (3, nbd): atom, atom, order")
         return
      end if
      if (mol%info%missing_hydrogen .or. mol%uhf /= 0) then
         call fatal_error(error, "OPLS-AA typing needs explicit hydrogens and a closed-shell structure")
         return
      end if
      allocate (graph%num(mol%nat), graph%degree(mol%nat), charge(mol%nat))
      graph%num = mol%num(mol%id)
      graph%degree = 0
      charge = 0
      if (allocated(mol%sdf)) then
         if (size(mol%sdf) /= mol%nat) then
            call fatal_error(error, "OPLS-AA formal charge annotations do not match the atom count")
            return
         end if
         charge = mol%sdf%charge
      end if
      if (abs(real(sum(charge), wp) - mol%charge) > 1.0e-8_wp) then
         call fatal_error(error, "OPLS-AA typing needs per-atom formal charges consistent with the total charge")
         return
      end if
      graph%charge = charge
      do ibond = 1, mol%nbd
         i = mol%bond(1, ibond)
         j = mol%bond(2, ibond)
         kind = mol%bond(3, ibond)
         if (i < 1 .or. i > mol%nat .or. j < 1 .or. j > mol%nat .or. i == j) then
            call fatal_error(error, "OPLS-AA graph contains an invalid bond endpoint")
            return
         end if
         if (kind < 1 .or. kind > 4) then
            call fatal_error(error, "OPLS-AA bond orders must be 1, 2, 3 or 4 (aromatic)")
            return
         end if
         graph%degree(i) = graph%degree(i) + 1
         graph%degree(j) = graph%degree(j) + 1
      end do
      allocate (graph%neighbor(maxval(graph%degree), mol%nat), source=0)
      allocate (graph%order(maxval(graph%degree), mol%nat), source=0)
      graph%degree = 0
      do ibond = 1, mol%nbd
         i = mol%bond(1, ibond)
         j = mol%bond(2, ibond)
         if (bonded(graph, i, j)) then
            call fatal_error(error, "OPLS-AA graph contains a duplicate bond")
            return
         end if
         graph%degree(i) = graph%degree(i) + 1
         graph%degree(j) = graph%degree(j) + 1
         graph%neighbor(graph%degree(i), i) = j
         graph%neighbor(graph%degree(j), j) = i
         graph%order(graph%degree(i), i) = mol%bond(3, ibond)
         graph%order(graph%degree(j), j) = mol%bond(3, ibond)
      end do
      do iat = 1, mol%nat
         ! Elements without rules: uncovered atoms for later terms
         if (.not. supported_element(graph%num(iat))) cycle
         if (valid_valence(graph, iat, charge(iat))) cycle
         write (label, "(i0)") iat
         call fatal_error(error, "Unsupported OPLS-AA valence at atom "//trim(label)// &
            & "; check explicit hydrogens, bond orders and formal charges")
         return
      end do
      allocate (graph%rings(mol%nat), source=0)
      allocate (graph%cyclic(mol%nat, 6), source=.false.)
      do iat = 1, mol%nat
         path = 0
         path(1) = iat
         depth = 1
         call find_cycles(graph, path, depth)
      end do
   end subroutine build_graph

   !> Basic closed-shell valence check including explicit hydrogens
   !>
   !> @param[in] graph Molecular graph
   !> @param[in] i Atom index
   !> @param[in] charge Formal charge
   pure function valid_valence(graph, i, charge) result(valid)
      !> Graph
      type(graph_type), intent(in) :: graph
      !> Atom index
      integer, intent(in) :: i
      !> Formal charge
      integer, intent(in) :: charge
      !> Whether the represented valence is supported
      logical :: valid
      integer :: valence, degree

      degree = graph%degree(i)
      valence = sum(graph%order(:degree, i))
      valid = .false.
      if (any(graph%order(:degree, i) == 4)) then
         if (count(graph%order(:degree, i) == 4) < 2) return
         if (any(graph%order(:degree, i) /= 1 .and. graph%order(:degree, i) /= 4)) return
         select case (graph%num(i))
         case (6)
            valid = degree == 3 .and. charge == 0
         case (7)
            valid = (degree == 2 .or. degree == 3) .and. charge == 0
         case (8, 16)
            valid = degree == 2 .and. charge == 0
         case default
            ! Unsupported aromatic-bond element
            valid = .false.
         end select
         return
      end if
      select case (graph%num(i))
      case (1)
         valid = degree == 1 .and. valence == 1 .and. charge == 0
      case (6)
         valid = valence == 4 .and. charge == 0
      case (7)
         valid = valence == 3 + charge .and. abs(charge) <= 1
      case (8)
         valid = valence == 2 + charge .and. abs(charge) <= 1
      case (9, 17, 35, 53)
         valid = (valence == 1 .and. charge == 0) .or. (degree == 0 .and. charge == -1)
      case (15)
         valid = (valence == 3 .or. valence == 5) .and. charge == 0
      case (16)
         valid = (valence == 2 .or. valence == 4 .or. valence == 6) .and. charge == 0
      case (3)
         valid = degree == 0 .and. charge == 1
      case default
         ! Unreachable: supported elements checked by build_graph
         valid = .false.
      end select
   end function valid_valence

   !> Whether the typing rules cover an element
   !>
   !> @param[in] num  atomic number
   pure function supported_element(num) result(supported)
      !> Atomic number
      integer, intent(in) :: num
      !> Whether a valence check and rules exist for the element
      logical :: supported
      select case (num)
      case (1, 3, 6, 7, 8, 9, 15, 16, 17, 35, 53)
         supported = .true.
      case default
         supported = .false.
      end select
   end function supported_element

   !> Restrict charged sites to environments explicitly identified by a rule
   !>
   !> @param[in] name Selected environment label
   !> @param[in] charge Formal charge
   pure function valid_rule_charge(name, charge) result(valid)
      !> Environment label
      character(len=*), intent(in) :: name
      !> Formal charge
      integer, intent(in) :: charge
      !> Whether the rule represents this formal charge
      logical :: valid

      valid = charge == 0
      if (valid) return
      select case (name)
      case ("opls_760", "opls_767", "opls_406")
         valid = charge == 1
      case ("opls_761", "opls_272", "opls_441", "opls_401")
         valid = charge == -1
      end select
   end function valid_rule_charge

   !> Enumerate each simple ring of length 3-6 once
   !>
   !> @param[in,out] graph Ring membership and counts
   !> @param[in,out] path Current simple path, atom indices
   !> @param[in] depth Current path length, 1-6
   recursive subroutine find_cycles(graph, path, depth)
      !> Graph
      type(graph_type), intent(inout) :: graph
      !> Current path
      integer, intent(inout) :: path(6)
      !> Current length
      integer, intent(in) :: depth
      integer :: k, next, current

      current = path(depth)
      do k = 1, graph%degree(current)
         next = graph%neighbor(k, current)
         if (next == path(1)) then
            if (depth < 3) cycle
            if (path(2) >= current) cycle
            graph%cyclic(path(:depth), depth) = .true.
            graph%rings(path(:depth)) = graph%rings(path(:depth)) + 1
         else
            if (depth == 6 .or. next <= path(1)) cycle
            if (any(path(:depth) == next)) cycle
            path(depth + 1) = next
            call find_cycles(graph, path, depth + 1)
         end if
      end do
   end subroutine find_cycles

   !> Compile the bundled connected SMARTS subset
   !>
   !> @param[in] smarts Pattern text
   !> @param[out] pattern Compiled atom predicates and edges
   !> @param[out] valid Whether every token and branch is structurally valid
   subroutine compile_pattern(smarts, pattern, valid)
      !> Pattern text
      character(len=*), intent(in) :: smarts
      !> Compiled pattern
      type(pattern_type), intent(out) :: pattern
      !> Parser status
      logical, intent(out) :: valid
      integer :: pos, last, previous, branch, stack(max_pattern_atoms), ring(0:9), digit
      character(len=96) :: token

      valid = .false.
      pos = 1
      previous = 0
      branch = 0
      ring = 0
      do while (pos <= len_trim(smarts))
         select case (smarts(pos:pos))
         case ("(")
            if (previous == 0 .or. branch == max_pattern_atoms) return
            branch = branch + 1
            stack(branch) = previous
            pos = pos + 1
            cycle
         case (")")
            if (branch == 0) return
            previous = stack(branch)
            branch = branch - 1
            pos = pos + 1
            cycle
         case ("0":"9")
            if (previous == 0) return
            digit = iachar(smarts(pos:pos)) - iachar("0")
            if (ring(digit) == 0) then
               ring(digit) = previous
            else
               pattern%edge(previous, ring(digit)) = .true.
               pattern%edge(ring(digit), previous) = .true.
               ring(digit) = 0
            end if
            pos = pos + 1
            cycle
         case ("[")
            last = index(smarts(pos:), "]")
            if (last <= 2) return
            last = pos + last - 1
            token = smarts(pos + 1:last - 1)
            pos = last + 1
         case ("*")
            token = "*"
            pos = pos + 1
         case ("A":"Z")
            last = pos
            if (pos < len_trim(smarts)) then
               if (smarts(pos + 1:pos + 1) >= "a" .and. smarts(pos + 1:pos + 1) <= "z") last = pos + 1
            end if
            token = smarts(pos:last)
            pos = last + 1
         case default
            return
         end select
         if (pattern%nat == max_pattern_atoms) return
         pattern%nat = pattern%nat + 1
         pattern%expression(pattern%nat) = token
         if (previous > 0) then
            pattern%edge(previous, pattern%nat) = .true.
            pattern%edge(pattern%nat, previous) = .true.
         end if
         previous = pattern%nat
      end do
      valid = pattern%nat > 0 .and. branch == 0 .and. all(ring == 0)
   end subroutine compile_pattern

   !> Match a connected atom environment with its first atom fixed
   !>
   !> @param[in] pattern Compiled environment
   !> @param[in] graph Molecular graph
   !> @param[in] root Structure atom assigned to the pattern's first atom
   !> @param[in] known Matching type-reference predicates, (nrules, nat)
   pure function matches_pattern(pattern, graph, root, known) result(matched)
      !> Environment
      type(pattern_type), intent(in) :: pattern
      !> Graph
      type(graph_type), intent(in) :: graph
      !> Root atom
      integer, intent(in) :: root
      !> Type-reference matches
      logical, intent(in) :: known(:, :)
      !> Whether an injective embedding exists
      logical :: matched
      integer :: map(max_pattern_atoms)

      map = 0
      map(1) = root
      matched = matches_expression(trim(pattern%expression(1)), graph, root, known)
      if (matched) matched = extend_match(pattern, graph, known, map, 2)
   end function matches_pattern

   !> Extend an injective environment embedding over adjacent atoms
   !>
   !> @param[in] pattern Compiled environment
   !> @param[in] graph Molecular graph
   !> @param[in] known Type-reference matches
   !> @param[in] map Structure atom per pattern atom
   !> @param[in] node Next pattern atom to assign
   recursive pure function extend_match(pattern, graph, known, map, node) result(matched)
      !> Environment
      type(pattern_type), intent(in) :: pattern
      !> Graph
      type(graph_type), intent(in) :: graph
      !> Type-reference matches
      logical, intent(in) :: known(:, :)
      !> Current embedding
      integer, intent(in) :: map(max_pattern_atoms)
      !> Next pattern atom
      integer, intent(in) :: node
      !> Whether a complete embedding exists
      logical :: matched
      integer :: parent, k, candidate, earlier
      integer :: next_map(max_pattern_atoms)
      logical :: valid

      matched = .true.
      if (node > pattern%nat) return
      matched = .false.
      parent = findloc(pattern%edge(:node - 1, node), .true., dim=1)
      if (parent == 0) return
      do k = 1, graph%degree(map(parent))
         candidate = graph%neighbor(k, map(parent))
         if (any(map(:node - 1) == candidate)) cycle
         if (.not. matches_expression(trim(pattern%expression(node)), graph, candidate, known)) cycle
         valid = .true.
         do earlier = 1, node - 1
            if (.not. pattern%edge(earlier, node)) cycle
            if (.not. bonded(graph, map(earlier), candidate)) valid = .false.
         end do
         if (.not. valid) cycle
         next_map = map
         next_map(node) = candidate
         if (extend_match(pattern, graph, known, next_map, node + 1)) then
            matched = .true.
            return
         end if
      end do
   end function extend_match

   !> Evaluate atom predicates with low-precedence AND and then OR
   !>
   !> @param[in] expression Predicate expression without brackets
   !> @param[in] graph Molecular graph
   !> @param[in] atom Structure atom
   !> @param[in] known Type-reference matches
   recursive pure function matches_expression(expression, graph, atom, known) result(matched)
      !> Predicate expression
      character(len=*), intent(in) :: expression
      !> Graph
      type(graph_type), intent(in) :: graph
      !> Atom
      integer, intent(in) :: atom
      !> Type-reference matches
      logical, intent(in) :: known(:, :)
      !> Predicate result
      logical :: matched
      integer :: split, value, stat

      matched = .false.
      if (len(expression) == 0) return
      split = index(expression, ";")
      if (split > 0) then
         matched = matches_expression(expression(:split - 1), graph, atom, known) .and. &
            & matches_expression(expression(split + 1:), graph, atom, known)
         return
      end if
      split = index(expression, ",")
      if (split > 0) then
         matched = matches_expression(expression(:split - 1), graph, atom, known) .or. &
            & matches_expression(expression(split + 1:), graph, atom, known)
         return
      end if
      select case (expression(1:1))
      case ("!")
         matched = .not. matches_expression(expression(2:), graph, atom, known)
      case ("*")
         matched = expression == "*"
      case ("%")
         value = rule_index(expression(2:))
         if (value > 0) matched = known(value, atom)
      case ("X", "r", "R", "#")
         read (expression(2:), *, iostat=stat) value
         if (stat /= 0) return
         select case (expression(1:1))
         case ("X")
            matched = graph%degree(atom) == value
         case ("r")
            if (value >= 3 .and. value <= 6) matched = graph%cyclic(atom, value)
         case ("R")
            matched = graph%rings(atom) == value
         case ("#")
            matched = graph%num(atom) == value
         end select
      case default
         matched = graph%num(atom) == to_number(expression)
      end select
   end function matches_expression

   !> Find a bundled rule by its exact name
   !>
   !> @param[in] name Rule label
   pure function rule_index(name) result(index)
      !> Rule label
      character(len=*), intent(in) :: name
      !> Rule index, zero when absent
      integer :: index
      integer :: i

      index = 0
      do i = 1, size(oplsaa_rules)
         if (oplsaa_rules(i)%name /= name) cycle
         index = i
         return
      end do
   end function rule_index

   !> Whether two structure atoms share an explicit bond
   !>
   !> @param[in] graph Molecular graph
   !> @param[in] i First atom
   !> @param[in] j Second atom
   pure function bonded(graph, i, j) result(found)
      !> Graph
      type(graph_type), intent(in) :: graph
      !> First atom
      integer, intent(in) :: i
      !> Second atom
      integer, intent(in) :: j
      !> Whether the bond exists
      logical :: found

      found = any(graph%neighbor(:graph%degree(i), i) == j)
   end function bonded

end module moist_model_moz_potential_lj_oplsaa_typing
