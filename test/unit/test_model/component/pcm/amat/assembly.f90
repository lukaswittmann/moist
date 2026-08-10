!> Assembly-level unit tests for the Gaussian PCM interaction matrix
module test_model_component_pcm_amat_assembly
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_model_component_pcm_amat, only: assemble_pcm_amat
   use test_helpers, only: get_test_structures, center_at_origin, &
                           get_test_cavity_iswig, rel_deviation
   use test_model_component_pcm_amat, only: count_branches, nmol, nleb_survey
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_model_component_pcm_amat_assembly

   !> Denser grid used for the single-structure refinement check
   integer, parameter :: nleb_dense = 110

   !> Relative tolerance for the entry-by-entry reference comparison
   real(wp), parameter :: ref_thr = 1.0e-14_wp

contains

   !> Collect the assembly test suite
   subroutine collect_model_component_pcm_amat_assembly(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("assembly_vs_reference", test_assembly_vs_reference), &
                  new_unittest("assembly_branch_coverage", test_assembly_branch_coverage), &
                  new_unittest("assembly_dense_grid", test_assembly_dense_grid), &
                  new_unittest("coincident_tesserae", test_coincident_tesserae) &
                  ]

   end subroutine collect_model_component_pcm_amat_assembly

   !* --------------------------------- Local helpers --------------------------------- *!

   !> Reference matrix entry evaluated straight from the closed form
   !>
   !> @param[in]  xi   Gaussian widths (ngrid)
   !> @param[in]  f    Gaussian switching factors (ngrid)
   !> @param[in]  xyz  Surface positions (3, ngrid)
   !> @param[in]  i    Row index
   !> @param[in]  j    Column index
   !> @return     a    Reference matrix entry
   pure function reference_entry(xi, f, xyz, i, j) result(a)
      !> Gaussian widths, switching factors and surface positions
      real(wp), intent(in) :: xi(:), f(:), xyz(:, :)
      !> Row and column indices
      integer, intent(in) :: i, j
      !> Reference matrix entry
      real(wp) :: a

      !> Pair width and separation
      real(wp) :: p, r

      if (i == j) then
         a = sqrt(2.0_wp/(4.0_wp*atan(1.0_wp)))*xi(i)/f(i)
      else
         p = xi(i)*xi(j)/sqrt(xi(i)*xi(i) + xi(j)*xi(j))
         r = norm2(xyz(:, i) - xyz(:, j))
         a = erf(p*r)/r
      end if

   end function reference_entry

   !> Largest elementwise relative deviation from the reference matrix
   !>
   !> @param[in]  xi    Gaussian widths (ngrid)
   !> @param[in]  f     Gaussian switching factors (ngrid)
   !> @param[in]  xyz   Surface positions (3, ngrid)
   !> @param[in]  amat  Assembled matrix (ngrid, ngrid)
   !> @return     dev   Largest |a - a_ref|/(1 + |a_ref|)
   function max_reference_deviation(xi, f, xyz, amat) result(dev)
      !> Gaussian widths, switching factors and surface positions
      real(wp), intent(in) :: xi(:), f(:), xyz(:, :)
      !> Assembled matrix
      real(wp), intent(in) :: amat(:, :)
      !> Largest elementwise relative deviation
      real(wp) :: dev

      !> Matrix indices and grid size
      integer :: i, j, ngrid
      !> Reference entry
      real(wp) :: a_ref

      ngrid = size(xi)
      dev = 0.0_wp
      do j = 1, ngrid
         do i = 1, ngrid
            a_ref = reference_entry(xi, f, xyz, i, j)
            dev = max(dev, rel_deviation(amat(i, j), a_ref))
         end do
      end do

   end function max_reference_deviation

   !* ------------------------------------- Tests ------------------------------------- *!

   !> Every assembled entry reproduces the closed-form kernel
   subroutine test_assembly_vs_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Surface of the current structure
      type(cavity_type_iswig) :: cavity
      !> Assembled matrix
      real(wp), allocatable :: amat(:, :)
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Structure and matrix indices, grid size
      integer :: imol, i, j, ngrid
      !> Largest relative deviation and largest asymmetry
      real(wp) :: dev, asym
      !> Failure context
      character(len=128) :: context

      call get_test_structures(mols, nmol)

      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb_survey)
         if (allocated(err)) then
            call test_failed(error, "surface construction failed: "//err%message)
            return
         end if

         ngrid = cavity%ngrid
         allocate (amat(ngrid, ngrid))
         call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat, err)
         if (allocated(err)) then
            call test_failed(error, "assembly failed: "//err%message)
            return
         end if

         dev = max_reference_deviation(cavity%xi0, cavity%f, cavity%xyz, amat)
         write (context, "(a,i0,a,i0,a)") "structure ", imol, " (ngrid = ", ngrid, &
            "): assembled matrix deviates from erf(p*r)/r"
         call check(error, dev, 0.0_wp, thr=ref_thr, more=trim(context))
         if (allocated(error)) return

         ! The lower triangle is written by the blocked mirror rather than by
         ! the kernel, so an indexing slip there shows up only as asymmetry.
         asym = 0.0_wp
         do j = 1, ngrid
            do i = j + 1, ngrid
               asym = max(asym, abs(amat(i, j) - amat(j, i)))
            end do
         end do
         write (context, "(a,i0,a)") "structure ", imol, &
            ": mirrored triangle does not match the assembled one"
         call check(error, asym, 0.0_wp, thr=0.0_wp, more=trim(context))
         if (allocated(error)) return

         deallocate (amat)
      end do

   end subroutine test_assembly_vs_reference

   !> Both kernel branches are actually populated by the survey surfaces
   subroutine test_assembly_branch_coverage(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Surface of the current structure
      type(cavity_type_iswig) :: cavity
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Structure index and per-branch pair counts
      integer :: imol, nfar, nnear
      !> Totals across all sampled structures
      integer :: total_far, total_near
      !> Failure context
      character(len=128) :: context

      call get_test_structures(mols, nmol)

      total_far = 0
      total_near = 0
      do imol = 1, size(mols)
         call center_at_origin(mols(imol))
         call get_test_cavity_iswig(mols(imol), cavity, err, nleb=nleb_survey)
         if (allocated(err)) then
            call test_failed(error, "surface construction failed: "//err%message)
            return
         end if

         call count_branches(cavity%xi0, cavity%xyz, nfar, nnear)
         write (context, "(a,i0,a)") "structure ", imol, &
            ": no unsaturated pairs, the transcendental branch is untested here"
         call check(error, nnear > 0, more=trim(context))
         if (allocated(error)) return

         total_far = total_far + nfar
         total_near = total_near + nnear
      end do

      call check(error, total_far > 0, &
                 more="no saturated pairs at all, the fast branch is untested here")

   end subroutine test_assembly_branch_coverage

   !> The reference agreement survives a much denser Lebedev grid
   subroutine test_assembly_dense_grid(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Dense surface
      type(cavity_type_iswig) :: cavity
      !> Assembled matrix
      real(wp), allocatable :: amat(:, :)
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Grid size
      integer :: ngrid
      !> Largest relative deviation
      real(wp) :: dev

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb_dense)
      if (allocated(err)) then
         call test_failed(error, "surface construction failed: "//err%message)
         return
      end if

      ngrid = cavity%ngrid
      allocate (amat(ngrid, ngrid))
      call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat, err)
      if (allocated(err)) then
         call test_failed(error, "assembly failed: "//err%message)
         return
      end if

      dev = max_reference_deviation(cavity%xi0, cavity%f, cavity%xyz, amat)
      call check(error, dev, 0.0_wp, thr=ref_thr, &
                 more="dense-grid assembly deviates from erf(p*r)/r")

   end subroutine test_assembly_dense_grid

   !> Coincident tesserae stay finite and hit the Gaussian overlap limit
   subroutine test_coincident_tesserae(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Synthetic surface with a duplicated point
      real(wp) :: xi(4), f(4), xyz(3, 4), amat(4, 4)
      !> Pair width, expected limit and matrix indices
      real(wp) :: p, limit
      integer :: i, j

      xi = [1.5_wp, 1.5_wp, 2.0_wp, 0.9_wp]
      f = [0.8_wp, 0.8_wp, 0.95_wp, 0.6_wp]
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [1.7_wp, -0.4_wp, 0.9_wp]
      xyz(:, 4) = [-3.1_wp, 2.2_wp, -1.4_wp]

      call assemble_pcm_amat(xi, f, xyz, amat, err)
      if (allocated(err)) then
         call test_failed(error, "assembly failed: "//err%message)
         return
      end if

      do j = 1, 4
         do i = 1, 4
            call check(error, amat(i, j) == amat(i, j), &
                       more="coincident tesserae produced a NaN")
            if (allocated(error)) return
            call check(error, abs(amat(i, j)) < huge(1.0_wp), &
                       more="coincident tesserae produced an infinity")
            if (allocated(error)) return
         end do
      end do

      p = xi(1)*xi(2)/sqrt(xi(1)*xi(1) + xi(2)*xi(2))
      limit = 2.0_wp*p/sqrt(4.0_wp*atan(1.0_wp))
      call check(error, amat(1, 2), limit, thr=ref_thr*limit, &
                 more="coincident pair does not reach the Gaussian overlap limit")

   end subroutine test_coincident_tesserae

end module test_model_component_pcm_amat_assembly
