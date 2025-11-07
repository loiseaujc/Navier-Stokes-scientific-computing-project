module TestLinalg
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private
   public :: collect_dense_linalg_tests
   public :: collect_tridiag_linalg_tests
contains

   !-----------------------------------------------------
   !-----     LINEAR ALGEBRA FOR DENSE MATRICES     -----
   !-----------------------------------------------------

   subroutine collect_dense_linalg_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Matrix-vector product", test_dense_matrix_vector_product)]
      testsuite = [testsuite, new_unittest("Gaussian elimination", test_gaussian_elimination)]
   end subroutine

   subroutine test_dense_matrix_vector_product(error)
      type(error_type), allocatable, intent(out) :: error
   end subroutine

   subroutine test_gaussian_elimination(error)
      type(error_type), allocatable, intent(out) :: error
   end subroutine

   !-----------------------------------------------------------
   !-----     LINEAR ALGEBRA FOR TRIDIAGONAL MATRICES     -----
   !-----------------------------------------------------------

   subroutine collect_tridiag_linalg_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Matrix-vector product", test_tridiag_matrix_vector_product)]
      testsuite = [testsuite, new_unittest("Thomas algorithm", test_tridiag_solver)]
   end subroutine

   subroutine test_tridiag_matrix_vector_product(error)
      type(error_type), allocatable, intent(out) :: error
   end subroutine

   subroutine test_tridiag_solver(error)
      type(error_type), allocatable, intent(out) :: error
   end subroutine

end module
