module TestLinalg
   use stdlib_math, only: all_close
   use stdlib_linalg, only: solve
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use constants
   use linalg
   implicit none(type, external)
   private
   public :: collect_dense_linalg_tests
   public :: collect_tridiag_linalg_tests

   integer, parameter :: n = 10
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
      real(wp) :: A(n, n), x(n)
      real(wp) :: y(n), z(n)

      ! Random matix and vector.
      call random_number(A); call random_number(x)

      ! Reference implementation.
      y = matmul(A, x)

      ! Student implementation.
      z = dense_matvec(A, x)

      ! Check that the two match.
      call check(error, all_close(y, z))
      if (allocated(error)) return
   end subroutine

   subroutine test_gaussian_elimination(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: A(n, n), b(n)
      real(wp) :: x(n), y(n)

      ! Random matix and vector.
      call random_number(A); call random_number(b)

      ! Reference implementation.
      x = solve(A, b)

      ! Student implementation.
      y = dense_solve(A, b)

      ! Check that the two match.
      call check(error, all_close(x, y))
      if (allocated(error)) return
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
      type(tridiagonal) :: A
      real(wp) :: coefs(3)
      real(wp) :: x(n), y(n), z(n)

      ! Random matix and vector.
      call random_number(coefs); call random_number(x)
      A = tridiagonal(coefs(1), coefs(2), coefs(3), n)

      ! Reference implementation.
      y = matmul(dense(A), x)

      ! Student implementation.
      z = tridiag_matvec(A, x)

      ! Check that the two match.
      call check(error, all_close(y, z))
      if (allocated(error)) return
   end subroutine

   subroutine test_tridiag_solver(error)
      type(error_type), allocatable, intent(out) :: error
      type(tridiagonal) :: A
      real(wp) :: coefs(3)
      real(wp) :: x(n), y(n), b(n)

      ! Random matix and vector.
      call random_number(coefs); call random_number(b)
      A = tridiagonal(coefs(1), coefs(2), coefs(3), n)

      ! Reference implementation.
      x = solve(dense(A), b)

      ! Student implementation.
      y = tridiag_solve(A, b)

      ! Check that the two match.
      call check(error, all_close(x, y))
      if (allocated(error)) return
   end subroutine

end module
