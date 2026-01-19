module TestPoisson
   use stdlib_math, only: all_close
   use stdlib_linalg, only: solve
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use constants
   use poisson
   implicit none(type, external)
   private
   public :: collect_poisson_tests

   integer, parameter :: nx = 4
   real(wp), parameter :: dx = 1.0_wp/(nx - 1)
   real(wp), parameter :: eps = epsilon(1.0_wp)
   real(wp), parameter :: tol = nx*sqrt(eps)
contains

   !-------------------------------------------
   !-----     TESTING POISSON SOLVERS     -----
   !-------------------------------------------

   subroutine collect_poisson_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Jacobi solver", test_jacobi_solver), &
                   new_unittest("Gauss-Seidel solver", test_gauss_seidel_solver), &
                   new_unittest("SOR solver", test_sor_solver), &
                   new_unittest("Conjugate Gradient solver", test_cg_solver)]
   end subroutine

   subroutine test_jacobi_solver(error)
      type(error_type), allocatable, intent(out) :: error

      !----- Internal variables -----
      real(wp), allocatable :: u(:, :), b(:, :), r(:, :)
      !> Solution and right-hand side of the Poisson equation.
      integer, parameter :: maxiter = nx**2
      !> Maximum number of iterations.
      character(len=*), parameter :: method = "jacobi"

      !----- Instructions -----
      !> Allocate data.
      allocate (b(nx, nx), source=1.0_wp)
      !> Zero-out boundary points.
      b(1, :) = 0.0_wp; b(nx, :) = 0.0_wp
      b(:, 1) = 0.0_wp; b(:, nx) = 0.0_wp
      !> Solve Poisson equation.
      u = poisson_solver(b, dx, dx, maxiter, method)
      !> Compute residual.
      r = apply_laplacian(u) - b
      !> Check correctness.
      call check(error, norm2(r*dx**2) <= tol)
   end subroutine

   subroutine test_gauss_seidel_solver(error)
      type(error_type), allocatable, intent(out) :: error

      !----- Internal variables -----
      real(wp), allocatable :: u(:, :), b(:, :), r(:, :)
      !> Solution and right-hand side of the Poisson equation.
      integer, parameter :: maxiter = nx**2
      !> Maximum number of iterations.
      character(len=*), parameter :: method = "gauss-seidel"

      !----- Instructions -----
      !> Allocate data.
      allocate (b(nx, nx), source=1.0_wp)
      !> Zero-out boundary points.
      b(1, :) = 0.0_wp; b(nx, :) = 0.0_wp
      b(:, 1) = 0.0_wp; b(:, nx) = 0.0_wp
      !> Solve Poisson equation.
      u = poisson_solver(b, dx, dx, maxiter, method)
      !> Compute residual.
      r = apply_laplacian(u) - b
      !> Check correctness.
      call check(error, norm2(r*dx**2) <= tol)
   end subroutine

   subroutine test_sor_solver(error)
      type(error_type), allocatable, intent(out) :: error

      !----- Internal variables -----
      real(wp), allocatable :: u(:, :), b(:, :), r(:, :)
      !> Solution and right-hand side of the Poisson equation.
      integer, parameter :: maxiter = nx**2
      !> Maximum number of iterations.
      character(len=*), parameter :: method = "sor"

      !----- Instructions -----
      !> Allocate data.
      allocate (b(nx, nx), source=1.0_wp)
      !> Zero-out boundary points.
      b(1, :) = 0.0_wp; b(nx, :) = 0.0_wp
      b(:, 1) = 0.0_wp; b(:, nx) = 0.0_wp
      !> Solve Poisson equation.
      u = poisson_solver(b, dx, dx, maxiter, method)
      !> Compute residual.
      r = apply_laplacian(u) - b
      !> Check correctness.
      call check(error, norm2(r*dx**2) <= tol)
   end subroutine

   subroutine test_cg_solver(error)
      type(error_type), allocatable, intent(out) :: error

      !----- Internal variables -----
      real(wp), allocatable :: u(:, :), b(:, :), r(:, :)
      !> Solution and right-hand side of the Poisson equation.
      integer, parameter :: maxiter = nx**2
      !> Maximum number of iterations.
      character(len=*), parameter :: method = "cg"

      !----- Instructions -----
      !> Allocate data.
      allocate (b(nx, nx), source=1.0_wp)
      !> Zero-out boundary points.
      b(1, :) = 0.0_wp; b(nx, :) = 0.0_wp
      b(:, 1) = 0.0_wp; b(:, nx) = 0.0_wp
      !> Solve Poisson equation.
      u = poisson_solver(b, dx, dx, maxiter, method)
      !> Compute residual.
      r = apply_laplacian(u) - b
      !> Check correctness.
      call check(error, norm2(r*dx**2) <= tol)
   end subroutine

   !-------------------------------------
   !-----     UTILITY FUNCTIONS     -----
   !-------------------------------------

   function apply_laplacian(u) result(v)
      real(wp), intent(in) :: u(:, :)
      real(wp), allocatable :: v(:, :)
      !----- Internal variables -----
      integer :: i, j
      !----- Instructions -----
      allocate (v(nx, nx), source=0.0_wp)
      do concurrent(i=2:nx - 1, j=2:nx - 1)
         v(i, j) = (u(i + 1, j) - 2*u(i, j) + u(i - 1, j))/dx**2 &
                   + (u(i, j + 1) - 2*u(i, j) + u(i, j - 1))/dx**2
      end do
   end function apply_laplacian

end module
