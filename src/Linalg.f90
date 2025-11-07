module Linalg
   use constants
   implicit none(type, external)
   private

   !----------------------------------------
   !-----     DENSE LINEAR ALGEBRA     -----
   !----------------------------------------

   interface
      !! Function for solving a square linear system of equation Ax = b
      !! using Gaussian elimination. The matrix A is an arbitrary
      !! square matrix of size n x n.
      module function dense_solver(A, b) result(x)
         real(dp), intent(in) :: A(:, :)
         !! n x n Matrix describing the system of equation.
         real(dp), intent(in) :: b(:)
         !! Right-hand side n-vector.
         real(dp), allocatable :: x(:)
         !! Solution n-vector.
      end function
   end interface
   public :: dense_solver

   !----------------------------------------
   !-----     TRIDIAGONAL MATRICES     -----
   !----------------------------------------

   ! Data structure.
   type, public :: tridiagonal_type
      private
      real(dp), allocatable :: dl(:), dv(:), du(:)
   end type

   ! Constructor.
   interface tridiagonal
      pure module function constant_tridiag(a, b, c, n) result(Amat)
         real(dp), intent(in) :: a, b, c
         integer, intent(in) :: n
         type(tridiagonal_type) :: Amat
      end function
   end interface
   public :: tridiagonal

   ! Matrix-vector product.
   interface matmul
      module function tridiag_matmul(A, x) result(y)
         type(tridiagonal_type), intent(in) :: A
         real(dp), intent(in) :: x(:)
         real(dp), allocatable :: y(:)
      end function
   end interface
   public :: matmul

   interface
      !! Function for solving a square linear system of equation Ax = b
      !! where A is a tridiagonal matrix of size n x n. It uses the
      !! Thomas algorithm.
      module function tridiag_solver(A, b) result(x)
         type(tridiagonal_type), intent(in) :: A
         !! n x n tridiagonal matrix describing the system of equation.
         real(dp), intent(in) :: b(:)
         !! Right-hand side n-vector.
         real(dp), allocatable :: x(:)
         !! Solution n-vector.
      end function
   end interface
   public :: tridiag_solver
contains

   !----------------------------------------
   !-----     DENSE LINEAR ALGEBRA     -----
   !----------------------------------------

   module procedure dense_solver
   end procedure

   !----------------------------------------
   !-----     TRIDIAGONAL MATRICES     -----
   !----------------------------------------

   module procedure constant_tridiag
   end procedure

   module procedure tridiag_matmul
   integer :: n
   n = size(x); allocate (y(n), source=0.0_dp)
   end procedure

   module procedure tridiag_solver
   end procedure

end module
