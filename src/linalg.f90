module linalg
   use constants
   implicit none(type, external)
   private

   public :: dense_matvec
   public :: dense_solve
   public :: tridiagonal, dense
   public :: tridiag_matvec
   public :: tridiag_solve

   !----------------------------------------
   !-----     DENSE LINEAR ALGEBRA     -----
   !----------------------------------------

   interface
      !! Function to perform the matrix-vector product y = Ax where
      !! A is an arbitrary m x n dense matrix.
      module function dense_matvec(A, x) result(y)
         real(wp), intent(in)  :: A(:, :)
         !! Arbitrary m x n dense matrix.
         real(wp), intent(in)  :: x(:)
         !! n-vector.
         real(wp), allocatable :: y(:)
         !! m-vector.
      end function
   end interface

   interface
      !! Function to solve the linear system Ax = b using Gaussian
      !! elimination.
      module function dense_solve(A, b) result(x)
         real(wp), intent(in)  :: A(:, :)
         !! n x n dense matrix describing the linear system.
         real(wp), intent(in)  :: b(:)
         !! Right-hand side n-vector.
         real(wp), allocatable :: x(:)
         !! Solution n-vector.
      end function
   end interface

   !----------------------------------------------
   !-----     TRIDIAGONAL MATRIX ALGEBRA     -----
   !----------------------------------------------

   ! Derived-type (i.e. object) to represent a tridiagonal matrix
   ! with constant coefficients.
   type :: tridiagonal
      real(wp) :: a, b, c
      !! Lower-diagonal, diagonal and upper-diagonal entries of
      !! the tridiagonal matrix.
      integer  :: n
      !! Dimension of the matrix.
   end type

   interface
      module function dense(A) result(B)
         type(tridiagonal), intent(in) :: A
         real(wp), allocatable :: B(:, :)
      end function
   end interface

   interface
      module function tridiag_matvec(A, x) result(y)
         type(tridiagonal), intent(in) :: A
         real(wp), intent(in)  :: x(:)
         real(wp), allocatable :: y(:)
      end function
   end interface

   interface
      module function tridiag_solve(A, b) result(x)
         type(tridiagonal), intent(in) :: A
         real(wp), intent(in)  :: b(:)
         real(wp), allocatable :: x(:)
      end function
   end interface

contains

   !----------------------------------------
   !-----     DENSE LINEAR ALGEBRA     -----
   !----------------------------------------

   module procedure dense_matvec
   integer :: n
   ! Vector dimension + memory allocation for the result.
   n = size(A, 1); allocate (y(n), source=0.0_wp)
   end procedure

   module procedure dense_solve
   integer :: n
   ! Vector dimension + memory allocation for the result.
   n = size(b); allocate (x(n), source=0.0_wp)
   end procedure

   !----------------------------------------------
   !-----     TRIDIAGONAL MATRIX ALGEBRA     -----
   !----------------------------------------------

   module procedure dense
   integer :: n, i
   n = A%n; allocate (B(n, n), source=0.0_wp)
   B(1, 1) = A%b; B(1, 2) = A%c
   do i = 2, n - 1
      B(i, i - 1) = A%a
      B(i, i) = A%b
      B(i, i + 1) = A%c
   end do
   B(n, n - 1) = A%a; B(n, n) = A%b
   end procedure

   module procedure tridiag_matvec
   integer :: n
   ! Vector dimension + memory allocation for the result.
   n = A%n; allocate (y(n), source=0.0_wp)
   end procedure

   module procedure tridiag_solve
   integer :: n
   ! Vector dimension + memory allocation for the result.
   n = size(b); allocate (x(n), source=0.0_wp)
   end procedure
end module linalg
