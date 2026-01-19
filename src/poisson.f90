module poisson
   use stdlib_optval, only: optval
   use constants
   implicit none
   private

   real(wp), parameter :: eps = epsilon(1.0_wp)
   real(wp), parameter :: tol = sqrt(eps)

   !----------------------------------------
   !-----     HIGH-LEVEL INTERFACE     -----
   !----------------------------------------

   interface
      module function poisson_solver(b, dx, dy, maxiter, method, omega) result(u)
         implicit none
         real(wp), intent(in) :: b(:, :)
         !> Right-hand side field.
         integer, intent(in) :: maxiter
         !> Maximum number of iterations before stopping.
         real(wp), intent(in) :: dx, dy
         !> Grid spacing in each direction.
         character(len=*), intent(in) :: method
         !> Which solver to use. Choose from
         !    - jacobi
         !    - gauss-seidel
         !    - sor (Gauss-Seidel with Successive Over-Relaxation)
         !    - cg  (Conjugate Gradient)
         real(wp), optional, intent(in) :: omega
         !> Over-relaxation weight for SOR. (optional)
         real(wp), allocatable :: u(:, :)
         !> Solution of the Poisson equation.
      end function poisson_solver
   end interface
   public :: poisson_solver

contains

   !-------------------------------------
   !-----     UTILITY FUNCTIONS     -----
   !-------------------------------------

   subroutine swap(u, v)
      real(wp), allocatable, intent(inout) :: u(:, :), v(:, :)
      real(wp), allocatable :: tmp(:, :)
      call move_alloc(from=u, to=tmp)
      call move_alloc(from=v, to=u)
      call move_alloc(from=tmp, to=v)
   end subroutine swap

   !-------------------------------------
   !-----     HIGH-LEVEL SOLVER     -----
   !-------------------------------------

   module procedure poisson_solver
   !----- Interval variables -----
   integer :: nx, ny
   !> Number of grid points in each direction.
   real(wp) :: omega_
   !> Over-relaxation weight.
   !----- Instructions -----
   !> Sanity checks.
   nx = size(b, 1); ny = size(b, 2)
   !> Dispatch to the appropriate solver.
   select case (method)
   case ("jacobi")
      call jacobi_solver(nx, ny, dx, dy, maxiter, b, u)
   case ("gauss-seidel")
      call gauss_seidel_solver(nx, ny, dx, dy, maxiter, b, u)
   case ("sor")
      omega_ = optval(omega, 1.0_wp)
      call sor(nx, ny, dx, dy, maxiter, omega_, b, u)
   case ("cg")
      call cg(nx, ny, dx, dy, maxiter, b, u)
   case default
      error stop "Selected algorithm is not available."
   end select
   end procedure poisson_solver

   !---------------------------------
   !-----     JACOBI METHOD     -----
   !---------------------------------

   subroutine jacobi_solver(nx, ny, dx, dy, maxiter, b, u)
      integer, intent(in) :: nx, ny
      !> Number of grid points in each direction.
      real(wp), intent(in) :: dx, dy
      !> Grid spacing in each direction.
      integer, intent(in) :: maxiter
      !> Maximum number of iteration.
      real(wp), dimension(nx, ny), intent(in) :: b
      !> Right-hand side field.
      real(wp), allocatable, dimension(:, :), intent(out) :: u
      !> Solution of the Poisson equation.

      !----- Internal variables -----
      real(wp), allocatable :: v(:, :)
      !> Temporary solution.
      real(wp) :: l2_norm
      !> Residual norm.
      integer :: iteration
      !> Iteration counter.

      !----- Instructions -----
      !> Initialize variables.
      allocate (u(nx, ny), v(nx, ny), source=0.0_wp)
      iteration = 0; l2_norm = huge(1.0_wp)
      !> Jacobi iterations.
      do while ((iteration < maxiter) .and. (l2_norm > tol))
         !> Jacobi update.
         call jacobi_kernel(nx, ny, dx, dy, b, u, v)
         !> Get the l2-norm of the residual.
         l2_norm = v(1, 1)
         !> Swap u and v for the next itertion.
         call swap(u, v)
      end do

      print *, "    Jacobi solver :"
      print *, "         - Number of iterations    :", iteration
      print *, "         - l2-norm of the residual :", l2_norm
   end subroutine jacobi_solver

   subroutine jacobi_kernel(nx, ny, dx, dy, b, u, v)
      integer, intent(in) :: nx, ny
      !> Number of grid points in each direction.
      real(wp), intent(in) :: dx, dy
      !> Grid spacing in each direction.
      real(wp), dimension(nx, ny), intent(in) :: b
      !> Right-hand side field.
      real(wp), dimension(nx, ny), intent(in) :: u
      !> Current approximate solution.
      real(wp), dimension(nx, ny), intent(out) :: v
      !> New approximate solution.

      !----- Internal variables -----
      integer :: i, j
      !> Loop counters.
      real(wp) :: l2_norm
      !> Residual norm.
      real(wp) :: sx, sy, sb
      !> Miscellaneous.

      !----- Instructions -----
      sx = 0.5_wp*dx**2/(dx**2 + dy**2)
      sy = 0.5_wp*dy**2/(dx**2 + dy**2)
      sb = 0.5_wp*dx**2*dy**2/(dx**2 + dy**2)
      l2_norm = 0.0_wp
      do j = 2, ny - 1
         do i = 2, nx - 1
            !> Jacobi update rule.
            v(i, j) = sy*(u(i + 1, j) + u(i - 1, j)) &
                      + sx*(u(i, j + 1) + u(i, j - 1)) &
                      - sb*b(i, j)
            !> Residual norm.
            l2_norm = l2_norm + dx*dy*(u(i, j) - v(i, j))**2
         end do
      end do
      !> l2-norm is stored in v(1, 1) for the sake of storage facility.
      v(1, 1) = sqrt(l2_norm)
   end subroutine jacobi_kernel

   !---------------------------------------
   !-----     GAUSS-SEIDEL METHOD     -----
   !---------------------------------------

   subroutine gauss_seidel_solver(nx, ny, dx, dy, maxiter, b, u)
      integer, intent(in) :: nx, ny
      !> Number of grid points in each direction.
      real(wp), intent(in) :: dx, dy
      !> Grid spacing in each direction.
      integer, intent(in) :: maxiter
      !> Maximum number of iterations.
      real(wp), dimension(nx, ny), intent(in) :: b
      !> Right-hand side field.
      real(wp), allocatable, dimension(:, :), intent(out) :: u
      !> Solution of the Poisson equation.

      !----- Internal variables -----

      !----- Instructions -----
      !> Initialize variables.
      allocate (u(nx, ny), source=0.0_wp)
   end subroutine gauss_seidel_solver

   !-----------------------------------------------------
   !-----     SUCCESSIVE OVER-RELAXATION METHOD     -----
   !-----------------------------------------------------

   subroutine sor(nx, ny, dx, dy, maxiter, omega, b, u)
      integer, intent(in) :: nx, ny
      !> Number of grid points in each direction.
      real(wp), intent(in) :: dx, dy
      !> Grid spacing in each direction.
      integer, intent(in) :: maxiter
      !> Maximum number of iterations.
      real(wp), intent(in) :: omega
      !> Over-relaxation weight.
      real(wp), dimension(nx, ny), intent(in) :: b
      !> Right-hand side field.
      real(wp), allocatable, dimension(:, :), intent(out) :: u
      !> Solution of the Poisson equation.

      !----- Internal variables -----

      !----- Instructions -----
      !> Initialize variables.
      allocate (u(nx, ny), source=0.0_wp)
   end subroutine sor

   !----------------------------------------------
   !-----      CONJUGATE GRADIENT METHOD     -----
   !----------------------------------------------

   subroutine cg(nx, ny, dx, dy, maxiter, b, u)
      integer, intent(in) :: nx, ny
      !> Number of grid points in each direction.
      real(wp), intent(in) :: dx, dy
      !> Grid spacing in each direction.
      integer, intent(in) :: maxiter
      !> Maximum number of iterations.
      real(wp), dimension(nx, ny), intent(in) :: b
      !> Right-hand side field.
      real(wp), allocatable, dimension(:, :), intent(out) :: u
      !> Solution of the Poisson equation.

      !----- Internal variables -----

      !----- Instructions -----
      !> Initialize variables.
      allocate (u(nx, ny), source=0.0_wp)
   end subroutine cg

end module
