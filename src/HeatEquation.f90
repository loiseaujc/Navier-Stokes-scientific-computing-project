module HeatEquation
   use stdlib_specialmatrices, only: tridiagonal, &         ! Constructor
                                     tridiagonal_dp_type, & ! Data type
                                     spmv                   ! Matrix-vector product
   use stdlib_linalg, only: norm
   use stdlib_math, only: linspace
   use stdlib_stats, only: mean, var
   use constants
   implicit none(type, external)
   private

   !------------------------------------------
   !-----     LINEAR ALGEBRA SOLVERS     -----
   !------------------------------------------

   interface matmul
      !! Utility wrapper for the tridiagonal matrix-vector product.
      module function tridiag_matmul(A, x) result(y)
         type(tridiagonal_dp_type), intent(in) :: A
         real(dp), intent(in) :: x(:)
         real(dp), allocatable :: y(:)
      end function
   end interface

   interface
      !! Function for solving a square linear system of equation Ax = b
      !! using Gaussian elimination. The matrix A is an arbitrary
      !! square matrix of size n x n.
      pure module function dense_solver(A, b) result(x)
         real(dp), intent(in) :: A(:, :)
         !! n x n Matrix describing the system of equation.
         real(dp), intent(in) :: b(:)
         !! Right-hand side n-vector.
         real(dp), allocatable :: x(:)
         !! Solution n-vector.
      end function
   end interface
   public :: dense_solver

   interface
      !! Function for solving a square linear system of equation Ax = b
      !! where A is a tridiagonal matrix of size n x n. It uses the
      !! Thomas algorithm.
      pure module function tridiag_solver(A, b) result(x)
         type(tridiagonal_dp_type), intent(in) :: A
         !! n x n tridiagonal matrix describing the system of equation.
         real(dp), intent(in) :: b(:)
         !! Right-hand side n-vector.
         real(dp), allocatable :: x(:)
         !! Solution n-vector.
      end function
   end interface
   public :: tridiag_solver

   !--------------------------------------------------------
   !-----     ONE-DIMENSIONAL STEADY HEAT EQUATION     -----
   !--------------------------------------------------------

   interface
      pure module function steadyheat_grid_convergence(ns) result(err)
         integer, intent(in) :: ns(:)
         !! Number of points for the list of problems to be tested.
         real(dp), allocatable :: err(:, :)
         !! Errors (in different norms) between the numerical and
         !! analytical solutions.
         !!    - err(:, 1) : Error in the L2 norm.
         !!    - err(:, 2) : Error in the L1 norm.
         !!    - err(:, 3) : Error in the L-infinity norm.
         !! Note: The first dimension of err is equal to the length of ns.
      end function
   end interface
   public :: steadyheat_grid_convergence

   interface
      pure module function steadyheat_computational_performances(ns, ntrials, solver) result(stats)
         integer, intent(in) :: ns(:)
         !! Number of points for the list of problems to be tested.
         integer, intent(in) :: ntrials
         !! Number of times each solve is done to obtain converged statisics.
         character(len=*), intent(in) :: solver
         !! Which linear solver needs to be tested:
         !!    - solver = "tridiag" : Thomas algorithm.
         !!    - solver = "dense"   : Dense linear solver.
         real(dp), allocatable :: stats(:, :)
         !! Statistics of the different solvers.
         !!    - stats(:, 1) : Average time needed to solve the problems.
         !!    - stats(:, 2) : Standard deviation for each problem size.
         !!    - stats(:, 3) : Minimum time for each problem size.
         !!    - stats(:, 4) : Maximum time for each problem size.
         !! Note: The first dimension of stats is equal to the length of ns.
      end function
   end interface
   public :: steadyheat_computational_performances

   !----------------------------------------------------------
   !-----     ONE-DIMENSIONAL UNSTEADY HEAT EQUATION     -----
   !----------------------------------------------------------

   interface
      pure module function unsteadyheat_dt_convergence(dts, Tmax, n) result(err)
         real(dp), intent(in) :: dts(:)
         !! Range of time steps considered for the analysis.
         real(dp), intent(in) :: Tmax
         !! End time for the temporal integration.
         integer, intent(in) :: n
         !! Number of grid points used to discretize the unit segment.
         real(dp), allocatable :: err(:, :)
         !! Errors (in different norms) between the numerical and
         !! analytical solutions.
         !!    - err(:, 1) : Error in the L2 norm.
         !!    - err(:, 2) : Error in the L1 norm.
         !!    - err(:, 3) : Error in the L-infinity norm.
         !! Note: The error needs to be evaluated at time t = Tmax only.
      end function
   end interface
   public :: unsteadyheat_dt_convergence

contains

   !------------------------------------------
   !-----     LINEAR ALGEBRA SOLVERS     -----
   !------------------------------------------

   module procedure tridiag_matmul
   integer :: n
   n = size(x); allocate (y(n), source=0.0_dp)
   call spmv(A, x, y)
   end procedure

   module procedure dense_solver
   end procedure

   module procedure tridiag_solver
   end procedure

   !--------------------------------------------------------
   !-----     ONE-DIMENSIONAL STEADY HEAT EQUATION     -----
   !--------------------------------------------------------

   module procedure steadyheat_grid_convergence
   end procedure

   module procedure steadyheat_computational_performances
   end procedure

   real(dp) elemental function analytic_rhs(x) result(f)
      real(dp), intent(in) :: x
        !! Spatial coordinate where the rhs vector needs to be evaluated.
   end function

   real(dp) elemental function analytic_steady_solution(x) result(u)
      real(dp), intent(in) :: x
      !! Spatial coordinate where the analytic solution needs to be evaluated.
   end function

   real(dp) function std(x) result(out)
      real(dp), intent(in) :: x(:)
      !! Array whose the standard deviation will be computed.
      out = sqrt(var(x))
   end function

   !----------------------------------------------------------
   !-----     ONE-DIMENSIONAL UNSTEADY HEAT EQUATION     -----
   !----------------------------------------------------------

   module procedure unsteadyheat_dt_convergence
   end procedure

   function explicit_euler_step(u, f, L, dt) result(v)
      real(dp), intent(in) :: u(:)
      !! Current solution at time t = k*dt.
      real(dp), intent(in) :: f(:)
      !! Constant forcing term.
      type(tridiagonal_dp_type), intent(in) :: L
      !! Tridiagonal matrix approximation of the 1D Laplacian operator.
      real(dp), intent(in) :: dt
      !! Time step.
      real(dp), allocatable :: v(:)
      !! New solution at time t = (k+1)*dt.
   end function

   function implicit_euler_step(u, f, H, dt) result(v)
      real(dp), intent(in) :: u(:)
      !! Current solution at time t = k*dt.
      real(dp), intent(in) :: f(:)
      !! Constant forcing term.
      type(tridiagonal_dp_type), intent(in) :: H
      !! 1D Helmholtz operator.
      real(dp), intent(in) :: dt
      !! Time step.
      real(dp), allocatable :: v(:)
      !! New solution at time t = (k+1)*dt.
   end function

   function crank_nicolson(u, f, H, L, dt) result(v)
      real(dp), intent(in) :: u(:)
      !! Current solution at time t = k*dt.
      real(dp), intent(in) :: f(:)
      !! Constant forcing term.
      type(tridiagonal_dp_type), intent(in) :: H
      !! 1D Helmholtz operator.
      type(tridiagonal_dp_type), intent(in) :: L
      !! 1D Laplace operator.
      real(dp), intent(in) :: dt
      !! Time step.
      real(dp), allocatable :: v(:)
      !! New solution at time t = (k+1)*dt.
   end function

   real(dp) elemental function analytic_unsteady_solution(x, t) result(u)
      real(dp), intent(in) :: x, t
      !! Spatio-temporal coordinates where to evaluate the analytical solution.
   end function
end module HeatEquation
