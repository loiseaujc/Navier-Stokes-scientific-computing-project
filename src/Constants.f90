module Constants
   implicit none(type, external)
   private

   ! Standard kind for double precision arithmetic.
   integer, parameter, public :: dp = selected_real_kind(15)
   ! Absolute tolerance.
   real(dp), parameter, public :: eps = epsilon(1.0_dp)
   ! Relative tolerance.
   real(dp), parameter, public :: tol = sqrt(eps)
end module Constants
