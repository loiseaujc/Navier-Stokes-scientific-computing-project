module constants
   implicit none(external)
   private
   ! Working precision for the floating point numbers.
   integer, parameter, public :: wp = selected_real_kind(15)
end module
