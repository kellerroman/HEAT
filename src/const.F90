module const
   implicit none
   public
   save
   integer, parameter :: dp = KIND(1.0D+0) ! DOUBLE PRECISION
   integer, parameter :: ip = KIND(1)      ! INTEGER PRECISION

   integer, parameter ::  WEST_FACE          = 1
   integer, parameter ::  EAST_FACE          = 2
   integer, parameter :: SOUTH_FACE          = 3
   integer, parameter :: NORTH_FACE          = 4
   integer, parameter ::  BACK_FACE          = 5
   integer, parameter :: FRONT_FACE          = 6

   integer, parameter ::    BC_ISOTHERMAL    = -1
   integer, parameter ::    BC_HEATFLUX      = -2
   integer, parameter ::    BC_SYMMETRY      = -3


end module const
