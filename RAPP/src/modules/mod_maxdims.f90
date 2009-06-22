!------------------------------------------------------------------------------------------!
!  Module maxdims: contains parameters with the maximum dimensions allowed.                !
!------------------------------------------------------------------------------------------!
module mod_maxdims
   implicit none

   !----- Maximum string length -----------------------------------------------------------!
   integer, parameter :: maxstr   = 256
   
   !----- Maximum number of files that RAPP can attempt to work ---------------------------!
   integer, parameter :: maxfiles = 2000

   !----- Maximum number of dimensions accepted by RAPP. ----------------------------------!
   integer, parameter :: maxrank  = 6

   !----- Maximum number of grids accepted by RAPP. ---------------------------------------!
   integer, parameter :: maxgrds  = 6

   !----- Maximum number of points --------------------------------------------------------!
   integer, parameter :: nxpmax   = 400 ! Maximum number of points in the horizontal
   integer, parameter :: nypmax   = 400 ! Maximum number of points in the horizontal
   integer, parameter :: nzpmax   = 1   ! Maximum number of points in the vertical

   !----- Maximum number of times per file ------------------------------------------------!
   integer, parameter :: maxtimes = 1489 ! Files every half an hour for a 31-day month.

end module mod_maxdims
