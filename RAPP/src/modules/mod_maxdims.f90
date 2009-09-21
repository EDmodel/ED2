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
   integer, parameter :: maxgrds  = 13

   !----- Maximum number of points --------------------------------------------------------!
   integer, parameter :: nxpmax   = 400 ! Maximum number of points in the horizontal
   integer, parameter :: nypmax   = 400 ! Maximum number of points in the horizontal
   integer, parameter :: nzpmax   = 1   ! Maximum number of points in the vertical

   !----- Maximum number of times per day and file ----------------------------------------!
   integer, parameter :: maxperday = 48             ! Data every ½ hour
   integer, parameter :: maxtimes  = maxperday * 31 ! Data every ½ hour, 31-day month.
   integer, parameter :: maxpdf    = maxperday * 12 ! Data every ½ hour, 12 months.

end module mod_maxdims
