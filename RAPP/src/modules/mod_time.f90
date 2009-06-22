!==========================================================================================!
!==========================================================================================!
! MEVI. Module time. This module contains the time structure and some time constants.      !
!==========================================================================================!
!==========================================================================================!
module mod_time
   type time_stt
      integer           :: year       ! Year
      integer           :: month      ! Month
      integer           :: day        ! Day
      integer           :: hour       ! Hour
      integer           :: minu       ! Minute
      integer           :: seco       ! Second
      real              :: fracday    ! Fraction of one day
      integer           :: doy        ! Day of year
      real(kind=8)      :: elapsed    ! Elapsed time since MEVI origin (Jan 1 1583, 00 GMT)
      character(len=3)  :: mmm        ! 3-letter month
      integer           :: yyyyddd    ! Vis5D date stamp (year//day of year)
      integer           :: hhmmss     ! Vis5D time stamp (hour//minute//second)
      character(len=15) :: gradsstamp ! Date/time stamp for GrADS
      character(len=17) :: timestr    ! Date/time stamp for MEVI output name
   end type time_stt

   character(len=3), dimension(12), parameter  :: monnames = (/'JAN','FEB','MAR','APR'     &
                                                              ,'MAY','JUN','JUL','AUG'     &
                                                              ,'SEP','OCT','NOV','DEC'/)


   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   type(time_stt) function no_time(nothing_int,nothing_real,nothing_dble,nothing_char)
      use mod_maxdims, only : maxstr
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=maxstr), intent(in) :: nothing_char
      integer              , intent(in) :: nothing_int
      real(kind=8)         , intent(in) :: nothing_dble
      real                 , intent(in) :: nothing_real
      !------------------------------------------------------------------------------------!

      no_time%year       = nothing_int
      no_time%month      = nothing_int
      no_time%day        = nothing_int
      no_time%hour       = nothing_int
      no_time%minu       = nothing_int
      no_time%seco       = nothing_int
      no_time%fracday    = nothing_real
      no_time%doy        = nothing_int
      no_time%elapsed    = nothing_dble
      no_time%mmm        = nothing_char(1:3)
      no_time%yyyyddd    = nothing_int
      no_time%hhmmss     = nothing_int
      no_time%gradsstamp = nothing_char(1:15)
      no_time%timestr    = nothing_char(1:17)

      return
   end function no_time
   !=======================================================================================!
   !=======================================================================================!
end module mod_time

