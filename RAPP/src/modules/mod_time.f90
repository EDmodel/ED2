!==========================================================================================!
!==========================================================================================!
! MEVI. Module time. This module contains the time structure.                              !
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
end module mod_time

