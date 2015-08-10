!==========================================================================================!
!==========================================================================================!
!      Function that determines the inverse of the number of days of the previous month.   !
!------------------------------------------------------------------------------------------!
subroutine lastmonthdate(thismonth,lastmonth,maxdaysi)
   use ed_misc_coms, only: simtime
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(simtime)               , intent(in)  :: thismonth
   type(simtime)               , intent(out) :: lastmonth
   real                        , intent(out) :: maxdaysi
   !----- Local variables. ----------------------------------------------------------------!
   integer      , dimension(12)              :: maxdays
   !----- External functions. -------------------------------------------------------------!
   logical                     , external    :: isleap
   !---------------------------------------------------------------------------------------!


   !----- Initialise the array assuming it is a regular year. -----------------------------!
   maxdays  = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   !---------------------------------------------------------------------------------------!


   !----- Copy the time structure to the previous month, and assign dummy date and time. --!
   lastmonth      = thismonth
   lastmonth%date = 1
   lastmonth%time = 0.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Shift one month back in time, and check if this would actually be December of the  !
   ! previous year.  If so, make it December of the previous year.                         !
   !---------------------------------------------------------------------------------------!
   lastmonth%month = lastmonth%month -1
   if (lastmonth%month == 0) then
      lastmonth%month = 12
      lastmonth%year  = lastmonth%year - 1
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Change the number of days in February if the previous month falls in a           !
   ! leap year.                                                                            !
   !---------------------------------------------------------------------------------------!
   if (isleap(lastmonth%year)) maxdays(2)=29
   !---------------------------------------------------------------------------------------!

   maxdaysi=1.0/real(maxdays(lastmonth%month))
   return
end subroutine lastmonthdate
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Function that determines the inverse of the number of days of the previous day.     !
!------------------------------------------------------------------------------------------!
subroutine yesterday_info(au_jour_d_hui,hier,maxjoursi)
   use ed_misc_coms, only: simtime
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(simtime)               , intent(in)  :: au_jour_d_hui
   type(simtime)               , intent(out) :: hier
   real                        , intent(out) :: maxjoursi
   !----- Local variables. ----------------------------------------------------------------!
   integer      , dimension(12)              :: maxdays
   !----- External functions. -------------------------------------------------------------!
   logical                     , external    :: isleap
   !---------------------------------------------------------------------------------------!


   !----- Initialise the array assuming it is a regular year. -----------------------------!
   maxdays  = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   !---------------------------------------------------------------------------------------!



   !----- Copy the time structure to the previous month, and assign dummy time. -----------!
   hier      = au_jour_d_hui
   hier%time = 0.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Call function "yesterday" to determine the correct date.                         !
   !---------------------------------------------------------------------------------------!
   call yesterday( au_jour_d_hui%year, au_jour_d_hui%month, au_jour_d_hui%date             &
                 , hier%year         , hier%month         , hier%date          )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Change the number of days in February if the previous month falls in a           !
   ! leap year.                                                                            !
   !---------------------------------------------------------------------------------------!
   if (isleap(hier%year)) maxdays(2)=29
   !---------------------------------------------------------------------------------------!

   maxjoursi = 1.0 / real(maxdays(hier%month))
   return
end subroutine yesterday_info
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Function that determines the which day was the day before.                          !
!------------------------------------------------------------------------------------------!
subroutine yesterday(inyear,inmonth,inday,outyear,outmonth,outday)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)  :: inyear
   integer               , intent(in)  :: inmonth
   integer               , intent(in)  :: inday
   integer               , intent(out) :: outyear
   integer               , intent(out) :: outmonth
   integer               , intent(out) :: outday
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(12)              :: maxdays
   !----- External functions. -------------------------------------------------------------!
   logical               , external    :: isleap
   !---------------------------------------------------------------------------------------!


   !----- Initialise the array assuming it is a regular year. -----------------------------!
   maxdays  = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Shift day by one day.  In case it becomes zero, it means that it should be the    !
   ! last day of the previous month.                                                       !
   !---------------------------------------------------------------------------------------!
   outday = inday - 1
   if (outday == 0) then
      !------------------------------------------------------------------------------------!
      !     Shift month by one month.  In case it becomes zero, it means that it should be !
      ! December of the previous year.                                                     !
      !------------------------------------------------------------------------------------!
      outmonth    = inmonth - 1
      if (outmonth == 0) then
         outyear  = inyear - 1
         outmonth = 12
      else
         outyear  = inyear
      end if

      !----- Fix the number of days in case the previous day falls in a leap year. --------!
      if (isleap(outyear)) maxdays(2)=29
      outday = maxdays(outmonth)
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !    Middle of the month, copy the month and year for current time.                  !
      !------------------------------------------------------------------------------------!
      outmonth = inmonth
      outyear  = inyear
      !------------------------------------------------------------------------------------!
   end if

   return
end subroutine yesterday
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function founds the number of days of a certain month and year.                  !
!------------------------------------------------------------------------------------------!
integer function num_days(month,year)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in) :: month
   integer               , intent(in) :: year
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(12)             :: maxdays
   !----- External functions. -------------------------------------------------------------!
   logical               , external   :: isleap
   !---------------------------------------------------------------------------------------!



   !----- Initialise the array assuming it is a regular year. -----------------------------!
   maxdays  = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   !---------------------------------------------------------------------------------------!



   !----- Correct the number of days in February in case this is a leap year. -------------!
   if (isleap(year)) maxdays(2) = 29
   !---------------------------------------------------------------------------------------!

   !------ Find the number of days of this month. -----------------------------------------!
   num_days = maxdays(month)
   !---------------------------------------------------------------------------------------!

   return
end function num_days
!==========================================================================================!
!==========================================================================================!
