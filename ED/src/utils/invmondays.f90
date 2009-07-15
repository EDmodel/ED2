!------------------------------------------------------------------------------------------!
! Function that determines the inverse of the number of days of the previous month.        !
!------------------------------------------------------------------------------------------!

subroutine lastmonthdate(time,lastmonth,ndaysi)
  use ed_misc_coms, only: simtime
  implicit none
  type(simtime), intent(in)  :: time
  type(simtime), intent(out) :: lastmonth
  real,          intent(out) :: ndaysi
  integer, dimension(12) :: ndays
  logical, external :: isleap

  ndays=(/31,28,31,30,31,30,31,31,30,31,30,31/)

  lastmonth=time
  lastmonth%date = 1
  lastmonth%time = 0.

  lastmonth%month = lastmonth%month -1
  if(lastmonth%month == 0)then
     lastmonth%month = 12
     lastmonth%year  = lastmonth%year - 1
  endif
!----- Changing the number of days in February for leap years -----------------------------!
  if(isleap(lastmonth%year)) ndays(2)=29

  ndaysi=1.0/real(ndays(lastmonth%month))
  return
end subroutine lastmonthdate
!------------------------------------------------------------------------------------------!



!------------------------------------------------------------------------------------------!
! Function that determines the which day was yesterday.                                    !
!------------------------------------------------------------------------------------------!
subroutine yesterday(inyear,inmonth,inday,outyear,outmonth,outday)
  implicit none
  integer, intent(in)  :: inyear,inmonth,inday
  integer, intent(out) :: outyear,outmonth,outday

  integer, dimension(12) :: ndays
  logical, external :: isleap

  ndays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  
  outday=inday - 1
  
  if(outday == 0)then
     outmonth    = inmonth - 1
     
     if(outmonth == 0)then
        outyear  = inyear - 1
        outmonth = 12
     else
        outyear  = inyear
     end if
     
     if(isleap(outyear)) ndays(2)=29
     outday =ndays(outmonth)
  else
     outmonth = inmonth
     outyear  = inyear
  end if
  
  return
end subroutine yesterday
!------------------------------------------------------------------------------------------!



!------------------------------------------------------------------------------------------!
! Function the number of days of a certain month                                           !
!------------------------------------------------------------------------------------------!
integer function num_days(month,year)
   implicit none
   integer, intent(in) :: month, year
   logical, external :: isleap
   integer, dimension(12) :: maxdays
   
   !----- Temporary array with # of days --------------------------------------------------!
   maxdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   if (isleap(year)) maxdays(2) = 29
   
   num_days = maxdays(month)
   return
end function num_days
!------------------------------------------------------------------------------------------!
