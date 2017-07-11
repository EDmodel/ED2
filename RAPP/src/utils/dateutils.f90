!==========================================================================================!
!==========================================================================================!
!   RAPP. dateutils.f90 library. These subroutines were adapted from ED, which were in     !
!         turn adapted from RAMS... They use January 1, 1000, 00 GMT as the time origin    !
!         and expects elapsed time to be in double precision.                              !
!------------------------------------------------------------------------------------------!
subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)
   use rconstants, only : day_sec,hr_sec,min_sec
   implicit none
   real(kind=8) :: seconds

   ! compute number of seconds past 1 January 1900 12:00 am

   real(kind=8) :: s1,s2,s3,s4
   integer :: year1,month1,date1,hour1,iy,ndays
   integer :: elapdays
   integer, external  :: dayofyear
   logical, external  :: isleap
   integer, parameter :: firstyear=1000

   !---------------------------------------------------------------------------------------!
   ! Counting the # of leap days between the reference and current year.                   !
   !---------------------------------------------------------------------------------------!
   elapdays = 0
   if (firstyear < year1) then
      do iy=firstyear,year1-1
         if (isleap(iy)) then
           elapdays=elapdays + 366
         else
           elapdays=elapdays + 365
         end if
      end do
   elseif (firstyear > year1) then
      do iy=firstyear-1,year1,-1
         if (isleap(iy)) then
           elapdays=elapdays - 366
         else
           elapdays=elapdays - 365
         end if
      end do
      elapdays = elapdays - 1
   end if
   
   ndays = elapdays + dayofyear(month1,date1,year1)
   s1= dble(ndays) * day_sec
   s2= dble(hour1/10000) * hr_sec
   s3= dble(mod(hour1,10000)/100)*min_sec
   s4= dble(mod(hour1,100))
   seconds= s1+s2+s3+s4

   return
end subroutine date_abs_secs2
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine date_2_seconds (iyearc,imonthc,idatec,itimec, &
    iyeara,imontha,idatea,itimea,total_seconds)
   implicit none

   ! Output model time = total_seconds
   real(kind=8) :: total_seconds,secs_a,secs_c

   ! Input the model initialization time
   ! NOT ITIMEA IS A SIX DIGIT INTEGER HHMMSS
   integer :: iyeara,imontha,idatea,itimea

   ! Input also the current integer times, and day-seconds
   integer :: iyearc,imonthc,idatec,itimec
   !real :: secondsc
   !integer :: itimec

   !itimec =  10000*int(secondsc/3600) + &     !hours
   !     100*int(secondsc/60) + &              !minutes
   !     mod(int(secondsc),60)                 !seconds

   call date_abs_secs2(iyearc,imonthc,idatec,itimec,secs_c)
   call date_abs_secs2(iyeara,imontha,idatea,itimea,secs_a)
   total_seconds = secs_c - secs_a

return
end subroutine date_2_seconds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine date_secs_ymdt (seconds,iyear1,imonth1,idate1,ihour1)
   use rconstants, only : day_sec,hr_sec,min_sec
   implicit none
   real(kind=8) :: seconds,s1
   integer :: iyear1,imonth1,idate1,ihour1

   integer,parameter :: firstyear=1000
   logical, external :: isleap

   ! compute real time given number of seconds past 1 January 1000 12:00 am  

   integer :: ny,nyr,ileap,nm,ihr,imn,isc

   integer :: mondays(12)
   data mondays/31,28,31,30,31,30,31,31,30,31,30,31/

   ! Get what year it is
   s1=seconds
   do ny=0,10000
      ileap=0
      if(isleap(firstyear+ny)) ileap=1
      s1=s1-(365.+ileap)* day_sec
      if(s1 < 0.) then
         nyr=ny
         s1=s1+(365.+ileap)* day_sec
         exit
      endif
   enddo
   iyear1=firstyear+nyr

   ! s1 is now number of secs into the year
   !   Get month
   do nm=1,12
      ileap=0
      if(isleap(firstyear+ny) .and. nm == 2) ileap=1
      s1=s1-(mondays(nm)+ileap)* day_sec
      if(s1 < 0.) then
         s1=s1+(mondays(nm)+ileap)* day_sec
         exit
      endif
   enddo
   imonth1=nm

   ! s1 is now number of secs into the month
   !   Get date and time

   idate1=int(s1/day_sec)
   s1=s1-idate1*day_sec
   !---- Don't know if this is right. idate1 should not be +1. Although days start
   ! at idate1, idate =0 means that it is still running the last day of previous month,
   ! just that it is past midnight so the previous test gave a negative. 
   ! idate1=idate1+1 ! Since date starts at 1
   if (idate1 == 0) then
     imonth1=imonth1-1
     if (imonth1 == 0) then
       imonth1 = 12
       iyear1  = iyear1 - 1
     end if
     ileap =0
     if (isleap(iyear1) .and. imonth1 == 2) ileap=1
     idate1=mondays(imonth1)+ileap
   end if

   ihr=int(s1/hr_sec)
   s1=s1-ihr*hr_sec
   imn=int(s1/min_sec)
   s1=s1-imn*min_sec
   isc=s1
   ihour1=ihr*10000+imn*100+isc

   return
end subroutine date_secs_ymdt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine date_add_to (inyear,inmonth,indate,inhour  &
                        ,tinc,tunits,outyear,outmonth,outdate,outhour)
   use rconstants, only : day_sec,hr_sec,min_sec
   implicit none

   integer inyear,inmonth,indate,inhour  &
          ,outyear,outmonth,outdate,outhour

   character(len=1) :: tunits

   ! adds/subtracts a time increment to a date and output new date
   ! -> uses hhmmss for hours, 4 digit year


   real(kind=8) :: tinc,ttinc,secs

   ! convert input time to seconds

   ttinc=tinc
   if(tunits.eq.'m') ttinc=tinc*min_sec
   if(tunits.eq.'h') ttinc=tinc*hr_sec
   if(tunits.eq.'d') ttinc=tinc*day_sec
   !print*,'inc:',tinc,tunits,ttinc


   call date_abs_secs2(inyear,inmonth,indate,inhour,secs)

   secs=secs+ttinc

   call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)

   !print*,'out stuff:',outyear,outmonth,outdate,outhour

   return
end subroutine date_add_to
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine date_unmake_big (inyear,inmonth,indate,inhour,outdate)
   implicit none
   integer :: inyear,inmonth,indate,inhour
   character(len=14) :: outdate

   read(outdate(1:4),10) inyear
   read(outdate(5:6),11) inmonth
   read(outdate(7:8),11) indate
   read(outdate(9:14),12) inhour
   10 format (i4)
   11 format (i2)
   12 format (i6)

   return
end subroutine date_unmake_big
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will find the day of year for a given date.                          !
!------------------------------------------------------------------------------------------!
integer function dayofyear(imonth,iday,iyear)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer            , intent(in) :: imonth
   integer            , intent(in) :: iday
   integer            , intent(in) :: iyear
   !----- Local variables. ----------------------------------------------------------------!
   real, dimension(12), parameter  :: elap_regu = (/  0,   31,  59,  90, 120, 151          &
                                                   , 181, 212, 243, 273, 304, 334 /)
   real, dimension(12), parameter  :: elap_leap = (/   0,  31,  60,  91, 121, 152          &
                                                   , 182, 213, 244, 274, 305, 335 /)
   !----- External functions. -------------------------------------------------------------!
   logical            , external   :: isleap
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Compute the day of year using one or other vector, depending on whether the year  !
   ! is leap or not.                                                                       !
   !---------------------------------------------------------------------------------------!

   if (isleap(iyear)) then
      dayofyear = iday + elap_leap(imonth)
   else
      dayofyear = iday + elap_regu(imonth)
   end if

   return
end function dayofyear
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function computes the Julian day (our reference is year 1000, and we assume       !
! Gregorian calendar, although it was implented in 1582).                                  !
!------------------------------------------------------------------------------------------!
integer function julday1000 (imonth,iday,iyear)
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: imonth
   integer, intent(in) :: iday
   integer, intent(in) :: iyear
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: iyy
   !----- External functions. -------------------------------------------------------------!
   logical, external   :: isleap
   integer, external   :: dayofyear
   !---------------------------------------------------------------------------------------!

   julday1000=0
   !---------------------------------------------------------------------------------------!
   !    If the year is at or after our time origin, we cycle over the years before our     !
   ! year, then add the day of year of this year.
   !---------------------------------------------------------------------------------------!
   if (iyear >= 1000) then
      do iyy=1000,iyear-1
         if (isleap(iyy)) then
            julday1000 = julday1000 + 366
         else
            julday1000 = julday1000 + 365
         end if
      end do

      julday1000 = julday1000 + dayofyear(imonth,iday,iyear)

   !---------------------------------------------------------------------------------------!
   !     Otherwise, we subtract year until our year, then add back the day of year of this !
   ! year.                                                                                 !
   !---------------------------------------------------------------------------------------!
   else
      do iyy=999,iyear,-1
         if (isleap(iyy)) then
            julday1000 = julday1000 - 366
         else
            julday1000 = julday1000 - 365
         end if
      end do

      julday1000=julday1000 + dayofyear(imonth,iday,iyear)

   end if 

   return
end function julday1000
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
logical function isleap(year)
   !This function runs a check on whether the year is leap or not, based 
   ! on Gregorian calendar
   integer, intent(in) :: year
   isleap = (mod(year,400) == 0) .or.  &
            (mod(year,4) == 0 .and. mod(year,100) /= 0)

   return
end function isleap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
character(len=3) function monchar(month)
   !----- This function simply gives the month in a 3-character letter --------------------!
   integer, intent(in) :: month
   character(len=3), dimension(12), parameter ::  m3letters=                               &
                 (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
   monchar=m3letters(month)
   return
end function monchar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function day_fraction(hour,minu,seco)
   use rconstants, only: day_sec,hr_sec,min_sec
   !----- This function returns hour/minute/second as the fraction of a full day ----------!
   implicit none
   integer, intent(in) :: hour,minu,seco
   day_fraction = (real(seco)+min_sec*real(minu)+hr_sec*real(hour))/day_sec
end function day_fraction
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function return the number of days in a given month, considering whether the     !
! year is leap or not.                                                                     !
!------------------------------------------------------------------------------------------!
integer function monndays(month,year)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: month
   integer, intent(in) :: year
   !----- Local constants. ----------------------------------------------------------------!
   integer, dimension(12), parameter :: monlength= (/ 31, 28, 31, 30, 31, 30               &
                                                    , 31, 31, 30, 31, 30, 31 /)
   !----- External functions. -------------------------------------------------------------!
   logical, external                 :: isleap
   !---------------------------------------------------------------------------------------!
   
   monndays = monlength(month)

   !----- If it is February in a leap year, add the leap day. -----------------------------!
   if (month == 2 .and. isleap(year)) monndays = monndays + 1

   return
end function monndays
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This function returns the month and day based on day of year (doy) and year.        !
!------------------------------------------------------------------------------------------!
subroutine doy_2_monday(doy,year,month,day,mmm)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer            , intent(in)  :: doy,year
   integer            , intent(out) :: month,day
   character(len=3)   , intent(out) :: mmm
   !----- Local constants. ----------------------------------------------------------------!
   real, dimension(12), parameter   :: elap_regu = (/   0,  31,  59,  90, 120, 151         &
                                                    , 181, 212, 243, 273, 304, 334 /)
   real, dimension(12), parameter   :: elap_leap = (/   0,  31,  60,  91, 121, 152         &
                                                    , 182, 213, 244, 274, 305, 335 /)
   !----- External functions. -------------------------------------------------------------!
   logical            , external    :: isleap
   character(len=3)   , external    :: monchar
   !---------------------------------------------------------------------------------------!

   if (isleap(year)) then
      month=minloc(doy-elap_leap,dim=1,mask=elap_leap < doy)
      day=doy-elap_leap(month)
   else
      month=minloc(doy-elap_regu,dim=1,mask=elap_regu < doy)
      day=doy-elap_regu(month)
   end if
   mmm=monchar(month)

   return
end subroutine doy_2_monday
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine gmt_2_hms(gmt,hour,minu,seco)
   !----- This function returns hour/minute/second as the fraction of a full day ----------!
   implicit none
   real,    intent(in)  :: gmt
   integer, intent(out) :: hour,minu,seco

   hour = int(gmt)
   minu = int(mod(gmt,1.) * 60)
   seco = int(mod(gmt,1.) * 3600) - 60.*minu

end subroutine gmt_2_hms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
integer function v5d_datestamp(year,doy)
   !----- This function simply concatenates the year with the day of year for vis5d -------!
   implicit none
   integer, intent(in) :: year, doy
   v5d_datestamp = 1000 * year + doy
   return
end function v5d_datestamp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
integer function v5d_timestamp(hour,minu,seco)
   !----- This function simply concatenates the hour, minute, and second for vis5d --------!
   implicit none
   integer, intent(in) :: hour,minu,seco

   v5d_timestamp = 10000*hour + 100 *minu + seco
   
   return
end function v5d_timestamp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
character(len=15) function grads_dtstamp(year,mmm,day,hour,minu)
   !----- This function simply concatenates time info for GrADS time stamp ----------------!
   implicit none
   integer         , intent(in) :: year,day,hour,minu
   character(len=3), intent(in) :: mmm
   
   if (minu == 0) then
      write(grads_dtstamp,fmt='(i2.2,a1,i2.2,a3,i4.4)') &
         hour,'z',day,mmm,year
   else
      write(grads_dtstamp,fmt='(2(i2.2,a1),i2.2,a3,i4.4)') &
         hour,':',minu,'z',day,mmm,year
   end if
   
   return
end function grads_dtstamp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
character(len=17) function rapp_dtstamp(year,month,day,hour,minu,seco)
   !----- This function simply concatenates time info for RAPP filename stamp -------------!
   implicit none
   integer, intent(in) :: year,month,day,hour,minu,seco

   write(rapp_dtstamp,fmt='(i4.4,3(a1,i2.2),2i2.2)') &
      year,'-',month,'-',day,'-',hour,minu,seco
   
   return
end function rapp_dtstamp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sort_time(nmytimes,mytimes)
   use mod_time, only : time_stt
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                             , intent(in)   :: nmytimes
   type(time_stt), dimension(nmytimes), intent(inout) :: mytimes
   !----- Local variables -----------------------------------------------------------------!
   type(time_stt)                                     :: timeholder
   real(kind=8)                                       :: elapholder
   real(kind=8) , dimension(nmytimes)                 :: elapsed
   logical      , dimension(nmytimes)                 :: unlocked
   integer                                            :: n,imin
   
   unlocked = .true.
   elapsed  = mytimes(:)%elapsed

   mainloop: do n = 1, nmytimes
      imin          = minloc(elapsed,dim=1,mask=unlocked)
      
      !----- If imin == n, then I don't need to sort anything this time, move on... -------!
      if (imin == n) then
        unlocked(n)   = .false. 
        cycle mainloop
      end if
      timeholder    = mytimes(imin)
      elapholder    = elapsed(imin)

      mytimes(imin) = mytimes(n)
      elapsed(imin) = elapsed(n)

      mytimes(n)    = timeholder
      elapsed(n)    = elapholder

      unlocked(n)   = .false. 
   end do mainloop
   
   return
end subroutine sort_time
!==========================================================================================!
!==========================================================================================!
