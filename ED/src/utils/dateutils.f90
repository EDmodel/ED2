!------------------------------------------------------------------------------------------!
! ED now uses January 1, 1000, 00 GMT as the origin.                                       !
!------------------------------------------------------------------------------------------!

subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)
   use consts_coms, only : day_sec,hr_sec,min_sec
   implicit none
   real(kind=8) :: seconds

   ! compute number of seconds past 1 January 1000 12:00 am

   real(kind=8) :: s1,s2,s3,s4
   integer :: year1,month1,date1,hour1,iy,ndays
   integer :: elapdays
   integer, external  :: julday
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
      do iy=firstyear-1,year1
         if (isleap(iy)) then
           elapdays=elapdays - 366
         else
           elapdays=elapdays - 365
         end if
      end do
   end if
   
   ndays = elapdays + julday(month1,date1,year1)
   s1= dble(ndays) * dble(day_sec)
   s2= dble(hour1/10000) * dble(hr_sec)
   s3= dble(mod(hour1,10000)/100)*dble(min_sec)
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
   use consts_coms, only : day_sec,hr_sec,min_sec
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
      s1=s1-(365.d0 + dble(ileap))* dble(day_sec)
      if(s1 < 0.d0) then
         nyr=ny
         s1=s1+(365.d0 +dble(ileap))* dble(day_sec)
         exit
      endif
   enddo
   iyear1=firstyear+nyr

   ! s1 is now number of secs into the year
   !   Get month
   do nm=1,12
      ileap=0
      if(isleap(firstyear+ny) .and. nm == 2) ileap=1
      s1=s1-dble(mondays(nm)+ileap)* dble(day_sec)
      if(s1 < 0.d0) then
         s1=s1+(dble(mondays(nm)+ileap))* dble(day_sec)
         exit
      endif
   enddo
   imonth1=nm

   ! s1 is now number of secs into the month
   !   Get date and time

   idate1=int(s1/dble(day_sec))
   s1=s1-dble(idate1)*dble(day_sec)
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

   ihr=int(s1/dble(hr_sec))
   s1=s1-dble(ihr)*dble(hr_sec)
   imn=int(s1/dble(min_sec))
   s1=s1-dble(imn)*dble(min_sec)
   isc=int(s1)
   ihour1=ihr*10000+imn*100+isc

   return
end subroutine date_secs_ymdt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine date_add_to (inyear,inmonth,indate,inhour  &
                        ,tinc,tunits,outyear,outmonth,outdate,outhour)
   use consts_coms, only : day_sec,hr_sec,min_sec
   implicit none

   integer inyear,inmonth,indate,inhour  &
          ,outyear,outmonth,outdate,outhour

   character(len=1) :: tunits

   ! adds/subtracts a time increment to a date and output new date
   ! -> uses hhmmss for hours, 4 digit year


   real(kind=8) :: tinc,ttinc,secs

   ! convert input time to seconds

   ttinc=tinc
   if(tunits.eq.'m') ttinc=tinc*dble(min_sec)
   if(tunits.eq.'h') ttinc=tinc*dble(hr_sec)
   if(tunits.eq.'d') ttinc=tinc*dble(day_sec)
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

   read(outdate(1:4),fmt='(i4)') inyear
   read(outdate(5:6),fmt='(i2)') inmonth
   read(outdate(7:8),fmt='(i2)') indate
   read(outdate(9:14),fmt='(i6)') inhour

   return
end subroutine date_unmake_big
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
integer function julday (imonth,iday,iyear)
   implicit none
   integer :: imonth,iday,iyear

   integer           :: febdays
   logical, external :: isleap

   ! compute the julian day from a normal date

   if (isleap(iyear)) then
      febdays=29
   else
      febdays=28
   end if
   
   julday= iday  &
         + min(1,max(0,imonth-1))*31  &
         + min(1,max(0,imonth-2))*febdays  &
         + min(1,max(0,imonth-3))*31  &
         + min(1,max(0,imonth-4))*30  &
         + min(1,max(0,imonth-5))*31  &
         + min(1,max(0,imonth-6))*30  &
         + min(1,max(0,imonth-7))*31  &
         + min(1,max(0,imonth-8))*31  &
         + min(1,max(0,imonth-9))*30  &
         + min(1,max(0,imonth-10))*31  &
         + min(1,max(0,imonth-11))*30  &
         + min(1,max(0,imonth-12))*31

   return
end function julday
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
integer function julday1000 (imonth,iday,iyear)
   implicit none
   integer :: imonth,iday,iyear

   integer :: i,imm,idd,jd
   
   integer           :: febdays
   logical, external :: isleap

   ! compute the julian day (from 1000) from a normal date w/4 digit yr
   ! 1583 is the first full year with Gregorian calendar, but we may want to do historical
   ! runs that predate the calandar so 1000 was chosen instead

   julday1000=0
   do i=1000,iyear

      imm=12
      idd=31
      if(i==iyear)then
         imm=imonth
         idd=iday
      endif
      
      if (isleap(i)) then
         febdays=29
      else
         febdays=28
      end if

      jd= idd  &
         + min(1,max(0,imm-1))*31  &
         + min(1,max(0,imm-2))*febdays  &
         + min(1,max(0,imm-3))*31  &
         + min(1,max(0,imm-4))*30  &
         + min(1,max(0,imm-5))*31  &
         + min(1,max(0,imm-6))*30  &
         + min(1,max(0,imm-7))*31  &
         + min(1,max(0,imm-8))*31  &
         + min(1,max(0,imm-9))*30  &
         + min(1,max(0,imm-10))*31  &
         + min(1,max(0,imm-11))*30  &
         + min(1,max(0,imm-12))*31

      julday1000=julday1000+jd

   enddo    

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
