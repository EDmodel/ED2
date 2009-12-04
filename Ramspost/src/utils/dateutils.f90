!############################# Change Log ##################################
! 1.0.0.4
!
! 010228 MJB date_abs_secs2 ##
!            Added version of date_abs_secs that accepts yyyy mm dd hhmmss
!            as arguments. ##
! 010111 MJB julday1970 ##
!            Added routine to calculate Julian day since 1970. ##
! 001019 CJT date_add_to ##
!            Fix for leap years and negative time increments. ##
! 000828 MJB date_add_to ##
!            Fix to incorrect minute calculation for output filenames. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine date_abs_secs (indate1,seconds)

! compute number of seconds past 1 January 1900 12:00 am

real*8 seconds,s1,s2,s3,s4
character*14 indate1
integer :: year1,month1,date1,hour1

call date_unmake_big(year1,month1,date1,hour1,indate1)

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1= dble(ndays) *86400.
s2= dble(hour1/10000)*3600.
s3= dble(mod(hour1,10000)/100)*60.
s4= dble(mod(hour1,100))
seconds= s1+s2+s3+s4

return
end

!***************************************************************************

subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)

! compute number of seconds past 1 January 1900 12:00 am

real*8 seconds,s1,s2,s3,s4
integer :: year1,month1,date1,hour1

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1= dble(ndays) *86400.
s2= dble(hour1/10000)*3600.
s3= dble(mod(hour1,10000)/100)*60.
s4= dble(mod(hour1,100))
seconds= s1+s2+s3+s4

return
end

!***************************************************************************

subroutine date_subtract (indate1,indate2,tinc,tunits)

! add (or subracts) a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real tinc
character*1 tunits
character*14 indate1, indate2
dimension mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/
integer year1,month1,date1,hour1,year2,month2,date2,hour2
real*8 secs1,secs2

call date_abs_secs(indate1,secs1)
call date_abs_secs(indate2,secs2)
!print*,'sub:',indate1,indate2,secs1,secs2

! convert time to requested unit

ttinc=secs2-secs1
if(tunits.eq.'s') tinc=ttinc
if(tunits.eq.'m') tinc=ttinc/60.
if(tunits.eq.'h') tinc=ttinc/3600.
if(tunits.eq.'d') tinc=ttinc/86400.
!print*,'sub:',secs1,secs2,tinc,tinc

return
end

!***************************************************************************

subroutine date_add_to (inyear,inmonth,indate,inhour  &
                        ,tinc,tunits,outyear,outmonth,outdate,outhour)

! add a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

integer inyear,inmonth,indate,inhour  &
       ,outyear,outmonth,outdate,outhour
real tinc
character*1 tunits
dimension mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/

! convert input time to seconds

ttinc=tinc
if(tunits.eq.'m') ttinc=tinc*60.
if(tunits.eq.'h') ttinc=tinc*3600.
if(tunits.eq.'d') ttinc=tinc*86400.
!print*,'inc:',tinc,tunits,ttinc

xhourin=inhour/10000
xminin=mod(inhour,10000)/100
xsecin=mod(inhour,100)
strtim=xhourin+xminin/60.+xsecin/3600.
!print*,'strtim=',strtim

izhours=int(mod(strtim+ttinc/3600.,24.)+.001)
izmin  =int(mod(strtim+ttinc/3600.,1.)*60+.001)
izsec  =int(mod(strtim*3600.+ttinc,60.)+.001)
!print*,'izs=',izhours,izmin,izsec

outhour= izhours*10000+izmin*100+izsec

iround=.001
if(ttinc<0.) iround=-.001
iadddays=int((strtim+ttinc/3600.)/24.+iround)

outyear=inyear
outdate=indate+iadddays
outmonth=inmonth

20 continue
   idays=mondays(outmonth)
   if(outmonth==2.and.mod(outyear,4)==0)  idays=29

   if(outdate.gt.idays) then
      outdate=outdate-idays
      outmonth=outmonth+1
      if(outmonth.gt.12) then
         outyear=outyear+1
         outmonth=1
      endif
   elseif(outdate.lt.1) then
      if(outmonth.eq.1)outmonth=13
      idays=mondays(outmonth-1)
      if(outmonth-1.eq.2.and.mod(outyear,4).eq.0)  idays=29
      outdate=idays+outdate
      outmonth=outmonth-1
      if(outmonth.eq.12)outyear=outyear-1
   else
      goto 21
   endif

   goto 20

21 continue

!print*,'out stuff:',outyear,outmonth,outdate,outhour

return
end

!***************************************************************************

subroutine date_make_big (inyear,inmonth,indate,inhour,outdate)

integer inyear,inmonth,indate,inhour
character*14 outdate

write(outdate(1:4),10) inyear
write(outdate(5:6),11) inmonth
write(outdate(7:8),11) indate
write(outdate(9:14),12) inhour
10 format (i4.4)
11 format (i2.2)
12 format (i6.6)

return
end

!***************************************************************************

subroutine date_unmake_big (inyear,inmonth,indate,inhour,outdate)

integer inyear,inmonth,indate,inhour
character*14 outdate

read(outdate(1:4),10) inyear
read(outdate(5:6),11) inmonth
read(outdate(7:8),11) indate
read(outdate(9:14),12) inhour
10 format (i4)
11 format (i2)
12 format (i6)

return
end

!***************************************************************************

subroutine RAMS_dintsort(ni,chnums,cstr)

! sort an array of character strings by an associated character field

character*14 chnums(*)
character cstr(*)*(*),cscr*200
character*14 mini,nscr

do n=1,ni
   mini='99999999999999'
   do nm=n,ni
      if(chnums(nm).lt.mini) then
         nmm=nm
         mini=chnums(nm)
      endif
   enddo
   nscr=chnums(n)
   chnums(n)=chnums(nmm)
   chnums(nmm)=nscr
   cscr=cstr(n)
   cstr(n)=cstr(nmm)
   cstr(nmm)=cscr
enddo

return
end

!***************************************************************************

subroutine RAMS_sort_dint3 (n1,ia1,n2,ia2,n3,ia3,nt,iall)

!     sort 3 arrays of char's, put back in 1 array
!     copy all to output array

character*14 ia1(*),ia2(*),ia3(*),iall(*)
character*14 mini,nscr

nt=0
do n=1,n1
   nt=nt+1
   iall(nt)=ia1(n)
enddo
do n=1,n2
   nt=nt+1
   iall(nt)=ia2(n)
enddo
do n=1,n3
   nt=nt+1
   iall(nt)=ia3(n)
enddo

do n=1,nt
   mini='99999999999999'
   do nm=n,nt
      if(iall(nm).lt.mini) then
         nmm=nm
         mini=iall(nm)
      endif
   enddo
   nscr=iall(n)
   iall(n)=iall(nmm)
   iall(nmm)=nscr
enddo

return
end

!***************************************************************************

subroutine RAMS_unique_dint (n1,ia1)

! reduce an array to get rid of duplicate entries

character*14 ia1(*)

nt=n1
10 continue
do n=2,nt
   if(ia1(n).eq.ia1(n-1)) then
      do nn=n,nt
         ia1(nn-1)=ia1(nn)
      enddo
      nt=nt-1
      goto 10
   endif
enddo
n1=nt

return
end

!***************************************************************************

function julday (imonth,iday,iyear)

! compute the julian day from a normal date

julday= iday  &
      + min(1,max(0,imonth-1))*31  &
      + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4))))  &
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
end

!************************************************************************

function julday1970 (imonth,iday,iyear)

! compute the julian day (from 1970) from a normal date w/4 digit yr

julday1970=0
do i=1970,iyear

   imm=12
   idd=31
   if(i==iyear)then
      imm=imonth
      idd=iday
   endif

   jd= idd  &
      + min(1,max(0,imm-1))*31  &
      + min(1,max(0,imm-2))*(28+(1-min(1,mod(i,4))))  &
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

   julday1970=julday1970+jd

enddo    

return
end
