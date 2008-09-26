!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

SUBROUTINE input_rawi (olat1,olat2,olon1,olon2)

use isan_coms
use rconstants

implicit none
real :: olat1,olat2,olon1,olon2

character(len=256) :: line,line2
character(len=32) :: cflags,tokens(40)

real, dimension(maxlev) :: p,t,z,h,d,f,pp,zp,tp,hp,vz,uz,zz,rp,dz,fz
integer :: iqflagsp(maxlev,4,3),iqflagsz(maxlev,3,3),iqfl(4)

character(len=4), dimension(maxlev) :: qp,q
character(len=1) :: iq1,iq2,iq3,iq4
character(len=8) :: idsta
character(len=6) :: asta
character(len=14) :: obsdate
real(kind=8) :: nasecs,obssecs

integer :: imarker,iver_up,nobs,k
integer :: ntok,jdate,jt,lp,lz,i,j,jyr,jmo,jdy,kl,iqf,iq
real :: xlat,xlon,elev

call date_abs_secs (iproc_dates(natime),nasecs)
call date_unmake_big (iyear,imonth,idate,ihour,iproc_dates(natime))
ihour=ihour/100

print 559,inrawi
559 format(/'   Acquiring rawindsonde input file    ',A)
call rams_f_open(4,inrawi,'FORMATTED','OLD','READ',0)

read(4,*) imarker
rewind 4

if(imarker.eq.999999) then
   read(4,*) imarker,iver_up
else
   iver_up=1
endif

do nobs=1,10000

   do k=1,maxlev
      p(k)=0.
      t(k)=0.
      z(k)=0.
      h(k)=0.
      f(k)=0.
      pp(k)=0.
      zp(k)=0.
      tp(k)=0.
      hp(k)=0.
      vz(k)=0.
      uz(k)=0.
      zz(k)=0.
      rp(k)=0.
      dz(k)=0.
      fz(k)=0.
   enddo
   
   if(iver_up.eq.1) then
   
      read(4,'(a)',end=100,err=101) line
      call parse(line,tokens,ntok)
      read(tokens(1),*) jdate
      read(tokens(2),*) jt
      read(tokens(3),'(a)') idsta
      read(tokens(4),*) lp
      read(tokens(5),*) lz
      read(tokens(6),*) xlat
      read(tokens(7),*) xlon
      read(tokens(8),*) elev

      IF(lp.GT.maxlev.or. lz.gt.maxlev) then
         print*,'Rawindsonde read error!'
         print*,'  Number of input levels greater than MAXLEV'
         print*,'  lp,lz,maxlev:',lp,lz,maxlev
         stop 'rawin-maxlev'
      endif

      do k=1,lp
         read(4,80)p(k),z(k),t(k),h(k),((iqflagsp(k,i,j),j=1,3),i=1,4)
      enddo
      80 format(2f12.3,f10.2,f10.4,2x,4(':',3i1))

      do k=1,lz
         read(4,85)zz(k),fz(k),dz(k),((iqflagsz(k,i,j),j=1,3),i=1,3)
      enddo
      85 format(f12.3,2f10.2,2x,3(':',3i1))

   elseif(iver_up.eq.2) then
   
      read(4,'(a)',end=100,err=101) line
      call parse(line,tokens,ntok)
      read(tokens(1),*) jyr
      read(tokens(2),*) jmo
      read(tokens(3),*) jdy
      read(tokens(4),*) jt
      read(tokens(5),'(a)') idsta
      read(tokens(6),*) lp
      read(tokens(7),*) lz
      read(tokens(8),*) xlat
      read(tokens(9),*) xlon
      read(tokens(10),*) elev

      IF(lp.GT.maxlev.or. lz.gt.maxlev) then
         print*,'Rawindsonde read error!'
         print*,'  Number of input levels greater than MAXLEV'
         print*,'  lp,lz,maxlev:',lp,lz,maxlev
         stop 'rawin-maxlev'
      endif

      jdate=jyr*10000+jmo*100+jdy

      do k=1,lp
         read(4,*,end=100,err=101) p(k),iqfl(1),z(k),iqfl(2)  &
                                  ,t(k),iqfl(3),h(k),iqfl(4)
         do i=1,4
            iqflagsp(k,i,1)=iqfl(i)/100
            iqflagsp(k,i,2)=mod(iqfl(i)/10,10)
            iqflagsp(k,i,3)=mod(iqfl(i),10)
         enddo
      enddo

      do k=1,lz
          read(4,*,end=100,err=101) zz(k),iqfl(1),fz(k),iqfl(2)  &
                                   ,dz(k),iqfl(3)
          do i=1,3
             iqflagsz(k,i,1)=iqfl(i)/100
             iqflagsz(k,i,2)=mod(iqfl(i)/10,10)
             iqflagsz(k,i,3)=mod(iqfl(i),10)
          enddo
      enddo

   endif
   
   ! Check date
   call date_make_big(jyr,jmo,jdy,jt*100,obsdate)
   call date_abs_secs(obsdate,obssecs)
   !print*,'===> ',natime,iobswin
   !print*,'===> ',iproc_dates(natime),' ',nasecs ,iyear,imonth,idate,ihour
   !print*,'===> ',obsdate            ,' ',obssecs,jyr  ,jmo   ,jdy  ,jt*100
   if((iobswin < 0 .and.  &
      (obssecs .lt. nasecs+iobswin .or. obssecs .gt. nasecs)).or.  &
      (iobswin > 0 .and.  &
      (obssecs .lt. nasecs-iobswin .or. obssecs .gt. nasecs+iobswin)).or.  &
      (iobswin .eq. 0 .and. nasecs .ne. obssecs)) then
      print*,'station ',idsta(1:len_trim(idsta))  &
            ,' excluded - out of obs window (s)',iobswin
      goto 987
   !else
   !   print*,'station ',idsta(1:len_trim(idsta))  &
   !         ,' within obs window (s)',iobswin
   endif

   ! ---------------------------------------------------------
   ! Convert input variables to SI units and appropriate form
   !    Set missing values to 1E30.
   !      PP - Pressure (Pascals)
   !      TP - Temperature (Kelvin)
   !      RP - Relative humidity (fraction)
   !      ZP - Height of pressure levels (meters)
   !      ZZ - Height of wind levels (meters)
   !      DZ - Wind direction (degrees)
   !      FZ - Wind speed (m/s)
   !      XLAT - Station latitude
   !      XLON - Station longitude
   !      ELEV - Station elevation (meters)
   !      LP   - Number of pressure levels
   !      LZ   - Number of wind levels in height
   !      IDSTA- Station ID
   ! ---------------------------------------------------------

   do kl=1,lp
      iqf=1
      do iq=1,3
         if((iqflagsp(kl,1,iq).ne.5.and.iqflagsp(kl,1,iq).ne.0  &
              .and.iqflagsp(kl,1,iq).ne.8).or.p(kl).lt.-998.) iqf=0
      enddo
      !print*,'kl:',kl,p(kl),iqf,iqflagsp(kl,1,1:3)
      if (iqf.eq.1) pp(kl)=p(kl)
      if (iqf.eq.0) pp(kl)=1.e30
      
      iqf=1
      do iq=1,3
         if((iqflagsp(kl,2,iq).ne.5.and.iqflagsp(kl,2,iq).ne.0  &
              .and.iqflagsp(kl,2,iq).ne.8).or.z(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) zp(kl)=z(kl)
      if (iqf.eq.0) zp(kl)=1.e30

      !MJW Make sure that either pressure or height is non-missing
      if(pp(kl) > 1.e20 .and. zp(kl) > 1.e20) then
         tp(kl)=1.e30
         rp(kl)=1.e30
         cycle
      endif

      iqf=1
      do iq=1,3
         if((iqflagsp(kl,3,iq).ne.5.and.iqflagsp(kl,3,iq).ne.0  &
              .and.iqflagsp(kl,3,iq).ne.8).or.t(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) tp(kl)=t(kl)+t00
      if (iqf.eq.0) tp(kl)=1.e30

      iqf=1
      do iq=1,3
         if((iqflagsp(kl,4,iq).ne.5.and.iqflagsp(kl,4,iq).ne.0  &
              .and.iqflagsp(kl,4,iq).ne.8).or.h(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) rp(kl)=h(kl)
      if (iqf.eq.0) rp(kl)=1.e30

   enddo

   do kl=1,lz
      iqf=1
      do iq=1,3
         if((iqflagsz(kl,1,iq).ne.5.and.iqflagsz(kl,1,iq).ne.0  &
              .and.iqflagsz(kl,1,iq).ne.8).or.zz(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) zz(kl)=zz(kl)
      if (iqf.eq.0) zz(kl)=1.e30

      !MJW Make sure that height is non-missing
      if(zz(kl) > 1.e20) then
         fz(kl)=1.e30
         dz(kl)=1.e30
         cycle
      endif

      iqf=1
      do iq=1,3
         if((iqflagsz(kl,2,iq).ne.5.and.iqflagsz(kl,2,iq).ne.0  &
              .and.iqflagsz(kl,2,iq).ne.8).or.fz(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) fz(kl)=fz(kl)
      if (iqf.eq.0) fz(kl)=1.e30

      iqf=1
      do iq=1,3
         if((iqflagsz(kl,3,iq).ne.5.and.iqflagsz(kl,3,iq).ne.0  &
              .and.iqflagsz(kl,3,iq).ne.8).or.dz(kl).lt.-998.) iqf=0
      enddo
      if (iqf.eq.1) dz(kl)=dz(kl)
      if (iqf.eq.0) dz(kl)=1.e30
      
   enddo
   
   !----------------------------------

   call sndproc (lp,lz,xlat,xlon,elev,idsta  &
                ,tp,zp,pp,rp,dz,fz,zz)
                
   987 continue

enddo

100 continue
print*,'End of rawindsonde file:',inrawi
goto 110

101 continue
print*,'Error in rawindsonde file:',inrawi
goto 110

102 continue
print*,'Error opening rawindsonde file:',inrawi

110 continue
close(4)

return
end

!***************************************************************************

subroutine sndproc (lp,lz,xlat,xlon,elev,idsta,tp,zp,pp,rp,dz,fz,zz)

use isan_coms
use rconstants
use therm_lib, only: ptrh2rvapl,virtt
implicit none
                   
real, dimension(*) :: tp,zp,pp,rp,dz,fz,zz
integer :: lp,lz
real :: xlat,xlon,elev
character(len=*) :: idsta

real, dimension(maxlev) :: uz,vz,thz

integer :: not,non,llp,k,llz,kk,lbc
real :: bchyd,rss,pio,tho,zso

do not=1,notsta
   if('r'//idsta(1:5).eq.notid(not)(1:6)) then
      print 31,idsta,xlat,xlon,elev,0,0
      return
   endif
enddo
31 format(' ----- Station omitted   ',1X,A8,' Lat:',F8.2,'  Lon:',F8.2  &
         ,'  Elev:',F10.2,'  P Levels:',I5,'  Wind levels:',I5)

nsta=nsta+1

!--------------------------------------------------------------------
!mjb added the following from 4a as the following fix does not handle
!    pressure levels flagged as suspect by qc (501) 

! Throw out any levels that have suspect pressure
! ALSO, don't read in any levels above the gridded data, since that leads
!  to bulls-eyes around all soundings near model top
      
llp=0
do k=1,lp
   if(pp(k).lt.1.e20) then
      if(pp(k) < levpr(nprz)*100.)then
         print 35,pp(k)/100.,levpr(nprz)
35       format('   discarding ',f6.1, &
                'mb level since it is above gridded data top at ',i4,'mb')
         cycle
      endif
      llp=llp+1
      tp(llp)=tp(k)
      zp(llp)=zp(k)
      pp(llp)=pp(k)
      rp(llp)=rp(k)
   endif
enddo
lp=llp

! Throw out any levels that have suspect winds
llz=0
do k=1,lz
   if(zz(k).lt.1.e20.and.dz(k).lt.1.e20.and.fz(k).lt.1.e20) then
      llz=llz+1
      zz(llz)=zz(K)
      fz(llz)=fz(K)
      dz(llz)=dz(K)
      !print*,'mjb w1',llz,zz(llz),fz(llz),dz(llz)
   endif
enddo
!print*,'mjb - wind levels reduced: ',lz,llz
lz=llz
!--------------------------------------------------------------------

! arrange winds on height levels
do k=1,lz
   if(dz(k).lt.1e19.and.fz(k).lt.1e19) call winduv (dz(k),fz(k),uz(k),vz(k))
enddo
call sort3(zz,uz,vz,lz)

up_lat(nsta)=xlat
up_lon(nsta)=xlon
up_top(nsta)=elev
up_lp(nsta)=lp
up_lz(nsta)=lz
up_chstid(nsta)=idsta(1:8)

! fix missing values if possible
do k=1,lp
   ! missing temperatures are:
   ! set to valid upper temperature if none below,
   ! interpolated linearly in logP if between valid values,
   ! left missing if no valid temperatures are above
   if(tp(k).gt.1e29.and.k.ne.lp) then
      !if(k > 1)then
      !   print*,'missing temp',k,lp,tp(k-1)
      !else
      !   print*,'missing temp',k,lp,tp(k)
      !endif
      kk=k+1
      do while (kk.lt.lp.and.tp(kk).gt.1e29)
         kk=kk+1
      enddo
      if (k.eq.1.and.tp(kk).lt.1e30) then
         tp(k)=tp(kk)
      elseif(tp(kk).lt.1e30) then
         tp(k)=tp(k-1)+(tp(kk)-tp(k-1))*log(pp(k)/pp(k-1))/log(pp(kk)/pp(k-1))
      endif
   endif
   ! same for humidities
   if(rp(k).gt.1e29.and.k.ne.lp) then
      !if(k > 1)then
      !   print*,'missing rh',k,lp,rp(k-1)
      !else
      !   print*,'missing rh',k,lp,rp(k)
      !endif
      kk=k+1
      do while (kk.lt.lp.and.rp(kk).gt.1e29)
         kk=kk+1
      enddo
      if (k.eq.1.and.rp(kk).lt.1e30) then
         rp(k)=rp(kk)
      elseif(rp(kk).lt.1e30) then
         rp(k)=rp(k-1)+(rp(kk)-rp(k-1))*log(pp(k)/pp(k-1))/log(pp(kk)/pp(k-1))
      endif
   endif

enddo

! Now that we filled in as many temps and rh's we could, we will recompute heights.
!        Find a boundary condition as first level at or below 500 mb. If we do 
!        not find a good level, we don't use the sounding.

bchyd=50000.
lbc=0
do k=lp,1,-1
   if(pp(k) < 1.e19 .and. pp(k) >= bchyd .and.  &
      tp(k) < 1.e19 .and. zp(k) < 1.e19) then
      lbc=k
      go to 52
   endif
enddo

! Didn't find a good level, so let's look up from 500mb...
do k=1,lp
   if(pp(k) < 1.e19 .and. pp(k) < bchyd .and.  &
      tp(k) < 1.e19 .and. zp(k) < 1.e19) then
      lbc=k
      go to 52
   endif
enddo
nsta=nsta-1
print*,' sndproc: Could not find good rawindsonde z boundary level.'
print*,'    Discarding sounding:',idsta(1:8)
return

52 continue

! Compute theta - if rh non-missing, compute vitual theta.

do k=1,lp
   thz(k)=1.e30
   if(tp(k) < 1.e19 .and. pp(k) < 1.e19) then
      thz(k)=tp(k)*(p00/pp(k))**rocp
      if(rp(k) < 1.e19) then
         !---------------------------------------------------------------------------------!
         !    Since this is radiosonde, always use liquid water to define vapour mixing    !
         ! ratio. This is the WMO standard, so we assume the data is in standard form.     !
         !---------------------------------------------------------------------------------!
         rss=ptrh2rvapl(rp(k),pp(k),tp(k))
         thz(k)=virtt(thz(k),rss)
      endif
   endif
enddo

! Recompute heights

pio=cp*(pp(lbc)/p00)**rocp
zso=zp(lbc)
tho=thz(lbc)
do k=lbc+1,lp
   zp(k)=1.e30
   if(pp(k) < 1e19 .and. thz(k) < 1e19)then
      zp(k)=zso+.5*(thz(k)+tho)*(pio-cp*(pp(k)/p00)**rocp)/g
      zso=zp(k)
      pio=cp*(pp(k)/p00)**rocp
      tho=thz(k)
   endif
enddo

zso=zp(lbc)
pio=cp*(pp(lbc)/p00)**rocp
tho=thz(lbc)
do k=lbc-1,1,-1
   zp(k)=1.e30
   if(pp(k) < 1e19 .and. thz(k) < 1e19)then
      zp(k)=zso+.5*(thz(k)+tho)*(pio-cp*(pp(k)/p00)**rocp)/g
      zso=zp(k)
      pio=cp*(pp(k)/p00)**rocp
      tho=thz(k)
   endif
enddo

! Do one last filter and fill final arrays. 
!    If any heights are missing, discard the level.

llp=0
do k=1,lp
   if(zp(k) < 1.e19) then
      llp=llp+1
      up_t(nsta,llp)=tp(k)
      up_z(nsta,llp)=zp(k)
      up_p(nsta,llp)=pp(k)
      up_r(nsta,llp)=rp(k)
   endif
enddo

up_lp(nsta)=llp

do k=1,lz
   up_uz(nsta,k)=uz(k)
   up_vz(nsta,k)=vz(k)
   up_zz(nsta,k)=zz(k)
enddo

print 1001,nsta,idsta,up_lat(nsta),up_lon(nsta),up_top(nsta),lp,lz
1001 format(' Rawindsonde station found ',I5,' ID:',1X,A8,' Lat:'  &
     ,F8.2,'  Lon:', F8.2,'  Elev:',0PF10.2,'  P levels:',I5  &
     ,'  Wind levels:',I5)

if(nsta.ge.maxsta) then
   print 92,maxsta
   92 format(' Number of stations found is greater than MAXSTA',I6)
   stop 'st3-maxsta'
endif

return
end

!***************************************************************************

subroutine staprt (nsta,lp,lz,m1,m2,usndz,vsndz,zsndz,psnd,zsnd,rsnd,tsnd)

use rconstants

implicit none
integer :: nsta,lp,lz,m1,m2
real :: usndz(m1,m2),vsndz(m1,m2),psnd(m1,m2),zsnd(m1,m2)  &
         ,rsnd(m1,m2),tsnd(m1,m2),zsndz(m1,m2)

integer :: k

print 1003,(psnd(nsta,k),zsnd(nsta,k)  &
     ,tsnd(nsta,k)*(p00/psnd(nsta,k))**rocp  &
     ,tsnd(nsta,k),rsnd(nsta,k) ,k=1,lp)
1003 format(//'  Pressure sounding'//,T5,'P(mb)',T15,'Height(m)'  &
     ,T25,'Theta(K)'  &
     ,T35,'Temp(K)',T45,'Rel hum(frac)',//  &
     (1x,-2pf10.2,0p3f10.2,0pf10.3))
print 1004,(zsndz(nsta,k),usndz(nsta,k),vsndz(nsta,k),k=1,lz)
1004 format(//'  Height levels',//,T5,'Z(m)',T15,'U(m/s)',  &
     T25,'V(m/s)',//,(1X, 0P3F10.2))
print '(///)'

return
end

!***************************************************************************

SUBROUTINE input_sfc (olat1,olat2,olon1,olon2)

use isan_coms
use rconstants

implicit none

real :: olat1,olat2,olon1,olon2

character*1 itflg,iddd*4,idsta*8,line*256,line2*256,cflags*64
character(len=32) :: tokens(40),var_string
integer :: iqflags(20,3),iqfl(20)
real*8 nasecs

integer :: iver_sfc,nvar,imarker
integer :: ic,nv,ntok,nsss,k,jyr,jmo,jdy,jd,jt,iqf,iq

real :: xlat,xlon,zx,ffx,ddx,tx,tdx,px
   
call date_abs_secs (iproc_dates(natime),nasecs)
call date_unmake_big (iyear,imonth,idate,ihour,iproc_dates(natime))
ihour=ihour/100

! Read surface observations from a file

print 559,insrfce
559 format(/'   Acquiring surface input file    ',A40)
call rams_f_open(4,insrfce,'FORMATTED','OLD','READ',0)

!read(4,*) imarker
!rewind 4

!if(imarker.eq.999999) then
!   read(4,*) imarker,iver_sfc
!else
!   iver_sfc=1
!endif

!read(4,*) nvar
!do nv=1,nvar
!   read(4,*) var_string
!enddo
!if (iver_sfc.eq.1) read(4,*)

do nsss=1,1000000

   read(4,'(a)',end=20,err=21) line
   
   ! Check to see if this is a header line. If so, (re)read header.
   if(line(1:6) == '999999') then
      !read(4,*) imarker
      ! rewind 4

     ! if(imarker.eq.999999) then
     !    read(4,*) imarker,iver_sfc
     ! else
     !    iver_sfc=1
     ! endif

      iver_sfc=2
      read(4,*) nvar
      do nv=1,nvar
         read(4,*) var_string
      enddo
      !if (iver_sfc.eq.1) read(4,*)
      
      read(4,'(a)',end=20,err=21) line
   endif
   
   ! Read a surface obs
   
   if(iver_sfc.eq.1) then
   
!      read(4,'(a)',end=20,err=21) line
      call parse(line,tokens,ntok)
      read(tokens(1),*) jd
      read(tokens(2),*) jt
      read(tokens(3),'(a)') idsta
      read(tokens(4),*) xlat
      read(tokens(5),*) xlon
      read(tokens(6),*) zx
      read(tokens(7),*) ffx
      read(tokens(8),*) ddx
      read(tokens(9),*) tx
      read(tokens(10),*) tdx
      read(tokens(11),*) px
      read(tokens(12),'(a)') cflags

      ic=1
      do nv=1,nvar
         read(cflags(ic:),'(1x,3i1)') (iqflags(nv,k),k=1,3)
         ic=ic+4
      enddo
      
      jyr=jd/10000
      jmo=mod(jd,10000)/100
      jdy=mod(jd,100)

   elseif (iver_sfc.eq.2) then
   
!      read(4,'(a)',end=20,err=21) line
      call parse(line,tokens,ntok)
      read(tokens(1),*) jyr
      read(tokens(2),*) jmo
      read(tokens(3),*) jdy
      read(tokens(4),*) jt
      read(tokens(5),'(a)') idsta
      read(tokens(6),*) xlat
      read(tokens(7),*) xlon
      read(tokens(8),*) zx
      read(tokens(9),*) ffx
      read(tokens(10),*) iqfl(1)
      read(tokens(11),*) ddx
      read(tokens(12),*) iqfl(2)
      read(tokens(13),*) tx
      read(tokens(14),*) iqfl(3)
      read(tokens(15),*) tdx
      read(tokens(16),*) iqfl(4)
      read(tokens(17),*) px
      read(tokens(18),*) iqfl(5)

      do nv=1,nvar
         iqflags(nv,1)=iqfl(nv)/100
         iqflags(nv,2)=mod(iqfl(nv)/10,10)
         iqflags(nv,3)=mod(iqfl(nv),10)
      enddo

   endif
 
   !-------------------------------------------------------------------
   !         Section for conversion to SI  units
   !-------------------------------------------------------------------
   !   At the end of this section, the following variables are
   !   required to be in SI  units.
   !
   !     XLAT   - STATION LATITUDE
   !     XLON   - STATION LONGITUDE
   !       ZX   - STATION ELEVATION
   !      DDX   - WIND DIRECTION
   !      FFX   - WIND SPEED
   !       TX   - TEMPERATURE IN KELVIN
   !      TDX   - DEWPOINT TEMPERATURE IN KELVIN
   !       PX   - SURFACE PRESSURE
   !    IDSTA   - Station ID (char*5)
   !
   !              Any missing values must be set to 1.E30
   !--------------------------------------------------------------------

   if(zx.lt.-998.) goto 987

   iqf=1
   do iq=1,3
      if((iqflags(1,iq).ne.5.and.iqflags(1,iq).ne.0  &
         .and.iqflags(1,iq).ne.8) .or. ffx.lt.-998.) iqf=0
   enddo
   if (iqf.eq.1) ffx=ffx
   if (iqf.eq.0) ffx=1.e30

   iqf=1
   do iq=1,3
      if((iqflags(2,iq).ne.5.and.iqflags(2,iq).ne.0  &
         .and.iqflags(2,iq).ne.8) .or. ddx.lt.-998.) iqf=0
   enddo
   if (iqf.eq.1) ddx=ddx
   if (iqf.eq.0) ddx=1.e30

   iqf=1
   do iq=1,3
      if((iqflags(3,iq).ne.5.and.iqflags(3,iq).ne.0  &
         .and.iqflags(3,iq).ne.8) .or. tx.lt.-998.) iqf=0
   enddo
   if (iqf.eq.1) tx=tx+t00
   if (iqf.eq.0) tx=1.e30

   iqf=1
   do iq=1,3
      if((iqflags(4,iq).ne.5.and.iqflags(4,iq).ne.0  &
         .and.iqflags(4,iq).ne.8) .or. tdx.lt.-998.) iqf=0
   enddo
   if (iqf.eq.1) tdx=tdx+t00
   if (iqf.eq.0) tdx=1.e30

   iqf=1
   do iq=1,3
      if((iqflags(5,iq).ne.5.and.iqflags(5,iq).ne.0  &
         .and.iqflags(5,iq).ne.8) .or. px.lt.-998.) iqf=0
   enddo
   if (iqf.eq.1) px=px
   if (iqf.eq.0) px=1.e30
   
   !-----------------------------------

   call sfcproc (nasecs,jyr,jmo,jdy,jt,idsta,xlat,xlon  &
                ,zx,ddx,ffx,tx,tdx,px)

   987 continue

ENDDO

20 continue
print*,'End of surface file:',insrfce
goto 110

21 continue
print*,'Error reading surface file:', insrfce
goto 110

102 continue
print*,'Error opening surface obs file',insrfce
goto 110

110 continue
close(4)

return
end

!***************************************************************************

subroutine sfcproc (nasecs,jyr,jmo,jdy,jt,idsta,xlat,xlon  &
                   ,zx,ddx,ffx,tx,tdx,px)
                   
use isan_coms
use rconstants
use therm_lib, only : rslf,rehul
implicit none

integer :: jyr,jmo,jdy,jt
real :: xlat,xlon,zx,ddx,ffx,tx,tdx,px
real*8 nasecs
character(len=*) :: idsta

real :: t1000,t900,t850,t700,t600,t500,p1000,p900,p850,p700,p600,p500  &
       ,z1000,z900,z850,z700,z600,z500,spd
real*8 obssecs,obssecsns
character(len=14) :: obsdate
character(len=8) :: csfc

integer :: nsu,ns,not,mis1,mis2,irepl,isfc
integer, save ::ncall=0

call date_make_big(jyr,jmo,jdy,jt*100,obsdate)
call date_abs_secs(obsdate,obssecs)

isfc=0
nssfc=nssfc+1
nsu=nssfc
sf_chstid(nsu)=idsta
sf_top(nsu)=zx
sf_lat(nsu)=xlat
sf_lon(nsu)=xlon
if(ddx.lt.1.e19.and.ffx.lt.1.e19) then
   call winduv(ddx,ffx,sf_u(nsu),sf_v(nsu))
else
   sf_u(nsu)=1.e30
   sf_v(nsu)=1.e30
endif
sf_date(nsu)=obsdate

! check whether to omit station

do not=1,notsta
   if('s'//idsta(1:5).eq.notid(not)(1:6)) then
      isfc=5
      nssfc=nssfc-1
      goto 110
   endif
enddo

! Check date against iobswin

if((iobswin < 0 .and.  &
   (obssecs .lt. nasecs+iobswin .or. obssecs .gt. nasecs)).or.  &
   (iobswin > 0 .and.  &
   (obssecs .lt. nasecs-iobswin .or. obssecs .gt. nasecs+iobswin)).or.  &
   (iobswin .eq. 0 .and. nasecs .ne. obssecs)) then
   isfc=3
   nssfc=nssfc-1
   goto 110
endif

t1000=287.6
t900=281.9
t850=278.6
t700=268.6
t600=260.8
t500=252.0
p1000=100000.
p900=90000.
p850=85000.
p700=70000.
p600=60000.
p500=50000.
z1000=111.
z900=988.
z850=1457.
z700=3012.
z600=4206.
z500=5574.
if(px.gt.1.e20.and.tx.lt.1.e20) then
   if(2.*zx.lt.z1000+z900) then
      px=p1000*exp(-g*(zx-z1000)/(rgas*.5*(tx+t1000)))
   elseif(2.*zx.lt.z900+z850) then
      px=p900*exp(-g*(zx-z900)/(rgas*.5*(tx+t900)))
   elseif(2.*zx.lt.z850+z700) then
      px=p850*exp(-g*(zx-z850)/(rgas*.5*(tx+t850)))
   elseif(2.*zx.lt.z700+z600) then
      px=p700*exp(-g*(zx-z700)/(rgas*.5*(tx+t700)))
   elseif(2.*zx.lt.z600+z500) then
      px=p600*exp(-g*(zx-z600)/(rgas*.5*(tx+t600)))
   endif
endif

if(tx.lt.1.e19.and.px.lt.1.e19) then
   if(tdx.lt.1.e19) then
      !------------------------------------------------------------------------------------!
      !     Assuming relative humidity with respect to the liquid phase. This is observed  !
      ! data, and following the WMO convention this should be in liquid phase.             !
      !------------------------------------------------------------------------------------!
      sf_r(nsu)=rehul(px,tx,rslf(px,tdx))
   else
      sf_r(nsu)=1.e30
   endif
   if(zx.lt.1.e19) then
      sf_s(nsu)=cp*tx+g*zx
   else
      sf_s(nsu)=1.e30
   endif
   sf_t(nsu)=tx*(p00/px)**rocp
   sf_p(nsu)=px
else
   sf_t(nsu)=1.e30
   sf_p(nsu)=1.e30
   sf_r(nsu)=1.e30
   sf_s(nsu)=1.e30
endif

do ns=1,nsu-1

   ! determine if current date is closer for idsta
   
   if(sf_chstid(ns)==idsta) then
      call date_abs_secs(sf_date(ns),obssecsns)
      if(abs(nasecs-obssecsns) > abs(nasecs-obssecs)) then
         isfc=7
         nssfc=nssfc-1
         mis1=0
         if(sf_u(ns).gt.1.e19) mis1=mis1+1
         if(sf_v(ns).gt.1.e19) mis1=mis1+1
         if(sf_t(ns).gt.1.e19) mis1=mis1+1
         if(sf_p(ns).gt.1.e19) mis1=mis1+1
         if(sf_r(ns).gt.1.e19) mis1=mis1+1
         if(sf_top(ns).gt.1.e19) mis1=mis1+1
         mis2=0
         if(sf_u(nsu).gt.1.e19) mis2=mis2+1
         if(sf_v(nsu).gt.1.e19) mis2=mis2+1
         if(sf_t(nsu).gt.1.e19) mis2=mis2+1
         if(sf_p(nsu).gt.1.e19) mis2=mis2+1
         if(sf_r(nsu).gt.1.e19) mis2=mis2+1
         if(sf_top(nsu).gt.1.e19) mis2=mis2+1
         if(mis2.ge.mis1) then
            sf_u(ns)=sf_u(nsu)
            sf_v(ns)=sf_v(nsu)
            sf_t(ns)=sf_t(nsu)
            sf_s(ns)=sf_s(nsu)
            sf_p(ns)=sf_p(nsu)
            sf_top(ns)=sf_top(nsu)
            sf_r(ns)=sf_r(nsu)
            sf_lat(ns)=sf_lat(nsu)
            sf_lon(ns)=sf_lon(nsu)
            csfc=sf_chstid(ns)
            sf_chstid(ns)=sf_chstid(nsu)
            sf_date(ns)=sf_date(nsu)
            isfc=4
            irepl=ns
            goto 110
         endif
      endif
   endif

enddo

do ns=1,nsu-1

   ! check for redundant surface obs

   if(abs(xlon-sf_lon(ns))<stasep.and.abs(xlat-sf_lat(ns))<stasep) then
      isfc=2
      nssfc=nssfc-1
      csfc=sf_chstid(ns)
      irepl=ns
      
      mis1=0
      if(sf_u(ns).gt.1.e19) mis1=mis1+1
      if(sf_v(ns).gt.1.e19) mis1=mis1+1
      if(sf_t(ns).gt.1.e19) mis1=mis1+1
      if(sf_p(ns).gt.1.e19) mis1=mis1+1
      if(sf_r(ns).gt.1.e19) mis1=mis1+1
      if(sf_top(ns).gt.1.e19) mis1=mis1+1
      mis2=0
      if(sf_u(nsu).gt.1.e19) mis2=mis2+1
      if(sf_v(nsu).gt.1.e19) mis2=mis2+1
      if(sf_t(nsu).gt.1.e19) mis2=mis2+1
      if(sf_p(nsu).gt.1.e19) mis2=mis2+1
      if(sf_r(nsu).gt.1.e19) mis2=mis2+1
      if(sf_top(nsu).gt.1.e19) mis2=mis2+1
      if(mis2.lt.mis1) then
         sf_u(ns)=sf_u(nsu)
         sf_v(ns)=sf_v(nsu)
         sf_t(ns)=sf_t(nsu)
         sf_s(ns)=sf_s(nsu)
         sf_p(ns)=sf_p(nsu)
         sf_top(ns)=sf_top(nsu)
         sf_r(ns)=sf_r(nsu)
         sf_lat(ns)=sf_lat(nsu)
         sf_lon(ns)=sf_lon(nsu)
         sf_chstid(ns)=sf_chstid(nsu)
         sf_date(ns)=sf_date(nsu)
         isfc=1
      endif
   endif

enddo

if(nssfc.ge.maxsfc.or.nssfc.gt.maxsname) then
   print 92,maxsfc
   92 format(' Number of sfc obs found greater','than MAXSFC or maxsname ',I6)
   stop 'st3-maxsfc'
endif

110 continue

if(isfc==0) print 100,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu)
100 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2)
   
if(isfc==1) print 101,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl        
101 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' miss - repl ',A8,1X,I5)
     
if(isfc==2) print 102,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl    
102 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' redundant w/',A8,1X,I5)
 
if(isfc==3) print 103,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),iobswin          
103 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' out of obswin',I5)
    
if(isfc==4) print 104,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl       
104 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' newer- repl ',A8,1X,I5)
     
if(isfc==5) print 105,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu)     
105 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' ommitted in namelist')
           
if(isfc==6) print 106,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl   
106 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' redundant w/',A8,1X,I5)
    
if(isfc==7) print 107,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl       
107 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' newer - miss ',A8,1X,I5)

!if(ncall==0) then
!   open(33,file='LANDMARKS-sfc',status='unknown')
!   ncall=1
!endif
!if(isfc==0) then
!   spd=sqrt(sf_v(nsu)**2+sf_u(nsu)**2)
!   if(abs(spd)<20) write(idsta(5:8),'(a1,f3.1)') '+',spd
!   write(33,*) idsta,sf_lat(nsu),sf_lon(nsu)
!endif

return
end
