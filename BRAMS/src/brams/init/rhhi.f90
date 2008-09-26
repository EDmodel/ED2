!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine inithh()
! +--------------------------------------------------------+
! _  Initialization of model grid, topography, and fields  _
! _    for horizontally homogeneous fields.                _
! +--------------------------------------------------------+

use mem_basic
use mem_grid

implicit none

integer :: ifm,icm

!     Arrange the input sounding.

call arrsnd

!     For GRID 1, and the second hemispheric grid if applicable,
!     compute the 1-D reference state variables, the 3-D
!     reference state variables, the 3-D model fields, the surface
!     layer parameters, and the initial soil model fields.


do ifm = 1,ngrids
   icm = nxtnest(ifm)
   if (icm .eq. 0) then

      call newgrid(ifm)
      call refs1d

      call refs3d(nzp,nxp,nyp  &
         ,basic_g(ifm)%pi0  (1,1,1)  ,basic_g(ifm)%dn0  (1,1,1)  &
         ,basic_g(ifm)%dn0u (1,1,1)  ,basic_g(ifm)%dn0v (1,1,1)  &
         ,basic_g(ifm)%th0  (1,1,1)  ,grid_g(ifm)%topt  (1,1)    &
         ,grid_g(ifm)%rtgt  (1,1)                                )

      if (if_adap == 0) then
 
         call flds3d(nzp,nxp,nyp  &
            ,basic_g(ifm)%uc  (1,1,1)  ,basic_g(ifm)%vc    (1,1,1)  &
            ,basic_g(ifm)%pi0 (1,1,1)  ,basic_g(ifm)%theta (1,1,1)  &
            ,basic_g(ifm)%thp (1,1,1)  ,basic_g(ifm)%rtp   (1,1,1)  &
            ,basic_g(ifm)%pc  (1,1,1)  ,basic_g(ifm)%rv    (1,1,1)  &
            ,grid_g(ifm)%topt (1,1)    ,grid_g(ifm)%topu   (1,1)    &
            ,grid_g(ifm)%topv (1,1)    ,grid_g(ifm)%rtgt   (1,1)    &
            ,grid_g(ifm)%rtgu (1,1)    ,grid_g(ifm)%rtgv   (1,1)    )

      else

         call flds3d_adap(nzp,nxp,nyp  ,grid_g(ifm)%flpu    (1,1)    &
            ,grid_g(ifm)%flpv  (1,1)    ,grid_g(ifm)%flpw    (1,1)    &
            ,basic_g(ifm)%uc  (1,1,1)  ,basic_g(ifm)%vc    (1,1,1)  &
            ,basic_g(ifm)%pi0 (1,1,1)  ,basic_g(ifm)%theta (1,1,1)  &
            ,basic_g(ifm)%thp (1,1,1)  ,basic_g(ifm)%rtp   (1,1,1)  &
            ,basic_g(ifm)%pc  (1,1,1)  ,basic_g(ifm)%rv    (1,1,1)  )

      endif

   endif
enddo
return
end

!     *****************************************************************

subroutine arrsnd

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use therm_lib, only: rslf,ptrh2rvapl,virtt

implicit none

integer :: nnns,k,kk,kkk
real :: toffset,dir,spd,zold1,zold2,tavg,wt
real, dimension(1) :: rtss
!     Arrange the input sounding

if (ps(1) .eq. 0.) then

   do nsndg=1,maxsndg
      ps(nsndg) = 0.
      ts(nsndg) = 0.
      rts(nsndg) = 0.
      us(nsndg) = 0.
      vs(nsndg) = 0.
   enddo

   open(1,file='SOUND_IN',status='old',form='formatted')
   do nsndg=1,maxsndg
      read(1,*,end=1999) ps(nsndg),ts(nsndg),rts(nsndg),us(nsndg),vs(nsndg)
      if(ps(nsndg).le.0.) go to 1999
   enddo
1999    continue
   close(1)
endif

if(itsflg.eq.0)then
   toffset=t00
elseif(itsflg.eq.1)then
   toffset=0.
endif

do nsndg=1,maxsndg
   nnns=nsndg
   if(ps(nsndg).eq.0.)go to 300
   if(us(nsndg).ne.9999.)then
      if(iusflg.ne.0)then
         dir=us(nsndg)
         spd=vs(nsndg)
         us(nsndg)=-spd*sin(pio180*dir)
         vs(nsndg)=-spd*cos(pio180*dir)
      endif
   endif

   if(ipsflg.eq.0)then
!     pressure given in millibars
      ps(nsndg)=ps(nsndg)*1.e2
   elseif(ipsflg.eq.1)then
!     pressure is height in meters with PS(1)=surface pressure

!  If sounding moisture is expressed as a mixing ratio (IRTSFLG=2),
!     take advantage of knowing virtual temperature effect when
!     integrating hydrostatically to get sounding pressure.
!
      if(irtsflg.eq.2)then
         vctr4(nsndg)=virtt(ts(nsndg+toffset),rts(nsndg)*1.e-3)
      else
         vctr4(nsndg)=ts(nsndg+toffset)
      endif
      if(nsndg.eq.1)then
         ps(nsndg)=ps(nsndg)*100.
         zold2=0.
      else
         zold1=zold2
         zold2=ps(nsndg)
         if(itsflg.eq.0.or.itsflg.eq.1)then
            tavg = 0.5*( vctr4(nsndg) + vctr4(nsndg-1) )
            ps(nsndg)=ps(nsndg-1)*exp(-g*(zold2-zold1)/(rgas*tavg))
         elseif(itsflg.eq.2)then
            tavg=(vctr4(nsndg) + vctr4(nsndg-1)*p00k/ps(nsndg-1)**rocp)*.5
            ps(nsndg)=(ps(nsndg-1)**rocp-g*(zold2-zold1)*p00k/(cp*tavg))**cpor
         endif
      endif
   else
      write(6,131)ipsflg
131        format(' PRESSURE TYPE (IPSFLG=',I2,') NOT IMPLEMENTED')
      stop
   endif

   if(itsflg.eq.0)then
!     Temperature in degrees celsius
      ts(nsndg)=ts(nsndg)+t00
   elseif(itsflg.eq.1)then
!     Temperature in Kelvin
   elseif(itsflg.eq.2)then
!     Temperature is potential temperature in kelvin
      ts(nsndg)=(ps(nsndg)*p00i)**rocp*ts(nsndg)
   else
      write(6,124)itsflg
124        format(' TEMPERATURE TYPE (ITSFLG=',I2,') NOT KNOWN')
      stop
   endif

   if(irtsflg.eq.0)then
!     Humidity given as dew point in degrees celsius
      rts(nsndg) = rslf(ps(nsndg),rts(nsndg)+t00)
   elseif(irtsflg.eq.1)then
!     Humidity given as dew point in degrees kelvin
      rts(nsndg) = rslf(ps(nsndg),rts(nsndg))
   elseif(irtsflg.eq.2)then
!     Humidity given as mixing ratio in g/kg
      rts(nsndg)=rts(nsndg)*1.e-3
   elseif(irtsflg.eq.3)then
!     Humidity given as relative humidity in percent
      rts(nsndg) = ptrh2rvapl(rts(nsndg)*.01,ps(nsndg),ts(nsndg))
   elseif(irtsflg.eq.4)then
!     Humidity given as dew point depression in kelvin
      rts(nsndg) = rslf(ps(nsndg),ts(nsndg)-rts(nsndg))
   else
      write(6,125) irtsflg
125        format(' HUMIDITY TYPE (IRTSFLG=',I2,') NOT KNOWN')
      stop
   endif
enddo
300  continue

nsndg=nnns-1
do k=1,nsndg
   if(us(k).eq.9999.)then
      do kk=k,1,-1
         if(us(kk).ne.9999.)go to 17
      enddo
17         continue
      do kkk=k,nsndg,1
         if(us(kkk).ne.9999.) go to 18
      enddo
18         continue
      wt=(ps(k)-ps(kk))/(ps(kkk)-ps(kk))
      us(k)=us(kk)+wt*(us(kkk)-us(kk))
      vs(k)=vs(kk)+wt*(vs(kkk)-vs(kk))
   endif
enddo

!     compute height levels of input sounding.

do k=2,nsndg
   hs(k)=hs(k-1)-rgas*.5  &
      *(virtt(ts(k),rts(k))+virtt(ts(k-1),rts(k-1)))  &
      *(log(ps(k))-log(ps(k-1)))/g
enddo

if(hs(nsndg).lt.zt(nzp)) then
   write(6,1) hs(nsndg),zt(nzp)
1    format('    Input sounding needs to go higher ! !', /,  &
    '      Sounding top (m) = ',F12.2,'  Model top (m) = ',F12.2)
   stop 'ARRSND'
endif

do k=1,nsndg
   thds(k)=ts(k)*(p00/ps(k))**rocp
enddo

return
end

!     ******************************************************************

subroutine refs1d

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use therm_lib, only: virtt,vapour_on
implicit none
! +---------------------------------------------------------------------
! \   This routine computes the reference state sounding on the model
! \     sigma-z levels from input sounding defined on pressure levels.
! +---------------------------------------------------------------------

integer :: k

if (ztn(nnzp(ngrid),ngrid) .gt. hs(nsndg)) then
   print*,' !!! Input sounding is not high enough !!!'
   print*,' !!! Sounding height: ',hs(nsndg)
   print*,' !!! Model top      : ',ztn(nnzp(ngrid),ngrid)
   stop 'refs1d'
endif

call htint(nsndg,thds,hs,nnzp(ngrid),vctr1,ztn(1,ngrid))
call htint(nsndg,us,hs,nnzp(ngrid),u01dn(1,ngrid),ztn(1,ngrid))
call htint(nsndg,vs,hs,nnzp(ngrid),v01dn(1,ngrid),ztn(1,ngrid))

if (vapour_on) then
   call htint(nsndg,rts,hs,nnzp(ngrid),rt01dn(1,ngrid),ztn(1,ngrid))
else
   do k = 1,nnzp(ngrid)
      rt01dn(k,ngrid) = 0.
   enddo
endif

do k = 1,nnzp(ngrid)
   th01dn(k,ngrid) = virtt(vctr1(k),rt01dn(k,ngrid))
enddo
u01dn(1,ngrid) = u01dn(2,ngrid)
v01dn(1,ngrid) = v01dn(2,ngrid)
rt01dn(1,ngrid) = rt01dn(2,ngrid)
th01dn(1,ngrid) = th01dn(2,ngrid)

pi01dn(1,ngrid) = cp * (ps(1) * p00i) ** rocp  &
   + g * (hs(1) - ztn(1,ngrid))  &
   / (.5 * (th01dn(1,ngrid) + virtt(thds(1),rts(1)) ) )
do k = 2,nnzp(ngrid)
   pi01dn(k,ngrid) = pi01dn(k-1,ngrid) - g / (dzmn(k-1,ngrid)  &
      * .5 * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
enddo

do k = 1,nnzp(ngrid)
   vctr4(k) = (pi01dn(k,ngrid) / cp) ** cpor * p00
   dn01dn(k,ngrid) = cp * vctr4(k)  &
      / (rgas * th01dn(k,ngrid) * pi01dn(k,ngrid))
enddo

return
end

!     ******************************************************************

subroutine flds3d(n1,n2,n3,uc,vc,pi0,theta,thp,rtp,pc,rv  &
   ,topt,topu,topv,rtgt,rtgu,rtgv)

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use therm_lib, only: mrsl,theta_iceliq,tv2temp,vapour_on,cloud_on

implicit none
integer :: n1,n2,n3
real :: pi0(n1,n2,n3),thp(n1,n2,n3),theta(n1,n2,n3)  &
   ,rtp(n1,n2,n3),pc(n1,n2,n3),rv(n1,n2,n3)  &
   ,uc(n1,n2,n3),vc(n1,n2,n3),topt(n2,n3),topu(n2,n3),topv(n2,n3)  &
   ,rtgt(n2,n3),rtgu(n2,n3),rtgv(n2,n3)

integer :: i,j,k
real :: qlatu,qlonu,qlatv,qlonv,dummy
real, dimension(nzpmax) :: p0,temp,rvls,rc

! +---------------------------------------------------------------------
! _    This routine initializes the 3-D velocity and thermodynamic
! _      fields from the 1-D reference state sounding.
! +---------------------------------------------------------------------

do j=1,nyp
   do i=1,nxp

      do k=1,nzp
         vctr11(k) = zt(k) * rtgt(i,j) + topt(i,j)
         vctr12(k) = zt(k) * rtgu(i,j) + topu(i,j)
         vctr13(k) = zt(k) * rtgv(i,j) + topv(i,j)
      enddo

      call htint(nzp, u01dn(1,ngrid),zt,nzp,vctr5,vctr12)
      call htint(nzp, v01dn(1,ngrid),zt,nzp,vctr6,vctr13)
      call htint(nzp,th01dn(1,ngrid),zt,nzp,vctr3,vctr11)
      if(vapour_on) call htint(nzp,rt01dn(1,ngrid),zt,nzp,rtp(1,i,j),vctr11)

! If sounding winds are to be interpreted as eastward (U) and
! northward (V) components, rotate winds from geographic to
! polar stereographic orientation

      if (ihtran .eq. 1) then

         call xy_ll(qlatu,qlonu,platn(ngrid),plonn(ngrid)  &
            ,xm(i),yt(j))
         call xy_ll(qlatv,qlonv,platn(ngrid),plonn(ngrid)  &
            ,xt(i),ym(j))

         do k = 1,nzp
            call uevetouv(uc(k,i,j),dummy,vctr5(k),vctr6(k)  &
                   ,qlatu,qlonu,platn(ngrid),plonn(ngrid))
            call uevetouv(dummy,vc(k,i,j),vctr5(k),vctr6(k)  &
                   ,qlatv,qlonv,platn(ngrid),plonn(ngrid))
         enddo
      else
         do k = 1,nzp
            uc(k,i,j)=vctr5(k)
            vc(k,i,j)=vctr6(k)
         enddo
      endif

      if(vapour_on)then
         do k=1,nzp
            thp(k,i,j)=tv2temp(vctr3(k),rtp(k,i,j))
         enddo
      else
         do k=1,nzp
            thp(k,i,j)=vctr3(k)
         enddo
      endif
      do k=1,nzp
         theta(k,i,j)=thp(k,i,j)
         pc(k,i,j)=0.
      enddo

      if(vapour_on)then
         do k=1,nzp
            rv(k,i,j)=rtp(k,i,j)
         enddo
      endif

      if(cloud_on)then
         do k=1,nzp
            p0(k)=(pi0(k,i,j)/cp)**cpor*p00
            temp(k)=pi0(k,i,j)*thp(k,i,j)/cp
         enddo
         call mrsl(nzp,p0(1),temp(1),rvls(1))
         do k=1,nzp
            rc(k)=max(0.,rtp(k,i,j)-rvls(k))
            thp(k,i,j)=theta_iceliq(pi0(k,i,j),temp(k),rc(k),0.)
            rv(k,i,j)=rtp(k,i,j)-rc(k)
         enddo
      endif

   enddo
enddo

return
end

! ***************************************************************

subroutine flds3d_adap(n1,n2,n3,flpu,flpv,flpw  &
   ,uc,vc,pi0,theta,thp,rtp,pc,rv)

use mem_grid
use ref_sounding
use rconstants
use therm_lib, only: mrsl,tv2temp,vapour_on,level
implicit none
integer :: n1,n2,n3
real, dimension(n2,n3) :: flpu,flpv,flpw
real, dimension(n1,n2,n3) :: uc,vc,pi0,thp,theta,rtp,pc,rv

real, dimension(nzpmax) :: p0,temp,rvls,rc
integer :: k,i,j
real :: qlatu,qlonu,qlatv,qlonv,dummy

! ---------------------------------------------------------------------
! _    This routine initializes the 3-D velocity and thermodynamic 
! _      fields from the 1-D reference state sounding.
! +---------------------------------------------------------------------

do j = 1,n3
   do i = 1,n2

! If sounding winds are to be interpreted as eastward (U) and 
! northward (V) components, rotate winds from geographic to
! polar stereographic orientation

      if (ihtran == 1) then

         call xy_ll(qlatu,qlonu,platn(ngrid),plonn(ngrid),xm(i),yt(j))
         call xy_ll(qlatv,qlonv,platn(ngrid),plonn(ngrid),xt(i),ym(j))


         do k = 1,n1
            call uevetouv(uc(k,i,j),dummy,u01dn(k,ngrid),v01dn(k,ngrid)  &
               ,qlatu,qlonu,platn(ngrid),plonn(ngrid))
         enddo

         do k = 1,n1
            call uevetouv(dummy,vc(k,i,j),u01dn(k,ngrid),v01dn(k,ngrid)  &
               ,qlatv,qlonv,platn(ngrid),plonn(ngrid))
         enddo

      else

         do k = 1,n1
            uc(k,i,j) = u01dn(k,ngrid)
         enddo
         do k = 1,n1
            vc(k,i,j) = v01dn(k,ngrid)
         enddo

      endif

      if (vapour_on) then
         do k = 1,n1
            rtp(k,i,j) = rt01dn(k,ngrid)
            thp(k,i,j) = tv2temp(th01dn(k,ngrid),rtp(k,i,j))
         enddo
      else
         do k = 1,n1
            thp(k,i,j) = th01dn(k,ngrid)
         enddo
      endif

      do k = 1,n1
         theta(k,i,j) = thp(k,i,j)
         pc(k,i,j) = 0.
      enddo

      if (level == 1) then
         do k = 1,n1
            rv(k,i,j) = rtp(k,i,j)
         enddo
      elseif (level >= 2) then
         do k = 1,n1
            p0(k) = (pi0(k,i,j)/cp) ** cpor * p00
            temp(k) = pi0(k,i,j) * thp(k,i,j) / cp
         enddo
         call mrsl(n1,p0(1),temp(1),rvls(1))
         do k = 1,n1
            rc(k) = max(0.,rtp(k,i,j) - rvls(k))
            thp(k,i,j) = theta(k,i,j)  &
               / (1.+(aklv*rc(k)) / max(temp(k),ttripoli))
            rv(k,i,j) = rtp(k,i,j) - rc(k)
         enddo
      endif

   enddo
enddo

return
end





