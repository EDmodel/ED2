!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine cuparm()

use mem_tend
use mem_cuparm
use mem_basic
use mem_grid
use node_mod

implicit none

real, save :: cptime=7200.


if (if_cuinv == 0) then

   !        Zero out tendencies initially

   if(time == 0.)then
      call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
      call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
      call azero(mxp*myp,cuparm_g(ngrid)%conprr(1,1))
   endif

   if(initial == 2 .and. time < cptime)return

   if(mod(time+dtlt+.001,dble(confrq)).le.dtlt.or.time.lt..01)then

!!$      print 90,time+dtlt,(time+dtlt)/3600.  &
!!$               +(itimea/100+mod(itimea,100)/60.)
!!$      90   format('  Convective tendencies updated    time =',f10.1,  &
!!$               '  Real time (hrs) =',f6.1)

      call azero(mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
      call azero(mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
      call azero(mxp*myp,cuparm_g(ngrid)%conprr(1,1))

      call conpar(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
          ,basic_g(ngrid)%up      (1,1,1)  &
          ,basic_g(ngrid)%vp      (1,1,1)  &
          ,basic_g(ngrid)%wp      (1,1,1)  &
          ,basic_g(ngrid)%theta   (1,1,1)  &
          ,basic_g(ngrid)%pp      (1,1,1)  &
          ,basic_g(ngrid)%pi0     (1,1,1)  &
          ,basic_g(ngrid)%dn0     (1,1,1)  &
          ,basic_g(ngrid)%rv      (1,1,1)  &
          ,cuparm_g(ngrid)%thsrc  (1,1,1)  &
          ,cuparm_g(ngrid)%rtsrc  (1,1,1)  &
          ,grid_g(ngrid)%rtgt     (1,1)    &
          ,cuparm_g(ngrid)%conprr (1,1), grid_g(ngrid)%flpw(1,1)  )

   endif
   
elseif (if_cuinv == 1) then
   ! Check cumulus inversion tendencies and see if they are usable. If so,
   !   put in thsrc,rtscr,conprr arrays.
   call cu_inv_tend(mzp,mxp,myp,ia,iz,ja,jz  &
                   ,cuparm_g(ngrid)%thsrc(1,1,1)  &
                   ,cuparm_g(ngrid)%thsrcp(1,1,1)  &
                   ,cuparm_g(ngrid)%thsrcf(1,1,1)  &
                   ,cuparm_g(ngrid)%rtsrc(1,1,1)  &
                   ,cuparm_g(ngrid)%rtsrcp(1,1,1)  &
                   ,cuparm_g(ngrid)%rtsrcf(1,1,1)  &
                   ,cuparm_g(ngrid)%conprr (1,1)  &
                   ,cuparm_g(ngrid)%conprrp (1,1)  &
                   ,cuparm_g(ngrid)%conprrf (1,1) )
endif
   

call accum(mxp*myp*mzp,tend%tht(1),cuparm_g(ngrid)%thsrc(1,1,1)    )
call accum(mxp*myp*mzp,tend%rtt(1),cuparm_g(ngrid)%rtsrc(1,1,1)    )

call update(mxp*myp,cuparm_g(ngrid)%aconpr(1,1)  &
                   ,cuparm_g(ngrid)%conprr(1,1)  ,dtlt)

return
end

!*************************************************************************

subroutine cu_inv_tend(m1,m2,m3,ia,iz,ja,jz  &
                      ,thsrc,thsrcp,thsrcf,rtsrc,rtsrcp,rtsrcf  &
                      ,conprr,conprrp,conprrf)

use mem_cuparm
use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m2,m3) :: conprr,conprrp,conprrf
real, dimension(m1,m2,m3) :: thsrc,thsrcp,thsrcf,rtsrc,rtsrcp,rtsrcf

integer :: k,i,j
real :: tfact,grwt

thsrc(1:m1,1:m2,1:m3) = 0.
rtsrc(1:m1,1:m2,1:m3) = 0.
conprr(1:m2,1:m3) = 0.

if (time < tcu_beg .or. time > tcu_end ) return

grwt = wt_cu_grid(ngrid)/tnudcu

if ( (cutime2-cutime1) <= CU_TIL ) then
   ! If past and future files are good for interpolation...
   
   tfact= (time-cutime1)/(cutime2-cutime1) * grwt
   do j=ja,jz
      do i=ia,iz
         do k=2,m1-1
            thsrc(k,i,j) = thsrcp(k,i,j) +  &
                           tfact * (thsrcf(k,i,j)-thsrcp(k,i,j))
            rtsrc(k,i,j) = rtsrcp(k,i,j) +  &
                           tfact * (rtsrcf(k,i,j)-rtsrcp(k,i,j))
         enddo
         conprr(i,j) = conprrp(i,j) +  &
                         tfact * (conprrf(i,j)-conprrp(i,j))
      enddo
   enddo
   
   return

elseif ( abs(time-cutime1) <= CU_TEL ) then
   ! Past file close enough...
   do j=ja,jz
      do i=ia,iz
         do k=2,m1-1
            thsrc(k,i,j) = thsrcp(k,i,j) * grwt
            rtsrc(k,i,j) = rtsrcp(k,i,j) * grwt
         enddo
         conprr(i,j) = conprrp(i,j) * grwt
      enddo
   enddo
   
   return

elseif ( abs(time-cutime2) <= CU_TEL ) then
   ! Future file close enough...
   do j=ja,jz
      do i=ia,iz
         do k=2,m1-1
            thsrc(k,i,j) = thsrcf(k,i,j) * grwt 
            rtsrc(k,i,j) = rtsrcf(k,i,j) * grwt
         enddo
         conprr(i,j) = conprrf(i,j) * grwt
      enddo
   enddo
   
   return
   
endif

return
end

!*************************************************************************

subroutine conpar(m1,m2,m3,ia,iz,ja,jz,ibcon  &
     ,up,vp,wp,theta,pp,pi0,dn0,rv  &
     ,thsrc,rtsrc,rtgt,conprr,flpw)

use conv_coms
use mem_grid
use mem_cuparm
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,lpw
real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3),theta(m1,m2,m3)  &
         ,pp(m1,m2,m3),pi0(m1,m2,m3),dn0(m1,m2,m3),rv(m1,m2,m3)  &
         ,thsrc(m1,m2,m3),rtsrc(m1,m2,m3),conprr(m2,m3)  &
         ,rtgt(m2,m3)
real, dimension(m2,m3) :: flpw

integer :: icpcnt=0,i1,i2,j1,j2,i,j,k,iprtfrq,iqmax,jqmax,kqmax
real :: dthmax

!
!        FLAG TO CONTROL PRINTOUT
!          ICPRTFL=0 - NO PRINTOUT
!                  1 - BRIEF INFO ON UP/DOWN DRAFT AND ITERATIONS
!                  2 - 1 PLUS MODEL TENDENCIES
!                  3 - 2 PLUS FINAL CONVECTIVE STRUCTURE
!                  4 - 3 PLUS UP/DOWN DRAFT AND ENVIRONMENT

icprtfl=0
iprtfrq=8
icpltfl=0
icpcnt=icpcnt+1
if(mod(icpcnt-iprtfrq+1,iprtfrq).eq.0) then
  icprtfl=1
endif

i1 = ia
i2 = iz
j1 = ja
j2 = jz

!  If variable initialization, on a coarse grid, not in a global simulation,
!  do not run convective parameterization in the lateral boundary region.

!  This is commented out for now since it is annoying to compute i1,i2,j1,j2
!    when in parallel. This means we are now computing conv tendencies in the
!    coarse grid nudging boundary regions.
!           We will see if this causes problems ...

!      IF (INITIAL .EQ. 2 .AND. nxtnest(ngrid) .EQ. 0 .and.
!     +    nhemgrd2 .le. 1) THEN

!         if (iand(ibcon,1) .gt. 0) i1 = 1 + nupts
!         if (iand(ibcon,2) .gt. 0) i2 = nxp - nupts
!         if (iand(ibcon,4) .gt. 0) j1 = 1 + nupts * jdim
!         if (iand(ibcon,8) .gt. 0) j2 = max (1,nyp - nupts)

!      ENDIF

dthmax=0.

do j = j1,j2
  do i = i1,i2

    do k = 1,m1
      ucon(k)=up(k,i,j)
      vcon(k)=vp(k,i,j)
      wcon(k)=wp(k,i,j)
      thtcon(k)=theta(k,i,j)
      picon(k)=(pp(k,i,j)+pi0(k,i,j))
      tmpcon(k)=thtcon(k)*picon(k)/cp
      dncon(k)=dn0(k,i,j)
      prcon(k)=(picon(k)/cp)**cpor*p00
      rvcon(k)=rv(k,i,j)
      zcon(k)=zt(k) *rtgt(i,j)
      zzcon(k)=zm(k) *rtgt(i,j)
      thsrc(k,i,j)=0.
      rtsrc(k,i,j)=0.
    enddo
    conprr(i,j)=0.
    wconmin=wcldbs
    contim=confrq

    lpw = nint(flpw(i,j))

    call cu_environ(lpw,m1-1)
    if(igo.ne.0) call kuocp

    if(igo.ne.0) then

      call cp2mod(lpw,m1)
      do k=lpw,m1-1
        thsrc(k,i,j)=ftcon(k)
        rtsrc(k,i,j)=frcon(k)

      enddo
      conprr(i,j)=cprecip

      do k= lpw,m1-1
        if(thsrc(k,i,j).gt.dthmax) then
          dthmax=thsrc(k,i,j)
          iqmax=i
          jqmax=j
          kqmax=k
        endif
      enddo

      if(icprtfl.gt.0) then
!!$        print 899,i,j,time
!!$899         format(' * CONVECTION AT',2I4,'  TIME = ',F10.2)
      endif

    endif

1000     continue

  enddo
enddo

!!$print 675,iqmax,jqmax,kqmax,dthmax*86400.
!!$675 format('   MAX CONVECTIVE HEATING RATE AT I,J,K,D/DAY -',3I4,F8.2)

return
end

!     ******************************************************************

subroutine cu_environ(k1,k2)

use conv_coms
use rconstants

implicit none
integer :: k1,k2

real :: hz(nkp),wcpmax,themax,tlll,plll,rlll,zlll,dzlll,dzdd,abe &
       ,thdu,tdu,rdsu,znz
integer :: k,nkmid,nk

!       Basic constants
dzlow=200.
dzhigh=500.
zmid=3000.

cdzmin=3000.

!         Compute moist static energy profile

do k=k1,k2
  hz(k)=cp*tmpcon(k)+g*zcon(k)+alvl*rvcon(k)
enddo

!         Check for conditional instability and any upward motion
!           greater than WCONMIN under ZMID

igo=0
do k=k1,k2
  if(hz(k).gt.hz(k+1))then
    igo=1
    go to 105
  endif
enddo
105 continue
if(igo.eq.0)return

igo=0
wcpmax=-1.e10
do k=k1,k2
  if(zcon(k).gt.zmid)go to 104
  wcpmax=max(wcpmax,wcon(k))
enddo
104 continue
if(wcpmax.gt.0.0.and.wcpmax.gt.wconmin)igo=1
if(igo.eq.0)return

!           INTERPOLATE MODEL SOUNDING (ENVIRONMENT) TO HIGHER
!             RESOLUTION GRID

nkmid=zmid/dzlow+1
zc(1)=0.
do k=2,nkmid
  zc(k)=zc(k-1)+dzlow
enddo
do k=nkmid+1,nkp
  zc(k)=zc(k-1)+dzhigh
enddo
ze(1)=0.
do k=2,nkp
  ze(k)=(zc(k)+zc(k-1))*.5
enddo
!                   FIND MODEL TOP ON CONVECTIVE GRID
znz=zcon(k2)
do k=nkp,1,-1
  if(ze(k).lt.znz)go to 13
enddo
stop ' envir stop 12'
13 continue
kmt=k
!                   DO ACTUAL INTERPOLATION

nk=k2-k1+1
call htint(nk,ucon,zcon(k1),kmt,upe,ze)
call htint(nk,vcon,zcon(k1),kmt,vpe,ze)
call htint(nk,wcon,zzcon(k1),kmt,wpe,ze)
call htint(nk,thtcon,zcon(k1),kmt,the,ze)
call htint(nk,rvcon,zcon(k1),kmt,rve,ze)
do k=1,kmt
  rve(k)=max(rve(k),1e-8)
enddo

!         COMPUTE THETA V, THETA E, AND GET PRESSURE PROFILE

pke(1)=picon(1)
do k=1,kmt
  thve(k)=the(k)*(1.+.61*rve(k))
enddo
do k=2,kmt
  pke(k)=pke(k-1)-g*2.*(ze(k)-ze(k-1))  &
        /(thve(k)+thve(k-1))
enddo
do k=1,kmt
  te(k)=the(k)*pke(k)/cp
  pe(k)=(pke(k)/cp)**cpor*p00
  rhoe(k)=pe(k)/(rgas*te(k)*(1.+.61*rve(k)))
enddo
do k=1,kmt
  call thetae(pe(k),te(k),rve(k),thee(k))
enddo

!         FIND THE MAIN SOURCE LEVEL OF THE UPDRAFT.

!           FIRST TEST - ANY INVERSION BELOW 1.2 KM

do k=3,nkmid
  if(te(k).gt.te(k-1).and.te(k).gt.te(k+1)  &
                     .and.ze(k).le.1200.)then
  kcon=k
  go to 77
  endif
enddo

!           IF THERE ISN'T AN INVERSION, USE THE LEVEL OF HIGHEST
!           THETA E .

themax=0.
do k=2,nkmid
  if(thee(k).gt.themax)then
    themax=thee(k)
    kcon=k
  endif
enddo

!         FIND THE LCL OF A LAYER AVERAGE AROUND THE SOURCE LEVEL

77 continue
tlll=(te(kcon)+te(kcon+1)+te(kcon-1))/3.
plll=pe(kcon)
rlll=(rve(kcon)+rve(kcon+1)+rve(kcon-1))/3.
zlll=ze(kcon)

call lcl(tlll,plll,rlll,tlcl,plcl,dzlcl)

!         FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL

dzlll=1e20
do k=1,kmt
  dzdd=abs(ze(k)-(zlll+dzlcl))
  if(dzdd.lt.dzlll)then
    dzlll=dzdd
    klcl=k
  endif
enddo

!         IF THERE IS NOT UPWARD MOTION AT THE LCL, NO CONVECTION
!           (MUST BE GREATER THAN WCONMIN )

if(wpe(klcl).lt.0.0.or.wpe(klcl).lt.wconmin)then
  igo=0
  return
endif

!         LOCATE EQUILIBRIUM TEMPERATURE LEVEL OF AN UNENTRAINED PARCEL.
!         COMPUTE INITIAL ABE.  IF ABE IS LESS THAN 0, NO CONVECTION.

theu(klcl)=the(kcon)*exp(alvl*rve(kcon)/(cp*tlcl))

do k=klcl,kmt
  if(theu(klcl).le.thve(k))go to 66
enddo
print*,'convection above model top:',klcl,theu(klcl),thve(k)
ketl=kmt-2
!cccccc      stop 65
66 continue
ketl=k
if(ze(ketl)-ze(klcl).lt.cdzmin)then
  igo=0
  return
endif

abe=0.
do k=klcl,ketl
  call the2t(theu(klcl),pe(k),thdu,tdu,rdsu)
  abe=abe+(thdu*(1.+.61*rdsu)-thve(k))/thve(k)*(zc(k)-zc(k-1))
enddo
if(abe.le.0.)then
  igo=0
  return
endif

!     if(icprtfl.gt.0)then
!       print 899
! 899   format(///,' * convection is activated * ')
!     endif

return
end

!     ******************************************************************

subroutine kuocp

use conv_coms
use rconstants

implicit none

integer :: k,idownd,klfs,kdiv,kdet,kover,kcoolh,kheat
real :: supplyw,anegl,apos,anegh,dddt,dzdiv,wtlfs,wtlcl,wtdiv,wtgnd &
       ,bkuo,zdetr,dzdet,vhint,vmint,vdint,avgmin,avtdiff,overmax,factr &
       ,heatmx,coolhi,c1 

!         Downdraft flag - 0 - no downdrafts
!                          1 - simple downdraft model
idownd=1

do k=1,nkp
  ftcon(k)=0.
  frcon(k)=0.
enddo

!         Compute vertical moisture convergence into the cloud layer.
!           Vertical flux out cloud top is assumed small.

supplyw=rhoe(klcl)*rve(klcl)*(wpe(klcl)+wpe(klcl-1))*.5

supply=supplyw

if(supply.le.0.) then
  igo=0
  return
endif

!         This is the cloud model.  Updraft is constant THETA e and
!           saturated with respect to water.  There is no ice.
!           Cloud top is one level above ETL.
!
!         THETA e of the updraft

theu(klcl)=the(kcon)*exp(alvl*rve(kcon)/(cp*tlcl))

!         Equilibrium Temperature Level of the source level air.

igo=0
do k=klcl,kmt
  call the2t(theu(klcl),pe(k),thu(k),tu(k),rsu(k))
  if(thu(k).gt.the(k).and.igo.eq.0) then
    igo=1
    klfc=k
  endif
  if(thu(k).le.the(k).and.igo.eq.1)go to 66
enddo
if(igo.eq.0) return
PRINT*,' Convection beyond model top - THup, THenv ',THU(KMT)  &
      ,THE(KMT)
k=kmt-1
66 continue
ketl=min(k,kmt)
kct=min(ketl+1,kmt)
call the2t(theu(klcl),pe(kct),thu(kct),tu(kct),rsu(kct))
do k=1,klcl-1
  thu(k)=the(k)
enddo

!         If the cloud is not at least CDZMIN deep or cloud top is
!           under 500 mb, no convection.

if(ze(ketl)-ze(klfc).lt.cdzmin.or.pe(kct).gt.50000.)then
  igo=0
  return
endif

!         Require the positive area be 50% greater than the negative
!           area below the LFC and  5% greater in total.

anegl=0.
do k=klcl,klfc-1
  anegl=anegl+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
apos=0.
do k=klfc,ketl-1
  apos=apos+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
anegh=0.
do k=ketl,kct
  anegh=anegh+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
if(apos.lt.abs(anegl)*1.5.or.apos.lt.abs(anegl+anegh)*1.05) then
  igo=0
  return
endif

if(idownd.eq.1) then

!         The downdraft model - starts at THETA e minimum (LFS).
!             Downdraft is 2 degrees colder than
!             environment at cloud base increasing to 5 degrees
!             colder at the ground.


!         Find LFS as THETA e minimum

do k=kct,2,-1
  if(thee(k).lt.thee(k+1).and.thee(k).lt.thee(k-1))go to 11
enddo
k=2
11 continue
klfs=k
if(klfs.le.klcl)klfs=klcl+1
thd(klfs)=the(klfs)

!        Limit dd deficit at the ground to the maximum of positive
!          temperature difference of updraft if less than 2.5 degrees.

dddt=0.
do k=klcl,kct
  dddt=max(dddt,thu(k)-the(k))
enddo
if(dddt.gt.2.5) dddt=5.

thd(2)=the(2)-dddt
thd(klcl)=the(klcl)-dddt*.2
do k=klcl,klfs
  thd(k)=thd(klcl)+(thd(klfs)-thd(klcl))/(ze(klfs)-ze(klcl))  &
    *(ze(k)-ze(klcl))
enddo
do k=3,klcl-1
  thd(k)=thd(2)+(thd(klcl)-thd(2))/(ze(klcl)-ze(2))  &
    *(ze(k)-ze(2))
enddo

!         Now we need to weight the downdraft relative to the updraft.
!           Assume that the dd weight is zero at the LFS, 1/2 of
!           updraft at cloud base, and equal to the updraft at cloud
!           base at the ground.


dzdiv=1e20
do k=1,kmt
  if(abs(ze(k)-800.).lt.dzdiv)then
    kdiv=k
    dzdiv=abs(ze(k)-800.)
  endif
enddo
kdiv=max(min(klcl,kdiv),2)
if(kdiv.eq.klcl) kdiv=klcl-1

do k=1,nkp
  wtd(k)=0.
enddo
wtlfs=0.
wtlcl=.1
wtdiv=.2
wtgnd=1.
do k=klcl+1,klfs
  wtd(k)=wtlcl+(wtlfs-wtlcl)/(ze(klfs)-ze(klcl))  &
    *(ze(k)-ze(klcl))
enddo
do k=kdiv,klcl
  wtd(k)=wtdiv+(wtlcl-wtdiv)/(ze(klcl)-ze(kdiv))  &
    *(ze(k)-ze(kdiv))
enddo
do k=2,kdiv-1
  wtd(k)=wtgnd+(wtdiv-wtgnd)/(ze(kdiv)-ze(2))  &
    *(ze(k)-ze(2))
enddo

else

  do k=1,nkp
    wtd(k)=0.
  enddo
  do k=2,klcl-1
    thu(k)=the(k)
  enddo

endif

!         Compute infamous b parameter.  Use Fritsch/Chappell's
!           precipitation efficiency.

envshr=sqrt((upe(kct)-upe(klfc))**2  &
           +(vpe(kct)-vpe(klfc))**2)  &
           /(ze(kct)-ze(klfc))*1e3
if(envshr.gt.1.35) then
  preff=1.591-.639*envshr+.0953*envshr**2-.00496*envshr**3
else
  preff=.9
endif
bkuo=1.-preff

!         Vertical profiles of convective heating and moistening

do k=2,kmt
  vheat(k)=0.
  vmois(k)=0.
  vmdry(k)=0.
enddo

!         Find the weighted THETA to use for the convection.

do k=2,kct
  thcon(k)=wtd(k)*thd(k)+(1.-wtd(k))*thu(k)
enddo

!         Heating profile is difference between convective THETAs and
!           environment.

do k=2,kct
  vheat(k)=thcon(k)-the(k)
enddo

!         Moisture profile is difference between vapor's of updraft and
!           environment in the cloud layer.  Below cloud base, air is
!           dried by SUPPLY.  Downdrafts are assumed to have no effect
!           on this.

zdetr=.66667*ze(kct)
dzdet=1000000.
do k=klcl,kct
  if(abs(ze(k)-zdetr).lt.dzdet)then
    dzdet=abs(ze(k)-zdetr)
    kdet=k
  endif
enddo

do k=kdet,kct
  vmois(k)=1.
enddo
!      do k=klcl,kct
!        vmois(k)=rsu(k)-rve(k)
!      enddo

do k=2,klcl-1
  vmdry(k)=rve(k)
enddo

vhint=0.
vmint=0.
vdint=0.
do k=2,kmt
  vhint=vhint+vheat(k)*(zc(k)-zc(k-1))
  vmint=vmint+vmois(k)*(zc(k)-zc(k-1))
  vdint=vdint+vmdry(k)*(zc(k)-zc(k-1))
enddo

!         If VHINT is less than 0, there is more negative area than
!           positive area.  No convection allowed.

if(vhint.le.0.) then
  igo=0
  return
endif

!         Also require that there is a minimum average
!           temperature difference between the updraft and environment
!           from the LFC to the ETL.  This eliminates the cases where
!           VHINT is very small and the heating and cooling rates get
!           astronomically large.

avgmin=.10
avtdiff=0.
do k=klfc,ketl-1
  avtdiff=avtdiff+(thcon(k)-the(k))
enddo
avtdiff=avtdiff/max(1,ketl-klfc)
if(avtdiff.lt.avgmin) then
  igo=0
  return
endif

!         Heating and moistening rates

3100 continue
do k=2,kmt
  ftcon(k)=alvl*preff*supply*vheat(k)  &
     /(pke(k)*rhoe(k)*vhint)
enddo
do k=klcl,kct
  frcon(k)=bkuo*supply*vmois(k)/(rhoe(k)*vmint)
enddo
do k=2,klcl-1
  frcon(k)=-supply*vmdry(k)/(rhoe(k)*vdint)
enddo

do k=klfc,ketl-1
  qvct1(k)=the(k)+contim*ftcon(k)
enddo
overmax=0.
do k=klfc,ketl-1
  if(qvct1(k)-thu(k).gt.overmax)then
    overmax=(qvct1(k)-thu(k))/(ftcon(k)*contim)
    kover=k
  endif
enddo

if(overmax.gt.0.) then
  factr=1.-overmax
  supply=factr*supply
!bob        if(icprtfl.ge.1)print*,' reducing supply ',kover,factr
!bob     +   ,qvct1(kover),thu(kover)
  go to 3100
endif

cprecip=preff*supply

if(icprtfl.gt.0) then
!bob        PRINT*,' ----------------------------------------------------'
!bob        PRINT 898,ZE(KCON),ZE(KLCL),ZE(KLFC),ZE(KETL),ZE(KCT)
!bob  898   FORMAT(' CLOUD LEVELS - SOURCE,LCL,LFC,ETL,TOP(KM) ',-3P,5F5.1)
!bobC       PRINT 896,SUPPLY/SUPPLYW
! 896   FORMAT(' SUPPLIES ' ,F8.4)
!bob        PRINT 897,THEU(KLCL),RSU(KLCL)*1E3
!897   FORMAT(' CLOUD PROPERTIES -  THETA E, RS AT LCL',F6.1,F8.2)
  coolhi=100000.
  heatmx=-10000.
  kcoolh=0
  kheat=0
  do k=ketl,kct
    if(ftcon(k).lt.coolhi) then
      coolhi=ftcon(k)
      kcoolh=k
    endif
  enddo
  do k=klcl,kct
    if(ftcon(k).gt.heatmx) then
      heatmx=ftcon(k)
      kheat=k
    endif
  enddo

  C1=86400.
!bob        PRINT 905,PE(KHEAT),HEATMX*C1,PE(KCOOLH),COOLHI*C1
!905   FORMAT(' MAX-MIN HEATING- P,(K/DAY)', 2(-2PF8.1,0PF7.1))
!bob        PRINT 906,PREFF,PREFF*SUPPLY*3600.*10.
!bob     +           ,(WPE(KLCL)+WPE(KLCL-1))*.5
!bob     +           ,RVE(KLCL)*1E3
!906   FORMAT(' PRECIPITATION - EFFICIENCY,RATE(MM/HR)',F5.2,F6.2,  &
!  '  LCL W,RV',F8.4,F7.2)
ENDIF
!
IF(ICPRTFL.GE.2) THEN
!bob      PRINT95,(K,ZE(K),PE(K),TE(K),THE(K),THEE(K),RVE(K),UCON(K)
!bob     +   ,VCON(K),WPE(K),K=1,KMT)
!95 FORMAT(//' ENVIRONMENT-K,Z,P,TE,THE,THEE,RVE,UP,VP,WP'/,  &
!   (I3, -2P,F8.1,-3P,F8.2,0P,3F7.2,3P,F6.2, -2P,3F8.2))
!bob      PRINT96,(K,THE(K),THU(K)   ,FTCON(K)*86400.,RVE(K),RSU(K),
!bob     +       FRCON(K)*86400.,VHEAT(K)*(ZC(K)-ZC(K-1))/VHINT,
!bob     +                       VMOIS(K)*(ZC(K)-ZC(K-1))/VMINT,
!bob     +                       VMDRY(K)*(ZC(K)-ZC(K-1))/VDINT,
!bob     +  K=2,KMT)
!96 FORMAT(//' HEATING/MOISTENING'  &
! ,'-K,THE,THU,THSRC,RVE,RSU,RTSRC,HEAT%,MOIST%,DRY%'/,  &
!   (I3,0P,3F7.1,3P,3F7.1,0P,3F6.2))
ENDIF
!
!      IF(ICPLTFL.EQ.1)THEN
!      FTCON(1)=FTCON(2)
!      FRCON(1)=FRCON(2)
!      FTMAX=-2000.
!      FRMAX=-2000.
!      FTMIN=2000.
!      FRMIN=2000.
!      DO 400 K=1,KMT
!      QVCT1(K)=FTCON(K)*86400.
!      QVCT2(K)=FRCON(K)*86400.
!      FTMAX=MAX(FTMAX,QVCT1(K))
!      FTMIN=MIN(FTMIN,QVCT1(K))
!      FRMAX=MAX(FRMAX,QVCT2(K))
!      FRMIN=MIN(FRMIN,QVCT2(K))
!      QVCT3(K)=ZE(K)*1E-5
!  400 CONTINUE
!
!        CALL DISPLA(2,1,1)
!        CALL SET (.05,.55,.15,.9,FTMIN-20.,FTMAX+20.
!     +            ,QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'dTH/dt$','height (km)$',1,4,1,'$$$$')
!c       CALL LINE(0.,0.,0.,QVCT3(KMT))
!       CALL EZXY(QVCT1,QVCT3,KMT,'Heating rates (K/day) $')
!C
!        CALL DISPLA(2,1,1)
!        CALL SET (.55,.99,.15,.9,FRMIN-.1,FRMAX+.1
!     +            ,QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'drT/dt$',' $',1,4,0,0)
!        CALL LINE(0.,0.,0.,QVCT3(KMT))
!        CALL EZXY(QVCT2,QVCT3,KMT,'Moistening rates (g/g/day)$')
!        CALL FRAME
!
!C
!      DO 410 K=2,KMT
!      THGRID=THE(K)+FTCON(K)*CONTIM
!      QVCT2(K)=RVE(K)+FRCON(K)*CONTIM
!      QVCT4(K)=RVE(K)
!      TGRID=THGRID*(PE(K)/P00)**ROCP
!      RVVV=RS(PE(K),TGRID)
!      QVCT1(K)=THGRID+AKLV*MAX(0.,QVCT2(K)-RVVV)
!  410 CONTINUE
!      QVCT1(1)=QVCT1(2)
!      QVCT2(1)=QVCT2(2)
!      QVCT4(1)=QVCT4(2)
!        CALL DISPLA(2,1,1)
!        CALL SET (.05,.55,.15,.8,290.,QVCT1(KMT),QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'THETA$','height (km)$',1,4,1,'$$$$')
!        CALL EZXY(THE,QVCT3,KMT,'INITIAL-SOLID  AFTER-DASHED$')
!        CALL ANOTAT(CHAR(0),CHAR(0),4,4,1,'$''$''')
!        CALL EZXY(QVCT1,QVCT3,KMT,CHAR(0))
!
!        CALL DISPLA(2,1,1)
!        CALL SET (.55,.97,.15,.8,  0.,QVCT4(2),QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'MIXING RATIO$',' $',1,4,1,'$$$$')
!        CALL EZXY(QVCT4,QVCT3,KMT,'INITIAL-SOLID  AFTER-DASHED$')
!        CALL ANOTAT(CHAR(0),' $',4,4,1,'$''$''')
!        CALL EZXY(QVCT2,QVCT3,KMT,CHAR(0))
!        CALL FRAME
!
!      ENDIF
!-----------------------------------------------------------------------
RETURN
END

!     ******************************************************************

subroutine cp2mod(k1,k2)

use conv_coms
use mem_scratch
use rconstants

implicit none

integer :: k1,k2

real, external :: ssum

integer :: k
real :: tftc,tftm,tfrc,tfrm,ftres,frres

!        Compute integrated heating and moistening tendencies


do k=2,kmt
  qvct1(k)=rhoe(k)*ftcon(k)*pke(k)
  qvct2(k)=rhoe(k)*alvl*frcon(k)
  qvct3(k)=(zc(k)-zc(k-1))*qvct1(k)
  qvct4(k)=(zc(k)-zc(k-1))*qvct2(k)
enddo
tftc=ssum(kmt-1,qvct3(2),1)
tfrc=ssum(kmt-1,qvct4(2),1)

!         Transfer tendencies to model grid

call vertmap2(qvct1,zc,kmt,vctr5,zzcon,k2)
call vertmap2(qvct2,zc,kmt,vctr6,zzcon,k2)

do k=k1,k2
  vctr5(k)=vctr5(k)*(zzcon(k)-zzcon(k-1))
  vctr6(k)=vctr6(k)*(zzcon(k)-zzcon(k-1))
enddo


!         Make sure the transfer from the convective grid to the model
!           grid happened correctly.

tftm=ssum(k2-k1+1,vctr5(k1),1)
tfrm=ssum(k2-k1+1,vctr6(k1),1)
!
ftres=tftm-tftc
frres=tfrm-tfrc
if(abs(ftres) > .01*abs(tftc)) then
  print*,' energy error in grid tranfser in convective param.'
  print*,' tftm,tftc ',tftm,tftc
endif

!         Change energy tendencies to temperature and mixing ratio
!           tendencies.

do k=k1,k2
  ftcon(k)=vctr5(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*picon(k))
  frcon(k)=vctr6(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*alvl)

enddo

return
end



subroutine vertmap2(datin,zin,n3in,datout,zout,n3out)

implicit none

integer :: n3in,n3out
real, dimension(n3in) :: datin,zin
real, dimension(n3out) :: datout,zout(n3out)

real, allocatable :: qvct(:),vctr(:)
integer :: k,l
real :: dzlft
!  This routine assumes that output vertical structure will be lower than
!   input vertical structure!!!!


!         Transfer quantity from input grid levels to output grid

allocate(qvct(n3in),vctr(n3out))

do k=1,n3out
   vctr(k)=0.
enddo

dzlft=0.
l=2
do k=2,n3out
   !     print *,'******************************* working on output layer ',k
   if(dzlft.ne.0.) then
      if(zin(l) .gt. zout(k)) then
         vctr(k)=vctr(k)+datin(l)*(zout(k)-zout(k-1))
         dzlft=zin(l)-zout(k)
          !  print*,'dzlft2 layer:',k,l,datin(l),dzlft
         go to 61
      else
         vctr(k)=vctr(k)+datin(l)*dzlft
         !   print*,'dzlft layer:',k,l,datin(l),dzlft
         l=l+1
         if(l > n3in) exit
         dzlft=0.
      endif
   endif
60   continue
   if(zin(l) <= zout(k)) then
      vctr(k)=vctr(k)+datin(l)*(zin(l)-zin(l-1))
      !   print*,'full layer:',k,l,datin(l),zin(L),zin(L-1)
      l=l+1
      if(l > n3in) exit
      dzlft=0.
      if(zin(l-1) == zout(k)) go to 61
      go to 60
   else
      vctr(k)=vctr(k)+datin(l)*(zout(k)-zin(l-1))
      !    print*,'part layer:',k,l,vctr(k),datin(l),zout(k),zin(l-1)
      dzlft=zin(l)-zout(k)
   endif
61    continue
      !  print *,'****************************** done with output layer ',k
enddo


!         Change energy tendencies to temperature and mixing ratio
!           tendencies.

do k=2,n3out
  datout(k)=vctr(k)/(zout(k)-zout(k-1))
enddo

deallocate(qvct,vctr)

return
end

!     ******************************************************************

subroutine lcl(t0,pp0,r0,tlcl,plcl,dzlcl)

use rconstants

implicit none
real :: t0,pp0,r0,tlcl,plcl,dzlcl

real, parameter :: cpg=102.45

integer :: nitt,ip
real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs
real, external :: td,rs

ip=0
11 continue

plcl=pp0
tlcl=t0
p0k=pp0**rocp
pi0i=p0k/p00k*cp
ttth0=t0*p00k/p0k
ttd=td(pp0,r0)
dz=cpg*(t0-ttd)
if(dz.le.0.)then
dzlcl=0.
return
endif
do 100 nitt=1,50
pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
pppi=(pki/cp)**cpor*p00
ti=ttth0*pki/cp
rvs=rs(pppi,ti)
if(abs(rvs-r0).lt..00003)go to 110
ttd=td(pppi,r0)
dz=dz+cp/g*(ti-ttd)
!print*,nitt,rvs-r0,ttd,ti,dz
100 continue
print*, 'no converge in LCL:',t0,pp0,r0
ip=ip+1
if(ip==1)go to 11
stop 'LCL no convergence'

110 continue
plcl=pppi
tlcl=ti
dzlcl=dz

return
end
