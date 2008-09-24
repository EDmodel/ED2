!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine oda_proc_obs(m1,m2,m3,i0,j0,pp,pi0,prs,ng,nobs)

  use mem_oda, only:  num_oda_sfc       &
                     ,oda_sfc_info      &
                     ,roda_sfc0         &
                     ,oda_sfc_obs       &
                     ,xkobs             &
                     ,ykobs             &
                     ,zkobs             &
                     ,ekobs             &
                     ,ikobs             &
                     ,jkobs             &
                     ,oda_sfc_til       &
                     ,oda_sfc_tel       &
                     ,ukobs             &
                     ,vkobs             &
                     ,tkobs             &
                     ,rkobs             &
                     ,num_oda_upa       &
                     ,oda_upa_info      &
                     ,roda_upa0         &
                     ,oda_upa_til       &
                     ,oda_upa_tel       &
                     ,oda_upa_obs

  use mem_grid, only: xtn        &
                     ,ytn        &
                     ,time       &
                     ,ztop       &
                     ,ztn        &
                     ,nnzp

  implicit none

  INTEGER, INTENT(IN)     :: m1
  INTEGER, INTENT(IN)     :: m2
  INTEGER, INTENT(IN)     :: m3
  INTEGER, INTENT(IN)     :: i0
  INTEGER, INTENT(IN)     :: j0
  INTEGER, INTENT(IN)     :: ng
  INTEGER, INTENT(OUT)    :: nobs
  REAL,    INTENT(IN)     :: pp(m1,m2,m3)
  REAL,    INTENT(IN)     :: pi0 (m1,m2,m3)
  REAL,    INTENT(INOUT)  :: prs(m2,m3)

  integer, external :: findtime
  integer :: ns,ntimes,ntpast,nobsp,ntfutr,k
  real :: ssx,ssy,oval,rtg

  real, allocatable :: dat1(:),dat2(:)

  allocate(dat1(m1),dat2(m1))


  nobs=0

  print*,'start proc_obs',ng,num_oda_sfc

  ! Do surface obs first

  do ns=1,num_oda_sfc

     ! First filter by horizontal location. If not within RODA_SFC0 
     !     of grid boundaries, we won't use it.

     ssx=oda_sfc_info(ns)%xsta
     ssy=oda_sfc_info(ns)%ysta
     ntimes=oda_sfc_info(ns)%ntimes

     !print*,'sx sy stime:',ssx,ssy,ntimes

     if (ssx < xtn(1+i0,ng) -roda_sfc0(ng) .or. &
          ssx > xtn(m2+i0,ng)+roda_sfc0(ng) .or. &
          ssy < ytn(1+j0,ng) -roda_sfc0(ng) .or. &
          ssy > ytn(m3+j0,ng)+roda_sfc0(ng) ) cycle

     ! We could filter by vertical here, but not now...

     ! Find the obs time for this station before 
     !    the current simulation time

     nobs = nobs + 1
     ntpast=findtime(time, oda_sfc_obs(ns)%time(1), ntimes)
     if(ntpast == 0 ) cycle

     if(nobs == 1) print*,'ntpast:',ntpast,time,ntimes
     xkobs(nobs)=ssx
     ykobs(nobs)=ssy
     zkobs(nobs)=0.
     ekobs(nobs)=oda_sfc_info(ns)%stopo
     ikobs(nobs)=oda_sfc_info(ns)%xista(ng)-i0
     jkobs(nobs)=oda_sfc_info(ns)%xjsta(ng)-j0

     call get_obs_value(ntpast,time, oda_sfc_obs(ns)%time(1), ntimes  &
          ,oda_sfc_obs(ns)%u(1), oda_sfc_til, oda_sfc_tel, oval)

     ukobs(nobs)=oval

     call get_obs_value(ntpast, time, oda_sfc_obs(ns)%time(1), ntimes  &
          ,oda_sfc_obs(ns)%v(1), oda_sfc_til, oda_sfc_tel, oval)

     vkobs(nobs)=oval

     call get_obs_value(ntpast, time, oda_sfc_obs(ns)%time(1), ntimes  &
          ,oda_sfc_obs(ns)%temp(1), oda_sfc_til, oda_sfc_tel, oval)  !(OUT:oval)

     tkobs(nobs)=oval

     !if(nobs==461) print*,'oval:',nobs,oval,oda_sfc_obs(ns)%temp(1:ntimes)
     !if(nobs==461) print*,'oval:',nobs,oval,oda_sfc_obs(ns)%time(1:ntimes)
     !if(nobs == 1) print*,'oval:',oval,oda_sfc_obs(ns)%time(1:ntimes)
     !if(nobs == 1) print*,'oval:',oval,oda_sfc_obs(ns)%temp(1:ntimes)

     call get_obs_value(ntpast, time,oda_sfc_obs(ns)%time(1), ntimes  &
          ,oda_sfc_obs(ns)%dewpt(1), oda_sfc_til, oda_sfc_tel, oval)

     rkobs(nobs)=oval

     ! call get_obs_value(ntpast,time,oda_sfc_obs(ns)%time(1),ntimes  &
     !       ,oda_sfc_obs(ns)%ps(1),oda_sfc_til,oda_sfc_tel,oval)
     ! pkobs(nobs)=oval

  enddo

  ! We have the set of values. Convert temp, dewpt to theta, mixrat.
  !    Most of the time, we need to interpolate sfc pressure from model fields.

  call sfc_obs_convert(m1,m2,m3,pp,pi0,prs,ng,nobs)


  ! Now do upper air

  print*,'start proc upa', ng,num_oda_upa

  do ns=1,num_oda_upa

     ! First filter by horizontal location. If not within RODA_UPA0 
     !     of grid boundaries, we won't use it.

     ssx=oda_upa_info(ns)%xsta
     ssy=oda_upa_info(ns)%ysta
     ntimes=oda_upa_info(ns)%ntimes
     print*,ns,ssx,ssy,ntimes
     if (ssx < xtn(1+i0,ng)  -roda_upa0(ng) .or. &
          ssx > xtn(m2+i0,ng) +roda_upa0(ng) .or. &
          ssy < ytn(1+j0,ng)  -roda_upa0(ng) .or. &
          ssy > ytn(m3+j0,ng) +roda_upa0(ng) ) cycle


     ! Find the obs time for this station before 
     !    the current simulation time

     nobsp = nobs + 1
     nobs  = nobs + m1
     ntpast=findtime(time, oda_upa_obs(ns)%time(1), ntimes)
     print*,'afterntp',ntpast,time,oda_upa_obs(ns)%time(1:ntimes)
     if(ntpast == 0 ) cycle

     xkobs(nobsp:nobs)=ssx
     ykobs(nobsp:nobs)=ssy
     ikobs(nobsp:nobs)=oda_upa_info(ns)%xista(ng)-i0
     jkobs(nobsp:nobs)=oda_upa_info(ns)%xjsta(ng)-j0

     ekobs(nobsp:nobs)=oda_upa_info(ns)%stopo

     ! Compute sigma-z heights based on sounding topo
     rtg=(1.-oda_upa_info(ns)%stopo/ztop)
     do k=1,m1
        zkobs(k+nobsp-1)=oda_upa_info(ns)%stopo+ztn(k,ng)*rtg
     enddo



     !  We could use upper air obs points directly in the krig, but that 
     !     wouldn't let us easily interpolate in time. So we'll first interpolate 
     !     to RAMS grid levels.

     !  We'll do this:
     !    - Interpolate ntpast and ntpast+1 soundings to RAMS vertical grid.
     !         We won't go beyond these times for now.
     !    - Assume that if a sounding is so incomplete as to be useless for
     !         this interpolation, it will have been removed on initialization.
     !    - Fill missing values above and below data
     !    - Handle each variable and level in the same manner as surface obs

     !  We have to fill the *kobs arrays assuming multiple variables at same 
     !    locations. So interpolate all variables first, then fill *kobs
     !    arrays correctly

     ntfutr=ntpast+1
     if(ntfutr > ntimes ) ntfutr=ntpast
     print*,'after',ntpast,ntfutr,ntimes

     ! Processing theta.....

     call upa_interp(oda_upa_obs(ns)%lp(ntpast)  &
          ,oda_upa_obs(ns)%theta(1,ntpast)  &
          ,oda_upa_obs(ns)%zgeo(1,ntpast)  &
          ,nnzp(ng), dat1(1), zkobs(nobsp) )

     call upa_interp(oda_upa_obs(ns)%lp(ntfutr)  &
          ,oda_upa_obs(ns)%theta(1,ntfutr)  &
          ,oda_upa_obs(ns)%zgeo(1,ntfutr)  &
          ,nnzp(ng), dat2(1), zkobs(nobsp) )

     call get_upaobs_value(time,oda_upa_obs(ns)%time(ntpast)  &
          ,oda_upa_obs(ns)%time(ntfutr)   &
          ,nnzp(ng), dat1(1), dat2(1),tkobs(nobsp)  &
          ,oda_upa_til, oda_upa_tel )

     ! Processing rv.....

     call upa_interp(oda_upa_obs(ns)%lp(ntpast)  &
          ,oda_upa_obs(ns)%rv(1,ntpast)  &
          ,oda_upa_obs(ns)%zgeo(1,ntpast)  &
          ,nnzp(ng), dat1(1), zkobs(nobsp) )

     call upa_interp(oda_upa_obs(ns)%lp(ntfutr)  &
          ,oda_upa_obs(ns)%rv(1,ntfutr)  &
          ,oda_upa_obs(ns)%zgeo(1,ntfutr)  &
          ,nnzp(ng), dat2(1), zkobs(nobsp) )

     call get_upaobs_value(time,oda_upa_obs(ns)%time(ntpast)  &
          ,oda_upa_obs(ns)%time(ntfutr)   &
          ,nnzp(ng), dat1(1), dat2(1),rkobs(nobsp)  &
          ,oda_upa_til, oda_upa_tel )

     ! Processing u.....

     call upa_interp(oda_upa_obs(ns)%lz(ntpast)  &
          ,oda_upa_obs(ns)%u(1,ntpast)  &
          ,oda_upa_obs(ns)%zz(1,ntpast)  &
          ,nnzp(ng), dat1(1), zkobs(nobsp) )

     call upa_interp(oda_upa_obs(ns)%lz(ntfutr)  &
          ,oda_upa_obs(ns)%u(1,ntfutr)  &
          ,oda_upa_obs(ns)%zz(1,ntfutr)  &
          ,nnzp(ng), dat2(1), zkobs(nobsp) )

     call get_upaobs_value(time,oda_upa_obs(ns)%time(ntpast)  &
          ,oda_upa_obs(ns)%time(ntfutr)   &
          ,nnzp(ng), dat1(1), dat2(1),ukobs(nobsp)  &
          ,oda_upa_til, oda_upa_tel )

     ! Processing v.....

     call upa_interp(oda_upa_obs(ns)%lz(ntpast)  &
          ,oda_upa_obs(ns)%v(1,ntpast)  &
          ,oda_upa_obs(ns)%zz(1,ntpast)  &
          ,nnzp(ng), dat1(1), zkobs(nobsp) )

     call upa_interp(oda_upa_obs(ns)%lz(ntfutr)  &
          ,oda_upa_obs(ns)%v(1,ntfutr)  &
          ,oda_upa_obs(ns)%zz(1,ntfutr)  &
          ,nnzp(ng), dat2(1), zkobs(nobsp) )

     call get_upaobs_value(time,oda_upa_obs(ns)%time(ntpast)  &
          ,oda_upa_obs(ns)%time(ntfutr)   &
          ,nnzp(ng), dat1(1), dat2(1),vkobs(nobsp)  &
          ,oda_upa_til, oda_upa_tel )


  enddo

  deallocate(dat1,dat2)

  !do k=1,nnzp(ng)
  !   print '(i3,5f12.4)',k,ukobs(k),vkobs(k),tkobs(k),rkobs(k),zkobs(k)
  !enddo

  return
end subroutine oda_proc_obs
!==============================

subroutine get_upaobs_value(time,timep,timef,nz,dat1,dat2,oval &
     ,oda_til,oda_tel)

  implicit none

  INTEGER, INTENT(IN)  :: nz
  REAL(KIND=8),    INTENT(IN)  :: time
  REAL(KIND=8),    INTENT(IN)  :: timep
  REAL(KIND=8),    INTENT(IN)  :: timef
  REAL,    INTENT(IN)  :: oda_til
  REAL,    INTENT(IN)  :: oda_tel
  REAL,    INTENT(OUT) :: oval(nz)
  REAL,    INTENT(IN)  :: dat1(nz)
  REAL,    INTENT(IN)  :: dat2(nz)

  integer :: k

  oval(1:nz)=-999.

  do k=1,nz

     ! - Time interpolate limit (ODA_TIL)- if the future-past obs
     !    is > this limit, do not use to interpolate
     !
     ! - Time extrapolate limit (ODA_TEL)- if past/future obs is greater than TIL,
     !    but less than TEL, use obs

     if (dat1(k) > -998. .and. dat2(k) > -998.) then
        ! We have both past and future good values. Check TIL:

        if (timef-timep <= dble(oda_til)) then
           ! If we got here, we can interpolate...

           oval(k)=dat1(k) + sngl((time-timep) *  &
                dble(dat2(k)-dat1(k))/(timef-timep))
           cycle
        endif
     endif

     ! Failed the TIL check, check TEL now
     if(dat1(k) > -998. .and. abs(time-timep) < dble(oda_tel)) then
        oval(k)=dat1(k)
     elseif(dat2(k) > -998. .and. abs(time-timef) < dble(oda_tel)) then
        oval(k)=dat2(k)
     endif

  enddo

  return
end subroutine get_upaobs_value


!==============================

subroutine upa_interp(np,varp,zp,nz,varz,zz)

  ! Interpolate vertical sounding data to RAMS sigma-z based on
  !    topography of the sounding

  implicit none

  INTEGER, INTENT(IN)  :: np
  INTEGER, INTENT(IN)  :: nz

  REAL,    INTENT(IN)  :: varp(np)
  REAL,    INTENT(IN)  :: zp(np)
  REAL,    INTENT(OUT) :: varz(nz) 
  REAL,    INTENT(IN)  :: zz(nz)

  integer :: np1,np2,k,kk,nb,nt
  real, allocatable :: v1(:)

  allocate(v1(max(np,nz)))


  varz(1:nz)=-999.

  ! First find actual top and bottom of good data. Return if problems.
  np1=0
  do k=1,np
     if (varp(k) > -998.) then
        np1=k
        exit
     endif
  enddo

  np2=0
  do k=np,1,-1
     if (varp(k) > -998.) then
        np2=k
        exit
     endif
  enddo

  if (np1 == 0 .or. np2 == 0 .or. np2-np1 <= 0) return

  ! Now fill in any missing values between np1 and np2. We're assuming
  !   that a single sounding does not have huge gaps and that all 
  !   zp (geopotential heights) are filled.

  nb=np1
  nt=nb+1
  v1(np1)=varp(np1)
  do k=np1+1,np2
     if(varp(k) > -998.) then
        v1(k)=varp(k)
        nb=k
        cycle
     else
        do kk=k,np2
           if(varp(k) > -998.) then
              nt=kk
              exit
           endif
        enddo
        v1(k) = v1(nb)+(varp(nt)-varp(nb))*(zp(k)-zp(nb))/(zp(nt)-zp(nb))
        nb=k
     endif
  enddo


  ! Finally, interpolate to sig-z levels

  kk=np1
  do k=1,nz
30   continue
     if(zz(k) < zp(1)) then
        cycle
     elseif(zz(k) > zp(np2)) then
        cycle
     elseif(zz(k) >= zp(kk).and.zz(k) <= zp(kk+1)) then
        varz(k)=varp(kk)+(varp(kk+1)-varp(kk))*(zz(k)-zp(kk))/(zp(kk+1)-zp(kk))
        cycle
     endif
     kk=kk+1
     if(kk == np2) then
        print *,'upa_interp:np',np
        do kk=1,np2
           print*,'kk,zp(kk),zz(kk)',zp(kk),zz(kk)
        enddo
        stop 'upa_interp'
     endif
     goto 30
  enddo

  deallocate(v1)

  return
end subroutine upa_interp

!==============================

subroutine sfc_obs_convert(n1,n2,n3,pp,pi0,prs,ng,nobs)

  use mem_oda
  use rconstants
  use therm_lib, only : rslif

  implicit none

  INTEGER, INTENT(IN) :: n1
  INTEGER, INTENT(IN) :: n2
  INTEGER, INTENT(IN) :: n3
  INTEGER, INTENT(IN) :: ng
  INTEGER, INTENT(IN) :: nobs
  
  REAL, INTENT(IN)    :: pp(n1,n2,n3)
  REAL, INTENT(IN)    :: pi0(n1,n2,n3)
  REAL, INTENT(INOUT) :: prs(n2,n3)

  integer :: i,j,ns


  ! Routine to convert input sfc temp and dewpt to theta and mixing ratio


  ! Compute surface pressure

  do j=1,n3
     do i=1,n2
        prs(i,j)=( (pp(1,i,j)+pp(2,i,j)+pi0(1,i,j)+pi0(2,i,j))*.5  &
             *cpi) ** cpor * p00
     enddo
  enddo

  ! Interpolate to station locations

!!!!!!! what to do if obs not on this grid or SUBDOMAIN ?????????
!!!!!!! use station pressure if reported ?????????

  do ns=1,nobs
     call gdtost(prs,n2,n3,ikobs(ns),jkobs(ns),pkobs(ns))
     if (pkobs(ns) > 1e10) pkobs(ns)=-999.
  enddo


  do ns=1,nobs
     if (rkobs(ns) > -998. .and. pkobs(ns) > -998.) then
        rkobs(ns) = rslif(pkobs(ns),rkobs(ns)+t00)
     else
        rkobs(ns)=-999.
     endif
     if (tkobs(ns) > -998. .and. pkobs(ns) > -998.) then
        !if(ns == 258) print*,'tconv:',ns,tkobs(ns),pkobs(ns)
        tkobs(ns) = (tkobs(ns)+t00)*(p00/pkobs(ns))**rocp
        !if(tkobs(ns) > 200. .and. tkobs(ns) < 273.) print*,'tconv:',ns,tkobs(ns),pkobs(ns)
     else
        tkobs(ns)=-999.
     endif
  enddo



  return
end subroutine sfc_obs_convert

!==============================

subroutine get_obs_value(ntpast,rtime,otimes,ntimes,odata  &
     ,oda_til,oda_tel,oval)

  implicit none

  INTEGER,         INTENT(IN)  :: ntpast
  INTEGER,         INTENT(IN)  :: ntimes
  REAL(KIND=8),    INTENT(IN)  :: rtime
  REAL,            INTENT(IN)  :: oda_til
  REAL,            INTENT(IN)  :: oda_tel
  REAL,            INTENT(OUT) :: oval

  REAL(KIND=8),    INTENT(IN)  :: otimes(ntimes)
  REAL,            INTENT(IN)  :: odata(ntimes)


  integer, external :: findgood

  integer :: nu1,nu2


  oval=-999.

  if (ntpast == 0) then
     if (abs(rtime-otimes(1)) < dble(oda_tel)) then
        oval=odata(1)
     endif
     return
  endif

  ! Find future time where there is a non-missing value,
  !   then go backward and find earlier non-missing value.

  nu2=findgood(ntpast+1,ntimes,odata(1),ntimes,1)
  nu1=findgood(1,nu2-1,odata(1),ntimes,2)
  ! - Time interpolate limit (ODA_TIL)- if the future-past obs
  !    is > this limit, do not use to interpolate
  !
  ! - Time extrapolate limit (ODA_TEL)- if past/future obs is greater than TIL,
  !    but less than TEL, use obs

  if (nu1 > 0 .and. nu2 > 0) then
     ! We have both past and future good values. Check TIL:

     if (otimes(nu2)-otimes(nu1) <= dble(oda_til)) then
        ! If we got here, we can interpolate...

        oval=odata(nu1) + sngl((rtime-otimes(nu1)) *  &
             dble(odata(nu2)-odata(nu1))/(otimes(nu2)-otimes(nu1)))
        return
     endif
  endif

  ! Failed the TIL check, check TEL now
  if(nu1 > 0 .and. abs(rtime-otimes(nu1)) < dble(oda_tel)) then
     oval=odata(nu1)
  elseif(nu2 > 0 .and. abs(rtime-otimes(nu2)) < dble(oda_tel)) then
     oval=odata(nu2)
  endif

  return
end subroutine get_obs_value

!==============================

integer function findgood(n1,n2,a,na,idir)
  implicit none

  INTEGER, INTENT(IN) :: n1
  INTEGER, INTENT(IN) :: n2
  INTEGER, INTENT(IN) :: idir
  INTEGER, INTENT(IN) :: na
  REAL,    INTENT(IN) :: a(na)

  integer :: nt, ntfound

  ! Assumes atime is sorted lowest to highest
  if(idir == 1) then
     do nt=n1,n2
        if (a(nt) > -998.) then
           findgood=nt
           return
        endif
     enddo
  elseif (idir == 2) then
     do nt=n2,n1,-1
        if (a(nt) > -998.) then
           findgood=nt
           return
        endif
     enddo
  endif

  findgood=0
  return
end function findgood

!==============================

integer function findtime(time,atime,ntime)
  implicit none

  INTEGER,         INTENT(IN) :: ntime
  REAL(KIND=8),            INTENT(IN) :: atime(ntime)
  REAL(KIND=8),    INTENT(IN) :: time
  
  integer :: nt,ntfound

  ! Assumes atime is sorted lowest to highest
  if(time < atime(1) ) then
     ntfound=0
     go to 100
  endif

  do nt=ntime,1,-1
     if (time >= atime(nt)) then
        ntfound=nt
        go to 100
     endif
  enddo
  ntfound=ntime

100 continue

  findtime=ntfound

  return
end function findtime

