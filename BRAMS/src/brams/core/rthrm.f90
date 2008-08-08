!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine thermo_boundary_driver(time, dtlong, f_thermo_e, f_thermo_w, f_thermo_s, &
     f_thermo_n, nzp, mxp, myp, jdim)

  implicit none
  ! Arguments
  real, intent(in)    :: dtlong
  ! RGK changed time
  real(kind=8),intent(in) :: time
  ! RGK changed time
  logical, intent(in) :: f_thermo_e, f_thermo_w, f_thermo_s, f_thermo_n
  integer, intent(in) :: nzp, mxp, myp, jdim

  ! Local Variables
  ! real, parameter :: frq_thermo_bd = 100. !in seconds


  ! Checkig frequency of calls
  ! Calling at each Long Timestep
  if (mod(time,dble(dtlong))/=0) return

  

  if (f_thermo_e) call thermo(nzp, mxp, myp, 1,   1,   1, myp, 'THRM_ONLY')

  if (f_thermo_w) call thermo(nzp, mxp, myp, mxp, mxp, 1, myp, 'THRM_ONLY')

  if (jdim==1) then

     if (f_thermo_s) call thermo(nzp, mxp, myp, 1, mxp, 1,   1,   'THRM_ONLY')

     if (f_thermo_n) call thermo(nzp, mxp, myp, 1, mxp, myp, myp, 'THRM_ONLY')

  endif

end subroutine thermo_boundary_driver

!--------------------------------------------------------------------------------

subroutine thermo(mzp,mxp,myp,ia,iz,ja,jz,micflg)

  use mem_grid, only: &
       ngrid !INTENT(IN)
  use mem_basic, only: &
       basic_g !INTENT(INOUT)
  use mem_micro, only: &
       micro_g
  use mem_scratch, only: &
       scratch, &
       vctr1,   &
       vctr2,   &
       vctr3,   &
       vctr4,   &
       vctr5,   &
       vctr6
  use micphys, only: &
       level, & !INTENT(IN)
       jnmb     !INTENT(IN)
  
  implicit none
  ! Arguments:
  integer, intent(in)          :: mzp, mxp, myp, ia, iz, ja, jz
  character(len=*), intent(in) :: micflg !Not used

  if (level .le. 1) then

     call drythrm(mzp,mxp,myp,ia,iz,ja,jz  &
          ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)   &
          ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv(1,1,1),level)

  elseif (level .eq. 2) then

     call satadjst(mzp,mxp,myp,ia,iz,ja,jz  &
          ,basic_g(ngrid)%pp(1,1,1)  ,scratch%scr1(1)             &
          ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1) &
          ,scratch%vt3db(1)          ,basic_g(ngrid)%pi0(1,1,1)   &
          ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv(1,1,1)    &
          ,micro_g(ngrid)%rcp(1,1,1) ,scratch%scr2(1))

  elseif (level .eq. 3) then

     call wetthrm3(mzp,mxp,myp,ia,iz,ja,jz,jnmb  &
          ,basic_g(ngrid)%pi0(1,1,1) ,basic_g(ngrid)%pp   (1,1,1)  &
          ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)  &
          ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv   (1,1,1)  &
          ,micro_g(ngrid)%rcp(1,1,1) ,micro_g(ngrid)%rrp  (1,1,1)  &
          ,micro_g(ngrid)%rpp(1,1,1) ,micro_g(ngrid)%rsp  (1,1,1)  &
          ,micro_g(ngrid)%rap(1,1,1) ,micro_g(ngrid)%rgp  (1,1,1)  &
          ,micro_g(ngrid)%rhp(1,1,1) ,micro_g(ngrid)%q6   (1,1,1)  &
          ,micro_g(ngrid)%q7(1,1,1)  &
          ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,ngrid)

  else

     stop 'Thermo option not supported...LEVEL out of bounds'

  endif

  return
end subroutine thermo
!
!     ***************************************************************
!
subroutine drythrm(m1,m2,m3,ia,iz,ja,jz,thil,theta,rt,rv,level)
  
  ! This routine calculates theta and rv for the case where no condensate is
  ! allowed.

  implicit none

  ! Arguments:
  integer, intent(in)   :: m1,m2,m3,ia,iz,ja,jz,level
  real, intent(in)      :: thil(m1,m2,m3),rt(m1,m2,m3)
  real, intent(inout)   :: theta(m1,m2,m3), rv(m1,m2,m3)
  
  ! Local Variables:
  
  integer               :: i,j,k
  
  do j = ja,jz
     do i = ia,iz
        do k = 1,m1
           theta(k,i,j) = thil(k,i,j)
        enddo
        if (level .eq. 1) then
           do k = 1,m1
              rv(k,i,j) = rt(k,i,j)
           enddo
        endif
     enddo
  enddo
  return
end subroutine drythrm
!
!     ***************************************************************
!
subroutine satadjst(m1,m2,m3,ia,iz,ja,jz  &
     ,pp,p,thil,theta,t,pi0,rtp,rv,rcp,rvls)

  ! This routine diagnoses theta, rv, and rcp using a saturation adjustment
  ! for the case when water is in the liquid phase only
  
  use rconstants, only: &
       cpi,  & ! INTENT(IN)
       p00,  & ! INTENT(IN)
       alvl, & ! INTENT(IN)
       cp      ! INTENT(IN)
  
  implicit none
  
  ! Arguments:
  integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz
  real, intent(in)    :: pp(m1,m2,m3), thil(m1,m2,m3)  &
       ,pi0(m1,m2,m3), rtp(m1,m2,m3)
  real, intent(inout) ::  p(m1,m2,m3), t(m1,m2,m3), rvls(m1,m2,m3), &
       rcp(m1,m2,m3), theta(m1,m2,m3), rv(m1,m2,m3)
  
  ! Local Variables:
  real, external      :: rslf
  integer             :: i,j,k,iterate
  real                :: picpi,til,tt
  integer             :: n
  
  do j = ja,jz
     do i = ia,iz
        do k = 1,m1
           picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
           p(k,i,j) = p00 * picpi ** 3.498
           til = thil(k,i,j) * picpi
           t(k,i,j) = til
           
           do iterate = 1,20
              rvls(k,i,j) = rslf(p(k,i,j),t(k,i,j))
              rcp(k,i,j) = max(rtp(k,i,j) - rvls(k,i,j),0.)
              tt = 0.7 * t(k,i,j) + 0.3 * til  &
                   * (1. + alvl * rcp(k,i,j)  &
                   / (cp * max(t(k,i,j),253.)))
              if (abs(tt - t(k,i,j)) .le. 0.001) go to 1
              t(k,i,j) = tt
           enddo
1          continue
           rv(k,i,j) = rtp(k,i,j) - rcp(k,i,j)
           theta(k,i,j) = t(k,i,j) / picpi
        enddo
     enddo
  enddo
  return
end subroutine satadjst
!
!     ***************************************************************
!
subroutine wetthrm3(m1,m2,m3,ia,iz,ja,jz,jnmb  &
   ,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap,rgp,rhp,q6,q7  &
   ,picpi,tair,til,rliq,rice,qhydm,ngrid)

  ! This routine calculates theta and rv for "level 3 microphysics"
  ! given prognosed theta_il, cloud, rain, pristine ice, snow, aggregates,
  ! graupel, hail, q6, and q7.

  use rconstants, only: &
       cpi,  & ! INTENT(IN)
       alvl, & ! INTENT(IN)
       alvi, & ! INTENT(IN)
       cpi4, & ! INTENT(IN)
       cp253i  ! INTENT(IN)

  implicit none

  ! Arguments:
  integer, intent(in)  :: m1, m2, m3, ia, iz, ja, jz, jnmb(*), ngrid
  real , intent(in)    :: pi0(m1,m2,m3), pp(m1,m2,m3), thp(m1,m2,m3)  &
       ,rtp(m1,m2,m3), rcp(m1,m2,m3), rrp(m1,m2,m3), rpp(m1,m2,m3),   &
       rsp(m1,m2,m3), rap(m1,m2,m3), rgp(m1,m2,m3), rhp(m1,m2,m3),    &
       q6(m1,m2,m3), q7(m1,m2,m3)
  real , intent(inout) :: picpi(*), tair(*), til(*), rliq(*), rice(*), &
       qhydm(*), rv(m1,m2,m3), theta(m1,m2,m3)

  ! Local Variables:
  integer :: i, j, k
  real    :: tcoal, fracliq, tairstr

  do j = ja,jz
     do i = ia,iz

        do k = 1,m1
           picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
           tair(k) = theta(k,i,j) * picpi(k)
           til(k) = thp(k,i,j) * picpi(k)
           rliq(k) = 0.
           rice(k) = 0.
        enddo

        if (jnmb(1) .ge. 1) then
           do k = 1,m1
              rliq(k) = rliq(k) + rcp(k,i,j)
           enddo
        endif

        if (jnmb(2) .ge. 1) then
           do k = 1,m1
              rliq(k) = rliq(k) + rrp(k,i,j)
           enddo
        endif

        if (jnmb(3) .ge. 1) then
           do k = 1,m1
              rice(k) = rice(k) + rpp(k,i,j)
           enddo
        endif

        if (jnmb(4) .ge. 1) then
           do k = 1,m1
              rice(k) = rice(k) + rsp(k,i,j)
           enddo
        endif

        if (jnmb(5) .ge. 1) then
           do k = 1,m1
              rice(k) = rice(k) + rap(k,i,j)
           enddo
        endif

        if (jnmb(6) .ge. 1) then
           do k = 1,m1
              call qtc(q6(k,i,j),tcoal,fracliq)
              rliq(k) = rliq(k) + rgp(k,i,j) * fracliq
              rice(k) = rice(k) + rgp(k,i,j) * (1. - fracliq)
           enddo
        endif

        if (jnmb(7) .ge. 1) then
           do k = 1,m1
              call qtc(q7(k,i,j),tcoal,fracliq)
              rliq(k) = rliq(k) + rhp(k,i,j) * fracliq
              rice(k) = rice(k) + rhp(k,i,j) * (1. - fracliq)
           enddo
        endif

        do k = 1,m1
           qhydm(k) = alvl * rliq(k) + alvi * rice(k)
           rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
        enddo

        do k = 1,m1
           if (tair(k) .gt. 253.) then
              tairstr = 0.5 * (til(k)  &
                   + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
           else
              tairstr = til(k) * (1. + qhydm(k) * cp253i)
           endif
           theta(k,i,j) = tairstr / picpi(k)
        enddo
        
     enddo
  enddo
  return
end subroutine wetthrm3


