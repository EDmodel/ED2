!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!  FOR ICNFLG=3, DEBATING WHETHER TO KEEP IT #/M4 OR CHANGE
!  PARM TO #/KG/M.  NEED TO DEFINE AVMIPSA, ETC. FOR ALL CATEGORIES.
!  MAY WANT TO DEFINE C1 TOO AND RENAME IT.
!  IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections
!  need to consider collection efficiencies for different habits?
!  collection efficiency for hail too high.  big hail should not
!  coallesce.

subroutine each_call(m1,dtlt)

  use rconstants, only : &
       pi4,              & ! INTENT(IN)
       alvl,             & ! INTENT(IN)
       alvi,             & ! INTENT(IN)
       alli,             & ! INTENT(IN)
       cliq,             & ! INTENT(IN)
       cice               ! INTENT(IN)

  use micphys, only :    &
       colf,             & ! INTENT(OUT)
       pi4dt,            & ! INTENT(OUT)
       sl,               & ! INTENT(OUT)
       sc,               & ! INTENT(OUT)
       sj,               & ! INTENT(OUT)
       sk,               & ! INTENT(OUT)
       jnmb,             & ! INTENT(IN)
       emb,              & ! INTENT(OUT)
       cfmas,            & ! INTENT(IN)
       parm,             & ! INTENT(IN)
       pwmas,            & ! INTENT(IN)
       jhcat,            & ! INTENT(OUT)
       sh,               & ! INTENT(OUT)
       sm                  ! INTENT(OUT)
       


  implicit none

  ! Arguments
  integer, intent(in) :: m1
  real, intent(in)    :: dtlt

  ! Local Variables
  integer :: lcat, k

  ! Initialize constants for vapor diffusion and, for fixed diameter cases, emb.

  colf = .785 * dtlt
  pi4dt = pi4 * dtlt
  sl(1) = alvl
  sl(2) = alvi
  sc(1) = cliq
  sc(2) = cice
  sj(1) = 0
  sj(2) = 1
  sj(3) = 0
  sj(4) = 0
  sj(5) = 0
  sj(6) = 1
  sj(7) = 1
  sk(1) = alli
  sk(2) = 0.

  do lcat = 1,7
     if (jnmb(lcat) == 2) then
        do k = 2,m1-1
           emb(k,lcat) = cfmas(lcat) * parm(lcat) ** pwmas(lcat)
        enddo
     endif
     do k = 2,m1-1
        jhcat(k,lcat) = lcat
     enddo
  enddo

  do k = 2,m1-1
     sh(k,1) = 0.
     sh(k,2) = 1.
     sh(k,6) = 1.
     sh(k,7) = 1.

     sm(k,1) = 1.
     sm(k,2) = 1.
  enddo

  return
end subroutine each_call

!******************************************************************************

subroutine range_check(m1,k1,k2,k3,i,j,flpw,micro)

  use mem_micro, only : micro_vars ! INTENT(IN) ! Only a type or structure

  use micphys, only : &
       ncat,          & ! INTENT(IN)
       rx,            & ! INTENT(OUT)
       cx,            & ! INTENT(OUT)
       qr,            & ! INTENT(OUT)
       qx,            & ! INTENT(OUT)
       vap,           & ! INTENT(OUT)
       tx,            & ! INTENT(OUT)
       jnmb,          & ! INTENT(IN)
       emb,           & ! INTENT(OUT)
       rxfer,         & ! INTENT(OUT)
       qrxfer,        & ! INTENT(OUT)
       enxfer,        & ! INTENT(OUT)
       cccnx,         & ! INTENT(OUT)
       cifnx            ! INTENT(OUT)

  implicit none

  ! Arguments
  type (micro_vars), intent(in)       :: micro 
  integer, intent(in)                 :: m1, i, j
  real, intent(in) :: flpw
  integer, dimension(10), intent(out) :: k1, k2, k3

  ! Local Variables
  integer                             :: k, lcatt, lcat, l, jcat
  integer                             :: lpw

  ! zero out microphysics scratch arrays for the present i,j column
  lpw=nint(flpw)
  do lcat = 1,ncat
     do k = 2,m1-1
        rx(k,lcat) = 0.
        cx(k,lcat) = 0.
        qr(k,lcat) = 0.
        qx(k,lcat) = 0.
        vap(k,lcat) = 0.
        tx(k,lcat) = 0.
     enddo

     if (jnmb(lcat) >= 3) then
        do k = 2,m1-1
           emb(k,lcat) = 0.
        enddo
     endif

     do jcat = 1,ncat
        do k = 2,m1-1
           rxfer(k,lcat,jcat) = 0.
           qrxfer(k,lcat,jcat) = 0.
           enxfer(k,lcat,jcat) = 0.
        enddo
     enddo
  enddo

  do l = 1,7
     k1(l) = lpw
     k2(l) = 1
  enddo

  ! fill scratch arrays for cloud water

  if (jnmb(1) >= 1) then
     do k = lpw,m1-1
        if (micro%rcp(k,i,j) >= 1.e-9) then
           k2(1) = k
           rx(k,1) = micro%rcp(k,i,j)
           if (jnmb(1) >= 5) cx(k,1) = micro%ccp(k,i,j)
           if (jnmb(1) == 7) cccnx(k) = micro%cccnp(k,i,j)
        else
           if (k2(1) == 1) k1(1) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for rain

  if (jnmb(2) >= 1) then
     do k = lpw,m1-1
        if (micro%rrp(k,i,j) >= 1.e-9) then
           k2(2) = k
           rx(k,2) = micro%rrp(k,i,j)
           qx(k,2) = micro%q2(k,i,j)
           qr(k,2) = qx(k,2) * rx(k,2)
           if (jnmb(2) >= 5) cx(k,2) = micro%crp(k,i,j)
        else
           if (k2(2) == 1) k1(2) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for pristine ice

  if (jnmb(3) >= 1) then
     do k = lpw,m1-1
        if (micro%rpp(k,i,j) >= 1.e-9) then
           k2(3) = k
           rx(k,3) = micro%rpp(k,i,j)
           cx(k,3) = micro%cpp(k,i,j)
           if (jnmb(3) == 7) cifnx(k) = micro%cifnp(k,i,j)
        else
           if (k2(3) == 1) k1(3) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for snow

  if (jnmb(4) >= 1) then
     do k = lpw,m1-1
        if (micro%rsp(k,i,j) >= 1.e-9) then
           k2(4) = k
           rx(k,4) = micro%rsp(k,i,j)
           if (jnmb(4) >= 5) cx(k,4) = micro%csp(k,i,j)
        else
           if (k2(4) == 1) k1(4) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for aggregates

  if (jnmb(5) >= 1) then
     do k = lpw,m1-1
        if (micro%rap(k,i,j) >= 1.e-9) then
           k2(5) = k
           rx(k,5) = micro%rap(k,i,j)
           if (jnmb(5) >= 5) cx(k,5) = micro%cap(k,i,j)
        else
           if (k2(5) == 1) k1(5) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for graupel

  if (jnmb(6) >= 1) then
     do k = lpw,m1-1
        if (micro%rgp(k,i,j) >= 1.e-9) then
           k2(6) = k
           rx(k,6) = micro%rgp(k,i,j)
           qx(k,6) = micro%q6(k,i,j)
           qr(k,6) = qx(k,6) * rx(k,6)
           if (jnmb(6) >= 5) cx(k,6) = micro%cgp(k,i,j)
        else
           if (k2(6) == 1) k1(6) = k + 1
        endif
     enddo
  endif

  ! fill scratch arrays for hail

  if (jnmb(7) >= 1) then
     do k = lpw,m1-1
        if (micro%rhp(k,i,j) >= 1.e-9) then
           k2(7) = k
           rx(k,7) = micro%rhp(k,i,j)
           qx(k,7) = micro%q7(k,i,j)
           qr(k,7) = qx(k,7) * rx(k,7)
           if (jnmb(7) >= 5) cx(k,7) = micro%chp(k,i,j)
        else
           if (k2(7) == 1) k1(7) = k + 1
        endif
     enddo
  endif

  k3(1) = k2(1)
  k3(3) = k2(3)

  k1(8) = min(k1(1),k1(2))
  k2(8) = max(k2(1),k2(2))
  k1(9) = min(k1(3),k1(4),k1(5),k1(6),k1(7))
  k2(9) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
  k1(10) = min(k1(8),k1(9))
  k2(10) = max(k2(8),k2(9))
  return
end subroutine range_check

!******************************************************************************

subroutine each_column(m1,k1,k2,i,j,flpw,rv,dn0)

  use rconstants, only : &
       alvl,             & ! INTENT(IN)
       alvi                ! INTENT(IN)

  use micphys, only :    &
       rvlsair,          & ! INTENT(OUT)
       press,            & ! INTENT(IN)
       tair,             & ! INTENT(IN)
       rvisair,          & ! INTENT(OUT)
       dn0i,             & ! INTENT(OUT)
       tairc,            & ! INTENT(OUT)
       tx,               & ! INTENT(OUT)
       thrmcon,          & ! INTENT(OUT)
       dynvisc,          & ! INTENT(OUT)
       jhcat,            & ! INTENT(OUT)
       jhabtab,          & ! INTENT(IN)
       vapdif,           & ! INTENT(OUT)
       rdynvsci,         & ! INTENT(OUT)
       denfac,           & ! INTENT(OUT)
       colfacr,          & ! INTENT(OUT)
       colf,             & ! INTENT(IN)
       colfacr2,         & ! INTENT(OUT)
       colfacc,          & ! INTENT(OUT)
       colfacc2,         & ! INTENT(OUT)
       tref,             & ! INTENT(OUT)
       sa,               & ! INTENT(OUT)
       tairstrc,         & ! INTENT(IN)
       rvstr,            & ! INTENT(IN)
       sumuy,            & ! INTENT(OUT)
       sumuz,            & ! INTENT(OUT)
       sumvr,            & ! INTENT(OUT)
       rvsref,           & ! INTENT(OUT)
       rvsrefp,          & ! INTENT(OUT)
       rvs0,             & ! INTENT(OUT)
       sh,               & ! INTENT(OUT)
       rxmin               ! INTENT(IN)
  implicit none

  ! Arguments
  integer, intent(in)                :: m1
  integer, intent(in)                :: i, j  ! Not used
  real, intent(in)                   :: flpw
  integer, dimension(10), intent(in) :: k1, k2
  real, dimension(m1), intent(in)    :: rv, dn0
  
  
  ! Local Variables
  integer :: k,nt,ns,lpw
  real :: ck1,ck2,ck3,elsref,elsrefp,dplinv,eisref,eisrefp,dpiinv,relhum

  real :: rslf,rsif,eslf,eslpf,esif,esipf

  data ck1,ck2,ck3/-4.818544e-3,1.407892e-4,-1.249986e-7/

  lpw = nint(flpw)
  
  do k = lpw,m1-1
     rvlsair(k) = rslf (press(k),tair(k))
     rvisair(k) = rsif (press(k),tair(k))
     dn0i(k) = 1. / dn0(k)
     tairc(k)   = tair(k) - 273.16
     tx(k,1) = tairc(k)
     thrmcon(k) = ck1 + (ck2 + ck3 * tair(k)) * tair(k)
     dynvisc(k) = .1718e-4 + .49e-7 * tairc(k)

     ! Diagnose habit of pristine ice and snow

     nt = max(1,min(31,-nint(tairc(k))))
     relhum = min(1.,rv(k) / rvlsair(k))
     ns = max(1,nint(100. * relhum))
     jhcat(k,3) = jhabtab(nt,ns,1)
     jhcat(k,4) = jhabtab(nt,ns,2)

  enddo

  do k = k1(10),k2(10)
     vapdif(k)     = 2.14 * (tair(k) / 273.15) ** 1.94 / press(k)
     rdynvsci(k) = sqrt(1. / dynvisc(k))
     denfac(k) = sqrt(dn0i(k))

     colfacr(k) = colf * denfac(k) * dn0(k)
     colfacr2(k) = 2. * colfacr(k)
     colfacc(k) = colfacr(k) * dn0(k)
     colfacc2(k) = 2. * colfacc(k)

     tref(k,1)   = tairc(k) - min(25.,700. * (rvlsair(k) - rv(k)))
     sa(k,2) = thrmcon(k) * sa(k,1)
     sa(k,3) = thrmcon(k) * (tairstrc(k) + sa(k,1) * rvstr(k))

     sumuy(k) = 0.
     sumuz(k) = 0.
     sumvr(k) = 0.
  enddo

  do k = k1(8),k2(8)
     elsref       = eslf(tref(k,1))
     elsrefp      = eslpf(tref(k,1))
     dplinv       = 1. / (press(k) - elsref)
     rvsref (k,1) = .622 * elsref * dplinv
     rvsrefp(k,1) = .622 * elsrefp * dplinv * (1. + elsref * dplinv)

     sa(k,4) = rvsrefp(k,1) * tref(k,1) - rvsref(k,1)
     sa(k,6) = alvl * rvsrefp(k,1)
     sa(k,8) = alvl * sa(k,4)
  enddo

  do k = k1(9),k2(9)
     tref(k,2)    = min(0.,tref(k,1))
     eisref       = esif (tref(k,2))
     eisrefp      = esipf(tref(k,2))
     dpiinv       = 1. / (press(k) - eisref)
     rvsref (k,2) = .622 * eisref * dpiinv
     rvsrefp(k,2) = .622 * eisrefp * dpiinv * (1. + eisref * dpiinv)
     rvs0(k)      = 379.4 / (press(k) - 610.)

     sa(k,5) = rvsrefp(k,2) * tref(k,2) - rvsref(k,2)
     sa(k,7) = alvi * rvsrefp(k,2)
     sa(k,9) = alvi * sa(k,5)
     sh(k,3) = 0.
     sh(k,4) = 0.
     sh(k,5) = 0.

  enddo

  return
end subroutine each_column

!******************************************************************************

subroutine enemb(m1,k1,k2,lcat,jflag,dn0,i,j)

  use micphys, only : &
       jnmb,          & ! INTENT(IN)
       emb,           & ! INTENT(INOUT)
       cx,            & ! INTENT(OUT)
       rx,            & ! INTENT(IN)
       rxmin,         & ! INTENT(IN)
       jhcat,         & ! INTENT(IN)
       cfemb0,        & ! INTENT(IN)
       pwemb0,        & ! INTENT(IN)
       cfen0,         & ! INTENT(IN)
       dn0i,          & ! INTENT(IN)
       pwen0,         & ! INTENT(IN)
       parm,          & ! INTENT(IN)
       emb0,          & ! INTENT(IN)
       emb1,          & ! INTENT(IN)
       vap,           & ! INTENT(IN)
       enmlttab         ! INTENT(IN)

  implicit none

  ! Arguments
  integer, intent(in)             :: m1, k1, k2, lcat, jflag
  integer, intent(in)             :: i, j ! Not used
  real, dimension(m1), intent(in) :: dn0

  ! Local Variables
  integer :: k,lhcat
  real :: embi,parmi,fracmass,cxloss

  if (jnmb(lcat) == 2) then
     embi = 1. / emb(2,lcat)
     do k = k1,k2
        cx(k,lcat) = rx(k,lcat) * embi
     enddo
  elseif (jnmb(lcat) == 3) then
     do k = k1,k2
        lhcat = jhcat(k,lcat)
        emb(k,lcat) = cfemb0(lhcat) * (dn0(k) * rx(k,lcat)) ** pwemb0(lhcat)
        cx(k,lcat) = cfen0(lhcat) * dn0i(k)  &
             * (dn0(k) * rx(k,lcat)) ** pwen0(lhcat)
     enddo
  elseif (jnmb(lcat) == 4) then
     parmi = 1. / parm(lcat)
     do k = k1,k2
        emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
        cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
     enddo
  elseif (jnmb(lcat) >= 5 .and. jflag == 1) then
     do k = k1,k2
        emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
             / max(1.e-9,cx(k,lcat))))
        cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
     enddo
  elseif (jnmb(lcat) >= 5 .and. jflag == 2) then
     do k = k1,k2

        if (rx(k,lcat) >= 1.e-9) then

           if (vap(k,lcat) < 0.) then
              fracmass = min(1.,-vap(k,lcat) / rx(k,lcat))
              cxloss = cx(k,lcat) * enmlttab(int(200. * fracmass) + 1  &
                   ,jhcat(k,lcat))
              cx(k,lcat) = cx(k,lcat) - cxloss
           endif
           emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
                / max(rxmin(lcat),cx(k,lcat))))
           cx(k,lcat) = rx(k,lcat) / emb(k,lcat)

        endif

     enddo
  endif

  return
end subroutine enemb

!******************************************************************************

subroutine x02(m1,k1,k2,lcat,dn0,i,j)

  use rconstants, only: alli   ! INTENT(IN)

  use micphys, only:  &
       rx,            & ! INTENT(INOUT)
       jhcat,         & ! INTENT(IN)
       vterm,         & ! INTENT(OUT)
       vtfac,         & ! INTENT(IN)
       emb,           & ! INTENT(IN)
       pwvtmasi,      & ! INTENT(IN)
       denfac,        & ! INTENT(IN)
       qx,            & ! INTENT(OUT)
       qr,            & ! INTENT(INOUT)
       cx,            & ! INTENT(INOUT)
       enmlttab,      & ! INTENT(IN)
       dnfac,         & ! INTENT(IN)
       pwmasi,        & ! INTENT(IN)
       gnu,           & ! INTENT(IN)
       shedtab          ! INTENT(IN)


  implicit none

  ! Arguments
  integer, intent(in)                :: m1, lcat
  integer, intent(in)                :: i, j  ! Not Used
  integer, dimension(10),intent(out) :: k1, k2
  real, dimension(m1), intent(in)    :: dn0

  ! Local Variables
  integer :: k,jflag,lhcat,inc,idns
  real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed  &
       ,qrmltshed,shedmass,fracmloss,dn


  k1(lcat) = k1(10)
  k2(lcat) = 1
  do k = k1(10),k2(10)
     if (rx(k,lcat) >= 1.e-9) k2(lcat) = k
     if (k2(lcat) == 1 .and. rx(k,lcat) < 1.e-9) k1(lcat) = k + 1
  enddo

  if (lcat == 2 .or. lcat >= 4) then
     jflag = 1

     call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0,i,j)

     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then

           lhcat = jhcat(k,lcat)
           vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) &
                * denfac(k)

        endif

     enddo
  endif

  if (lcat == 2) then

     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then

           rxinv = 1. / rx(k,lcat)
           qx(k,lcat) = qr(k,lcat) * rxinv
           ! limit rain to under 48C and over -80C
           qx(k,lcat) = max(0.,min(1.6*alli,qx(k,lcat)))

        endif

     enddo

  elseif (lcat == 3) then

     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then


           rinv = 1. / rx(k,lcat)
           qx(k,lcat) = qr(k,lcat) * rinv

           call qtc(qx(k,lcat),tcoal,fracliq)

           rmelt = rx(k,lcat) * fracliq
           cmelt = cx(k,lcat) * fracliq

           rx(k,lcat) = rx(k,lcat) - rmelt
           rx(k,1) = rx(k,1) + rmelt
           cx(k,lcat) = cx(k,lcat) - cmelt
           cx(k,1) = cx(k,1) + cmelt


        endif

     enddo
     !
     ! meyers - source for cloud aerosol number here?
     !
  elseif (lcat == 4 .or. lcat == 5) then

     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then

           rinv = 1. / rx(k,lcat)
           qx(k,lcat) = qr(k,lcat) * rinv
           call qtc(qx(k,lcat),tcoal,fracliq)

           if (fracliq > 1.e-6) then
              rmelt = rx(k,lcat) * fracliq

              ! change this??? move to rain instead ??? look at melting decisions in col2

              ricetor6 = min(rx(k,lcat) - rmelt,rmelt)
              rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
              rx(k,6) = rx(k,6) + rmelt + ricetor6
              qr(k,6) = qr(k,6) + rmelt * alli
              qx(k,lcat) = 0.

              ! keep the above the same with ricetor6
              ! meyers - use sa melt table here? yes
              !
              fracmloss = (rmelt + ricetor6) * rinv
              closs = enmlttab(int(200. * fracmloss) + 1,jhcat(k,lcat)) * cx(k,lcat)
              cx(k,lcat) = cx(k,lcat) - closs
              cx(k,6) = cx(k,6) + closs
           endif


        endif

     enddo


  elseif (lcat == 6) then

     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then

           rxinv = 1. / rx(k,lcat)
           qx(k,lcat) = qr(k,lcat) * rxinv
           call qtc(qx(k,lcat),tcoal,fracliq)

           if (fracliq > 0.95) then
              rx(k,2) = rx(k,2) + rx(k,6)
              qr(k,2) = qr(k,2) + rx(k,6) * alli
              cx(k,2) = cx(k,2) + cx(k,6)
              rx(k,6) = 0.
              qr(k,6) = 0.
              cx(k,6) = 0.
           endif

        endif

     enddo

  elseif (lcat == 7) then

     shedmass = 5.236e-7
     do k = k1(lcat),k2(lcat)

        if (rx(k,lcat) >= 1.e-9) then

           rxinv = 1. / rx(k,lcat)
           qx(k,lcat) = qr(k,lcat) * rxinv
           !c          qx(k,lcat) = max(-50.,qx(k,lcat))
           call qtc(qx(k,lcat),tcoal,fracliq)

           if (fracliq > 0.95) then
              rx(k,2) = rx(k,2) + rx(k,7)
              qr(k,2) = qr(k,2) + rx(k,7) * alli
              cx(k,2) = cx(k,2) + cx(k,7)
              rx(k,7) = 0.
              qr(k,7) = 0.
              cx(k,7) = 0.
              !
              !  take out following IF statement?
              !

           elseif (fracliq > 0.3) then

              lhcat = jhcat(k,lcat)
              inc = nint(200. * fracliq) + 1
              dn = dnfac(lhcat) * emb(k,lcat) ** pwmasi(lhcat)
              idns = max(1,nint(1.e3 * dn * gnu(lcat)))
              rshed = rx(k,lcat) * shedtab(inc,idns)
              !cc               rmltshed = rx(k,lcat) * rmlttab(inc) + rshed
              rmltshed = rshed
              qrmltshed = rmltshed * alli

              rx(k,2) = rx(k,2) + rmltshed
              qr(k,2) = qr(k,2) + qrmltshed
              rx(k,lcat) = rx(k,lcat) - rmltshed
              qr(k,lcat) = qr(k,lcat) - qrmltshed
              !               closs = cx(k,lcat) * enmlttab(inc,lhcat)
              !               cx(k,lcat) = cx(k,lcat) - closs
              !               cx(k,2) = cx(k,2) + closs + rshed / shedmass
              cx(k,2) = cx(k,2) + rshed / shedmass
           endif

        endif

     enddo

  endif
  
  return
end subroutine x02

!******************************************************************************

subroutine c03(m1,k1,k2,lcat,dn0,i,j)

  use micphys, only: jnmb  ! INTENT(IN)

  implicit none

  ! Arguments
  integer, intent(in)             :: m1, k1, k2, lcat, i, j
  real, dimension(m1), intent(in) :: dn0

  ! Local Variables
  integer :: jflag

  jflag = 1
  if (jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,jflag,dn0,i,j)

  return
end subroutine c03

!******************************************************************************

subroutine pc03(m1,k1,k2,lcat,dn0,i,j)

  use micphys, only: &
       jnmb,         & ! INTENT(IN)
       rx,           & ! INTENT(IN)
       jhcat,        & ! INTENT(IN)
       vterm,        & ! INTENT(OUT)
       vtfac,        & ! INTENT(IN)
       emb,          & ! INTENT(IN)
       pwvtmasi,     & ! INTENT(IN)
       denfac          ! INTENT(IN)

  implicit none

  ! Arguments
  integer, intent(in)             :: m1, k1, k2, lcat, i, j
  real, dimension(m1), intent(in) :: dn0

  ! Local Variables
  integer :: k, lhcat, jflag

  jflag = 1
  if (jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,jflag,dn0,i,j)

  do k = k1,k2

     if (rx(k,lcat) >= 1.e-9) then

        lhcat = jhcat(k,lcat)
        vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) * denfac(k)

     endif

  enddo

  return
end subroutine pc03

!******************************************************************************

subroutine sedim(m1,lcat,ngr,nembfall,maxkfall,k1,k2,flpw,i,j  &
     ,rtp,thp,theta,dn0,alphasfc  &
     ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew,pcpfillc,pcpfillr,sfcpcp, &
     dzt, if_adap)

  use rconstants, only: cpi  ! INTENT(IN)

  use micphys, only: &
       nhcat,        & ! INTENT(IN)
       jhcat,        & ! INTENT(IN)
       pcprx,        & ! INTENT(IN)
       rx,           & ! INTENT(INOUT)
       cx,           & ! INTENT(INOUT)
       qx,           & ! INTENT(INOUT)
       ch1,          & ! INTENT(IN)
       emb,          & ! INTENT(IN)
       cfmas,        & ! INTENT(IN)
       ch3,          & ! INTENT(IN)
       dn0i,         & ! INTENT(IN)
       ch2,          & ! INTENT(IN)
       dispemb0,     & ! INTENT(IN)
       accpx,        & ! INTENT(OUT)
       tairc,        & ! INTENT(OUT)
       tair            ! INTENT(IN)

  implicit none

  ! Arguments
  integer, intent(in) :: m1, lcat, ngr, nembfall, maxkfall, k1, k2, i, j
  real, intent(in) :: flpw
  real, dimension(m1), intent(out) :: rtp
  real, dimension(m1), intent(in)  :: thp, theta, dn0
  real, intent(in) :: alphasfc, dtlti
  real, intent(inout) :: pcpg, qpcpg, dpcpg
  real, dimension(m1), intent(out) :: cnew, rnew, qrnew
  real, dimension(m1,maxkfall,nembfall,nhcat), intent(in) :: pcpfillc, pcpfillr
  real, dimension(maxkfall,nembfall,nhcat), intent(in) :: sfcpcp
  real, dimension(m1), intent(in) :: dzt          ! From RAMS 6.0
  integer, intent(in) :: if_adap                  ! From RAMS 6.0

  ! Local Variables
  integer :: k,lhcat,iemb,iemb2,kkf,kk,lpw
  real :: colddn0,rolddn0,qrolddn0,dispemb,riemb,wt2,psfc,qnew

  pcprx(lcat) = 0.
  !do k = 2,k2
  do k = 1, m1 ! From RAMS 6.0
     rnew(k) = 0.
     cnew(k) = 0.
     qrnew(k) = 0.
  enddo

  do k = k1,k2
     lhcat = jhcat(k,lcat)

     if (rx(k,lcat) > 1.e-9) then
        colddn0 = cx(k,lcat) * dn0(k)
        rolddn0 = rx(k,lcat) * dn0(k)
        qrolddn0 = qx(k,lcat) * rolddn0

        dispemb = ch1(lhcat)  &
             * (emb(k,lcat)/cfmas(lhcat)) ** ch3(lhcat) * sqrt(dn0i(k))
        riemb = 1. + ch2(lhcat,ngr) * log10(dispemb / dispemb0(lhcat,ngr))

        !Bob (10/24/00):  Now, limiting iemb to max of nembfall

        iemb = min(nint(riemb),nembfall)

        if (k <= maxkfall) then
           psfc = rolddn0 * sfcpcp(k,iemb,lhcat)
        endif

        do kkf = 1,min(maxkfall,k-1)
           kk = k + 1 - kkf

           cnew(kk)  = cnew(kk)   &
                +  colddn0 * dn0i(kk) * pcpfillc(k,kkf,iemb,lhcat)
           rnew(kk)  = rnew(kk)   &
                +  rolddn0 * dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)
           qrnew(kk) = qrnew(kk)  &
                + qrolddn0 * dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)

           !---------------------------------------------------------------
           ! adjustment for underground and partial grid cells for ada grid

           !if (ada_flag == 1) then
           !   do k = 2,nint(flpw)-1

           !must consider volt correction to pcpfillc and pcpfillr tables

           !---------------------------------------------------------------

        enddo

        if (k <= maxkfall) then
           qpcpg = qpcpg + psfc * qx(k,lcat)
           pcprx(lcat) = pcprx(lcat) + psfc
        endif

     endif
  enddo

  lpw = nint(flpw)

  ! From RAMS 6.0
  if (if_adap == 1) then
     do k = 2,lpw-1
        pcprx(lcat) = pcprx(lcat) + rnew(k) * dn0(k) / dzt(k)
        qpcpg = qpcpg + qrnew(k) * dn0(k) / dzt(k)
        
        cnew(k) = 0.
        rnew(k) = 0.
        qrnew(k) = 0.
     enddo
  endif

  pcpg = pcpg + pcprx(lcat)
  accpx(lcat) = pcprx(lcat)
  dpcpg = dpcpg + pcprx(lcat) * alphasfc
  pcprx(lcat) = pcprx(lcat) * dtlti

  !do k = 2,k2
  do k = lpw,k2 ! From RAMS 6.0
     rtp(k) = rtp(k) + rnew(k) - rx(k,lcat)
     qnew = qrnew(k) / max(1.e-20, rnew(k))

     !         if (iqflag == 1) then
     tairc(k) = tairc(k) - thp(k) * thp(k)  &
          * (2820. * (rnew(k) - rx(k,lcat))  &
          - cpi * (qrnew(k) - qx(k,lcat) * rx(k,lcat)))  &
          / (max(tair(k), 253.) * theta(k))
     !         else
     !            tairc(k) = tairc(k) - thp(k) * thp(k) * 2820.
     !     +          * (rnew(k) - rx(k,lcat)) / (max(tair(k), 253.) * theta(k))
     !         endif

     rx(k,lcat) = rnew(k)
     cx(k,lcat) = cnew(k)
     qx(k,lcat) = qnew

     if (rx(k,lcat) < 1.e-9) then
        rx(k,lcat) = 0.
        cx(k,lcat) = 0.
        qx(k,lcat) = 0.
     endif

  enddo
  return
end subroutine sedim

!******************************************************************************

subroutine negadj1(m1,m2,m3)

  use mem_basic, only: basic_g ! INTENT(OUT)

  use mem_micro, only: micro_g ! INTENT(OUT)

  use mem_grid,  only: &
       grid_g,         & ! INTENT(IN)
       ngrid             ! INTENT(IN)

  use micphys, only:   level ! INTENT(IN)

  use mem_scratch, only: vctr9 ! INTENT(OUT)

  implicit none

  ! Arguments
  integer, intent(in) :: m1, m2, m3

  if (level == 0) return

  call adj1(m1,m2,m3,grid_g(ngrid)%flpw(1,1),basic_g(ngrid)%rtp(1,1,1)  &
       ,basic_g(ngrid)%thp(1,1,1),micro_g(ngrid),vctr9)

  return
end subroutine negadj1

!******************************************************************************

subroutine adj1(m1,m2,m3,flpw,rtp,thp,micro,vctr9)

  use mem_micro, only: &
       micro_vars        ! INTENT(IN)


  use micphys, only:   &
       level,          & ! INTENT(IN)
       ncat,           & ! INTENT(IN)
       rx,             & ! INTENT(OUT)
       cx,             & ! INTENT(OUT)
       jnmb              ! INTENT(IN)

  implicit none

  ! Arguments
  integer, intent(in)                      :: m1, m2, m3
  real, dimension(m2,m3), intent(in)       :: flpw
  type (micro_vars), intent(inout)         :: micro
  real, dimension(m1,m2,m3), intent(inout) :: rtp, thp
  real, dimension(m1), intent(out)         :: vctr9

  ! Local Variables
  integer :: i,j,k,lcat,ka
  real :: frac




  if (level .eq. 0) return

  do lcat = 1,ncat
     do k = 1,m1
        rx(k,lcat) = 0.
        cx(k,lcat) = 0.
     enddo
  enddo

  do j = 1,m3
     do i = 1,m2
        ka = nint(flpw(i,j))

        if (jnmb(1) > 0) call ae1kmic(ka,m1,rx(1,1),micro%rcp(1,i,j))
        if (jnmb(2) > 0) call ae1kmic(ka,m1,rx(1,2),micro%rrp(1,i,j))
        if (jnmb(3) > 0) call ae1kmic(ka,m1,rx(1,3),micro%rpp(1,i,j))
        if (jnmb(4) > 0) call ae1kmic(ka,m1,rx(1,4),micro%rsp(1,i,j))
        if (jnmb(5) > 0) call ae1kmic(ka,m1,rx(1,5),micro%rap(1,i,j))
        if (jnmb(6) > 0) call ae1kmic(ka,m1,rx(1,6),micro%rgp(1,i,j))
        if (jnmb(7) > 0) call ae1kmic(ka,m1,rx(1,7),micro%rhp(1,i,j))

        do lcat = 1,ncat
           do k = ka,m1
              if (rx(k,lcat) < 1.e-9) rx(k,lcat) = 0.
           enddo
        enddo

        do k = ka,m1
           rtp(k,i,j) = max(0.,rtp(k,i,j))
           vctr9(k) = 1.001 * (rx(k,1)+ rx(k,2) + rx(k,3)  &
                + rx(k,4) + rx(k,5) + rx(k,6) + rx(k,7))
        enddo

        do k = ka,m1
           if (vctr9(k) > rtp(k,i,j)) then
              frac = rtp(k,i,j) / (1.e-9 + vctr9(k))
              do lcat = 1,ncat
                 rx(k,lcat) = rx(k,lcat) * frac
              enddo
           endif
        enddo

        if (jnmb(1) > 0) call ae1kmic(ka,m1,micro%rcp(1,i,j),rx(1,1))
        if (jnmb(2) > 0) call ae1kmic(ka,m1,micro%rrp(1,i,j),rx(1,2))
        if (jnmb(3) > 0) call ae1kmic(ka,m1,micro%rpp(1,i,j),rx(1,3))
        if (jnmb(4) > 0) call ae1kmic(ka,m1,micro%rsp(1,i,j),rx(1,4))
        if (jnmb(5) > 0) call ae1kmic(ka,m1,micro%rap(1,i,j),rx(1,5))
        if (jnmb(6) > 0) call ae1kmic(ka,m1,micro%rgp(1,i,j),rx(1,6))
        if (jnmb(7) > 0) call ae1kmic(ka,m1,micro%rhp(1,i,j),rx(1,7))

        if (jnmb(1) >= 5)  &
             call ae1mic(ka,m1,micro%ccp(1,i,j),micro%rcp(1,i,j),rx(1,1))
        if (jnmb(2) >= 5)  &
             call ae1mic(ka,m1,micro%crp(1,i,j),micro%rrp(1,i,j),rx(1,2))
        if (jnmb(3) >= 5)  &
             call ae1mic(ka,m1,micro%cpp(1,i,j),micro%rpp(1,i,j),rx(1,3))
        if (jnmb(4) >= 5)  &
             call ae1mic(ka,m1,micro%csp(1,i,j),micro%rsp(1,i,j),rx(1,4))
        if (jnmb(5) >= 5)  &
             call ae1mic(ka,m1,micro%cap(1,i,j),micro%rap(1,i,j),rx(1,5))
        if (jnmb(6) >= 5)  &
             call ae1mic(ka,m1,micro%cgp(1,i,j),micro%rgp(1,i,j),rx(1,6))
        if (jnmb(7) >= 5)  &
             call ae1mic(ka,m1,micro%chp(1,i,j),micro%rhp(1,i,j),rx(1,7))

        !
        !  Think about how thp should change here - should it be due to a change in
        !     rtp or to a change in the condensate?
        !
        !               vctr10(k) = rrp(k,i,j +rpp(k,i,j)  + rsp(k,i,j) + rap(k,i,j)
        !     +                   + rgp(k,i,j) + rhp(k,i,j)
        !               thp(k,i,j) = thp(k,i,j)
        !     +                    * (1. - aklv * (vctr8(k) - rtp(k,i,j))
        !c or +                    * (1. - aklv * (vctr10(k) - vctr9(k,i,j))
        !     +                    /(max(temp, 253.)))

     enddo
  enddo
  return
end subroutine adj1

!---------------------------------------------------------------------------

subroutine ae1mic(ka,m1,c3,r3,r1)

  implicit none

  ! Arguments
  integer, intent(in)                :: m1,ka
  real, dimension(m1), intent(inout) :: c3
  real, dimension(m1), intent(in)    :: r3, r1

  ! Local Variables
  integer :: k

  do k = ka,m1
     c3(k) = c3(k) * r1(k) / (1.e-9 + r3(k))
     if (c3(k) < 0.) c3(k) = 0.
  enddo

  return
end subroutine ae1mic

!---------------------------------------------------------------------------

subroutine ae1kmic(ka,kb,cr3,cr1)

  implicit none

  ! Arguments
  integer, intent(in) :: ka,kb
  real, dimension(kb), intent(out) :: cr3
  real, dimension(kb), intent(in)  :: cr1

  ! Local Variables
  integer :: k

  do k = ka,kb
     cr3(k) = cr1(k)
  enddo

  return
end subroutine ae1kmic
