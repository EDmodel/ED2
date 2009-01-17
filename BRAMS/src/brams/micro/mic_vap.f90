!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!   08/31/08 - MLO - Switching micro by OLAM equivalent, which uses densities rather than  !
!                    mixing ratio.                                                         !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes some thermodynamic properties (such as temperature and its   !
! derivative with respect to the mixing ratio) for a moist process.                        !
!------------------------------------------------------------------------------------------!
subroutine thrmstr(m1)

   use rconstants
   use micphys
   use therm_lib, only: qtk,thil2temp,dtempdrs

   implicit none
   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in) :: m1
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,lcat
   real :: fracliq,tcoal,t1stguess
   !---------------------------------------------------------------------------------------!

   do k = lpw,m1
      tair(k) = pottemp(k) * exner(k) * cpi
   end do

   !----- Outside the levels in which condensation happens, use non-condensation thermo ---!
   do k = 1,k1(10)-1
      pottemp(k) = thil(k)
      til(k)     = tair(k)
      rvap(k)    = rtot(k)
   end do
   do k = k2(10)+1,m1
      pottemp(k) = thil(k)
      til(k)     = tair(k)
      rvap(k)    = rtot(k)
   end do
   !---------------------------------------------------------------------------------------!

   !----- Initialise til, rliq and rice at the levels in which condensates exist ----------!
   do k = k1(10),k2(10)
      til(k) = thil(k) * exner(k) * cpi
      rliq(k) = 0.
      rice(k) = 0.
   end do

   !----- Cloud and rain are all liquid ---------------------------------------------------!
   do lcat = 1,2
      do k = k1(lcat),k2(lcat)
         rliq(k) = rliq(k) + rx(k,lcat)
      end do
   end do

   !----- Pristine ice, snow, and aggregates are all ice ----------------------------------!
   do lcat = 3,5
      do k = k1(lcat),k2(lcat)
         rice(k) = rice(k) + rx(k,lcat)
      end do
   end do

   !----- Graupel and Hail have both phases, figure out how much of each ------------------!
   do lcat = 6,7
      do k = k1(lcat),k2(lcat)
         call qtk(qx(k,lcat),tcoal,fracliq)
         rliq(k) = rliq(k) + rx(k,lcat) * fracliq
         rice(k) = rice(k) + rx(k,lcat) * (1. - fracliq)
      end do
   end do

   !----- Find the vapour pressure and the total enthalpy due to phase change -------------!
   do k = k1(10),k2(10)
      qhydm(k) = alvl * rliq(k) + alvi * rice(k)
      rvstr(k) = rtot(k) - rliq(k) - rice(k)
   end do

   !----- Finding the air temperature and -dT/drs for this temperature --------------------!
   do k = k1(10),k2(10)
      tairstr(k) = thil2temp(thil(k),exner(k),press(k),rliq(k),rice(k),tair(k))
      sa(k,1) = (-1) * dtempdrs(exner(k),thil(k),tairstr(k),rliq(k),rice(k),1.e-12)
   end do

   return
end subroutine thrmstr
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine diffprep(lcat)

   use rconstants
   use micphys

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: lcat
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,if1,if4,if6,if8,lhcat
   real    :: fre,scdei
   !---------------------------------------------------------------------------------------!

   !----- Selecting whether to use liquid or ice related stuff ----------------------------!
   select case (lcat)
   case (1,2)
      if1 = 1
      if4 = 4
      if6 = 6
      if8 = 8
   case default
      if1 = 2
      if4 = 5
      if6 = 7
      if8 = 9
   end select

   !---------------------------------------------------------------------------------------!
   !     This part is solving the set of equations around Walko et al.(2000), equation 14. !
   !  All the s?(k,lcat) correspond to Walko et al. (2000) ?j, that are defined at their   !
   ! main table.                                                                           !
   !---------------------------------------------------------------------------------------!
   mainloop: do k = k1(lcat),k2(lcat)
      lhcat = jhcat(k,lcat)

      if (rx(k,lcat) < rxmin(lcat)) cycle mainloop

      fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat)                                  &
          + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

      sb(k,lcat)    = cx(k,lcat) * rhoa(k) * fre * pi4dt
      su(k,lcat)    = vapdif(k)  * sb(k,lcat)
      sd(k,lcat)    = sh(k,lcat) * rx(k,lcat)
      se(k,lcat)    = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k)
      sf(k,lcat)    = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)
      sg(k,lcat)    = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3) + sj(lcat) * qr(k,lcat)
      
      scdei         = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
      ss(k,lcat)    = sf(k,lcat) * scdei
      sw(k,lcat)    = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      ttest(k,lcat) = ss(k,lcat) * rvap(k) + sw(k,lcat)
   end do mainloop


   !---------------------------------------------------------------------------------------!
   !    Iced categories. The first trial assumed Mj = 1 and assuming rvap to be the        !
   ! previous rvap. If that gives T above the triple point, then force it to be at the     !
   ! triple point, and Mj becomes 0. Otherwise, happy with the temperature and Mj is set   !
   ! to 1.                                                                                 !
   !---------------------------------------------------------------------------------------!
   select case (lcat)
   case(3:5)
      iceloop: do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) < rxmin(lcat)) cycle iceloop
         !---------------------------------------------------------------------------------!
         !     For all-ice categories, the test is done by assuming Hj to be zero, and the !
         ! variable is switched to 1 in case temperature is greater than triple point.     !
         ! (as explained in Walko et al. (2000) between equations 16 and 17).              !
         !---------------------------------------------------------------------------------!
         if (ttest(k,lcat) >= t3ple-t00 ) then
            sm(k,lcat) = 0.
            sh(k,lcat) = 1.
            sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
            scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
            ss(k,lcat) = sf(k,lcat) * scdei
            sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
         else
            sm(k,lcat) = 1.
         end if
      end do iceloop
   case (6,7)
      mixloop: do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) < rxmin(lcat)) cycle mixloop

         if (ttest(k,lcat) >= t3ple-t00) then
            sm(k,lcat) = 0.
         else
            sm(k,lcat) = 1.
         end if
      end do mixloop

   end select

   lastloop: do k = k1(lcat),k2(lcat)
      if (rx(k,lcat) < rxmin(lcat)) cycle lastloop

      sy(k,lcat) = rvsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
      sz(k,lcat) = 1. - rvsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
      sumuy(k)   = sumuy(k) + su(k,lcat) * sy(k,lcat)
      sumuz(k)   = sumuz(k) + su(k,lcat) * sz(k,lcat)
      
   end do lastloop

   return
end subroutine diffprep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    The subroutine computes the first attempt to find a closed solution for water vapour  !
! in which the fast-response phenomena are implicitly balanced while slow-response         !
! phenomena are represented as explicit forcing terms, and all the effects of heat storage !
! and latent heat are included where appropriate. This is solving equation (17) of Walko   !
! et al. (2000) paper.                                                                     !
!------------------------------------------------------------------------------------------!
subroutine vapdiff ()
   use micphys, only : &
           k1          & ! intent(in)
          ,k2          & ! intent(in)
          ,rvstr       & ! intent(in)
          ,sumuy       & ! intent(in)
          ,sumuz       & ! intent(in)
          ,rvap        ! ! intent(out)
 
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer :: k
   !---------------------------------------------------------------------------------------!


   do k = k1(10),k2(10)
      rvap(k) = (rvstr(k) + sumuy(k)) / (1.0 + sumuz(k))
   end do

   return
end subroutine vapdiff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the Walko et al (2000)'s equation (17) again. depending on   !
! the outcome of the first evaluation. At that point it was assumed that no hydrometeor    !
! had completely evaporated. If equation (9) stated that the category completely evaporat- !
! ed, then Uj and Vj are switched to their alternate values, rvap will be corrected from   !
! equation (17) and rx will be set to zero in that category. The order in which the        !
! hydrometeors are called does matter here.                                                !
!------------------------------------------------------------------------------------------!
subroutine vapflux(lcat)

   use micphys
   use rconstants, only : alli,t3ple,t00,cliqt3
   use micro_coms, only : qmixedmin,qmixedmax,qprismax
   implicit none

   integer, intent(in) :: lcat
   integer             :: k,if1,if4
   real                :: rxx


   !----- Selecting whether to use liquid or ice related stuff ----------------------------!
   select case (lcat)
   case (1:2)
      if1 = 1
      if4 = 4
   case (3:7)
      if1 = 2
      if4 = 5
   end select

   mainloop: do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) < rxmin(lcat)) cycle mainloop
      !------------------------------------------------------------------------------------!
      !    We now solve equation (16) to get the temperature [tx] with the known Mj, and   !
      ! we then compute the first term of the RHS of equation (9) [vap(k,lcat)]            !
      ! sa(k,if4) has rs'(Ref)*T(Ref) - rs(Ref).                                           !
      !------------------------------------------------------------------------------------!
      tx(k,lcat)  = (ss(k,lcat) * rvap(k) + sw(k,lcat)) * sm(k,lcat)                       &
                  + t3ple * (1.-sm(k,lcat))
      vap(k,lcat) = su(k,lcat) * (rvap(k) + sa(k,if4) - rvsrefp(k,if1) * tx(k,lcat))


      !------------------------------------------------------------------------------------!
      !    In case it wasn't completely evaporated.                                        !
      !------------------------------------------------------------------------------------!
      if (vap(k,lcat) > -rx(k,lcat)) then
         rxx = rx(k,lcat) + vap(k,lcat) ! New mixing ratio.
         !----- 
         if (sm(k,lcat) > .5) then
            qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
            qr(k,lcat) = qx(k,lcat) * rxx
         else
            qx(k,lcat) = (rvap(k)*sf(k,lcat)+sg(k,lcat)-tx(k,lcat)*se(k,lcat))       &
                       / sd(k,lcat)
            !----- It seems a bound sanity check, not sure what 350000 and -100000 mean. --!
            !qx(k,lcat) = min(qmixedmax,max(qmixedmin,qx(k,lcat)))
            qr(k,lcat) = qx(k,lcat) * rxx
         end if
      end if

      !------------------------------------------------------------------------------------!
      !     For all hydrometeors, eliminate them in case the solution of equation (9) led  !
      ! to complete evaporation.                                                           !
      !                                                                                    ! 
      ! Bob - Now also do the following section if pristine ice totally melts:evaporate it !
      !       too.                                                                         !
      !------------------------------------------------------------------------------------!
      if ((lcat == 3 .and. qx(k,lcat) > qprismax) .or. vap(k,lcat) <= -rx(k,lcat)) then
         sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
         sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
         sumvr(k) = sumvr(k) + rx(k,lcat)
         rvap(k)  = (rvstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k))

         vap(k,lcat) = - rx(k,lcat)
         tx(k,lcat)  = 0.
         rx(k,lcat)  = 0.
         qx(k,lcat)  = 0.
         qr(k,lcat)  = 0.
         !---- Not sure about this one, but if rx is 0, I guess that cx should also be 0. -!
         cx(k,lcat)  = 0.
      !----- Just copy rxx to rx and we are all set ---------------------------------------!
      elseif (vap(k,lcat) > -rx(k,lcat)) then
         rx(k,lcat) = rxx

      end if

   end do mainloop

   return
end subroutine vapflux
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    I think this subroutine is based on Walko et al. (1995), section 2.7 (other           !
! conversions). This is the first part, conversion between pristine ice and snow due to    !
! secondary ice production based on riming.                                                !
!                                                                                          !
!    In any case, we had to switch most of the internal calculation to double precision,   !
! to avoid underflow/overflow.                                                             !
!------------------------------------------------------------------------------------------!
subroutine psxfer()

   use micphys

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,lhcat,it,uuu
   real(kind=8)        :: embx,dn,xlim,dvap,dqr,dnum
   real(kind=8)        :: rx8,cx8,qx8,vapx8,rxmin8,cxmin8,pwmas8,pwmasi8,dnfac8,dpsmi8
   real(kind=8)        :: gamx8,gam18,dps28,dps8,gnux8,gamn1x8
   !---------------------------------------------------------------------------------------!

   mainloop: do k = k1(3),k2(3)  ! These are now equal to k1(4),k2(4)

      if (vap(k,3) > 0. .or. vap(k,4) < 0.) then
         !----- Assigning a value so in case nothing happens here, nothing will happen... -!
         dvap = 0.
         dqr  = 0.
         dnum = 0.

         if (vap(k,3) > 0. .and. rx(k,3) >= rxmin(3)) then
            lhcat = jhcat(k,3)

            !----- Copying data to double precision scratch variables ---------------------!
            dps8    = dble(dps)
            dps28   = dble(dps2)
            rx8     = dble(rx(k,3))
            cx8     = dble(cx(k,3))
            qx8     = dble(qx(k,3))
            vapx8   = dble(vap(k,3))
            gnux8   = dble(gnu(3))
            gamn1x8 = dble(gamn1(3))
            rxmin8  = dble(rxmin(3))
            cxmin8  = 1.e-3 ! dble(cxmin(3))
            pwmas8  = dble(pwmas(lhcat))
            pwmasi8 = dble(pwmasi(lhcat))
            dnfac8  = dble(dnfac(lhcat))
            dpsmi8  = dble(dpsmi(lhcat))
            !------------------------------------------------------------------------------!
            
            embx    = max(rxmin8,rx8) / max(cxmin8,cx8)
            dn      = dnfac8 * embx**pwmasi8
            it      = nint(dn*1.e6)

            gamx8   = dble(gam(it,3))
            gam18   = dble(gam(it,1))

            xlim = gamx8 * dps28 * (dps8/dn)**(gnux8-1.) / (gamn1x8 * pwmas8 * dn * dn)
            dvap = min(rx8,vapx8 * (xlim + gam18/gamn1x8))
            dqr  = dvap * qx8
            dnum = dvap * min(dpsmi8,1./embx)
         elseif (vap(k,4) < 0. .and. rx(k,4) >= rxmin(4)) then
         
            lhcat = jhcat(k,4)

            !----- Copying data to double precision scratch variables ---------------------!
            dps8    = dble(dps)
            dps28   = dble(dps2)
            rx8     = dble(rx(k,4))
            cx8     = dble(cx(k,4))
            qx8     = dble(qx(k,4))
            vapx8   = dble(vap(k,4))
            gnux8   = dble(gnu(4))
            rxmin8  = dble(rxmin(4))
            cxmin8  = 1.e-3 ! dble(cxmin(4))
            gamn1x8 = dble(gamn1(4))
            pwmas8  = dble(pwmas(lhcat))
            pwmasi8 = dble(pwmasi(lhcat))
            dnfac8  = dble(dnfac(lhcat))
            dpsmi8  = dble(dpsmi(lhcat))
            !------------------------------------------------------------------------------!

            embx    = max(rxmin8,rx8) / max(cxmin8,cx8)
            dn      = dnfac8 * embx**pwmasi8
            it      = nint(dn*1.e6)

            gamx8   = dble(gam(it,3)) 

            xlim = gamx8 * dps28 * (dps8/dn)**(gnux8-1.)/ (gamn1x8 * pwmas8 * dn*dn)
            dvap = max(-rx8,vapx8* xlim)
            dqr = dvap * qx8
            dnum = dvap * max(dpsmi8,1./embx)
         end if

         rx(k,3) = rx(k,3) - sngl(dvap)
         cx(k,3) = cx(k,3) - sngl(dnum)
         qr(k,3) = qr(k,3) - sngl(dqr)
         rx(k,4) = rx(k,4) + sngl(dvap)
         cx(k,4) = cx(k,4) + sngl(dnum)
         qr(k,4) = qr(k,4) + sngl(dqr)
         if (rx(k,3) > rxmin(3)) qx(k,3) = qr(k,3) / rx(k,3)
         if (rx(k,4) > rxmin(4)) qx(k,4) = qr(k,4) / rx(k,4)
      end if

   end do mainloop

   return
end subroutine psxfer
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine updates the temperature and corresponding saturation mixing ratios.   !
!------------------------------------------------------------------------------------------!
subroutine newtemp()

   use rconstants
   use micphys
   use therm_lib, only: rslf,rsif

   implicit none

   !----- Local variables -----------------------------------------------------------------!
   integer             :: k
   !---------------------------------------------------------------------------------------!

   do k = k1(10),k2(10)
      tair(k)    = tairstr(k) + sa(k,1) * (rvstr(k) - rvap(k))
      tairc(k)   = tair(k)    - t00
      pottemp(k) = tair(k) * cp / exner(k)

      rvlsair(k) = rslf(press(k),tair(k))
      rvisair(k) = rsif (press(k),tair(k))
   end do

   return
end subroutine newtemp
!==========================================================================================!
!==========================================================================================!
