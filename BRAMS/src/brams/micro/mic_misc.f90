!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!

subroutine each_call(m1,dtlt)
   !---------------------------------------------------------------------------------------!
   !     FOR ICNFLG=3, DEBATING WHETHER TO KEEP IT #/M4 OR CHANGE PARM TO #/KG/M.  NEED TO !
   ! DEFINE AVMIPSA, ETC. FOR ALL CATEGORIES. MAY WANT TO DEFINE C1 TOO AND RENAME IT.     !
   ! IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections need to consider collection !
   ! efficiencies for different habits? Collection efficiency for hail too high. Big hail  !
   ! should not coallesce.                                                                 !
   !---------------------------------------------------------------------------------------!

   use rconstants, only : &
           pi4            & ! intent(in)
          ,alvl           & ! intent(in)
          ,alvi           & ! intent(in)
          ,alli           & ! intent(in)
          ,cliq           & ! intent(in)
          ,cice           & ! intent(in)
          ,tsupercool     ! ! intent(in)

   use micphys, only :    &
           jnmb           & ! intent(in)
          ,cfmas          & ! intent(in)
          ,parm           & ! intent(in)
          ,pwmas          & ! intent(in)
          ,emb            & ! intent(out)
          ,colf           & ! intent(out)
          ,pi4dt          & ! intent(out)
          ,sl             & ! intent(out)
          ,sc             & ! intent(out)
          ,sj             & ! intent(out)
          ,jhcat          & ! intent(out)
          ,sh             & ! intent(out)
          ,sm             & ! intent(out)
          ,sq             ! ! intent(out)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: m1
   real   , intent(in) :: dtlt
   !----- Local Variables -----------------------------------------------------------------!
   integer :: lcat, lhcat, k
   !---------------------------------------------------------------------------------------!


   !----- Initialize constants for vapor diffusion and, for fixed diameter cases, emb. ----!
   colf  = .785 * dtlt
   pi4dt = pi4  * dtlt

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
   sq(1) = tsupercool
   sq(2) = 0.

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
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_thermovars(m1,i,j,flpw,thp,btheta,pp,rtp,rv,wp,dn0,pi0,micro)

   use mem_micro, only : micro_vars ! INTENT(IN) - micro structure

   use micphys, only : &
        ncat           & ! intent(in)
       ,availcat       & ! intent(in)
       ,progncat       & ! intent(in)
       ,jnmb           & ! intent(in)
       ,rxmin          & ! intent(in)
       ,cxmin          & ! intent(in)
       ,jhcat          & ! intent(in)
       ,k1             & ! intent(out)
       ,k2             & ! intent(out)
       ,k3             & ! intent(out)
       ,lpw            & ! intent(out)
       ,rx             & ! intent(out)
       ,cx             & ! intent(out)
       ,qr             & ! intent(out)
       ,qx             & ! intent(out)
       ,sa             & ! intent(out)
       ,thil           & ! intent(out)
       ,pottemp        & ! intent(out)
       ,til            & ! intent(out)
       ,theiv          & ! intent(out)
       ,rvstr          & ! intent(out)
       ,tair           & ! intent(out)
       ,tairstr        & ! intent(out)
       ,pottemp        & ! intent(out)
       ,qhydm          & ! intent(out)
       ,press          & ! intent(out)
       ,exner          & ! intent(out)
       ,vertvelo       & ! intent(out)
       ,rhoa           & ! intent(out)
       ,rhoi           & ! intent(out)
       ,rvap           & ! intent(out)
       ,rtot           & ! intent(out)
       ,rliq           & ! intent(out)
       ,rice           & ! intent(out)
       ,totcond        & ! intent(out)
       ,vap            & ! intent(out)
       ,tx             & ! intent(out)
       ,emb            & ! intent(out)
       ,rxfer          & ! intent(out)
       ,qrxfer         & ! intent(out)
       ,enxfer         & ! intent(out)
       ,cccnx          & ! intent(out)
       ,cifnx          ! ! intent(out)

   use rconstants, only : p00, cpi, cp, cpor,alvi,alvl, tsupercool, cliqi
   use therm_lib , only : qtk,thil2temp,dtempdrs,thetaeiv

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                         , intent(in ) :: m1, i, j
   real                            , intent(in ) :: flpw
   real             , dimension(m1), intent(in ) :: thp,btheta,pp,rtp,rv,wp,dn0,pi0
   type (micro_vars)               , intent(in ) :: micro 
   !----- Local Variables -----------------------------------------------------------------!
   integer                             :: k, lcatt, lcat, jcat,lhcat
   real                                :: rhomin, frac, tcoal, fracliq
   !---------------------------------------------------------------------------------------!

   !----- Finding the lowest level above ground -------------------------------------------!
   lpw=nint(flpw)

   !----- Initialising the scratch structures ---------------------------------------------!
   do k=1,m1
      thil     (k) = thp   (k)
      exner    (k) = pi0   (k) + pp (k)
      press    (k) = p00 * (cpi * exner(k))**cpor
      vertvelo (k) = wp    (k)
      rhoa     (k) = dn0   (k)
      rhoi     (k) = 1./ rhoa(k)
      rtot     (k) = rtp   (k)
      rliq     (k) = 0.
      rice     (k) = 0.
   end do

   !----- Zero out microphysics scratch arrays for the present i,j column -----------------!
   do lcat = 1,ncat
      do k = 1,m1
         rx(k,lcat)  = 0.
         cx(k,lcat)  = 0.
         qr(k,lcat)  = 0.
         qx(k,lcat)  = 0.
         vap(k,lcat) = 0.
      end do

      if (jnmb(lcat) >= 3) then
         do k = 1,m1-1
            emb(k,lcat) = 0.
         end do
      end if

      do jcat = 1,ncat
         do k = 1,m1-1
            rxfer (k,lcat,jcat) = 0.
            qrxfer(k,lcat,jcat) = 0.
            enxfer(k,lcat,jcat) = 0.
         end do
      end do
   end do

   do lcat = 1,ncat
      k1(lcat) = lpw
      k2(lcat) = 1
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! 1. Fill scratch arrays for cloud water.                                               !
   !---------------------------------------------------------------------------------------!
   if (availcat(1)) then
      do k = lpw,m1-1
         if (micro%rcp(k,i,j) >= rxmin(1)) then
            k2(1) = k
            rx(k,1) = micro%rcp(k,i,j)
            rliq(k) = rliq(k) + rx(k,1)
            if (progncat(1))  cx(k,1)  = micro%ccp(k,i,j)
            if (jnmb(1) == 7) cccnx(k) = micro%cccnp(k,i,j)
         else
            if (k2(1) == 1) k1(1) = k + 1
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 2. Fill scratch arrays for rain.                                                      !
   !---------------------------------------------------------------------------------------!
   if (availcat(2)) then
      do k = lpw,m1-1
         if (micro%rrp(k,i,j) >= rxmin(2)) then
            k2(2) = k
            rx(k,2)    = micro%rrp(k,i,j)
            rliq(k)    = rliq(k) + rx(k,2)
            qx(k,2)    = micro%q2(k,i,j)
            tx(k,2)    = tsupercool + cliqi * qx(k,2)
            qr(k,2)    = qx(k,2) * rx(k,2)
            if (progncat(2)) cx(k,2) = micro%crp(k,i,j)
         else
            if (k2(2) == 1) k1(2) = k + 1
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 3. Fill scratch arrays for pristine ice.                                              !
   !---------------------------------------------------------------------------------------!
   if (availcat(3)) then
      do k = lpw,m1-1
         if (micro%rpp(k,i,j) >= rxmin(3)) then
            k2(3)      = k
            rx(k,3)    = micro%rpp(k,i,j)
            rice(k)    = rice(k) + rx(k,3)
            cx(k,3)    = micro%cpp(k,i,j)
            if (jnmb(3) == 7) cifnx(k) = micro%cifnp(k,i,j)
         else
            if (k2(3) == 1) k1(3) = k + 1
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 4. Fill scratch arrays for snow.                                                      !
   !---------------------------------------------------------------------------------------!
   if (availcat(4)) then
      do k = lpw,m1-1
         if (micro%rsp(k,i,j) >= rxmin(4)) then
            k2(4)      = k
            rx(k,4)    = micro%rsp(k,i,j)
            rice(k)    = rice(k) + rx(k,4)
            if (progncat(4)) cx(k,4) = micro%csp(k,i,j)
         else
            if (k2(4) == 1) k1(4) = k + 1
         endif
      enddo
   endif
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 5. Fill scratch arrays for aggregates.                                                !
   !---------------------------------------------------------------------------------------!
   if (availcat(5)) then
      do k = lpw,m1-1
         if (micro%rap(k,i,j) >= rxmin(5)) then
            k2(5)      = k
            rx(k,5)    = micro%rap(k,i,j)
            rice(k)    = rice(k) + rx(k,5)
            if (progncat(5)) cx(k,5) = micro%cap(k,i,j)
         else
            if (k2(5) == 1) k1(5) = k + 1
         endif
      enddo
   endif
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 6. Fill scratch arrays for graupel.                                                   !
   !---------------------------------------------------------------------------------------!
   if (availcat(6)) then
      do k = lpw,m1-1
         if (micro%rgp(k,i,j) >= rxmin(6)) then
            k2(6)      = k
            rx(k,6)    = micro%rgp(k,i,j)
            qx(k,6)    = micro%q6(k,i,j)
            call qtk(qx(k,6),tx(k,6),fracliq)
            rliq(k)    = rliq(k) + rx(k,6)*fracliq
            rice(k)    = rice(k) + rx(k,6)*(1.-fracliq)
            qr(k,6)    = qx(k,6)          * rx(k,6)
            if (progncat(6)) cx(k,6) = micro%cgp(k,i,j)
         else
            if (k2(6) == 1) k1(6) = k + 1
         endif
      enddo
   endif
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 7. Fill scratch arrays for hail.                                                      !
   !---------------------------------------------------------------------------------------!
   if (availcat(7)) then
      do k = lpw,m1-1
         if (micro%rhp(k,i,j) >= rxmin(7)) then
            k2(7)      = k
            rx(k,7)    = micro%rhp(k,i,j)
            qx(k,7)    = micro%q7(k,i,j)
            call qtk(qx(k,7),tx(k,7),fracliq)
            rliq(k)    = rliq(k) + rx(k,7)*fracliq
            rice(k)    = rice(k) + rx(k,7)*(1.-fracliq)
            qr(k,7)    = qx(k,7)          * rx(k,7)
            if (progncat(7)) cx(k,7) = micro%chp(k,i,j)
         else
            if (k2(7) == 1) k1(7) = k + 1
         endif
      enddo
   endif
   !---------------------------------------------------------------------------------------!

   !----- Setting the limits in which hydrometeors exist ----------------------------------!
   k3(1)  = k2(1) ! k3 saves this initial value for copyback
   k3(3)  = k2(3) ! k3 saves this initial value for copyback
   
   k1(8)  = min(k1(1),k1(2))
   k2(8)  = max(k2(1),k2(2))
   k1(9)  = min(k1(3),k1(4),k1(5),k1(6),k1(7))
   k2(9)  = max(k2(3),k2(4),k2(5),k2(6),k2(7))
   k1(10) = min(k1(8),k1(9))
   k2(10) = max(k2(8),k2(9))
   
   !---------------------------------------------------------------------------------------!
   !    Doing a thermodynamic adjustment. Theta-il and rtp were recently updated, but the  !
   ! other thermodynamic variables (such as theta, rvap) were not and may be off by this   !
   ! time. Defining this values, and also the ice-vapour equivalent potential temperature, !
   ! which will be used to update thil after precipitation (sedimentation) takes place.    !
   !---------------------------------------------------------------------------------------!
   do k=1,m1
      til(k)     = thil(k) * exner(k) * cpi
      totcond(k) = rliq(k) + rice(k)
      rvap(k)    = rtot(k)-totcond(k)
      rvstr(k)   = rvap(k)
      qhydm(k)   = alvl*rliq(k) + alvi *rice(k)
      !----- 1st guess for temperature, then temperature from theta-il and condensates. ---!
      tcoal      = btheta(k) * exner(k) * cpi
      tair(k)    = thil2temp(thil(k),exner(k),press(k),rliq(k),rice(k),tcoal)
      tairstr(k) = tair(k)
      pottemp(k) = cp * tair(k) / exner(k)
      !----- The A1 term in Walko et al. (2000) according to the new thermodynamics. ------!
      sa(k,1) = (-1) * dtempdrs(exner(k),thil(k),tairstr(k),rliq(k),rice(k),1.e-12)

      !------------------------------------------------------------------------------------!
      !    Ice-vapour equivalent potential temperature. This variable is conserved even    !
      ! when sedimentation (aka rain, snow, hail, etc.) occurs. This will be used to       !
      ! compute the new theta-il (ice-liquid potential temperature) in the end of this     !
      ! subroutine.                                                                        !
      !------------------------------------------------------------------------------------!
      theiv(k)   = thetaeiv(thil(k),press(k),tair(k),rvap(k),rtot(k),.true.)
   end do
   
   !---------------------------------------------------------------------------------------!
   !    Initialise temperature for cloud and pristine ice, because they can nucleate and   !
   ! that will require an initial temperature. Snow and aggregates would not need this     !
   ! in principle, but we will play it safe anyway. For rain, graupel, and hail, we use    !
   ! their actual temperature if they already exist, otherwise assume environment          !
   ! temperature. Again, this is probably unecessary but we are just trying to play it     !
   ! safe.                                                                                 !
   !---------------------------------------------------------------------------------------!
   do lcat=1,7
      select case (lcat)
      case (1,3,4,5)
         !----- Hydrometeors with no known temperature. Assume tair. ----------------------!
         do k=1,m1
            tx(k,lcat) = tair(k)
         end do
      case (2,6,7)
         !---------------------------------------------------------------------------------!
         !    The levels with hydrometeors were already initialised using the internal     !
         ! energy. Therefore we only to fill the other levels with tair.                   !
         !---------------------------------------------------------------------------------!
         do k=1,m1
            if (rx(k,lcat) < rxmin(lcat)) then
               tx(k,lcat) = tair(k)
            end if
         end do
      end select
   end do

   return
end subroutine fill_thermovars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will update the temperature- and water-related variables, as the      !
! "thermo" would do. Here we also double check the mixing ratio of all hydrometeors to     !
! ensure they will not lead to negative vapour mixing ratio, and theta-il will be also     !
! updated to account sedimentation.                                                        !
!------------------------------------------------------------------------------------------!
subroutine update_thermo(m1)
   use micphys    , only : theiv         & ! intent(in)    
                         , press         & ! intent(in)    
                         , exner         & ! intent(in)    
                         , ncat          & ! intent(in)    
                         , lpw           & ! intent(in)    
                         , qx            & ! intent(in)    
                         , rxmin         & ! intent(in)    
                         , rtot          & ! intent(inout) 
                         , rx            & ! intent(inout) 
                         , rvap          & ! intent(inout) 
                         , totcond       & ! intent(out)   
                         , rliq          & ! intent(out)   
                         , rice          & ! intent(out)   
                         , thil          & ! intent(out)   
                         , pottemp       & ! intent(out)   
                         , tair          ! ! intent(out)   
   use therm_lib  , only : qtk           & ! subroutine    
                         , thetaeiv2thil & ! function
                         , thil2temp     ! ! function
   use rconstants , only : cp            ! ! intent(in)
   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in) :: m1
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,lcat
   real                :: tcoal, fracliq
   !---------------------------------------------------------------------------------------!
   !     Make sure that vapour and total water are above a minimum, and rescaling mixing   !
   ! ratio and count number if needed.                                                     !
   !---------------------------------------------------------------------------------------!
   totcond(:) = 0.
   do lcat=1,ncat
      do k = lpw,m1
         totcond(k) = totcond(k) + rx(k,lcat)
      end do
   end do
   call minvap_1d_adj(m1,lpw,m1,rtot(1:m1),rvap(1:m1),totcond(1:m1))

   !---------------------------------------------------------------------------------------!
   !    Finding the amount of ice and liquid 
   !---------------------------------------------------------------------------------------!
   rliq(:) = 0.
   rice(:) = 0.
   !----- All liquid hydrometeors ---------------------------------------------------------!
   do lcat=1,2
      do k = lpw,m1
         if (rx(k,lcat) > rxmin(lcat)) rliq(k) = rliq(k) + rx(k,lcat)
      end do
   end do
   !----- All ice hydrometeors ------------------------------------------------------------!
   do lcat=3,5
      do k = lpw,m1
         if (rx(k,lcat) > rxmin(lcat)) rice(k) = rice(k) + rx(k,lcat)
      end do
   end do
   !----- Mixed phase hydrometeors --------------------------------------------------------!
   do lcat=6,7
      do k = lpw,m1
         if (rx(k,lcat) > rxmin(lcat)) then
            call qtk(qx(k,lcat),tcoal,fracliq)
            rliq(k) = rliq(k) + rx(k,lcat)*fracliq
            rice(k) = rice(k) + rx(k,lcat)*(1.-fracliq)
         end if
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !    Theta-E_iv is conserved even with sedimentation. Note this is not exactly the same !
   ! variable as Tripoli and Cotton (1981). Their theta-E_iv is more like a saturation     !
   ! value, so it is not conserved in saturated processes in case sedimentation occurs.    !
   ! This one is based on the saturation point, so in principle this should be conserved   !
   ! even in precipitation events. The only possible caveat is with mixed phases and       !
   ! supercooled water, but I think it should be okay as long as they exist only in excess !
   ! of saturation.                                                                        !
   !    After thil is updated, we also update the potential temperature because it is used !
   ! in the acoustic time step.                                                            !
   !---------------------------------------------------------------------------------------!
   do k = lpw,m1
      thil(k)    = thetaeiv2thil(theiv(k),press(k),rtot(k),.true.)
      tair(k)    = thil2temp(thil(k),exner(k),press(k),rliq(k),rice(k),tair(k))
      pottemp(k) = cp * tair(k) / exner(k)
   end do

   return
end subroutine update_thermo
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine each_column(m1,dtlt)

  use rconstants, only :   &
          pi4              & ! intent(in)
         ,t3ple            & ! intent(in)
         ,es3ple           & ! intent(in)
         ,t00              & ! intent(in)
         ,ep               & ! intent(in)
         ,epes3ple         & ! intent(in)
         ,alvl             & ! intent(in)
         ,alvi               ! intent(in)

  use micphys, only :      &
          k1               & ! intent(in)
         ,k2               & ! intent(in)
         ,lpw              & ! intent(in)
         ,press            & ! intent(in)
         ,tair             & ! intent(in)
         ,jhabtab          & ! intent(in)
         ,colf             & ! intent(in)
         ,tairstr          & ! intent(in)
         ,rvstr            & ! intent(in)
         ,rxmin            & ! intent(in)
         ,rvap             & ! intent(in)
         ,rvlsair          & ! intent(out)
         ,rvisair          & ! intent(out)
         ,rhoa             & ! intent(in)
         ,rhoi             & ! intent(in)
         ,press            & ! intent(in)
         ,colf             & ! intent(out)
         ,pi4dt            & ! intent(out)
         ,tairc            & ! intent(out)
         ,thrmcon          & ! intent(out)
         ,dynvisc          & ! intent(out)
         ,jhcat            & ! intent(out)
         ,vapdif           & ! intent(out)
         ,rdynvsci         & ! intent(out)
         ,denfac           & ! intent(out)
         ,colfacr          & ! intent(out)
         ,colfacr2         & ! intent(out)
         ,colfacc          & ! intent(out)
         ,colfacc2         & ! intent(out)
         ,tref             & ! intent(out)
         ,sa               & ! intent(out)
         ,sumuy            & ! intent(out)
         ,sumuz            & ! intent(out)
         ,sumvr            & ! intent(out)
         ,rvsref           & ! intent(out)
         ,rvsrefp          & ! intent(out)
         ,sh               ! ! intent(out)


    use micro_coms, only : &
           ckcoeff         & ! intent(in)
          ,dvcoeff         & ! intent(in)
          ,vdcoeff         & ! intent(in)
          ,dtempmax        ! ! intent(in)

   use therm_lib  , only : &
           rslf            & ! Function
          ,rsif            & ! Function
          ,rslfp           & ! Function
          ,rsifp           & ! Function
          ,rehuil          ! ! Function



  
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: m1
   real   , intent(in) :: dtlt
   !----- Local Variables -----------------------------------------------------------------!
   integer :: k,nt,ns,irh
   real    :: relhum,pvap
   real    :: tauxkelvin
   !---------------------------------------------------------------------------------------!

   do k = lpw,m1-1
      rvlsair(k) = rslf(press(k),tair(k))
      rvisair(k) = rsif(press(k),tair(k))

      tairc(k)   = tair(k) - t00
      thrmcon(k) = ckcoeff(1) + (ckcoeff(2) + ckcoeff(3) * tair(k)) * tair(k)
      dynvisc(k) = dvcoeff(1) + dvcoeff(2) * tairc(k)
      denfac(k)  = sqrt(rhoi(k))


      !----- Diagnose habit of pristine ice and snow --------------------------------------!
      nt = max(1,min(31,-nint(tairc(k))))
      !----- THERMO DILEMMA : Shouldn't it be rehuil here? --------------------------------!
      relhum = min(1.,rehuil(press(k),tair(k),rvap(k)))

      irh = nint(100.*relhum)
      ns = max(1,irh)
      jhcat(k,3) = jhabtab(nt,ns,1)
      jhcat(k,4) = jhabtab(nt,ns,2)

   end do


   do k = k1(10),k2(10)
      vapdif(k)   = vdcoeff(1) * (tair(k) / t00) ** vdcoeff(2) / press(k)
      rdynvsci(k) = sqrt(1. / dynvisc(k))

      colfacr(k)  = colf * denfac(k) * rhoa(k)
      colfacr2(k) = 2. * colfacr(k)
      colfacc(k)  = colfacr(k) * rhoa(k)
      colfacc2(k) = 2. * colfacc(k)

      tref(k,1)   = tair(k) - min(dtempmax,700. * (rvlsair(k) - rvap(k)))
      sa(k,2)     = thrmcon(k) * sa(k,1)
      sa(k,3)     = thrmcon(k) * (tairstr(k) + sa(k,1) * rvstr(k))

      sumuy(k) = 0.
      sumuz(k) = 0.
      sumvr(k) = 0.
   end do

   !----- Liquid water properties ---------------------------------------------------------!
   do k = k1(8),k2(8)
      rvsref (k,1) = rslf(press(k),tref(k,1))
      rvsrefp(k,1) = rslfp(press(k),tref(k,1))

      sa(k,4)      = rvsrefp(k,1) * tref(k,1) - rvsref(k,1)
      sa(k,6)      = alvl * rvsrefp(k,1)
      sa(k,8)      = alvl * sa(k,4)
   enddo
   !----- Ice properties ------------------------------------------------------------------!
   do k = k1(9),k2(9)
      tref(k,2)    = min(t3ple,tref(k,1))
      rvsref (k,2) = rsif(press(k),tref(k,2))
      rvsrefp(k,2) = rsifp(press(k),tref(k,2))

      sa(k,5)      = rvsrefp(k,2) * tref(k,2) - rvsref(k,2)
      sa(k,7)      = alvi * rvsrefp(k,2)
      sa(k,9)      = alvi * sa(k,5)
      sh(k,3)      = 0.
      sh(k,4)      = 0.
      sh(k,5)      = 0.
   end do


   return
end subroutine each_column
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine enemb(lcat,jflag)

   use micphys, only : &
           k1          & ! intent(in)
          ,k2          & ! intent(in)
          ,rhoa        & ! intent(in)
          ,jnmb        & ! intent(in)
          ,emb         & ! intent(inout)
          ,cx          & ! intent(out)
          ,rx          & ! intent(in)
          ,cxmin       & ! intent(in)
          ,rxmin       & ! intent(in)
          ,jhcat       & ! intent(in)
          ,cfemb0      & ! intent(in)
          ,pwemb0      & ! intent(in)
          ,cfen0       & ! intent(in)
          ,rhoi        & ! intent(in)
          ,pwen0       & ! intent(in)
          ,parm        & ! intent(in)
          ,emb0        & ! intent(in)
          ,emb1        & ! intent(in)
          ,vap         & ! intent(in)
          ,enmlttab      ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)             :: lcat, jflag
   !----- Local Variables -----------------------------------------------------------------!
   integer                         :: k,lhcat,lk1,lk2
   real                            :: embi,parmi,fracmass,cxloss
   !---------------------------------------------------------------------------------------!


   lk1=k1(lcat)
   lk2=k2(lcat)

   select case (jnmb(lcat))
   case(2)
      embi = 1. / emb(2,lcat)
      do k = lk1,lk2
         cx(k,lcat) = rx(k,lcat) * embi
      end do

   case(3)
      do k = lk1,lk2
         lhcat       = jhcat(k,lcat)
         emb(k,lcat) = cfemb0(lhcat)        * (rhoa(k) * rx(k,lcat)) ** pwemb0(lhcat)
         cx(k,lcat)  = cfen0(lhcat)*rhoi(k) * (rhoa(k) * rx(k,lcat)) ** pwen0(lhcat)
      end do

   case (4)
      parmi = 1. / parm(lcat)
      do k = lk1,lk2
         emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
         cx(k,lcat)  = rx(k,lcat) / emb(k,lcat)
      end do

   case (5:7)
      select case (jflag)
      case (1)
         do k = lk1,lk2
             lhcat = jhcat(k,lcat)
             emb(k,lcat) = max(emb0(lcat),min(emb1(lcat)                                   &
                                             ,rx(k,lcat)/max(1.e-9,cx(k,lcat)) ) )
             cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
         end do

      case(2)
         do k = lk1,lk2
            lhcat = jhcat(k,lcat)
            if (rx(k,lcat) >= rxmin(lcat)) then

               if (vap(k,lcat) < 0.) then
                  fracmass = min(1.,-vap(k,lcat) / rx(k,lcat))
                  cxloss = cx(k,lcat) * enmlttab(int(200. * fracmass) + 1,lhcat)
                  cx(k,lcat) = cx(k,lcat) - cxloss
               end if
               !---------------------------------------------------------------------------!
               ! Leaving cxmin commented out for now, I will check the impact later.       !
               !---------------------------------------------------------------------------!
               ! emb(k,lcat) = max(emb0(lcat),min(emb1(lcat)                                 &
               !                                 ,rx(k,lcat)/max(cxmin(lhcat),cx(k,lcat)) ) )
               emb(k,lcat) = max(emb0(lcat),min(emb1(lcat)                                 &
                                               ,rx(k,lcat)/max(1.e-9,cx(k,lcat)) ) )
               cx(k,lcat) = rx(k,lcat) / emb(k,lcat)

            end if

         end do
      end select
   end select

   if (jflag == 2) call getict(lcat)

   return
end subroutine enemb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine x02(lcat)

   use rconstants, only: &
        qicet3,          & ! INTENT(IN)
        qliqt3,          & ! INTENT(IN)
        alli               ! INTENT(IN)

   use micphys, only:  &
        rx,            & ! INTENT(INOUT)
        rxmin,         & ! intent(in)
        jhcat,         & ! INTENT(IN)
        k1,            & ! intent(in)
        k2,            & ! intent(in)
        rhoa,          & ! intent(in)
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

   use therm_lib, only : qtk
   
   use micro_coms, only : &
        qrainmin,         & ! INTENT(IN)
        qrainmax          ! ! INTENT(IN)

   implicit none

   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in)                :: lcat
   !----- Local Variables -----------------------------------------------------------------!
   integer :: k,jflag,lhcat,inc,idns
   real    :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed,qrmelt
   real    :: qrmltshed,shedmass,fracmloss,dn
   !---------------------------------------------------------------------------------------!

   k1(lcat) = k1(10)
   k2(lcat) = 1
   do k = k1(10),k2(10)
      if (rx(k,lcat) >= rxmin(lcat)) k2(lcat) = k
      if (k2(lcat) == 1 .and. rx(k,lcat) < rxmin(lcat)) k1(lcat) = k + 1
   end do

   if (k1(lcat) > k2(lcat)) return


   !----- Finding the terminal velocity ---------------------------------------------------!
   select case (lcat)
   case (2,4:7)
      jflag = 1
      call enemb(lcat,jflag)
      terminalloop: do k = k1(lcat),k2(lcat)

         if (rx(k,lcat) < rxmin(lcat) ) cycle terminalloop

         lhcat = jhcat(k,lcat)
         vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) * denfac(k)

      end do terminalloop
   end select
   !---------------------------------------------------------------------------------------!


   select case (lcat) 
   case (2) !----- Rain -------------------------------------------------------------------!
      do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) >= rxmin(lcat)) then
            rxinv = 1. / rx(k,lcat)
            qx(k,lcat) = qr(k,lcat) * rxinv
            !----- Limit rain to under 48C and over -80C ----------------------------------!
            qx(k,lcat) = max(qrainmin,min(qrainmax,qx(k,lcat)))
            !----- To avoid mass and total energy inconsistencies -------------------------!
            qr(k,lcat) = qx(k,lcat) * rx(k,lcat)
         end if
      end do

   case (3) !----- Pristine Ice -----------------------------------------------------------!
      do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) >= rxmin(lcat)) then
            rinv       = 1. / rx(k,lcat)
            qx(k,lcat) = qr(k,lcat) * rinv

            call qtk(qx(k,lcat),tcoal,fracliq)

            rmelt  = rx(k,lcat)    * fracliq
            cmelt  = cx(k,lcat)    * fracliq
            qrmelt = qliqt3        * rmelt

            rx(k,lcat) = rx(k,lcat) - rmelt
            rx(k,1)    = rx(k,1)    + rmelt
            cx(k,lcat) = cx(k,lcat) - cmelt
            cx(k,1)    = cx(k,1)    + cmelt
            qr(k,lcat) = qr(k,lcat) - qrmelt
            qr(k,1)    = qr(k,1)    + qrmelt
            if (rx(k,lcat) > rxmin(lcat)) qx(k,lcat) = qr(k,lcat) / rx(k,lcat)
            if (rx(k,1)    > rxmin(1)   ) qx(k,1)    = qr(k,1)    / rx(k,1)

         end if
      end do


   case (4,5) !----- Snow, aggregates -----------------------------------------------------!
      do k = k1(lcat),k2(lcat)

         if (rx(k,lcat) >= rxmin(lcat)) then

            rinv = 1. / rx(k,lcat)
            qx(k,lcat) = qr(k,lcat) * rinv
            call qtk(qx(k,lcat),tcoal,fracliq)

            if (fracliq > 1.e-6) then
               rmelt = rx(k,lcat) * fracliq

               !---------------------------------------------------------------------------!
               !    Change this??? Move to rain instead ??? Look at melting decisions in   !
               ! col2.                                                                     !
               !---------------------------------------------------------------------------!
               ricetor6   = min(rx(k,lcat) - rmelt,rmelt)
               rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
               rx(k,6)    = rx(k,6) + rmelt + ricetor6
               qr(k,6)    = qr(k,6) + rmelt * qliqt3 + ricetor6 * qicet3
               qx(k,lcat) = qicet3 !---- All water is gone, leave ice at 0°C only ---------!
               qr(k,lcat) = qx(k,lcat) * rx(k,lcat)
               if (rx(k,6) > rxmin(6)) qx(k,6)    = qr(k,6)    / rx(k,6)

               !----- Keep the above the same with ricetor6 -------------------------------!
               fracmloss  = (rmelt + ricetor6) * rinv
               closs      = enmlttab(int(200. * fracmloss) + 1,jhcat(k,lcat)) * cx(k,lcat)
               cx(k,lcat) = cx(k,lcat) - closs
               cx(k,6)    = cx(k,6) + closs
            end if
         end if
      end do

   case (6) !----- Graupel ----------------------------------------------------------------!
      do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) >= rxmin(lcat)) then
            rxinv = 1. / rx(k,lcat)
            qx(k,lcat) = qr(k,lcat) * rxinv
            call qtk(qx(k,lcat),tcoal,fracliq)

            if (fracliq > 0.95) then
               rx(k,2) = rx(k,2) + rx(k,6)
               qr(k,2) = qr(k,2) + rx(k,6) * qliqt3
               cx(k,2) = cx(k,2) + cx(k,6)
               qx(k,2) = qr(k,2) / rx(k,2)
               rx(k,6) = 0.
               qr(k,6) = 0.
               cx(k,6) = 0.
               qx(k,6) = 0.
            end if
         end if
      end do

   case (7) !----- Hail -------------------------------------------------------------------!
      shedmass = 5.236e-7
      do k = k1(lcat),k2(lcat)
         if (rx(k,lcat) >= rxmin(lcat)) then
            rxinv = 1. / rx(k,lcat)
            qx(k,lcat) = qr(k,lcat) * rxinv
            call qtk(qx(k,lcat),tcoal,fracliq)

            if (fracliq > 0.95) then
               rx(k,2) = rx(k,2) + rx(k,7)
               qr(k,2) = qr(k,2) + rx(k,7) * qliqt3
               cx(k,2) = cx(k,2) + cx(k,7)
               qx(k,2) = qr(k,2) / rx(k,2)
               rx(k,7) = 0.
               qr(k,7) = 0.
               cx(k,7) = 0.
               qx(k,7) = 0.
            !----- Take out following IF statement? ---------------------------------------!
            elseif (fracliq > 0.3) then

               lhcat      = jhcat(k,lcat)
               inc        = nint(200. * fracliq) + 1
               dn         = dnfac(lhcat) * emb(k,lcat) ** pwmasi(lhcat)
               idns       = max(1,nint(1.e3 * dn * gnu(lcat)))
               rshed      = rx(k,lcat) * shedtab(inc,idns)
               rmltshed   = rshed
               qrmltshed  = rmltshed * qliqt3

               rx(k,2)    = rx(k,2)    + rmltshed
               qr(k,2)    = qr(k,2)    + qrmltshed
               rx(k,lcat) = rx(k,lcat) - rmltshed
               qr(k,lcat) = qr(k,lcat) - qrmltshed
               cx(k,2)    = cx(k,2)    + rshed / shedmass
               if (rx(k,lcat) > rxmin(lcat)) qx(k,lcat) = qr(k,lcat) / rx(k,lcat)
               if (rx(k,2)    > rxmin(2)   ) qx(k,2)    = qr(k,2)    / rx(k,2)
            end if
         end if
      end do

   end select

   return
end subroutine x02
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine determines the terminal velocity.                                    !
!------------------------------------------------------------------------------------------!
subroutine pc03(lcat,jflag)

   use micphys, only: &
           k1         & ! intent(in)
          ,k2         & ! intent(in)
          ,jnmb       & ! intent(in)
          ,rx         & ! intent(in)
          ,jhcat      & ! intent(in)
          ,vtfac      & ! intent(in)
          ,emb        & ! intent(in)
          ,pwvtmasi   & ! intent(in)
          ,denfac     & ! intent(in)
          ,rxmin      & ! intent(in)
          ,vterm      ! ! intent(out)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)             :: lcat, jflag
   !----- Local Variables -----------------------------------------------------------------!
   integer :: k, lhcat
   !---------------------------------------------------------------------------------------!

   if (jnmb(lcat) >= 3) call enemb(lcat,jflag)

   mainloop: do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) < rxmin(lcat)) cycle mainloop

      lhcat = jhcat(k,lcat)
      vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) * denfac(k)

   end do mainloop

   return
end subroutine pc03
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sedim(m1,lcat,if_adap,mynum,pcpg,qpcpg,dpcpg,dtlti,pcpfillc,pcpfillr,sfcpcp,dzt)

   use rconstants, only : cpi,ttripoli,alvl,alvi,alli,cp  ! intent(in)
   use micphys   , only : &
           k1             & ! intent(in   )
          ,k2             & ! intent(in   )
          ,nembfall       & ! intent(in   )
          ,maxkfall       & ! intent(in   )
          ,nhcat          & ! intent(in   )
          ,jhcat          & ! intent(in   )
          ,progncat       & ! intent(in   )
          ,pcprx          & ! intent(in   )
          ,ch1            & ! intent(in   )
          ,emb            & ! intent(in   )
          ,cfmasi         & ! intent(in   )
          ,ch3            & ! intent(in   )
          ,rhoi           & ! intent(in   )
          ,ch2            & ! intent(in   )
          ,dispemb0i      & ! intent(in   )
          ,tair           & ! intent(in   )
          ,rhoa           & ! intent(in   )
          ,rhoi           & ! intent(in   )
          ,rtot           & ! intent(in   )
          ,thil           & ! intent(in   )
          ,pottemp        & ! intent(in   )
          ,lpw            & ! intent(in   )
          ,denfac         & ! intent(in   )
          ,rx             & ! intent(inout)
          ,rxmin          & ! intent(inout)
          ,cx             & ! intent(inout)
          ,cxmin          & ! intent(inout)
          ,qx             & ! intent(inout)
          ,qr             & ! intent(inout)
          ,cfall          & ! intent(  out)
          ,rfall          & ! intent(  out)
          ,qrfall         & ! intent(  out)
          ,accpx          & ! intent(  out)
          ,tairc          ! ! intent(  out)
   use micro_coms, only : &
           alphasfc       ! ! intent(in   )

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                                    , intent(in)    :: m1,lcat,if_adap,mynum
   real                                       , intent(in)    :: dtlti
   real, dimension(m1                        ), intent(in)    :: dzt
   real, dimension(m1,maxkfall,nembfall,nhcat), intent(in)    :: pcpfillc, pcpfillr
   real, dimension(   maxkfall,nembfall,nhcat), intent(in)    :: sfcpcp
   real                                       , intent(inout) :: pcpg, qpcpg, dpcpg
   !----- Local Variables -----------------------------------------------------------------!
   integer :: k,lhcat,iemb,iemb2,kkf,kk,jcat,ee
   real    :: dispemb,riemb,wt2,psfc
   real    :: tcoal,fracliq,dqlat,coldrhoa,roldrhoa,qroldrhoa
   !---------------------------------------------------------------------------------------!

   !----- Zero out any "fall" cells that might accumulate precipitation -------------------!
   pcprx(lcat) = 0.
   do k = 1, m1
      rfall(k)  = 0.
      cfall(k)  = 0.
      qrfall(k) = 0.
   end do

   !----- Loop over potential donor cells -------------------------------------------------!
   mainloop: do k = k1(lcat),k2(lcat)
      lhcat = jhcat(k,lcat)

      !----- Jump to end of loop if current cell has little or no hydrometeor mass --------!
      if (rx(k,lcat) < rxmin(lcat)) cycle mainloop

      coldrhoa  = cx(k,lcat) * rhoa(k)
      roldrhoa  = rx(k,lcat) * rhoa(k)
      qroldrhoa = qx(k,lcat) * roldrhoa

      dispemb = ch1(lhcat) * (emb(k,lcat) * cfmasi(lhcat)) ** ch3(lhcat) * denfac(k)
      riemb   = 1. + ch2(lhcat) * log10(dispemb * dispemb0i(lhcat))

      !----- Bob (10/24/00):  Now, limiting iemb to max of nembfall -----------------------!
      iemb = min(nint(riemb),nembfall)

      do kkf = 1,min(maxkfall,k-1)
         kk = k + 1 - kkf

         cfall(kk) = cfall(kk) + coldrhoa  * rhoi(kk) * pcpfillc(k,kkf,iemb,lhcat)
         rfall(kk) = rfall(kk) + roldrhoa  * rhoi(kk) * pcpfillr(k,kkf,iemb,lhcat)
         qrfall(kk) = qrfall(kk) + qroldrhoa * rhoi(kk) * pcpfillr(k,kkf,iemb,lhcat)
      end do

      if (k <= maxkfall) then
         psfc        = roldrhoa * sfcpcp(k,iemb,lhcat)
         qpcpg       = qpcpg + psfc * qx(k,lcat)
         pcprx(lcat) = pcprx(lcat) + psfc
      end if
   end do mainloop


   !----- Copy accumulated precip in "below-ground" cells to surface precip ---------------!
   if (if_adap == 1) then
      do k=2,lpw-1
         pcprx(lcat) = pcprx(lcat) + rfall(k)  * rhoa(k) / dzt(k)
         qpcpg       = qpcpg       + qrfall(k) * rhoa(k) / dzt(k)
         
         cfall(k)  = 0.
         rfall(k)  = 0.
         qrfall(k) = 0.
      end do
   end if

   pcpg        = pcpg  + pcprx(lcat)
   accpx(lcat) = pcprx(lcat)
   dpcpg       = dpcpg + pcprx(lcat) * alphasfc(lcat)
   pcprx(lcat) = pcprx(lcat)  * dtlti

   do k = lpw,k2(lcat) ! From RAMS 6.0
      rtot(k) = rtot(k) + rfall(k) - rx(k,lcat)

      if (rfall(k) > rxmin(lcat)) then
         rx(k,lcat) = rfall(k)
         cx(k,lcat) = cfall(k)
         qx(k,lcat) = qrfall(k) / rfall(k)
      else
         rx(k,lcat) = 0.
         cx(k,lcat) = 0.
         qx(k,lcat) = 0.
      end if
   end do

  return
end subroutine sedim
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine is a wrapper for the microphysics sanity check of hydrometeors mixing  !
! ratio. This is to ensure that all species are either above a minimum value or zero, and  !
! that both vapour and total mixing ratio are above a safe minimum to prevent singular-    !
! ities, memory invasion and other horrible things that overflow and underflow can cause.  !
!------------------------------------------------------------------------------------------!
subroutine negadj1(m1,m2,m3,ia,iz,ja,jz)

   use mem_basic  , only:  &
           basic_g         ! intent(out)

   use mem_micro  , only : &
           micro_g         ! ! intent(out)

   use mem_grid   , only : &
           grid_g          & ! intent(in)
          ,ngrid           ! ! intent(in)

   use therm_lib  , only : &
           vapour_on       ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz

   if (.not. vapour_on) return

   call adj1(m1,m2,m3,ia,iz,ja,jz                                                          &
            , grid_g(ngrid)%flpw                  , basic_g(ngrid)%rv                      &
            , basic_g(ngrid)%rtp                  , basic_g(ngrid)%dn0                     &
            , micro_g(ngrid)                      )

   return
end subroutine negadj1
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine adj1(m1,m2,m3,ia,iz,ja,jz,flpw,rv,rtp,dn0,micro)

   use mem_micro  , only : &
           micro_vars      ! ! INTENT(IN)

   use micphys    , only : &
           ncat            & ! intent(in)
          ,rx              & ! intent(out)
          ,cx              & ! intent(out)
          ,rxmin           & ! intent(in)
          ,cxmin           & ! intent(in)
          ,jnmb            & ! intent(in)
          ,progncat        & ! intent(in)
          ,availcat        ! ! intent(in)

   use mem_scratch, only : &
           vctr6           ! ! intent(out)

   use rconstants , only : &
           toodry          ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                              , intent(in)    :: m1, m2, m3,ia,iz,ja,jz
   real            , dimension(m2,m3)   , intent(in)    :: flpw
   type(micro_vars)                     , intent(inout) :: micro
   real            , dimension(m1,m2,m3), intent(inout) :: rv,rtp,dn0
   !----- Local Variables -----------------------------------------------------------------!
   integer :: i,j,k,lcat,ka
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 1. Initialise all r and c arrays.                                                     !
   !---------------------------------------------------------------------------------------!
   do lcat = 1,ncat
      do k = 1,m1
         rx(k,lcat) = 0.
         cx(k,lcat) = 0.
      enddo
   enddo


   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))

         !---------------------------------------------------------------------------------!
         ! 2. Copying the arrays to scratch structures.                                    !
         !---------------------------------------------------------------------------------!
         do k=ka,m1
            !----- 1. Cloud ---------------------------------------------------------------!
            if (availcat(1)) rx(k,1) = micro%rcp(k,i,j)
            if (progncat(1)) cx(k,1) = micro%ccp(k,i,j)
            !----- 2. Rain ----------------------------------------------------------------!
            if (availcat(2)) rx(k,2) = micro%rrp(k,i,j)
            if (progncat(2)) cx(k,2) = micro%crp(k,i,j)
            !----- 3. Pristine ice --------------------------------------------------------!
            if (availcat(3)) rx(k,3) = micro%rpp(k,i,j)
            if (progncat(3)) cx(k,3) = micro%cpp(k,i,j)
            !----- 4. Snow ----------------------------------------------------------------!
            if (availcat(4)) rx(k,4) = micro%rsp(k,i,j)
            if (progncat(4)) cx(k,4) = micro%csp(k,i,j)
            !----- 5. Aggregates ----------------------------------------------------------!
            if (availcat(5)) rx(k,5) = micro%rap(k,i,j)
            if (progncat(5)) cx(k,5) = micro%cap(k,i,j)
            !----- 6. Graupel -------------------------------------------------------------!
            if (availcat(6)) rx(k,6) = micro%rgp(k,i,j)
            if (progncat(6)) cx(k,6) = micro%cgp(k,i,j)
            !----- 7. Hail ----------------------------------------------------------------!
            if (availcat(7)) rx(k,7) = micro%rhp(k,i,j)
            if (progncat(7)) cx(k,7) = micro%chp(k,i,j)
         end do


         !---------------------------------------------------------------------------------!
         ! 3. Now I will check for very low mixing ratios and flush them to zero.          !
         !---------------------------------------------------------------------------------!
         do lcat = 1,ncat
            do k=ka,m1
               if (rx(k,lcat)  < rxmin(lcat)) then
                  rx(k,lcat) = 0.
                  cx(k,lcat) = 0.
               end if
            end do
         end do
         !---------------------------------------------------------------------------------!
         ! 4. Finding the total condensate mixing ratio.                                   !
         !---------------------------------------------------------------------------------!
         vctr6(:) = 0.
         do lcat=1,ncat
            do k = ka,m1
               vctr6(k) = vctr6(k) + rx(k,lcat)
            end do
         end do

         !---------------------------------------------------------------------------------!
         ! 5. Check whether the condensates would make either rtp or rv fall below toodry. !
         !    If needed, rescale the condensate mixing ratio. rx and cx will be updated    !
         !    as needed inside this subroutine.                                            !
         !---------------------------------------------------------------------------------!
         call minvap_1d_adj(m1,ka,m1,rtp(1:m1,i,j),rv(1:m1,i,j),vctr6(1:m1))

         !---------------------------------------------------------------------------------!
         ! 6. Copy back to the original structures.                                        !
         !---------------------------------------------------------------------------------!
         do k=ka,m1
            !----- 1. Cloud ---------------------------------------------------------------!
            if (availcat(1)) micro%rcp(k,i,j) = rx(k,1)
            if (progncat(1)) micro%ccp(k,i,j) = cx(k,1)
            !----- 2. Rain ----------------------------------------------------------------!
            if (availcat(2)) micro%rrp(k,i,j) = rx(k,2)
            if (progncat(2)) micro%crp(k,i,j) = cx(k,2)
            !----- 3. Pristine ice --------------------------------------------------------!
            if (availcat(3)) micro%rpp(k,i,j) = rx(k,3)
            if (progncat(3)) micro%cpp(k,i,j) = cx(k,3)
            !----- 4. Snow ----------------------------------------------------------------!
            if (availcat(4)) micro%rsp(k,i,j) = rx(k,4)
            if (progncat(4)) micro%csp(k,i,j) = cx(k,4)
            !----- 5. Aggregates ----------------------------------------------------------!
            if (availcat(5)) micro%rap(k,i,j) = rx(k,5)
            if (progncat(5)) micro%cap(k,i,j) = cx(k,5)
            !----- 6. Graupel -------------------------------------------------------------!
            if (availcat(6)) micro%rgp(k,i,j) = rx(k,6)
            if (progncat(6)) micro%cgp(k,i,j) = cx(k,6)
            !----- 7. Hail ----------------------------------------------------------------!
            if (availcat(7)) micro%rhp(k,i,j) = rx(k,7)
            if (progncat(7)) micro%chp(k,i,j) = cx(k,7)
         end do
     end do
  end do
  return
end subroutine adj1
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This is somewhat a sanity check. We don't want rtp to be zero or too small, so we fix !
! the minimum possible to be the "toodry" variable. Here we check this and also sum all    !
! hydrometeors, and check whether this amount of hydrometeors can exist. Not only rtp must !
! be >= toodry, vapour mixing ratio should also be. If these criteria are not met, then we !
! scale down the hydrometeor mixing ratio and count number.                                !
!------------------------------------------------------------------------------------------!
subroutine minvap_1d_adj(m1,ka,kz,rtot,rvap,totcond)
   use micphys    , only : rx     & ! intent(inout)
                         , cx     & ! intent(inout)
                         , rxmin  & ! intent(inout)
                         , ncat   ! ! intent(in)
   use rconstants , only : toodry ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)    :: m1,ka,kz
   real   , dimension(m1), intent(inout) :: rtot,totcond,rvap
   !----- Local Variables -----------------------------------------------------------------!
   integer :: k,lcat
   real    :: frac
   !---------------------------------------------------------------------------------------!

   negadjloop: do k = ka,kz
      rtot(k) = max(toodry,rtot(k))
      if (rtot(k) == toodry) then
         frac    = 0.
      elseif (rtot(k) - totcond(k) <= toodry) then
         frac    = max(0., min(1.,(rtot(k)-toodry)/max(1.e-20,totcond(k)) ) )
      else
         rvap (k) = rtot(k) - totcond(k)
         cycle negadjloop
      end if
      !----- Rescaling the condensates, always checking whether they are above minimum. ---!
      totcond(k) = 0.
      do lcat = 1,ncat
         rx(k,lcat) = rx(k,lcat) * frac
         if (rx(k,lcat) < rxmin(lcat)) then
           rx(k,lcat) = 0.
           cx(k,lcat) = 0.
         else
            cx(k,lcat) = cx(k,lcat) * frac
         end if
         totcond(k) = totcond(k) + rx(k,lcat)
      end do
      rvap (k) = rtot(k) - totcond(k)
   end do negadjloop

   return
end subroutine minvap_1d_adj
!==========================================================================================!
!==========================================================================================!
