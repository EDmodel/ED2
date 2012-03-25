!==========================================================================================!
! grell_cupar_ensemble.f90                                                                 !
!                                                                                          !
!    This file contains subroutines that will calculate ensemble-related stuff.            !
!    Some references that may be useful to follow what is going on here:                   !
!                                                                                          !
!    Grell, G.A., 1993: Prognostic evaluation assumptions used by cumulus parameter-       !
!        izations. Mon. Wea. Rev., vol. 121, 764-787.                                      !
!                                                                                          !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the rates of change per unit of mass at the bottom level. There !
! is an optional mass check here in which the run will halt should mass conservation be    !
! violated. This subroutine generalizes equations B.18/19 from Grell (1993) paper and is   !
! fine for any thermodynamic variable, here represented by "this".                         !
!------------------------------------------------------------------------------------------!
subroutine grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,this,p_cup,this_cup         &
                                  ,mentrd_rate,cdd,dzd_cld,etad_cld,thisd_cld,dellathis)
   use rconstants, only: grav
   implicit none

   integer               , intent(in)    :: mgmzp       ! Number of levels
   logical               , intent(in)    :: checkmass   ! Flag for mass balance check
   real                  , intent(in)    :: masstol     ! Mass tolerance. 
   real                  , intent(in)    :: edt         ! Efficiency, epsilon
   real, dimension(mgmzp), intent(in)    :: this        ! Thermo variable @ model levels
   real, dimension(mgmzp), intent(in)    :: p_cup       ! Pressure at cloud levels [Pa]
   real, dimension(mgmzp), intent(in)    :: this_cup    ! Thermo variable @ cloud levels
   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Downdraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: cdd         ! Downdraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Delta-z for downdrafts;
   real, dimension(mgmzp), intent(in)    :: etad_cld    ! Normalized dndraft mass flux;
   real, dimension(mgmzp), intent(in)    :: thisd_cld   ! Thermo variable at downdraft;
   real                  , intent(inout) :: dellathis   ! Change of thermo per unit of mass

   real                                :: subin       ! Subsidence from level aloft;
   real                                :: detdo1      ! 1st downdraft detrainment term
   real                                :: detdo2      ! 2nd downdraft detrainment term
   real                                :: entdo       ! Downdraft entrainment term
   real                                :: totmass     ! Total mass balance

   detdo1  = edt*etad_cld(2)*cdd(1)*dzd_cld(1)
   detdo2  = edt*etad_cld(1)
   entdo   = edt*etad_cld(2)*mentrd_rate(1)*dzd_cld(1)
   subin   = -edt*etad_cld(2)
   
   if (checkmass) then
      totmass = detdo1+detdo2-entdo+subin
      if (abs(totmass) > masstol) then
         write (unit=*,fmt='(a)')                 '--------- Mass check failed!!! ---------'
         write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
         write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
         write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmass=',  totmass
         write(unit=*, fmt='(a)')                 '----------------------------------------'
         call abort_run('The model will stop since it is not conserving mass...'           &
                       ,'grell_dellabot_ensemble','grell_cupar_ensemble.f90')
      end if
   end if
   
   dellathis = (detdo1 * .5*(thisd_cld(1)+thisd_cld(2)) + detdo2 * thisd_cld(1)            &
             + subin * this_cup(2) - entdo * this(1) ) * grav /(p_cup(1)-p_cup(2))
   
   return
end subroutine grell_dellabot_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the rates of change per unit of mass at the levels other than   !
! the surface. There is an optional mass check here in which the run will halt should mass !
! conservation be violated. This subroutine generalizes equations B.14/15 (entire column)  !
! and B.16/17 (cloud top) from Grell (1993) paper and is fine for any thermodynamic        !
! variable, here represented by "this".                                                    !
!------------------------------------------------------------------------------------------!
subroutine grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,klou,klfc,klod,ktop,this &
                                ,p_cup,this_cup,mentrd_rate,mentru_rate,cdd,cdu,dzd_cld    &
                                ,etad_cld,etau_cld,thisd_cld,thisu_cld,dellathis)
   use rconstants, only : grav

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)    :: mgmzp       ! Number of levels
   logical               , intent(in)    :: checkmass   ! Flag for mass balance check
   real                  , intent(in)    :: masstol     ! Mass tolerance.
   real                  , intent(in)    :: edt         ! Efficiency, epsilon
   integer               , intent(in)    :: kdet        ! Detrainment level
   integer               , intent(in)    :: klou        ! Level of origin of updraft
   integer               , intent(in)    :: klfc        ! Level of free convection
   integer               , intent(in)    :: klod        ! Level of origin of downdraft
   integer               , intent(in)    :: ktop        ! Cloud top
   real, dimension(mgmzp), intent(in)    :: this        ! Thermo variable @ model levels
   real, dimension(mgmzp), intent(in)    :: p_cup       ! Pressure at cloud levels [Pa]
   real, dimension(mgmzp), intent(in)    :: this_cup    ! Thermo variable @ cloud levels
   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Downdraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: cdd         ! Downdraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Delta-z for downdrafts;
   real, dimension(mgmzp), intent(in)    :: etad_cld    ! Normalized dndraft mass flux;
   real, dimension(mgmzp), intent(in)    :: etau_cld    ! Normalized updraft mass flux;
   real, dimension(mgmzp), intent(in)    :: thisd_cld   ! Thermo variable at downdraft;
   real, dimension(mgmzp), intent(in)    :: thisu_cld   ! Thermo variable at updraft;
   real, dimension(mgmzp), intent(inout) :: dellathis   ! Change of thermo per unit of mass
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: k           ! Counter
   real                                  :: subin       ! Subsidence from level aloft;
   real                                  :: subout      ! Subsidence to level below;
   real                                  :: detdo       ! Downdraft detrainment term
   real                                  :: entdo       ! Downdraft entrainment term
   real                                  :: detup       ! Updraft detrainment term
   real                                  :: entup       ! Updraft entrainment term
   real                                  :: totmass     ! Total mass balance
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    The vertical loop is done between the 2nd level and cloud top. The first level was !
   ! solved by grell_dellabot_ensemble, and the cumulus parameterization has nothing to    !
   ! do above the cloud top. Because of the grid staggering dz is always the downdraft dz. !
   !---------------------------------------------------------------------------------------!
   vertloop: do k=2,ktop
      !------------------------------------------------------------------------------------!
      ! Computing the generic terms                                                        !
      !------------------------------------------------------------------------------------!
      subin  = etau_cld(k+1) - edt * etad_cld(k+1)
      subout = etau_cld(k)   - edt * etad_cld(k)
      
      !------------------------------------------------------------------------------------!
      ! Computing updraft terms, depending on where I am.                                  !
      !------------------------------------------------------------------------------------!
      !----- Below the level of free convection, no entrainment or detrainment ------------!
      if (k < klfc .and. k /= klou-1) then 
         entup = 0.
         detup = 0.

      !------------------------------------------------------------------------------------!
      !    Where the updrafts begin, entrainment only. You may ask yourself why only at    !
      ! klou-1, and the levels between klou and klfc-1 are all zero? This is because the   !
      ! net value is zero, since the rates cancel out in this layer.                       !
      !------------------------------------------------------------------------------------!
      elseif (k == klou-1) then
         entup = etau_cld(klou)
         detup = 0.

      !----- In-cloud, both entrainment and detrainment -----------------------------------!
      elseif (k >= klfc .and. k < ktop) then
         entup = mentru_rate(k) * dzd_cld(k) * etau_cld(k)
         detup =       cdu(k+1) * dzd_cld(k) * etau_cld(k)
      
      !----- At the cloud top, detrainment only -------------------------------------------!
      else
         entup = 0.
         subin = 0.          !---- NOTHING enters through the top, not even from above. ---!
         detup = etau_cld(k)
      end if
      !------------------------------------------------------------------------------------!
      
      
      !------------------------------------------------------------------------------------!
      !    Compute downdraft terms, depending on where I am. Note that it's safe to use    !
      ! this for shallow clouds, because it kdet and klod will be both zero, so it will    !
      ! always fall in the "nothing happens" case.                                         !
      !------------------------------------------------------------------------------------!
      !----- Below detrainment level, both entrainment and detrainment happen -------------!
      if (k <= kdet) then
         detdo  = edt *         cdd(k) * dzd_cld(k) * etad_cld(k+1)
         entdo  = edt * mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- Within the downdraft layer, but above kdet, only entrainment happens ---------!
      elseif (k > kdet .and. k < klod) then
         detdo  = 0.
         entdo  = edt * mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- klod requires special assumption otherwise the entrainment would be zero -----!
      elseif (k == klod) then 
         detdo  = 0.
         entdo  = edt * etad_cld(k)
      !----- Outside the downdraft layer, nothing happens ---------------------------------!
      else 
         detdo  = 0.
         entdo  = 0.
      end if
      !------------------------------------------------------------------------------------!
      
      
      !------------------------------------------------------------------------------------!
      !    Checking mass balance if the user wants so. If the mass balance doesn't close   !
      ! within a good tolerance, print the terms and abort the run.                        ! 
      !------------------------------------------------------------------------------------!
      if (checkmass) then
         totmass=subin-subout-entup-entdo+detup+detdo
         if (totmass > masstol) then
            write(unit=*,fmt='(a)')                 '-------- Mass check failed!!! --------'
            write(unit=*,fmt='(2(a,1x,i3,1x))')     ' k      =',     k,'kdet   =',kdet
            write(unit=*,fmt='(2(a,1x,i3,1x))')     ' klou   =',  klou,'klfc   =',klfc
            write(unit=*,fmt='(2(a,1x,i3,1x))')     ' klod   =',  klod,'ktop   =',ktop
            write(unit=*,fmt='(2(a,1x,es10.3,1x))') ' subin  =', subin,'subout =',subout
            write(unit=*,fmt='(2(a,1x,es10.3,1x))') ' detup  =', detup,'entup  =',entup
            write(unit=*,fmt='(2(a,1x,es10.3,1x))') ' entdo  =', entdo,'detdo  =',detdo
            write(unit=*,fmt='(2(a,1x,es10.3,1x))') ' edt    =',   edt,'totmas =',totmass
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' etau_cld(k) = ',etau_cld(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' etad_cld(k) = ',etad_cld(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' etau_cld(kr)= ',etau_cld(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' etad_cld(kr)= ',etad_cld(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' cdu(k)      = ',cdu(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' cdd(k)      = ',cdd(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' cdu(kr)     = ',cdu(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' cdd(kr)     = ',cdd(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' mentru_rate(k)  = ',mentru_rate(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' mentrd_rate(k)  = ',mentrd_rate(k)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' mentru_rate(kr) = ',mentru_rate(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' mentrd_rate(kr) = ',mentrd_rate(k+1)
            write(unit=*,fmt='(1(a,1x,es10.3,1x))') ' dzd_cld(k)      = ',dzd_cld(k)

            write(unit=*,fmt='(a)')                 '--------------------------------------'
            call abort_run('The model will stop since it is not conserving mass...'           &
                          ,'grell_dellas_ensemble','grell_cupar_ensemble.f90')
         end if
      end if
      !------------------------------------------------------------------------------------!
      
      !----- Finding the tendency for our thermodynamic variable --------------------------!
      dellathis(k) = (  subin * this_cup(k+1) - subout * this_cup(k)                       &
                      + detup * .5 * (thisu_cld(k+1) + thisu_cld(k))                       &
                      + detdo * .5 * (thisd_cld(k+1) + thisd_cld(k))                       &
                      - entup * this(k)       - entdo * this(k)      )                     &
                      * grav /(p_cup(k)-p_cup(k+1)) !---- p decreases with height... ------!

   end do vertloop

   return
end subroutine grell_dellas_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes epsilon, the fraction between downdraft and updraft mass     !
! fluxes function associated with downdrafts. This considers that downdraft strength       !
! depends on the wind shear, and the bulk of this subroutine is based on:                  !
!                                                                                          !
! Fritsch, J.M; Chappell, C. F., 1980: Numerical prediction of convectively driven         !
!      mesoscale pressure systems. Part I: Convective parameterization. J. Atmos. Sci.,    !
!      vol. 37(8), 1722-1733. (fc or FC80).                                                !
!                                                                                          !
! Zhang, D.-L.; Fritsch, J.M., 1986: Numerical simulation of the meso-beta scale structure !
!      and evolution of the 1977 Johnstown flood. Part I: Model description and            !
!      verification. J. Atm. Sci., vol. 43(18). 1913-1943. (zf or ZF86).                   !
!------------------------------------------------------------------------------------------!
subroutine grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,klou,klfc,ktop,edtmin,edtmax     &
                                    ,pwav,pwev,z_cup,uwind,vwind,dzd_cld,edt_eff,icld,icap)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)  :: mkx        ! Grid dimesnsion
   integer                    , intent(in)  :: mgmzp      ! Grid dimesnsion
   integer                    , intent(in)  :: maxens_eff ! # of prec. efficiency members
   integer                    , intent(in)  :: klou       ! Updraft origin
   integer                    , intent(in)  :: klfc       ! Level of free convection
   integer                    , intent(in)  :: ktop       ! Cloud top
   integer                    , intent(in)  :: icld       ! Cloud type
   integer                    , intent(in)  :: icap       ! Static control type
   real                       , intent(in)  :: edtmin     ! Minimum efficiency
   real                       , intent(in)  :: edtmax     ! Maximum efficiency
   real                       , intent(in)  :: pwav       ! Integ. updraft cond.       (I1)
   real                       , intent(in)  :: pwev       ! Integ. downdraft evap.    (-I2)
   real, dimension(mgmzp)     , intent(in)  :: z_cup      ! Height @ cloud levels;
   real, dimension(mgmzp)     , intent(in)  :: uwind      ! Zonal wind @ model levels;
   real, dimension(mgmzp)     , intent(in)  :: vwind      ! Meridional wind @ model levels;
   real, dimension(mgmzp)     , intent(in)  :: dzd_cld    ! Layer thickness
   real, dimension(maxens_eff), intent(out) :: edt_eff    ! Prec. efficiency for ensemble.
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: k          ! Counter
   integer                                  :: e          ! Counter
   integer                                  :: iun        ! File unit
   integer                                  :: incr       ! Increment
   real                                     :: botwind    ! Wind magnitude at cloud bottom
   real                                     :: topwind    ! Wind magnitude at cloud top
   real                                     :: vshear     ! Wind shear
   real                                     :: edt        ! Standard prec. Efficiency 
   real                                     :: delta_edt  ! Efficiency "step" for ensemble.
   real                                     :: pef_fc     ! Precipitation efficiency (FC80).
   real                                     :: er_zf      ! Temp. variable to compute...
   real                                     :: pef_zf     ! Prec. efficiency (ZF86)
   real                                     :: zkbc_kft   ! Level of free convection height
   !---------------------------------------------------------------------------------------!
   !      Local constants.                                                                 !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    The coefficients below were extracted from the papers, and they are named after    !
   ! the equation they came from. For example fc58 is the vector with all coefficients     !
   ! from Fritsch-Chappell equation 58.                                                    !
   !---------------------------------------------------------------------------------------!
   real, dimension(0:3)  , parameter     :: fc58=(/1.591,-0.639,0.0953,-0.00496/)
   !----- Factor to convert meters to thoushands of feet, following ZF86 ------------------!
   real                  , parameter     :: m2thfeet=3.2808399e-3
   !---------------------------------------------------------------------------------------!
   !    The paper goes up to the third order only. This however, would lead to values out- !
   ! side the allowable range (0-1) at "normal" LCLs.                                      !
   !---------------------------------------------------------------------------------------!
   real, dimension(0:5)  , parameter     :: zf13=(/    .96729352,-.70034167,.162179896     &
                                                  ,-1.2569798E-2, 4.2772E-4,  -5.44E-6 /)
   !---------------------------------------------------------------------------------------!
   !    Variable that allow extra information to be printed.                               !
   !---------------------------------------------------------------------------------------!
   logical               , parameter     :: debug = .true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    First, estimate the average wind shear between the cloud bottom and cloud top.     !
   ! This is done somewhat differently from the original code, because I think wind shear  !
   ! is not an additive operation, but rather a vector magnitude (considering that it is   !
   ! horizontal vorticity when we neglect the vertical wind gradient):                     !
   !    Notation, based on the closest symbols I could find :)                             !
   !    § means integral; ß means sum; _ means from; ^ means to or exponent                !
   !    Ð?k means ?(k)-?(k-1); U=[u²+v²]^½                                                 !
   !                                                                                       !
   !       Vsm = 1/(ztop-zlou) × §_zlou^ztop [ dU/dz dz]                                   !
   !    This is approximatedly                                                             !
   !                                                                                       !
   !       Vsm = 1/(ztop-zlou) × ß_(k=klou)^ktop [(ÐUk / Ðzk) Ðzk ]                        !
   !    Since Ðzk is never zero, and the intermediate points cancel out, we can define     !
   ! the mean wind shear simply as:                                                        !
   !                                                                                       !
   !       Vsm = [U(ktop)-U(klou)]/(ztop-zlou)                                             !
   !                                                                                       !
   !    Since the wind is not at the cloud level, it needs to be interpolated. Here we'll  !
   ! simply take the average between k and k+1.                                            !
   !---------------------------------------------------------------------------------------!


   !vshear = 0.
   !do k=klou,ktop
   !   vshear = vshear                                                                      &
   !          + abs(sqrt(uwind(k+1)*uwind(k+1)+vwind(k+1)*vwind(k+1))                       &
   !               -sqrt(uwind(k  )*uwind(k  )+vwind(k  )*vwind(k  )))
   !end do
   !vshear  = 1000.*vshear/(z_cup(ktop)-z_cup(klou))

   topwind = 0.5*( sqrt(uwind(ktop     )**2 + vwind(ktop     )**2)                         &
                 + sqrt(uwind(ktop  + 1)**2 + vwind(ktop  + 1)**2))
   botwind = 0.5*( sqrt(uwind(klou     )**2 + vwind(klou     )**2)                         &
                 + sqrt(uwind(klou  + 1)**2 + vwind(klou  + 1)**2))

   !----- Fritsch-Chappell used wind shear in units of 1/(1000*s) -------------------------!
   vshear = 1000.*abs(topwind-botwind)/(z_cup(ktop)-z_cup(klou))
   
   !---------------------------------------------------------------------------------------!
   ! Computing the precipitation efficiency based on FC80 equation (58):                   !
   !---------------------------------------------------------------------------------------!
   if (vshear <= 1.2) then
      pef_fc = 0.95
   elseif (vshear >= 4.0) then
      pef_fc = 0.24
   else
      pef_fc= fc58(0)+ vshear*(fc58(1)+ vshear*(fc58(2)+ fc58(3)*vshear) )
   end if
   !---------------------------------------------------------------------------------------!
   !    Alternative precipitation efficiency, based on ZF86 equations 12 and 13.           !
   !---------------------------------------------------------------------------------------!
   !----- Converting the height at LFC from m to thousands of feet (kft) ------------------!
   zkbc_kft = z_cup(klfc) * m2thfeet
   !----- Using equation (13) to find Er --------------------------------------------------!
   if (zkbc_kft <= 3.) then !----- 5th order polynomial screws here, bound it. ------------!
     er_zf = 0.02
   elseif (zkbc_kft > 24.) then !----- 5th order polynomial screws here, bound it ---------!
     er_zf = 2.40
   else
     er_zf = zf13(0) + zkbc_kft * (zf13(1)+zkbc_kft * (zf13(2)+zkbc_kft *                  &
                                  (zf13(3)+zkbc_kft * (zf13(4) + zkbc_kft * zf13(5)) ) ) )
   end if

   !---------------------------------------------------------------------------------------!
   ! Computing the precipitation efficiency based on ZF86 equation (12):                   !
   !---------------------------------------------------------------------------------------!
   pef_zf = 1./(1. + er_zf)


   !---------------------------------------------------------------------------------------!
   !    The standard precipitation efficiency is simply the average between both methods.  !
   !---------------------------------------------------------------------------------------!
   edt = max(edtmin,min(edtmax,(1. - .5 * (pef_fc+pef_zf))*pwav/max(1.e-20,-pwev)))

   !---------------------------------------------------------------------------------------!
   !    The ensemble will be done with a range of precipitation efficiencies around the    !
   ! main one. Note that there are two ways to define. The first one sounds more intuitive !
   ! to me, but it may give poorer results, thus the second method. I left both here, but  !
   ! that should not be totally controlable by the user.                                   !
   !---------------------------------------------------------------------------------------!
   delta_edt = edt / (float(maxens_eff+1))
   do e=1,maxens_eff
      incr    = (-1)**e * e / 2 !---- This gives 0, 1, -1, 2, -2, etc. --------------------!
      !----- Centred around edt -----------------------------------------------------------!
      edt_eff(e) = max(edtmin,min(edtmax, edt  + incr * delta_edt))
      !----- Shifted from edt -------------------------------------------------------------!
      ! edt_eff(e) = max(edtmin,min(edtmax,edt-float(k-1)*edtinc))
   end do

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   if (debug) then
      iun = 50 + icld
      write (unit=iun,fmt='(a)')          '-----------------------------------------------'
      write (unit=iun,fmt='(a,1x,i5)'   ) ' ICAP    =',icap
      write (unit=iun,fmt='(a)')          ' '
      write (unit=iun,fmt='(a,1x,f11.4)') ' BOTWIND [   m/s] =',botwind
      write (unit=iun,fmt='(a,1x,f11.4)') ' TOPWIND [   m/s] =',topwind
      write (unit=iun,fmt='(a,1x,f11.4)') ' ZKLOU   [     m] =',z_cup(klou)
      write (unit=iun,fmt='(a,1x,f11.4)') ' ZKLFC   [     m] =',z_cup(klfc)
      write (unit=iun,fmt='(a,1x,f11.4)') ' ZKTOP   [     m] =',z_cup(ktop)
      write (unit=iun,fmt='(a,1x,f11.4)') ' VSHEAR  [m/s/km] =',vshear*1000.
      write (unit=iun,fmt='(a,1x,f11.4)') ' PEF_FC  [   ---] =',pef_fc
      write (unit=iun,fmt='(a,1x,f11.4)') ' PEF_ZF  [   ---] =',pef_zf
      write (unit=iun,fmt='(a,1x,f11.4)') ' PWAV    [  g/kg] =',pwav*1000.
      write (unit=iun,fmt='(a,1x,f11.4)') ' PWEV    [  g/kg] =',pwev*1000.
      write (unit=iun,fmt='(a,1x,f11.4)') ' EDT     [   ---] =',edt
      write (unit=iun,fmt='(a)')          '-----------------------------------------------'
      write (unit=iun,fmt='(a)') ' '
      write (unit=iun,fmt='(a)') ' '
   end if
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   
   return
end subroutine grell_efficiency_ensemble
!==========================================================================================!
!==========================================================================================!
