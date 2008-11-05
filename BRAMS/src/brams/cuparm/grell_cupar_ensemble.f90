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
   use rconstants, only: g
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
   real, dimension(mgmzp), intent(inout) :: dellathis   ! Change of thermo per unit of mass

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
   
   dellathis(1) = (detdo1 * .5*(thisd_cld(1)+thisd_cld(2)) + detdo2 * thisd_cld(1)         &
                  + subin * this_cup(2) - entdo * this(1) ) * g /(p_cup(1)-p_cup(2))
   
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
subroutine grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop,this &
                                ,p_cup,this_cup,mentrd_rate,mentru_rate,cdd,cdu,dzd_cld    &
                                ,etad_cld,etau_cld,thisd_cld,thisu_cld,dellathis)
   use rconstants, only : g

   implicit none

   integer               , intent(in)    :: mgmzp       ! Number of levels
   logical               , intent(in)    :: checkmass   ! Flag for mass balance check
   real                  , intent(in)    :: masstol     ! Mass tolerance.
   real                  , intent(in)    :: edt         ! Efficiency, epsilon
   integer               , intent(in)    :: kdet        ! Detrainment level
   integer               , intent(in)    :: k22         ! Updraft origin level
   integer               , intent(in)    :: kbcon       ! Level of free convection
   integer               , intent(in)    :: jmin        ! Downdraft origin level
   integer               , intent(in)    :: ktop        ! Cloud top level
   
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

   integer                               :: k           ! Counter
   real                                  :: subin       ! Subsidence from level aloft;
   real                                  :: subout      ! Subsidence to level below;
   real                                  :: detdo       ! Downdraft detrainment term
   real                                  :: entdo       ! Downdraft entrainment term
   real                                  :: detup       ! Updraft detrainment term
   real                                  :: entup       ! Updraft entrainment term
   real                                  :: totmass     ! Total mass balance

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
      if (k < kbcon .and. k /= k22-1) then 
         entup = 0.
         detup = 0.

      !------------------------------------------------------------------------------------!
      !    Where the updrafts begin, entrainment only. You may ask yourself why only at    !
      ! k22-1, and the levels between k22 and kbcon-1 are all zero? This is because the    !
      ! net value is zero, since the rates cancel out in this layer.                       !
      !------------------------------------------------------------------------------------!
      elseif (k == k22-1) then
         entup = etau_cld(k22)
         detup = 0.

      !----- In-cloud, both entrainment and detrainment -----------------------------------!
      elseif (k >= kbcon .and. k < ktop) then
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
      ! this for shallow clouds, because it kdet and jmin will be both zero, so it will    !
      ! always fall in the "nothing happens" case.                                         !
      !------------------------------------------------------------------------------------!
      !----- Below detrainment level, both entrainment and detrainment happen -------------!
      if (k <= kdet) then
         detdo  = edt *         cdd(k) * dzd_cld(k) * etad_cld(k+1)
         entdo  = edt * mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- Within the downdraft layer, but above kdet, only entrainment happens ---------!
      elseif (k > kdet .and. k < jmin) then
         detdo  = 0.
         entdo  = edt * mentrd_rate(k) * dzd_cld(k) * etad_cld(k+1)
      !----- Jmin requires special assumption otherwise the entrainment would be zero -----!
      elseif (k == jmin) then 
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
            write(unit=*,fmt='(2(a,1x,i3,1x))')     ' k22    =',   k22,'kbcon  =',kbcon
            write(unit=*,fmt='(2(a,1x,i3,1x))')     ' jmin   =',  jmin,'ktop   =',ktop
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
                      * g /(p_cup(k)-p_cup(k+1)) !---- p decreases with height... ---------!

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
subroutine grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,k22,kbcon,ktop,edtmin,edtmax     &
                                    ,pwav,pwev,z_cup,uwind,vwind,dzd_cld,edt_eff)

   implicit none

   integer                 , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer                 , intent(in)  :: maxens_eff  ! # of prec. efficiency members
   integer                 , intent(in)  :: k22         ! Updraft origin
   integer                 , intent(in)  :: kbcon       ! Level of free convection
   integer                 , intent(in)  :: ktop        ! Cloud top


   real                       , intent(in)  :: edtmin   ! Minimum efficiency
   real                       , intent(in)  :: edtmax   ! Maximum efficiency
   real                       , intent(in)  :: pwav     ! Integrated updraft cond.    (I1)
   real                       , intent(in)  :: pwev     ! Integrated downdraft evap. (-I2)
   real, dimension(mgmzp)     , intent(in)  :: z_cup    ! Height @ cloud levels;
   real, dimension(mgmzp)     , intent(in)  :: uwind    ! Zonal wind @ model levels;
   real, dimension(mgmzp)     , intent(in)  :: vwind    ! Meridional wind @ model levels;
   real, dimension(mgmzp)     , intent(in)  :: dzd_cld  ! Layer thickness
   real, dimension(maxens_eff), intent(out) :: edt_eff  ! Prec. efficiency for ensemble.
   
   integer                               :: k,e         ! Counters
   integer                               :: incr        ! Increment
   real                                  :: botwind     ! Wind magnitude at cloud bottom;
   real                                  :: topwind     ! Wind magnitude at cloud top;
   real                                  :: vshear      ! Wind shear
   real                                  :: edt         ! Standard prec. Efficiency 
   real                                  :: delta_edt   ! Efficiency "step" for ensemble.
   real                                  :: pef_fc      ! Precipitation efficiency (FC80).
   real                                  :: er_zf       ! Temporary variable to compute...
   real                                  :: pef_zf      ! Prec. efficiency (ZF86)
   real                                  :: zkbc_kft    ! Level of free convection height
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


   !---------------------------------------------------------------------------------------!
   !    First, estimate the average wind shear between the cloud bottom and cloud top.     !
   ! This is done somewhat differently from the original code, because I think wind shear  !
   ! is not an additive operation, but rather a vector magnitude (considering that it is   !
   ! horizontal vorticity when we neglect the vertical wind gradient):                     !
   !    Notation, based on the closest symbols I could find :)                             !
   !    § means integral; ß means sum; _ means from; ^ means to or exponent                !
   !    Ð?k means ?(k)-?(k-1); U=[u²+v²]^½                                                 !
   !                                                                                       !
   !       Vsm = 1/(ztop-z22) × §_z22^ztop [ dU/dz dz]                                     !
   !    This is approximatedly                                                             !
   !                                                                                       !
   !       Vsm = 1/(ztop-z22) × ß_(k=k22)^ktop [(ÐUk / Ðzk) Ðzk ]                          !
   !    Since Ðzk is never zero, and the intermediate points cancel out, we can define     !
   ! the mean wind shear simply as:                                                        !
   !                                                                                       !
   !       Vsm = [U(ktop)-U(k22)]/(ztop-z22)                                               !
   !                                                                                       !
   !    Since the wind is not at the cloud level, it needs to be interpolated. Here we'll  !
   ! simply take the average between k and k+1.                                            !
   !---------------------------------------------------------------------------------------!


   !vshear = 0.
   !do k=k22,ktop
   !   vshear = vshear                                                                      &
   !          + abs(sqrt(uwind(k+1)*uwind(k+1)+vwind(k+1)*vwind(k+1))                       &
   !               -sqrt(uwind(k  )*uwind(k  )+vwind(k  )*vwind(k  )))
   !end do
   !vshear  = 1000.*vshear/(z_cup(ktop)-z_cup(k22))

   topwind = 0.5*( sqrt(uwind(ktop    )**2 + vwind(ktop    )**2)                           &
                 + sqrt(uwind(ktop + 1)**2 + vwind(ktop + 1)**2))
   botwind = 0.5*( sqrt(uwind(k22     )**2 + vwind(k22     )**2)                           &
                 + sqrt(uwind(k22  + 1)**2 + vwind(k22  + 1)**2))

   !----- Fritsch-Chappell used wind shear in units of 1/(1000*s) -------------------------!
   vshear = 1000.*abs(topwind-botwind)/(z_cup(ktop)-z_cup(k22))
   
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
   zkbc_kft = z_cup(kbcon) * m2thfeet
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

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=57,fmt='(a)') '-------------------------------------------------------------'
   !write (unit=57,fmt='(5(a9,1x),a12,2(1x,a8))')                                           &
   !     '  botwind','  topwind','   z(k22)','z(kbcon)','  z(ktop)','  wind_shear'          &
   !     ,'  pef_fc','  pef_zf'
   !write (unit=57,fmt='(5(f9.3,1x),es12.5,2(1x,f8.4))')                                    &
   !        botwind ,topwind ,z_cup(k22),z_cup(kbcon),z_cup(ktop),vshear,pef_fc ,pef_zf
   !write (unit=57,fmt='(a)') ' '
   !write (unit=57,fmt='(2(a,1x,es11.4,1x))') 'pwav=',pwav,'pwev=',pwev
   !write (unit=57,fmt='(a)') '-------------------------------------------------------------'
   !write (unit=57,fmt='(a)') ' '
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


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
   
   return
end subroutine grell_efficiency_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the mass flux associated with the dynamic control. This      !
! should be called inside the other loops, so the ensemble variables that should have      !
! the permutation of all ensemble dimensions is fully filled.                              !
!------------------------------------------------------------------------------------------!
subroutine grell_dyncontrol_ensemble(closure_type,comp_down,mgmzp,maxens_dyn,dtime         &
                                    ,tscal_kf,dens_curr,prev_dnmf,mconv,k22,kbcon,ktop     &
                                    ,omeg,p_cup,edt,mbprime,one_b,aatot0,aatot,x_aatot     &
                                    ,pwav,pwev,dnmf_dyn,upmf_dyn)
   use rconstants, only: g
   use grell_coms, only: pclim   & ! Levels with available climatological work functions
                        ,aclim1  & ! Standard cloud work function climatology
                        ,aclim2    ! Alternative cloud work function climatology
   implicit none
   !----- Input flags ---------------------------------------------------------------------!
   character(len=2)           , intent(in)  :: closure_type ! My dynamic control
   logical                    , intent(in)  :: comp_down    ! I have downdrafts and rain.
   !----- Input dimensions ----------------------------------------------------------------!
   integer                    , intent(in)  :: mgmzp        ! # of levels
   integer                    , intent(in)  :: maxens_dyn   ! # of dynamic control members
   real                       , intent(in)  :: dtime        ! Current grid time step
   real                       , intent(in)  :: tscal_kf     ! Kain-Fritsch time scale
   real                       , intent(in)  :: dens_curr    ! Density current
   real                       , intent(in)  :: prev_dnmf    ! Previous dndraft mass flux
   real                       , intent(in)  :: mconv        ! Integ. moisture convergence
   !----- Input levels --------------------------------------------------------------------!
   integer                    , intent(in)  :: k22          ! Updraft origin
   integer                    , intent(in)  :: kbcon        ! Level of free convection
   integer                    , intent(in)  :: ktop         ! Cloud top
   !----- Input large scale variables -----------------------------------------------------!
   real, dimension(mgmzp)     , intent(in)  :: omeg         ! Vertical pressure velocity
   real, dimension(mgmzp)     , intent(in)  :: p_cup        ! Pressure with forcing
   !----- Input ensemble-dependent variables ----------------------------------------------!
   real                       , intent(in)  :: edt          ! Epsilon, I1/I2
   real                       , intent(in)  :: mbprime      ! Arbitrary mass flux
   real                       , intent(in)  :: one_b        ! (1-b) from Krishnamurti '83.
   real                       , intent(in)  :: pwav         ! Available condensed water
   real                       , intent(in)  :: pwev         ! Water evaporated by ddrfts
   real                       , intent(in)  :: aatot0       ! Current cloud work function
   real                       , intent(in)  :: aatot        ! Cloud work with forcing
   real                       , intent(in)  :: x_aatot      ! Perturbed cloud work fctn.
   !----- Downdraft is inout because I may not compute it ---------------------------------!
   real, dimension(maxens_dyn), intent(inout) :: dnmf_dyn     ! Ref. downdraft mass flux
   real, dimension(maxens_dyn), intent(out)   :: upmf_dyn     ! Ref. updraft mass flux
   !----- Local variables -----------------------------------------------------------------!
   integer           :: kclim        ! Closest climatological level
   integer           :: kclim2nd     ! Next-to-closest one...
   real              :: divisor            ! Scratch, with the divisor
   real, parameter   :: tinyden=epsilon(1.)! Small number for division
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   select case (closure_type)
   !---------------------------------------------------------------------------------------!
   ! 1. Standard Grell (1993), modified quasi-equilibrium buoyant energy.                  !
   !---------------------------------------------------------------------------------------!
   case ('gr')
      divisor     = sign(max(tinyden,abs(aatot-x_aatot)),aatot-x_aatot)
      upmf_dyn(1) = max(0.,mbprime*(aatot-aatot0)/divisor) + dens_curr
   !---------------------------------------------------------------------------------------!
   ! 2. Standard Arakawa and Schubert (1974), quasi-equilibrium buoyant energy.            !
   !---------------------------------------------------------------------------------------!
   case ('as')
      !----- Finding the climatological height which is the closest to my cloud -----------!
      divisor     = sign(max(tinyden,abs(aatot-x_aatot)),aatot-x_aatot)
      kclim       = minloc(abs(pclim-p_cup(ktop)),dim=1)
      upmf_dyn(1) = max(0.,mbprime * (aatot-aclim1(kclim))/divisor) + dens_curr

   !---------------------------------------------------------------------------------------!
   ! 3. Standard Kain and Fritsch (1990), instability removal.                             !
   !---------------------------------------------------------------------------------------!
   case ('kf')
      divisor     = sign(max(tinyden,abs(aatot-x_aatot)),aatot-x_aatot)
      upmf_dyn(1) = max(0., mbprime * dtime * aatot0 / (divisor * tscal_kf)) + dens_curr
   !---------------------------------------------------------------------------------------!
   ! 4. Standard Frank and Cohen (1987), low-level environment mass flux                   !
   !---------------------------------------------------------------------------------------!
   case ('lo')
      upmf_dyn(1) = max(0.,-omeg(k22)/g - prev_dnmf) + dens_curr

   !---------------------------------------------------------------------------------------!
   ! 5. Standard Krishnamurti et al. (1983), moisture convergence.                         !
   !---------------------------------------------------------------------------------------!
   case ('mc')
      divisor     = sign(max(tinyden,abs(pwav- edt * pwev)),pwav-edt*pwev)
      upmf_dyn(1) = max(0.,mconv * 1.4 * one_b / divisor) + dens_curr

   !---------------------------------------------------------------------------------------!
   ! 6. Ensemble, based on Grell and Dévényi (2002)                                        !
   !    Here for each different dynamic control style, the first element is the standard,  !
   !    like in the previous cases, and the others are perturbations.                      !
   !---------------------------------------------------------------------------------------!
   case ('en','nc')
      !----- Divisor is the same for Grell, Arakawa-Schubert and Kain-Fritsch. ------------!
      divisor     = sign(max(tinyden,abs(aatot-x_aatot)),aatot-x_aatot)

      !------------------------------------------------------------------------------------!
      ! 6a. Grell (1993), modified quasi-equilibrium buoyant energy.                       !
      !------------------------------------------------------------------------------------!
      upmf_dyn(1) = max(0.,mbprime*(aatot-aatot0)/divisor)
      upmf_dyn(2) = 0.9 * upmf_dyn(1) + dens_curr
      upmf_dyn(3) = 1.1 * upmf_dyn(1) + dens_curr
      upmf_dyn(1) = upmf_dyn(1)       + dens_curr
      
      !------------------------------------------------------------------------------------!
      ! 6b. Arakawa and Schubert (1974), quasi-equilibrium buoyant energy.                 !
      !------------------------------------------------------------------------------------!
      !------ Finding the closest climatological level ------------------------------------!
      kclim = minloc(abs(pclim-p_cup(ktop)),dim=1)
      !------ Finding the next-to-closest climatological level ----------------------------!
      kclim2nd = minloc(abs(pclim-p_cup(ktop)),dim=1,mask=(pclim /= pclim(kclim)))
      !------ Computing the upward mass flux ----------------------------------------------!
      upmf_dyn(4)= max(0., mbprime*(aatot-aclim1(kclim))/divisor) + dens_curr
      upmf_dyn(5)= max(0., mbprime*(aatot-aclim2(kclim))/divisor) + dens_curr
      upmf_dyn(6)= max(0., mbprime*(aatot-aclim1(kclim2nd))/divisor) + dens_curr
      upmf_dyn(7)= max(0., mbprime*(aatot-aclim2(kclim2nd))/divisor) + dens_curr
      
      !------------------------------------------------------------------------------------!
      ! 6c. Kain and Fritsch (1990), instability removal.                                  !
      !------------------------------------------------------------------------------------!
      upmf_dyn(8)  = max(0., mbprime * dtime * aatot0 / (divisor * tscal_kf))
      upmf_dyn(9)  = upmf_dyn(8) * 1.111111111 + dens_curr !---- Equiv. to tscal_kf * 0.9 -!
      upmf_dyn(10) = upmf_dyn(8) * 0.909090909 + dens_curr !---- Equiv. to tscal_kf * 1.1 -!
      upmf_dyn(8)  = upmf_dyn(8) + dens_curr

      if (closure_type == 'en') then
         !---------------------------------------------------------------------------------!
         ! 6d. Frank and Cohen (1987), low-level environment mass flux.                    !
         !---------------------------------------------------------------------------------!
         upmf_dyn(11) = max(0.,-omeg(k22)/g   - prev_dnmf) + dens_curr
         upmf_dyn(12) = max(0.,-omeg(kbcon)/g - prev_dnmf) + dens_curr
         !----- Picking up the strongest mass flux below the LFC (except k22) -------------!
         upmf_dyn(13) = max(0.,-minval(omeg(1:(kbcon-1)),omeg(1:(kbcon-1)) /= omeg(k22))/g &
                               -prev_dnmf) + dens_curr
         
         !---------------------------------------------------------------------------------!
         ! 6e. Krishnamurti et al. (1983), moisture convergence. Here I am using some      !
         !     different values of (1+eta) (what Grell and Devenyi (2002) called f_emp).   !
         !     I'm oscillating it between (1+0.25) and (1+0.55), which is respectively     !
         !     roughly 0.9 and 1.1 times (1+0.4), their average value                      !
         !---------------------------------------------------------------------------------!
         divisor      = sign(max(tinyden,abs(pwav- edt * pwev)),pwav-edt*pwev)
         upmf_dyn(14) = max(0.,mconv * 1.4 * one_b / divisor)
         upmf_dyn(15) = 1.1 * upmf_dyn(14) + dens_curr
         upmf_dyn(16) = 0.9 * upmf_dyn(14) + dens_curr
         upmf_dyn(14) = upmf_dyn(14)       + dens_curr
      end if
   end select
   
   if (comp_down) dnmf_dyn(:) = edt * upmf_dyn(:)
   
   return
end subroutine grell_dyncontrol_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes takes the ensemble averages, integrate them and then compute !
! should be called inside the other loops, so the ensemble variables that should have      !
! the permutation of all ensemble dimensions is fully filled.                              !
!------------------------------------------------------------------------------------------!
subroutine grell_feedback(comp_down,mgmzp,maxens_cap,maxens_eff,maxens_lsf,maxens_dyn      &
                         ,inv_ensdim,max_heat,ktop,edt_eff,dellathil_eff,dellaqtot_eff     &
                         ,pw_eff,dnmf_ens,upmf_ens,ierr,upmf,dnmf,edt,outthil,outqtot      &
                         ,precip)

   use rconstants,   only: day_sec
   implicit none
   logical               , intent(in)    :: comp_down  ! I am computing downdrafts
   integer               , intent(in)    :: mgmzp      ! # of levels
   integer               , intent(in)    :: maxens_cap ! # of static controls
   integer               , intent(in)    :: maxens_eff ! # of prec. efficiency members
   integer               , intent(in)    :: maxens_lsf ! # of mb members
   integer               , intent(in)    :: maxens_dyn ! # of dynamic control members
   real                  , intent(in)    :: inv_ensdim ! 1 / # of members 
   real                  , intent(in)    :: max_heat   ! Maximum heat allowed
   integer               , intent(in)    :: ktop       ! Cloud top

   !---------------------------------------------------------------------------------------!
   ! Ensemble variables                                                                    !
   !---------------------------------------------------------------------------------------!
   !----- Fraction between reference downdraft and updraft --------------------------------!
   real, dimension(maxens_eff,maxens_cap)      , intent(in) ::  &
                                            edt_eff        ! ! Grell's epsilon
   !----- Changes of the thermodynamic properties per unit of mass. -----------------------!
   real, dimension(mgmzp,maxens_eff,maxens_cap), intent(in) ::  & 
                                            dellathil_eff  & ! Ice-liquid potential temp.
                                           ,dellaqtot_eff  & ! Total mixing ratio
                                           ,pw_eff         ! ! Fall-out water
   !----- Reference mass fluxes, in [kg/m²/s] ---------------------------------------------!
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(in) :: &
                                            dnmf_ens    & ! Downdraft reference
                                           ,upmf_ens    ! ! Updraft reference
   !---------------------------------------------------------------------------------------!
   
   integer               , intent(inout) :: ierr        ! Error flag
   real                  , intent(out)   :: upmf        ! Ref. upward mass flx (Grell's mb)
   real                  , intent(out)   :: dnmf        ! Ref. dnward mass flx (Grell's m0)
   real                  , intent(out)   :: edt         ! m0/mb, Grell's epsilon
   real, dimension(mgmzp), intent(out)   :: outthil     ! Change in temperature profile
   real, dimension(mgmzp), intent(out)   :: outqtot     ! Change in total mixing ratio
   real                  , intent(out)   :: precip      ! Precipitation rate

   integer                               :: k,l               ! Counter
   integer                               :: kmin              ! Minimum location index
   integer                               :: kmax              ! Maximum location index
   integer                               :: iedt              ! efficiency counter
   integer                               :: icap              ! Cap_max counter
   real                                  :: rescale           ! Rescaling factor
   real                                  :: inv_maxens_effcap ! 1/( maxens_eff*maxens_cap)

   inv_maxens_effcap = 1./(maxens_eff*maxens_cap)

   !---------------------------------------------------------------------------------------!
   !    Initialize all output variables. I may need to break in case upmf is a bad one.    !
   !---------------------------------------------------------------------------------------!
   dnmf    = 0.
   edt     = 0.
   outthil = 0.
   outqtot = 0.
   precip  = 0.
      
   !---------------------------------------------------------------------------------------!
   !    Here I am also doing something different from the original code. upmf (the mb term !
   ! using Grell's notation) is simply the average accross all members, not only the       !
   ! positive ones. My claim is that if most terms are zero, then little convection        !
   ! should happen. If it turns out to be too little rain, I switch this back to the old   !
   ! style.                                                                                ! 
   !---------------------------------------------------------------------------------------!
   upmf=sum(upmf_ens) * inv_ensdim

   !---------------------------------------------------------------------------------------!
   !  If the reference upward mass flux is zero, that means that there is no cloud...      !
   !---------------------------------------------------------------------------------------!
   if (upmf == 0.) then
      ierr = 10
      return
   end if

   !---------------------------------------------------------------------------------------!
   !   Average the temperature tendency among the precipitation efficiency ensemble. If it !
   ! is heating/cooling too much, rescale the reference upward mass flux. Since max_heat   !
   ! is given in K/day, I will compute outt in K/day just to test the value and check      !
   ! whether it is outside the allowed range or not. After I'm done, I will rescale it     !
   ! back to K/s.                                                                          !
   !---------------------------------------------------------------------------------------!
   do k=1,ktop
      outthil(k) = upmf * sum(dellathil_eff(k,:,:)) * inv_maxens_effcap * day_sec
   end do
   !----- Get minimum and maximum outt, and where they happen -----------------------------!
   kmin = minloc(outthil,dim=1)
   kmax = maxloc(outthil,dim=1)
   
   !----- If excessive heat happens, scale down -------------------------------------------!
   if (kmax > 2 .and. outthil(kmax) > 2.0 * max_heat) then
      rescale = 2.0*max_heat / outthil(kmax)
      upmf = upmf * rescale
      do k=1,ktop
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- If excessive cooling happens, scale down ----------------------------------------!
   if (outthil(kmin)  < - max_heat) then
      rescale = - max_heat/ outthil(kmin)
      upmf = upmf * rescale
      do k=1,ktop
         outthil(k) = outthil(k) * rescale
      end do
   end if
   !----- Heating close to the surface needs to be smaller, being strict there ------------!
   do k=1,2
      if (outthil(k) > 0.5 * max_heat) then
         rescale = 0.5 * max_heat / outthil(k)
         upmf = upmf * rescale
         do l=1,ktop
            outthil(l) = outthil(l) * rescale
         end do
      end if
   end do
   !----- Converting outt to K/s ----------------------------------------------------------!
   do k=1,ktop
      outthil(k) = outthil(k) / day_sec
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    With the mass flux and heating on the track, compute other sources/sinks           !
   !---------------------------------------------------------------------------------------!
   do k=1,ktop
      outqtot(k)  = upmf * sum(dellaqtot_eff(k,:,:) ) * inv_maxens_effcap
   end do

   !---------------------------------------------------------------------------------------!
   !   Computing precipitation. It should never be negative, so making sure that this      !
   ! never happens. I will skip this in case this cloud is too shallow.                    !
   !---------------------------------------------------------------------------------------!
   precip = 0.
   if (comp_down) then
      do icap=1,maxens_cap
         do iedt=1,maxens_eff
            precip        = precip + max(0.,upmf*sum(pw_eff(1:ktop,iedt,icap)))
         end do
      end do
      precip = precip * inv_maxens_effcap
   end if
   
   !---------------------------------------------------------------------------------------!
   !    Compute downdraft mass flux, and redefining epsilon.                               !
   !---------------------------------------------------------------------------------------!
   if (comp_down .and. upmf > 0) then
      dnmf = sum(dnmf_ens) * inv_ensdim
      edt  = dnmf/upmf
   end if

   return
end subroutine grell_feedback
!==========================================================================================!
!==========================================================================================!
