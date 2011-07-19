!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the variables that depend on heat, water, and (eventually)  !
! carbon exchanges with the canopy air space and vegetation.  This solves only the         !
! prognostic variables, the diagnostic ones will be updated in a separate sub-routine.     !
!------------------------------------------------------------------------------------------!
subroutine leaf3_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass         &
                       ,ustar,tstar,rstar,cstar,zeta,ribulk,soil_rough,veg_rough           &
                       ,patch_rough,veg_height,veg_displace,veg_lai,veg_tai,veg_water      &
                       ,veg_hcap,veg_energy,leaf_class,veg_fracarea,stom_condct,can_prss   &
                       ,can_rvap,can_co2,sensible_gc,sensible_vc,evap_gc,evap_vc           &
                       ,transp,gpp,plresp,resphet,ground_rsat,ground_rvap,ground_temp      &
                       ,ground_fliq,available_water,rshort)
   use leaf_coms
   use rconstants
   use therm_lib , only : eslif          & ! function
                        , rslif          & ! function
                        , thetaeiv       & ! function
                        , thrhsh2temp    & ! function
                        , tslif          & ! function
                        , qwtk           ! ! subroutine
   use catt_start, only : CATT
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer             , intent(in)    :: mzg
   integer             , intent(in)    :: mzs
   integer             , intent(in)    :: ksn
   real, dimension(mzg), intent(in)    :: soil_energy
   real, dimension(mzg), intent(in)    :: soil_water
   real, dimension(mzg), intent(in)    :: soil_text
   real, dimension(mzs), intent(in)    :: sfcwater_mass
   real                , intent(in)    :: ustar
   real                , intent(in)    :: tstar
   real                , intent(in)    :: rstar
   real                , intent(in)    :: cstar
   real                , intent(in)    :: zeta
   real                , intent(in)    :: ribulk
   real                , intent(in)    :: soil_rough
   real                , intent(in)    :: veg_rough
   real                , intent(in)    :: patch_rough
   real                , intent(in)    :: veg_height
   real                , intent(in)    :: veg_displace
   real                , intent(in)    :: veg_lai
   real                , intent(in)    :: veg_tai
   real                , intent(in)    :: leaf_class
   real                , intent(in)    :: veg_fracarea
   real                , intent(in)    :: ground_rsat
   real                , intent(in)    :: ground_rvap
   real                , intent(in)    :: ground_temp
   real                , intent(in)    :: ground_fliq
   real                , intent(in)    :: available_water
   real                , intent(in)    :: rshort
   real                , intent(inout) :: veg_water
   real                , intent(inout) :: veg_hcap
   real                , intent(inout) :: veg_energy
   real                , intent(inout) :: stom_condct
   real                , intent(inout) :: can_prss
   real                , intent(inout) :: can_rvap
   real                , intent(inout) :: can_co2
   real                , intent(inout) :: sensible_gc
   real                , intent(inout) :: sensible_vc
   real                , intent(inout) :: evap_gc
   real                , intent(inout) :: evap_vc
   real                , intent(inout) :: transp
   real                , intent(inout) :: gpp
   real                , intent(inout) :: plresp
   real                , intent(inout) :: resphet
   !----- Local arguments. ----------------------------------------------------------------!
   integer                             :: k
   integer                             :: kk
   integer                             :: nveg
   integer                             :: nsoil
   real                                :: aux
   real                                :: dsm
   real                                :: fac
   real                                :: factv
   real                                :: fthi
   real                                :: ftlo
   real                                :: frad
   real                                :: fswp
   real                                :: fvpd
   real                                :: qwtot
   real                                :: lsfc_rvap
   real                                :: lint_rvap
   real                                :: lsfc_pvap
   real                                :: lint_pvap
   real                                :: sigmaw 
   real                                :: slai
   real                                :: stai
   real                                :: slpotv
   real                                :: swp
   real                                :: vpd
   real                                :: wtemp
   real                                :: wtroot
   real                                :: x
   real                                :: zognd
   real                                :: zoveg
   real                                :: zdisp
   real                                :: zveg
   real                                :: wflx
   real                                :: dewgndflx
   real                                :: qdewgndflx
   real                                :: ddewgndflx
   real                                :: transp_test
   real                                :: transp_wilt
   real                                :: gsw_wilt
   real                                :: gsw_inf
   real                                :: wshed
   real                                :: qwshed
   real                                :: dwshed
   real                                :: old_veg_water
   real                                :: old_veg_energy
   !----- Local constants. ----------------------------------------------------------------!
   character(len=9)      , parameter   :: fmti='(a,1x,i6)'
   character(len=13)     , parameter   :: fmtf='(a,1x,es12.5)'
   character(len=3)      , parameter   :: fmtc='(a)'
   character(len=9)      , parameter   :: fmtl='(a,1x,l1)'
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: leaf3_reduced_wind
   !---------------------------------------------------------------------------------------!



   !----- Find the time step auxiliary variables. -----------------------------------------!
   dtllowcc = dtll / (can_depth * can_rhos)
   dtllohcc = dtll / (can_depth * can_rhos * cp * can_temp)
   dtlloccc = dtllowcc * mmdry
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the atmosphere -> canopy fluxes.                                             !
   !---------------------------------------------------------------------------------------!
   rho_ustar  = can_rhos  * ustar
   eflxac     = rho_ustar * estar * cp * can_temp ! Enthalpy exchange
   hflxac     = rho_ustar * tstar * can_exner     ! Sensible heat exchange
   wflxac     = rho_ustar * rstar                 ! Water vapour exchange
   cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute ground-canopy resistance rgnd.  Assume zognd not affected by snow.        !
   ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining the values          !
   ! (zognd,0) when vegetation is fully buried in snow.                                    !
   !---------------------------------------------------------------------------------------!
   zognd = soil_rough
   zoveg = veg_rough    * (1. - snowfac) + zognd * snowfac
   zdisp = veg_displace * (1. - snowfac)
   zveg  = veg_height   * (1. - snowfac)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the ground conductance due to vegetation and find the total conductance      !
   ! (the inverse of the sum due to bare ground plus vegetation).                          !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then

      !------------------------------------------------------------------------------------!
      !    If vegetation is sufficiently abundant and not covered by snow, compute  heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      factv = log((geoht-zdisp) / zoveg) / (vonk * vonk * atm_vels)
      aux   = exp(exar * (1. - (zdisp + zoveg) / zveg))
      ggveg = (exar * (zveg - zdisp)) / (factv * zveg  * (exp(exar) - aux))
   else
      ggveg = huge_num
   end if
   ggnet    = ggfact * ggbare * ggveg / (ggbare + ggveg)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Split the total precipitation into intercepted water and through-fall, according  !
   ! to the fraction of open canopy.                                                       !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      intercepted_tot  = pcpgl  * veg_fracarea * (1. - snowfac)
      qintercepted_tot = qpcpgl * veg_fracarea * (1. - snowfac)
      dintercepted_tot = dpcpgl * veg_fracarea * (1. - snowfac)
   else
      intercepted_tot  = 0.0
      qintercepted_tot = 0.0
      dintercepted_tot = 0.0
   end if
   !----- The throughfall is what is not intercepted.  Shedding is zero at this point. ----!
   throughfall_tot  = pcpgl  - intercepted_tot
   qthroughfall_tot = qpcpgl - qintercepted_tot
   dthroughfall_tot = dpcpgl - dintercepted_tot
   wshed_tot        = 0.0
   qwshed_tot       = 0.0
   dwshed_tot       = 0.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between the ground and the canopy air   !
   ! space.  The ground may be either the top soil layer or the temporary surface water    !
   ! snow surface.  wflx is the possible vapour flux, and can be either positive or        !
   ! negative.  After we compute it, we check its sign, because the value is good only     !
   ! when it is negative or when the surface is water/snow.  The ground calculation must   !
   ! be done differently.  Later, two variables will define the flux between ground and    !
   ! canopy air space:                                                                     !
   ! - wflxgc [kg/m2/s] is going to be the evaporation flux from ground to canopy.         !
   ! - dewgnd [kg/m2/s] is going to be the dew/frost flux from canopy to ground.           !
   !     Both are defined as positive quantities.  Sensible heat is defined by only one    !
   ! variable, hflxgc [J/m2/s], which can be either positive or negative.                  !
   !---------------------------------------------------------------------------------------!
   hflxgc = ggnet * can_rhos * cp * (ground_temp - can_temp)
   wflx   = ggnet * can_rhos      * (ground_rsat - can_rvap)
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !   Here we will decide how to compute the evaporation and condensation fluxes based on !
   ! the sign of wflx, the number of temporary surface water/snow layers, and whether the  !
   ! canopy air is saturation and whether the user is concerned about super-saturation.    !
   !---------------------------------------------------------------------------------------!
   if (wflx <= 0.0) then
      !------------------------------------------------------------------------------------!
      !     Flux is negative, which means that dew/frost is going to form.  We must also   !
      ! decide whether dew or frost will form, and this is based on the current partition  !
      ! of the surface water phase.  The depth gain uses frost density rather than ice     !
      ! density based on MCD suggestion on 11/16/2009.                                     !
      !------------------------------------------------------------------------------------!
      dewgndflx  = min(-wflx,(can_rvap - toodry) / dtllowcc)
      qdewgndflx = dewgndflx * (alvi - ground_fliq * alli)
      ddewgndflx = dewgndflx * (ground_fliq * wdnsi + (1.0-ground_fliq) * fdnsi)
      !----- Set evaporation fluxes to zero. ----------------------------------------------!
      wflxgc     = 0.0
      qwflxgc    = 0.0

   elseif (can_rhv >= 1.0 .and. (.not. supersat_ok)) then
      !------------------------------------------------------------------------------------!
      !     In principle evaporation could happen, but the canopy air space is saturated.  !
      ! The user doesn't want too much super-saturation, so we will not let any            !
      ! evaporation to occur.                                                              !
      !------------------------------------------------------------------------------------!
      dewgndflx  = 0.0
      qdewgndflx = 0.0
      ddewgndflx = 0.0
      wflxgc     = 0.0
      qwflxgc    = 0.0

   elseif (ksn > 0) then
      !------------------------------------------------------------------------------------!
      !     Evaporation will happen, and there is a temporary surface water/snow layer     !
      ! above the soil, so wflx is still good, we simply use it, and make the dew/frost    !
      ! fluxes to be zero.                                                                 !
      !------------------------------------------------------------------------------------!
      dewgndflx  = 0.0
      qdewgndflx = 0.0
      ddewgndflx = 0.0
      wflxgc     = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
      qwflxgc    = wflx * (alvi - ground_fliq * alli)

   else
      !------------------------------------------------------------------------------------!
      !    Evaporation will happen, but water will come from the top most layer.  Wflx     !
      ! cannot be used because the ground specific humidity is not the saturation specific !
      ! humidity at the soil temperature only, it depends on the canopy specific humidity  !
      ! itself and the soil moisture.  ground_rvap cannot be less than can_rvap, so we no  !
      ! longer need to check for sign.                                                     !
      !------------------------------------------------------------------------------------!
      wflxgc  = can_rhos * ggnet * ggsoil * (ground_rvap - can_rvap) / (ggnet + ggsoil)
      !----- Adjust the flux according to the surface fraction (no phase bias). -----------!
      qwflxgc = wflxgc * ( alvi - ground_fliq * alli)
      !----- Set condensation fluxes to zero. ---------------------------------------------!
      dewgndflx  = 0.0
      qdewgndflx = 0.0
      ddewgndflx = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !----- Find the total input of water due to dew/frost formation. -----------------------!
   dewgnd_tot  = dewgndflx  * dtll
   qdewgnd_tot = qdewgndflx * dtll
   ddewgnd_tot = ddewgndflx * dtll
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If vegetation is sufficiently abundant and not covered by snow, compute heat and  !
   ! moisture fluxes from vegetation to canopy, and flux resistance from soil or snow to   !
   ! canopy.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      !----- Save the vegetation class in an alias. ---------------------------------------!
      nveg = nint(leaf_class)
   
      !------------------------------------------------------------------------------------!
      !     Find the wind speed at the top of the canopy, using the similarity theory.     !
      !------------------------------------------------------------------------------------!
      veg_wind = leaf3_reduced_wind(ustar,zeta,ribulk,geoht,veg_displace,veg_height        &
                                  ,patch_rough)


      !------------------------------------------------------------------------------------!
      !    Find the aerodynamic conductances for heat and water at the leaf boundary       !
      ! layer.                                                                             !
      !------------------------------------------------------------------------------------!
      call leaf3_aerodynamic_conductances(nveg,veg_wind,veg_temp,can_temp,can_rvap,can_rhos)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the "perceived" leaf area and tree area indices, dependending on snow     !
      ! cover.                                                                             !
      !------------------------------------------------------------------------------------!
      stai = veg_tai * (1. - snowfac)
      slai = veg_lai * (1. - snowfac)
      !------------------------------------------------------------------------------------!


      !----- Soil water potential factor for stomatal control. ----------------------------!
      swp = -200.0
      do k = kroot(nveg),mzg
         nsoil  = nint(soil_text(k))
         slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
         if (slpotv > swp) swp = slpotv
      end do
      swp = swp * grav * wdns
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the intercellular vapour pressure and mixing ratio; they are both assumed !
      ! to be the saturation value at the leaf temperature.                                !
      !------------------------------------------------------------------------------------!
      lint_pvap = eslif(veg_temp)
      lint_rvap = rslif(can_prss,veg_temp)
      !------------------------------------------------------------------------------------!


      !----- Compute mixing ratio at leaf surface using previous gsw ----------------------!
      gsw       = stom_condct
      lsfc_rvap = (gbw * can_rvap + gsw * lint_rvap) / (gbw + gsw)
      lsfc_pvap = lsfc_rvap * can_prss / (ep + lsfc_rvap)
      vpd       = max(lint_pvap - lsfc_pvap,0.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Use a combination of 5 environmental factors that may control stomatal        !
      ! conductance to estimate the target value for this new time step.                   !
      !------------------------------------------------------------------------------------!
      ftlo    = 1. + exp(-stlo * (veg_temp - btlo))
      fthi    = 1. + exp(-sthi * (veg_temp - bthi))
      frad    = 1. + exp(-srad * (rshort   - brad))
      fswp    = 1. + exp(-sswp * (swp      - bswp))
      fvpd    = 1. + exp(-svpd * (vpd      - bvpd))
      if (slai > tai_min) then
         gsw_inf = can_rhos * gsw_max(nveg) / (ftlo * fthi * frad * fvpd * fswp) / slai
      else
         gsw_inf = 0.0
      end if
      !------------------------------------------------------------------------------------!

      !----- 15-minute response time for stomatal conductance (must have dtll <= 900.). ---!
      gsw     = gsw + dtll * (gsw_inf - gsw) / max(dtll,900.)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Transpiration can only happen if there are leaves and water.                    !
      !------------------------------------------------------------------------------------!
      if (available_water > 0. .and. gsw > 0.0) then

         !----- Find the maximum transpiration allowed. -----------------------------------!
         transp_wilt  = min(available_water / dtll, transp_max)

         !------ Find the stomatal conductance. -------------------------------------------!
         if (transp_wilt >= gbw * (lint_rvap - can_rvap)) then
            !------------------------------------------------------------------------------!
            !   In this case, wilting would only happen if gsw could be negative, which is !
            ! physically unrealistic.  No need to worry, use the actual conductance.       !
            !------------------------------------------------------------------------------!
            transp_test  = stom_side(nveg) * gbw * gsw * slai * (lint_rvap - can_rvap)     &
                         / (gbw + gsw)
         else
            !------------------------------------------------------------------------------!
            !    Limit maximum transpiration to be <= transp_max and less than the maximum !
            ! water available by decreasing stomatal conductance if necessary.             !
            !------------------------------------------------------------------------------!
            gsw_wilt     = gbw * transp_wilt / ((lint_rvap - can_rvap) * gbw - transp_wilt)
            gsw          = min(gsw,gsw_wilt)
            transp_test  = stom_side(nveg) * gbw * gsw * slai * (lint_rvap - can_rvap)     &
                         / (gbw + gsw)
         end if
      else
         transp_test = 0.
      end if
      stom_condct = gsw
      !------------------------------------------------------------------------------------!


      !----- Sensible heat flux from leaf to canopy. --------------------------------------!
      hflxvc = 2. * stai * gbh * (veg_temp - can_temp)

      !------------------------------------------------------------------------------------!
      !     Check whether evaporation or condensation should happen.                       !
      !------------------------------------------------------------------------------------!
      if (lint_rvap >= can_rvap) then
         !---------------------------------------------------------------------------------!
         !     Saturation at leaf temperature is greater than canopy air space mixing      !
         ! ratio.  If some water is sitting on the leaf surface, let it evaporate.         !
         ! Evaporation should happen from only one side of the leaf.                       !
         !---------------------------------------------------------------------------------!
         !----- Compute the fraction of leaves covered in water. --------------------------!
         sigmaw = min(1.0, (veg_water / (.2 * stai)) ** twothirds)
         !----- Find the evaporation rate. ------------------------------------------------!
         wflxvc = min(stai * gbw * (lint_rvap - can_rvap) * sigmaw,veg_water/dtll)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Also, if there is available soil water and let it transpire.                !
         !---------------------------------------------------------------------------------!
         transp_loc = max(0., transp_test)
      else
         !---------------------------------------------------------------------------------!
         !     Flow is towards the leaf, dew/frost formation on one side of the leaf.  We  !
         ! don't let the dew formation completely eliminate water from the canopy air      !
         ! space, and impose 0. transpiration.                                             !
         !---------------------------------------------------------------------------------!
         wflxvc     = max( stai * gbw * (lint_rvap - can_rvap)                             &
                         , min(0.,(toodry-can_rvap)/ dtllowcc) )
         transp_loc = 0.
      end if
      transp_tot = transp_tot + transp_loc
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Update vegetation moisture from precipitation, vapor flux with canopy,         !
      ! and shedding of excess moisture.                                                   !
      !------------------------------------------------------------------------------------!
      old_veg_water = veg_water
      wtemp         = veg_water - wflxvc * dtll + intercepted_tot
      !----- If this result will lead to negative solution, then reduce evaporation. ------!
      if (wtemp < 0.005 * leaf_maxwhc * stai) then
         wflxvc     = (veg_water + intercepted_tot) / dtll
         veg_water  = 0.
         wshed_tot  = 0.
         qwshed_tot = 0.
         dwshed_tot = 0.
      elseif (wtemp > leaf_maxwhc * stai) then
         !----- Shed water in excess of the leaf holding water capacity. ------------------!
         wshed_tot  = wtemp - leaf_maxwhc * stai
         qwshed_tot = wshed_tot * (     veg_fliq  * cliq * (veg_temp - tsupercool)         &
                                  + (1.-veg_fliq) * cice *  veg_temp )
         dwshed_tot = wshed_tot * (veg_fliq * wdnsi + (1.0 - veg_fliq) * fdnsi)
         veg_water  = leaf_maxwhc * stai
      else
         wshed_tot  = 0.
         qwshed_tot = 0.
         dwshed_tot = 0.
         veg_water  = wtemp
      end if
      !------------------------------------------------------------------------------------!



      !------ Find the associated latent heat flux from vegetation to canopy. -------------!
      qwflxvc     = wflxvc     * (alvi - alli * veg_fliq)
      qtransp_loc = transp_loc *  alvl !----- Liquid phase only in transpiration. ---------!
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Update vegetation internal energy from radiative fluxes and sensible and      !
      ! latent heat transfer with canopy air.                                              !
      !------------------------------------------------------------------------------------!
      old_veg_energy = veg_energy
      veg_energy     = veg_energy                                                          &
                     + dtll * ( rshort_v + rlonga_v + rlonggs_v - rlongv_gs - rlongv_a     &
                              - hflxvc   - qwflxvc   - qtransp_loc)                        &
                     + qintercepted_tot - qwshed_tot




      !------------------------------------------------------------------------------------!
      !      Update enthalpy, CO2, and canopy mixing ratio.                                !
      !------------------------------------------------------------------------------------!
      can_lntheta      = can_lntheta                                                       &
                       + dtllohcc * ( hflxgc + hflxvc + hflxac)
      can_rvap         = can_rvap                                                          &
                       + dtllowcc * (wflxgc - dewgndflx + wflxvc + transp_loc + wflxac)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc   = 0.
      cflxvc   = 0.
      cflxcv   = 0.
      can_co2  = can_co2  + dtllowcc * mmdry * (cflxgc + cflxvc + cflxcv + cflxac)
      !------------------------------------------------------------------------------------!



      !----- Update the fluxes. -----------------------------------------------------------!
      sensible_gc = sensible_gc +  hflxgc                * dtll_factor
      sensible_vc = sensible_vc +  hflxvc                * dtll_factor
      evap_gc     = evap_gc     + (qwflxgc - qdewgndflx) * dtll_factor
      evap_vc     = evap_vc     +  qwflxvc               * dtll_factor
      transp      = transp      +  qtransp_loc           * dtll_factor

   else
      !------------------------------------------------------------------------------------!
      !    No vegetation, or vegetation buried in snow...                                  !
      !------------------------------------------------------------------------------------!

      !----- Update the canopy prognostic variables. --------------------------------------!
      can_lntheta  = can_lntheta  + dtllohcc * (hflxgc + qwflxgc   + hflxac)
      can_rvap     = can_rvap     + dtllowcc * (wflxgc - dewgndflx + wflxac)

      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc       = 0.
      can_co2      = can_co2      + dtlloccc * (cflxgc + cflxac)
      !------------------------------------------------------------------------------------!

      !----- Update the fluxes. -----------------------------------------------------------!
      sensible_gc = sensible_gc +  hflxgc                * dtll_factor
      evap_gc     = evap_gc     + (qwflxgc - qdewgndflx) * dtll_factor

   end if

   return
end subroutine leaf3_canopy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of canopy air space variables for a given     !
! patch, before the time step iteration loop.                                              !
!------------------------------------------------------------------------------------------!
subroutine leaf3_can_diag(ip,can_theta,can_theiv,can_rvap,leaf_class,can_prss,initial)
   use leaf_coms , only : atm_prss      & ! intent(in)
                        , atm_theta     & ! intent(in)
                        , atm_shv       & ! intent(in)
                        , geoht         & ! intent(in)
                        , veg_ht        & ! intent(in)
                        , can_shv       & ! intent(out)
                        , can_rsat      & ! intent(out)
                        , can_rhv       & ! intent(out)
                        , can_depth     & ! intent(inout)
                        , can_depth_min & ! intent(in)
                        , can_temp      & ! intent(out)
                        , can_exner     & ! intent(inout)
                        , can_lntheta   & ! intent(inout)
                        , can_rhos      ! ! intent(inout)
   use rconstants, only : cp            & ! intent(in)
                        , cpi           & ! intent(in)
                        , ep            & ! intent(in)
                        , p00i          & ! intent(in)
                        , rocp          ! ! intent(in)
   use therm_lib , only : reducedpress  & ! function
                        , rslif         & ! function
                        , idealdenssh   & ! function
                        , thetaeiv      ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)      :: ip
   real   , intent(inout)   :: can_theta
   real   , intent(out)     :: can_theiv
   real   , intent(in)      :: can_rvap
   real   , intent(in)      :: leaf_class
   real   , intent(inout)   :: can_prss
   logical, intent(in)      :: initial
   !----- Local variables. ----------------------------------------------------------------!
   integer                  :: nveg
   !---------------------------------------------------------------------------------------!




   !----- Define the canopy depth. --------------------------------------------------------!
   if (initial) then
      nveg      = nint(leaf_class)
      if (nveg == 0 .or. ip == 1) then
         can_depth = can_depth_min
      else
         can_depth = max(can_depth_min,veg_ht(nveg))
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Set initial value of specific humidity ------------------------------------------!
   can_shv     = can_rvap / (can_rvap + 1.)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     If this is the initial step, we initialise the canopy air "dry enthropy" (log of  !
   ! potential temperature).  If not, we update the actual potential temperature from the  !
   ! logarithm.                                                                            !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_lntheta = log(can_theta)
   else
      can_theta   = exp(can_lntheta)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update canopy air pressure and the Exner function here.  Canopy air pressure is    !
   ! assumed to remain constant during one LEAF full timestep, which means that heat       !
   ! equals to change in enthalpy.  Between time steps, atmospheric pressure changes so    !
   ! enthalpy is no longer conserved.  Potential temperature and equivalent potential      !
   ! temperature, on the other hand, are conserved if no heat happens, which is the case   !
   ! between successive calls.  Therefore, we track theta and theta_eiv.                   !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_prss  = reducedpress(atm_prss,atm_theta,atm_shv,geoht,can_theta,can_shv          &
                              ,can_depth)
      can_exner = cp  * (p00i * can_prss) ** rocp
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Find the canopy air temperature.
   !---------------------------------------------------------------------------------------!
   can_temp = cpi * can_theta * can_exner
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Even though this leads to small violations of the ideal gas law, we are also      !
   ! going to keep the canopy air density constant over one LEAF time step.  This is to    !
   ! avoid that non-linearities sprout and cause problems to the code.                     !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_rhos = idealdenssh (can_prss,can_temp,can_shv)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the canopy air saturation mixing ratio and the relative humidity.  These are !
   ! used to decide whether the canopy air space can hold more water or not.               !
   !---------------------------------------------------------------------------------------!
   can_rsat = rslif(can_prss,can_temp)
   can_rhv  = can_rvap * (ep + can_rsat) / (can_rsat * (ep + can_rvap))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the canopy air space ice-vapour equivalent potential temperature.  This is   !
   ! currently a diagnostic variable only, but it should become the main variable if we    !
   ! ever switch to foggy canopy air space.                                                !
   !---------------------------------------------------------------------------------------!
   can_theiv = thetaeiv(can_theta,can_prss,can_temp,can_rvap,can_rvap,-84)
   !---------------------------------------------------------------------------------------!



   return
end subroutine leaf3_can_diag
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of leaf-level variables for a given patch,    !
! before the time step iteration loop.                                                     !
!------------------------------------------------------------------------------------------!
subroutine leaf3_veg_diag(veg_energy,veg_water,veg_hcap)
   use leaf_coms , only : can_temp         & ! intent(in)
                        , min_patch_area   & ! intent(in)
                        , veg_temp         & ! intent(out)
                        , veg_fliq         & ! intent(out)
                        , resolvable       ! ! intent(in)
   use rconstants, only : t3ple            ! ! intent(in)
   use therm_lib , only : qwtk             ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(inout) :: veg_energy
   real   , intent(inout) :: veg_water
   real   , intent(in)    :: veg_hcap
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find leaf temperature and liquid fraction of water, but only if we are solving    !
   ! leaves.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)
   else
      veg_temp = can_temp
      if (veg_temp > t3ple) then
         veg_fliq = 1.0
      elseif (veg_temp < t3ple) then
         veg_fliq = 0.0
      else
         veg_fliq = 0.5
      end if
      veg_energy = veg_hcap * veg_temp
      veg_water  = 0.0
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_veg_diag
!==========================================================================================!
!==========================================================================================!
