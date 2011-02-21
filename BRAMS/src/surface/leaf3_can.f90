!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the variables that depend on heat, water, and (eventually)  !
! carbon exchanges with the canopy air space and vegetation.                               !
!------------------------------------------------------------------------------------------!
subroutine leaf_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass          &
                      ,ustar,tstar,rstar,cstar,zeta,ribulk,soil_rough,veg_rough            &
                      ,patch_rough,veg_height,veg_lai,veg_tai,veg_water,veg_hcap           &
                      ,veg_energy,leaf_class,veg_fracarea,stom_condct,can_prss,can_theiv   &
                      ,can_theta,can_rvap,can_co2,sensible_gc,sensible_vc,evap_gc,evap_vc  &
                      ,transp,gpp,plresp,resphet,ground_rsat,ground_rvap,ground_temp       &
                      ,ground_fliq,available_water,rshort,i,j,ip)
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
   real                , intent(inout) :: can_theiv
   real                , intent(inout) :: can_theta
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
   integer                             :: i
   integer                             :: j
   integer                             :: ip
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
   real                                :: transp_test
   real                                :: transp_wilt
   real                                :: gsw_wilt
   real                                :: gsw_inf
   real                                :: wshed_loc
   real                                :: qwshed_loc
   real                                :: old_veg_water
   real                                :: old_veg_energy
   real                                :: old_can_rvap
   real                                :: old_can_theiv
   real                                :: old_can_shv
   real                                :: old_can_theta
   real                                :: old_can_temp
   real                                :: old_can_prss
   !----- Local constants. ----------------------------------------------------------------!
   character(len=9)      , parameter   :: fmti='(a,1x,i6)'
   character(len=13)     , parameter   :: fmtf='(a,1x,es12.5)'
   character(len=3)      , parameter   :: fmtc='(a)'
   character(len=9)      , parameter   :: fmtl='(a,1x,l1)'
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: leaf_reduced_wind
   !---------------------------------------------------------------------------------------!



   !----- Find the time step auxiliary variables. -----------------------------------------!
   dtllowcc = dtll / (can_depth * can_rhos)
   dtllohcc = dtll / (can_depth * can_rhos * cp * can_temp)
   dtlloccc = dtllowcc * mmdry
   !---------------------------------------------------------------------------------------!


   !----- Save the previous values for debugging. -----------------------------------------!
   old_can_theiv    = can_theiv
   old_can_rvap     = can_rvap
   old_can_shv      = can_shv
   old_can_theta    = can_theta
   old_can_temp     = can_temp
   old_can_prss     = can_prss



   !---------------------------------------------------------------------------------------!
   !     Compute ground-canopy resistance rgnd.  Assume zognd not affected by snow.        !
   ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining the values          !
   ! (zognd,0) when vegetation is fully buried in snow.                                    !
   !---------------------------------------------------------------------------------------!

   zognd = soil_rough
   zoveg = veg_rough * (1.-snowfac) + zognd * snowfac
   zdisp = veg_height * (1.-snowfac)
   zveg = zdisp / 0.63

   if (veg_tai >= .1 .and. snowfac < .9) then

      !------------------------------------------------------------------------------------!
      !    If vegetation is sufficiently abundant and not covered by snow, compute  heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      factv       = log(geoht / zoveg) / (vonk * vonk * atm_vels)
      aux         = exp(exar * (1. - (zdisp + zoveg) / zveg))
      ggveg       = (exar * (zveg - zdisp)) / (factv * zveg  * (exp(exar) - aux))
      !------------------------------------------------------------------------------------!
      !    The net ground conductance is computed based on the weighted average of the     !
      ! bare ground and vegetated ground _resistances_.                                    !
      !------------------------------------------------------------------------------------!
      wshed_tot   = 0.
      qwshed_tot  = 0.
      transp_tot  = 0.
   else
      !------------------------------------------------------------------------------------!
      !     If the TAI is very small or if snow mostly covers the vegetation, bypass       !
      ! vegetation computations.  Set heat and moisture flux resistance rgnd between the   !
      ! "canopy" and snow or soil surface to its bare soil value.  Set shed precipitation  !
      ! heat and moisture to unintercepted values.                                         !
      !------------------------------------------------------------------------------------!
      wshed_tot  = pcpgl
      qwshed_tot = qpcpgl
      transp_tot = 0.
      ggveg      = 0.0
   end if
   ggnet = ggfact * ggbare

   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between top soil or snow surface and    !
   ! canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow evaporation  !
   ! and dewgnd is the mass of dew that forms on the snow/soil surface this timestep; both !
   ! are defined as always positive or zero.                                               !
   !---------------------------------------------------------------------------------------!
   hflxgc     = cp * can_rhos * ggnet * (ground_temp - can_temp)
   wflx       =      can_rhos * ggnet * (ground_rsat - can_rvap)
   
   if (wflx >= 0.) then
      dewgndflx = 0.
   else
      dewgndflx = min(-wflx,(can_rvap - toodry) / dtllowcc)
   end if
   dewgnd_tot = dewgndflx * dtll

   if (ksn == 0) then
      wflxgc = max(0.,can_rhos * ggnet * (ground_rvap - can_rvap))
   else
      wflxgc = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
   end if

   qdewgnd_tot = (alvi - alli * ground_fliq) * dewgnd_tot
   qdewgndflx  = (alvi - alli * ground_fliq) * dewgndflx
   qwflxgc     = (alvi - alli * ground_fliq) * wflxgc

   !---------------------------------------------------------------------------------------!
   !     If vegetation is sufficiently abundant and not covered by snow, compute heat and  !
   ! moisture fluxes from vegetation to canopy, and flux resistance from soil or snow to   !
   ! canopy.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (veg_lai >= .1 .and. snowfac < .9) then
      !----- Save the vegetation class in an alias. ---------------------------------------!
      nveg = nint(leaf_class)
   
      !------------------------------------------------------------------------------------!
      !     Find the wind speed at the top of the canopy, using the similarity theory.     !
      !------------------------------------------------------------------------------------!
      veg_wind = leaf_reduced_wind(ustar,zeta,ribulk,geoht,0.0,veg_height,patch_rough)


      !------------------------------------------------------------------------------------!
      !    Find the aerodynamic conductances for heat and water at the leaf boundary       !
      ! layer.                                                                             !
      !------------------------------------------------------------------------------------!
      call leaf_aerodynamic_conductances(nveg,veg_wind,veg_temp,can_temp,can_rvap,can_rhos)
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
      gsw_inf = can_rhos * gsw_max(nveg) / (ftlo * fthi * frad * fvpd * fswp) / slai
      !------------------------------------------------------------------------------------!

      !----- 15-minute response time for stomatal conductance (must have dtll <= 900.). ---!
      gsw     = gsw + dtll * (gsw_inf - gsw) / 900.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Transpiration can only happen if there is water.                                !
      !------------------------------------------------------------------------------------!
      if (available_water > 0.) then

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
      wtemp         = veg_water - wflxvc * dtll + pcpgl * vf
      !----- If this result will lead to negative solution, then reduce evaporation. ------!
      if (wtemp < 0.0005 * stai) then
         wshed_loc = 0.
         wflxvc    = (veg_water + pcpgl * vf) / dtll
         veg_water = 0.
      else
         wshed_loc = max(0.,wtemp - 0.2 * stai)
         wshed_tot = wshed_tot + wshed_loc
         veg_water = wtemp - wshed_loc
      end if
      !------------------------------------------------------------------------------------!

      !------ Find the associated latent heat flux from vegetation to canopy. -------------!
      qwflxvc     = wflxvc     * (alvi - alli * veg_fliq)
      qtransp_loc = transp_loc *  alvl !----- Liquid phase only in transpiration. ---------!
      !------------------------------------------------------------------------------------!


      !----- Exchange heat between vegetation and precipitation in implicit scheme --------!
      qwshed_loc = wshed_loc * (     veg_fliq  * cliq * (veg_temp - tsupercool)            &
                               + (1.-veg_fliq) * cice *  veg_temp )
      qwshed_tot = qwshed_tot + qwshed_loc

      !------------------------------------------------------------------------------------!
      !      Update vegetation internal energy from radiative fluxes and sensible and      !
      ! latent heat transfer with canopy air.                                              !
      !------------------------------------------------------------------------------------!
      old_veg_energy = veg_energy
      veg_energy     = veg_energy                                                          &
                     + dtll * ( rshort_v + rlonga_v + rlonggs_v - rlongv_gs - rlongv_a     &
                              - hflxvc   - qwflxvc   - qtransp_loc)                        &
                     + qpcpgl * vf - qwshed_loc

      !------------------------------------------------------------------------------------!
      !     Find the atmosphere -> canopy fluxes.                                          !
      !------------------------------------------------------------------------------------!
      rho_ustar  = can_rhos  * ustar
      eflxac     = rho_ustar * estar * cp * can_temp ! Enthalpy exchange
      hflxac     = rho_ustar * tstar * can_exner     ! Sensible heat exchange
      wflxac     = rho_ustar * rstar                 ! Water vapour exchange
      cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
      !------------------------------------------------------------------------------------!




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

      !----- Find the vegetation diagnostic variables for next time step. -----------------!
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)

      !----- Find the canopy air diagnostic variables for next time step. -----------------!
      can_theta = exp(can_lntheta)
      can_shv   = can_rvap / (can_rvap + 1.)
      can_temp  = thrhsh2temp(can_theta,can_rhos,can_shv)
      can_prss  = can_rhos * rdry * can_temp * (1. + epim1 * can_shv)
      can_exner = cp  * (p00i * can_prss) ** rocp
      can_theiv = thetaeiv(can_theta,can_prss,can_temp,can_rvap,can_rvap,-8)
      !------------------------------------------------------------------------------------!


      if (can_lntheta /= can_lntheta .or. can_theiv /= can_theiv) then
         nsoil = nint(soil_text(mzg))
         write (unit=*,fmt=fmtc) '------------ Canopy theta_Eiv is screwed. ------------'
         write (unit=*,fmt=fmtc) ' - Vegetated patch. '
         write (unit=*,fmt=fmti) ' - SOIL_TEXT        = ',nsoil
         write (unit=*,fmt=fmti) ' - LEAF_CLASS       = ',nveg
         write (unit=*,fmt=fmti) ' - KSN              = ',ksn
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - CAN_THEIV        = ',can_theiv
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THEIV    = ',old_can_theiv
         write (unit=*,fmt=fmtf) ' - CAN_RVAP         = ',can_rvap
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RVAP     = ',old_can_rvap
         write (unit=*,fmt=fmtf) ' - CAN_CO2          = ',can_co2
         write (unit=*,fmt=fmtf) ' - CAN_PRSS         = ',can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_PRSS     = ',old_can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THETA    = ',old_can_theta
         write (unit=*,fmt=fmtf) ' - OLD_CAN_TEMP     = ',old_can_temp
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - VEG_ENERGY       = ',veg_energy
         write (unit=*,fmt=fmtf) ' - OLD_VEG_ENERGY   = ',old_veg_energy
         write (unit=*,fmt=fmtf) ' - VEG_WATER        = ',veg_water
         write (unit=*,fmt=fmtf) ' - OLD_VEG_WATER    = ',old_veg_water
         write (unit=*,fmt=fmtf) ' - VEG_HCAP         = ',veg_hcap
         write (unit=*,fmt=fmtf) ' - VEG_TEMP         = ',veg_temp
         write (unit=*,fmt=fmtf) ' - VEG_FLIQ         = ',veg_fliq
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - GROUND_RVAP      = ',ground_rvap
         write (unit=*,fmt=fmtf) ' - GROUND_RSAT      = ',ground_rsat
         write (unit=*,fmt=fmtf) ' - GROUND_TEMP      = ',ground_temp
         write (unit=*,fmt=fmtf) ' - GROUND_FLIQ      = ',ground_fliq
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - HFLXGC           = ',hflxgc
         write (unit=*,fmt=fmtf) ' - HFLXVC           = ',hflxvc
         write (unit=*,fmt=fmtf) ' - HFLXAC           = ',hflxac
         write (unit=*,fmt=fmtf) ' - EFLXAC           = ',eflxac
         write (unit=*,fmt=fmtf) ' - QWFLXGC          = ',qwflxgc
         write (unit=*,fmt=fmtf) ' - QDEWGNDFLX       = ',qdewgndflx
         write (unit=*,fmt=fmtf) ' - QWFLXVC          = ',qwflxvc
         write (unit=*,fmt=fmtf) ' - QTRANSP_LOC      = ',qtransp_loc
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - RSHORT_V         = ',rshort_v
         write (unit=*,fmt=fmtf) ' - RLONGA_V         = ',rlonga_v
         write (unit=*,fmt=fmtf) ' - RLONGGS_V        = ',rlonggs_v
         write (unit=*,fmt=fmtf) ' - RLONGV_GS        = ',rlongv_gs
         write (unit=*,fmt=fmtf) ' - RLONGV_A         = ',rlongv_a
         write (unit=*,fmt=fmtf) ' - QPCPGL           = ',qpcpgl
         write (unit=*,fmt=fmtf) ' - VF               = ',vf
         write (unit=*,fmt=fmtf) ' - QWSHED_LOC       = ',qwshed_loc
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - WFLXGC           = ',wflxgc
         write (unit=*,fmt=fmtf) ' - DEWGNDFLX        = ',dewgndflx
         write (unit=*,fmt=fmtf) ' - TRANSP_LOC       = ',transp_loc
         write (unit=*,fmt=fmtf) ' - WFLXAC           = ',wflxac
         write (unit=*,fmt=fmtf) ' - WFLXVC           = ',wflxvc
         write (unit=*,fmt=fmtf) ' - PCPGL            = ',pcpgl
         write (unit=*,fmt=fmtf) ' - WSHED_LOC        = ',wshed_loc
         write (unit=*,fmt=fmtf) ' - VF               = ',vf
         write (unit=*,fmt=fmtc) '-----------------------------------------------------'
      end if

   else
      !------------------------------------------------------------------------------------!
      !    No vegetation, or vegetation buried in snow...                                  !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the atmosphere -> canopy fluxes.                                          !
      !------------------------------------------------------------------------------------!
      rho_ustar  = can_rhos  * ustar
      eflxac     = rho_ustar * estar * cp * can_temp ! Enthalpy exchange
      hflxac     = rho_ustar * tstar * can_exner     ! Sensible heat exchange
      wflxac     = rho_ustar * rstar                 ! Water vapour exchange
      cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
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

      !----- Find the diagnostic canopy air variables for next time step. -----------------!

      !----- Find the canopy air diagnostic variables for next time step. -----------------!
      can_theta = exp(can_lntheta)
      can_shv   = can_rvap / (can_rvap + 1.)
      can_temp  = thrhsh2temp(can_theta,can_rhos,can_shv)
      can_prss  = can_rhos * rdry * can_temp * (1. + epim1 * can_shv)
      can_exner = cp  * (p00i * can_prss) ** rocp
      can_theiv = thetaeiv(can_theta,can_prss,can_temp,can_rvap,can_rvap,-8)
      !------------------------------------------------------------------------------------!


      !----- Update the vegetation temperature. -------------------------------------------!
      veg_temp = can_temp
      if (veg_temp > t3ple) then
         veg_fliq = 1.0
      elseif (veg_temp < t3ple) then
         veg_fliq = 0.0
      else
         veg_fliq = 0.5
      end if
      veg_energy = 0.0
      veg_water  = 0.0
      !------------------------------------------------------------------------------------!




      if (can_lntheta /= can_lntheta .or. can_theiv /= can_theiv) then
         write (unit=*,fmt=fmtc) '------------ Canopy theta_Eiv is screwed. ------------'
         write (unit=*,fmt=fmtc) ' - Non-vegetated patch. '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - CAN_THEIV        = ',can_theiv
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THEIV    = ',old_can_theiv
         write (unit=*,fmt=fmtf) ' - CAN_RVAP         = ',can_rvap
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RVAP     = ',old_can_rvap
         write (unit=*,fmt=fmtf) ' - CAN_CO2          = ',can_co2
         write (unit=*,fmt=fmtf) ' - CAN_PRSS         = ',can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_PRSS     = ',old_can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THETA    = ',old_can_theta
         write (unit=*,fmt=fmtf) ' - OLD_CAN_TEMP     = ',old_can_temp
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - HFLXGC           = ',hflxgc
         write (unit=*,fmt=fmtf) ' - EFLXAC           = ',eflxac
         write (unit=*,fmt=fmtf) ' - QWFLXGC          = ',qwflxgc
         write (unit=*,fmt=fmtf) ' - QDEWGNDFLX       = ',qdewgndflx
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - WFLXGC           = ',wflxgc
         write (unit=*,fmt=fmtf) ' - DEWGNDFLX        = ',dewgndflx
         write (unit=*,fmt=fmtf) ' - WFLXAC           = ',wflxac
         write (unit=*,fmt=fmtc) '-----------------------------------------------------'
      end if
   end if

   return
end subroutine leaf_canopy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of canopy air space variables for a given     !
! patch, before the time step iteration loop.                                              !
!------------------------------------------------------------------------------------------!
subroutine leaf_can_init(ip,can_theta,can_rvap,leaf_class,can_prss)
   use leaf_coms , only : atm_prss      & ! intent(in)
                        , atm_theta     & ! intent(in)
                        , atm_shv       & ! intent(in)
                        , geoht         & ! intent(in)
                        , veg_ht        & ! intent(in)
                        , can_exner     & ! intent(in)
                        , can_lntheta   & ! intent(in)
                        , can_shv       & ! intent(out)
                        , can_depth     & ! intent(out)
                        , can_depth_min & ! intent(out)
                        , can_temp      & ! intent(out)
                        , can_rhos      ! ! intent(out)
   use rconstants, only : cp            & ! intent(in)
                        , cpi           & ! intent(in)
                        , p00i          & ! intent(in)
                        , rocp          ! ! intent(in)
   use therm_lib , only : reducedpress  & ! function
                        , idealdenssh   ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)    :: ip
   real   , intent(in)    :: can_theta
   real   , intent(in)    :: can_rvap
   real   , intent(in)    :: leaf_class
   real   , intent(out)   :: can_prss
   !----- Local variables. ----------------------------------------------------------------!
   integer                :: nveg
   !---------------------------------------------------------------------------------------!




   !----- Define the canopy depth. --------------------------------------------------------!
   nveg      = nint(leaf_class)
   if (nveg == 0 .or. ip == 1) then
      can_depth = can_depth_min
   else
      can_depth = max(can_depth_min,veg_ht(nveg))
   end if
   !---------------------------------------------------------------------------------------!



   !----- Set initial value of specific humidity ------------------------------------------!
   can_shv     = can_rvap / (can_rvap + 1.)
   !---------------------------------------------------------------------------------------!




   !----- Update canopy air "dry enthropy" (log of theta). --------------------------------!
   can_lntheta            = log(can_theta)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update canopy air pressure and the Exner function here.  Canopy air pressure is    !
   ! assumed to remain constant during one LEAF full timestep, which means that heat       !
   ! equals to change in enthalpy.  Between time steps, atmospheric pressure changes so    !
   ! enthalpy is no longer conserved.  Potential temperature and equivalent potential      !
   ! temperature, on the other hand, are conserved if no heat happens, which is the case   !
   ! between successive calls.  Therefore, we track theta and theta_eiv.                   !
   !---------------------------------------------------------------------------------------!
   can_prss  = reducedpress(atm_prss,atm_theta,atm_shv,geoht,can_theta,can_shv,can_depth)
   can_exner = cp  * (p00i * can_prss) ** rocp
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Find the canopy air temperature, and compute the canopy air density.  Even though  !
   ! this leads to small violations of the ideal gas law, we are also going to keep the    !
   ! canopy air density constant over one LEAF time step.  This is to avoid that non-      !
   ! linearities sprout and cause problems to the code.                                    !
   !---------------------------------------------------------------------------------------!
   can_temp = cpi * can_theta * can_exner
   can_rhos = idealdenssh (can_prss,can_temp,can_shv)
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaf_can_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of leaf-level variables for a given patch,    !
! before the time step iteration loop.                                                     !
!------------------------------------------------------------------------------------------!
subroutine leaf_veg_init(ip,veg_energy,veg_water,veg_hcap,leaf_class,veg_lai,patch_area)
   use leaf_coms , only : can_temp     & ! intent(in)
                        , tiny_parea   & ! intent(in)
                        , veg_temp     & ! intent(out)
                        , veg_fliq     ! ! intent(out)
   use rconstants, only : t3ple        ! ! intent(in)
   use therm_lib , only : qwtk         ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)    :: ip
   real   , intent(inout) :: veg_energy
   real   , intent(inout) :: veg_water
   real   , intent(in)    :: veg_hcap
   real   , intent(in)    :: leaf_class
   real   , intent(in)    :: veg_lai
   real   , intent(in)    :: patch_area
   !----- Local variables. ----------------------------------------------------------------!
   integer                :: nveg
   !---------------------------------------------------------------------------------------!




   !----- Define the canopy depth. --------------------------------------------------------!
   nveg      = nint(leaf_class)
   if (nveg == 0 .or. ip == 1 .or. patch_area < tiny_parea) then
      veg_temp = can_temp
      if (veg_temp > t3ple) then
         veg_fliq = 1.0
      elseif (veg_temp < t3ple) then
         veg_fliq = 0.0
      else
         veg_fliq = 0.5
      end if
      veg_energy = 0.0
      veg_water  = 0.0
   else
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf_veg_init
!==========================================================================================!
!==========================================================================================!
