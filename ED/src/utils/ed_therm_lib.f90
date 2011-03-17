!==========================================================================================!
!==========================================================================================!
module ed_therm_lib
  
   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! This function calculates the total heat capacity in (J/K/m2) of the cohort            !
   ! leaf biomass.  This function is primarily used to calculate leaf temperatures         !
   ! based on leaf energy.                                                                 !
   !                                                                                       !
   ! Inputs:                                                                               !
   ! + BLEAF        - the leaf biomass of the cohort in kgC/plant.  NOTE: it will be       !
   !                  converted to kg/m2, by the product of leaf_carbon in units kg/plant  !
   !                  and number of plants per square meter, and the conversion of carbon  !
   !                  to total biomass.  The right hand side of the main equation accounts !
   !                  for the mass of insterstitial water and its ability to hold energy.  !
   ! + BDEAD        - the structural stem biomass of the cohort in kgC/plant.              !
   ! + BALIVE       - the live tissue biomass of the cohort, in kgC/plant.                 !
   ! + NPLANTS      - the number of plants per m2.                                         !
   ! + PFT          - the plant functional type of the current cohort, which may serve     !
   !                  for defining different parameterizations of specific heat capacity   !
   ! + PHEN_STATUS  - this is probably redundant with the LAI check, but if phen_status is !
   !                  2, it means no leaves exist so we set hcapveg to 0.                  !
   !                                                                                       !
   ! These methods follow the ways of Gu et al. 2007, with the only difference that for    !
   ! non-green biomass we dropped the temperature dependence and assumed T=T3ple, just to  !
   ! make it simpler.  See the module pft_coms.f90 for a description of the parameters,    !
   ! and see ed_params.f90 for the setting of these parameters.                            !
   !                                                                                       !
   ! Reference:                                                                            !
   !                                                                                       !
   ! Gu, L., T. Meyers, S. G. Pallardy, 2007: Influences of biomass heat and biochemical   !
   !      energy storages on the land surface fluxes and radiative temperature.            !
   !      J. Geophys. Res., v. 112, doi: 10.1029/2006JD007425.                             !
   !---------------------------------------------------------------------------------------!
   real function calc_hcapveg(bleaf,bdead,balive,nplant,hite,pft,phen_status, bsapwood)
      use consts_coms          , only : cliq                ! ! intent(in)
      use pft_coms             , only : c_grn_leaf_dry      & ! intent(in)
                                      , wat_dry_ratio_grn   & ! intent(in)
                                      , c_ngrn_biom_dry     & ! intent(in)
                                      , wat_dry_ratio_ngrn  & ! intent(in)
                                      , delta_c             & ! intent(in)
                                      , C2B                 ! ! intent(in)
      use rk4_coms             , only : ibranch_thermo      ! ! intent(in)
      use allometry            , only : wood_biomass        ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in)    :: bleaf         ! Biomass of leaves              [kgC/plant]
      real    , intent(in)    :: bdead         ! Biomass of structural stem     [kgC/plant]
      real    , intent(in)    :: balive        ! Biomass of live tissue         [kgC/plant]
      real    , intent(in)    :: bsapwood      ! Biomass of sapwood             [kgC/plant]
      real    , intent(in)    :: nplant        ! Number of plants               [ plant/m2]
      real    , intent(in)    :: hite          ! Cohort mean height             [        m]
      integer , intent(in)    :: pft           ! Plant functional type          [     ----]
      integer , intent(in)    :: phen_status   ! Phenology status               [     ----]
      !----- Local variables --------------------------------------------------------------!
      real                    :: bwood         ! Wood biomass                   [kgC/plant]
      real                    :: spheat_leaf   ! Leaf specific heat             [   J/kg/K]
      real                    :: spheat_wood   ! Wood specific heat             [   J/kg/K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Here we decide whether we compute the branch heat capacity or not.              !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0)
         !----- Skip it, the user doesn't want to solve for branches. ---------------------!
         spheat_wood = 0.
         bwood       = 0.
      case default
         !----- Finding branch/twig specific heat and biomass. ----------------------------!
         spheat_wood = (c_ngrn_biom_dry(pft) + wat_dry_ratio_ngrn(pft) * cliq)             &
                     / (1. + wat_dry_ratio_ngrn(pft)) + delta_c(pft)
         bwood = wood_biomass(bdead,balive,pft,hite, bsapwood)
      end select

      select case (phen_status)
      case (2)
         !---------------------------------------------------------------------------------!
         !     If phenology is 2 (i.e., no leaves), then the heat capacity is only due to  !
         ! the wood (it will be 0 in case wood is excluded).                               !
         !---------------------------------------------------------------------------------!
         calc_hcapveg = nplant * C2B * bwood * spheat_wood
      case default
         !---------------------------------------------------------------------------------!
         !     If phenology is either 0 or 1 (with leaves), then the heat capacity is the  !
         ! sum of the heat capacity due to wood and leaves.                                !
         !---------------------------------------------------------------------------------!
         !----- Finding the leaf specific heat. -------------------------------------------!
         spheat_leaf = (c_grn_leaf_dry(pft) + wat_dry_ratio_grn(pft) * cliq)               &
                     / (1. + wat_dry_ratio_grn(pft))
         !----- Then the heat capacity. ---------------------------------------------------!
         calc_hcapveg = nplant * C2B                                                       &
                      * ( bwood * spheat_wood * (1. + wat_dry_ratio_ngrn(pft))             &
                        + bleaf * spheat_leaf * (1. + wat_dry_ratio_grn(pft) ) )
      end select

      return
   end function calc_hcapveg
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine updates the vegetation energy when the plant heat capacity has     !
   ! changed. This routine should be used only when leaf or structural biomass has         !
   ! changed, it should never be used in fast time steps.                                  !
   !                                                                                       !
   !    The "cweh" mean "consistent water&energy&hcap" assumption                          !
   !---------------------------------------------------------------------------------------!
   subroutine update_veg_energy_cweh(csite,ipa,ico,old_hcapveg)
      use ed_state_vars, only : sitetype   & ! Structure
                              , patchtype  ! ! Structure
      use therm_lib    , only : qwtk       ! ! subroutine
      use consts_coms  , only : cliq       & ! intent(in)
                              , cice       & ! intent(in)
                              , tsupercool ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      integer        , intent(in) :: ico
      real           , intent(in) :: old_hcapveg
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      real(kind=8)                :: new_energy
      real                        :: new_temp
      real                        :: new_fliq
      integer                     :: kclosest
      integer                     :: k
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !    When heat capacity is zero (i.e., no leaves and we are ignoring branches), we   !
      ! cannot find the temperature using qwtk because it is a singularity.  Notice that   !
      ! this is different skipping when cohorts are not solvable... If the cohort is not   !
      ! solvable but still has some heat capacity, we should update internal energy using  !
      ! the traditional method, and NEVER force the heat capacity to be zero, otherwise we !
      ! violate the fact that heat capacity is a linear function of mass and this will     !
      ! cause problems during the fusion/splitting process.                                !
      !------------------------------------------------------------------------------------!
      if (cpatch%hcapveg(ico) == 0. ) then
         cpatch%veg_energy(ico) = 0.
         cpatch%veg_water(ico)  = 0.
         cpatch%veg_fliq(ico)   = 0.
         if (cpatch%hite(ico) > csite%total_snow_depth(ipa)) then
            !----- Plant is exposed, set temperature to the canopy temperature. -----------!
            cpatch%veg_temp(ico) = csite%can_temp(ipa)
         else
            !----- Find the snow layer that is the closest to where the leaves would be. --!
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) then
                  kclosest = k
               end if
            end do
            cpatch%veg_temp(ico) = csite%sfcwater_tempk(kclosest,ipa)
         end if

      !------------------------------------------------------------------------------------!
      !    Cohort heat capacity is not zero. Since we track vegetation temperature and     !
      ! liquid fraction of vegetation coating water, we can recalculate the internal       !
      ! energy by just switching the old heat capacity by the new one.                     !
      !------------------------------------------------------------------------------------!
      else
         cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)               &
                                + cpatch%veg_water(ico)                                    &
                                * ( cliq * cpatch%veg_fliq(ico)                            &
                                  * (cpatch%veg_temp(ico) - tsupercool)                    &
                                  + cice * (1.-cpatch%veg_fliq(ico)) * cpatch%veg_temp(ico))

         !----- This is a sanity check, it can be removed if it doesn't crash. ------------!
         call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)        &
                  ,new_temp,new_fliq)
      
         if (abs(new_temp - cpatch%veg_temp(ico)) > 0.1) then
            write (unit=*,fmt='(a)') ' ENERGY CONSERVATION FAILED!:'
            write (unit=*,fmt='(a,1x,es12.5)') ' Old temperature:  ',cpatch%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.5)') ' New temperature:  ',new_temp
            write (unit=*,fmt='(a,1x,es12.5)') ' Old heat capacity:',old_hcapveg
            write (unit=*,fmt='(a,1x,es12.5)') ' New heat capacity:',cpatch%hcapveg(ico)
            write (unit=*,fmt='(a,1x,es12.5)') ' Veg energy:       ',cpatch%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.5)') ' Veg water:        ',cpatch%veg_water(ico)
            call fatal_error('Energy is leaking!!!','update_veg_energy_cweh'               &
                            &,'ed_therm_lib.f90')
         end if
      end if

      return
   end subroutine update_veg_energy_cweh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine computes ground_shv, which is is the saturation mixing ratio of     !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   ! References:                                                                           !
   !                                                                                       !
   ! NP89 - Noilhan, J., S. Planton, 1989: A simple parameterization of land surface       !
   !        processes for meteorological models. Mon. Wea. Rev., 117, 536-549.             !
   !                                                                                       !
   ! LP92 - Lee, T. J., R. A. Pielke, 1992: Estimating the soil surface specific humidity  !
   !        J. Appl. Meteorol., 31, 480-484.                                               !
   !                                                                                       !
   ! LP93 - Lee, T. J., R. A. Pielke, 1993: CORRIGENDUM, J. Appl. Meteorol., 32, 580.      !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap(ksn,nsoil,topsoil_water,topsoil_temp,topsoil_fliq,sfcwater_temp   &
                        ,sfcwater_fliq,can_prss,can_shv,ground_shv,ground_ssh              &
                        ,ground_temp,ground_fliq)

      use soil_coms   , only : soil      & ! intent(in)
                             , betapower ! ! intent(in)
      use consts_coms , only : pi1       & ! intent(in)
                             , wdns      & ! intent(in)
                             , gorh2o    & ! intent(in)
                             , lnexp_min ! ! intent(in)
      use therm_lib   , only : rslif     ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer     , intent(in)  :: ksn           ! # of surface water layers    [     ----]
      integer     , intent(in)  :: nsoil         ! Soil type                    [     ----]
      real(kind=4), intent(in)  :: topsoil_water ! Top soil water               [m³_h2o/m³]
      real(kind=4), intent(in)  :: topsoil_temp  ! Top soil temperature         [        K]
      real(kind=4), intent(in)  :: topsoil_fliq  ! Top soil liquid water frac.  [       --]
      real(kind=4), intent(in)  :: sfcwater_temp ! Snow/water temperature       [        K]
      real(kind=4), intent(in)  :: sfcwater_fliq ! Snow/water liq. water frac.  [       --]
      real(kind=4), intent(in)  :: can_prss      ! canopy pressure              [       Pa]
      real(kind=4), intent(in)  :: can_shv       ! canopy vapour spec humidity  [kg_vap/kg]
      real(kind=4), intent(out) :: ground_shv    ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=4), intent(out) :: ground_ssh    ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=4), intent(out) :: ground_temp   ! Surface temperature          [        K]
      real(kind=4), intent(out) :: ground_fliq   ! Surface liquid water frac.   [       --]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)              :: slpotvn       ! soil water potential         [        m]
      real(kind=4)              :: alpha         ! alpha term (Lee-Pielke,1992) [     ----]
      real(kind=4)              :: beta          ! beta term  (Lee-Pielke,1992) [     ----]
      real(kind=4)              :: lnalpha       ! ln(alpha)                    [     ----]
      real(kind=4)              :: smterm        ! soil moisture term           [     ----]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide what are we going to call surface.                                      !
      !------------------------------------------------------------------------------------!
      select case (ksn)
      case (0)
         !---------------------------------------------------------------------------------!
         !      Without snowcover or water ponding, ground_shv is the effective specific   !
         ! humidity of soil and is used for soil evaporation.  This value is a combination !
         ! of the canopy air specific humidity, the saturation specific humidity at the    !
         ! soil temperature.  When the soil tends to dry air soil moisture, ground_shv     !
         ! tends to the canopy air space specific humidity, whereas it tends to the        !
         ! saturation value when the soil moisture is near or above field capacity.  These !
         ! tendencies will be determined by the alpha and beta parameters.                 !
         !---------------------------------------------------------------------------------!
         ground_temp = topsoil_temp
         ground_fliq = topsoil_fliq
         !----- Compute the saturation specific humidity at ground temperature. -----------!
         ground_ssh  = rslif(can_prss,ground_temp)
         ground_ssh  = ground_ssh / (1.0 + ground_ssh)
         !----- Determine alpha. ----------------------------------------------------------!
         slpotvn      = soil(nsoil)%slpots                                                 &
                      * (soil(nsoil)%slmsts / topsoil_water) ** soil(nsoil)%slbs
         lnalpha     = gorh2o * slpotvn / ground_temp
         if (lnalpha > lnexp_min) then
            alpha   = exp(lnalpha)
         else
            alpha   = 0.0
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Determine Beta, following NP89.  However, because we want evaporation to be !
         ! shut down when the soil approaches the dry air soil moisture, we offset both    !
         ! the soil moisture and field capacity to the soil moisture above dry air soil.   !
         ! This is necessary to avoid evaporation to be large just slightly above the dry  !
         ! air soil, which was happening especially for those soil types rich in clay.     !
         !---------------------------------------------------------------------------------!
         smterm     = (topsoil_water - soil(nsoil)%soilcp)                                 &
                    / (soil(nsoil)%sfldcap - soil(nsoil)%soilcp)
         beta       = (.5 * (1. - cos (min(1.,smterm) * pi1))) ** betapower
         !----- Use the expression from LP92 to determine the specific humidity. ----------!
         ground_shv = ground_ssh * alpha * beta + (1. - beta) * can_shv
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !    If a temporary layer exists, we use the top layer as the surface.  Since     !
         ! this is "pure" water or snow, we let it evaporate freely.  We can understand    !
         ! this as the limit of alpha and beta tending to one.                             !
         !---------------------------------------------------------------------------------!
         ground_temp = sfcwater_temp
         ground_fliq = sfcwater_fliq
         !----- Compute the saturation specific humidity at ground temperature. -----------!
         ground_ssh = rslif(can_prss,ground_temp)
         ground_ssh = ground_ssh / (1.0 + ground_ssh)
         !----- The ground specific humidity in this case is just the saturation value. ---!
         ground_shv  = ground_ssh
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine ed_grndvap
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine computes ground_shv, which is is the saturation mixing ratio of     !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   ! References:                                                                           !
   !                                                                                       !
   ! NP89 - Noilhan, J., S. Planton, 1989: A simple parameterization of land surface       !
   !        processes for meteorological models. Mon. Wea. Rev., 117, 536-549.             !
   !                                                                                       !
   ! LP92 - Lee, T. J., R. A. Pielke, 1992: Estimating the soil surface specific humidity  !
   !        J. Appl. Meteorol., 31, 480-484.                                               !
   !                                                                                       !
   ! LP93 - Lee, T. J., R. A. Pielke, 1993: CORRIGENDUM, J. Appl. Meteorol., 32, 580.      !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap8(ksn,nsoil,topsoil_water,topsoil_temp,topsoil_fliq,sfcwater_temp  &
                         ,sfcwater_fliq,can_prss,can_shv,ground_shv,ground_ssh             &
                         ,ground_temp,ground_fliq)
      use soil_coms   , only : soil8      & ! intent(in)
                             , betapower8 ! ! intent(in)
      use consts_coms , only : pi18       & ! intent(in)
                             , wdns8      & ! intent(in)
                             , gorh2o8    & ! intent(in)
                             , lnexp_min8 ! ! intent(in) 
      use therm_lib8  , only : rslif8     ! ! function

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer     , intent(in)  :: ksn           ! # of surface water layers    [     ----]
      integer     , intent(in)  :: nsoil         ! Soil type                    [     ----]
      real(kind=8), intent(in)  :: topsoil_water ! Top soil water               [m³_h2o/m³]
      real(kind=8), intent(in)  :: topsoil_temp  ! Top soil temperature         [        K]
      real(kind=8), intent(in)  :: topsoil_fliq  ! Top soil liquid water frac.  [       --]
      real(kind=8), intent(in)  :: sfcwater_temp ! Snow/water temperature       [        K]
      real(kind=8), intent(in)  :: sfcwater_fliq ! Snow/water liq. water frac.  [       --]
      real(kind=8), intent(in)  :: can_prss      ! canopy pressure              [       Pa]
      real(kind=8), intent(in)  :: can_shv       ! canopy vapour spec humidity  [kg_vap/kg]
      real(kind=8), intent(out) :: ground_shv    ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=8), intent(out) :: ground_ssh    ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=8), intent(out) :: ground_temp   ! Surface temperature          [        K]
      real(kind=8), intent(out) :: ground_fliq   ! Surface liquid water frac.   [       --]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)              :: slpotvn       ! soil water potential         [        m]
      real(kind=8)              :: alpha         ! alpha term (Lee-Pielke,1992) [     ----]
      real(kind=8)              :: beta          ! beta term  (Lee-Pielke,1992) [     ----]
      real(kind=8)              :: lnalpha       ! ln(alpha)                    [     ----]
      real(kind=8)              :: smterm        ! soil moisture term           [     ----]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide what are we going to call surface.                                      !
      !------------------------------------------------------------------------------------!
      select case (ksn)
      case (0)
         !---------------------------------------------------------------------------------!
         !      Without snowcover or water ponding, ground_shv is the effective specific   !
         ! humidity of soil and is used for soil evaporation.  This value is a combination !
         ! of the canopy air specific humidity, the saturation specific humidity at the    !
         ! soil temperature.  When the soil tends to dry air soil moisture, ground_shv     !
         ! tends to the canopy air space specific humidity, whereas it tends to the        !
         ! saturation value when the soil moisture is near or above field capacity.  These !
         ! tendencies will be determined by the alpha and beta parameters.                 !
         !---------------------------------------------------------------------------------!
         ground_temp = topsoil_temp
         ground_fliq = topsoil_fliq
         !----- Compute the saturation specific humidity at ground temperature. -----------!
         ground_ssh  = rslif8(can_prss,ground_temp)
         ground_ssh  = ground_ssh / (1.d0 + ground_ssh)
         !----- Determine alpha. ----------------------------------------------------------!
         slpotvn      = soil8(nsoil)%slpots                                                &
                      * (soil8(nsoil)%slmsts / topsoil_water) ** soil8(nsoil)%slbs
         lnalpha     = gorh2o8 * slpotvn / ground_temp
         if (lnalpha > lnexp_min8) then
            alpha   = exp(lnalpha)
         else
            alpha   = 0.d0
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Determine Beta, following NP89.  However, because we want evaporation to be !
         ! shut down when the soil approaches the dry air soil moisture, we offset both    !
         ! the soil moisture and field capacity to the soil moisture above dry air soil.   !
         ! This is necessary to avoid evaporation to be large just slightly above the dry  !
         ! air soil, which was happening especially for those soil types rich in clay.     !
         !---------------------------------------------------------------------------------!
         smterm     = (topsoil_water - soil8(nsoil)%soilcp)                                &
                    / (soil8(nsoil)%sfldcap - soil8(nsoil)%soilcp)
         beta       = (5.d-1 * (1.d0 - cos (min(1.d0,smterm) * pi18))) ** betapower8
         !----- Use the expression from LP92 to determine the specific humidity. ----------!
         ground_shv = ground_ssh * alpha * beta + (1.d0 - beta) * can_shv
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !    If a temporary layer exists, we use the top layer as the surface.  Since     !
         ! this is "pure" water or snow, we let it evaporate freely.  We can understand    !
         ! this as the limit of alpha and beta tending to one.                             !
         !---------------------------------------------------------------------------------!
         ground_temp = sfcwater_temp
         ground_fliq = sfcwater_fliq
         !----- Compute the saturation specific humidity at ground temperature. -----------!
         ground_ssh = rslif8(can_prss,ground_temp)
         ground_ssh = ground_ssh / (1.d0 + ground_ssh)
         !----- The ground specific humidity in this case is just the saturation value. ---!
         ground_shv  = ground_ssh
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine ed_grndvap8
   !=======================================================================================!
   !=======================================================================================!
 end module ed_therm_lib
!==========================================================================================!
!==========================================================================================!
