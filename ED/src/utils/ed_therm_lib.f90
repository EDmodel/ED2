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
   real function calc_hcapveg(bleaf,bdead,balive,nplant,hite,pft,phen_status)
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
         bwood = wood_biomass(bdead,balive,pft,hite)
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
   !      This routine computes surface_ssh, which is is the saturation mixing ratio of    !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap(nlev_sfcwater,nts,soil_water,soil_energy,sfcwater_energy,rhos     &
                        ,can_shv,ground_shv,surface_ssh,surface_tempk,surface_fracliq)
     
      use soil_coms   , only: ed_nstyp,soil
      use grid_coms   , only: nzg
      use consts_coms , only: pi1,wdns,gorvap
      use therm_lib   , only: rhovsil,qtk,qwtk,qwtk
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
      integer, intent(in)  :: nts             ! soil textural class (local name)
     
      real   , intent(in)  :: soil_water      ! soil water content              [m訛h2o/m設
      real   , intent(in)  :: soil_energy     ! Soil internal energy            [     J/m設
      real   , intent(in)  :: sfcwater_energy ! Snow/water internal energy      [     J/kg]
      real   , intent(in)  :: rhos            ! air density                     [    kg/m設
      real   , intent(in)  :: can_shv         ! canopy vapor spec humidity      [kg_vap/kg]
      real   , intent(out) :: ground_shv      ! ground equilibrium spec hum     [kg_vap/kg]
      real   , intent(out) :: surface_ssh     ! surface (saturation) spec hum   [kg_vap/kg]
      real   , intent(out) :: surface_tempk   ! Surface water temperature       [        K]
      real   , intent(out) :: surface_fracliq ! fraction of surface water in liquid phase
      !----- Local variables --------------------------------------------------------------!
      real                 :: slpotvn         ! soil water potential            [        m]
      real                 :: alpha           ! "alpha" term in Lee and Pielke (1993)
      real                 :: beta            ! "beta" term in Lee and Pielke (1993)
      real                 :: lnalpha         ! ln(alpha)
      !------------------------------------------------------------------------------------!
      
      if (nlev_sfcwater > 0 .and. sfcwater_energy > 0.) then
         !----- If a temporary layer exists, this is the surface. -------------------------!
         call qtk(sfcwater_energy,surface_tempk,surface_fracliq)
         surface_ssh = rhovsil(surface_tempk) / rhos
      else
         !---------------------------------------------------------------------------------!
         !      Without snowcover, ground_shv is the effective saturation mixing ratio of  !
         ! soil and is used for soil evaporation.  First, compute the "alpha" term or soil !
         ! "relative humidity" and the "beta" term.                                        !
         !---------------------------------------------------------------------------------!
         call qwtk(soil_energy,soil_water*wdns,soil(nts)%slcpd,surface_tempk               &
                  ,surface_fracliq)
         surface_ssh = rhovsil(surface_tempk) / rhos
         slpotvn     = soil(nts)%slpots * (soil(nts)%slmsts / soil_water) ** soil(nts)%slbs
         lnalpha     = gorvap * slpotvn / surface_tempk
         if (lnalpha > -38.) then
            alpha       = exp(lnalpha)
         else
            alpha       = 0.0
         end if
         beta        = .25 * (1. - cos (min(1.,soil_water / soil(nts)%sfldcap) * pi1)) ** 2
         ground_shv  = surface_ssh * alpha * beta + (1. - beta) * can_shv
      end if
      return
   end subroutine ed_grndvap
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine computes surface_ssh, which is is the saturation mixing ratio of    !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   ! This uses double precision arguments for real numbers.                                !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap8(nlev_sfcwater,nts,soil_water,soil_energy,sfcwater_energy,rhos    &
                        ,can_shv,ground_shv,surface_ssh,surface_tempk,surface_fracliq)
      use soil_coms   , only: ed_nstyp,soil8
      use grid_coms   , only: nzg
      use consts_coms , only: pi18,wdns8,gorvap8
      use therm_lib   , only: rhovsil,qtk8,qwtk8
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
      integer, intent(in)  :: nts             ! soil textural class (local name)
     
      real(kind=8), intent(in)  :: soil_water      ! soil water content         [m訛h2o/m設
      real(kind=8), intent(in)  :: soil_energy     ! Soil internal energy       [     J/m設
      real(kind=8), intent(in)  :: sfcwater_energy ! Snow/water internal energy [     J/kg]
      real(kind=8), intent(in)  :: rhos            ! air density                [    kg/m設
      real(kind=8), intent(in)  :: can_shv         ! canopy vapor spec humidity [kg_vap/kg]
      real(kind=8), intent(out) :: ground_shv      ! Gnd. equilibrium spec hum  [kg_vap/kg]
      real(kind=8), intent(out) :: surface_ssh     ! Sfc. sat. spec hum         [kg_vap/kg]
      real(kind=8), intent(out) :: surface_tempk   ! Surface water temperature  [        K]
      real(kind=8), intent(out) :: surface_fracliq ! liquid fraction of surface water
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)              :: slpotvn         ! soil water potential       [        m]
      real(kind=8)              :: alpha           ! "alpha" term in Lee and Pielke (1993)
      real(kind=8)              :: beta            ! "beta" term in Lee and Pielke (1993)
      real(kind=8)              :: lnalpha         ! ln(alpha)
      !------------------------------------------------------------------------------------!
      
      if (nlev_sfcwater > 0 .and. sfcwater_energy > 0.d0) then
         !----- If a temporary layer exists, this is the surface. -------------------------!
         call qtk8(sfcwater_energy,surface_tempk,surface_fracliq)
         surface_ssh = dble(rhovsil(sngl(surface_tempk))) / rhos
      else
         !---------------------------------------------------------------------------------!
         !      Without snowcover, ground_shv is the effective saturation mixing ratio of  !
         ! soil and is used for soil evaporation.  First, compute the "alpha" term or soil !
         ! "relative humidity" and the "beta" term.                                        !
         !---------------------------------------------------------------------------------!
         call qwtk8(soil_energy,soil_water*wdns8,soil8(nts)%slcpd,surface_tempk            &
                  ,surface_fracliq)
         surface_ssh = dble(rhovsil(sngl(surface_tempk))) / rhos
         slpotvn     = soil8(nts)%slpots*(soil8(nts)%slmsts / soil_water)                  &
                     ** soil8(nts)%slbs
         lnalpha     = gorvap8 * slpotvn / surface_tempk
         if (lnalpha > -38.) then
            alpha       = exp(lnalpha)
         else
            alpha       = 0.d0
         end if
         beta        = 2.5d-1                                                              &
                     * (1.d0 - cos (min(1.d0,soil_water / soil8(nts)%sfldcap) * pi18)) ** 2
         ground_shv  = surface_ssh * alpha * beta + (1.d0 - beta) * can_shv
      end if
      return
   end subroutine ed_grndvap8
   !=======================================================================================!
   !=======================================================================================!
 end module ed_therm_lib
!==========================================================================================!
!==========================================================================================!
