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
   ! + LEAF_BIOMASS - the leaf biomass of the cohort in kg/plant.  NOTE: it will be        !
   !                  converted to kg/m2, by the product of leaf_carbon in units kg/plant  !
   !                  and number of plants per square meter, and the conversion of carbon  !
   !                  to total biomass.  The right hand side of the main equation accounts !
   !                  for the mass of insterstitial water and its ability to hold energy.  !
   ! + NPLANTS      - the number of plants per m2.                                         !
   ! + PFT          - the plant functional type of the current cohort, which may serve     !
   !                  for defining different parameterizations of specific heat capacity   !
   ! + LAI          - LAI, used to check whether we will keep track of this cohort heat    !
   !                  capacity or simply set it to zero (when LAI < LAI_MIN and we don't   !
   !                  prognose the leaf energy).                                           !
   ! + PHEN_STATUS  - this is probably redundant with the LAI check, but if phen_status is !
   !                  2, it means no leaves exist so we set hcapveg to 0.                  !
   !                                                                                       !
   ! These methods follow the ways of Gu et al. 2007. See the module pft_coms.f90          !
   ! for a description of the parameters, and see ed_params.f90 for the setting            !
   ! of these parameters.                                                                  !
   !---------------------------------------------------------------------------------------!
   real function calc_hcapveg(leaf_carbon,nplants,lai,pft,phen_status)
      use consts_coms          , only : cliq                ! ! intent(in)
      use pft_coms             , only : c_grn_leaf_dry      & ! intent(in)
                                      , wat_dry_ratio_grn   & ! intent(in)
                                      , C2B                 ! ! intent(in)
      use canopy_radiation_coms, only : lai_min             ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in)    :: leaf_carbon         ! Leaf biomass              [kg/plant]
      real    , intent(in)    :: nplants             ! Number of plants          [plant/m2]
      real    , intent(in)    :: lai                 ! Leaf area index           [   m2/m2]
      integer , intent(in)    :: pft                 ! Plant functional type     [    ----]
      integer , intent(in)    :: phen_status         ! Phenology status          [    ----]
      !------------------------------------------------------------------------------------!

      !----- Leaf heat capacity should be zero if there is no leaf ------------------------!
      if (phen_status == 2) then
         calc_hcapveg = 0.
      else
         !---------------------------------------------------------------------------------!
         !    Including an offset in the specific heat capacity. This is as bad as set-    !
         ! ting the minimum heat capacity in the sense that this is a kluge, but it has a  !
         ! few advantages...                                                               !
         ! 1. If there is no plant, the heat capacity will be 0.                           !
         ! 2. If patches are fused or split, then hcapveg can be really treated as an      !
         !    extensive quantity;                                                          !
         ! 3. The heat capacity becomes a smooth curve, continuous and differentiable.     !
         !---------------------------------------------------------------------------------!
         calc_hcapveg = nplants * leaf_carbon * C2B                                        &
                      * ( c_grn_leaf_dry(pft) + wat_dry_ratio_grn(pft)*cliq)
      end if

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
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      integer        , intent(in) :: ico
      real           , intent(in) :: old_hcapveg
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      real                        :: new_temp
      real                        :: new_fliq
      integer                     :: kclosest
      integer                     :: k
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !    If the heat capacity is zero we cannot use qwtk to find the temperature.  This  !
      ! will happen only when the plant has no leaves, so we set up the leaf temperature   !
      ! to be the canopy air temperature, or the snow temperature in case it is buried in  !
      ! snow.                                                                              !
      !------------------------------------------------------------------------------------!
      if (cpatch%hcapveg(ico) == 0. .or. cpatch%phenology_status(ico) == 2) then
         cpatch%hcapveg(ico)    = 0.
         cpatch%veg_energy(ico) = 0.
         cpatch%veg_water(ico)  = 0.
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
      elseif(cpatch%phenology_status(ico) < 2) then
         cpatch%veg_energy(ico) = cpatch%veg_energy(ico)                                   &
                                + (cpatch%hcapveg(ico)-old_hcapveg) * cpatch%veg_temp(ico)
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
      else 
            write (unit=*,fmt='(a)') ' WHAT AM I DOING HERE???:'
            call fatal_error('Entered in a forbidden place!!!','update_veg_energy_cweh'    &
                            &,'ed_therm_lib.f90')
      end if

      return
   end subroutine update_veg_energy_cweh
   !=======================================================================================!
   !=======================================================================================!






   !==========================================================================================!
   !==========================================================================================!
   
   subroutine ed_grndvap(nlev_sfcwater, nts, soil_water, soil_energy,    &
        sfcwater_energy,rhos, can_shv, ground_shv, surface_ssh, &
        surface_tempk,surface_fracliq)
     
     use soil_coms,   only: ed_nstyp, soil
     use grid_coms,   only: nzg
     use consts_coms, only:  pi1, wdns, gorvap
     use therm_lib  , only: rhovsil, qtk, qwtk, qwtk8
     
     implicit none
     integer, intent(in) :: nlev_sfcwater ! # active levels of surface water
     integer, intent(in) :: nts           ! soil textural class (local name)
     
     real(kind=8), intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
     real, intent(in)  :: soil_energy     ! [J/m³]
     real, intent(in)  :: sfcwater_energy ! [J/kg]
     real, intent(in)  :: rhos            ! air density [kg/m³]
     real, intent(in)  :: can_shv         ! canopy vapor spec hum [kg_vap/kg_air]
     real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
     real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]
     real, intent(out) :: surface_tempk   ! Surface water temperature [K]
     real, intent(out) :: surface_fracliq ! fraction of surface water in liquid phase
     
     ! Local variables

     real :: slpotvn ! soil water potential [m]
     real :: alpha   ! "alpha" term in Lee and Pielke (1993)
     real :: beta    ! "beta" term in Lee and Pielke (1993)
     
     ! surface_ssh is the saturation mixing ratio of the top soil or snow surface
     ! and is used for dew formation and snow evaporation.
     
     if (nlev_sfcwater > 0 .and. sfcwater_energy > 0.) then
        call qtk(sfcwater_energy,surface_tempk,surface_fracliq)
        surface_ssh = rhovsil(surface_tempk) / rhos
     else
        
        ! Without snowcover, ground_shv is the effective saturation mixing
        ! ratio of soil and is used for soil evaporation.  First, compute the
        ! "alpha" term or soil "relative humidity" and the "beta" term.
        
        call qwtk8(soil_energy,soil_water*dble(wdns),soil(nts)%slcpd &
                  ,surface_tempk,surface_fracliq)
        surface_ssh = rhovsil(surface_tempk) / rhos
        
        slpotvn = soil(nts)%slpots * (soil(nts)%slmsts / sngl(soil_water)) ** soil(nts)%slbs
        alpha = exp(gorvap * slpotvn / surface_tempk)
        beta = .25 * (1. - cos (min(1.,sngl(soil_water) / soil(nts)%sfldcap) * pi1)) ** 2
        ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv
        
     endif
     
     return
   end subroutine ed_grndvap


 end module ed_therm_lib
