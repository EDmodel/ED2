module ed_therm_lib
  
contains

   !=======================================================================================!
   !=======================================================================================!
   real function calc_hcapveg(leaf_carbon,structural_carbon,nplants,pft)
     
     ! This function calculates the total heat capacity in (J/K/m2) of the cohort
     ! leaf biomass.  This function is primarily used to calculate leaf temperatures
     ! based on leaf energy.

     ! Inputs:
     ! leaf_biomass - the leaf biomass of the cohort in kg/m2 - NOTE
     ! this is constructed by the product of leaf_carbon in units kg/plant
     ! and number of plants per square meter, and the conversion of carbon to 
     ! total biomass.  The right hand side of the main equation accounts for 
     ! the mass of insterstitial water and its ability to hold energy.
     ! pft - the plant functional type of the current cohort, which may serve
     ! for defining different parameterizations of specific heat capacity

     ! These methods follow the ways of Gu et al. 2007. See the module pft_coms.f90
     ! for a description of the parameters, and see ed_params.f90 for the setting
     ! of these parameters.

     use consts_coms,only : cliq

     use pft_coms,only:       &
          c_grn_leaf_dry,     &
          c_ngrn_biom_dry,    &
          wat_dry_ratio_grn,  &
          wat_dry_ratio_ngrn, &
          hcap_stem_fraction, &
          C2B


     implicit none
     
     real :: leaf_carbon         ! leaf biomass per plant
     real :: structural_carbon   ! structural biomass per plant
     real :: nplants             ! Number of plants per square meter
     integer :: pft
     real,parameter :: min_hcapveg = 200.0    

     
     calc_hcapveg = nplants*leaf_carbon*C2B*( c_grn_leaf_dry(pft) + wat_dry_ratio_grn(pft)*cliq)

     calc_hcapveg = max(calc_hcapveg,min_hcapveg)
     
     
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
   subroutine update_veg_energy_cweh(veg_energy,veg_temp,old_hcapveg,new_hcapveg)     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(inout) :: veg_energy
      real, intent(in)    :: veg_temp
      real, intent(in)    :: old_hcapveg
      real, intent(in)    :: new_hcapveg
      !------------------------------------------------------------------------------------!
      veg_energy = veg_energy + (new_hcapveg-old_hcapveg) * veg_temp
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
