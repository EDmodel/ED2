module ed_therm_lib
  
contains

   !=======================================================================================!
   !=======================================================================================!
   real function calc_hcapveg(leaf_carbon,structural_carbon,nplants,pft)
     
     ! This function calculates the total heat capacity in (J/K) of the cohort
     ! biomass.  This function is primarily used to calculate leaf temperatures
     ! based on leaf energy.  At present, heat capacity is primarily based off
     ! of leaf biomass, but is also accounting for a variable fraction of stems.
     ! This is considered partially for stability.

     ! Inputs:
     ! leaf_biomass - the leaf biomass of the cohort in kg/m2 - NOTE
     ! you are most likely going to be passing the variables bleaf and bdead
     ! , but you must scale (multiply) them by the plant density (nplant) first.
     ! structural_biomass - the biomass of live and dead hardwood in kg/m2
     ! pft - the plant functional type of the current cohort, which may serve
     ! for defining different parameterizations of specific heat capacity

     ! These methods follow the ways of Gu et al. 2007. See the module pft_coms.f90
     ! for a description of the parameters, and see ed_params.f90 for the setting
     ! of these parameters.  The parameters are derived from secondary sources
     ! used by Gu et al. 2007 and the references are listed where parameters
     ! are defined.

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
!     real :: spec_hcap_leaf
!     real :: spec_hcap_stem
     

     real,parameter :: biomass_factor = 1.00 ! This is a biomass kluge factor
                                             ! The model is much faster and more stable
                                             ! When heat capacity is high.
                                             ! It was also found that when net-radiation
                                             ! matched tower observations, the dynamic
                                             ! range of the 

     real,parameter :: min_hcapveg = 30.0    ! This is roughly 10 g of biomass at 3000 J/kg/K
                                             ! Dont be fooled, this is quite high
     integer :: pft     
     real,parameter :: veg_temp = 285.0      ! RIght now we are using a nominal vegetation
                                             ! temperature, but 

     ! Old Method
     ! calc_hcapveg = hcapveg_ref * max(canopy_height,heathite_min)*cohort_lai/patch_lai

     ! Specific heat capacity of leaf biomass (J/kg/K)
     ! The calculation of leaf heat capacity follows Gu et al. 2007
!     spec_hcap_leaf  = (c_grn_leaf_dry(pft) + cliq*wat_dry_ratio_grn(pft))/(1.+wat_dry_ratio_grn(pft))
     
     
     ! Specific heat capacity of the stems and structural biomass.
     ! Also following Gu et al. 2007
     
!     spec_hcap_stem  = (c_ngrn_biom_dry(pft) + cliq*wat_dry_ratio_ngrn(pft))/(1+wat_dry_ratio_ngrn(pft))&
!          + 100. * wat_dry_ratio_ngrn(pft) * &
!          (-0.06191 + 2.36*1.0e-4*veg_temp - 1.33*1.0e-2*wat_dry_ratio_ngrn(pft))
     
     

     calc_hcapveg = nplants *                                                              &
                  (  leaf_carbon*C2B*c_grn_leaf_dry(pft)                                   &
                   + wat_dry_ratio_grn(pft)*leaf_carbon*C2B*cliq                           &
                  )



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
        sfcwater_energy, rhos, can_shv, ground_shv, surface_ssh)
     
     use soil_coms,   only: ed_nstyp, soil
     use grid_coms,   only: nzg
     use consts_coms, only:  pi1, grav, rvap, wdns
     use therm_lib  , only: rhovsil, qtk, qwtk8
     
     implicit none
     
     integer, intent(in) :: nlev_sfcwater ! # active levels of surface water
     integer, intent(in) :: nts           ! soil textural class (local name)
     
     real(kind=8), intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
     real, intent(in)  :: soil_energy     ! [J/m^3]
     real, intent(in)  :: sfcwater_energy ! [J/kg]
     real, intent(in)  :: rhos            ! air density [kg/m^3]
     real, intent(in)  :: can_shv         ! canopy vapor spec hum [kg_vap/kg_air]
     real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
     real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]
     
     
     real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

     
     ! Local variables

     real :: slpotvn ! soil water potential [m]
     real :: alpha   ! "alpha" term in Lee and Pielke (1993)
     real :: beta    ! "beta" term in Lee and Pielke (1993)
     real :: tempk   ! surface water temp [K]
     real :: fracliq ! fraction of surface water in liquid phase
     
     ! surface_ssh is the saturation mixing ratio of the top soil or snow surface
     ! and is used for dew formation and snow evaporation.
     
     if (nlev_sfcwater > 0) then
        call qtk(sfcwater_energy,tempk,fracliq)
        surface_ssh = rhovsil(tempk) / rhos
     else
        
        ! Without snowcover, ground_shv is the effective saturation mixing
        ! ratio of soil and is used for soil evaporation.  First, compute the
        ! "alpha" term or soil "relative humidity" and the "beta" term.
        
        call qwtk8(soil_energy,soil_water*dble(wdns),soil(nts)%slcpd,tempk,fracliq)
        surface_ssh = rhovsil(tempk) / rhos
        
        slpotvn = soil(nts)%slpots * (soil(nts)%slmsts / sngl(soil_water)) ** soil(nts)%slbs
        alpha = exp(gorvap * slpotvn / tempk)
        beta = .25 * (1. - cos (min(1.,sngl(soil_water) / soil(nts)%sfldcap) * pi1)) ** 2
        ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv
        
     endif
     
     return
   end subroutine ed_grndvap


 end module ed_therm_lib
