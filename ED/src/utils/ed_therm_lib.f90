module ed_therm_lib

  
  
  !---------------------------------------------------------------------------------------!
  !     These constants came from the paper in which the saturation vapour pressure is    !
  ! based on:                                                                             !
  !                                                                                       !
  !  Murphy, D. M.; Koop, T., 2005: Review of the vapour pressures of ice and supercooled !
  !     water for atmospheric applications. Q. J. Royal Meteor. Soc., vol. 31, pp. 1539-  !
  !     1565 (hereafter MK05).                                                            !
  !                                                                                       !
  !  These equations give the triple point at t3ple, with vapour pressure being es3ple.   !
  !---------------------------------------------------------------------------------------!
  !----- Coefficients based on equation (7): ---------------------------------------------!
  real, dimension(0:3), parameter :: iii_7 = (/ 9.550426,-5723.265, 3.53068,-0.00728332 /)
  !----- Coefficients based on equation (10), first fit ----------------------------------!
  real, dimension(0:3), parameter :: l01_10= (/54.842763,-6763.22 ,-4.210  , 0.000367   /)
  !----- Coefficients based on equation (10), second fit ---------------------------------!
  real, dimension(0:3), parameter :: l02_10= (/53.878   ,-1331.22 ,-9.44523, 0.014025   /)
  !----- Coefficients based on the hyperbolic tangent ------------------------------------!
  real, dimension(2)  , parameter :: ttt_10= (/0.0415,218.8/)
  !---------------------------------------------------------------------------------------!
  
contains

   !=======================================================================================!
   !=======================================================================================!
   real function calc_hcapveg(leaf_carbon,structural_carbon,nplants,pft)
     
     ! This function calculates the total heat capacity in (J/K) of the cohort
     ! biomass.  This function is primarily used to calculate leaf temperatures
     ! based on leaf energy.  At present, heat capacity is primarily based off
     ! of leaf biomass, but is also can account for a variable fraction of stems.
     ! This is considered partially for stability.

     ! Inputs:
     ! leaf_biomass - the leaf biomass of the cohort in kg/m2 - NOTE
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
     

     real,parameter :: biomass_factor = 1.0 ! This is a biomass kluge factor
                                             ! The model is much faster and more stable
                                             ! When heat capacity is high.

     real,parameter :: min_hcapveg = 30.0    

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
     
     
     ! New Method
     !     calc_hcapveg = biomass_factor * nplants * (leaf_carbon*C2B*spec_hcap_leaf + &
     !          structural_carbon *C2B * hcap_stem_fraction * spec_hcap_stem )
     
     calc_hcapveg = nplants * (leaf_carbon*C2B*c_grn_leaf_dry(pft) + wat_dry_ratio_grn(pft)*leaf_carbon*C2B*cliq)




     calc_hcapveg = max(calc_hcapveg,min_hcapveg)
     
     
     return
     
   end function calc_hcapveg
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine update_veg_energy_cweh(cpatch,ico)

     ! Update the vegetation energy and heat capacity
     ! This version is usefull if the previous leaf water leaf energy
     ! and leaf heat capacities have not been updated or changed, ie
     ! they are still in balance. This is usefull, because we can use
     ! the previous quanitities to diagnose the leaf temperature, and 
     ! especially the leaf liquid fraction; and then use those quantities
     ! to re-calculate the leaf energy and heat capacity after
     ! we incorporate a change in biomass.
     ! The "cwe" mean "consistent water&energy&hcap" assumption
     
     use ed_state_vars,only: patchtype
     use consts_coms,only:cliq,cice,alli,t3ple
     use therm_lib, only : qwtk
     
     implicit none
     
     type(patchtype),target :: cpatch
     integer,intent(in) :: ico
     real :: fracliq

     ! First lets use the existing vegetation energy, to calculate the
     ! liquid fraction of leaf water
     
     
     call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico), &
          cpatch%hcapveg(ico),cpatch%veg_temp(ico),fracliq)
     
     
     ! Next, we assume that the leaf's biophysical quantities have changed, so
     ! we update the vegetation heat capacity based on these new values
     
     cpatch%hcapveg(ico) = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico), &
          cpatch%nplant(ico),cpatch%pft(ico))
     
     ! Now, calculate the new energy, based on the updated quantities
     
     cpatch%veg_energy(ico) = ( cpatch%veg_temp(ico)-t3ple)*cpatch%hcapveg(ico) + & ! U of the leaf tissue
          fracliq*cpatch%veg_water(ico)*alli + &                                    ! latent heat of fusion
          fracliq*cpatch%veg_water(ico)*cliq*( cpatch%veg_temp(ico) -t3ple) + &     ! thermal energy of any liquid
          (1.-fracliq)*cpatch%veg_water(ico)*cice*( cpatch%veg_temp(ico) -t3ple)     ! thermal energy of any ice
     
     return
   end subroutine update_veg_energy_cweh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine update_veg_energy_ct(cpatch,ico)

     ! Update the vegetation energy and heat capacity
     ! This version is usefull if man cohort level quanities
     ! have changed. If the biomass quanities have been updated, and the
     ! leaf water has changed, then we must use this version.
     ! This method is not as robust as the "cweh" version, and makes the following
     ! two assumptions:
     ! 1) That the temperature of the cohort is still representative of the actual
     ! termperature given that there may have been all sorts of changes to the cohort
     ! 2) That if the temperature is above or equal to the triple point, 
     ! than the liquid fraction is 1.0, if the temperature is below the triple point
     ! then there is no liquid at all.
     ! Again the other method is more robust, because it handles energies continuously
     ! whereas this scheme can not produce vegetation energies with mixed liquid
     ! fractions
     ! the "ct" means it makes a "consistent temperature" assumption
     
     
     use ed_state_vars,only:patchtype
     use consts_coms,only:cliq,cice,alli,t3ple
     
     implicit none
     
     type(patchtype),target :: cpatch
     integer,intent(in) :: ico
     
     ! Update the vegetation heat capacity based on these new values
     
     cpatch%hcapveg(ico) = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico), &
          cpatch%nplant(ico),cpatch%pft(ico))
     
     ! Now, calculate the new energy, based on the updated quantities
     
     if(cpatch%veg_temp(ico)>=t3ple) then
        
        cpatch%veg_energy(ico) = ( cpatch%veg_temp(ico)-t3ple)*cpatch%hcapveg(ico) + & ! U of the leaf tissue
             cpatch%veg_water(ico)*alli + &                                    ! latent heat of fusion
             cpatch%veg_water(ico)*cliq*( cpatch%veg_temp(ico) -t3ple)       ! thermal energy of any liquid
        
     else
        
        cpatch%veg_energy(ico) = ( cpatch%veg_temp(ico)-t3ple)*cpatch%hcapveg(ico) + & ! U of the leaf tissue
             cpatch%veg_water(ico)*cice*( cpatch%veg_temp(ico) -t3ple)     ! thermal energy of any ice
        
     endif
     
     return
   end subroutine update_veg_energy_ct

   !==========================================================================================!
   !==========================================================================================!
   
   subroutine ed_grndvap(nlev_sfcwater, nts, soil_water, soil_energy,    &
        sfcwater_energy, rhos, can_shv, ground_shv, surface_ssh)
     
     use soil_coms,   only: ed_nstyp, soil
     use grid_coms,   only: nzg
     use consts_coms, only:  pi1, grav, rvap
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
        
        call qwtk8(soil_energy,soil_water*1.d3,soil(nts)%slcpd,tempk,fracliq)
        surface_ssh = rhovsil(tempk) / rhos
        
        slpotvn = soil(nts)%slpots * (soil(nts)%slmsts / sngl(soil_water)) ** soil(nts)%slbs
        alpha = exp(gorvap * slpotvn / tempk)
        beta = .25 * (1. - cos (min(1.,sngl(soil_water) / soil(nts)%sfldcap) * pi1)) ** 2
        ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv
        
     endif
     
     return
   end subroutine ed_grndvap
   
   !================================================================================

  real function fast_svp(pres,temp)
     
     use consts_coms, only: t3ple,ep
     implicit none
     real,intent(in) :: temp
     real,intent(in) :: pres
     real            :: ttfun,fun2,fun1
     real :: esz

     if(temp<t3ple) then
        !----- Updated method, using MK05 ------------------------------------------------!
        fun1 = iii_7(0) + iii_7(1)/temp + iii_7(2) * log(temp) + iii_7(3) * temp
        esz  = exp(fun1)
     else
        fun1 = l01_10(0) + l01_10(1)/temp + l01_10(2)*log(temp) + l01_10(3) * temp
        fun2 = l02_10(0) + l02_10(1)/temp + l02_10(2)*log(temp) + l02_10(3) * temp
        ttfun = tanh(ttt_10(1) * (temp - ttt_10(2)))
        esz  = exp(fun1 + ttfun*fun2)
     endif
     
     fast_svp = ep * esz / (pres - esz)
     
   end function fast_svp


 end module ed_therm_lib
