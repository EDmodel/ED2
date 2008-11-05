! ===================================================

subroutine canopy_photosynthesis_ar(csite, ipa, vels, rhos, prss,   &
     ed_ktrans, ntext_soil, soil_water, soil_fracliq, lsl, sum_lai_rbi,  &
     leaf_aging_factor, green_leaf_factor)

  use ed_state_vars,only:sitetype,patchtype
  use max_dims, only: n_pft
  use pft_coms, only: leaf_width, water_conductance, q, qsw, include_pft
  use canopy_radiation_coms, only: lai_min
  use grid_coms, only: nzg
  use soil_coms, only: soil, dslz
  use consts_coms, only : t00,mmdov

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  real, intent(in) :: vels
  real, intent(in) :: rhos
  real, intent(in) :: prss
  integer, dimension(nzg), intent(out) :: ed_ktrans
  integer, dimension(nzg), intent(in) :: ntext_soil
  real(kind=8), dimension(nzg), intent(in) :: soil_water
  real, dimension(nzg), intent(in) :: soil_fracliq
  integer, intent(in) :: lsl
  real, intent(in), dimension(n_pft) :: leaf_aging_factor
  real, intent(in), dimension(n_pft) :: green_leaf_factor

  integer, dimension(nzg) :: root_depth_indices

  integer :: ipft
  real :: P_op
  real :: P_cl
  real :: leaf_resp
  real :: cumulative_lai
  real :: water_demand
  real :: water_supply
  integer :: k1
  integer :: k2
  real :: swp
  integer :: nts
  real :: slpotv
  real, dimension(nzg) :: available_liquid_water
  real, intent(out) :: sum_lai_rbi
  logical :: las

  las = .false.

  cpatch => csite%patch(ipa)

  ! calculate liquid water available for transpiration
  available_liquid_water(nzg) = max(0.0, 1.0e3 * dslz(nzg) *   &
       soil_fracliq(nzg) * (soil_water(nzg) - soil(ntext_soil(nzg))%soilcp))

  do k1 = nzg-1, lsl, -1
     available_liquid_water(k1) = available_liquid_water(k1+1) +  &
          dslz(k1) * 1.0e3 * soil_fracliq(k1) * max(0.0, soil_water(k1) -  &
          soil(ntext_soil(k1))%soilcp)
  enddo

  ! Initialize the array of maximum photosynthesis rates used in the 
  ! mortality function.

  csite%A_o_max(1:n_pft,ipa) = 0.0
  csite%A_c_max(1:n_pft,ipa) = 0.0

  ! Find the first cohort with leaves above snow
  
  cohort_with_leaves:do ico = 1,cpatch%ncohorts
     if(cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > csite%total_snow_depth(ipa)) then
        las = .true.
        exit cohort_with_leaves
        
     endif
  enddo cohort_with_leaves
  
  if (las) then

     ! Compute maximum photosynthetic rates,ie, the rate the cohort would have 
     ! if it were at the top of the canopy (used for the mortality function)
     do ipft = 1, n_pft
        
        if(include_pft(ipft) == 1)then
           ! Compute aerodynamic resistance between leaf and canopy air space
           
           !           if (data%icanopy == 3) then
           !              cpatch%rb = (1.0 + 5.5 * pss%lai) /   &
           !                   (0.01 * sqrt(max(0.1,min(pss%ustar,4.)) * const1))
           !           else
           cpatch%rb(ico) = min(25.0 * csite%lai(ipa), 1.0/(  &
                0.003*sqrt(vels/leaf_width(cpatch%pft(ico)))  &
                + 1.03e-5   &
                * ( 1.6e8 * abs(cpatch%veg_temp(ico)-csite%can_temp(ipa))  &
                * leaf_width(cpatch%pft(ico))**3 )**0.25 &
                / leaf_width(cpatch%pft(ico))))

              call lphysiol_full(  &
                   cpatch%veg_temp(ico)-t00  &  ! Vegetation temperature (C)
                   ,csite%can_shv(ipa)*mmdov  &  ! canopy specific humidity (mol/mol)
                   ,csite%can_co2(ipa)*1e-6 &            ! effective CO2 mixing ratio (umol/mol)*(mol/1e6 umol)=(mol/mol)
!                   ,3.6e-4  &                      ! effective CO2 mixing ratio (mol/mol)
                   ,cpatch%par_v(ico)/cpatch%lai(ico)  &             ! absorbed PAR (Ein/m2 leaf/s)
                   ,cpatch%rb(ico)  &                       ! aerodynamic resistance (s/m)
                   ,rhos  &                        ! air density (kg/m3)
                   ,csite%A_o_max(ipft,ipa)   &     ! maximum open photosynthetic rate (umol/m2 leaf/s)
                   
                   ,csite%A_c_max(ipft,ipa)  &      ! maximum closed photosynthetic rate (umol/m2 leaf/s)
                   ,P_op  & ! open stomata resistance for water [s/m]
                   ,P_cl  & ! closed stomata resistance for water [s/m]
                   ,ipft  & ! PFT
                   ,prss  & ! pressure (kg/m/s2)
                   ,leaf_resp  &  ! leaf respiration rate (umol/m2 leaf/s)
                   ,green_leaf_factor(ipft)  &  ! fraction of actual green leaves relative to on-allometry value.
                   ,leaf_aging_factor(ipft)  &
                   ,csite%old_stoma_data_max(ipft,ipa))
              ! type containing the exact stomatal derivatives and met info
              !           else
              !              call lphysiol_sell(cpatch%veg_temp-t00  &
              !                   ,csite%can_shv*mmdov    &
              !                   ,3.6e-4  &
              !                   ,cpatch%par_v/cpatch%lai,cpatch%rb  &
              !                   ,rhos,csite%A_o_max(ipft),csite%A_c_max(ipft)  &
              !                   ,P_op,P_cl,dumarg  &
              !                   ,ipft,prss,data,success_flag,leaf_resp  &
              !                   ,cs%elongation_factor(ipft),cs%gee_phen_delay(ipft))
              !           endif


           endif
        enddo
        
     else

        csite%A_o_max(1:n_pft,ipa) = 0.0
        csite%A_c_max(1:n_pft,ipa) = 0.0

     endif
     

  

     ! cumulative LAI
     cumulative_lai = 0.0
     
     ! LAI/rb, summed over all cohorts. Used in the Euler scheme.
     sum_lai_rbi = 0.0
     
     ! Initialize variables for transpiration calculation
     root_depth_indices(:) = 0
     
     ! Loop over all cohorts
     
     do ico = 1,cpatch%ncohorts
        
        ! Only need to worry about photosyn if radiative transfer has been 
        ! done for this cohort
        
        if(cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > csite%total_snow_depth(ipa))then
        
        ! Aerodynamic resistance [s/m]
        cpatch%rb(ico) = min(25.0 * csite%lai(ipa), 1.0/ &
             (0.003*sqrt(vels*exp(-0.5*cumulative_lai)  &
             /leaf_width(cpatch%pft(ico)))  &
             + 0.5 * 2.06e-5   &
             * ( 1.6e8*abs(cpatch%veg_temp(ico)-csite%can_temp(ipa))  &
             *leaf_width(cpatch%pft(ico))**3 )**0.25 &
             /leaf_width(cpatch%pft(ico))))

!        if (iphoto == 1) then
           call lphysiol_full(  &
                cpatch%veg_temp(ico)-t00  &  ! Vegetation temperature (C)
                ,csite%can_shv(ipa)*mmdov  &  ! canopy specific humidity (mol/mol)
                ,csite%can_co2(ipa)*1e-6 &            ! effective CO2 mixing ratio (umol/mol)*(mol/1e6 umol)=(mol/mol)
               ! ,3.6e-4  &          ! effective CO2 mixing ratio (mol/mol)
                ,cpatch%par_v(ico)/cpatch%lai(ico)  &  ! absorbed PAR (Ein/m2 leaf/s)
                ,cpatch%rb(ico)  &       ! aerodynamic resistance (s/m)
                ,rhos  &      ! air density (kg/m3)
                ,cpatch%A_open(ico) & ! maximum open photosynthetic rate (umol/m2 leaf/s)
                ,cpatch%A_closed(ico) & !maximum closed photosynthetic rate (umol/m2 leaf/s)
                ,cpatch%rsw_open(ico)  & ! open stomata resistance for water [s/m]
                ,cpatch%rsw_closed(ico)  & ! closed stomata resistance for water [s/m]
                ,cpatch%pft(ico)  & ! PFT
                ,prss  & ! pressure (kg/m/s2)
                ,cpatch%leaf_respiration(ico)  &! leaf respiration rate (umol/m2 leaf/s)
                ,green_leaf_factor(cpatch%pft(ico))  & ! fraction of actual green leaves relative to on-allometry value.
                ,leaf_aging_factor(cpatch%pft(ico))  &
                ,cpatch%old_stoma_data(ico)) ! type containing the exact stomatal derivatives and met info
!        else
           ! TRYING SELLERS AND COLLATZ PHOTOSYNTHESIS... 
           ! ---------------------------------------------------------
!           call lphysiol_sell( &
!                cpatch%veg_temp-t00,                          &
!                csite%can_shv*mmdov,   &
!                3.6e-4,                               &
!                cpatch%par_v/cpatch%lai,                    &
!                cpatch%rb,                               &
!                rhos,                               &
!                cpatch%A_op,                             &
!                cpatch%A_cl,                             &
!                cpatch%P_op,                             &
!                cpatch%P_cl,                             &
!                dumarg,                               &
!                cpatch%pft,                               &
!                prss,                              &
!                data,                                 &
!                success_flag,                         &
!                leaf_resp,                            &
!                elongation_factor(cpatch%pft),         &
!                gee_phen_delay(cpatch%pft))
!        endif
        sum_lai_rbi = sum_lai_rbi + cpatch%lai(ico) / cpatch%rb(ico)

        ! Leaf respiration
        cpatch%leaf_respiration(ico) = cpatch%leaf_respiration(ico)  * cpatch%lai(ico)
        cpatch%mean_leaf_resp(ico)   = cpatch%mean_leaf_resp(ico)    + cpatch%leaf_respiration(ico)
        cpatch%dmean_leaf_resp(ico)  = cpatch%dmean_leaf_resp(ico)   + cpatch%leaf_respiration(ico)

        ! Demand for water
        water_demand = cpatch%Psi_open(ico) ! kg/m2/s; Psi_open is from last time step

        ! Supply of water
        water_supply = water_conductance(cpatch%pft(ico)) *   &
             available_liquid_water(cpatch%krdepth(ico)) * 1.0e-3 *   &
             q(cpatch%pft(ico)) * cpatch%balive(ico) / (1.0 + q(cpatch%pft(ico)) + cpatch%hite(ico) *   &
             qsw(cpatch%pft(ico))) * cpatch%nplant(ico)

        root_depth_indices(cpatch%krdepth(ico)) = 1

        ! Weighting between open/closed stomata
        cpatch%fsw(ico) = water_supply / max(1.0e-30,water_supply + water_demand)

        ! Account for nitrogen limitation
        cpatch%fs_open(ico) = cpatch%fsw(ico) * cpatch%fsn(ico)

        ! Photorespiration can become important at high temperatures.  If so,
        ! close down the stomata.
        if(cpatch%A_open(ico) < cpatch%A_closed(ico)) cpatch%fs_open(ico) = 0.0

        ! Net stomatal resistance
        cpatch%stomatal_resistance(ico) = 1.0 / (cpatch%fs_open(ico) / cpatch%rsw_open(ico) +   &
             (1.0 - cpatch%fs_open(ico)) / cpatch%rsw_closed(ico))

        ! GPP, averaged over frqstate
        cpatch%gpp(ico) = cpatch%lai(ico) * (cpatch%fs_open(ico) * cpatch%A_open(ico) &
             + (1.0 - cpatch%fs_open(ico)) *   &
             cpatch%A_closed(ico)) + cpatch%leaf_respiration(ico)

        cpatch%mean_gpp(ico) = cpatch%mean_gpp(ico) + cpatch%gpp(ico)

        ! GPP, summed over 1 day [umol/m2]
        cpatch%dmean_gpp(ico) = cpatch%dmean_gpp(ico) + cpatch%gpp(ico)

        ! Potential GPP if no N limitation
        cpatch%dmean_gpp_pot(ico) = cpatch%dmean_gpp_pot(ico) + cpatch%lai(ico) * &
             (cpatch%fsw(ico) * cpatch%A_open(ico) +  &
             (1.0 - cpatch%fsw(ico)) * cpatch%A_closed(ico)) + cpatch%leaf_respiration(ico)

        ! Maximum GPP if at the top of the canopy
        cpatch%dmean_gpp_max(ico) = cpatch%dmean_gpp_max(ico) + cpatch%lai(ico) * (cpatch%fs_open(ico) *   &
             csite%A_o_max(cpatch%pft(ico),ipa) + (1.0 - cpatch%fs_open(ico)) *   &
             csite%A_c_max(cpatch%pft(ico),ipa)) + cpatch%leaf_respiration(ico)

        ! Update cumulative LAI
        cumulative_lai = cumulative_lai + cpatch%lai(ico)

     else

        ! This is the case where a cohort does not have leaves or is 
        ! snow-covered
        cpatch%A_open(ico)              = 0.0
        cpatch%A_closed(ico)            = 0.0
        cpatch%Psi_open(ico)            = 0.0
        cpatch%Psi_closed(ico)          = 0.0
        cpatch%rsw_open(ico)            = 0.0
        cpatch%rsw_closed(ico)          = 0.0
        cpatch%rb(ico)                  = 0.0
        cpatch%stomatal_resistance(ico) = 0.0
        cpatch%gpp(ico)                 = 0.0
        cpatch%leaf_respiration(ico)    = 0.0
     endif

  enddo

  ! For plants of a given rooting depth, determine soil level from which 
  ! transpired water is to be extracted.
  ed_ktrans(:) = 0
  do k1 = lsl, nzg
     swp = -1.e10
     if(root_depth_indices(k1) == 1)then
        do k2 = k1, nzg
           nts = ntext_soil(k2)
           slpotv = soil(nts)%slpots * (soil(nts)%slmsts / soil_water(k2)) ** soil(nts)%slbs
           ! Multiply by liquid fraction (ice is unavailable for transpiration)
           slpotv = slpotv * soil_fracliq(k2)
           ! Find layer in root zone with highest slpotv AND soil_water 
           ! above minimum soilcp.  Set ktrans to this layer.
           if (slpotv > swp .and. soil_water(k2) > soil(nts)%soilcp) then
              swp = slpotv
              ed_ktrans(k1) = k2
           endif
        enddo
     endif
  enddo

  return
end subroutine canopy_photosynthesis_ar
