module growth_balive_ar


contains

subroutine dbalive_dt_ar(cgrid, tfact)

  ! Do not change the order of operations below unless you know what you
  ! are doing.  Changing the order can affect the C/N budgets.

  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
  use pft_coms, only: q, qsw, plant_N_supply_scale, c2n_storage, &
       growth_resp_factor, storage_turnover_rate
  use physiology_coms, only: N_plant_lim
  use grid_coms, only: nzg

  implicit none

  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy,isi,ipa,ico

  real :: salloc
  real :: salloci
  real :: bl
  real :: br
  real, intent(in) :: tfact
  real :: daily_C_gain
  real :: carbon_balance
  real :: carbon_balance_pot
  real :: carbon_balance_max
  real :: balive_in
  real :: nitrogen_supply
  real, external :: mortality_rates_ar
  real :: dndt
  real :: nitrogen_uptake
  real :: N_uptake_pot

  do ipy = 1,cgrid%npolygons
     
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)
        
        do ipa = 1,csite%npatches

           cpatch => csite%patch(ipa)
        
     
           ! Reset averaged variables
           csite%total_plant_nitrogen_uptake(ipa) = 0.0

           ! Loop over cohorts
           do ico = 1,cpatch%ncohorts

              ! Initialize cohort nitrogen uptake
              nitrogen_uptake = 0.0
              N_uptake_pot = 0.0
              
              ! Set allocation factors
              salloc = 1.0 + qsw(cpatch%pft(ico)) * cpatch%hite(ico) + q(cpatch%pft(ico))
              salloci = 1.0 / salloc
              
              ! Subtract preceding day's storage respiration.
              cpatch%bstorage(ico) = cpatch%bstorage(ico) - cpatch%storage_respiration(ico)
!              cpatch%cb(13,ico) = cpatch%cb(13,ico) - cpatch%storage_respiration(ico)
!              cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%storage_respiration(ico)
              
              ! When you lose storage carbon, allow the associated nitrogen 
              ! to go to litter in order to maintain prescribed C2N ratio.
              csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%storage_respiration(ico) /   &
                   c2n_storage * cpatch%nplant(ico)

              ! If plants are off allometry, move carbon from bstorage
              ! to balive
              call transfer_C_from_storage_ar(cpatch,ico, salloc, nitrogen_uptake,   &
                   N_uptake_pot)
              
              ! Calculate leaf, fine root biomass
              if(cpatch%phenology_status(ico) /= 2)then
                 bl = cpoly%green_leaf_factor(cpatch%pft(ico),isi) * cpatch%balive(ico) * salloci
              else
                 bl = 0.0
              endif
              br = q(cpatch%pft(ico)) * cpatch%balive(ico) * salloci 

              ! Compute maintenance costs and growth/storage/vleaf respiration
              ! for the coming day.
              
              call plant_maintenance_and_resp_ar(cpatch,ico, br, bl, tfact,   &
                   daily_C_gain, csite%avg_daily_temp(ipa))
              
              ! Subtract maintenance costs from balive.
              cpatch%balive(ico) = cpatch%balive(ico) - cpatch%maintenance_costs(ico)
              cpatch%cb(13,ico) = cpatch%cb(13,ico) - cpatch%maintenance_costs(ico)
              cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%maintenance_costs(ico)

              ! Calculate actual, potential and maximum carbon balances
              call plant_carbon_balances_ar(cpatch,ico, daily_C_gain, carbon_balance,  &
                   carbon_balance_pot, carbon_balance_max)

              ! Compute respiration rates for coming day [kgC/plant/day]
              cpatch%growth_respiration(ico) = max(0.0, daily_C_gain *   &
                   growth_resp_factor(cpatch%pft(ico)))
              cpatch%storage_respiration(ico) = cpatch%bstorage(ico) *   &
                   storage_turnover_rate(cpatch%pft(ico)) * tfact
              cpatch%vleaf_respiration(ico) = (1.0 - cpoly%green_leaf_factor(cpatch%pft(ico),isi))  &
                   / (1.0 + q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) * cpatch%hite(ico))  &
                   * cpatch%balive(ico) * storage_turnover_rate(cpatch%pft(ico)) * tfact


              ! Allocate plant carbon balance to balive and bstorage
              balive_in = cpatch%balive(ico)
              call alloc_plant_c_balance_ar(cpatch,ico, salloc, salloci,   &
                   carbon_balance, nitrogen_uptake,   &
                   cpoly%green_leaf_factor(cpatch%pft(ico),isi), csite%fsn_in(ipa))
             
              ! Do a shadow calculation to see what would have happened if 
              ! stomata were open.  This is used to calculate potential 
              ! nitrogen uptake, N_uptake_pot.

              if(N_plant_lim == 1)call potential_N_uptake_ar(cpatch,ico, salloc,   &
                   salloci, balive_in, carbon_balance_pot, N_uptake_pot,   &
                   cpoly%green_leaf_factor(cpatch%pft(ico),isi))
              
              !  Increment the [kgN/m2] taken up during previous day.
              csite%total_plant_nitrogen_uptake(ipa) =   &
                   csite%total_plant_nitrogen_uptake(ipa)   &
                   + nitrogen_uptake * cpatch%nplant(ico)
              
              ! Calculate plant N limitation factor
              if(n_plant_lim == 0 .or. N_uptake_pot <= 0.0)then
                 cpatch%fsn(ico) = 1.0
              else
                 br = q(cpatch%pft(ico)) * cpatch%balive(ico) * salloci 
                 nitrogen_supply = plant_N_supply_scale * br   &
                      * csite%mineralized_soil_N(ipa)
                 cpatch%fsn(ico) = nitrogen_supply / (nitrogen_supply + N_uptake_pot)
              endif
              
              ! Do mortality --- note that only frost mortality changes daily.
              dndt = - mortality_rates_ar(cpatch,ico, csite%avg_daily_temp(ipa)) *   &
                   cpatch%nplant(ico) * tfact
              
              ! Update monthly mortality rate [plants/m2/month]
              cpatch%monthly_dndt(ico) = cpatch%monthly_dndt(ico) + dndt

           enddo
           
           ! Update litter
           call litter_ar(csite,ipa)
           
           ! Recompute patch LAI
           call update_patch_derived_props_ar(csite,cpoly%lsl(isi),   &
                cpoly%met(isi)%rhos,ipa)

           !reset average daily temperature
           csite%avg_daily_temp(ipa) = 0.0 
           
        enddo

     enddo

  enddo


  return
end subroutine dbalive_dt_ar

!====================================================================

subroutine transfer_C_from_storage_ar(cpatch,ico, salloc, nitrogen_uptake, N_uptake_pot)

  use ed_state_vars,only:patchtype

  use pft_coms, only: c2n_leaf, c2n_storage, c2n_stem
  use decomp_coms, only: f_labile

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: salloc
  real, intent(inout) :: nitrogen_uptake
  real, intent(inout) :: N_uptake_pot
  real :: off_allometry_cb
  real, external :: dbh2bl
  real :: increment

  ! Only do the transfer there are supposed to be leaves

  if(cpatch%phenology_status(ico) < 2)then
     
     ! If plants have storage, transfer it to balive
     off_allometry_cb = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * salloc - cpatch%balive(ico)
     increment = max(0.0,min( max(0.0, off_allometry_cb), cpatch%bstorage(ico)))
     cpatch%balive(ico) = cpatch%balive(ico) + increment
     cpatch%bstorage(ico) = cpatch%bstorage(ico) - increment

     ! N uptake is required since c2n_leaf < c2n_storage.
     ! Units are kgN/plant/day.
     nitrogen_uptake = increment * ( f_labile (cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) +   &
          (1.0 - f_labile(cpatch%pft(ico))) / c2n_stem - 1.0 / c2n_storage)
     N_uptake_pot = nitrogen_uptake
     
  endif

  return
end subroutine transfer_C_from_storage_ar

!====================================================================

subroutine plant_maintenance_and_resp_ar(cpatch,ico, br, bl, tfact, daily_C_gain, tempk)
  
  use ed_state_vars,only:patchtype
  use pft_coms, only: phenology, root_turnover_rate, leaf_turnover_rate
  use consts_coms, only: umol_2_kgC,day_sec
  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: br
  real, intent(in) :: bl
  real, intent(in) :: tfact
  real, intent(in) :: tempk
  real, intent(out) :: daily_C_gain
  real :: maintenance_temp_dep

  ! Get the temperature dependence
  if(phenology(cpatch%pft(ico)) /= 1)then
     maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
  else
     maintenance_temp_dep = 1.0
  endif

  ! Calculate maintenance demand (kgC/plant/year)
  
  cpatch%maintenance_costs(ico) = (root_turnover_rate(cpatch%pft(ico)) * br +  &
       leaf_turnover_rate(cpatch%pft(ico)) * bl) * maintenance_temp_dep

  ! Convert units of maintenance to [kgC/plant/day]
  
  cpatch%maintenance_costs(ico) = cpatch%maintenance_costs(ico) * tfact
  
        
  ! Compute daily C uptake [kgC/plant/day]

  if(cpatch%nplant(ico) .gt. tiny(1.0)) then
     daily_C_gain = umol_2_kgC * day_sec * (cpatch%dmean_gpp(ico) -   &
          cpatch%dmean_leaf_resp(ico) - cpatch%dmean_root_resp(ico)) / cpatch%nplant(ico)
  else
     daily_C_gain = 0.0
  endif

  return
end subroutine plant_maintenance_and_resp_ar

!===================================================================

subroutine plant_carbon_balances_ar(cpatch,ico, daily_C_gain, carbon_balance,  &
     carbon_balance_pot, carbon_balance_max)
 
  use ed_state_vars,only:patchtype
  use pft_coms, only: growth_resp_factor
  use consts_coms, only: umol_2_kgC,day_sec

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: daily_C_gain
  real, intent(out) :: carbon_balance
  real, intent(out) :: carbon_balance_pot
  real, intent(out) :: carbon_balance_max
  real :: daily_C_gain_pot
  real :: daily_C_gain_max
  real :: growth_respiration_pot
  real :: growth_respiration_max

  ! Calculate actual daily carbon balance: kgC/plant/day.

  carbon_balance = daily_C_gain - cpatch%growth_respiration(ico) - cpatch%vleaf_respiration(ico)

  if(cpatch%nplant(ico) .gt. tiny(1.0)) then

     ! Calculate potential carbon balance (used for nitrogen 
     ! demand function).  [kgC/plant/day]
     
     daily_C_gain_pot = umol_2_kgC * day_sec * (cpatch%dmean_gpp_pot(ico) -   &
          cpatch%dmean_leaf_resp(ico) - cpatch%dmean_root_resp(ico)) / cpatch%nplant(ico)
     growth_respiration_pot = max(0.0, daily_C_gain_pot *   &
          growth_resp_factor(cpatch%pft(ico)))
     carbon_balance_pot = daily_C_gain_pot - growth_respiration_pot -   &
          cpatch%vleaf_respiration(ico)
     
     ! Calculate maximum carbon balance (used for mortality) 
     daily_C_gain_max = umol_2_kgC * day_sec * (cpatch%dmean_gpp_max(ico) -   &
          cpatch%dmean_leaf_resp(ico) - cpatch%dmean_root_resp(ico)) / cpatch%nplant(ico)
     growth_respiration_max = max(0.0, daily_C_gain_max *   &
          growth_resp_factor(cpatch%pft(ico)))
     carbon_balance_max = daily_C_gain_max - growth_respiration_max -   &
          cpatch%vleaf_respiration(ico)
     
  else
     carbon_balance_max = 0.0
     carbon_balance_pot = 0.0
  end if
  ! Carbon balances for mortality
  cpatch%cb(13,ico) = cpatch%cb(13,ico) + carbon_balance
  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) + carbon_balance_max

  return
end subroutine plant_carbon_balances_ar

!====================================================================

subroutine alloc_plant_c_balance_ar(cpatch,ico, salloc, salloci, carbon_balance,   &
     nitrogen_uptake, green_leaf_factor, fsn_in)

  use ed_state_vars,only:patchtype
  use pft_coms, only: c2n_storage, c2n_leaf, sla, c2n_stem
  use decomp_coms, only: f_labile
  use therm_lib,only : update_veg_energy_cweh

  implicit none
  
  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: salloc
  real, intent(in) :: salloci
  real, intent(in) :: carbon_balance
  real, intent(inout) :: nitrogen_uptake
  real, intent(in) :: green_leaf_factor
  real, intent(inout) :: fsn_in
  real :: bl_max
  real :: bl_pot
  real, external :: dbh2bl
  real :: increment

  if(cpatch%phenology_status(ico) == 0 .and. carbon_balance > 0.0 )then

     ! Simply update monthly carbon gain.  This will be 
     ! used for structural growth at the end of the month.
     cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
     nitrogen_uptake = nitrogen_uptake +  &
          carbon_balance / c2n_storage
     cpatch%bleaf(ico) = cpatch%balive(ico) * salloci * green_leaf_factor

  else

     ! are there leaves?
     if(cpatch%phenology_status(ico) < 2)then

        bl_max = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * green_leaf_factor
        bl_pot = green_leaf_factor  &
             * (cpatch%balive(ico) + carbon_balance) * salloci
        ! will this increment take us over the limit?
        if(bl_pot > bl_max)then

           ! if so, put remainder in storage
           increment = carbon_balance   &
                - (dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * salloc - cpatch%balive(ico))
           cpatch%bstorage(ico) = cpatch%bstorage(ico) + increment
           nitrogen_uptake = nitrogen_uptake +  &
                increment / c2n_storage
           increment = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * salloc - cpatch%balive(ico)
           cpatch%balive(ico) = cpatch%balive(ico) + increment
           nitrogen_uptake = nitrogen_uptake +  &
                increment * (f_labile(cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) +  &
                (1.0 - f_labile(cpatch%pft(ico))) / c2n_stem)
           cpatch%bleaf(ico) = bl_max

           cpatch%phenology_status(ico) = 0

        else
           

           ! it will not exceed limit, so just add to balive
           cpatch%balive(ico) = max(0.0,cpatch%balive(ico) + carbon_balance)
           cpatch%phenology_status(ico) = 1
           cpatch%bleaf(ico) = cpatch%balive(ico) * salloci *   &
                green_leaf_factor


           if(carbon_balance < 0.0)then
              fsn_in = fsn_in - carbon_balance *  &
                   (f_labile(cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) + (1.0 -  &
                   f_labile(cpatch%pft(ico))) / c2n_stem ) * cpatch%nplant(ico)
           else
              nitrogen_uptake = nitrogen_uptake +  &
                   carbon_balance * (f_labile(cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) +  &
                   (1.0 - f_labile(cpatch%pft(ico))) / c2n_stem)
           endif

        endif

     else

        ! in this case, carbon balance in negative so just subtract
        cpatch%balive(ico) = max(0.0,cpatch%balive(ico) + carbon_balance)
        cpatch%bleaf(ico) = 0.0
        cpatch%lai(ico) = 0.0
        fsn_in = fsn_in - carbon_balance *  &
             (f_labile(cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) + (1.0 -   &
             f_labile(cpatch%pft(ico))) / c2n_stem) * cpatch%nplant(ico)
        
     endif
  endif
  
  cpatch%lai(ico) = cpatch%bleaf(ico) * cpatch%nplant(ico) * sla(cpatch%pft(ico))
  
  ! It is likely that the leaf biomass has changed, therefore, update vegetation
  ! energy and heat capacity.
  
  call update_veg_energy_cweh(cpatch,ico)
  

  return
end subroutine alloc_plant_c_balance_ar

!====================================================================

subroutine potential_N_uptake_ar(cpatch,ico, salloc, salloci, balive_in,   &
     carbon_balance_pot, N_uptake_pot, green_leaf_factor)
 
  use ed_state_vars,only:patchtype
  use pft_coms, only: c2n_storage, c2n_leaf, c2n_stem
  use decomp_coms, only: f_labile

  implicit none
  type(patchtype),target :: cpatch
  integer :: ico
  real, intent(in) :: salloc
  real, intent(in) :: salloci
  real, intent(in) :: balive_in
  real, intent(in) :: carbon_balance_pot
  real, intent(inout) :: N_uptake_pot
  real, intent(in) :: green_leaf_factor
  real :: bl_max
  real :: bl_pot
  real, external :: dbh2bl
  real :: increment

  if(cpatch%phenology_status(ico) == 0 .and. carbon_balance_pot > 0.0 )then

     ! Positive carbon balance with plants fully flushed
     N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage

  else

     if(cpatch%phenology_status(ico) < 2)then

        ! There are at least some leaves
        bl_max = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * green_leaf_factor
        bl_pot = green_leaf_factor  &
             * (balive_in + carbon_balance_pot) * salloci

        if(bl_pot > bl_max)then

           ! this increment took us over the limit, so remainder is 
           ! put in storage
           increment = carbon_balance_pot   &
                - (dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * salloc - balive_in)
           N_uptake_pot = N_uptake_pot + increment / c2n_storage
           increment = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * salloc - balive_in
           N_uptake_pot = N_uptake_pot + increment * (f_labile(cpatch%pft(ico)) /   &
                c2n_leaf(cpatch%pft(ico)) + (1.0 - f_labile(cpatch%pft(ico))) / c2n_stem)
        else

           ! this increment did not exceed the limit.

           if(carbon_balance_pot > 0.0)then

              ! There is uptake if carbon_balance_pot is positive

              N_uptake_pot = N_uptake_pot + carbon_balance_pot *   &
                   (f_labile(cpatch%pft(ico)) / c2n_leaf(cpatch%pft(ico)) + (1.0 -   &
                   f_labile(cpatch%pft(ico))) / c2n_stem)
           endif

        endif

     endif

  endif

  return
end subroutine potential_N_uptake_ar

!====================================================================

subroutine litter_ar(csite,ipa)

  use ed_state_vars,only:patchtype,sitetype
  use pft_coms, only: c2n_leaf, c2n_stem, l2n_stem
  use decomp_coms, only: f_labile

  implicit none
  type(sitetype),target   :: csite
  type(patchtype),pointer :: cpatch
  integer :: ico,ipa
  real :: plant_litter
  real :: plant_litter_f
  real :: plant_litter_s

  cpatch => csite%patch(ipa)

  ! Add fine root and leaf turnover to the litter

  ! Loop over cohorts

  do ico=1,cpatch%ncohorts

     plant_litter = cpatch%maintenance_costs(ico) * cpatch%nplant(ico)
     plant_litter_f = plant_litter * f_labile(cpatch%pft(ico))
     plant_litter_s = plant_litter - plant_litter_f

     csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f
     csite%fsn_in(ipa) = csite%fsn_in(ipa) + plant_litter_f / c2n_leaf(cpatch%pft(ico))

     csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
     csite%ssl_in(ipa) = csite%ssl_in(ipa) + plant_litter_s * l2n_stem / c2n_stem

  enddo

  return
end subroutine litter_ar

end module growth_balive_ar
