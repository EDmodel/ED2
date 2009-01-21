!====================================================================
! ============================================

subroutine init_ed_cohort_vars_array(cpatch,ico, lsl)
  
  use ed_state_vars,only : patchtype
  use allometry, only: calc_root_depth, assign_root_depth

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real :: root_depth
  integer, intent(in) :: lsl

  cpatch%mean_gpp(ico) = 0.0
  cpatch%mean_leaf_resp(ico) = 0.0
  cpatch%mean_root_resp(ico) = 0.0
  
  cpatch%dmean_leaf_resp(ico) = 0.0
  cpatch%dmean_root_resp(ico) = 0.0
  cpatch%dmean_gpp(ico) = 0.0
  cpatch%dmean_gpp_pot(ico) = 0.0
  cpatch%dmean_gpp_max(ico) = 0.0
  
  cpatch%leaf_respiration(ico) = 0.0
  cpatch%root_respiration(ico) = 0.0
  cpatch%growth_respiration(ico) = 0.0
  cpatch%storage_respiration(ico) = 0.0
  cpatch%vleaf_respiration(ico) = 0.0

  cpatch%fsw(ico) = 1.0
  cpatch%fsn(ico) = 1.0
  
  !----- This variable would never be assigned for low LAI cohorts 
  cpatch%fs_open(ico) = cpatch%fsw(ico)*cpatch%fsn(ico)

  cpatch%monthly_dndt(ico) = 0.0

  cpatch%par_v(ico)            = 0.0
  cpatch%par_v_beam(ico)       = 0.0
  cpatch%par_v_diffuse(ico)    = 0.0
  cpatch%rshort_v(ico)         = 0.0
  cpatch%rshort_v_beam(ico)    = 0.0
  cpatch%rshort_v_diffuse(ico) = 0.0
  cpatch%rlong_v(ico)          = 0.0
  cpatch%rlong_v_surf(ico)     = 0.0
  cpatch%rlong_v_incid(ico)    = 0.0
       
  cpatch%rb(ico)               = 0.0
  cpatch%A_open(ico)           = 0.0
  cpatch%A_closed(ico)         = 0.0
  cpatch%Psi_closed(ico)       = 0.0
  cpatch%rsw_open(ico)         = 0.0
  cpatch%rsw_closed(ico)       = 0.0
      
      
  cpatch%stomatal_resistance(ico) = 0.0
  cpatch%maintenance_costs(ico)   = 0.0
  cpatch%paw_avg10d(ico)          = 0.0

  ! From the ed_state_vars comment, these variables are now deprecated, commenting
  ! them so it will compile...
  ! cpatch%co_srad_h(ico)        = 1.0
  ! cpatch%co_lrad_h(ico)        = 1.0
  ! cpatch%co_sens_h(ico)        = 1.0
  ! cpatch%co_evap_h(ico)        = 1.0
  ! cpatch%co_liqr_h(ico)        = 1.0

  ! cpatch%co_srad_h(ico)  =   cpatch%co_srad_h(ico)  -      1.0
  ! cpatch%co_lrad_h(ico)  =   cpatch%co_lrad_h(ico)  -      1.0
  ! cpatch%co_sens_h(ico)  =   cpatch%co_sens_h(ico)  -      1.0
  ! cpatch%co_evap_h(ico)  =   cpatch%co_evap_h(ico)  -      1.0
  ! cpatch%co_liqr_h(ico)  =   cpatch%co_liqr_h(ico)  -      1.0


  cpatch%Psi_open(ico) = 0.0

  cpatch%cb(1:13,ico) = 1.0
  cpatch%cb_max(1:13,ico) = 1.0
  cpatch%cbr_bar(ico) = 1.0

  cpatch%old_stoma_data(ico)%recalc = 1
  root_depth = calc_root_depth(cpatch%hite(ico),cpatch%dbh(ico), cpatch%pft(ico))
  cpatch%krdepth(ico) = assign_root_depth(root_depth, lsl)
  
  cpatch%first_census(ico) = 1
  cpatch%new_recruit_flag(ico) = 0
  cpatch%bseeds(ico) = 0.0

  return
end subroutine init_ed_cohort_vars_array

! ==========================================

subroutine init_ed_patch_vars_array(csite,ip1,ip2)
  
  use ed_state_vars,only:sitetype
  use max_dims, only: n_pft
!  use fuse_fiss_utils, only: count_cohorts
  use grid_coms,     only: nzs, nzg

  implicit none

  type(sitetype), target :: csite
  integer :: ipft,ncohorts,ipa
  integer :: ip1,ip2
  
  do ipa = ip1,ip2
     do ipft = 1,n_pft
        csite%old_stoma_data_max(ipft,ipa)%recalc = 1
     enddo
  enddo


  ! Initialize sfcwater state variables
  csite%sfcwater_mass(1:nzs,ip1:ip2) = 0.0
  csite%sfcwater_energy(1:nzs,ip1:ip2) = 0.0
  csite%sfcwater_depth(1:nzs,ip1:ip2) = 0.0
  csite%rshort_s(:,ip1:ip2) = 0.0
  csite%rshort_s_beam(:,ip1:ip2) = 0.0
  csite%rshort_s_diffuse(:,ip1:ip2) = 0.0
  csite%nlev_sfcwater(ip1:ip2) = 0
  
  csite%rlong_s(ip1:ip2) = 0.0
  
  csite%avg_daily_temp(ip1:ip2) = 0.0

  csite%mean_rh(ip1:ip2) = 0.0
  csite%mean_nep(ip1:ip2) = 0.0

  csite%mean_runoff(ip1:ip2) = 0.0
  csite%mean_wflux(ip1:ip2) = 0.0
  csite%mean_latflux(ip1:ip2) = 0.0
  csite%mean_qrunoff(ip1:ip2) = 0.0
  csite%mean_hflux(ip1:ip2) = 0.0

  csite%dmean_A_decomp(ip1:ip2) = 0.0
  csite%dmean_Af_decomp(ip1:ip2) = 0.0

  csite%repro(1:n_pft,ip1:ip2) = 0.0

  csite%htry(ip1:ip2) = 1.0

  csite%can_co2(ip1:ip2) = 370.0

  csite%wbudget_loss2atm(ip1:ip2) = 0.0
  csite%wbudget_loss2runoff(ip1:ip2) = 0.0
  csite%co2budget_loss2atm(ip1:ip2) = 0.0
  csite%ebudget_loss2atm(ip1:ip2) = 0.0
  csite%ebudget_latent(ip1:ip2) = 0.0

  csite%wbudget_precipgain(ip1:ip2) = 0.0
  csite%ebudget_precipgain(ip1:ip2) = 0.0
  csite%ebudget_netrad(ip1:ip2) = 0.0
  csite%co2budget_gpp(ip1:ip2) = 0.0
  csite%co2budget_gpp_dbh(:,ip1:ip2) = 0.0
  csite%co2budget_plresp(ip1:ip2) = 0.0
  csite%co2budget_rh(ip1:ip2) = 0.0

  !----------------------------------------------------------------------------------------!
  !    These variables need to be initialized here otherwise it will fail when new patches !
  ! are created.                                                                           !
  !----------------------------------------------------------------------------------------!
  csite%avg_carbon_ac(ip1:ip2)    = 0.0
  csite%avg_vapor_vc(ip1:ip2)     = 0.0
  csite%avg_dew_cg(ip1:ip2)       = 0.0
  csite%avg_vapor_gc(ip1:ip2)     = 0.0
  csite%avg_wshed_vg(ip1:ip2)     = 0.0
  csite%avg_vapor_ac(ip1:ip2)     = 0.0
  csite%avg_transp(ip1:ip2)       = 0.0
  csite%avg_evap(ip1:ip2)         = 0.0
  csite%avg_smoist_gg(:,ip1:ip2)    = 0.0
  csite%avg_smoist_gc(:,ip1:ip2)    = 0.0
  csite%avg_runoff(ip1:ip2)       = 0.0
  csite%avg_sensible_vc(ip1:ip2)  = 0.0
  csite%avg_sensible_2cas(ip1:ip2)= 0.0
  csite%avg_qwshed_vg(ip1:ip2)    = 0.0
  csite%avg_sensible_gc(ip1:ip2)  = 0.0
  csite%avg_sensible_ac(ip1:ip2)  = 0.0
  csite%avg_sensible_tot(ip1:ip2) = 0.0
  csite%avg_sensible_gg(:,ip1:ip2)  = 0.0
  csite%avg_runoff_heat(ip1:ip2)  = 0.0
  csite%aux(ip1:ip2)              = 0.0
  csite%aux_s(:,ip1:ip2)            = 0.0
  
  csite%avg_heatstor_veg(ip1:ip2) = 0.0  !SHOULD THIS BE ZERO'D ALSO?

  csite%rshort_g(ip1:ip2) = 0.0
  csite%rshort_g_beam(ip1:ip2) = 0.0
  csite%rshort_g_diffuse(ip1:ip2) = 0.0
  csite%rlong_g(ip1:ip2) = 0.0
  csite%rlong_g_surf(ip1:ip2) = 0.0
  csite%rlong_g_incid(ip1:ip2) = 0.0
  csite%rlong_s(ip1:ip2) = 0.0
  csite%rlong_s_surf(ip1:ip2) = 0.0
  csite%rlong_s_incid(ip1:ip2) = 0.0
  csite%albedt(ip1:ip2) = 0.0
  csite%albedo_beam(ip1:ip2) = 0.0
  csite%albedo_diffuse(ip1:ip2) = 0.0
  csite%rlongup(ip1:ip2) = 0.0
  csite%rlong_albedo(ip1:ip2) = 0.0

  csite%fsc_in(ip1:ip2)                      = 0.0
  csite%ssc_in(ip1:ip2)                      = 0.0
  csite%ssl_in(ip1:ip2)                      = 0.0
  csite%fsn_in(ip1:ip2)                      = 0.0
  csite%total_plant_nitrogen_uptake(ip1:ip2) = 0.0

  ncohorts = 0
  do ipa=1,csite%npatches
     ncohorts = ncohorts + csite%patch(ipa)%ncohorts
  enddo

  csite%cohort_count = ncohorts

  return
end subroutine init_ed_patch_vars_array

!======================================================================


subroutine init_ed_site_vars_array(cpoly, lat)

  use ed_state_vars,only:polygontype
  use max_dims, only: n_pft, n_dbh, n_dist_types 
  use grid_coms,     only: nzs, nzg

  implicit none

  type(polygontype),target :: cpoly
  real, intent(in) :: lat
  integer, external :: julday

  cpoly%basal_area(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%basal_area_growth(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%basal_area_mort(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%basal_area_cut(1:n_pft, 1:n_dbh,:) = 0.0
!  cpoly%basal_area_recruit(1:n_pft, 1:n_dbh,:) = 0.0

  cpoly%agb(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%agb_growth(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%agb_mort(1:n_pft, 1:n_dbh,:) = 0.0
  cpoly%agb_cut(1:n_pft, 1:n_dbh,:) = 0.0
!  cpoly%agb_recruit(1:n_pft, 1:n_dbh,:) = 0.0

  cpoly%green_leaf_factor(1:n_pft,:) = 1.0
  cpoly%leaf_aging_factor(1:n_pft,:) = 1.0

  cpoly%min_monthly_temp(:) = 0.0

 ! cpoly%mm_gpp(:) = 0.0
 ! cpoly%mm_plresp(:) = 0.0
 ! cpoly%mm_rh(:) = 0.0
 ! cpoly%mm_nep(:) = 0.0

!  cpoly%mean_precip(:) = 0.0
!  cpoly%mean_qprecip(:) = 0.0
!  cpoly%mean_netrad(:) = 0.0

  cpoly%lambda_fire(1:12,:) = 0.0
  
  cpoly%disturbance_memory(1:n_dist_types, 1:n_dist_types,:) = 0.0

  cpoly%agri_stocking_density(:) = 10.0

!KIM - using latitude for c3/c4 crops may not always work.  
!    - e.g., a lot of corn (C4 crop) grows in the US Midwest.
!    - agri_stocking_pft would better be sepcified according to the need.
!    - possibly in ED2IN?
!    - anyway, this part should be more elaborate for the case 
!    - that we have different crops/pastures.
  if(lat > 40.0)then
     cpoly%agri_stocking_pft(:) = 12
     cpoly%plantation_stocking_pft(:) = 6
  else
     cpoly%agri_stocking_pft(:) = 14
     cpoly%plantation_stocking_pft(:) = 7
  endif
  cpoly%plantation_stocking_density(:) = 4.0

  cpoly%primary_harvest_memory(:) = 0.0
  cpoly%secondary_harvest_memory(:) = 0.0

!  cpoly%soil_water(1:nzg,:)     = 0.0
!  cpoly%soil_tempk(1:nzg,:)     = 0.0
!  cpoly%soil_energy(1:nzg,:)    = 0.0
!  cpoly%sfcwater_tempk(1:nzs,:) = 0.0
!  cpoly%sfcwater_mass(1:nzs,:)  = 0.0
!  cpoly%sfcwater_depth(1:nzs,:) = 0.0

!  cpoly%avg_smoist_gc(1:nzg)  = 0.0
!  cpoly%avg_smoist_gg(1:nzg)  = 0.0
!  cpoly%avg_sensible_gg(1:nzg)= 0.0
!  cpoly%aux_s(1:nzg)          = 0.0
  
  return
end subroutine init_ed_site_vars_array

!======================================================================
subroutine init_ed_poly_vars_array(cgrid)
  
  use ed_state_vars,only:edtype
  
  implicit none
  
  type(edtype),target :: cgrid
  integer :: ipy
  
  real :: soil_C
  real :: soil_N
  real :: veg_C
  real :: veg_N
  
  

  do ipy = 1,cgrid%npolygons
     !Moved inside the loop for the cases in which npolygons is 0
     cgrid%mean_precip(ipy)  = 0.0
     cgrid%mean_qprecip(ipy) = 0.0
     cgrid%mean_netrad(ipy)  = 0.0
     call compute_C_and_N_storage(cgrid,ipy,soil_C, soil_N, veg_C, veg_N)
     cgrid%cbudget_initialstorage(ipy) = soil_C + veg_C
     cgrid%nbudget_initialstorage(ipy) = soil_N + veg_N
     cgrid%cbudget_nep(ipy) = 0.0
  enddo

  return
end subroutine init_ed_poly_vars_array


!===================================================

subroutine new_patch_sfc_props_ar(csite,ipa, rhos)
  
  use ed_state_vars,only:sitetype,patchtype
  use grid_coms, only: nzg, nzs
  use soil_coms, only: soil,slz
  use therm_lib, only: qwtk8,qtk
  use ed_therm_lib,only:ed_grndvap
  use consts_coms, only: wdns
  
  implicit none
  integer :: ipa,ico
  type(sitetype), target :: csite
  type(patchtype), pointer :: cpatch
  integer :: k
  real :: surface_temp, surface_fliq
  real, intent(in) :: rhos
  
  do k = 1, nzg
     call qwtk8(csite%soil_energy(k,ipa), csite%soil_water(k,ipa)*dble(wdns),  &
          soil(csite%ntext_soil(k,ipa))%slcpd, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa))
  enddo

  csite%nlev_sfcwater(ipa) = 0
  k = 1
  do while(csite%sfcwater_mass(min(k,nzs),ipa) > 1.0e-6 .and. k <= nzs)
     csite%nlev_sfcwater(ipa) = k
     csite%sfcwater_energy(k,ipa) = csite%sfcwater_energy(k,ipa) / csite%sfcwater_mass(k,ipa)
     call qtk(csite%sfcwater_energy(k,ipa), csite%sfcwater_tempk(k,ipa),   &
          csite%sfcwater_fracliq(k,ipa))
     k = k+1
  enddo
  
  call ed_grndvap(csite%nlev_sfcwater(ipa), csite%ntext_soil(nzg,ipa),   &
       csite%soil_water(nzg,ipa), csite%soil_energy(nzg,ipa),    &
       csite%sfcwater_energy(max(1,csite%nlev_sfcwater(ipa)),ipa), rhos,   &
       csite%can_shv(ipa), csite%ground_shv(ipa), csite%surface_ssh(ipa), &
       surface_temp, surface_fliq)
  
  
  !---- paw_avg10d is a 10-day running average of available water. Initialize it with the current value.
  cpatch => csite%patch(ipa)
  do ico = 1, cpatch%ncohorts
     cpatch%paw_avg10d(ico) = 0.0
     do k = cpatch%krdepth(ico), nzg - 1
        cpatch%paw_avg10d(ico) = cpatch%paw_avg10d(ico)                           &
             + real(csite%soil_water(k,ipa)-dble(soil(csite%ntext_soil(k,ipa))%soilcp))     &
             * (slz(k+1)-slz(k))/(soil(csite%ntext_soil(k,ipa))%slmsts            &
             - soil(csite%ntext_soil(k,ipa))%soilcp) 
     end do
     cpatch%paw_avg10d(ico)= cpatch%paw_avg10d(ico) + real(csite%soil_water(nzg,ipa)  &
          -dble(soil(csite%ntext_soil(nzg,ipa))%soilcp)) &
          *(-1.0*slz(nzg))              &
          /(soil(csite%ntext_soil(nzg,ipa))%slmsts &
          -soil(csite%ntext_soil(nzg,ipa))%soilcp) 
     cpatch%paw_avg10d(ico) = cpatch%paw_avg10d(ico)/(-1.0*slz(cpatch%krdepth(ico)))
  end do
  
  return
end subroutine new_patch_sfc_props_ar
