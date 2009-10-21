!====================================================================
! ============================================

subroutine init_ed_cohort_vars(cpatch,ico, lsl)
  
  use ed_state_vars,only : patchtype
  use allometry, only: calc_root_depth, assign_root_depth
  use pft_coms, only : leaf_turnover_rate, Vm0, sla

  implicit none

  type(patchtype),target :: cpatch
  integer :: ico
  real :: root_depth
  integer, intent(in) :: lsl
  
  cpatch%solvable(ico) = .false.

  cpatch%mean_gpp(ico)        = 0.0
  cpatch%mean_leaf_resp(ico)  = 0.0
  cpatch%mean_root_resp(ico)  = 0.0
  
  cpatch%dmean_leaf_resp(ico) = 0.0
  cpatch%dmean_root_resp(ico) = 0.0
  cpatch%dmean_gpp(ico)       = 0.0
  cpatch%dmean_gpp_pot(ico)   = 0.0
  cpatch%dmean_gpp_max(ico)   = 0.0

  cpatch%mmean_gpp(ico)          = 0.0
  cpatch%mmean_leaf_resp(ico)    = 0.0
  cpatch%mmean_root_resp(ico)    = 0.0
  cpatch%mmean_growth_resp(ico)  = 0.0
  cpatch%mmean_storage_resp(ico) = 0.0
  cpatch%mmean_vleaf_resp(ico)   = 0.0
  
  cpatch%dmean_light_level(ico) = 0.0
  cpatch%mmean_light_level(ico) = 0.0
  
  cpatch%dmean_fs_open(ico)     = 0.0
  cpatch%dmean_fsw(ico)         = 0.0
  cpatch%dmean_fsn(ico)         = 0.0
  
  cpatch%mmean_fs_open(ico)     = 0.0
  cpatch%mmean_fsw(ico)         = 0.0
  cpatch%mmean_fsn(ico)         = 0.0


  cpatch%gpp(ico) = 0.0
  cpatch%leaf_respiration(ico) = 0.0
  cpatch%root_respiration(ico) = 0.0
  cpatch%growth_respiration(ico) = 0.0
  cpatch%storage_respiration(ico) = 0.0
  cpatch%vleaf_respiration(ico) = 0.0

  cpatch%fsw(ico) = 1.0
  cpatch%fsn(ico) = 1.0
  
  !----- This variable would never be assigned for low LAI cohorts 
  cpatch%fs_open(ico) = cpatch%fsw(ico)*cpatch%fsn(ico)

  cpatch%monthly_dndt(ico)      = 0.0
  cpatch%mort_rate(:,ico)       = 0.0
  cpatch%mmean_mort_rate(:,ico) = 0.0

  cpatch%dagb_dt(ico)          = 0.0
  cpatch%dba_dt(ico)           = 0.0
  cpatch%ddbh_dt(ico)          = 0.0


  cpatch%light_level(ico)      = 0.0
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
  cpatch%paw_avg(ico)             = 0.5 !0.0 - [KIM] starting from the mid point.  if starting from the driest point, plants'll drop leaves initially due to the water stress


  cpatch%Psi_open(ico) = 0.0

  cpatch%cb(1:13,ico) = 1.0
  cpatch%cb_max(1:13,ico) = 1.0
  cpatch%cbr_bar(ico) = 1.0

  cpatch%old_stoma_data(ico)%recalc = 1
  cpatch%old_stoma_data(ico)%T_L = 0.0
  cpatch%old_stoma_data(ico)%e_A              = 0.0
  cpatch%old_stoma_data(ico)%PAR              = 0.0
  cpatch%old_stoma_data(ico)%rb_factor        = 0.0
  cpatch%old_stoma_data(ico)%prss             = 0.0
  cpatch%old_stoma_data(ico)%phenology_factor = 0.0
  cpatch%old_stoma_data(ico)%gsw_open         = 0.0
  cpatch%old_stoma_data(ico)%ilimit           = 0
  cpatch%old_stoma_data(ico)%T_L_residual     = 0.0
  cpatch%old_stoma_data(ico)%e_a_residual     = 0.0
  cpatch%old_stoma_data(ico)%par_residual     = 0.0
  cpatch%old_stoma_data(ico)%rb_residual      = 0.0
  cpatch%old_stoma_data(ico)%leaf_residual    = 0.0
  cpatch%old_stoma_data(ico)%gsw_residual     = 0.0
  
  
  cpatch%old_stoma_vector(:,ico) = 0.
  cpatch%old_stoma_vector(1,ico) = 1.


  root_depth = calc_root_depth(cpatch%hite(ico),cpatch%dbh(ico), cpatch%pft(ico))
  cpatch%krdepth(ico) = assign_root_depth(root_depth, lsl)
  
  cpatch%first_census(ico) = 1
  cpatch%new_recruit_flag(ico) = 0
  cpatch%bseeds(ico) = 0.0

  cpatch%hcapveg(ico)    = 0.
  cpatch%veg_energy(ico) = 0.
  cpatch%veg_temp(ico)   = 0.
  cpatch%veg_water(ico)  = 0.
  cpatch%veg_fliq(ico)   = 0.

  cpatch%turnover_amp(ico) = 1.0
  
  if (leaf_turnover_rate(cpatch%pft(ico)) > 0.0) then
     cpatch%llspan(ico) = 12.0/leaf_turnover_rate(cpatch%pft(ico)) !in month
  else
     cpatch%llspan(ico) = 9999.
  end if
  cpatch%vm_bar(ico) = Vm0(cpatch%pft(ico))
  cpatch%sla(ico) = sla(cpatch%pft(ico))

  return
end subroutine init_ed_cohort_vars

! ==========================================

subroutine init_ed_patch_vars(csite,ip1,ip2,lsl)
  
  use ed_state_vars,only:sitetype
  use ed_max_dims, only: n_pft
!  use fuse_fiss_utils, only: count_cohorts
  use grid_coms,     only: nzs, nzg
  use soil_coms, only: slz
  use canopy_air_coms, only : veg_height_min, minimum_canopy_depth

  implicit none

  type(sitetype), target :: csite
  integer :: ipft,ncohorts,ipa
  integer, intent(in) :: ip1,ip2,lsl
  
  do ipa = ip1,ip2
     do ipft = 1,n_pft
        csite%old_stoma_data_max(ipft,ipa)%recalc = 1
     enddo
  enddo


  ! Initialize sfcwater state variables
  csite%sfcwater_mass(1:nzs,ip1:ip2)   = 0.0
  csite%sfcwater_energy(1:nzs,ip1:ip2) = 0.0
  csite%sfcwater_depth(1:nzs,ip1:ip2)  = 0.0
  csite%total_snow_depth(ip1:ip2)      = 0.0
  csite%snowfac(ip1:ip2)               = 0.0

  csite%rshort_s(:,ip1:ip2) = 0.0
  csite%rshort_s_beam(:,ip1:ip2) = 0.0
  csite%rshort_s_diffuse(:,ip1:ip2) = 0.0
  csite%nlev_sfcwater(ip1:ip2) = 0
  
  csite%rlong_s(ip1:ip2) = 0.0
  
  csite%avg_daily_temp(ip1:ip2) = 0.0

  csite%mean_rh(ip1:ip2) = 0.0
  csite%mean_nep(ip1:ip2) = 0.0
    
  csite%A_decomp(ip1:ip2)             = 0.0
  csite%f_decomp(ip1:ip2)             = 0.0
  csite%rh(ip1:ip2)                   = 0.0
  csite%cwd_rh(ip1:ip2)               = 0.0
  csite%fuse_flag(ip1:ip2)            = 0.0
  csite%plant_ag_biomass(ip1:ip2)     = 0.0

  csite%mean_runoff(ip1:ip2) = 0.0
  csite%mean_wflux(ip1:ip2) = 0.0
  csite%mean_latflux(ip1:ip2) = 0.0
  csite%mean_qrunoff(ip1:ip2) = 0.0
  csite%mean_hflux(ip1:ip2) = 0.0

  csite%dmean_rh(ip1:ip2) = 0.0
  csite%mmean_rh(ip1:ip2) = 0.0
  
  csite%dmean_A_decomp(ip1:ip2) = 0.0
  csite%dmean_Af_decomp(ip1:ip2) = 0.0

  csite%repro(1:n_pft,ip1:ip2) = 0.0
  csite%A_o_max(1:n_pft,ip1:ip2) = 0.0
  csite%A_c_max(1:n_pft,ip1:ip2) = 0.0

  csite%htry(ip1:ip2) = 1.0

  csite%co2budget_gpp(ip1:ip2)            = 0.0
  csite%co2budget_gpp_dbh(:,ip1:ip2)      = 0.0
  csite%co2budget_rh(ip1:ip2)             = 0.0
  csite%co2budget_plresp(ip1:ip2)         = 0.0
  csite%co2budget_initialstorage(ip1:ip2) = 0.0
  csite%co2budget_loss2atm(ip1:ip2)       = 0.0
  csite%co2budget_residual(ip1:ip2)       = 0.0
  csite%wbudget_precipgain(ip1:ip2)       = 0.0
  csite%wbudget_loss2atm(ip1:ip2)         = 0.0
  csite%wbudget_loss2runoff(ip1:ip2)      = 0.0
  csite%wbudget_loss2drainage(ip1:ip2)    = 0.0
  csite%wbudget_initialstorage(ip1:ip2)   = 0.0
  csite%wbudget_residual(ip1:ip2)         = 0.0
  csite%ebudget_precipgain(ip1:ip2)       = 0.0
  csite%ebudget_netrad(ip1:ip2)           = 0.0
  csite%ebudget_loss2atm(ip1:ip2)         = 0.0
  csite%ebudget_loss2runoff(ip1:ip2)      = 0.0
  csite%ebudget_loss2drainage(ip1:ip2)    = 0.0
  csite%ebudget_latent(ip1:ip2)           = 0.0
  csite%ebudget_initialstorage(ip1:ip2)   = 0.0
  csite%ebudget_residual(ip1:ip2)         = 0.0
  csite%dmean_co2_residual(ip1:ip2)       = 0.0
  csite%dmean_energy_residual(ip1:ip2)    = 0.0
  csite%dmean_water_residual(ip1:ip2)     = 0.0
  csite%mmean_co2_residual   (ip1:ip2)     = 0.0
  csite%mmean_energy_residual(ip1:ip2)     = 0.0
  csite%mmean_water_residual (ip1:ip2)     = 0.0

  !----------------------------------------------------------------------------------------!
  !    These variables need to be initialized here otherwise it will fail when new patches !
  ! are created.                                                                           !
  !----------------------------------------------------------------------------------------!
  csite%avg_carbon_ac    (ip1:ip2)  = 0.0
  csite%avg_vapor_vc     (ip1:ip2)  = 0.0
  csite%avg_dew_cg       (ip1:ip2)  = 0.0
  csite%avg_vapor_gc     (ip1:ip2)  = 0.0
  csite%avg_wshed_vg     (ip1:ip2)  = 0.0
  csite%avg_vapor_ac     (ip1:ip2)  = 0.0
  csite%avg_transp       (ip1:ip2)  = 0.0
  csite%avg_evap         (ip1:ip2)  = 0.0
  csite%avg_netrad       (ip1:ip2)  = 0.0
  csite%avg_runoff       (ip1:ip2)  = 0.0
  csite%avg_drainage     (ip1:ip2)  = 0.0
  csite%aux              (ip1:ip2)  = 0.0
  csite%avg_sensible_vc  (ip1:ip2)  = 0.0
  csite%avg_qwshed_vg    (ip1:ip2)  = 0.0
  csite%avg_sensible_gc  (ip1:ip2)  = 0.0
  csite%avg_sensible_ac  (ip1:ip2)  = 0.0
  csite%avg_runoff_heat  (ip1:ip2)  = 0.0
  csite%avg_sensible_gg(:,ip1:ip2)  = 0.0
  csite%avg_smoist_gg  (:,ip1:ip2)  = 0.0
  csite%avg_smoist_gc  (:,ip1:ip2)  = 0.0
  csite%aux_s          (:,ip1:ip2)  = 0.0

  csite%avg_veg_energy(ip1:ip2)     = 0.0
  csite%avg_veg_temp(ip1:ip2)       = 0.0
  csite%avg_veg_fliq(ip1:ip2)       = 0.0
  csite%avg_veg_water(ip1:ip2)      = 0.0

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
  
  csite%watertable(ip1:ip2)                  = slz(lsl)
  csite%ustar(ip1:ip2) = 0.0
  csite%tstar(ip1:ip2) = 0.0
  csite%qstar(ip1:ip2) = 0.0
  csite%cstar(ip1:ip2) = 0.0
  csite%upwp(ip1:ip2)  = 0.0
  csite%tpwp(ip1:ip2)  = 0.0
  csite%qpwp(ip1:ip2)  = 0.0
  csite%cpwp(ip1:ip2)  = 0.0
  csite%wpwp(ip1:ip2)  = 0.0

  csite%can_enthalpy(ip1:ip2) = 0.0
  csite%can_temp    (ip1:ip2) = 0.0
  csite%can_rhos    (ip1:ip2) = 0.0
  csite%ground_shv  (ip1:ip2) = 0.0
  csite%surface_ssh (ip1:ip2) = 0.0

  csite%old_stoma_data_max(:,ip1:ip2)%recalc = 1
  csite%old_stoma_data_max(:,ip1:ip2)%T_L = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%e_A              = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%PAR              = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%rb_factor        = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%prss             = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%phenology_factor = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%gsw_open         = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%ilimit           = 0
  csite%old_stoma_data_max(:,ip1:ip2)%T_L_residual     = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%e_a_residual     = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%par_residual     = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%rb_residual      = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%leaf_residual    = 0.0
  csite%old_stoma_data_max(:,ip1:ip2)%gsw_residual     = 0.0
  
  
  csite%old_stoma_vector_max(:,:,ip1:ip2) = 0.
  csite%old_stoma_vector_max(1,:,ip1:ip2) = real(csite%old_stoma_data_max(:,ip1:ip2)%recalc)

  ncohorts = 0
  do ipa=1,csite%npatches
     ncohorts = ncohorts + csite%patch(ipa)%ncohorts
  enddo

  csite%cohort_count = ncohorts

  return
end subroutine init_ed_patch_vars

!======================================================================

subroutine init_ed_site_vars(cpoly, lat)

  use ed_state_vars,only:polygontype
  use ed_max_dims, only: n_pft, n_dbh, n_age, n_dist_types 
  use pft_coms, only: agri_stock,plantation_stock
  use grid_coms, only: nzs, nzg

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

  ! Initialising minimum monthly temperature with a very large value, to be reduced
  ! as the canopy temperature is updated.
  cpoly%min_monthly_temp(:) = huge(1.)

 ! cpoly%mm_gpp(:) = 0.0
 ! cpoly%mm_plresp(:) = 0.0
 ! cpoly%mm_rh(:) = 0.0
 ! cpoly%mm_nep(:) = 0.0

!  cpoly%mean_precip(:) = 0.0
!  cpoly%mean_qprecip(:) = 0.0
!  cpoly%mean_netrad(:) = 0.0

  cpoly%lambda_fire(1:12,:) = 0.0
  
  cpoly%disturbance_memory(1:n_dist_types, 1:n_dist_types,:) = 0.0
  cpoly%disturbance_rates(1:n_dist_types, 1:n_dist_types,:) = 0.0

  cpoly%agri_stocking_density(:) = 10.0

!KIM - 
!    - anyway, this part should be more elaborate for the case 
!    - that we have different crops/pastures.
! It's now defined in ED2IN, but it assumes only one PFT. It is probably not the
! ideal solution for regional runs...
  cpoly%agri_stocking_pft(:) = agri_stock
  cpoly%plantation_stocking_pft(:) = plantation_stock

  cpoly%plantation_stocking_density(:) = 4.0

  cpoly%primary_harvest_memory(:) = 0.0
  cpoly%secondary_harvest_memory(:) = 0.0
  
  cpoly%rad_avg(:) = 0.0
  
  
  return
end subroutine init_ed_site_vars

!======================================================================
subroutine init_ed_poly_vars(cgrid)
  
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
   end do

   return
end subroutine init_ed_poly_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will assign the values of some diagnostic variables, such as soil    !
! and temporary layer temperature and liquid fraction, and the surface properties.         !
!------------------------------------------------------------------------------------------!
subroutine new_patch_sfc_props(csite,ipa)
   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         ! ! structure
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use soil_coms     , only : soil              & ! intent(in), look-up table
                            , slz               & ! intent(in)
                            , min_sfcwater_mass ! ! intent(in)
   use consts_coms   , only : wdns              ! ! intent(in)
   use therm_lib     , only : qwtk              & ! subroutine
                            , qtk               ! ! subroutine
   use ed_therm_lib  , only : ed_grndvap        ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype) , target     :: csite          ! Current site
   integer        , intent(in) :: ipa            ! Number of the current patch
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype), pointer    :: cpatch         ! Current patch
   integer                     :: k              ! Layer counter
   integer                     :: ico            ! Cohort counter
   integer                     :: nsoil          ! Alias for soil texture class
   real                        :: surface_temp   ! Scratch variable for ed_grndvap
   real                        :: surface_fliq   ! Scratch variable for ed_grndvap
   !---------------------------------------------------------------------------------------!
  
   !----- Finding soil temperature and liquid water fraction. -----------------------------!
   do k = 1, nzg
      nsoil = csite%ntext_soil(k,ipa)
      call qwtk(csite%soil_energy(k,ipa), csite%soil_water(k,ipa)*wdns                     &
                ,soil(nsoil)%slcpd, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa))
   end do
   !---------------------------------------------------------------------------------------! 



   !---------------------------------------------------------------------------------------! 
   !   Determining number of temporary snow/surface water layers.  This is done by check-  !
   ! ing the mass.  If it is a layer, then we convert sfcwater_energy from J/m2 to J/kg,   !
   ! and compute the temperature and liquid fraction.                                      !
   !---------------------------------------------------------------------------------------! 
   csite%nlev_sfcwater(ipa) = 0
   snowloop: do k=1,nzs
      !----- Leave the loop if there is not enough mass in this layer... ------------------!
      if (csite%sfcwater_mass(k,ipa) <= min_sfcwater_mass)  exit snowloop
      csite%nlev_sfcwater(ipa) = k
      csite%sfcwater_energy(k,ipa) = csite%sfcwater_energy(k,ipa)                          &
                                   / csite%sfcwater_mass(k,ipa)
      call qtk(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa)                    &
              ,csite%sfcwater_fracliq(k,ipa))
   end do snowloop
   !---------------------------------------------------------------------------------------!
   !     Now, just to be safe, we will assign zeroes to layers above.                      !
   !---------------------------------------------------------------------------------------!
   do k=csite%nlev_sfcwater(ipa)+1,nzs
      csite%sfcwater_mass(k,ipa)   = 0.
      csite%sfcwater_energy(k,ipa) = 0.
      csite%sfcwater_depth(k,ipa)  = 0.
      if (k == 1) then
         csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(nzg,ipa)
         csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa) 
      else
         csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
         csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
      end if
   end do
   !---------------------------------------------------------------------------------------! 


   
   !----- Now we can compute the surface properties. --------------------------------------!
   k=max(1,csite%nlev_sfcwater(ipa))
   call ed_grndvap(csite%nlev_sfcwater(ipa),csite%ntext_soil(nzg,ipa)                      &
                  ,csite%soil_water(nzg,ipa),csite%soil_energy(nzg,ipa)                    &
                  ,csite%sfcwater_energy(k,ipa), csite%can_rhos(ipa),csite%can_shv(ipa)    &
                  ,csite%ground_shv(ipa),csite%surface_ssh(ipa),surface_temp,surface_fliq)
   !---------------------------------------------------------------------------------------! 



   !---------------------------------------------------------------------------------------!
   !      paw_avg is a 10-day running average of available water. Initialize it with    !
   ! the current value.                                                                    !
   !  MLO: I don't think this is currently used, but if we do use this, then we should add !
   !       a flag to initialise like this only if it is a new patch.  This is also called  !
   !       during patch fusion, when the 10-day running average is available...            !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   do ico = 1, cpatch%ncohorts
      cpatch%paw_avg(ico) = 0.0

      do k = cpatch%krdepth(ico), nzg - 1
         nsoil = csite%ntext_soil(k,ipa)
         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)                                      &
                             + (csite%soil_water(k,ipa) - soil(nsoil)%soilcp)           &
                             * (slz(k+1)-slz(k))                                        &
                             / (soil(nsoil)%slmsts - soil(nsoil)%soilcp) 
      end do
      nsoil = csite%ntext_soil(nzg,ipa)
      cpatch%paw_avg(ico) = cpatch%paw_avg(ico)                                         &
                          + (csite%soil_water(nzg,ipa) - soil(nsoil)%soilcp)            &
                          * (-1.0*slz(nzg)) / (soil(nsoil)%slmsts - soil(nsoil)%soilcp) 
      cpatch%paw_avg(ico) = cpatch%paw_avg(ico)/(-1.0*slz(cpatch%krdepth(ico)))
   end do
   !---------------------------------------------------------------------------------------! 
 
   return
end subroutine new_patch_sfc_props
!==========================================================================================!
!==========================================================================================!
