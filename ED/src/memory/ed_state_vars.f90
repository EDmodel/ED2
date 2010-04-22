module ed_state_vars

  use grid_coms, only: nzg,nzs,ngrids
  use c34constants, only : stoma_data,n_stoma_atts
  use ed_max_dims, only: max_site,n_pft,n_dbh,n_age,n_mort,n_dist_types      &
                       ,maxmach,maxgrds
  use disturb_coms, only : lutime,num_lu_trans,max_lu_years
  use met_driver_coms, only: met_driv_data,met_driv_state
  use fusion_fission_coms, only: ff_ndbh
  use phenology_coms, only: prescribed_phen
  use ed_misc_coms, only: idoutput, imoutput

  implicit none
!============================================================================!
!============================================================================!




!============================================================================!
!============================================================================!
  !------------------------------------------------------------------------!
  !    Some constants. All ED-structure variables are now initialized with !
  ! these huge numbers, so the pointers will point to numbers they can     !
  ! represent. They are not the largest to avoid overflow, but they are    !
  ! well off the acceptable range for most variables. That way all of them !
  ! must be initialised properly afterwards.                               !
  !------------------------------------------------------------------------!
  real    , parameter     :: large_real=1.e34
  real(kind=8), parameter :: large_double=1.d+34
  integer , parameter     :: large_integer=huge(1)-huge(1)/10
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!

  !--------------------------------------------------------------------------!  
  ! Patch type:
  ! The following are arrays of cohorts that populate the current patch.
  ! This is the lowest (most infantile) part of the memory structure.
  !--------------------------------------------------------------------------!  

  type patchtype

     integer :: ncohorts

     ! Global index of the first cohort across all cohorts

     integer :: coglob_id

     !================================================================
     ! Plant functional type:
     ! 1 - C4 grass
     ! 2 - tropical early successional tree
     ! 3 - tropical mid successional tree
     ! 4 - tropical late successional tree
     ! 5 - C3 grass
     ! 6 - northern pines (temperate)
     ! 7 - southern pines (temperate)
     ! 8 - late successional conifers
     ! 9 - early successional cold-deciduous hardwood
     ! 10 - mid successional cold-deciduous hardwood
     ! 11 - late successional cold-deciduous hardwood
     ! 12 - c3 pasture
     ! 13 - c3 crop (e.g.,wheat, rice, soybean) 
     ! 14 - c4 pasture
     ! 15 - c4 crop (e.g.,corn/maize)

     integer ,pointer,dimension(:) :: pft

     ! Density of plants (number/m2)
     real ,pointer,dimension(:) :: nplant

     ! Plant height (m)
     real ,pointer,dimension(:) :: hite

     ! Above-ground biomass (kgC/plant)
     real, pointer, dimension(:) :: agb

     ! Basal area (cm2)
     real, pointer, dimension(:) :: basarea
     
     ! Rate of change in agb (kgC/plant/yr)
     real, pointer, dimension(:) :: dagb_dt
     
     ! Rate of change in basal area (cm2/yr)
     real, pointer, dimension(:) :: dba_dt
     
     ! Rate of change in dbh (cm/yr)
     real, pointer, dimension(:) :: ddbh_dt

     ! Plant diameter at breast height (cm)
     real ,pointer,dimension(:) :: dbh

     ! Biomass of the structural stem (kgC/plant)
     real ,pointer,dimension(:) :: bdead

     ! Biomass of leaves (kgC/plant)
     real ,pointer,dimension(:) :: bleaf
     

     ! phenology_status codes:
     ! 0 - plant has the maximum LAI, given its size
     ! 1 - plant has an LAI between 0 and its maximum
     ! 2 - plant has no leaves
     integer ,pointer,dimension(:) :: phenology_status

     ! Biomass of live tissue (kgC/plant)
     real ,pointer,dimension(:) :: balive

     ! Biomass of fine roots (kgC/plant)
     real, pointer, dimension(:) :: broot
     
     ! Biomass of sapwood (kgC/plant)
     real, pointer, dimension(:) :: bsapwood

     ! Leaf area index (m2 leaf / m2 ground)
     real ,pointer,dimension(:) :: lai

     ! Wood projected area (m2 wood / m2 ground)
     real ,pointer,dimension(:) :: wpa

     ! Wood area index (m2 wood / m2 ground)
     real ,pointer,dimension(:) :: wai

     ! Logical test to check whether this cohort can be solved...
     logical, pointer, dimension(:) :: solvable

     ! Plant storage pool of carbon [kgC/plant]
     real ,pointer,dimension(:) :: bstorage
     
     ! Monthly carbon balance for past 12 months and the current month
     ! (kgC/plant/yr) - 13th column will have the partial month integration
     real, pointer,dimension(:,:) :: cb       !(13,ncohorts)
     
     ! Maximum monthly carbon balance for past 12 months and the current 
     ! month  if cohort were at the top of the canopy (kgC/plant, over one month)
     ! 13th column will have the partial month integration.
     real, pointer,dimension(:,:) :: cb_max  !(13,ncohorts)
     
     ! Annual average ratio of cb/cb_max
     real ,pointer,dimension(:) :: cbr_bar
     
     ! Monthly mean of cb. The only difference from cb is the way it is scaled, so
     ! all months have the same "weight"
     real, pointer, dimension(:) :: mmean_cb

     ! Vegetation internal energy (J/m2 ground)
     real ,pointer,dimension(:) :: veg_energy

     ! Vegetation temperature (K)
     real ,pointer,dimension(:) :: veg_temp

     ! Fraction of liquid water (dimensionless)
     real ,pointer,dimension(:) :: veg_fliq

     ! Vegetation surface water (kg/m2 ground)
     real ,pointer,dimension(:) :: veg_water

     ! Leaf (vegetation) heat capacity [ J/m2/K]
     real, pointer,dimension(:) :: hcapveg

     ! Gross primary productivity (GPP) [umol/m2/s], averaged over the 
     ! output frequency (FRQSTATE)
     real ,pointer,dimension(:) :: mean_gpp

     ! Mean leaf respiration rate (umol/m2 ground/s), averaged over FRQSTATE
     real ,pointer,dimension(:) :: mean_leaf_resp

     ! Mean root respiration rate (umol/m2 ground/s), averaged over FRQSTATE
     real ,pointer,dimension(:) :: mean_root_resp

     ! Mean growth respiration rate (umol/m2 ground/s), averaged over FRQSTATE
     real ,pointer,dimension(:) :: mean_growth_resp

     ! Mean storage respiration rate (umol/m2 ground/s), averaged over FRQSTATE
     real ,pointer,dimension(:) :: mean_storage_resp

     ! Mean vleaf respiration rate (umol/m2 ground/s), averaged over FRQSTATE
     real ,pointer,dimension(:) :: mean_vleaf_resp


     ! Mean leaf respiration rate (umol/m2 ground/s), averaged over 1 day
     real ,pointer,dimension(:) :: today_leaf_resp

     ! Mean root respiration rate (umol/m2 ground/s), averaged over 1 day
     real ,pointer,dimension(:) :: today_root_resp

     ! Gross primary productivity (GPP) [umol/m2 ground/s], averaged over 1 day
     real ,pointer,dimension(:) :: today_gpp

     ! Potential GPP in the absence of N limitation [umol/m2 ground/s],
     ! averaged over 1 day
     real ,pointer,dimension(:) :: today_gpp_pot

     ! Maximum GPP if cohort were at the top of the canopy 
     ! [umol/m2 ground/s], averaged over 1 day
     real ,pointer,dimension(:) :: today_gpp_max

     ! Plant growth respiration (kgC/plant/day)
     real ,pointer,dimension(:) :: growth_respiration

     ! Plant storage respiration (kgC/plant/day)
     real ,pointer,dimension(:) :: storage_respiration

     ! Plant virtual leaf respiration (kgC/plant/day)
     real ,pointer,dimension(:) :: vleaf_respiration

     ! Daily and monthly means of the plant productivity terms, all in kgC/m2/yr
     real, pointer, dimension(:) :: dmean_gpp
     real, pointer, dimension(:) :: dmean_leaf_resp
     real, pointer, dimension(:) :: dmean_root_resp
     real, pointer, dimension(:) :: mmean_gpp
     real, pointer, dimension(:) :: mmean_leaf_resp
     real, pointer, dimension(:) :: mmean_root_resp
     real, pointer, dimension(:) :: mmean_growth_resp
     real, pointer, dimension(:) :: mmean_storage_resp
     real, pointer, dimension(:) :: mmean_vleaf_resp

     ! Weighting factor for open and closed stomata due to N limitation
     real ,pointer,dimension(:) :: fsn

     ! Plant density tendency [plants/m2/month]
     real ,pointer,dimension(:) :: monthly_dndt

     ! Mortality rate [yr-1]
     real , pointer, dimension(:,:) :: mort_rate

     ! Monthly mean Mortality rate [yr-1] (only frost mortality changes...)
     real , pointer, dimension(:,:) :: mmean_mort_rate

     ! This is where you keep the derivatives of the 
     ! stomatal conductance residual and the old met conditions.
     type(stoma_data) ,pointer,dimension(:) :: old_stoma_data

     ! This vector is just for transporting the data into and out
     ! of the HDF5 file
     real ,pointer,dimension(:,:)           :: old_stoma_vector

     ! Transpiration rate, open stomata (mm/s)
     real ,pointer,dimension(:) :: Psi_open

     ! This specifies the index of the deepest soil layer of which the 
     ! cohort can access water.
     integer ,pointer,dimension(:) :: krdepth

     ! The model reports annual growth, mortality, cut and recruitment.  
     ! These rates are calculated with respect to two censuses.  This 
     ! is the flag specifying if a cohort was present at the first
     ! census.
     integer ,pointer,dimension(:) :: first_census

     ! Flag specifying if this cohort is a new recruit with respect
     ! to the first census.
     integer ,pointer,dimension(:) :: new_recruit_flag

     ! Light level received by this cohort, its diurnal and monthly means
     real, pointer, dimension(:) :: light_level
     real, pointer, dimension(:) :: dmean_light_level
     real, pointer, dimension(:) :: mmean_light_level
     real, pointer, dimension(:) :: light_level_beam
     real, pointer, dimension(:) :: dmean_light_level_beam
     real, pointer, dimension(:) :: mmean_light_level_beam
     real, pointer, dimension(:) :: light_level_diff
     real, pointer, dimension(:) :: dmean_light_level_diff
     real, pointer, dimension(:) :: mmean_light_level_diff
     real, pointer, dimension(:) :: beamext_level
     real, pointer, dimension(:) :: dmean_beamext_level
     real, pointer, dimension(:) :: mmean_beamext_level
     real, pointer, dimension(:) :: diffext_level
     real, pointer, dimension(:) :: dmean_diffext_level
     real, pointer, dimension(:) :: mmean_diffext_level
     real, pointer, dimension(:) :: norm_par_beam
     real, pointer, dimension(:) :: dmean_norm_par_beam
     real, pointer, dimension(:) :: mmean_norm_par_beam
     real, pointer, dimension(:) :: norm_par_diff
     real, pointer, dimension(:) :: dmean_norm_par_diff
     real, pointer, dimension(:) :: mmean_norm_par_diff

     ! Light extinction of this cohort, its diurnal and monthly means
     real, pointer, dimension(:) :: lambda_light
     real, pointer, dimension(:) :: dmean_lambda_light
     real, pointer, dimension(:) :: mmean_lambda_light

     ! Photosynthetically active radiation (PAR) absorbed by the 
     ! cohort (units are Einsteins/m2/s)
     real ,pointer,dimension(:) :: par_v
     real ,pointer,dimension(:) :: dmean_par_v
     real ,pointer,dimension(:) :: mmean_par_v

     ! Photosynthetically active radiation (PAR) absorbed by the 
     ! cohort (units are Einsteins/m2/s), beam component
     real ,pointer,dimension(:) :: par_v_beam
     real ,pointer,dimension(:) :: dmean_par_v_beam
     real ,pointer,dimension(:) :: mmean_par_v_beam

     ! Photosynthetically active radiation (PAR) absorbed by the 
     ! cohort (units are Einsteins/m2/s), diffuse component
     real ,pointer,dimension(:) :: par_v_diffuse
     real ,pointer,dimension(:) :: dmean_par_v_diff
     real ,pointer,dimension(:) :: mmean_par_v_diff

     ! Total short wave radiation absorbed by the cohort, W/m2
     real ,pointer,dimension(:) :: rshort_v

     ! Total short wave radiation absorbed by the cohort, W/m2, beam component
     real ,pointer,dimension(:) :: rshort_v_beam

     ! Total short wave radiation absorbed by the cohort, W/m2, diffuse 
     ! component
     real ,pointer,dimension(:) :: rshort_v_diffuse

     ! Total long wave radiation absorbed by the cohort (W/m2)
     real ,pointer,dimension(:) :: rlong_v

     ! Total long wave radiation absorbed by the cohort (W/m2), due to 
     ! the temperature of the vegetation and surface alone
     real ,pointer,dimension(:) :: rlong_v_surf

     ! Total long wave radiation absorbed by the cohort (W/m2), due to 
     ! the incident atmospheric long wave alone
     real ,pointer,dimension(:) :: rlong_v_incid

     ! Leaf aerodynamic resistance (s/m)
     real ,pointer,dimension(:) :: rb

     ! Photosynthesis rate, open stomata (umol/m2 leaf/s)
     real ,pointer,dimension(:) :: A_open

     ! Photosynthesis rate, closed stomata (umol/m2 leaf/s)
     real ,pointer,dimension(:) :: A_closed

     ! Transpiration rate, closed stomata (mm/s)
     real ,pointer,dimension(:) :: Psi_closed

     ! Stomatal resistance for water, open stomata (s/m)
     real ,pointer,dimension(:) :: rsw_open

     ! Stomatal resistance for water, closed stomata (s/m)
     real ,pointer,dimension(:) :: rsw_closed

     ! Weighting factor for open and closed stomata (fsw=1 => fully open)
     real ,pointer,dimension(:) :: fsw

     ! Product of fsw and fsn
     real ,pointer,dimension(:) :: fs_open

     ! Daily and monthly means of the weighting factors for open and closed stomata
     real ,pointer,dimension(:) :: dmean_fs_open
     real ,pointer,dimension(:) :: dmean_fsw
     real ,pointer,dimension(:) :: dmean_fsn
     real ,pointer,dimension(:) :: mmean_fs_open
     real ,pointer,dimension(:) :: mmean_fsw
     real ,pointer,dimension(:) :: mmean_fsn

     ! Net stomatal resistance [s/m]
     real ,pointer,dimension(:) :: stomatal_resistance

     ! Plant maintenance costs due to turnover of leaves and fine 
     ! roots [kgC/plant]
     real ,pointer,dimension(:) :: leaf_maintenance
     real ,pointer,dimension(:) :: root_maintenance
     ! Monthly mean, in kgC/plant/year
     real ,pointer,dimension(:) :: mmean_leaf_maintenance
     real ,pointer,dimension(:) :: mmean_root_maintenance
     
     ! Leaf loss to litter layer due to phenology [kgC/plant]
     real , pointer, dimension(:) :: leaf_drop
     ! Monthly mean, in kgC/plant/year
     real , pointer, dimension(:) :: mmean_leaf_drop

     ! Amount of seeds produced for dispersal [kgC/plant]
     real ,pointer,dimension(:) :: bseeds

     ! Instantaneous values of leaf and root respiration [umol/m2/s]
     real ,pointer,dimension(:) :: leaf_respiration
     real ,pointer,dimension(:) :: root_respiration

     ! Gross Primary Productivity [umol/m2/s]
     real, pointer,dimension(:) :: gpp

     ! Plant available water, (stress function)  0 to 1 [unitless]
     real, pointer, dimension(:) :: paw_avg

     ! Phenology-related
     real, pointer, dimension(:) :: turnover_amp
     real, pointer, dimension(:) :: llspan
     real, pointer, dimension(:) :: vm_bar
     real, pointer, dimension(:) :: sla

  end type patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !---------------------------------------------------------------------------!  
  ! Site type:
  ! The following are the patch level arrays that populate the current site.
  !---------------------------------------------------------------------------!  
  
  type sitetype
     

     integer :: npatches

     !  The global index of the first cohort in all patches
     integer,pointer,dimension(:) :: paco_id

     ! The number of cohorts in each patch
     integer,pointer,dimension(:) :: paco_n

     ! Global index of the first patch in this vector, across all patches
     ! on the grid

     integer :: paglob_id

     ! The patches containing the cohort arrays
     type(patchtype),pointer,dimension(:) :: patch


     ! Leaf area index (m2 leaf / m2 ground)
     real, pointer,dimension(:) :: lai

     ! Wood projected area (m2 wood / m2 ground)
     real, pointer,dimension(:) :: wpa

     ! Wood area index (m2 wood / m2 ground)
     real, pointer,dimension(:) :: wai

     ! Patch type index:
     !       Agriculture = 1
     !       Secondary Forest = 2
     !       Primary Forest = 3
     integer , pointer,dimension(:) :: dist_type

     ! Time since last disturbance (years)
     real , pointer,dimension(:) :: age

     ! Fractional area of the patch
     real , pointer,dimension(:) :: area

     ! Fractional area of the patch (considers lai weighting)
     real, pointer,dimension(:) :: laiarea

     ! Total patch heat capacity (J / m2 ground / K)
     real, pointer,dimension(:) :: hcapveg

     ! Soil carbon concentration, fast pool (kg/m2)
     real , pointer,dimension(:) :: fast_soil_C 

     ! Soil carbon concentration, slow pool (kg/m2)
     real , pointer,dimension(:) :: slow_soil_C 

     ! Soil carbon concentration, structural pool (kg/m2)
     real , pointer,dimension(:) :: structural_soil_C

     ! Soil lignin concentration, structural pool (kg/m2)
     real , pointer,dimension(:) :: structural_soil_L

     ! Soil nitrogen concentration, mineralized pool (kg/m2)
     real , pointer,dimension(:) :: mineralized_soil_N  

     ! Soil nitrogen concentration, fast pool (kg/m2)
     real , pointer,dimension(:) :: fast_soil_N 

     ! Number of degree days
     ! Degree days --- sum of daily average temperatures above 278.15 K 
     real , pointer,dimension(:) :: sum_dgd

     ! Chill days --- number of days with average temperatures below 278.15 K
     real , pointer,dimension(:) :: sum_chd  

     ! Flag specifying whether (1) or not (0) this patch is a plantation
     integer , pointer,dimension(:) :: plantation

     ! Enthalpy (J/kg) of canopy air
     real , pointer,dimension(:) :: can_enthalpy

     ! Temperature (K) of canopy air
     real , pointer,dimension(:) :: can_temp

     ! Water vapor specific humidity (kg/kg) of canopy air
     real , pointer,dimension(:) :: can_shv

     ! CO2 mixing ratio [umol/mol] of canopy air
     real , pointer,dimension(:) :: can_co2

     ! Density [kg/m3] of canopy air
     real , pointer,dimension(:) :: can_rhos

     ! Canopy air pressure.
     real , pointer,dimension(:) :: can_prss

     ! Canopy air potential temperature
     real , pointer,dimension(:) :: can_theta

     ! Canopy depth (m)
     real , pointer,dimension(:) :: can_depth

     ! Mean light extinction coefficient, its diurnal and monthly cycle.
     real , pointer, dimension(:) :: lambda_light
     real , pointer, dimension(:) :: dmean_lambda_light
     real , pointer, dimension(:) :: mmean_lambda_light


     ! Number of cohorts in the patch
     integer , pointer,dimension(:) :: cohort_count

     ! Patch name (only used when restarting from ED1)
     character(len=256), pointer,dimension(:) :: pname
     
     ! Number of surface water layers
     integer,pointer,dimension(:) :: nlev_sfcwater

     ! Surface water mass (kg/m2)
     real, pointer,dimension(:,:) :: sfcwater_mass

     ! Surface water internal energy (J/kg)
     real, pointer,dimension(:,:) :: sfcwater_energy

     ! Depth of surface water (m)
     real, pointer,dimension(:,:) :: sfcwater_depth

     ! Temperature of surface water (K)
     real, pointer,dimension(:,:) :: sfcwater_tempk
    
     ! Liquid fraction of surface water
     real, pointer,dimension(:,:) :: sfcwater_fracliq

     ! Short wave radiation absorbed by the surface water (W/m2)
     real, pointer,dimension(:,:) :: rshort_s

     ! Short wave radiation absorbed by the surface water, 
     ! beam component (W/m2)
     real, pointer,dimension(:,:) :: rshort_s_beam

     ! Short wave radiation absorbed by the surface water, 
     ! diffuse component (W/m2)
     real, pointer,dimension(:,:) :: rshort_s_diffuse

     ! Soil textural class index
     integer, pointer,dimension(:,:) :: ntext_soil
     
     ! Soil internal energy (J/m3)
     real,    pointer,dimension(:,:) :: soil_energy
     
     ! Soil water (m3 water / m3 soil)
     real,    pointer,dimension(:,:) :: soil_water
     
     ! Temperature of soil (K)
     real,    pointer,dimension(:,:) :: soil_tempk
     
     ! Liquid fraction of soil
     real,    pointer,dimension(:,:) :: soil_fracliq
     
     ! Effective specific humidity (kg/kg) just above soil
     real,    pointer,dimension(:) :: ground_shv
     
     ! Surface saturation specific humidity (kg/kg)
     real,    pointer,dimension(:) :: surface_ssh

     ! Net roughness length (m)
     real,    pointer,dimension(:)  :: rough

     ! Photosynthetic rates for different PFTs, if they were at the top 
     ! of the canopy (umol/m2 leaf/s).  Used by mortality function.
     real, pointer,dimension(:,:) :: A_o_max  ! open stomata
     real, pointer,dimension(:,:) :: A_c_max  ! closed stomata
     
     ! This will hold the stomatal conductance data from the previous 
     ! time step corresponding to A_o_max
     type(stoma_data), pointer,dimension(:,:)  :: old_stoma_data_max
     real, pointer,dimension(:,:,:)            :: old_stoma_vector_max

     ! Average daily temperature [K]
     real , pointer,dimension(:) :: avg_daily_temp  

     ! average of rh [umol/m2/s] over FRQSTATE
     real , pointer,dimension(:) :: mean_rh

     ! average of rh [umol/m2/s] over a day and a month, respectively
     real , pointer,dimension(:) :: dmean_rh
     real , pointer,dimension(:) :: mmean_rh

     ! average of net ecosystem productivity (NEP) [umol/m2/s] over FRQSTATE
     real , pointer,dimension(:) :: mean_nep

     !=======================================================================!
     !      These variables are terms of the water, carbon, and energy       !
     ! budget, and will compute what enters, what stays, and what leaves ED. !
     ! They are averaged over FRQSTATE, so the units are now per unit of     !
     ! time for all fluxes, and the total storage at the beginning of the    !
     ! next time step (or the final storage between the previous analysis    !
     ! and this one).                                                        !
     !=======================================================================!
     ! Mean moisture transfer from the canopy air to the atmosphere 
     ! (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2atm

     ! Contribution of density change on change in the final storage
     ! (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_denseffect

     ! Precipitation [kg_H2O/m2/s]
     real , pointer,dimension(:) :: wbudget_precipgain

     ! Mean runoff (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2runoff

     ! Mean drainage (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2drainage

     ! Total water (soil, sfcwater, can_shv, veg_water) at the beginning
     ! of budget-averaging time. [kg_H2O/m2]
     real , pointer,dimension(:) :: wbudget_initialstorage

     ! Residual of the water budget (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_residual

     ! Mean sensible heat transfer from the canopy air to the atmosphere
     ! (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2atm

     ! Mean change in storage due to change in density
     ! (J/m2/s)
     real , pointer,dimension(:) :: ebudget_denseffect

     ! Energy associated with runoff (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2runoff

     ! Energy associated with drainage (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2drainage

     ! Net absorbed radiation by soil, sfcwater, vegetation [J/m2/s]
     real , pointer,dimension(:) :: ebudget_netrad

     ! Mean latent heat flux (J/m2/s)
     real , pointer,dimension(:) :: ebudget_latent

     ! Energy associated with precipitation (J/m2/s)
     real , pointer,dimension(:) :: ebudget_precipgain

     ! Total energy (soil, sfcwater, can_shv, veg_water) at the beginning
     ! of budget-averaging time. [J/m2]
     real , pointer,dimension(:) :: ebudget_initialstorage

     ! Residual of the energy budget [J/m2/s]
     real, pointer, dimension(:) :: ebudget_residual

     ! Total CO2 (can_shv) at the beginning of budget-averaging 
     ! time. [umol_CO2/m2]
     real , pointer,dimension(:) :: co2budget_initialstorage

     ! Residual of the CO2 budget [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_residual

     ! Flux of CO2 from the canopy air to the atmosphere [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_loss2atm

     ! Change in CO2 total storage due to change in density [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_denseffect

     ! Average GPP [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_gpp

     ! Average GPP by GPP class[umol_CO2/m2/s]
     real , pointer,dimension(:,:) :: co2budget_gpp_dbh

     ! Average plant respiration [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_plresp

     ! Average heterotrophic respiration [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_rh

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

     ! Daily average of A_decomp, the temperature and moisture dependence
     ! of heterotrophic respiration.  The "today" variable is used in the
     ! model, whereas dmean and mmean are for output only.
     real , pointer,dimension(:) :: today_A_decomp
     real , pointer,dimension(:) :: dmean_A_decomp
     real , pointer,dimension(:) :: mmean_A_decomp

     ! Daily average of the product A_decomp * f_decomp, which incorporates
     ! temperature, moisture, and N dependence of decomposition.  The "today" 
     ! variable is used in the model, whereas dmean and mmean are for output only.
     real , pointer,dimension(:) :: today_Af_decomp
     real , pointer,dimension(:) :: dmean_Af_decomp
     real , pointer,dimension(:) :: mmean_Af_decomp

     ! Carbon available to establish recruits [kgC/m2]
     real, pointer,dimension(:,:) :: repro    !(n_pft,npatches)  

     ! Vegetation roughness length (m)
     real , pointer,dimension(:) :: veg_rough

     ! Vegetation height (m)
     real , pointer,dimension(:) :: veg_height

     ! Input to fast soil carbon pool [kgC/m2/day]
     real , pointer,dimension(:) :: fsc_in

     ! Input to structural soil carbon pool [kgC/m2/day]
     real , pointer,dimension(:) :: ssc_in

     ! Input to soil lignin pool [kgC/m2/day]
     real , pointer,dimension(:) :: ssl_in  

     ! Input to fast soil nitrogen pool [kgN/m2/day]
     real , pointer,dimension(:) :: fsn_in  

     ! Plant nitrogen update summed over all cohorts [kgN/m2/day]
     real , pointer,dimension(:) :: total_plant_nitrogen_uptake

     !real, pointer,dimension(:) :: Nnet_min
     !real, pointer,dimension(:) :: Ngross_min
     real, pointer,dimension(:) :: mineralized_N_input
     real, pointer,dimension(:) :: mineralized_N_loss

     ! Short wave radiation absorbed by the ground (W/m2)
     real , pointer,dimension(:) :: rshort_g
     
     ! Short wave radiation absorbed by the ground, beam component (W/m2)
     real , pointer,dimension(:) :: rshort_g_beam
     
     ! Short wave radiation absorbed by the ground, diffuse component (W/m2)
     real , pointer,dimension(:) :: rshort_g_diffuse
     
     ! Long wave radiation absorbed by the ground (W/m2)
     real , pointer,dimension(:) :: rlong_g

     ! Long wave radiation absorbed by the ground (W/m2), due to the 
     ! surface and vegetation alone
     real , pointer,dimension(:) :: rlong_g_surf

     ! Long wave radiation absorbed by the ground (W/m2), due to the 
     ! incident long wave alone
     real , pointer,dimension(:) :: rlong_g_incid

     ! Long wave radiation absorbed by the surface water (W/m2)
     real , pointer,dimension(:) :: rlong_s

     ! Long wave radiation absorbed by the surface water (W/m2), due to 
     ! the surface and vegetation alone
     real , pointer,dimension(:) :: rlong_s_surf

     ! Long wave radiation absorbed by the surface water (W/m2), due 
     ! to the incident atmospheric long wave alone
     real , pointer,dimension(:) :: rlong_s_incid

     ! Patch albedo
     real , pointer,dimension(:) :: albedt

     ! Patch albedo, beam component
     real , pointer,dimension(:) :: albedo_beam

     ! Patch albedo, diffuse component
     real , pointer,dimension(:) :: albedo_diffuse

     ! Upward long wave radiation at the top of the canopy (W/m2)
     real , pointer,dimension(:) :: rlongup

     ! Albedo for long wave radiation
     real , pointer,dimension(:) :: rlong_albedo

     ! Total snow depth as calculated in the radiation scheme.  Used for 
     ! checking if cohorts are buried.  [m]
     real , pointer,dimension(:) :: total_snow_depth

     ! Fraction of vegetation covered with snow.  Used for computing 
     ! surface roughness.
     real , pointer,dimension(:) :: snowfac

     ! limitation of heterotrophic respiration due to physical 
     ! environmental factors (0-1 coefficient)
     real , pointer,dimension(:) :: A_decomp

     ! damping of decomposition due to nitrogen immobilization 
     ! (0-1 coefficient)
     real , pointer,dimension(:) :: f_decomp

     ! total heterotrophic respiration (umol/m2/s)
     real , pointer,dimension(:) :: rh

     ! coarse woody debris contribution to rh (umol/m2/s)
     real , pointer,dimension(:) :: cwd_rh

     ! Integer flag specifying whether this patch is to be fused
     integer , pointer,dimension(:) :: fuse_flag

     ! Plant density broken down into size and PFT bins.  Used in patch fusion
     real, pointer,dimension(:,:,:) :: pft_density_profile !(n_pft,ff_ndbh,npatches)

     ! Above ground biomass in this patch [kgC/m2]
     real , pointer,dimension(:) :: plant_ag_biomass

     ! Mean water flux from the canopy air to the atmosphere (kg_H2O/mÂ²/s).
     real, pointer,dimension(:) :: mean_wflux

     ! Mean latent heat flux from the canopy air to the atmosphere (W/mÂ²).
     real, pointer,dimension(:) :: mean_latflux

     ! Mean sensible heat flux from the canopy air to the atmosphere (W/mÂ²).
     real, pointer,dimension(:) :: mean_hflux
     
     ! Mean runoff (kg_H2O/mÂ²/s).
     real, pointer,dimension(:) :: mean_runoff

     ! Energy associated with runoff (W/mÂ²), averaged over FRQSTATE
     real, pointer,dimension(:) :: mean_qrunoff

     ! Last time step successfully completed by integrator.
     real, pointer,dimension(:)  :: htry

     ! Average time step used by the integrator, its daily and monthly mean
     real, pointer, dimension(:) :: avg_rk4step
     real, pointer, dimension(:) :: dmean_rk4step
     real, pointer, dimension(:) :: mmean_rk4step

     ! Average available water for transpiration.
     real, pointer,dimension(:)  :: avg_available_water

     real, pointer,dimension(:)  :: ustar ! Friction velocity [m/s]

     real, pointer,dimension(:)  :: tstar ! Characteristic temperature fluctuation scale [K]

     real, pointer,dimension(:)  :: qstar ! Characteristic specific humidity fluct. scale [kg/kg]

     real, pointer,dimension(:)  :: cstar ! Characteristic CO2 mix. ratio fluct. scale [ppm]

     real, pointer,dimension(:)  :: upwp !eddy mom. flux u-prime w-prime

     real, pointer,dimension(:)  :: tpwp !eddy heat flux t-prime w-prime

     real, pointer,dimension(:)  :: qpwp !eddy moist. flux q-prime w-prime

     real, pointer,dimension(:)  :: cpwp !eddy CO2 flux flux c-prime w-prime

     real, pointer,dimension(:)  :: wpwp !eddy v. mom. flux w-prime w-prime

     ! Marcos: are we still leaving questions for each other in the code?  Or should I
     ! just send you an email?  RGK 8-4-09
     ! Ryan: I prefer e-mail because I can always add some new spam rules :)...

     ! -----------------------------------------------------------------
     ! Fast time flux diagnostic variables
     ! The following variables are all averaged over 
     ! the fast time flux period frqfast
     !------------------------------------------------------------------
     
     !----- Radiation --------------------------------------------------
     
     real,pointer,dimension(:) :: avg_netrad        ! Net radiation of soil,surface water and vegetation sum
                                                    ! [W/m2]

     ! Notation
     ! ATM - atmosphere, CAS - canopy air space

     !----- Moisture ---------------------------------------------------

     real,pointer,dimension(:)   :: avg_carbon_ac   ! Average carbon flux, ATM -> CAS,       [umol/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_vc    ! Average water vapor flux, VEG -> CAS   [kg/m2/s]
     real,pointer,dimension(:)   :: avg_dew_cg      ! Average dew flux, CAS -> GND           [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Average water vapor flux, GND -> CAS   [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Average flux of precip that bypasses
                                                    !  leaves at water holding capacity (throughfall)
     real,pointer,dimension(:)   :: avg_intercepted ! Average flux of precip that is intercepted
                                                    !                                        [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Average water vapor flux, ATM-CAS      [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Average transpiration of water vapor   [kg/m2/s]
     real,pointer,dimension(:)   :: avg_evap        ! Average evaporation from leaf and
                                                    !  and ground surfaces -> CAS            [kg/m2/s]
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between soil layers
                                                    !  where layer (nzg) is the flux from the
                                                    !  surface water into the top layer      [kg/m2/s]
     real,pointer,dimension(:,:) :: avg_smoist_gc   ! Transpired soil moisture sink in each layer
                                                    !                                        [kg/m2/s]
     real,pointer,dimension(:)   :: avg_runoff      ! Average surface water runoff           [kg/m2/s]
     real,pointer,dimension(:)   :: avg_drainage    ! Average water drainage through the lower 
                                                    !  soil layer                            [kg/m2/s]
     real,pointer,dimension(:)   :: avg_drainage_heat! Average internal energy loss due to water 
                                                     ! drainage through the lower soil layer [kg/m2/s]

     !----- Auxillary Variables (user can modify to view any variable ----------------------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable
     
     !----- Sensible heat ------------------------------------------------------------------!

     real,pointer,dimension(:) :: avg_sensible_vc   ! Sensible heat flux, VEG -> CAS         [W/m2]
     real,pointer,dimension(:) :: avg_qwshed_vg     ! Internal energy flux of precipitation
                                                    !  througfall                            [W/m2]
     real,pointer,dimension(:) :: avg_qintercepted  ! Internal energy flux of intercepted
                                                    !  precipitation                         [W/m2]
     real,pointer,dimension(:) :: avg_sensible_gc   ! Sensible heat flux, GND -> CAS         [W/m2]
     real,pointer,dimension(:) :: avg_sensible_ac   ! Sensible heat flux, ATM -> CAS         [W/m2]
     real,pointer,dimension(:,:) :: avg_sensible_gg ! Net soil heat flux between layers      [W/m2]
     real,pointer,dimension(:) :: avg_runoff_heat   ! Surface runoff internal energy flux    [W/m2]

     !----- Mass and Energy ----------------------------------------------------------------!

     real,pointer,dimension(:) :: avg_veg_energy    ! The sum of vegetation energy of cohorts in patch
                                                    ! [J/m2]   
     real,pointer,dimension(:) :: avg_veg_temp      ! The mean patch level vegetation temp  [K]
                                                    ! 
     real,pointer,dimension(:) :: avg_veg_fliq      ! The mean patch level leaf surf liquid frac 
     real,pointer,dimension(:) :: avg_veg_water     ! The sum of water residing on all leaf surface
                                                    ! [kg/m2]
     !----- Hydrology variables ------------------------------------------------------------!

     real,pointer,dimension(:)   :: watertable
     real,pointer,dimension(:)   :: moist_dz
     real,pointer,dimension(:)   :: ksat
     real,pointer,dimension(:)   :: soil_sat_energy
     real,pointer,dimension(:)   :: soil_sat_water
     real,pointer,dimension(:)   :: soil_sat_heat
     real,pointer,dimension(:,:) :: runoff_A ! Runoff parameters, 1st dimension is 3.
     real,pointer,dimension(:)   :: runoff_rate
     real,pointer,dimension(:)   :: runoff
          
  end type sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !----------------------------------------------------------------------!
  ! Polygon Type:
  ! The following are the arrays of site level variables
  ! that populate the current polygon.
  ! ----------------------------------------------------
  type polygontype

     integer :: nsites

     !  The global index of the first patch in all each site
     integer,pointer,dimension(:) :: sipa_id

     ! The number of patches in each site
     integer,pointer,dimension(:) :: sipa_n

     ! Global index of the first site in this vector across all sites
     ! on the grid
     integer :: siglob_id     ! Global index of the first patch across all cohorts

     integer,pointer,dimension(:) :: patch_count    ! number of patches per site
     
     ! Pointer to the patch vectors
     type(sitetype),pointer,dimension(:) :: site

     integer,pointer,dimension(:) :: sitenum

     ! The vectorized landuse matrix is allocated in landuse_init
     type(lutime), pointer,dimension(:,:) :: clutimes !(luyears,nsites)

     integer,pointer,dimension(:) :: num_landuse_years  !The number of years of landuse data
                                                        ! calculated at read-in of data during 
                                                        ! landuse_init

     real,   pointer,dimension(:,:) :: lai_pft    ! Site level mean LAI, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)
     real,   pointer,dimension(:,:) :: wpa_pft    ! Woody Projected Area, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)
     real,   pointer,dimension(:,:) :: wai_pft    ! Woody area index, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)

     !-----------------------------------
     ! BASIC INFO
     !-----------------------------------

     real,pointer,dimension(:) :: area            ! Proportion of a grid cell occupied by the site
     real,pointer,dimension(:) :: patch_area      ! Un-normalized sum of patch areas
     real,pointer,dimension(:) :: elevation       ! mean site elevation (meters)
     real,pointer,dimension(:) :: slope           ! mean site slope (degrees)
     real,pointer,dimension(:) :: aspect          ! mean site aspect (degrees)

     
     !--------------------------
     !  hydrologic routing info
     !--------------------------
     real,pointer,dimension(:) :: TCI             ! topographic convergence index

     integer,pointer,dimension(:) :: lsl          ! Index of the lowest soil layer

     real,pointer,dimension(:) :: pptweight       ! lapse weight for precip

     ! The Following variables used to point to structures
     ! Now they point to the array index of the site
     !   type(site), pointer :: hydro_next  !site run-off goes to
     !   type(site), pointer :: hydro_prev  !site run-on comes from
     integer, pointer,dimension(:) :: hydro_next  !site run-off goes to
     integer, pointer,dimension(:) :: hydro_prev  !site run-on comes from
     
     real,pointer,dimension(:)  :: moist_W   !Wetness index
     real,pointer,dimension(:)  :: moist_f   ! decay of soil conductance w/ depth
     real,pointer,dimension(:)  :: moist_tau ! tau characteristic time scale for water redistribution
     real,pointer,dimension(:)  :: moist_zi  ! TOPMODEL "equilibrium" water table depth 
     real,pointer,dimension(:)  :: baseflow  ! loss of water from site to watershed discharge (kg/m2/s)

     integer,pointer,dimension(:,:) :: ntext_soil  ! The soil classifications index of the soil layer
                                                   ! refer to soil_coms to view the soil parameters
                                                   ! associated with each of these classes

     !-----------------------------------
     ! TEMPERATURE VARIABLES
     !-----------------------------------
     ! minimum daily-averaged temperature in this month -- determines whether
     ! or not recruits come up this month.
     real,pointer,dimension(:)  :: min_monthly_temp

     !-----------------------------------
     ! FORESTRY
     !-----------------------------------

     ! is the site under plantation management? (1=yes, 0=no)
     integer,pointer,dimension(:) :: plantation  ! initialized to zero in site creation

     ! Upon creating an agriculture patch in this site, stock it with this 
     ! PFT.  Set, along with the other stocking parameters, in 
     ! init_ed_site_vars().
     integer,pointer,dimension(:) :: agri_stocking_pft
     
     ! Upon creating an agriculture patch in this site, stock it with 
     ! this density of plants [plants/m2]
     real,pointer,dimension(:) :: agri_stocking_density
     
     ! Upon creating a plantation patch in this site, stock it with this PFT
     integer,pointer,dimension(:) :: plantation_stocking_pft
     
     ! Upon creating an plantation patch in this site, stock it with 
     ! this density of plants [plants/m2]
     real,pointer,dimension(:) :: plantation_stocking_density

     ! Unapplied primary forest harvest from previous years (save until 
     ! harvest is above minimum threshold.) [kgC/m2].  Initialized 
     ! together with secondary memory in init_ed_site_vars().
     real,pointer,dimension(:) :: primary_harvest_memory
     
     ! Unapplied secondary forest harvest from previous years (save until 
     ! harvest is above minimum threshold.) [kgC/m2]
     real ,pointer,dimension(:):: secondary_harvest_memory


     !-----------------------------------
     ! FIRE
     !-----------------------------------

     ! site average fire disturbance rate
     real,pointer,dimension(:) :: fire_disturbance_rate
     ! total fuel in the dry patches
     real,pointer,dimension(:) :: ignition_rate

     real,pointer, dimension(:,:) :: lambda_fire ! initialized in create_site !(12,nsites)

     type(prescribed_phen),pointer, dimension(:) :: phen_pars

     !-----------------------------------
     ! TREEFALL
     !-----------------------------------
     ! site average treefall disturbance rate
     real,pointer, dimension(:) :: treefall_disturbance_rate

     !-----------------------------------
     ! DISTURBANCE
     !-----------------------------------
     ! rate of natural disturbance
     real,pointer, dimension(:) :: nat_disturbance_rate

     ! disturbance dist_type id: 
     !        dist_type = 0 if treefall was last disturbance
     !        dist_type = 1 if fire was last disturbance
     integer,pointer, dimension(:) :: nat_dist_type

     ! if new patch is less than min size, store information in the memory
     real,pointer, dimension(:,:,:) :: disturbance_memory !(n_dist_types,n_dist_types,nsites)

     ! the disturbance matrix (to,from)
     real,pointer,dimension(:,:,:) :: disturbance_rates !(n_dist_types,n_dist_types,nsites)

     ! fraction of non-surviving a.g. material removed from patch
     real,pointer,dimension(:,:) :: loss_fraction !(n_dist_types,nsites)

     real,pointer,dimension(:,:):: green_leaf_factor !(n_pft,nsites)
     real,pointer,dimension(:,:) :: leaf_aging_factor !(n_pft,nsites)

     type(met_driv_state),pointer,dimension(:) :: met

     real,pointer, dimension(:,:,:) :: basal_area !(n_pft,n_dbh,nsites), cm2/m2
     real,pointer, dimension(:,:,:) :: agb        !(n_pft,n_dbh,nsites), kgC/m2
     real,pointer, dimension(:,:,:) :: pldens     !(n_pft,n_dbh,nsites)  plant/m2
     real,pointer, dimension(:,:,:) :: bseeds     !(n_pft,n_dbh,nsites)  kgC/m2

     real,pointer, dimension(:,:,:) :: basal_area_growth ! cm2/m2/yr
     real,pointer, dimension(:,:,:) :: basal_area_mort   ! c2/m2/yr
     real,pointer, dimension(:,:,:) :: basal_area_cut    ! c2/m2/yr

     real,pointer, dimension(:,:,:) :: agb_growth        ! kgC/m2/yr
     real,pointer, dimension(:,:,:) :: agb_mort          ! kgC/m2/yr
     real,pointer, dimension(:,:,:) :: agb_cut           ! kgC/m2/yr


     real,pointer,dimension(:) :: cosaoi
     real,pointer,dimension(:) :: albedo_beam
     real,pointer,dimension(:) :: albedo_diffuse
     real,pointer,dimension(:) :: rlong_albedo     

     real,pointer,dimension(:) :: albedt
     real,pointer,dimension(:) :: rlongup
     !----- Length of day light, used to average light levels properly. -------------------!
     real, pointer, dimension(:) :: daylight

     ! ------------------------------------------
     ! Fast time flux diagnostic variables == Polygon Level
     !-------------------------------------------
     
     !----- Moisture ----------------------------
     !                                              | Description
     real,pointer,dimension(:)   :: avg_vapor_vc    ! Vegetation to canopy air latent heat flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_dew_cg      ! Dew to ground flux
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Ground to canopy air latent heat flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Throughfall
     real,pointer,dimension(:)   :: avg_intercepted ! Intercepted
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Canopy to atmosphere water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Transpiration
     real,pointer,dimension(:)   :: avg_evap        ! Evaporation
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between layers
     real,pointer,dimension(:,:) :: avg_smoist_gc   ! Trabspired soil moisture sink
     real,pointer,dimension(:)   :: avg_runoff      ! Total runoff
     real,pointer,dimension(:)   :: avg_drainage    ! Total drainage
     real,pointer,dimension(:)   :: avg_drainage_heat! Total drainage heat flux

     !----- Auxiliary Variables (user can modify to view any variable -----------------------------------------------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable

     !---- Polygon LAI, WPA, and WAI ------------------------------------------------------!
     real, pointer, dimension(:) :: lai
     real, pointer, dimension(:) :: avg_lma
     real, pointer, dimension(:) :: wpa
     real, pointer, dimension(:) :: wai

     !----- Sensible heat -------------------------------------------------------------------------------------------!
     !                                              | Description
     real,pointer,dimension(:) :: avg_sensible_vc   ! Vegetation to Canopy sensible heat flux
     real,pointer,dimension(:) :: avg_qwshed_vg     ! Internal energy of throughfall precipitation
     real,pointer,dimension(:) :: avg_qintercepted  ! Internal energy of intercepted precipitation
     real,pointer,dimension(:) :: avg_sensible_gc   ! Ground to canopy air sensible heat flux
     real,pointer,dimension(:) :: avg_sensible_ac   ! Canopy to atmosphere sensible heat flux
     real,pointer,dimension(:,:) :: avg_sensible_gg ! Net soil heat flux between layers
     real,pointer,dimension(:) :: avg_runoff_heat   ! Total runoff internal energy flux

     !----- Mass and Energy --------------------------------------------------!

     real,pointer,dimension(:) :: avg_veg_energy
     real,pointer,dimension(:) :: avg_veg_temp
     real,pointer,dimension(:) :: avg_veg_hcap
     real,pointer,dimension(:) :: avg_veg_fliq
     real,pointer,dimension(:) :: avg_veg_water
     
     real,pointer,dimension(:) :: avg_can_temp
     real,pointer,dimension(:) :: avg_can_shv
     real,pointer,dimension(:) :: avg_can_co2
     real,pointer,dimension(:) :: avg_can_rhos
     real,pointer,dimension(:) :: avg_can_prss
     real,pointer,dimension(:) :: avg_can_theta
     real,pointer,dimension(:) :: avg_can_enthalpy
     real,pointer,dimension(:) :: avg_can_depth

     real,pointer,dimension(:,:) :: avg_soil_energy
     real,pointer,dimension(:,:) :: avg_soil_water
     real,pointer,dimension(:,:) :: avg_soil_temp
     real,pointer,dimension(:,:) :: avg_soil_fracliq
     real,pointer,dimension(:)   :: avg_soil_wetness	!relative to wilting point
     
     real,pointer,dimension(:) :: avg_skin_temp
     real,pointer,dimension(:) :: avg_available_water

     real,pointer,dimension(:) :: runoff

     !-----  Phenology
     real, pointer, dimension(:) :: rad_avg

     !----- Meteorological forcing
     real,pointer,dimension(:) :: avg_atm_tmp
     real,pointer,dimension(:) :: avg_atm_shv
     real,pointer,dimension(:) :: avg_atm_prss


     !----- NACP intercomparison ---------------------------------------------!
     real,pointer,dimension(:) :: avg_sfcw_depth
     real,pointer,dimension(:) :: avg_sfcw_energy
     real,pointer,dimension(:) :: avg_sfcw_mass
     real,pointer,dimension(:) :: avg_sfcw_tempk
     real,pointer,dimension(:) :: avg_sfcw_fracliq
     real,pointer,dimension(:) :: avg_fsc
     real,pointer,dimension(:) :: avg_ssc
     real,pointer,dimension(:) :: avg_stsc
     real,pointer,dimension(:) :: avg_balive
     real,pointer,dimension(:) :: avg_bleaf
     real,pointer,dimension(:) :: avg_broot
     real,pointer,dimension(:) :: avg_bsapwood
     real,pointer,dimension(:) :: avg_bdead
     real,pointer,dimension(:) :: avg_fsn
     real,pointer,dimension(:) :: avg_msn
     real,pointer,dimension(:) :: avg_bstorage
     real,pointer,dimension(:) :: avg_bseeds

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

  end type polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !---------------------------------------------------
  ! Ed Type:
  ! The following are the arrays of polygon level variables
  ! that populate the current grid.
  !---------------------------------------------------

  type edtype

     !---- Data space indexing variables ---------------
     !
     ! These variables are the total number of each
     ! respective hierarchical type, summed across
     ! all compute nodes if it is a parallel process
     ! These values are mostly needed for output
     ! dataspace dimensioning.

     integer :: npolygons_global
     integer :: nsites_global
     integer :: npatches_global
     integer :: ncohorts_global

     ! Index offsets for total counts of cohorts, patches
     ! and sites.  If this is the nth machine writing to
     ! a file, it needs to know how many of these types
     ! had come before, so it writes data contiguously
     ! into the same file with the others

     integer :: mach_cohort_offset_index
     integer :: mach_patch_offset_index
     integer :: mach_site_offset_index
     integer :: mach_polygon_offset_index

     !
     !---------------------------------------------------


     integer :: npolygons

     integer :: pyglob_id


     !  The global index of the first site in each polygon
     integer,pointer,dimension(:) :: pysi_id

     ! The number of sites in each polygon
     integer,pointer,dimension(:) :: pysi_n

     real(kind=8),pointer,dimension(:) :: walltime_py

     real,pointer,dimension(:) :: lat

     real,pointer,dimension(:) :: lon

     integer,pointer,dimension(:) :: xatm
     
     integer,pointer,dimension(:) :: yatm

     integer,pointer,dimension(:) :: natm

     ! Sensible heat flux from the canopy to atmosphere, averaged across sites
   
     real,pointer,dimension(:) :: sensflux_py       ! of dimension npolys

     integer,pointer,dimension(:,:) :: ntext_soil   ! Soil texture classification
                                                    ! (nzg,polygon)
     
     integer,pointer,dimension(:) :: lsl            ! Layer of lowest soil
     
     ! matrix of site hydrologic adjacency
     real,pointer,dimension(:,:,:) :: site_adjacency

     real,pointer,dimension(:) :: wbar
     real,pointer,dimension(:) :: Te
     real,pointer,dimension(:) :: zbar
     real,pointer,dimension(:) :: tau
     real,pointer,dimension(:) :: sheat
     real,pointer,dimension(:) :: baseflow
     real,pointer,dimension(:) :: runoff
     real,pointer,dimension(:) :: swliq

     type(polygontype),pointer,dimension(:) :: polygon

     integer,pointer,dimension(:) :: ilon   ! index for matching met. data
     integer,pointer,dimension(:) :: ilat   ! index for matching met. data

     ! Polygon AGB (kgC/m2)
     real,pointer,dimension(:) :: total_agb

     ! Polygon basal area (cm2/m2)
     real,pointer,dimension(:) :: total_basal_area

     ! AGB accruing due to growth (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_growth

     ! AGB lost due to mortality (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_mort

     ! AGB used to generate recruits (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_recruit



     ! CHANGED THE FOLLOWING UNIT DESCRIPTORS: FROM (cm2/m2/yr) RGK 6-13-08
     !------------
     ! BASAL_AREA accruing due to growth (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_growth

     ! BASAL_AREA lost due to mortality (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_mort

     ! BASAL_AREA used to generate recruits (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_recruit
     !------------


     integer,pointer,dimension(:) :: nsites
     ! list of site numbers
     integer,pointer,dimension(:,:) :: sitenums !(max_site,npolygons)
 
     ! specification if the adjacency table was loaded from a file
     integer,pointer,dimension(:) :: load_adjacency

     !! Meteorological driver data
     type(met_driv_data),pointer,dimension(:) :: metinput
     
     !! Lapse rate transfer data
     type(met_driv_state),pointer,dimension(:) :: met, lapse

     real,pointer,dimension(:) :: cosz
     real,pointer,dimension(:) :: mean_gpp
     real,pointer,dimension(:) :: mean_precip
     real,pointer,dimension(:) :: mean_qprecip
     real,pointer,dimension(:) :: mean_netrad

     ! Total carbon (vegetation plus soil) at the beginning of budget-averaging
     ! time [kgC/m2]
     real,pointer,dimension(:) :: cbudget_initialstorage

     ! Average NEP (GPP - plant respiration - heterotrophic respiration)
     ! [kgC/m2/day], used for evaluating daily carbon budget.
     real,pointer,dimension(:) :: cbudget_nep

     ! Total nitrogen (vegetation plus soil) at the beginning of  
     ! budget-averaging time [kgN/m2]
     real,pointer,dimension(:) :: nbudget_initialstorage

     ! Polygon basal area profile [cm2/m2]
     real,pointer,dimension(:,:,:) :: basal_area !(n_pft,n_dbh,npolygons)

     ! Polygon above-ground biomass [kgC/m2]
     real,pointer,dimension(:,:,:) :: agb  !(n_pft,n_dbh,npolygons)
     ! Polygon plant density [plant/m2]
     real,pointer, dimension(:,:,:) :: pldens     !(n_pft,n_dbh,npolygons)

     ! Seed biomass [kgC/m2]
     real,pointer, dimension(:,:,:) :: bseeds     !(n_pft,n_dbh,npolygons)

     ! ------------------------------------------
     ! Diagnostic variables
     !-------------------------------------------
     
     !----- Moisture Flux ----------------------------

     real,pointer,dimension(:)   :: avg_vapor_vc    ! Vegetation to canopy air latent heat flux
     real,pointer,dimension(:)   :: avg_dew_cg      ! Dew to ground flux
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Ground to canopy air latent heat flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Throughfall
     real,pointer,dimension(:)   :: avg_intercepted ! Intercepted
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Canopy to atmosphere water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Transpiration
     real,pointer,dimension(:)   :: avg_evap        ! Evaporation
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between layers
     real,pointer,dimension(:,:) :: avg_smoist_gc   ! Trabspired soil moisture sink
     real,pointer,dimension(:)   :: avg_runoff      ! Total runoff
     real,pointer,dimension(:)   :: avg_drainage      ! Total drainage through the soil bottom
     real,pointer,dimension(:)   :: avg_drainage_heat ! Total drainage internal heat loss

     !----- Auxillary Variables (user can modify to view any variable --------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable

     !----- LAI, WPA, and WAI ------------------------------------------------!
     real,pointer,dimension(:) :: lai
     real,pointer,dimension(:) :: avg_lma
     real,pointer,dimension(:) :: wpa
     real,pointer,dimension(:) :: wai


     !----- Sensible heat Flux -----------------------------------------------!

     real,pointer,dimension(:) :: avg_sensible_vc   ! Vegetation to Canopy sensible heat flux
     real,pointer,dimension(:) :: avg_qwshed_vg     ! Internal energy of throughfall precipitation
     real,pointer,dimension(:) :: avg_qintercepted  ! Internal energy of intercepted precipitation
     real,pointer,dimension(:) :: avg_sensible_gc   ! Ground to canopy air sensible heat flux
     real,pointer,dimension(:) :: avg_sensible_ac   ! Canopy to atmosphere sensible heat flux
     real,pointer,dimension(:,:) :: avg_sensible_gg ! Net soil heat flux between layers
     real,pointer,dimension(:) :: avg_runoff_heat   ! Total runoff internal energy flux

     !----- Mass and Energy --------------------------------------------------!

     ! New variables for testing the model stability
     real,pointer,dimension(:) :: max_veg_temp
     real,pointer,dimension(:) :: min_veg_temp
     real,pointer,dimension(:) :: max_soil_temp
     real,pointer,dimension(:) :: min_soil_temp

     real,pointer,dimension(:) :: avg_veg_energy
     real,pointer,dimension(:) :: avg_veg_temp
     real,pointer,dimension(:) :: avg_veg_hcap
     real,pointer,dimension(:) :: avg_veg_fliq
     real,pointer,dimension(:) :: avg_veg_water
     
     real,pointer,dimension(:) :: avg_can_temp
     real,pointer,dimension(:) :: avg_can_shv
     real,pointer,dimension(:) :: avg_can_co2
     real,pointer,dimension(:) :: avg_can_rhos
     real,pointer,dimension(:) :: avg_can_prss
     real,pointer,dimension(:) :: avg_can_theta
     real,pointer,dimension(:) :: avg_can_enthalpy
     real,pointer,dimension(:) :: avg_can_depth

     real,pointer,dimension(:,:) :: avg_soil_energy
     real,pointer,dimension(:,:) :: avg_soil_water
     real,pointer,dimension(:,:) :: avg_soil_temp
     real,pointer,dimension(:,:) :: avg_soil_fracliq
     
     real,pointer,dimension(:) :: avg_soil_wetness
     real,pointer,dimension(:) :: avg_skin_temp
     real,pointer,dimension(:) :: avg_available_water


     real,pointer,dimension(:) :: avg_sfcw_depth     !sfcwater_depth
     real,pointer,dimension(:) :: avg_sfcw_energy    !sfcwater_mass
     real,pointer,dimension(:) :: avg_sfcw_mass      !sfcwater_mass
     real,pointer,dimension(:) :: avg_sfcw_tempk     !sfcwater_tempk
     real,pointer,dimension(:) :: avg_sfcw_fracliq   !sfcwater_fracliq
     real,pointer,dimension(:) :: avg_bdead
     real,pointer,dimension(:) :: avg_balive
     real,pointer,dimension(:) :: avg_bleaf
     real,pointer,dimension(:) :: avg_broot
     real,pointer,dimension(:) :: avg_bsapwood
     real,pointer,dimension(:) :: avg_bseeds
     real,pointer,dimension(:) :: avg_bstorage
     real,pointer,dimension(:) :: avg_fsc
     real,pointer,dimension(:) :: avg_ssc
     real,pointer,dimension(:) :: avg_stsc

     real,pointer,dimension(:) :: avg_fsn
     real,pointer,dimension(:) :: avg_msn

     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     real,pointer,dimension(:) :: Cleaf
     real,pointer,dimension(:) :: Croot
     real,pointer,dimension(:) :: Cstore
     real,pointer,dimension(:) :: Ccwd
     real,pointer,dimension(:) :: Nleaf
     real,pointer,dimension(:) :: Ndead
     real,pointer,dimension(:) :: Nroot
     real,pointer,dimension(:) :: Nstore
     real,pointer,dimension(:) :: Ncwd
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     real,pointer,dimension(:) :: Cleaf_grow
     real,pointer,dimension(:) :: Croot_grow
     real,pointer,dimension(:) :: Cdead_grow
     real,pointer,dimension(:) :: Cstore_grow
     real,pointer,dimension(:) :: Cleaf_litter_flux
     real,pointer,dimension(:) :: Croot_litter_flux
     real,pointer,dimension(:) :: Ccwd_flux
     real,pointer,dimension(:) :: Nleaf_grow
     real,pointer,dimension(:) :: Ndead_grow
     real,pointer,dimension(:) :: Nroot_grow
     real,pointer,dimension(:) :: Nstore_grow
     real,pointer,dimension(:) :: Nleaf_litter_flux
     real,pointer,dimension(:) :: Nroot_litter_flux
     real,pointer,dimension(:) :: Ncwd_flux
     real,pointer,dimension(:) :: Nbiomass_uptake
     real,pointer,dimension(:) :: Ngross_min
     real,pointer,dimension(:) :: Nnet_min

     !----- Meteorologic Conditions ----------------------------------------------------!
     real,pointer,dimension(:) :: avg_nir_beam
     real,pointer,dimension(:) :: avg_nir_diffuse
     real,pointer,dimension(:) :: avg_par_beam
     real,pointer,dimension(:) :: avg_par_diffuse
     real,pointer,dimension(:) :: avg_atm_tmp
     real,pointer,dimension(:) :: avg_atm_shv
     real,pointer,dimension(:) :: avg_rshort
     real,pointer,dimension(:) :: avg_rshort_diffuse
     real,pointer,dimension(:) :: avg_rlong
     real,pointer,dimension(:) :: avg_pcpg
     real,pointer,dimension(:) :: avg_qpcpg
     real,pointer,dimension(:) :: avg_dpcpg
     real,pointer,dimension(:) :: avg_vels
     real,pointer,dimension(:) :: avg_atm_prss
     real,pointer,dimension(:) :: avg_exner
     real,pointer,dimension(:) :: avg_geoht
     real,pointer,dimension(:) :: avg_atm_co2
     real,pointer,dimension(:) :: avg_albedt
     real,pointer,dimension(:) :: avg_rlongup
     real,pointer,dimension(:,:,:) :: avg_lai_ebalvars   ! This diagnostic partitions energy flux
                                                         ! variables into LAI regimes.
                                                         ! The matrix is arranged as follows
                                                         ! (LAI,VARIABLE,POLYGON)
                                                         ! Where LAI is (<2.0,2-4,>4)
                                                         ! Where Variable is (RNET,LHF,SHF,CANTEMP)
     real,pointer,dimension(:) :: avg_gpp
     real,pointer,dimension(:) :: avg_leaf_resp
     real,pointer,dimension(:) :: avg_root_resp
     real,pointer,dimension(:) :: avg_growth_resp
     real,pointer,dimension(:) :: avg_storage_resp
     real,pointer,dimension(:) :: avg_vleaf_resp
     real,pointer,dimension(:) :: avg_plant_resp
     real,pointer,dimension(:) :: avg_htroph_resp
     real,pointer,dimension(:) :: avg_leaf_drop
     real,pointer,dimension(:) :: avg_leaf_maintenance
     real,pointer,dimension(:) :: avg_root_maintenance


     !-------------------------------------------------------------------!
     ! These variables carry the daily mean, and are allocated only when !
     ! daily or monthly means are requested by the user. For now only    !
     ! polygon-level averages are available, and these are (almost) the  !
     ! same as the old structure, now only at the polygon level. This is !
     ! not a requirement at all, if you feel like looking at daily means !
     ! of site/patch/cohort level variables, feel free to include them.  !
     !-------------------------------------------------------------------!

     real, pointer, dimension(:)   :: dmean_pcpg           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_runoff         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_drainage       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_evap           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_transp         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_vc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_gc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_ac       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_sensible_vc    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_sensible_gc    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_sensible_ac    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_gpp            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_nep            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_plresp         ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_rh             ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_leaf_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_root_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_growth_resp    ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_storage_resp   ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_vleaf_resp     ! [ kgC/m²/yr]
     
     real, pointer, dimension(:,:) :: dmean_soil_temp      !(nzg,npolygons)
     real, pointer, dimension(:,:) :: dmean_soil_water     !(nzg,npolygons)
     real, pointer, dimension(:)   :: dmean_fs_open
     real, pointer, dimension(:)   :: dmean_fsw
     real, pointer, dimension(:)   :: dmean_fsn
     
     real, pointer, dimension(:)   :: dmean_can_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_shv       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_co2       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_rhos      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_prss      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_theta     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_veg_energy    ! (npolygons)
     real, pointer, dimension(:)   :: dmean_veg_water     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_veg_hcap      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_veg_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_shv       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_prss      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_vels      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rshort        ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rlong         ! (npolygons)

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

     !-------------------------------------------------------------------!
     !   These are variables updated at a daily basis, so they are not   !
     ! averages but they are written at the daily analysis               !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: lai_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: wpa_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: wai_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: dmean_gpp_dbh    !(n_dbh       ,npolygons)


     !-------------------------------------------------------------------!
     !   These are variables updated at a monthly basis, so they are not !
     ! averages but they are written at the monthly analysis             !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: agb_pft    !(n_pft       ,npolygons), kgC/m2
     real, pointer, dimension(:,:) :: ba_pft     !(n_pft       ,npolygons), cm2/m2
     real, pointer, dimension(:,:) :: bseeds_pft !(n_pft      ,npolygons), kgC/m2

     real, pointer, dimension(:,:) :: mmean_gpp_dbh  !(n_dbh       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_lai_pft  !(n_pft       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_wpa_pft  !(n_pft       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_wai_pft  !(n_pft       ,npolygons)

     !-------------------------------------------------------------------!
     ! These variables carry the montlhly mean, and are allocated only   !
     ! when monthly means are requested by the user. For now only        !
     ! polygon-level averages are available, and these are (almost) the  !
     ! same as the old structure, now only at the polygon level. This is !
     ! not a requirement at all, if you feel like looking at daily means !
     ! of site/patch/cohort level variables, feel free to include them.  !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:)   :: mmean_pcpg           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_runoff         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_drainage       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_evap           ! [   kg/m2/s]
     real, pointer, dimension(:)   :: mmean_transp         ! [   kg/m2/s]
     real, pointer, dimension(:)   :: mmean_vapor_vc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_gc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_ac       ! [   kg/m2/s]
     real, pointer, dimension(:)   :: mmean_sensible_vc    ! [      W/m2]
     real, pointer, dimension(:)   :: mmean_sensible_gc    ! [      W/m2]
     real, pointer, dimension(:)   :: mmean_sensible_ac    ! [      W/m2]
     real, pointer, dimension(:)   :: mmean_gpp            ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_nep            ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_plresp         ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_rh             ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_leaf_resp      ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_root_resp      ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_growth_resp    ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_storage_resp   ! [ kgC/m2/yr]
     real, pointer, dimension(:)   :: mmean_vleaf_resp     ! [ kgC/m2/yr]
     real, pointer, dimension(:,:) :: mmean_soil_temp      !(nzg,npolygons)
     real, pointer, dimension(:,:) :: mmean_soil_water     !(nzg,npolygons)
     real, pointer, dimension(:)   :: mmean_fs_open
     real, pointer, dimension(:)   :: mmean_fsw
     real, pointer, dimension(:)   :: mmean_fsn

     real, pointer, dimension(:)   :: mmean_can_temp      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_shv       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_co2       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_rhos      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_prss      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_theta     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_veg_energy    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_veg_water     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_veg_temp      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_veg_hcap      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_temp      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rlong     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rshort    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_shv       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_prss      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_vels      ! (npolygons)

     !-------------------------------------------------------------------!
     ! These variables serve two purposes. During the run they carry the !
     ! monthly square sum, and at the time of the monthly output, it     !
     ! is converted to standard deviation.                               !
     ! This needs to be stored at the history run, and thus it will have !
     ! different meaning if one looks at the -S- and -E- files.          !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:) :: stdev_gpp
     real, pointer, dimension(:) :: stdev_evap
     real, pointer, dimension(:) :: stdev_transp
     real, pointer, dimension(:) :: stdev_sensible
     real, pointer, dimension(:) :: stdev_nep
     real, pointer, dimension(:) :: stdev_rh

     
     !----- Disturbance rates. ----------------------------------------------!
     real, pointer, dimension(:,:,:) :: disturbance_rates

     !-----------------------------------------------------------------------!
     !     This variable storages the workload of this polygon, and this is  !
     ! "measured" by the total number of Runge-Kutta time steps this polygon !
     ! takes on average every day.  This is not averaged by patch, instead   ! 
     ! this is added (because having two patches is computationally more     !
     ! expensive than having one).  Since the number of patches changes by   !
     ! season, and this may be dramatically different, we store the previous !
     ! 12 months.                                                            !
     !-----------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: workload

  end type edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !----------------------------------------------!  
  !  GLOBAL VARIABLES



  !------------------------------------------------------------------------------------!
  ! These variables are allocated and assigned before the parallel distribution, so we !
  ! don't want to keep it inside the structure. In case of serial runs, we should have !
  ! offset full of zeroes and the polygons numbered from 1 to the number of polygons.  !
  ! In a non-POI parallel run, we want to allocate "mpolygons" in each node, but then  !
  ! we need to keep track of the actual polygon "ID" to write the output correctly.    !
  !------------------------------------------------------------------------------------!
  
  ! Number of polygons in each grid, for each machine. so this has a ngrids size. 
  integer, dimension(maxmach,maxgrds) :: gdpy

  ! Number of sites in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdsi

  ! Number of patches in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdpa

  ! Number of cohorts in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdco
  

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: py_off 
  
  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: si_off 

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: pa_off 

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: co_off 


  type(edtype),pointer,dimension(:) :: edgrid_g

  ! The following are swap variables used during
  ! deallocation-reallocation

  type(edtype)      :: edswap_g

  type(polygontype) :: polyswap_g

  type(sitetype)    :: siteswap_g

  type(patchtype)   :: patchswap_g
  
!------------------------------------------------------------------------------------------!
!   The following variables are for tracking the number of variables written in the output,!
! this way we avoid having the same ID used twice.                                         !
!------------------------------------------------------------------------------------------!
  integer :: nioglobal, niogrid, niopoly, niosite

!------------------------------------------------------------------------------------------!
! Logical switch that decides if the pointer tables for IO need to be updated
! The number and allocation of cohorts and patches dictates this, and they 
! change at a monthly frequency typically.
!------------------------------------------------------------------------------------------!
  logical :: filltables

  
  
contains
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  ! ===================================================
  !  Allocation subroutines of ed state variables
  ! ===================================================

  subroutine allocate_edglobals(ngrids)

    use ed_var_tables,only:num_var,vt_info,maxvars

    implicit none
    integer :: ngrids
    
    if (associated(edgrid_g))  then
       print*,"SHOULD NOT HAVE ASSOCIATED GLOBALS"
    else
       nullify(edgrid_g)
       allocate(edgrid_g(ngrids))
    endif

    !  Allocate the basic var-table structures

    allocate(num_var(ngrids))
    num_var = 0
    
    allocate(vt_info(maxvars,ngrids))
    vt_info(:,:)%first=.true.


    !  Initialize the global offsets

    edgrid_g(:)%mach_cohort_offset_index = 0
    edgrid_g(:)%mach_patch_offset_index = 0
    edgrid_g(:)%mach_site_offset_index = 0
    edgrid_g(:)%mach_polygon_offset_index = 0


    return
  end subroutine allocate_edglobals
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_edtype(cgrid,npolygons)
    
    implicit none

    integer :: npolygons
    type(edtype),target :: cgrid

    call nullify_edtype(cgrid)
    cgrid%npolygons = npolygons
    ! This if is needed for coupled runs: some nodes may receive areas exclusively over oceans, and 
    ! npolygons will be 0 in these nodes. Nothing should be allocated in these circumstances. 
    if (npolygons > 0) then
       allocate(cgrid%polygon(npolygons))
       allocate(cgrid%lat(npolygons))
       allocate(cgrid%lon(npolygons))
       allocate(cgrid%natm(npolygons))
       allocate(cgrid%xatm(npolygons))
       allocate(cgrid%yatm(npolygons))
       allocate(cgrid%ntext_soil(nzg,npolygons))
       allocate(cgrid%lsl(npolygons))
       
       allocate(cgrid%pysi_id(npolygons))
       allocate(cgrid%pysi_n(npolygons))
       allocate(cgrid%walltime_py(npolygons))
       allocate(cgrid%sensflux_py(npolygons))
       allocate(cgrid%site_adjacency(max_site,(max_site+1),npolygons))
       allocate(cgrid%wbar(npolygons))
       allocate(cgrid%Te(npolygons))
       allocate(cgrid%zbar(npolygons))
!! NOT USED       allocate(cgrid%tau(npolygons))
       allocate(cgrid%sheat(npolygons))
       allocate(cgrid%baseflow(npolygons))
       allocate(cgrid%runoff(npolygons))
       allocate(cgrid%swliq(npolygons))
       
       allocate(cgrid%ilon(npolygons))
       allocate(cgrid%ilat(npolygons))
       allocate(cgrid%total_agb(npolygons))
       allocate(cgrid%total_basal_area(npolygons))
       allocate(cgrid%total_agb_growth(npolygons))
       allocate(cgrid%total_agb_mort(npolygons))
       allocate(cgrid%total_agb_recruit(npolygons))
       allocate(cgrid%total_basal_area_growth(npolygons))
       allocate(cgrid%total_basal_area_mort(npolygons))
       allocate(cgrid%total_basal_area_recruit(npolygons))
       allocate(cgrid%nsites(npolygons))
       allocate(cgrid%sitenums(max_site,npolygons))
       allocate(cgrid%load_adjacency(npolygons))
       allocate(cgrid%cosz(npolygons))
       allocate(cgrid%mean_gpp(npolygons))
       allocate(cgrid%mean_precip(npolygons))
       allocate(cgrid%mean_qprecip(npolygons))
       allocate(cgrid%mean_netrad(npolygons))
       allocate(cgrid%cbudget_initialstorage(npolygons))
       allocate(cgrid%cbudget_nep(npolygons))
       allocate(cgrid%nbudget_initialstorage(npolygons))
       allocate(cgrid%basal_area(n_pft,n_dbh,npolygons))
       allocate(cgrid%agb(n_pft,n_dbh,npolygons))
       allocate(cgrid%pldens(n_pft,n_dbh,npolygons))
       allocate(cgrid%bseeds(n_pft,n_dbh,npolygons))
       
       allocate(cgrid%metinput(npolygons))
       allocate(cgrid%met(npolygons))

       allocate(cgrid%lapse(npolygons))

       allocate(cgrid%lai  (npolygons))
       allocate(cgrid%avg_lma(npolygons))
       allocate(cgrid%wpa  (npolygons))
       allocate(cgrid%wai  (npolygons))

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       allocate(cgrid%avg_vapor_vc  (npolygons))
       allocate(cgrid%avg_dew_cg    (npolygons))
       allocate(cgrid%avg_vapor_gc  (npolygons))
       allocate(cgrid%avg_wshed_vg  (npolygons))
       allocate(cgrid%avg_intercepted (npolygons))
       allocate(cgrid%avg_vapor_ac  (npolygons))
       allocate(cgrid%avg_transp    (npolygons))
       allocate(cgrid%avg_evap      (npolygons))
       allocate(cgrid%avg_smoist_gg (nzg,npolygons))
       allocate(cgrid%avg_smoist_gc (nzg,npolygons))
       allocate(cgrid%avg_runoff        (npolygons))
       allocate(cgrid%avg_drainage       (npolygons))
       allocate(cgrid%avg_drainage_heat  (npolygons))
       allocate(cgrid%aux           (npolygons))
       allocate(cgrid%aux_s         (nzg,npolygons))
       allocate(cgrid%avg_sensible_vc  (npolygons))
       allocate(cgrid%avg_qwshed_vg    (npolygons))
       allocate(cgrid%avg_qintercepted (npolygons))
       allocate(cgrid%avg_sensible_gc  (npolygons))
       allocate(cgrid%avg_sensible_ac  (npolygons))
       allocate(cgrid%avg_sensible_gg  (nzg,npolygons))
       allocate(cgrid%avg_runoff_heat  (npolygons))

       allocate(cgrid%max_veg_temp(npolygons))
       allocate(cgrid%min_veg_temp(npolygons))
       allocate(cgrid%max_soil_temp(npolygons))
       allocate(cgrid%min_soil_temp(npolygons))

       ! Fast time state diagnostics
       allocate(cgrid%avg_veg_energy  (npolygons))
       allocate(cgrid%avg_veg_hcap    (npolygons))
       allocate(cgrid%avg_veg_temp    (npolygons))
       allocate(cgrid%avg_veg_fliq    (npolygons))
       allocate(cgrid%avg_veg_water   (npolygons))
       allocate(cgrid%avg_can_temp    (npolygons))
       allocate(cgrid%avg_can_shv     (npolygons))
       allocate(cgrid%avg_can_co2     (npolygons))
       allocate(cgrid%avg_can_rhos    (npolygons))
       allocate(cgrid%avg_can_prss    (npolygons))
       allocate(cgrid%avg_can_theta   (npolygons))
       allocate(cgrid%avg_can_enthalpy(npolygons))
       allocate(cgrid%avg_can_depth   (npolygons))
       allocate(cgrid%avg_soil_energy(nzg,npolygons))
       allocate(cgrid%avg_soil_water(nzg,npolygons))
       allocate(cgrid%avg_soil_temp (nzg,npolygons))
       allocate(cgrid%avg_soil_fracliq (nzg,npolygons))

       allocate(cgrid%avg_soil_wetness   (npolygons))
       allocate(cgrid%avg_skin_temp      (npolygons))
       allocate(cgrid%avg_available_water(npolygons))

       allocate(cgrid%avg_lai_ebalvars (3,4,npolygons))

       allocate(cgrid%avg_gpp  (npolygons))
       allocate(cgrid%avg_leaf_resp  (npolygons))
       allocate(cgrid%avg_root_resp  (npolygons))
       allocate(cgrid%avg_growth_resp  (npolygons))
       allocate(cgrid%avg_storage_resp  (npolygons))
       allocate(cgrid%avg_vleaf_resp  (npolygons))
       allocate(cgrid%avg_plant_resp  (npolygons))
       allocate(cgrid%avg_htroph_resp  (npolygons))
       allocate(cgrid%avg_leaf_drop       (npolygons))
       allocate(cgrid%avg_leaf_maintenance(npolygons))
       allocate(cgrid%avg_root_maintenance(npolygons))

       !! added MCD for NACP intercomparison
       allocate(cgrid%avg_sfcw_depth   (npolygons))
       allocate(cgrid%avg_sfcw_energy  (npolygons))
       allocate(cgrid%avg_sfcw_mass    (npolygons))
       allocate(cgrid%avg_sfcw_tempk   (npolygons))
       allocate(cgrid%avg_sfcw_fracliq (npolygons))
       allocate(cgrid%avg_bdead       (npolygons))
       allocate(cgrid%avg_balive      (npolygons))
       allocate(cgrid%avg_bleaf       (npolygons))
       allocate(cgrid%avg_broot       (npolygons))
       allocate(cgrid%avg_bsapwood    (npolygons))
       allocate(cgrid%avg_bstorage    (npolygons))
       allocate(cgrid%avg_bseeds      (npolygons))
       allocate(cgrid%avg_fsc         (npolygons))
       allocate(cgrid%avg_ssc         (npolygons))
       allocate(cgrid%avg_stsc        (npolygons))
       allocate(cgrid%avg_fsn         (npolygons))
       allocate(cgrid%avg_msn         (npolygons))

       !! C & N fluxes and pools for NCEAS/FACE
       allocate(cgrid%Cleaf  (npolygons))
     allocate(cgrid%Croot  (npolygons))
     allocate(cgrid%Cstore (npolygons))
     allocate(cgrid%Ccwd   (npolygons))
     allocate(cgrid%Nleaf  (npolygons))
     allocate(cgrid%Ndead  (npolygons))
     allocate(cgrid%Nroot  (npolygons))
     allocate(cgrid%Nstore (npolygons))
     allocate(cgrid%Ncwd   (npolygons))
     allocate(cgrid%Cleaf_grow       (npolygons))
     allocate(cgrid%Croot_grow       (npolygons))
     allocate(cgrid%Cdead_grow       (npolygons))
     allocate(cgrid%Cstore_grow      (npolygons))
     allocate(cgrid%Cleaf_litter_flux(npolygons))
     allocate(cgrid%Croot_litter_flux(npolygons))
     allocate(cgrid%Ccwd_flux        (npolygons))
     allocate(cgrid%Nleaf_grow       (npolygons))
     allocate(cgrid%Ndead_grow       (npolygons))
     allocate(cgrid%Nroot_grow       (npolygons))
     allocate(cgrid%Nstore_grow      (npolygons))
     allocate(cgrid%Nleaf_litter_flux(npolygons))
     allocate(cgrid%Nroot_litter_flux(npolygons))
     allocate(cgrid%Ncwd_flux        (npolygons))
     allocate(cgrid%Nbiomass_uptake  (npolygons))
     allocate(cgrid%Ngross_min       (npolygons))
     allocate(cgrid%Nnet_min         (npolygons))



       ! Meteorologic conditions (forcing)
       allocate(cgrid%avg_nir_beam     (npolygons))
       allocate(cgrid%avg_nir_diffuse  (npolygons))
       allocate(cgrid%avg_par_beam     (npolygons))
       allocate(cgrid%avg_par_diffuse  (npolygons))
       allocate(cgrid%avg_atm_tmp      (npolygons))
       allocate(cgrid%avg_atm_shv      (npolygons))
       allocate(cgrid%avg_rshort       (npolygons))
       allocate(cgrid%avg_rshort_diffuse (npolygons))
       allocate(cgrid%avg_rlong        (npolygons))
       allocate(cgrid%avg_pcpg         (npolygons))
       allocate(cgrid%avg_qpcpg        (npolygons))
       allocate(cgrid%avg_dpcpg        (npolygons))
       allocate(cgrid%avg_vels         (npolygons))
       allocate(cgrid%avg_atm_prss     (npolygons))
       allocate(cgrid%avg_exner        (npolygons))
       allocate(cgrid%avg_geoht        (npolygons))
       allocate(cgrid%avg_atm_co2      (npolygons))
       allocate(cgrid%avg_albedt       (npolygons))
       allocate(cgrid%avg_rlongup      (npolygons))
       
       allocate(cgrid%lai_pft            (n_pft       ,npolygons))
       allocate(cgrid%wpa_pft            (n_pft       ,npolygons))
       allocate(cgrid%wai_pft            (n_pft       ,npolygons))

          allocate(cgrid%workload             (13          ,npolygons))

       !-----------------------------------------------------------------!
       ! Allocating the daily means, only if daily or monthly means were !
       ! requested by the user.                                          !
       !-----------------------------------------------------------------!
       if (idoutput > 0 .or. imoutput > 0) then
          allocate(cgrid%dmean_pcpg           (             npolygons))
          allocate(cgrid%dmean_runoff         (             npolygons))
          allocate(cgrid%dmean_drainage       (             npolygons))
          allocate(cgrid%dmean_vapor_ac       (             npolygons))
          allocate(cgrid%dmean_vapor_gc       (             npolygons))
          allocate(cgrid%dmean_vapor_vc       (             npolygons))
          allocate(cgrid%dmean_gpp            (             npolygons))
          allocate(cgrid%dmean_evap           (             npolygons))
          allocate(cgrid%dmean_transp         (             npolygons))
          allocate(cgrid%dmean_sensible_ac    (             npolygons))
          allocate(cgrid%dmean_sensible_gc    (             npolygons))
          allocate(cgrid%dmean_sensible_vc    (             npolygons))
          allocate(cgrid%dmean_plresp         (             npolygons))
          allocate(cgrid%dmean_rh             (             npolygons))
          allocate(cgrid%dmean_leaf_resp      (             npolygons))
          allocate(cgrid%dmean_root_resp      (             npolygons))
          allocate(cgrid%dmean_growth_resp    (             npolygons))
          allocate(cgrid%dmean_storage_resp   (             npolygons))
          allocate(cgrid%dmean_vleaf_resp     (             npolygons))
          allocate(cgrid%dmean_nep            (             npolygons))
          allocate(cgrid%dmean_soil_temp      (nzg,         npolygons))
          allocate(cgrid%dmean_soil_water     (nzg,         npolygons))
          allocate(cgrid%dmean_fs_open        (             npolygons))
          allocate(cgrid%dmean_fsw            (             npolygons))
          allocate(cgrid%dmean_fsn            (             npolygons))
          allocate(cgrid%dmean_gpp_dbh        (n_dbh       ,npolygons))
          allocate(cgrid%dmean_can_temp       (             npolygons))
          allocate(cgrid%dmean_can_shv        (             npolygons))
          allocate(cgrid%dmean_can_co2        (             npolygons))
          allocate(cgrid%dmean_can_rhos       (             npolygons))
          allocate(cgrid%dmean_can_prss       (             npolygons))
          allocate(cgrid%dmean_can_theta      (             npolygons))
          allocate(cgrid%dmean_veg_energy     (             npolygons))
          allocate(cgrid%dmean_veg_water      (             npolygons))
          allocate(cgrid%dmean_veg_hcap       (             npolygons))
          allocate(cgrid%dmean_veg_temp       (             npolygons))
          allocate(cgrid%dmean_atm_temp       (             npolygons))
          allocate(cgrid%dmean_atm_shv        (             npolygons))
          allocate(cgrid%dmean_atm_prss       (             npolygons))
          allocate(cgrid%dmean_atm_vels       (             npolygons))
          allocate(cgrid%dmean_rshort         (             npolygons))
          allocate(cgrid%dmean_rlong          (             npolygons))

          allocate(cgrid%dmean_co2_residual   (             npolygons))
          allocate(cgrid%dmean_energy_residual(             npolygons))
          allocate(cgrid%dmean_water_residual (             npolygons))

       end if
       !-------------------------------------------------------------------!
       ! Allocating the monthly means, only if monthly means were          !
       ! requested by the user.                                            !
       !-------------------------------------------------------------------!
       if (imoutput > 0) then
          allocate(cgrid%mmean_pcpg           (             npolygons))
          allocate(cgrid%mmean_runoff         (             npolygons))
          allocate(cgrid%mmean_drainage       (             npolygons))
          allocate(cgrid%mmean_evap           (             npolygons))
          allocate(cgrid%mmean_transp         (             npolygons))
          allocate(cgrid%mmean_vapor_ac       (             npolygons))
          allocate(cgrid%mmean_vapor_gc       (             npolygons))
          allocate(cgrid%mmean_vapor_vc       (             npolygons))
          allocate(cgrid%mmean_sensible_ac    (             npolygons))
          allocate(cgrid%mmean_sensible_gc    (             npolygons))
          allocate(cgrid%mmean_sensible_vc    (             npolygons))
          allocate(cgrid%mmean_gpp            (             npolygons))
          allocate(cgrid%mmean_nep            (             npolygons))
          allocate(cgrid%mmean_plresp         (             npolygons))
          allocate(cgrid%mmean_rh             (             npolygons))
          allocate(cgrid%mmean_leaf_resp      (             npolygons))
          allocate(cgrid%mmean_root_resp      (             npolygons))
          allocate(cgrid%mmean_growth_resp    (             npolygons))
          allocate(cgrid%mmean_storage_resp   (             npolygons))
          allocate(cgrid%mmean_vleaf_resp     (             npolygons))
          allocate(cgrid%mmean_soil_temp      (nzg,         npolygons))
          allocate(cgrid%mmean_soil_water     (nzg,         npolygons))
          allocate(cgrid%mmean_fs_open        (             npolygons))
          allocate(cgrid%mmean_fsw            (             npolygons))
          allocate(cgrid%mmean_fsn            (             npolygons))
          allocate(cgrid%mmean_gpp_dbh        (n_dbh       ,npolygons))
          allocate(cgrid%mmean_lai_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_wpa_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_wai_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_can_temp       (             npolygons))
          allocate(cgrid%mmean_can_shv        (             npolygons))
          allocate(cgrid%mmean_can_co2        (             npolygons))
          allocate(cgrid%mmean_can_rhos       (             npolygons))
          allocate(cgrid%mmean_can_prss       (             npolygons))
          allocate(cgrid%mmean_can_theta      (             npolygons))
          allocate(cgrid%mmean_veg_energy     (             npolygons))
          allocate(cgrid%mmean_veg_water      (             npolygons))
          allocate(cgrid%mmean_veg_temp       (             npolygons))
          allocate(cgrid%mmean_veg_hcap       (             npolygons))
          allocate(cgrid%mmean_atm_temp       (             npolygons))
          allocate(cgrid%mmean_rshort         (             npolygons))
          allocate(cgrid%mmean_rlong          (             npolygons))
          allocate(cgrid%mmean_atm_shv        (             npolygons))
          allocate(cgrid%mmean_atm_prss       (             npolygons))
          allocate(cgrid%mmean_atm_vels       (             npolygons))
          allocate(cgrid%mmean_co2_residual   (             npolygons))
          allocate(cgrid%mmean_energy_residual(             npolygons))
          allocate(cgrid%mmean_water_residual (             npolygons))
          allocate(cgrid%bseeds_pft           (n_pft       ,npolygons))
          allocate(cgrid%agb_pft              (n_pft       ,npolygons))
          allocate(cgrid%ba_pft               (n_pft       ,npolygons))
          allocate(cgrid%stdev_gpp            (             npolygons))
          allocate(cgrid%stdev_evap           (             npolygons))
          allocate(cgrid%stdev_transp         (             npolygons))
          allocate(cgrid%stdev_sensible       (             npolygons))
          allocate(cgrid%stdev_nep            (             npolygons))
          allocate(cgrid%stdev_rh             (             npolygons))
          allocate(cgrid%disturbance_rates    (n_dist_types,n_dist_types,npolygons))

       end if



       ! Initialize the variables with a non-sense number.
       !call huge_edtype(cgrid)
    end if
    return
  end subroutine allocate_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_polygontype(cpoly,nsites)

    implicit none

    integer :: nsites
    type(polygontype),target :: cpoly

    call nullify_polygontype(cpoly)

    cpoly%nsites = nsites

    allocate(cpoly%sipa_id(nsites))
    allocate(cpoly%sipa_n(nsites))

    allocate(cpoly%patch_count(nsites))  
    allocate(cpoly%site(nsites))
    allocate(cpoly%sitenum(nsites))

    allocate(cpoly%lsl(nsites))   

    allocate(cpoly%area(nsites))
    allocate(cpoly%patch_area(nsites))
    allocate(cpoly%elevation(nsites))
    allocate(cpoly%slope(nsites))
    allocate(cpoly%aspect(nsites))
    
    allocate(cpoly%num_landuse_years(nsites))

    allocate(cpoly%lai_pft(n_pft,nsites))
    allocate(cpoly%wpa_pft(n_pft,nsites))
    allocate(cpoly%wai_pft(n_pft,nsites))

    allocate(cpoly%TCI(nsites))      
    allocate(cpoly%pptweight(nsites))      
    allocate(cpoly%hydro_next(nsites))
    allocate(cpoly%hydro_prev(nsites))
    allocate(cpoly%moist_W(nsites))
    allocate(cpoly%moist_f(nsites))  
    allocate(cpoly%moist_tau(nsites))
    allocate(cpoly%moist_zi(nsites)) 
    allocate(cpoly%baseflow(nsites)) 
    allocate(cpoly%ntext_soil(nzg,nsites))
    allocate(cpoly%min_monthly_temp(nsites))
    allocate(cpoly%plantation(nsites)) 
    allocate(cpoly%agri_stocking_pft(nsites))
    allocate(cpoly%agri_stocking_density(nsites))
    allocate(cpoly%plantation_stocking_pft(nsites))
    allocate(cpoly%plantation_stocking_density(nsites))
    allocate(cpoly%primary_harvest_memory(nsites))
    allocate(cpoly%secondary_harvest_memory(nsites))
    allocate(cpoly%fire_disturbance_rate(nsites))
    allocate(cpoly%ignition_rate(nsites))
    allocate(cpoly%lambda_fire(12,nsites))
    allocate(cpoly%phen_pars(nsites))!THIS PTR IS ALLOCATED IN PHENOLOGY_INIT
    allocate(cpoly%treefall_disturbance_rate(nsites))
    allocate(cpoly%nat_disturbance_rate(nsites))
    allocate(cpoly%nat_dist_type(nsites))
    allocate(cpoly%disturbance_memory(n_dist_types,n_dist_types,nsites))
    allocate(cpoly%disturbance_rates(n_dist_types,n_dist_types,nsites))
    allocate(cpoly%loss_fraction(n_dist_types,nsites))

    allocate(cpoly%green_leaf_factor(n_pft,nsites))
    allocate(cpoly%leaf_aging_factor(n_pft,nsites))
    
    allocate(cpoly%met(nsites))

    allocate(cpoly%basal_area  (n_pft,n_dbh,nsites))
    allocate(cpoly%agb         (n_pft,n_dbh,nsites))
    allocate(cpoly%pldens      (n_pft,n_dbh,nsites))
    allocate(cpoly%bseeds      (n_pft,n_dbh,nsites))

    allocate(cpoly%basal_area_growth (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_growth        (n_pft,n_dbh,nsites))
    allocate(cpoly%basal_area_mort   (n_pft,n_dbh,nsites))
    allocate(cpoly%basal_area_cut   (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_mort          (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_cut          (n_pft,n_dbh,nsites))

    allocate(cpoly%cosaoi(nsites))
    allocate(cpoly%albedo_beam(nsites))
    allocate(cpoly%albedo_diffuse(nsites))
    allocate(cpoly%rlong_albedo(nsites))
    
    allocate(cpoly%albedt(nsites))
    allocate(cpoly%rlongup(nsites))
    allocate(cpoly%daylight(nsites))

    allocate(cpoly%lai  (nsites)) 
    allocate(cpoly%avg_lma(nsites))
    allocate(cpoly%wpa  (nsites))
    allocate(cpoly%wai  (nsites))
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    allocate(cpoly%avg_vapor_vc  (nsites))
    allocate(cpoly%avg_dew_cg    (nsites))
    allocate(cpoly%avg_vapor_gc  (nsites))
    allocate(cpoly%avg_wshed_vg  (nsites))
    allocate(cpoly%avg_intercepted (nsites))
    allocate(cpoly%avg_vapor_ac  (nsites))
    allocate(cpoly%avg_transp    (nsites))
    allocate(cpoly%avg_evap      (nsites))
    
    allocate(cpoly%avg_smoist_gg (nzg,nsites))
    allocate(cpoly%avg_smoist_gc (nzg,nsites))
    allocate(cpoly%avg_runoff        (nsites))
    allocate(cpoly%avg_drainage  (nsites))
    allocate(cpoly%avg_drainage_heat  (nsites))
    allocate(cpoly%aux           (nsites))
    allocate(cpoly%aux_s         (nzg,nsites))
    allocate(cpoly%avg_sensible_vc  (nsites))
    allocate(cpoly%avg_qwshed_vg    (nsites))
    allocate(cpoly%avg_qintercepted (nsites))
    allocate(cpoly%avg_sensible_gc  (nsites))
    allocate(cpoly%avg_sensible_ac  (nsites))
    allocate(cpoly%avg_sensible_gg  (nzg,nsites))
    allocate(cpoly%avg_runoff_heat         (nsites))

    ! Fast time state diagnostics
    allocate(cpoly%avg_veg_energy          (nsites))
    allocate(cpoly%avg_veg_hcap            (nsites))
    allocate(cpoly%avg_veg_temp            (nsites))
    allocate(cpoly%avg_veg_fliq            (nsites))
    allocate(cpoly%avg_veg_water           (nsites))
    allocate(cpoly%avg_can_temp            (nsites))
    allocate(cpoly%avg_can_shv             (nsites))
    allocate(cpoly%avg_can_co2             (nsites))
    allocate(cpoly%avg_can_rhos            (nsites))
    allocate(cpoly%avg_can_prss            (nsites))
    allocate(cpoly%avg_can_theta           (nsites))
    allocate(cpoly%avg_can_enthalpy        (nsites))
    allocate(cpoly%avg_can_depth           (nsites))
    allocate(cpoly%avg_soil_energy     (nzg,nsites))
    allocate(cpoly%avg_soil_water      (nzg,nsites))
    allocate(cpoly%avg_soil_temp       (nzg,nsites))
    allocate(cpoly%avg_soil_fracliq    (nzg,nsites))
    allocate(cpoly%avg_soil_wetness        (nsites))
    allocate(cpoly%avg_skin_temp           (nsites))
    allocate(cpoly%avg_available_water     (nsites))
    allocate(cpoly%runoff                  (nsites))
    ! Phenology-related
    allocate(cpoly%rad_avg                 (nsites))
    ! Meteorological data
    allocate(cpoly%avg_atm_tmp             (nsites))
    allocate(cpoly%avg_atm_shv             (nsites))
    allocate(cpoly%avg_atm_prss            (nsites))

    !!!NACP
    allocate(cpoly%avg_sfcw_depth           (nsites))
    allocate(cpoly%avg_sfcw_energy          (nsites))
    allocate(cpoly%avg_sfcw_mass            (nsites))
    allocate(cpoly%avg_sfcw_fracliq         (nsites))
    allocate(cpoly%avg_sfcw_tempk           (nsites))
    allocate(cpoly%avg_fsc                 (nsites))
    allocate(cpoly%avg_stsc                (nsites))
    allocate(cpoly%avg_ssc                 (nsites))
    allocate(cpoly%avg_balive              (nsites))
    allocate(cpoly%avg_bleaf               (nsites))
    allocate(cpoly%avg_broot               (nsites))
    allocate(cpoly%avg_bsapwood            (nsites))
    allocate(cpoly%avg_bdead               (nsites))
    allocate(cpoly%avg_fsn                 (nsites))
    allocate(cpoly%avg_msn                 (nsites))


    allocate(cpoly%avg_bstorage            (nsites))
    allocate(cpoly%avg_bseeds              (nsites))


    allocate(cpoly%dmean_co2_residual      (nsites))
    allocate(cpoly%dmean_energy_residual   (nsites))
    allocate(cpoly%dmean_water_residual    (nsites))
    allocate(cpoly%mmean_co2_residual      (nsites))
    allocate(cpoly%mmean_energy_residual   (nsites))
    allocate(cpoly%mmean_water_residual    (nsites))


    ! Initialize the variables with a non-sense number.
    !call huge_polygontype(cpoly)
    
    return
  end subroutine allocate_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_sitetype(csite,npatches)

    implicit none

    integer :: npatches,ipa
    type(sitetype),target :: csite
    
    call nullify_sitetype(csite)
    csite%npatches = npatches
    allocate(csite%paco_id(npatches))
    allocate(csite%paco_n(npatches))
    allocate(csite%patch(npatches))

    ! Initialize zero cohorts
    do ipa=1,npatches
       csite%patch(ipa)%ncohorts = 0
    enddo

    allocate(csite%lai(npatches))
    allocate(csite%wpa(npatches))
    allocate(csite%wai(npatches))
    allocate(csite%dist_type(npatches))
    allocate(csite%age(npatches))
    allocate(csite%area(npatches))
    allocate(csite%laiarea(npatches))
    allocate(csite%hcapveg(npatches))
    allocate(csite%fast_soil_C(npatches))
    allocate(csite%slow_soil_C(npatches))
    allocate(csite%structural_soil_C(npatches))
    allocate(csite%structural_soil_L(npatches))
    allocate(csite%mineralized_soil_N(npatches))
    allocate(csite%fast_soil_N(npatches))
    allocate(csite%sum_dgd(npatches))
    allocate(csite%sum_chd(npatches))
    allocate(csite%plantation(npatches))
    allocate(csite%can_enthalpy(npatches))
    allocate(csite%can_temp(npatches))
    allocate(csite%can_shv(npatches))
    allocate(csite%can_co2(npatches))
    allocate(csite%can_rhos(npatches))
    allocate(csite%can_prss(npatches))
    allocate(csite%can_theta(npatches))
    allocate(csite%can_depth(npatches))
    allocate(csite%lambda_light(npatches))
    allocate(csite%cohort_count(npatches))
    allocate(csite%pname(npatches))
    
    allocate(csite%sfcwater_mass(nzs,npatches))
    allocate(csite%sfcwater_energy(nzs,npatches))
    allocate(csite%sfcwater_depth(nzs,npatches))
    allocate(csite%rshort_s(nzs,npatches))
    allocate(csite%rshort_s_beam(nzs,npatches))
    allocate(csite%rshort_s_diffuse(nzs,npatches))
    allocate(csite%sfcwater_tempk(nzs,npatches))
    allocate(csite%sfcwater_fracliq(nzs,npatches))
    allocate(csite%nlev_sfcwater(npatches))
    allocate(csite%ntext_soil(nzg,npatches))
    allocate(csite%soil_energy(nzg,npatches))
    allocate(csite%soil_water(nzg,npatches))
    allocate(csite%soil_tempk(nzg,npatches))
    allocate(csite%soil_fracliq(nzg,npatches))
    allocate(csite%ground_shv(npatches))
    allocate(csite%surface_ssh(npatches))
    allocate(csite%rough(npatches))
    allocate(csite%A_o_max(n_pft,npatches)) 
    allocate(csite%A_c_max(n_pft,npatches)) 

    allocate(csite%old_stoma_data_max(n_pft,npatches))
    allocate(csite%old_stoma_vector_max(n_stoma_atts,n_pft,npatches))

    allocate(csite%avg_daily_temp(npatches))  
    allocate(csite%mean_rh(npatches))
    allocate(csite%mean_nep(npatches))
    allocate(csite%wbudget_loss2atm(npatches))
    allocate(csite%wbudget_denseffect(npatches))
    allocate(csite%wbudget_precipgain(npatches))
    allocate(csite%wbudget_loss2runoff(npatches))
    allocate(csite%wbudget_loss2drainage(npatches))
    allocate(csite%wbudget_initialstorage(npatches))
    allocate(csite%wbudget_residual(npatches))
    allocate(csite%ebudget_latent(npatches))
    allocate(csite%ebudget_loss2atm(npatches))
    allocate(csite%ebudget_denseffect(npatches))
    allocate(csite%ebudget_loss2runoff(npatches))
    allocate(csite%ebudget_loss2drainage(npatches))
    allocate(csite%ebudget_netrad(npatches))
    allocate(csite%ebudget_precipgain(npatches))
    allocate(csite%ebudget_initialstorage(npatches))
    allocate(csite%ebudget_residual(npatches))
    allocate(csite%co2budget_initialstorage(npatches))
    allocate(csite%co2budget_residual(npatches))
    allocate(csite%co2budget_loss2atm(npatches))
    allocate(csite%co2budget_denseffect(npatches))
    allocate(csite%co2budget_gpp(npatches))
    allocate(csite%co2budget_gpp_dbh(n_dbh,npatches))
    allocate(csite%co2budget_plresp(npatches))
    allocate(csite%co2budget_rh(npatches))
    allocate(csite%today_A_decomp(npatches))
    allocate(csite%today_Af_decomp(npatches))
    allocate(csite%repro(n_pft,npatches))
    allocate(csite%veg_rough(npatches))
    allocate(csite%veg_height (npatches))
    allocate(csite%fsc_in(npatches))
    allocate(csite%ssc_in(npatches))
    allocate(csite%ssl_in(npatches))
    allocate(csite%fsn_in(npatches))
    allocate(csite%total_plant_nitrogen_uptake(npatches))
    allocate(csite%mineralized_N_loss(npatches))
    allocate(csite%mineralized_N_input(npatches))
    allocate(csite%rshort_g(npatches))
    allocate(csite%rshort_g_beam(npatches))
    allocate(csite%rshort_g_diffuse(npatches))
    allocate(csite%rlong_g(npatches))
    allocate(csite%rlong_g_surf(npatches))
    allocate(csite%rlong_g_incid(npatches))
    allocate(csite%rlong_s(npatches))
    allocate(csite%rlong_s_surf(npatches))
    allocate(csite%rlong_s_incid(npatches))
    allocate(csite%albedt(npatches))
    allocate(csite%albedo_beam(npatches))
    allocate(csite%albedo_diffuse(npatches))
    allocate(csite%rlongup(npatches))
    allocate(csite%rlong_albedo(npatches))
    allocate(csite%total_snow_depth(npatches))
    allocate(csite%snowfac(npatches))
    allocate(csite%A_decomp(npatches))
    allocate(csite%f_decomp(npatches))
    allocate(csite%rh(npatches))
    allocate(csite%cwd_rh(npatches))
    allocate(csite%fuse_flag(npatches))
    allocate(csite%pft_density_profile(n_pft,ff_ndbh,npatches))
    allocate(csite%plant_ag_biomass(npatches))

    allocate(csite%mean_wflux(npatches))
    allocate(csite%mean_latflux(npatches))
    allocate(csite%mean_hflux(npatches))
    allocate(csite%mean_runoff(npatches))
    allocate(csite%mean_qrunoff(npatches))

    allocate(csite%htry       (npatches))
    allocate(csite%avg_rk4step(npatches))

    allocate(csite%avg_available_water(npatches))

    allocate(csite%ustar(npatches))
    allocate(csite%tstar(npatches))
    allocate(csite%qstar(npatches))
    allocate(csite%cstar(npatches))
    
    allocate(csite%upwp(npatches))
    allocate(csite%qpwp(npatches))
    allocate(csite%cpwp(npatches))
    allocate(csite%tpwp(npatches))
    allocate(csite%wpwp(npatches))

    allocate(csite%avg_carbon_ac(npatches))

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    allocate(csite%avg_vapor_vc  (npatches))
    allocate(csite%avg_dew_cg    (npatches))
    allocate(csite%avg_vapor_gc  (npatches))
    allocate(csite%avg_wshed_vg  (npatches))
    allocate(csite%avg_intercepted (npatches))
    allocate(csite%avg_vapor_ac  (npatches))
    allocate(csite%avg_transp    (npatches))
    allocate(csite%avg_evap      (npatches))
    allocate(csite%avg_netrad    (npatches))
    allocate(csite%avg_smoist_gg (nzg,npatches))
    allocate(csite%avg_smoist_gc (nzg,npatches))
    allocate(csite%avg_runoff    (npatches))
    allocate(csite%avg_drainage  (npatches))
    allocate(csite%avg_drainage_heat  (npatches))
    allocate(csite%aux           (npatches))
    allocate(csite%aux_s         (nzg,npatches))
    allocate(csite%avg_sensible_vc  (npatches))
    allocate(csite%avg_qwshed_vg    (npatches))
    allocate(csite%avg_qintercepted (npatches))
    allocate(csite%avg_sensible_gc  (npatches))
    allocate(csite%avg_sensible_ac  (npatches))
    allocate(csite%avg_sensible_gg  (nzg,npatches))
    allocate(csite%avg_runoff_heat  (npatches))

    ! ----------------------------------------------

    allocate(csite%avg_veg_energy(npatches))
    allocate(csite%avg_veg_temp  (npatches))
    allocate(csite%avg_veg_fliq  (npatches))
    allocate(csite%avg_veg_water (npatches))


    allocate(csite%watertable      (npatches))
    allocate(csite%moist_dz        (npatches))
    allocate(csite%ksat            (npatches))
    allocate(csite%soil_sat_energy (npatches))
    allocate(csite%soil_sat_water  (npatches))
    allocate(csite%soil_sat_heat   (npatches))

    allocate(csite%runoff_A        (3,npatches))
    allocate(csite%runoff_rate     (npatches))
    allocate(csite%runoff          (npatches))

    if (imoutput > 0 .or. idoutput > 0) then
       allocate(csite%dmean_rk4step           (npatches))
       allocate(csite%dmean_lambda_light      (npatches))
       allocate(csite%dmean_co2_residual      (npatches))
       allocate(csite%dmean_energy_residual   (npatches))
       allocate(csite%dmean_water_residual    (npatches))
       allocate(csite%dmean_rh                (npatches))
       allocate(csite%dmean_A_decomp          (npatches))
       allocate(csite%dmean_Af_decomp         (npatches))
    end if
    if (imoutput > 0 ) then
       allocate(csite%mmean_rk4step           (npatches))
       allocate(csite%mmean_lambda_light      (npatches))
       allocate(csite%mmean_co2_residual      (npatches))
       allocate(csite%mmean_energy_residual   (npatches))
       allocate(csite%mmean_water_residual    (npatches))
       allocate(csite%mmean_rh                (npatches))
       allocate(csite%mmean_A_decomp          (npatches))
       allocate(csite%mmean_Af_decomp         (npatches))
    end if

    ! Initialize the variables with a non-sense number.
    !call huge_sitetype(csite)

    return
  end subroutine allocate_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_patchtype(cpatch,ncohorts)

    implicit none

    integer :: ncohorts
    type(patchtype),target :: cpatch
    
    call nullify_patchtype(cpatch)
    cpatch%ncohorts = ncohorts

    !----- We need to leave the subroutine in case this patch is empty. -------------------!
    if (ncohorts == 0) return
    
    allocate(cpatch%pft(ncohorts))
    allocate(cpatch%nplant(ncohorts))
    allocate(cpatch%hite(ncohorts))
    allocate(cpatch%agb(ncohorts))
    allocate(cpatch%basarea(ncohorts))
    allocate(cpatch%dagb_dt(ncohorts))
    allocate(cpatch%dba_dt(ncohorts))
    allocate(cpatch%ddbh_dt(ncohorts))
    allocate(cpatch%dbh(ncohorts))
    allocate(cpatch%bdead(ncohorts))
    allocate(cpatch%bleaf(ncohorts))
    allocate(cpatch%phenology_status(ncohorts))
    allocate(cpatch%balive(ncohorts))
    allocate(cpatch%broot(ncohorts))
    allocate(cpatch%bsapwood(ncohorts))
    allocate(cpatch%lai(ncohorts))
    allocate(cpatch%wpa(ncohorts))
    allocate(cpatch%wai(ncohorts))
    allocate(cpatch%solvable(ncohorts))
    allocate(cpatch%bstorage(ncohorts))
    allocate(cpatch%cb(13,ncohorts))
    allocate(cpatch%cb_max(13,ncohorts))
    allocate(cpatch%cbr_bar(ncohorts))
    allocate(cpatch%veg_energy(ncohorts))
    allocate(cpatch%veg_temp(ncohorts))
    allocate(cpatch%veg_fliq(ncohorts))
    allocate(cpatch%veg_water(ncohorts))
    allocate(cpatch%mean_gpp(ncohorts))
    allocate(cpatch%mean_leaf_resp(ncohorts))
    allocate(cpatch%mean_root_resp(ncohorts))
    allocate(cpatch%mean_growth_resp(ncohorts))
    allocate(cpatch%mean_storage_resp(ncohorts))
    allocate(cpatch%mean_vleaf_resp(ncohorts))
    allocate(cpatch%today_leaf_resp(ncohorts))
    allocate(cpatch%today_root_resp(ncohorts))
    allocate(cpatch%today_gpp(ncohorts))
    allocate(cpatch%today_gpp_pot(ncohorts))
    allocate(cpatch%today_gpp_max(ncohorts))
    allocate(cpatch%growth_respiration(ncohorts))
    allocate(cpatch%storage_respiration(ncohorts))
    allocate(cpatch%vleaf_respiration(ncohorts))
    allocate(cpatch%fsn(ncohorts))
    allocate(cpatch%monthly_dndt(ncohorts))
    allocate(cpatch%mort_rate(n_mort,ncohorts))

    ! Soon to be replaced
    allocate(cpatch%old_stoma_data(ncohorts))
    allocate(cpatch%old_stoma_vector(n_stoma_atts,ncohorts))

    allocate(cpatch%Psi_open(ncohorts))
    allocate(cpatch%krdepth(ncohorts))
    allocate(cpatch%first_census(ncohorts))
    allocate(cpatch%new_recruit_flag(ncohorts))
    allocate(cpatch%light_level(ncohorts))
    allocate(cpatch%light_level_beam(ncohorts))
    allocate(cpatch%light_level_diff(ncohorts))
    allocate(cpatch%beamext_level(ncohorts))
    allocate(cpatch%diffext_level(ncohorts))
    allocate(cpatch%norm_par_beam(ncohorts))
    allocate(cpatch%norm_par_diff(ncohorts))
    allocate(cpatch%lambda_light(ncohorts))
    allocate(cpatch%par_v(ncohorts))
    allocate(cpatch%par_v_beam(ncohorts))
    allocate(cpatch%par_v_diffuse(ncohorts))
    allocate(cpatch%rshort_v(ncohorts))
    allocate(cpatch%rshort_v_beam(ncohorts))
    allocate(cpatch%rshort_v_diffuse(ncohorts))
    allocate(cpatch%rlong_v(ncohorts))
    allocate(cpatch%rlong_v_surf(ncohorts))
    allocate(cpatch%rlong_v_incid(ncohorts))
    allocate(cpatch%rb(ncohorts))
    allocate(cpatch%A_open(ncohorts))
    allocate(cpatch%A_closed(ncohorts))
    allocate(cpatch%Psi_closed(ncohorts))
    allocate(cpatch%rsw_open(ncohorts))
    allocate(cpatch%rsw_closed(ncohorts))
    allocate(cpatch%fsw(ncohorts))
    allocate(cpatch%fs_open(ncohorts))
    allocate(cpatch%stomatal_resistance(ncohorts))
    allocate(cpatch%leaf_maintenance(ncohorts))
    allocate(cpatch%root_maintenance(ncohorts))
    allocate(cpatch%leaf_drop(ncohorts))
    allocate(cpatch%bseeds(ncohorts))
    allocate(cpatch%leaf_respiration(ncohorts))
    allocate(cpatch%root_respiration(ncohorts))
    allocate(cpatch%hcapveg(ncohorts))
    allocate(cpatch%gpp(ncohorts))
    allocate(cpatch%paw_avg(ncohorts))
    allocate(cpatch%turnover_amp(ncohorts))
    allocate(cpatch%llspan(ncohorts))
    allocate(cpatch%vm_bar(ncohorts))
    allocate(cpatch%sla(ncohorts))

    if (idoutput > 0 .or. imoutput > 0) then
       allocate(cpatch%dmean_par_v(ncohorts))
       allocate(cpatch%dmean_par_v_beam(ncohorts))
       allocate(cpatch%dmean_par_v_diff(ncohorts))
       allocate(cpatch%dmean_light_level(ncohorts))
       allocate(cpatch%dmean_light_level_beam(ncohorts))
       allocate(cpatch%dmean_light_level_diff(ncohorts))
       allocate(cpatch%dmean_beamext_level(ncohorts))
       allocate(cpatch%dmean_diffext_level(ncohorts))
       allocate(cpatch%dmean_norm_par_beam(ncohorts))
       allocate(cpatch%dmean_norm_par_diff(ncohorts))
       allocate(cpatch%dmean_lambda_light(ncohorts))
       allocate(cpatch%dmean_fs_open(ncohorts))
       allocate(cpatch%dmean_fsw(ncohorts))
       allocate(cpatch%dmean_fsn(ncohorts))
       allocate(cpatch%dmean_gpp(ncohorts))
       allocate(cpatch%dmean_leaf_resp(ncohorts))
       allocate(cpatch%dmean_root_resp(ncohorts))
    end if
    
    if (imoutput > 0) then
       allocate(cpatch%mmean_par_v(ncohorts))
       allocate(cpatch%mmean_par_v_beam(ncohorts))
       allocate(cpatch%mmean_par_v_diff(ncohorts))
       allocate(cpatch%mmean_light_level(ncohorts))
       allocate(cpatch%mmean_light_level_beam(ncohorts))
       allocate(cpatch%mmean_light_level_diff(ncohorts))
       allocate(cpatch%mmean_beamext_level(ncohorts))
       allocate(cpatch%mmean_diffext_level(ncohorts))
       allocate(cpatch%mmean_norm_par_beam(ncohorts))
       allocate(cpatch%mmean_norm_par_diff(ncohorts))
       allocate(cpatch%mmean_lambda_light(ncohorts))
       allocate(cpatch%mmean_fs_open(ncohorts))
       allocate(cpatch%mmean_fsw(ncohorts))
       allocate(cpatch%mmean_fsn(ncohorts))
       allocate(cpatch%mmean_leaf_maintenance(ncohorts))
       allocate(cpatch%mmean_root_maintenance(ncohorts))
       allocate(cpatch%mmean_leaf_drop(ncohorts))
       allocate(cpatch%mmean_cb(ncohorts))
       allocate(cpatch%mmean_gpp(ncohorts))
       allocate(cpatch%mmean_leaf_resp(ncohorts))
       allocate(cpatch%mmean_root_resp(ncohorts))
       allocate(cpatch%mmean_growth_resp(ncohorts))
       allocate(cpatch%mmean_storage_resp(ncohorts))
       allocate(cpatch%mmean_vleaf_resp(ncohorts))
       allocate(cpatch%mmean_mort_rate(n_mort,ncohorts))
    end if


    ! Initialize the variables with a non-sense number.
    !call huge_patchtype(cpatch)

    return
  end subroutine allocate_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_edtype(cgrid)
    
    implicit none
    type(edtype),target :: cgrid


       nullify(cgrid%polygon                 )
       nullify(cgrid%lat                     )
       nullify(cgrid%lon                     )
       nullify(cgrid%natm                    )
       nullify(cgrid%xatm                    )
       nullify(cgrid%yatm                    )
       nullify(cgrid%ntext_soil              )
       nullify(cgrid%lsl                     )
       
       nullify(cgrid%pysi_id                 )
       nullify(cgrid%pysi_n                  )
       nullify(cgrid%walltime_py             )
       nullify(cgrid%sensflux_py             )
       nullify(cgrid%site_adjacency          )
       nullify(cgrid%wbar                    )
       nullify(cgrid%Te                      )
       nullify(cgrid%zbar                    )
!! NOT USED       nullify(cgrid%tau                     )
       nullify(cgrid%sheat                   )
       nullify(cgrid%baseflow                )
       nullify(cgrid%runoff                  )
       nullify(cgrid%swliq                   )
       
       nullify(cgrid%ilon                    )
       nullify(cgrid%ilat                    )
       nullify(cgrid%total_agb               )
       nullify(cgrid%total_basal_area        )
       nullify(cgrid%total_agb_growth        )
       nullify(cgrid%total_agb_mort          )
       nullify(cgrid%total_agb_recruit       )
       nullify(cgrid%total_basal_area_growth )
       nullify(cgrid%total_basal_area_mort   )
       nullify(cgrid%total_basal_area_recruit)
       nullify(cgrid%nsites                  )
       nullify(cgrid%sitenums                )
       nullify(cgrid%load_adjacency          )
       nullify(cgrid%cosz                    )
       nullify(cgrid%mean_gpp                )
       nullify(cgrid%mean_precip             )
       nullify(cgrid%mean_qprecip            )
       nullify(cgrid%mean_netrad             )
       nullify(cgrid%cbudget_initialstorage  )
       nullify(cgrid%cbudget_nep             )
       nullify(cgrid%nbudget_initialstorage  )
       nullify(cgrid%basal_area              )
       nullify(cgrid%agb                     )
       nullify(cgrid%pldens                  )
       nullify(cgrid%bseeds                  )
       
       nullify(cgrid%metinput                )
       nullify(cgrid%met                     )

       nullify(cgrid%lapse                   )

       nullify(cgrid%lai                     )
       nullify(cgrid%avg_lma                 )
       nullify(cgrid%wpa                     )
       nullify(cgrid%wai                     )

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       nullify(cgrid%avg_vapor_vc            )
       nullify(cgrid%avg_dew_cg              )
       nullify(cgrid%avg_vapor_gc            )
       nullify(cgrid%avg_wshed_vg            )
       nullify(cgrid%avg_intercepted         )
       nullify(cgrid%avg_vapor_ac            )
       nullify(cgrid%avg_transp              )
       nullify(cgrid%avg_evap                )
       nullify(cgrid%avg_smoist_gg           )
       nullify(cgrid%avg_smoist_gc           )
       nullify(cgrid%avg_runoff              )
       nullify(cgrid%avg_drainage            )
       nullify(cgrid%avg_drainage_heat       )
       nullify(cgrid%aux                     )
       nullify(cgrid%aux_s                   )
       nullify(cgrid%avg_sensible_vc         )
       nullify(cgrid%avg_qwshed_vg           )
       nullify(cgrid%avg_qintercepted        )
       nullify(cgrid%avg_sensible_gc         )
       nullify(cgrid%avg_sensible_ac         )
       nullify(cgrid%avg_sensible_gg         )
       nullify(cgrid%avg_runoff_heat         )

       ! Fast time state diagnostics 
       nullify(cgrid%avg_veg_energy          )
       nullify(cgrid%avg_veg_hcap            )
       nullify(cgrid%avg_veg_temp            )
       nullify(cgrid%avg_veg_fliq            )
       nullify(cgrid%avg_veg_water           )
       nullify(cgrid%avg_can_temp            )
       nullify(cgrid%avg_can_shv             )
       nullify(cgrid%avg_can_co2             )
       nullify(cgrid%avg_can_rhos            )
       nullify(cgrid%avg_can_prss            )
       nullify(cgrid%avg_can_theta           )
       nullify(cgrid%avg_can_enthalpy        )
       nullify(cgrid%avg_can_depth           )
       nullify(cgrid%avg_soil_energy         )
       nullify(cgrid%avg_soil_water          )
       nullify(cgrid%avg_soil_temp           )
       nullify(cgrid%avg_soil_fracliq        )
       nullify(cgrid%avg_soil_wetness        )
       nullify(cgrid%avg_skin_temp           )
       nullify(cgrid%avg_available_water     )

       nullify(cgrid%avg_lai_ebalvars)

       nullify(cgrid%avg_gpp             )
       nullify(cgrid%avg_leaf_resp       )
       nullify(cgrid%avg_root_resp       )
       nullify(cgrid%avg_growth_resp     )
       nullify(cgrid%avg_storage_resp    )
       nullify(cgrid%avg_vleaf_resp      )
       nullify(cgrid%avg_plant_resp      )
       nullify(cgrid%avg_htroph_resp     )
       nullify(cgrid%avg_leaf_drop       )
       nullify(cgrid%avg_leaf_maintenance)
       nullify(cgrid%avg_root_maintenance)

       !!! added for NACP intercomparison (MCD)
       nullify(cgrid%avg_sfcw_depth     )
       nullify(cgrid%avg_sfcw_energy    )
       nullify(cgrid%avg_sfcw_mass      )
       nullify(cgrid%avg_sfcw_tempk     )
       nullify(cgrid%avg_sfcw_fracliq   )
       nullify(cgrid%avg_bdead         )
       nullify(cgrid%avg_balive        )
       nullify(cgrid%avg_bleaf          )
       nullify(cgrid%avg_broot          )
       nullify(cgrid%avg_bsapwood       )
       nullify(cgrid%avg_bstorage       )
       nullify(cgrid%avg_bseeds         )

       nullify(cgrid%avg_fsc           )
       nullify(cgrid%avg_ssc           )
       nullify(cgrid%avg_stsc          )
       nullify(cgrid%avg_fsn           )
       nullify(cgrid%avg_msn           )


       
     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     nullify(cgrid%Cleaf  )
     nullify(cgrid%Croot  )
     nullify(cgrid%Cstore )
     nullify(cgrid%Ccwd   )
     nullify(cgrid%Nleaf  )
     nullify(cgrid%Ndead  )
     nullify(cgrid%Nroot  )
     nullify(cgrid%Nstore )
     nullify(cgrid%Ncwd   )
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     nullify(cgrid%Cleaf_grow       )
     nullify(cgrid%Croot_grow       )
     nullify(cgrid%Cdead_grow       )
     nullify(cgrid%Cstore_grow      )
     nullify(cgrid%Cleaf_litter_flux)
     nullify(cgrid%Croot_litter_flux)
     nullify(cgrid%Ccwd_flux        )
     nullify(cgrid%Nleaf_grow       )
     nullify(cgrid%Ndead_grow       )
     nullify(cgrid%Nroot_grow       )
     nullify(cgrid%Nstore_grow      )
     nullify(cgrid%Nleaf_litter_flux)
     nullify(cgrid%Nroot_litter_flux)
     nullify(cgrid%Ncwd_flux        )
     nullify(cgrid%Nbiomass_uptake  )
     nullify(cgrid%Ngross_min       )
     nullify(cgrid%Nnet_min         )


       ! Meteorologic conditions (forcing)
       nullify(cgrid%avg_nir_beam            )
       nullify(cgrid%avg_nir_diffuse         )
       nullify(cgrid%avg_par_beam            )
       nullify(cgrid%avg_par_diffuse         )
       nullify(cgrid%avg_atm_tmp             )
       nullify(cgrid%avg_atm_shv             )
       nullify(cgrid%avg_rshort              )
       nullify(cgrid%avg_rshort_diffuse      )
       nullify(cgrid%avg_rlong               )
       nullify(cgrid%avg_pcpg                )
       nullify(cgrid%avg_qpcpg               )
       nullify(cgrid%avg_dpcpg               )
       nullify(cgrid%avg_vels                )
       nullify(cgrid%avg_atm_prss            )
       nullify(cgrid%avg_exner               )
       nullify(cgrid%avg_geoht               )
       nullify(cgrid%avg_atm_co2             )
       nullify(cgrid%avg_albedt              )
       nullify(cgrid%avg_rlongup             )
       nullify(cgrid%dmean_pcpg              )
       nullify(cgrid%dmean_runoff            )
       nullify(cgrid%dmean_drainage          )
       nullify(cgrid%dmean_vapor_ac          )
       nullify(cgrid%dmean_vapor_gc          )
       nullify(cgrid%dmean_vapor_vc          )
       nullify(cgrid%dmean_gpp               )
       nullify(cgrid%dmean_evap              )
       nullify(cgrid%dmean_transp            )
       nullify(cgrid%dmean_sensible_vc       )
       nullify(cgrid%dmean_sensible_gc       )
       nullify(cgrid%dmean_sensible_ac       )
       nullify(cgrid%dmean_plresp            )
       nullify(cgrid%dmean_rh                )
       nullify(cgrid%dmean_leaf_resp         )
       nullify(cgrid%dmean_root_resp         )
       nullify(cgrid%dmean_growth_resp       )
       nullify(cgrid%dmean_storage_resp      )
       nullify(cgrid%dmean_vleaf_resp        )
       nullify(cgrid%dmean_nep               )
       nullify(cgrid%dmean_soil_temp         )
       nullify(cgrid%dmean_soil_water        )
       nullify(cgrid%dmean_fs_open           )
       nullify(cgrid%dmean_fsw               )
       nullify(cgrid%dmean_fsn               )
       nullify(cgrid%dmean_gpp_dbh           )
       nullify(cgrid%dmean_can_temp          )
       nullify(cgrid%dmean_can_shv           )
       nullify(cgrid%dmean_can_co2           )
       nullify(cgrid%dmean_can_rhos          )
       nullify(cgrid%dmean_can_prss          )
       nullify(cgrid%dmean_can_theta         )
       nullify(cgrid%dmean_veg_energy        )
       nullify(cgrid%dmean_veg_water         )
       nullify(cgrid%dmean_veg_hcap          )
       nullify(cgrid%dmean_veg_temp          )
       nullify(cgrid%dmean_atm_temp          )
       nullify(cgrid%dmean_atm_shv           )
       nullify(cgrid%dmean_atm_prss          )
       nullify(cgrid%dmean_atm_vels          )
       nullify(cgrid%dmean_rshort            )
       nullify(cgrid%dmean_rlong             )
       nullify(cgrid%dmean_co2_residual      )
       nullify(cgrid%dmean_energy_residual   )
       nullify(cgrid%dmean_water_residual    )
       nullify(cgrid%lai_pft                 )
       nullify(cgrid%wpa_pft                 )
       nullify(cgrid%wai_pft                 )
       nullify(cgrid%mmean_gpp               )
       nullify(cgrid%mmean_evap              )
       nullify(cgrid%mmean_transp            )
       nullify(cgrid%mmean_sensible_vc       )
       nullify(cgrid%mmean_sensible_gc       )
       nullify(cgrid%mmean_sensible_ac       )
       nullify(cgrid%mmean_vapor_vc          )
       nullify(cgrid%mmean_vapor_gc          )
       nullify(cgrid%mmean_vapor_ac          )
       nullify(cgrid%mmean_nep               )
       nullify(cgrid%mmean_soil_temp         )
       nullify(cgrid%mmean_soil_water        )
       nullify(cgrid%mmean_fs_open           )
       nullify(cgrid%mmean_fsw               )
       nullify(cgrid%mmean_fsn               )
       nullify(cgrid%mmean_plresp            )
       nullify(cgrid%mmean_rh                )
       nullify(cgrid%mmean_leaf_resp         )
       nullify(cgrid%mmean_root_resp         )
       nullify(cgrid%mmean_growth_resp       )
       nullify(cgrid%mmean_storage_resp      )
       nullify(cgrid%mmean_vleaf_resp        )
       nullify(cgrid%mmean_gpp_dbh           )
       nullify(cgrid%mmean_lai_pft           )
       nullify(cgrid%mmean_wpa_pft           )
       nullify(cgrid%mmean_wai_pft           )
       nullify(cgrid%mmean_can_temp          )
       nullify(cgrid%mmean_can_shv           )
       nullify(cgrid%mmean_can_co2           )
       nullify(cgrid%mmean_can_rhos          )
       nullify(cgrid%mmean_can_prss          )
       nullify(cgrid%mmean_can_theta         )
       nullify(cgrid%mmean_veg_energy        )
       nullify(cgrid%mmean_veg_water         )
       nullify(cgrid%mmean_veg_hcap          )
       nullify(cgrid%mmean_veg_temp          )
       nullify(cgrid%mmean_atm_temp          )
       nullify(cgrid%mmean_rshort            )
       nullify(cgrid%mmean_rlong             )
       nullify(cgrid%mmean_atm_shv           )
       nullify(cgrid%mmean_atm_prss          )
       nullify(cgrid%mmean_atm_vels          )
       nullify(cgrid%mmean_pcpg              )
       nullify(cgrid%mmean_runoff            )
       nullify(cgrid%mmean_drainage          )
       nullify(cgrid%mmean_co2_residual      )
       nullify(cgrid%mmean_energy_residual   )
       nullify(cgrid%mmean_water_residual    )
       nullify(cgrid%bseeds_pft              )
       nullify(cgrid%agb_pft                 )
       nullify(cgrid%ba_pft                  )
       nullify(cgrid%stdev_gpp               )
       nullify(cgrid%stdev_evap              )
       nullify(cgrid%stdev_transp            )
       nullify(cgrid%stdev_sensible          )
       nullify(cgrid%stdev_nep               )
       nullify(cgrid%stdev_rh                )
       nullify(cgrid%disturbance_rates       )
       nullify(cgrid%workload                )

    return
  end subroutine nullify_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_polygontype(cpoly)

    implicit none

    type(polygontype),target :: cpoly

    nullify(cpoly%sipa_id)
    nullify(cpoly%sipa_n)
    nullify(cpoly%patch_count)  
    nullify(cpoly%site)
    nullify(cpoly%sitenum)

    nullify(cpoly%area)
    nullify(cpoly%patch_area)
    nullify(cpoly%elevation)
    nullify(cpoly%slope)
    nullify(cpoly%aspect)

    nullify(cpoly%num_landuse_years)
    nullify(cpoly%lai_pft)
    nullify(cpoly%wpa_pft)
    nullify(cpoly%wai_pft)

    nullify(cpoly%TCI)   
    nullify(cpoly%pptweight)   
    nullify(cpoly%lsl)   
    nullify(cpoly%hydro_next)
    nullify(cpoly%hydro_prev)
    nullify(cpoly%moist_W)
    nullify(cpoly%moist_f)  
    nullify(cpoly%moist_tau)
    nullify(cpoly%moist_zi) 
    nullify(cpoly%baseflow) 
    nullify(cpoly%ntext_soil)
    nullify(cpoly%min_monthly_temp)
    nullify(cpoly%plantation) 
    nullify(cpoly%agri_stocking_pft)
    nullify(cpoly%agri_stocking_density)
    nullify(cpoly%plantation_stocking_pft)
    nullify(cpoly%plantation_stocking_density)
    nullify(cpoly%primary_harvest_memory)
    nullify(cpoly%secondary_harvest_memory)
    nullify(cpoly%fire_disturbance_rate)
    nullify(cpoly%ignition_rate)
    nullify(cpoly%lambda_fire)
    nullify(cpoly%phen_pars)
    nullify(cpoly%treefall_disturbance_rate)
    nullify(cpoly%nat_disturbance_rate)
    nullify(cpoly%nat_dist_type)
    nullify(cpoly%disturbance_memory)
    nullify(cpoly%disturbance_rates)
    nullify(cpoly%loss_fraction)
    nullify(cpoly%green_leaf_factor)
    nullify(cpoly%leaf_aging_factor)
    nullify(cpoly%met)
    nullify(cpoly%basal_area)
    nullify(cpoly%agb)    
    nullify(cpoly%pldens)    
    nullify(cpoly%bseeds)    
    
    nullify(cpoly%basal_area_growth)
    nullify(cpoly%agb_growth       )
    nullify(cpoly%basal_area_mort  )
    nullify(cpoly%agb_mort         )
    nullify(cpoly%basal_area_cut   ) !NOT IN REGISTRY
    nullify(cpoly%agb_cut          )

    nullify(cpoly%cosaoi)
    nullify(cpoly%albedo_beam)
    nullify(cpoly%albedo_diffuse)
    nullify(cpoly%rlong_albedo)
    
    nullify(cpoly%albedt)
    nullify(cpoly%rlongup)
    nullify(cpoly%avg_lma    )
    nullify(cpoly%daylight)
    nullify(cpoly%lai    )
    nullify(cpoly%wpa    )
    nullify(cpoly%wai    )
    
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    nullify(cpoly%avg_vapor_vc  )
    nullify(cpoly%avg_dew_cg    )
    nullify(cpoly%avg_vapor_gc  )
    nullify(cpoly%avg_wshed_vg  )
    nullify(cpoly%avg_intercepted)
    nullify(cpoly%avg_vapor_ac  )
    nullify(cpoly%avg_transp    )
    nullify(cpoly%avg_evap      )
    nullify(cpoly%avg_smoist_gg )
    nullify(cpoly%avg_smoist_gc )
    nullify(cpoly%avg_runoff    )
    nullify(cpoly%avg_drainage  )
    nullify(cpoly%avg_drainage_heat  )
    nullify(cpoly%aux           )
    nullify(cpoly%aux_s         )
    nullify(cpoly%avg_sensible_vc  )
    nullify(cpoly%avg_qwshed_vg    )
    nullify(cpoly%avg_qintercepted )
    nullify(cpoly%avg_sensible_gc  )
    nullify(cpoly%avg_sensible_ac  )
    nullify(cpoly%avg_sensible_gg  )
    nullify(cpoly%avg_runoff_heat  )

    ! ----------------------------------------------
    nullify(cpoly%avg_veg_energy)
    nullify(cpoly%avg_veg_hcap  )
    nullify(cpoly%avg_veg_temp  )
    nullify(cpoly%avg_veg_fliq  )
    nullify(cpoly%avg_veg_water )
    nullify(cpoly%avg_can_temp  )
    nullify(cpoly%avg_can_shv   )
    nullify(cpoly%avg_can_co2   )
    nullify(cpoly%avg_can_rhos  )
    nullify(cpoly%avg_can_prss  )
    nullify(cpoly%avg_can_theta )
    nullify(cpoly%avg_can_enthalpy)
    nullify(cpoly%avg_can_depth)
    nullify(cpoly%avg_soil_energy)
    nullify(cpoly%avg_soil_water)
    nullify(cpoly%avg_soil_temp )
    nullify(cpoly%avg_soil_fracliq )
    nullify(cpoly%avg_soil_wetness )
    nullify(cpoly%avg_skin_temp  )
    nullify(cpoly%avg_available_water)
    nullify(cpoly%runoff        )
    ! Phenology
    nullify(cpoly%rad_avg          )
    ! Meteorological conditions
    nullify(cpoly%avg_atm_tmp      )
    nullify(cpoly%avg_atm_shv      )
    nullify(cpoly%avg_atm_prss     )
    ! NACP
    nullify(cpoly%avg_sfcw_depth   )
    nullify(cpoly%avg_sfcw_energy  )
    nullify(cpoly%avg_sfcw_mass    )
    nullify(cpoly%avg_sfcw_tempk   )
    nullify(cpoly%avg_sfcw_fracliq )
    nullify(cpoly%avg_fsc          )
    nullify(cpoly%avg_stsc         )
    nullify(cpoly%avg_ssc          )
    nullify(cpoly%avg_balive       )
    nullify(cpoly%avg_bleaf        )
    nullify(cpoly%avg_broot        )
    nullify(cpoly%avg_bsapwood     )
    nullify(cpoly%avg_bdead        )
    nullify(cpoly%avg_fsn          )
    nullify(cpoly%avg_msn          )
    nullify(cpoly%avg_bstorage     )
    nullify(cpoly%avg_bseeds       )

    nullify(cpoly%dmean_co2_residual      )
    nullify(cpoly%dmean_energy_residual   )
    nullify(cpoly%dmean_water_residual    )
    nullify(cpoly%mmean_co2_residual      )
    nullify(cpoly%mmean_energy_residual   )
    nullify(cpoly%mmean_water_residual    )

    return
  end subroutine nullify_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_sitetype(csite)

    implicit none

    type(sitetype),target :: csite
    
    nullify(csite%paco_id)
    nullify(csite%paco_n)

    nullify(csite%dist_type)
    nullify(csite%age)
    nullify(csite%area)
    nullify(csite%laiarea)
    nullify(csite%hcapveg)
    nullify(csite%fast_soil_C)
    nullify(csite%slow_soil_C)
    nullify(csite%structural_soil_C)
    nullify(csite%structural_soil_L)
    nullify(csite%mineralized_soil_N)
    nullify(csite%fast_soil_N)
    nullify(csite%pname)
    nullify(csite%sum_dgd)
    nullify(csite%sum_chd)
    nullify(csite%plantation) 
    nullify(csite%cohort_count)
    nullify(csite%can_enthalpy)
    nullify(csite%can_temp)
    nullify(csite%can_shv)
    nullify(csite%can_co2)
    nullify(csite%can_rhos)
    nullify(csite%can_prss)
    nullify(csite%can_theta)
    nullify(csite%can_depth)
    nullify(csite%lambda_light)
    nullify(csite%dmean_lambda_light)
    nullify(csite%mmean_lambda_light)
    nullify(csite%lai)
    nullify(csite%wpa)
    nullify(csite%wai)

    nullify(csite%sfcwater_mass)
    nullify(csite%sfcwater_energy)
    nullify(csite%sfcwater_depth)
    nullify(csite%rshort_s)
    nullify(csite%rshort_s_beam)
    nullify(csite%rshort_s_diffuse)
    nullify(csite%sfcwater_tempk)
    nullify(csite%sfcwater_fracliq)
    nullify(csite%nlev_sfcwater)
    nullify(csite%ntext_soil)
    nullify(csite%soil_energy)
    nullify(csite%soil_water)
    nullify(csite%soil_tempk)
    nullify(csite%soil_fracliq)
    nullify(csite%ground_shv)
    nullify(csite%surface_ssh)
    nullify(csite%rough)
    nullify(csite%A_o_max) 
    nullify(csite%A_c_max) 
    nullify(csite%old_stoma_data_max)
    nullify(csite%old_stoma_vector_max)
    nullify(csite%avg_daily_temp)  
    nullify(csite%mean_rh)
    nullify(csite%dmean_rh)
    nullify(csite%mmean_rh)
    nullify(csite%mean_nep)
    nullify(csite%wbudget_loss2atm)
    nullify(csite%wbudget_denseffect)
    nullify(csite%wbudget_precipgain)
    nullify(csite%wbudget_loss2runoff)
    nullify(csite%wbudget_loss2drainage)
    nullify(csite%wbudget_initialstorage)
    nullify(csite%wbudget_residual)
    nullify(csite%ebudget_latent)
    nullify(csite%ebudget_loss2atm)
    nullify(csite%ebudget_denseffect)
    nullify(csite%ebudget_loss2runoff)
    nullify(csite%ebudget_loss2drainage)
    nullify(csite%ebudget_netrad)
    nullify(csite%ebudget_precipgain)
    nullify(csite%ebudget_initialstorage)
    nullify(csite%ebudget_residual)
    nullify(csite%co2budget_initialstorage)
    nullify(csite%co2budget_residual)
    nullify(csite%co2budget_loss2atm)
    nullify(csite%co2budget_denseffect)
    nullify(csite%co2budget_gpp)
    nullify(csite%co2budget_gpp_dbh)
    nullify(csite%co2budget_plresp)
    nullify(csite%co2budget_rh)
    nullify(csite%today_A_decomp)
    nullify(csite%today_Af_decomp)
    nullify(csite%dmean_A_decomp)
    nullify(csite%dmean_Af_decomp)
    nullify(csite%mmean_A_decomp)
    nullify(csite%mmean_Af_decomp)
    nullify(csite%repro)
    nullify(csite%veg_rough)
    nullify(csite%veg_height)
    nullify(csite%fsc_in)
    nullify(csite%ssc_in)
    nullify(csite%ssl_in)
    nullify(csite%fsn_in)
    nullify(csite%total_plant_nitrogen_uptake)
    nullify(csite%mineralized_N_loss)
    nullify(csite%mineralized_N_input)
    nullify(csite%rshort_g)
    nullify(csite%rshort_g_beam)
    nullify(csite%rshort_g_diffuse)
    nullify(csite%rlong_g)
    nullify(csite%rlong_g_surf)
    nullify(csite%rlong_g_incid)
    nullify(csite%rlong_s)
    nullify(csite%rlong_s_surf)
    nullify(csite%rlong_s_incid)
    nullify(csite%albedt)
    nullify(csite%albedo_beam)
    nullify(csite%albedo_diffuse)
    nullify(csite%rlongup)
    nullify(csite%rlong_albedo)
    nullify(csite%total_snow_depth)
    nullify(csite%snowfac)
    nullify(csite%A_decomp)
    nullify(csite%f_decomp)
    nullify(csite%rh)
    nullify(csite%cwd_rh)
    nullify(csite%fuse_flag)
    nullify(csite%pft_density_profile)
    nullify(csite%plant_ag_biomass)

    nullify(csite%mean_wflux)
    nullify(csite%mean_latflux)
    nullify(csite%mean_hflux)
    nullify(csite%mean_runoff)
    nullify(csite%mean_qrunoff)

    nullify(csite%htry)
    nullify(csite%avg_rk4step)
    nullify(csite%dmean_rk4step)
    nullify(csite%mmean_rk4step)

    nullify(csite%avg_available_water)

    nullify(csite%ustar)
    nullify(csite%tstar)
    nullify(csite%qstar)
    nullify(csite%cstar)

    nullify(csite%upwp)
    nullify(csite%qpwp)
    nullify(csite%cpwp)
    nullify(csite%tpwp)
    nullify(csite%wpwp)


    nullify(csite%avg_carbon_ac)

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    nullify(csite%avg_vapor_vc  )
    nullify(csite%avg_dew_cg    )
    nullify(csite%avg_vapor_gc  )
    nullify(csite%avg_wshed_vg  )
    nullify(csite%avg_intercepted)
    nullify(csite%avg_vapor_ac  )
    nullify(csite%avg_transp    )
    nullify(csite%avg_evap      )
    nullify(csite%avg_netrad    )
    nullify(csite%avg_smoist_gg )
    nullify(csite%avg_smoist_gc )
    nullify(csite%avg_runoff    )
    nullify(csite%avg_drainage  )
    nullify(csite%avg_drainage_heat  )
    nullify(csite%aux           )
    nullify(csite%aux_s         )
    nullify(csite%avg_sensible_vc  )
    nullify(csite%avg_qwshed_vg    )
    nullify(csite%avg_qintercepted )
    nullify(csite%avg_sensible_gc  )
    nullify(csite%avg_sensible_ac  )
    nullify(csite%avg_sensible_gg  )
    nullify(csite%avg_runoff_heat  )

    ! ----------------------------------------------
    nullify(csite%avg_veg_energy) 
    nullify(csite%avg_veg_temp) 
    nullify(csite%avg_veg_fliq) 
    nullify(csite%avg_veg_water)

    nullify(csite%watertable      )
    nullify(csite%moist_dz        )
    nullify(csite%ksat            )
    nullify(csite%soil_sat_energy )
    nullify(csite%soil_sat_water  )
    nullify(csite%soil_sat_heat   )
    nullify(csite%runoff_A        )
    nullify(csite%runoff_rate     )
    nullify(csite%runoff          )

    nullify(csite%patch)

    nullify(csite%dmean_co2_residual      )
    nullify(csite%dmean_energy_residual   )
    nullify(csite%dmean_water_residual    )
    nullify(csite%mmean_co2_residual      )
    nullify(csite%mmean_energy_residual   )
    nullify(csite%mmean_water_residual    )

    return
  end subroutine nullify_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_patchtype(cpatch)

    implicit none

    type(patchtype),target :: cpatch
    
    nullify(cpatch%pft)
    nullify(cpatch%nplant)
    nullify(cpatch%hite)
    nullify(cpatch%agb)
    nullify(cpatch%basarea)
    nullify(cpatch%dagb_dt)
    nullify(cpatch%dba_dt)
    nullify(cpatch%ddbh_dt)
    nullify(cpatch%dbh)
    nullify(cpatch%bdead)
    nullify(cpatch%bleaf)
    nullify(cpatch%phenology_status)
    nullify(cpatch%balive)
    nullify(cpatch%broot)
    nullify(cpatch%bsapwood)
    nullify(cpatch%lai)
    nullify(cpatch%wpa)
    nullify(cpatch%wai)
    nullify(cpatch%solvable)
    nullify(cpatch%bstorage)
    nullify(cpatch%cb)
    nullify(cpatch%cb_max)
    nullify(cpatch%cbr_bar)
    nullify(cpatch%mmean_cb)
    nullify(cpatch%veg_energy)
    nullify(cpatch%veg_temp)
    nullify(cpatch%veg_fliq)
    nullify(cpatch%veg_water)
    nullify(cpatch%mean_gpp)
    nullify(cpatch%mean_leaf_resp)
    nullify(cpatch%mean_root_resp)
    nullify(cpatch%mean_storage_resp)
    nullify(cpatch%mean_growth_resp)
    nullify(cpatch%mean_vleaf_resp)
    nullify(cpatch%today_leaf_resp)
    nullify(cpatch%today_root_resp)
    nullify(cpatch%today_gpp)
    nullify(cpatch%today_gpp_pot)
    nullify(cpatch%today_gpp_max)
    nullify(cpatch%dmean_gpp         )
    nullify(cpatch%dmean_leaf_resp   )
    nullify(cpatch%dmean_root_resp   )
    nullify(cpatch%mmean_gpp         )
    nullify(cpatch%mmean_leaf_resp   )
    nullify(cpatch%mmean_root_resp   )
    nullify(cpatch%mmean_growth_resp )
    nullify(cpatch%mmean_storage_resp)
    nullify(cpatch%mmean_vleaf_resp  )
    nullify(cpatch%growth_respiration)
    nullify(cpatch%storage_respiration)
    nullify(cpatch%vleaf_respiration)
    nullify(cpatch%fsn)
    nullify(cpatch%monthly_dndt)
    nullify(cpatch%mort_rate)
    nullify(cpatch%mmean_mort_rate)
    nullify(cpatch%old_stoma_data)
    nullify(cpatch%old_stoma_vector)
    nullify(cpatch%Psi_open)
    nullify(cpatch%krdepth)
    nullify(cpatch%first_census)
    nullify(cpatch%new_recruit_flag)
    nullify(cpatch%light_level)
    nullify(cpatch%dmean_light_level)
    nullify(cpatch%mmean_light_level)
    nullify(cpatch%light_level_beam)
    nullify(cpatch%dmean_light_level_beam)
    nullify(cpatch%mmean_light_level_beam)
    nullify(cpatch%light_level_diff)
    nullify(cpatch%dmean_light_level_diff)
    nullify(cpatch%mmean_light_level_diff)
    nullify(cpatch%beamext_level)
    nullify(cpatch%dmean_beamext_level)
    nullify(cpatch%mmean_beamext_level)
    nullify(cpatch%diffext_level)
    nullify(cpatch%dmean_diffext_level)
    nullify(cpatch%mmean_diffext_level)
    nullify(cpatch%lambda_light)
    nullify(cpatch%dmean_lambda_light)
    nullify(cpatch%mmean_lambda_light)
    nullify(cpatch%norm_par_beam)
    nullify(cpatch%dmean_norm_par_beam)
    nullify(cpatch%mmean_norm_par_beam)
    nullify(cpatch%norm_par_diff)
    nullify(cpatch%dmean_norm_par_diff)
    nullify(cpatch%mmean_norm_par_diff)
    nullify(cpatch%dmean_par_v)
    nullify(cpatch%dmean_par_v_beam)
    nullify(cpatch%dmean_par_v_diff)
    nullify(cpatch%mmean_par_v)
    nullify(cpatch%mmean_par_v_beam)
    nullify(cpatch%mmean_par_v_diff)
    nullify(cpatch%par_v)
    nullify(cpatch%par_v_beam)
    nullify(cpatch%par_v_diffuse)
    nullify(cpatch%rshort_v)
    nullify(cpatch%rshort_v_beam)
    nullify(cpatch%rshort_v_diffuse)
    nullify(cpatch%rlong_v)
    nullify(cpatch%rlong_v_surf)
    nullify(cpatch%rlong_v_incid)
    nullify(cpatch%rb)
    nullify(cpatch%A_open)    
    nullify(cpatch%A_closed)
    nullify(cpatch%Psi_closed)
    nullify(cpatch%rsw_open)
    nullify(cpatch%rsw_closed)
    nullify(cpatch%fsw)
    nullify(cpatch%fs_open)
    nullify(cpatch%dmean_fs_open)
    nullify(cpatch%dmean_fsw)
    nullify(cpatch%dmean_fsn)
    nullify(cpatch%mmean_fs_open)
    nullify(cpatch%mmean_fsw)
    nullify(cpatch%mmean_fsn)
    nullify(cpatch%stomatal_resistance)
    nullify(cpatch%leaf_maintenance)
    nullify(cpatch%root_maintenance)
    nullify(cpatch%mmean_leaf_maintenance)
    nullify(cpatch%mmean_root_maintenance)
    nullify(cpatch%leaf_drop)
    nullify(cpatch%mmean_leaf_drop)
    nullify(cpatch%bseeds)
    nullify(cpatch%leaf_respiration)
    nullify(cpatch%root_respiration)
    nullify(cpatch%hcapveg)
    nullify(cpatch%gpp)
    nullify(cpatch%paw_avg)
    nullify(cpatch%turnover_amp)
    nullify(cpatch%llspan)
    nullify(cpatch%vm_bar)
    nullify(cpatch%sla)

    return
  end subroutine nullify_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_edtype(cgrid)
  ! ===================================================
  !  Deallocation subroutines of ed state variables
  ! ===================================================
    
    implicit none
    
    type(edtype),target :: cgrid

       if(associated(cgrid%polygon                 )) deallocate(cgrid%polygon                 )
       if(associated(cgrid%lat                     )) deallocate(cgrid%lat                     )
       if(associated(cgrid%lon                     )) deallocate(cgrid%lon                     )
       if(associated(cgrid%natm                    )) deallocate(cgrid%natm                    )
       if(associated(cgrid%xatm                    )) deallocate(cgrid%xatm                    )
       if(associated(cgrid%yatm                    )) deallocate(cgrid%yatm                    )
       if(associated(cgrid%ntext_soil              )) deallocate(cgrid%ntext_soil              )
       if(associated(cgrid%lsl                     )) deallocate(cgrid%lsl                     )
       
       if(associated(cgrid%pysi_n                  )) deallocate(cgrid%pysi_n                  )
       if(associated(cgrid%pysi_id                 )) deallocate(cgrid%pysi_id                 )
       if(associated(cgrid%walltime_py             )) deallocate(cgrid%walltime_py             )
       if(associated(cgrid%sensflux_py             )) deallocate(cgrid%sensflux_py             )
       if(associated(cgrid%site_adjacency          )) deallocate(cgrid%site_adjacency          )
       if(associated(cgrid%wbar                    )) deallocate(cgrid%wbar                    )
       if(associated(cgrid%Te                      )) deallocate(cgrid%Te                      )
       if(associated(cgrid%zbar                    )) deallocate(cgrid%zbar                    )
!! NOT USED       if(associated(cgrid%tau          )) deallocate(cgrid%tau                     )
       if(associated(cgrid%sheat                   )) deallocate(cgrid%sheat                   )
       if(associated(cgrid%baseflow                )) deallocate(cgrid%baseflow                )
       if(associated(cgrid%runoff                  )) deallocate(cgrid%runoff                  )
       if(associated(cgrid%swliq                   )) deallocate(cgrid%swliq                   )
       
       if(associated(cgrid%ilon                    )) deallocate(cgrid%ilon                    )
       if(associated(cgrid%ilat                    )) deallocate(cgrid%ilat                    )
       if(associated(cgrid%total_agb               )) deallocate(cgrid%total_agb               )
       if(associated(cgrid%total_basal_area        )) deallocate(cgrid%total_basal_area        )
       if(associated(cgrid%total_agb_growth        )) deallocate(cgrid%total_agb_growth        )
       if(associated(cgrid%total_agb_mort          )) deallocate(cgrid%total_agb_mort          )
       if(associated(cgrid%total_agb_recruit       )) deallocate(cgrid%total_agb_recruit       )
       if(associated(cgrid%total_basal_area_growth )) deallocate(cgrid%total_basal_area_growth )
       if(associated(cgrid%total_basal_area_mort   )) deallocate(cgrid%total_basal_area_mort   )
       if(associated(cgrid%total_basal_area_recruit)) deallocate(cgrid%total_basal_area_recruit)
       if(associated(cgrid%nsites                  )) deallocate(cgrid%nsites                  )
       if(associated(cgrid%sitenums                )) deallocate(cgrid%sitenums                )
       if(associated(cgrid%load_adjacency          )) deallocate(cgrid%load_adjacency          )
       if(associated(cgrid%cosz                    )) deallocate(cgrid%cosz                    )
       if(associated(cgrid%mean_gpp                )) deallocate(cgrid%mean_gpp                )
       if(associated(cgrid%mean_precip             )) deallocate(cgrid%mean_precip             )
       if(associated(cgrid%mean_qprecip            )) deallocate(cgrid%mean_qprecip            )
       if(associated(cgrid%mean_netrad             )) deallocate(cgrid%mean_netrad             )
       if(associated(cgrid%cbudget_initialstorage  )) deallocate(cgrid%cbudget_initialstorage  )
       if(associated(cgrid%cbudget_nep             )) deallocate(cgrid%cbudget_nep             )
       if(associated(cgrid%nbudget_initialstorage  )) deallocate(cgrid%nbudget_initialstorage  )
       if(associated(cgrid%basal_area              )) deallocate(cgrid%basal_area              )
       if(associated(cgrid%agb                     )) deallocate(cgrid%agb                     )
       if(associated(cgrid%pldens                  )) deallocate(cgrid%pldens                  )
       if(associated(cgrid%bseeds                  )) deallocate(cgrid%bseeds                  )
       
       if(associated(cgrid%metinput                )) deallocate(cgrid%metinput                )
       if(associated(cgrid%met                     )) deallocate(cgrid%met                     )
       if(associated(cgrid%lapse                   )) deallocate(cgrid%lapse                   )

       if(associated(cgrid%lai                     )) deallocate(cgrid%lai                     ) 
       if(associated(cgrid%avg_lma                 )) deallocate(cgrid%avg_lma                    )
       if(associated(cgrid%wpa                     )) deallocate(cgrid%wpa                     )
       if(associated(cgrid%wai                     )) deallocate(cgrid%wai                     )

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       if(associated(cgrid%avg_vapor_vc            )) deallocate(cgrid%avg_vapor_vc            )
       if(associated(cgrid%avg_dew_cg              )) deallocate(cgrid%avg_dew_cg              )
       if(associated(cgrid%avg_vapor_gc            )) deallocate(cgrid%avg_vapor_gc            )
       if(associated(cgrid%avg_wshed_vg            )) deallocate(cgrid%avg_wshed_vg            )
       if(associated(cgrid%avg_intercepted         )) deallocate(cgrid%avg_intercepted         )
       if(associated(cgrid%avg_vapor_ac            )) deallocate(cgrid%avg_vapor_ac            )
       if(associated(cgrid%avg_transp              )) deallocate(cgrid%avg_transp              )
       if(associated(cgrid%avg_evap                )) deallocate(cgrid%avg_evap                )
       if(associated(cgrid%avg_smoist_gg           )) deallocate(cgrid%avg_smoist_gg           )
       if(associated(cgrid%avg_smoist_gc           )) deallocate(cgrid%avg_smoist_gc           )
       if(associated(cgrid%avg_runoff              )) deallocate(cgrid%avg_runoff              )
       if(associated(cgrid%avg_drainage            )) deallocate(cgrid%avg_drainage            )
       if(associated(cgrid%avg_drainage_heat       )) deallocate(cgrid%avg_drainage_heat       )
       if(associated(cgrid%aux                     )) deallocate(cgrid%aux                     )
       if(associated(cgrid%aux_s                   )) deallocate(cgrid%aux_s                   )
       if(associated(cgrid%avg_sensible_vc         )) deallocate(cgrid%avg_sensible_vc         )
       if(associated(cgrid%avg_qwshed_vg           )) deallocate(cgrid%avg_qwshed_vg           )
       if(associated(cgrid%avg_qintercepted        )) deallocate(cgrid%avg_qintercepted        )
       if(associated(cgrid%avg_sensible_gc         )) deallocate(cgrid%avg_sensible_gc         )
       if(associated(cgrid%avg_sensible_ac         )) deallocate(cgrid%avg_sensible_ac         )
       if(associated(cgrid%avg_sensible_gg         )) deallocate(cgrid%avg_sensible_gg         )
       if(associated(cgrid%avg_runoff_heat         )) deallocate(cgrid%avg_runoff_heat         )

       if(associated(cgrid%max_veg_temp          )) deallocate(cgrid%max_veg_temp          )
       if(associated(cgrid%min_veg_temp          )) deallocate(cgrid%min_veg_temp          )
       if(associated(cgrid%max_soil_temp          )) deallocate(cgrid%max_soil_temp          )
       if(associated(cgrid%min_soil_temp          )) deallocate(cgrid%min_soil_temp          )

       ! Fast time state diagnostics
       if(associated(cgrid%avg_veg_energy          )) deallocate(cgrid%avg_veg_energy          )
       if(associated(cgrid%avg_veg_hcap            )) deallocate(cgrid%avg_veg_hcap            )
       if(associated(cgrid%avg_veg_temp            )) deallocate(cgrid%avg_veg_temp            )
       if(associated(cgrid%avg_veg_fliq            )) deallocate(cgrid%avg_veg_fliq            )
       if(associated(cgrid%avg_veg_water           )) deallocate(cgrid%avg_veg_water           )
       if(associated(cgrid%avg_can_temp            )) deallocate(cgrid%avg_can_temp            )
       if(associated(cgrid%avg_can_shv             )) deallocate(cgrid%avg_can_shv             )
       if(associated(cgrid%avg_can_co2             )) deallocate(cgrid%avg_can_co2             )
       if(associated(cgrid%avg_can_rhos            )) deallocate(cgrid%avg_can_rhos            )
       if(associated(cgrid%avg_can_prss            )) deallocate(cgrid%avg_can_prss            )
       if(associated(cgrid%avg_can_theta           )) deallocate(cgrid%avg_can_theta           )
       if(associated(cgrid%avg_can_enthalpy        )) deallocate(cgrid%avg_can_enthalpy        )
       if(associated(cgrid%avg_can_depth           )) deallocate(cgrid%avg_can_depth           )
       if(associated(cgrid%avg_soil_energy         )) deallocate(cgrid%avg_soil_energy         )
       if(associated(cgrid%avg_soil_water          )) deallocate(cgrid%avg_soil_water          )
       if(associated(cgrid%avg_soil_temp           )) deallocate(cgrid%avg_soil_temp           )
       if(associated(cgrid%avg_soil_fracliq        )) deallocate(cgrid%avg_soil_fracliq        )
       if(associated(cgrid%avg_soil_wetness        )) deallocate(cgrid%avg_soil_wetness        )
       if(associated(cgrid%avg_skin_temp           )) deallocate(cgrid%avg_skin_temp           )
       if(associated(cgrid%avg_available_water     )) deallocate(cgrid%avg_available_water     )
      
       if(associated(cgrid%avg_lai_ebalvars        )) deallocate(cgrid%avg_lai_ebalvars        )

       if(associated(cgrid%avg_gpp                 )) deallocate(cgrid%avg_gpp                 )
       if(associated(cgrid%avg_leaf_resp           )) deallocate(cgrid%avg_leaf_resp           )
       if(associated(cgrid%avg_root_resp           )) deallocate(cgrid%avg_root_resp           )
       if(associated(cgrid%avg_growth_resp         )) deallocate(cgrid%avg_growth_resp         )
       if(associated(cgrid%avg_storage_resp        )) deallocate(cgrid%avg_storage_resp        )
       if(associated(cgrid%avg_vleaf_resp          )) deallocate(cgrid%avg_vleaf_resp          )
       if(associated(cgrid%avg_plant_resp          )) deallocate(cgrid%avg_plant_resp          )
       if(associated(cgrid%avg_htroph_resp         )) deallocate(cgrid%avg_htroph_resp         )
       if(associated(cgrid%avg_leaf_drop           )) deallocate(cgrid%avg_leaf_drop           )
       if(associated(cgrid%avg_leaf_maintenance    )) deallocate(cgrid%avg_leaf_maintenance    )
       if(associated(cgrid%avg_root_maintenance    )) deallocate(cgrid%avg_root_maintenance    )


       !!! added for NACP intercomparison (MCD)
       if(associated(cgrid%avg_sfcw_depth          )) deallocate(cgrid%avg_sfcw_depth      )
       if(associated(cgrid%avg_sfcw_energy         )) deallocate(cgrid%avg_sfcw_energy     )
       if(associated(cgrid%avg_sfcw_mass           )) deallocate(cgrid%avg_sfcw_mass       )
       if(associated(cgrid%avg_sfcw_tempk          )) deallocate(cgrid%avg_sfcw_tempk      )
       if(associated(cgrid%avg_sfcw_fracliq        )) deallocate(cgrid%avg_sfcw_fracliq    )
       if(associated(cgrid%avg_bdead               )) deallocate(cgrid%avg_bdead           )
       if(associated(cgrid%avg_balive              )) deallocate(cgrid%avg_balive          )
       if(associated(cgrid%avg_broot               )) deallocate(cgrid%avg_broot           )
       if(associated(cgrid%avg_bleaf               )) deallocate(cgrid%avg_bleaf           )
       if(associated(cgrid%avg_bsapwood            )) deallocate(cgrid%avg_bsapwood        )
       if(associated(cgrid%avg_bstorage            )) deallocate(cgrid%avg_bstorage        )
       if(associated(cgrid%avg_bseeds              )) deallocate(cgrid%avg_bseeds          )
       if(associated(cgrid%avg_fsc                 )) deallocate(cgrid%avg_fsc             )
       if(associated(cgrid%avg_ssc                 )) deallocate(cgrid%avg_ssc             )
       if(associated(cgrid%avg_stsc                )) deallocate(cgrid%avg_stsc            )

       if(associated(cgrid%avg_fsn                 )) deallocate(cgrid%avg_fsn             )
       if(associated(cgrid%avg_msn                 )) deallocate(cgrid%avg_msn             )


     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf  )) deallocate(cgrid%Cleaf)
     if(associated(cgrid%Croot  )) deallocate(cgrid%Croot)
     if(associated(cgrid%Cstore )) deallocate(cgrid%Cstore)
     if(associated(cgrid%Ccwd   )) deallocate(cgrid%Ccwd)
     if(associated(cgrid%Nleaf  )) deallocate(cgrid%Nleaf)
     if(associated(cgrid%Ndead  )) deallocate(cgrid%Ndead)
     if(associated(cgrid%Nroot  )) deallocate(cgrid%Nroot)
     if(associated(cgrid%Nstore )) deallocate(cgrid%Nstore)
     if(associated(cgrid%Ncwd   )) deallocate(cgrid%Ncwd)
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf_grow       )) deallocate(cgrid%Cleaf_grow)
     if(associated(cgrid%Croot_grow       )) deallocate(cgrid%Croot_grow)
     if(associated(cgrid%Cdead_grow       )) deallocate(cgrid%Cdead_grow)
     if(associated(cgrid%Cstore_grow      )) deallocate(cgrid%Cstore_grow)
     if(associated(cgrid%Cleaf_litter_flux)) deallocate(cgrid%Cleaf_litter_flux)
     if(associated(cgrid%Croot_litter_flux)) deallocate(cgrid%Croot_litter_flux)
     if(associated(cgrid%Ccwd_flux        )) deallocate(cgrid%Ccwd_flux)
     if(associated(cgrid%Nleaf_grow       )) deallocate(cgrid%Nleaf_grow)
     if(associated(cgrid%Ndead_grow       )) deallocate(cgrid%Ndead_grow)
     if(associated(cgrid%Nroot_grow       )) deallocate(cgrid%Nroot_grow)
     if(associated(cgrid%Nstore_grow      )) deallocate(cgrid%Nstore_grow)
     if(associated(cgrid%Nleaf_litter_flux)) deallocate(cgrid%Nleaf_litter_flux)
     if(associated(cgrid%Nroot_litter_flux)) deallocate(cgrid%Nroot_litter_flux)
     if(associated(cgrid%Ncwd_flux        )) deallocate(cgrid%Ncwd_flux)
     if(associated(cgrid%Nbiomass_uptake  )) deallocate(cgrid%Nbiomass_uptake)
     if(associated(cgrid%Ngross_min       )) deallocate(cgrid%Ngross_min)
     if(associated(cgrid%Nnet_min         )) deallocate(cgrid%Nnet_min)


 


       ! ----------------------------------------------

       if(associated(cgrid%avg_nir_beam            )) deallocate(cgrid%avg_nir_beam            )
       if(associated(cgrid%avg_nir_diffuse         )) deallocate(cgrid%avg_nir_diffuse         )
       if(associated(cgrid%avg_par_beam            )) deallocate(cgrid%avg_par_beam            )
       if(associated(cgrid%avg_par_diffuse         )) deallocate(cgrid%avg_par_diffuse         )
       if(associated(cgrid%avg_atm_tmp             )) deallocate(cgrid%avg_atm_tmp             )
       if(associated(cgrid%avg_atm_shv             )) deallocate(cgrid%avg_atm_shv             )
       if(associated(cgrid%avg_rshort              )) deallocate(cgrid%avg_rshort              )
       if(associated(cgrid%avg_rshort_diffuse      )) deallocate(cgrid%avg_rshort_diffuse      )
       if(associated(cgrid%avg_rlong               )) deallocate(cgrid%avg_rlong               )
       if(associated(cgrid%avg_pcpg                )) deallocate(cgrid%avg_pcpg                )
       if(associated(cgrid%avg_qpcpg               )) deallocate(cgrid%avg_qpcpg               )
       if(associated(cgrid%avg_dpcpg               )) deallocate(cgrid%avg_dpcpg               )
       if(associated(cgrid%avg_vels                )) deallocate(cgrid%avg_vels                )
       if(associated(cgrid%avg_atm_prss            )) deallocate(cgrid%avg_atm_prss            )
       if(associated(cgrid%avg_exner               )) deallocate(cgrid%avg_exner               )
       if(associated(cgrid%avg_geoht               )) deallocate(cgrid%avg_geoht               )
       if(associated(cgrid%avg_atm_co2             )) deallocate(cgrid%avg_atm_co2             )
       if(associated(cgrid%avg_albedt              )) deallocate(cgrid%avg_albedt              )
       if(associated(cgrid%avg_rlongup             )) deallocate(cgrid%avg_rlongup             )


       if(associated(cgrid%runoff                  )) deallocate(cgrid%runoff                  )
       
       if(associated(cgrid%dmean_pcpg              )) deallocate(cgrid%dmean_pcpg              )
       if(associated(cgrid%dmean_runoff            )) deallocate(cgrid%dmean_runoff            )
       if(associated(cgrid%dmean_drainage          )) deallocate(cgrid%dmean_drainage          )
       if(associated(cgrid%dmean_vapor_ac          )) deallocate(cgrid%dmean_vapor_ac          )
       if(associated(cgrid%dmean_vapor_gc          )) deallocate(cgrid%dmean_vapor_gc          )
       if(associated(cgrid%dmean_vapor_vc          )) deallocate(cgrid%dmean_vapor_vc          )
       if(associated(cgrid%dmean_gpp               )) deallocate(cgrid%dmean_gpp               )
       if(associated(cgrid%dmean_evap              )) deallocate(cgrid%dmean_evap              )
       if(associated(cgrid%dmean_transp            )) deallocate(cgrid%dmean_transp            )
       if(associated(cgrid%dmean_sensible_vc       )) deallocate(cgrid%dmean_sensible_vc       )
       if(associated(cgrid%dmean_sensible_gc       )) deallocate(cgrid%dmean_sensible_gc       )
       if(associated(cgrid%dmean_sensible_ac       )) deallocate(cgrid%dmean_sensible_ac       )
       if(associated(cgrid%dmean_plresp            )) deallocate(cgrid%dmean_plresp            )
       if(associated(cgrid%dmean_rh                )) deallocate(cgrid%dmean_rh                )
       if(associated(cgrid%dmean_leaf_resp         )) deallocate(cgrid%dmean_leaf_resp         )
       if(associated(cgrid%dmean_root_resp         )) deallocate(cgrid%dmean_root_resp         )
       if(associated(cgrid%dmean_growth_resp       )) deallocate(cgrid%dmean_growth_resp       )
       if(associated(cgrid%dmean_storage_resp      )) deallocate(cgrid%dmean_storage_resp      )
       if(associated(cgrid%dmean_vleaf_resp        )) deallocate(cgrid%dmean_vleaf_resp        )
       if(associated(cgrid%dmean_nep               )) deallocate(cgrid%dmean_nep               )
       if(associated(cgrid%dmean_soil_temp         )) deallocate(cgrid%dmean_soil_temp         )
       if(associated(cgrid%dmean_soil_water        )) deallocate(cgrid%dmean_soil_water        )
       if(associated(cgrid%dmean_fs_open           )) deallocate(cgrid%dmean_fs_open           )
       if(associated(cgrid%dmean_fsw               )) deallocate(cgrid%dmean_fsw               )
       if(associated(cgrid%dmean_fsn               )) deallocate(cgrid%dmean_fsn               )
       if(associated(cgrid%dmean_gpp_dbh           )) deallocate(cgrid%dmean_gpp_dbh           )
       if(associated(cgrid%dmean_can_temp          )) deallocate(cgrid%dmean_can_temp          )
       if(associated(cgrid%dmean_can_shv           )) deallocate(cgrid%dmean_can_shv           )
       if(associated(cgrid%dmean_can_co2           )) deallocate(cgrid%dmean_can_co2           )
       if(associated(cgrid%dmean_can_rhos          )) deallocate(cgrid%dmean_can_rhos          )
       if(associated(cgrid%dmean_can_prss          )) deallocate(cgrid%dmean_can_prss          )
       if(associated(cgrid%dmean_can_theta         )) deallocate(cgrid%dmean_can_theta         )
       if(associated(cgrid%dmean_veg_energy        )) deallocate(cgrid%dmean_veg_energy        )
       if(associated(cgrid%dmean_veg_water         )) deallocate(cgrid%dmean_veg_water         )
       if(associated(cgrid%dmean_veg_hcap          )) deallocate(cgrid%dmean_veg_hcap          )
       if(associated(cgrid%dmean_veg_temp          )) deallocate(cgrid%dmean_veg_temp          )
       if(associated(cgrid%dmean_atm_temp          )) deallocate(cgrid%dmean_atm_temp          )
       if(associated(cgrid%dmean_atm_shv           )) deallocate(cgrid%dmean_atm_shv           )
       if(associated(cgrid%dmean_atm_prss          )) deallocate(cgrid%dmean_atm_prss          )
       if(associated(cgrid%dmean_atm_vels          )) deallocate(cgrid%dmean_atm_vels          )
       if(associated(cgrid%dmean_rshort            )) deallocate(cgrid%dmean_rshort            )
       if(associated(cgrid%dmean_rlong             )) deallocate(cgrid%dmean_rlong             )
       if(associated(cgrid%dmean_co2_residual      )) deallocate(cgrid%dmean_co2_residual      )
       if(associated(cgrid%dmean_energy_residual   )) deallocate(cgrid%dmean_energy_residual   )
       if(associated(cgrid%dmean_water_residual    )) deallocate(cgrid%dmean_water_residual    )

       if(associated(cgrid%lai_pft                 )) deallocate(cgrid%lai_pft                 )
       if(associated(cgrid%wpa_pft                 )) deallocate(cgrid%wpa_pft                 )
       if(associated(cgrid%wai_pft                 )) deallocate(cgrid%wai_pft                 )
       if(associated(cgrid%mmean_gpp               )) deallocate(cgrid%mmean_gpp               )
       if(associated(cgrid%mmean_evap              )) deallocate(cgrid%mmean_evap              )
       if(associated(cgrid%mmean_transp            )) deallocate(cgrid%mmean_transp            )
       if(associated(cgrid%mmean_vapor_vc          )) deallocate(cgrid%mmean_vapor_vc          )
       if(associated(cgrid%mmean_vapor_gc          )) deallocate(cgrid%mmean_vapor_gc          )
       if(associated(cgrid%mmean_vapor_ac          )) deallocate(cgrid%mmean_vapor_ac          )
       if(associated(cgrid%mmean_sensible_vc       )) deallocate(cgrid%mmean_sensible_vc       )
       if(associated(cgrid%mmean_sensible_gc       )) deallocate(cgrid%mmean_sensible_gc       )
       if(associated(cgrid%mmean_sensible_ac       )) deallocate(cgrid%mmean_sensible_ac       )
       if(associated(cgrid%mmean_nep               )) deallocate(cgrid%mmean_nep               )
       if(associated(cgrid%mmean_soil_temp         )) deallocate(cgrid%mmean_soil_temp         )
       if(associated(cgrid%mmean_soil_water        )) deallocate(cgrid%mmean_soil_water        )
       if(associated(cgrid%mmean_fs_open           )) deallocate(cgrid%mmean_fs_open           )
       if(associated(cgrid%mmean_fsw               )) deallocate(cgrid%mmean_fsw               )
       if(associated(cgrid%mmean_fsn               )) deallocate(cgrid%mmean_fsn               )
       if(associated(cgrid%mmean_plresp            )) deallocate(cgrid%mmean_plresp            )
       if(associated(cgrid%mmean_rh                )) deallocate(cgrid%mmean_rh                )
       if(associated(cgrid%mmean_leaf_resp         )) deallocate(cgrid%mmean_leaf_resp         )
       if(associated(cgrid%mmean_root_resp         )) deallocate(cgrid%mmean_root_resp         )
       if(associated(cgrid%mmean_growth_resp       )) deallocate(cgrid%mmean_growth_resp       )
       if(associated(cgrid%mmean_storage_resp      )) deallocate(cgrid%mmean_storage_resp      )
       if(associated(cgrid%mmean_vleaf_resp        )) deallocate(cgrid%mmean_vleaf_resp        )
       if(associated(cgrid%mmean_gpp_dbh           )) deallocate(cgrid%mmean_gpp_dbh           )
       if(associated(cgrid%mmean_lai_pft           )) deallocate(cgrid%mmean_lai_pft           )
       if(associated(cgrid%mmean_wpa_pft           )) deallocate(cgrid%mmean_wpa_pft           )
       if(associated(cgrid%mmean_wai_pft           )) deallocate(cgrid%mmean_wai_pft           )
       if(associated(cgrid%mmean_can_temp          )) deallocate(cgrid%mmean_can_temp          )
       if(associated(cgrid%mmean_can_shv           )) deallocate(cgrid%mmean_can_shv           )
       if(associated(cgrid%mmean_can_co2           )) deallocate(cgrid%mmean_can_co2           )
       if(associated(cgrid%mmean_can_rhos          )) deallocate(cgrid%mmean_can_rhos          )
       if(associated(cgrid%mmean_can_prss          )) deallocate(cgrid%mmean_can_prss          )
       if(associated(cgrid%mmean_can_theta         )) deallocate(cgrid%mmean_can_theta         )
       if(associated(cgrid%mmean_veg_energy        )) deallocate(cgrid%mmean_veg_energy        )
       if(associated(cgrid%mmean_veg_water         )) deallocate(cgrid%mmean_veg_water         )
       if(associated(cgrid%mmean_veg_hcap          )) deallocate(cgrid%mmean_veg_hcap          )
       if(associated(cgrid%mmean_veg_temp          )) deallocate(cgrid%mmean_veg_temp          )
       if(associated(cgrid%mmean_atm_temp          )) deallocate(cgrid%mmean_atm_temp          )
       if(associated(cgrid%mmean_rshort            )) deallocate(cgrid%mmean_rshort            )
       if(associated(cgrid%mmean_rlong             )) deallocate(cgrid%mmean_rlong             )
       if(associated(cgrid%mmean_atm_shv           )) deallocate(cgrid%mmean_atm_shv           )
       if(associated(cgrid%mmean_atm_prss          )) deallocate(cgrid%mmean_atm_prss          )
       if(associated(cgrid%mmean_atm_vels          )) deallocate(cgrid%mmean_atm_vels          )
       if(associated(cgrid%mmean_pcpg              )) deallocate(cgrid%mmean_pcpg              )
       if(associated(cgrid%mmean_runoff            )) deallocate(cgrid%mmean_runoff            )
       if(associated(cgrid%mmean_drainage          )) deallocate(cgrid%mmean_drainage          )
       if(associated(cgrid%mmean_co2_residual      )) deallocate(cgrid%mmean_co2_residual      )
       if(associated(cgrid%mmean_energy_residual   )) deallocate(cgrid%mmean_energy_residual   )
       if(associated(cgrid%mmean_water_residual    )) deallocate(cgrid%mmean_water_residual    )
       if(associated(cgrid%bseeds_pft              )) deallocate(cgrid%bseeds_pft              )
       if(associated(cgrid%agb_pft                 )) deallocate(cgrid%agb_pft                 )
       if(associated(cgrid%ba_pft                  )) deallocate(cgrid%ba_pft                  )
       if(associated(cgrid%stdev_gpp               )) deallocate(cgrid%stdev_gpp               )
       if(associated(cgrid%stdev_evap              )) deallocate(cgrid%stdev_evap              )
       if(associated(cgrid%stdev_transp            )) deallocate(cgrid%stdev_transp            )
       if(associated(cgrid%stdev_sensible          )) deallocate(cgrid%stdev_sensible          )
       if(associated(cgrid%stdev_nep               )) deallocate(cgrid%stdev_nep               )
       if(associated(cgrid%stdev_rh                )) deallocate(cgrid%stdev_rh                )
       if(associated(cgrid%disturbance_rates       )) deallocate(cgrid%disturbance_rates       )
       if(associated(cgrid%workload                )) deallocate(cgrid%workload                )

    return
  end subroutine deallocate_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_polygontype(cpoly)

    implicit none

    type(polygontype),target :: cpoly

    if(associated(cpoly%sipa_id                     )) deallocate(cpoly%sipa_id                     )
    if(associated(cpoly%sipa_n                      )) deallocate(cpoly%sipa_n                      )
    if(associated(cpoly%patch_count                 )) deallocate(cpoly%patch_count                 )
    if(associated(cpoly%site                        )) deallocate(cpoly%site                        )
    if(associated(cpoly%sitenum                     )) deallocate(cpoly%sitenum                     )

    if(associated(cpoly%lsl                         )) deallocate(cpoly%lsl                         )

    if(associated(cpoly%area                        )) deallocate(cpoly%area                        )
    if(associated(cpoly%patch_area                  )) deallocate(cpoly%patch_area                  )
    if(associated(cpoly%elevation                   )) deallocate(cpoly%elevation                   )
    if(associated(cpoly%slope                       )) deallocate(cpoly%slope                       )
    if(associated(cpoly%aspect                      )) deallocate(cpoly%aspect                      )

    if(associated(cpoly%num_landuse_years           )) deallocate(cpoly%num_landuse_years           )
    if(associated(cpoly%lai_pft                     )) deallocate(cpoly%lai_pft                     )
    if(associated(cpoly%wpa_pft                     )) deallocate(cpoly%wpa_pft                     )
    if(associated(cpoly%wai_pft                     )) deallocate(cpoly%wai_pft                     )

    if(associated(cpoly%TCI                         )) deallocate(cpoly%TCI                         )
    if(associated(cpoly%pptweight                   )) deallocate(cpoly%pptweight                         )
    if(associated(cpoly%lsl                         )) deallocate(cpoly%lsl                         )
    if(associated(cpoly%hydro_next                  )) deallocate(cpoly%hydro_next                  )
    if(associated(cpoly%hydro_prev                  )) deallocate(cpoly%hydro_prev                  )
    if(associated(cpoly%moist_W                     )) deallocate(cpoly%moist_W                     )
    if(associated(cpoly%moist_f                     )) deallocate(cpoly%moist_f                     )
    if(associated(cpoly%moist_tau                   )) deallocate(cpoly%moist_tau                   )
    if(associated(cpoly%moist_zi                    )) deallocate(cpoly%moist_zi                    )
    if(associated(cpoly%baseflow                    )) deallocate(cpoly%baseflow                    )
    if(associated(cpoly%ntext_soil                  )) deallocate(cpoly%ntext_soil                  )
    if(associated(cpoly%min_monthly_temp            )) deallocate(cpoly%min_monthly_temp            )
    if(associated(cpoly%plantation                  )) deallocate(cpoly%plantation                  )
    if(associated(cpoly%agri_stocking_pft           )) deallocate(cpoly%agri_stocking_pft           )
    if(associated(cpoly%agri_stocking_density       )) deallocate(cpoly%agri_stocking_density       )
    if(associated(cpoly%plantation_stocking_pft     )) deallocate(cpoly%plantation_stocking_pft     )
    if(associated(cpoly%plantation_stocking_density )) deallocate(cpoly%plantation_stocking_density )
    if(associated(cpoly%primary_harvest_memory      )) deallocate(cpoly%primary_harvest_memory      )
    if(associated(cpoly%secondary_harvest_memory    )) deallocate(cpoly%secondary_harvest_memory    )
    if(associated(cpoly%fire_disturbance_rate       )) deallocate(cpoly%fire_disturbance_rate       )
    if(associated(cpoly%ignition_rate               )) deallocate(cpoly%ignition_rate               )
    if(associated(cpoly%lambda_fire                 )) deallocate(cpoly%lambda_fire                 )
    if(associated(cpoly%phen_pars                   )) deallocate(cpoly%phen_pars                   )
    if(associated(cpoly%treefall_disturbance_rate   )) deallocate(cpoly%treefall_disturbance_rate   )
    if(associated(cpoly%nat_disturbance_rate        )) deallocate(cpoly%nat_disturbance_rate        )
    if(associated(cpoly%nat_dist_type               )) deallocate(cpoly%nat_dist_type               )
    if(associated(cpoly%disturbance_memory          )) deallocate(cpoly%disturbance_memory          )
    if(associated(cpoly%disturbance_rates           )) deallocate(cpoly%disturbance_rates           )
    if(associated(cpoly%loss_fraction               )) deallocate(cpoly%loss_fraction               )
    if(associated(cpoly%green_leaf_factor           )) deallocate(cpoly%green_leaf_factor           )
    if(associated(cpoly%leaf_aging_factor           )) deallocate(cpoly%leaf_aging_factor           )
    if(associated(cpoly%met                         )) deallocate(cpoly%met                         )
    if(associated(cpoly%basal_area                  )) deallocate(cpoly%basal_area                  )
    if(associated(cpoly%agb                         )) deallocate(cpoly%agb                         )
    if(associated(cpoly%pldens                      )) deallocate(cpoly%pldens                      )
    if(associated(cpoly%bseeds                      )) deallocate(cpoly%bseeds                      )
    
    if(associated(cpoly%basal_area_growth           )) deallocate(cpoly%basal_area_growth           )
    if(associated(cpoly%agb_growth                  )) deallocate(cpoly%agb_growth                  )
    if(associated(cpoly%basal_area_mort             )) deallocate(cpoly%basal_area_mort             )
    if(associated(cpoly%agb_mort                    )) deallocate(cpoly%agb_mort                    )
    if(associated(cpoly%basal_area_cut              )) deallocate(cpoly%basal_area_cut              )!NOT IN REGISTRY
    if(associated(cpoly%agb_cut                     )) deallocate(cpoly%agb_cut                     )

    if(associated(cpoly%cosaoi                      )) deallocate(cpoly%cosaoi                      )
    if(associated(cpoly%albedo_beam                 )) deallocate(cpoly%albedo_beam                 )
    if(associated(cpoly%albedo_diffuse              )) deallocate(cpoly%albedo_diffuse              )
    if(associated(cpoly%rlong_albedo                )) deallocate(cpoly%rlong_albedo                )
    
    if(associated(cpoly%albedt                      )) deallocate(cpoly%albedt                      )
    if(associated(cpoly%rlongup                     )) deallocate(cpoly%rlongup                     )
    if(associated(cpoly%daylight                    )) deallocate(cpoly%daylight                    )
    
    if(associated(cpoly%lai                         )) deallocate(cpoly%lai                         )
    if(associated(cpoly%avg_lma                     )) deallocate(cpoly%avg_lma                         )
    if(associated(cpoly%wpa                         )) deallocate(cpoly%wpa                         )
    if(associated(cpoly%wai                         )) deallocate(cpoly%wai                         )

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    if(associated(cpoly%avg_vapor_vc                )) deallocate(cpoly%avg_vapor_vc                )
    if(associated(cpoly%avg_dew_cg                  )) deallocate(cpoly%avg_dew_cg                  )
    if(associated(cpoly%avg_vapor_gc                )) deallocate(cpoly%avg_vapor_gc                )
    if(associated(cpoly%avg_wshed_vg                )) deallocate(cpoly%avg_wshed_vg                )
    if(associated(cpoly%avg_intercepted             )) deallocate(cpoly%avg_intercepted             )
    if(associated(cpoly%avg_vapor_ac                )) deallocate(cpoly%avg_vapor_ac                )
    if(associated(cpoly%avg_transp                  )) deallocate(cpoly%avg_transp                  )
    if(associated(cpoly%avg_evap                    )) deallocate(cpoly%avg_evap                    )
    if(associated(cpoly%avg_smoist_gg               )) deallocate(cpoly%avg_smoist_gg               )
    if(associated(cpoly%avg_smoist_gc               )) deallocate(cpoly%avg_smoist_gc               )
    if(associated(cpoly%avg_runoff                  )) deallocate(cpoly%avg_runoff                  )
    if(associated(cpoly%avg_drainage                )) deallocate(cpoly%avg_drainage                )
    if(associated(cpoly%avg_drainage_heat           )) deallocate(cpoly%avg_drainage_heat           )
    if(associated(cpoly%aux                         )) deallocate(cpoly%aux                         )
    if(associated(cpoly%aux_s                       )) deallocate(cpoly%aux_s                       )
    if(associated(cpoly%avg_sensible_vc             )) deallocate(cpoly%avg_sensible_vc             )
    if(associated(cpoly%avg_qwshed_vg               )) deallocate(cpoly%avg_qwshed_vg               )
    if(associated(cpoly%avg_qintercepted            )) deallocate(cpoly%avg_qintercepted            )
    if(associated(cpoly%avg_sensible_gc             )) deallocate(cpoly%avg_sensible_gc             )
    if(associated(cpoly%avg_sensible_ac             )) deallocate(cpoly%avg_sensible_ac             )
    if(associated(cpoly%avg_sensible_gg             )) deallocate(cpoly%avg_sensible_gg             )
    if(associated(cpoly%avg_runoff_heat             )) deallocate(cpoly%avg_runoff_heat             )
    if(associated(cpoly%avg_veg_energy              )) deallocate(cpoly%avg_veg_energy              )
    if(associated(cpoly%avg_veg_hcap                )) deallocate(cpoly%avg_veg_hcap                )
    if(associated(cpoly%avg_veg_temp                )) deallocate(cpoly%avg_veg_temp                )
    if(associated(cpoly%avg_veg_fliq                )) deallocate(cpoly%avg_veg_fliq                )
    if(associated(cpoly%avg_veg_water               )) deallocate(cpoly%avg_veg_water               )
    if(associated(cpoly%avg_can_temp                )) deallocate(cpoly%avg_can_temp                )
    if(associated(cpoly%avg_can_shv                 )) deallocate(cpoly%avg_can_shv                 )
    if(associated(cpoly%avg_can_co2                 )) deallocate(cpoly%avg_can_co2                 )
    if(associated(cpoly%avg_can_rhos                )) deallocate(cpoly%avg_can_rhos                )
    if(associated(cpoly%avg_can_prss                )) deallocate(cpoly%avg_can_prss                )
    if(associated(cpoly%avg_can_theta               )) deallocate(cpoly%avg_can_theta               )
    if(associated(cpoly%avg_can_enthalpy            )) deallocate(cpoly%avg_can_enthalpy            )
    if(associated(cpoly%avg_can_depth               )) deallocate(cpoly%avg_can_depth               )
    if(associated(cpoly%avg_soil_energy             )) deallocate(cpoly%avg_soil_energy             )
    if(associated(cpoly%avg_soil_water              )) deallocate(cpoly%avg_soil_water              )
    if(associated(cpoly%avg_soil_temp               )) deallocate(cpoly%avg_soil_temp               )
    if(associated(cpoly%avg_soil_fracliq            )) deallocate(cpoly%avg_soil_fracliq            )
    if(associated(cpoly%avg_soil_wetness            )) deallocate(cpoly%avg_soil_wetness            )
    if(associated(cpoly%avg_skin_temp               )) deallocate(cpoly%avg_skin_temp               )
    if(associated(cpoly%avg_available_water         )) deallocate(cpoly%avg_available_water         )
    if(associated(cpoly%runoff                      )) deallocate(cpoly%runoff                      )
    if(associated(cpoly%rad_avg                     )) deallocate(cpoly%rad_avg                     )
    ! Meteorological information
    if(associated(cpoly%avg_atm_tmp                 )) deallocate(cpoly%avg_atm_tmp                 )
    if(associated(cpoly%avg_atm_shv                 )) deallocate(cpoly%avg_atm_shv                 )
    if(associated(cpoly%avg_atm_prss                )) deallocate(cpoly%avg_atm_prss                )
    ! NACP
    if(associated(cpoly%avg_sfcw_depth            )) deallocate(cpoly%avg_sfcw_depth            )
    if(associated(cpoly%avg_sfcw_energy           )) deallocate(cpoly%avg_sfcw_energy           )
    if(associated(cpoly%avg_sfcw_mass             )) deallocate(cpoly%avg_sfcw_mass             )
    if(associated(cpoly%avg_sfcw_fracliq          )) deallocate(cpoly%avg_sfcw_fracliq          )
    if(associated(cpoly%avg_sfcw_tempk            )) deallocate(cpoly%avg_sfcw_tempk            )
    if(associated(cpoly%avg_fsc                   )) deallocate(cpoly%avg_fsc                   )
    if(associated(cpoly%avg_stsc                  )) deallocate(cpoly%avg_stsc                  )
    if(associated(cpoly%avg_ssc                   )) deallocate(cpoly%avg_ssc                   )
    if(associated(cpoly%avg_bdead                 )) deallocate(cpoly%avg_bdead                 )
    if(associated(cpoly%avg_balive                )) deallocate(cpoly%avg_balive                )
    if(associated(cpoly%avg_fsn                   )) deallocate(cpoly%avg_fsn                   )
    if(associated(cpoly%avg_msn                   )) deallocate(cpoly%avg_msn                   )

    if(associated(cpoly%avg_bleaf                 )) deallocate(cpoly%avg_bleaf                 )
    if(associated(cpoly%avg_broot                 )) deallocate(cpoly%avg_broot                 )
    if(associated(cpoly%avg_bsapwood              )) deallocate(cpoly%avg_bsapwood              )
    if(associated(cpoly%avg_bstorage              )) deallocate(cpoly%avg_bstorage              )
    if(associated(cpoly%avg_bseeds                )) deallocate(cpoly%avg_bseeds                )

    if(associated(cpoly%dmean_co2_residual        )) deallocate(cpoly%dmean_co2_residual        )
    if(associated(cpoly%dmean_energy_residual     )) deallocate(cpoly%dmean_energy_residual     )
    if(associated(cpoly%dmean_water_residual      )) deallocate(cpoly%dmean_water_residual      )
    if(associated(cpoly%mmean_co2_residual        )) deallocate(cpoly%mmean_co2_residual        )
    if(associated(cpoly%mmean_energy_residual     )) deallocate(cpoly%mmean_energy_residual     )
    if(associated(cpoly%mmean_water_residual      )) deallocate(cpoly%mmean_water_residual      )

    return
  end subroutine deallocate_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_sitetype(csite)

    implicit none

    type(sitetype),target :: csite
    integer :: ipa
    
    if(associated(csite%paco_id                      )) deallocate(csite%paco_id                      )
    if(associated(csite%paco_n                       )) deallocate(csite%paco_n                       )

    if(associated(csite%dist_type                    )) deallocate(csite%dist_type                    )
    if(associated(csite%age                          )) deallocate(csite%age                          )
    if(associated(csite%area                         )) deallocate(csite%area                         )
    if(associated(csite%laiarea                      )) deallocate(csite%laiarea                      )
    if(associated(csite%hcapveg                      )) deallocate(csite%hcapveg                      )
    if(associated(csite%fast_soil_C                  )) deallocate(csite%fast_soil_C                  )
    if(associated(csite%slow_soil_C                  )) deallocate(csite%slow_soil_C                  )
    if(associated(csite%structural_soil_C            )) deallocate(csite%structural_soil_C            )
    if(associated(csite%structural_soil_L            )) deallocate(csite%structural_soil_L            )
    if(associated(csite%mineralized_soil_N           )) deallocate(csite%mineralized_soil_N           )
    if(associated(csite%fast_soil_N                  )) deallocate(csite%fast_soil_N                  )
    if(associated(csite%pname                        )) deallocate(csite%pname                        )
    if(associated(csite%sum_dgd                      )) deallocate(csite%sum_dgd                      )
    if(associated(csite%sum_chd                      )) deallocate(csite%sum_chd                      )
    if(associated(csite%plantation                   )) deallocate(csite%plantation                   )
    if(associated(csite%cohort_count                 )) deallocate(csite%cohort_count                 )
    if(associated(csite%can_enthalpy                 )) deallocate(csite%can_enthalpy                 )
    if(associated(csite%can_temp                     )) deallocate(csite%can_temp                     )
    if(associated(csite%can_shv                      )) deallocate(csite%can_shv                      )
    if(associated(csite%can_co2                      )) deallocate(csite%can_co2                      )
    if(associated(csite%can_rhos                     )) deallocate(csite%can_rhos                     )
    if(associated(csite%can_prss                     )) deallocate(csite%can_prss                     )
    if(associated(csite%can_theta                    )) deallocate(csite%can_theta                    )
    if(associated(csite%can_depth                    )) deallocate(csite%can_depth                    )
    if(associated(csite%lambda_light                 )) deallocate(csite%lambda_light                 )
    if(associated(csite%dmean_lambda_light           )) deallocate(csite%dmean_lambda_light           )
    if(associated(csite%mmean_lambda_light           )) deallocate(csite%mmean_lambda_light           )
    if(associated(csite%lai                          )) deallocate(csite%lai                          )
    if(associated(csite%wpa                          )) deallocate(csite%wpa                          )
    if(associated(csite%wai                          )) deallocate(csite%wai                          )

    if(associated(csite%sfcwater_mass                )) deallocate(csite%sfcwater_mass                )
    if(associated(csite%sfcwater_energy              )) deallocate(csite%sfcwater_energy              )
    if(associated(csite%sfcwater_depth               )) deallocate(csite%sfcwater_depth               )
    if(associated(csite%rshort_s                     )) deallocate(csite%rshort_s                     )
    if(associated(csite%rshort_s_beam                )) deallocate(csite%rshort_s_beam                )
    if(associated(csite%rshort_s_diffuse             )) deallocate(csite%rshort_s_diffuse             )
    if(associated(csite%sfcwater_tempk               )) deallocate(csite%sfcwater_tempk               )
    if(associated(csite%sfcwater_fracliq             )) deallocate(csite%sfcwater_fracliq             )
    if(associated(csite%nlev_sfcwater                )) deallocate(csite%nlev_sfcwater                )
    if(associated(csite%ntext_soil                   )) deallocate(csite%ntext_soil                   )
    if(associated(csite%soil_energy                  )) deallocate(csite%soil_energy                  )
    if(associated(csite%soil_water                   )) deallocate(csite%soil_water                   )
    if(associated(csite%soil_tempk                   )) deallocate(csite%soil_tempk                   )
    if(associated(csite%soil_fracliq                 )) deallocate(csite%soil_fracliq                 )
    if(associated(csite%ground_shv                   )) deallocate(csite%ground_shv                   )
    if(associated(csite%surface_ssh                  )) deallocate(csite%surface_ssh                  )
    if(associated(csite%rough                        )) deallocate(csite%rough                        )
    if(associated(csite%A_o_max                      )) deallocate(csite%A_o_max                      )
    if(associated(csite%A_c_max                      )) deallocate(csite%A_c_max                      )

    if(associated(csite%old_stoma_data_max           )) deallocate(csite%old_stoma_data_max           )
    if(associated(csite%old_stoma_vector_max         )) deallocate(csite%old_stoma_vector_max         )

    if(associated(csite%avg_daily_temp               )) deallocate(csite%avg_daily_temp               )
    if(associated(csite%mean_rh                      )) deallocate(csite%mean_rh                      )
    if(associated(csite%dmean_rh                     )) deallocate(csite%dmean_rh                     )
    if(associated(csite%mmean_rh                     )) deallocate(csite%mmean_rh                     )
    if(associated(csite%mean_nep                     )) deallocate(csite%mean_nep                     )
    if(associated(csite%wbudget_loss2atm             )) deallocate(csite%wbudget_loss2atm             )
    if(associated(csite%wbudget_denseffect           )) deallocate(csite%wbudget_denseffect           )
    if(associated(csite%wbudget_precipgain           )) deallocate(csite%wbudget_precipgain           )
    if(associated(csite%wbudget_loss2runoff          )) deallocate(csite%wbudget_loss2runoff          )
    if(associated(csite%wbudget_loss2drainage        )) deallocate(csite%wbudget_loss2drainage        )
    if(associated(csite%wbudget_initialstorage       )) deallocate(csite%wbudget_initialstorage       )
    if(associated(csite%wbudget_residual             )) deallocate(csite%wbudget_residual             )
    if(associated(csite%ebudget_latent               )) deallocate(csite%ebudget_latent               )
    if(associated(csite%ebudget_loss2atm             )) deallocate(csite%ebudget_loss2atm             )
    if(associated(csite%ebudget_denseffect           )) deallocate(csite%ebudget_denseffect           )
    if(associated(csite%ebudget_loss2runoff          )) deallocate(csite%ebudget_loss2runoff          )
    if(associated(csite%ebudget_loss2drainage        )) deallocate(csite%ebudget_loss2drainage        )
    if(associated(csite%ebudget_netrad               )) deallocate(csite%ebudget_netrad               )
    if(associated(csite%ebudget_precipgain           )) deallocate(csite%ebudget_precipgain           )
    if(associated(csite%ebudget_initialstorage       )) deallocate(csite%ebudget_initialstorage       )
    if(associated(csite%ebudget_residual             )) deallocate(csite%ebudget_residual             )
    if(associated(csite%co2budget_initialstorage     )) deallocate(csite%co2budget_initialstorage     )
    if(associated(csite%co2budget_residual           )) deallocate(csite%co2budget_residual           )
    if(associated(csite%co2budget_loss2atm           )) deallocate(csite%co2budget_loss2atm           )
    if(associated(csite%co2budget_denseffect         )) deallocate(csite%co2budget_denseffect         )
    if(associated(csite%co2budget_gpp                )) deallocate(csite%co2budget_gpp                )
    if(associated(csite%co2budget_gpp_dbh            )) deallocate(csite%co2budget_gpp_dbh            )
    if(associated(csite%co2budget_plresp             )) deallocate(csite%co2budget_plresp             )
    if(associated(csite%co2budget_rh                 )) deallocate(csite%co2budget_rh                 )
    if(associated(csite%today_A_decomp               )) deallocate(csite%today_A_decomp               )
    if(associated(csite%today_Af_decomp              )) deallocate(csite%today_Af_decomp              )
    if(associated(csite%dmean_A_decomp               )) deallocate(csite%dmean_A_decomp               )
    if(associated(csite%dmean_Af_decomp              )) deallocate(csite%dmean_Af_decomp              )
    if(associated(csite%mmean_A_decomp               )) deallocate(csite%mmean_A_decomp               )
    if(associated(csite%mmean_Af_decomp              )) deallocate(csite%mmean_Af_decomp              )
    if(associated(csite%repro                        )) deallocate(csite%repro                        )
    if(associated(csite%veg_rough                    )) deallocate(csite%veg_rough                    )
    if(associated(csite%veg_height                   )) deallocate(csite%veg_height                   )
    if(associated(csite%fsc_in                       )) deallocate(csite%fsc_in                       )
    if(associated(csite%ssc_in                       )) deallocate(csite%ssc_in                       )
    if(associated(csite%ssl_in                       )) deallocate(csite%ssl_in                       )
    if(associated(csite%fsn_in                       )) deallocate(csite%fsn_in                       )
    if(associated(csite%total_plant_nitrogen_uptake  )) deallocate(csite%total_plant_nitrogen_uptake  )
    if(associated(csite%mineralized_N_input )) deallocate(csite%mineralized_N_input  )
    if(associated(csite%mineralized_N_loss  )) deallocate(csite%mineralized_N_loss  )
    if(associated(csite%rshort_g                     )) deallocate(csite%rshort_g                     )
    if(associated(csite%rshort_g_beam                )) deallocate(csite%rshort_g_beam                )
    if(associated(csite%rshort_g_diffuse             )) deallocate(csite%rshort_g_diffuse             )
    if(associated(csite%rlong_g                      )) deallocate(csite%rlong_g                      )
    if(associated(csite%rlong_g_surf                 )) deallocate(csite%rlong_g_surf                 )
    if(associated(csite%rlong_g_incid                )) deallocate(csite%rlong_g_incid                )
    if(associated(csite%rlong_s                      )) deallocate(csite%rlong_s                      )
    if(associated(csite%rlong_s_surf                 )) deallocate(csite%rlong_s_surf                 )
    if(associated(csite%rlong_s_incid                )) deallocate(csite%rlong_s_incid                )
    if(associated(csite%albedt                       )) deallocate(csite%albedt                       )
    if(associated(csite%albedo_beam                  )) deallocate(csite%albedo_beam                  )
    if(associated(csite%albedo_diffuse               )) deallocate(csite%albedo_diffuse               )
    if(associated(csite%rlongup                      )) deallocate(csite%rlongup                      )
    if(associated(csite%rlong_albedo                 )) deallocate(csite%rlong_albedo                 )
    if(associated(csite%total_snow_depth             )) deallocate(csite%total_snow_depth             )
    if(associated(csite%snowfac                      )) deallocate(csite%snowfac                      )
    if(associated(csite%A_decomp                     )) deallocate(csite%A_decomp                     )
    if(associated(csite%f_decomp                     )) deallocate(csite%f_decomp                     )
    if(associated(csite%rh                           )) deallocate(csite%rh                           )
    if(associated(csite%cwd_rh                       )) deallocate(csite%cwd_rh                       )
    if(associated(csite%fuse_flag                    )) deallocate(csite%fuse_flag                    )
    if(associated(csite%pft_density_profile          )) deallocate(csite%pft_density_profile          )
    if(associated(csite%plant_ag_biomass             )) deallocate(csite%plant_ag_biomass             )

    if(associated(csite%mean_wflux                   )) deallocate(csite%mean_wflux                   )
    if(associated(csite%mean_latflux                 )) deallocate(csite%mean_latflux                 )
    if(associated(csite%mean_hflux                   )) deallocate(csite%mean_hflux                   )
    if(associated(csite%mean_runoff                  )) deallocate(csite%mean_runoff                  )
    if(associated(csite%mean_qrunoff                 )) deallocate(csite%mean_qrunoff                 )

    if(associated(csite%htry                         )) deallocate(csite%htry                         )
    if(associated(csite%avg_rk4step                  )) deallocate(csite%avg_rk4step                  )
    if(associated(csite%dmean_rk4step                )) deallocate(csite%dmean_rk4step                )
    if(associated(csite%mmean_rk4step                )) deallocate(csite%mmean_rk4step                )

    if(associated(csite%avg_available_water          )) deallocate(csite%avg_available_water          )

    if(associated(csite%ustar                        )) deallocate(csite%ustar                        )
    if(associated(csite%tstar                        )) deallocate(csite%tstar                        )
    if(associated(csite%qstar                        )) deallocate(csite%qstar                        )
    if(associated(csite%cstar                        )) deallocate(csite%cstar                        )

    if(associated(csite%upwp                         )) deallocate(csite%upwp                         )
    if(associated(csite%qpwp                         )) deallocate(csite%qpwp                         )
    if(associated(csite%cpwp                         )) deallocate(csite%cpwp                         )
    if(associated(csite%tpwp                         )) deallocate(csite%tpwp                         )
    if(associated(csite%wpwp                         )) deallocate(csite%wpwp                         )

    if(associated(csite%avg_carbon_ac                )) deallocate(csite%avg_carbon_ac                )

    if(associated(csite%avg_vapor_vc                 )) deallocate(csite%avg_vapor_vc                 )
    if(associated(csite%avg_dew_cg                   )) deallocate(csite%avg_dew_cg                   )
    if(associated(csite%avg_vapor_gc                 )) deallocate(csite%avg_vapor_gc                 )
    if(associated(csite%avg_wshed_vg                 )) deallocate(csite%avg_wshed_vg                 )
    if(associated(csite%avg_intercepted              )) deallocate(csite%avg_intercepted              )
    if(associated(csite%avg_vapor_ac                 )) deallocate(csite%avg_vapor_ac                 )
    if(associated(csite%avg_transp                   )) deallocate(csite%avg_transp                   )
    if(associated(csite%avg_evap                     )) deallocate(csite%avg_evap                     )
    if(associated(csite%avg_netrad                   )) deallocate(csite%avg_netrad                   )
    if(associated(csite%avg_smoist_gg                )) deallocate(csite%avg_smoist_gg                )
    if(associated(csite%avg_smoist_gc                )) deallocate(csite%avg_smoist_gc                )
    if(associated(csite%avg_runoff                   )) deallocate(csite%avg_runoff                   )
    if(associated(csite%avg_drainage                 )) deallocate(csite%avg_drainage                 )
    if(associated(csite%avg_drainage_heat            )) deallocate(csite%avg_drainage_heat            )
    if(associated(csite%aux                          )) deallocate(csite%aux                          )
    if(associated(csite%aux_s                        )) deallocate(csite%aux_s                        )
    if(associated(csite%avg_sensible_vc              )) deallocate(csite%avg_sensible_vc              )
    if(associated(csite%avg_qwshed_vg                )) deallocate(csite%avg_qwshed_vg                )
    if(associated(csite%avg_qintercepted             )) deallocate(csite%avg_qintercepted             )
    if(associated(csite%avg_sensible_gc              )) deallocate(csite%avg_sensible_gc              )
    if(associated(csite%avg_sensible_ac              )) deallocate(csite%avg_sensible_ac              )
    if(associated(csite%avg_sensible_gg              )) deallocate(csite%avg_sensible_gg              )
    if(associated(csite%avg_runoff_heat              )) deallocate(csite%avg_runoff_heat              )
    if(associated(csite%avg_veg_energy               )) deallocate(csite%avg_veg_energy               )
    if(associated(csite%avg_veg_temp                 )) deallocate(csite%avg_veg_temp                 )
    if(associated(csite%avg_veg_fliq                 )) deallocate(csite%avg_veg_fliq                 )
    if(associated(csite%avg_veg_water                )) deallocate(csite%avg_veg_water                )

    if(associated(csite%watertable                   )) deallocate(csite%watertable                   )
    if(associated(csite%moist_dz                     )) deallocate(csite%moist_dz                     )
    if(associated(csite%ksat                         )) deallocate(csite%ksat                         )
    if(associated(csite%soil_sat_energy              )) deallocate(csite%soil_sat_energy              )
    if(associated(csite%soil_sat_water               )) deallocate(csite%soil_sat_water               )
    if(associated(csite%soil_sat_heat                )) deallocate(csite%soil_sat_heat                )
    if(associated(csite%runoff_A                     )) deallocate(csite%runoff_A                     )
    if(associated(csite%runoff_rate                  )) deallocate(csite%runoff_rate                  )
    if(associated(csite%runoff                       )) deallocate(csite%runoff                       )

    if(associated(csite%dmean_co2_residual        )) deallocate(csite%dmean_co2_residual        )
    if(associated(csite%dmean_energy_residual     )) deallocate(csite%dmean_energy_residual     )
    if(associated(csite%dmean_water_residual      )) deallocate(csite%dmean_water_residual      )
    if(associated(csite%mmean_co2_residual        )) deallocate(csite%mmean_co2_residual        )
    if(associated(csite%mmean_energy_residual     )) deallocate(csite%mmean_energy_residual     )
    if(associated(csite%mmean_water_residual      )) deallocate(csite%mmean_water_residual      )

    do ipa=1,csite%npatches
       if (csite%patch(ipa)%ncohorts > 0) call deallocate_patchtype(csite%patch(ipa))
    enddo

    if(associated(csite%patch                        )) deallocate(csite%patch                        )

    return
  end subroutine deallocate_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_patchtype(cpatch)

    implicit none

    type(patchtype),target :: cpatch
    
    if (cpatch%ncohorts == 0) return
    
    if(associated(cpatch%pft))                 deallocate(cpatch%pft)
    if(associated(cpatch%nplant))              deallocate(cpatch%nplant)
    if(associated(cpatch%hite))                deallocate(cpatch%hite)
    if(associated(cpatch%agb))                 deallocate(cpatch%agb)
    if(associated(cpatch%basarea))             deallocate(cpatch%basarea)
    if(associated(cpatch%dagb_dt))             deallocate(cpatch%dagb_dt)
    if(associated(cpatch%dba_dt))              deallocate(cpatch%dba_dt)
    if(associated(cpatch%ddbh_dt))             deallocate(cpatch%ddbh_dt)
    if(associated(cpatch%dbh))                 deallocate(cpatch%dbh)
    if(associated(cpatch%bdead))               deallocate(cpatch%bdead)
    if(associated(cpatch%bleaf))               deallocate(cpatch%bleaf)
    if(associated(cpatch%phenology_status))    deallocate(cpatch%phenology_status)
    if(associated(cpatch%balive))              deallocate(cpatch%balive)
    if(associated(cpatch%broot))               deallocate(cpatch%broot)
    if(associated(cpatch%bsapwood))            deallocate(cpatch%bsapwood)
    if(associated(cpatch%lai))                 deallocate(cpatch%lai)
    if(associated(cpatch%wpa))                 deallocate(cpatch%wpa)
    if(associated(cpatch%wai))                 deallocate(cpatch%wai)
    if(associated(cpatch%solvable))            deallocate(cpatch%solvable)
    if(associated(cpatch%bstorage))            deallocate(cpatch%bstorage)
    if(associated(cpatch%cb))                  deallocate(cpatch%cb)
    if(associated(cpatch%cb_max))              deallocate(cpatch%cb_max)
    if(associated(cpatch%cbr_bar))             deallocate(cpatch%cbr_bar)
    if(associated(cpatch%mmean_cb))            deallocate(cpatch%mmean_cb)
    if(associated(cpatch%veg_energy))          deallocate(cpatch%veg_energy)
    if(associated(cpatch%veg_temp))            deallocate(cpatch%veg_temp)
    if(associated(cpatch%veg_fliq))            deallocate(cpatch%veg_fliq)
    if(associated(cpatch%veg_water))           deallocate(cpatch%veg_water)
    if(associated(cpatch%mean_gpp))            deallocate(cpatch%mean_gpp)
    if(associated(cpatch%mean_leaf_resp))      deallocate(cpatch%mean_leaf_resp)
    if(associated(cpatch%mean_root_resp))      deallocate(cpatch%mean_root_resp)
    if(associated(cpatch%mean_growth_resp))    deallocate(cpatch%mean_growth_resp)
    if(associated(cpatch%mean_storage_resp))   deallocate(cpatch%mean_storage_resp)
    if(associated(cpatch%mean_vleaf_resp))     deallocate(cpatch%mean_vleaf_resp)
    if(associated(cpatch%today_leaf_resp))     deallocate(cpatch%today_leaf_resp)
    if(associated(cpatch%today_root_resp))     deallocate(cpatch%today_root_resp)
    if(associated(cpatch%today_gpp))           deallocate(cpatch%today_gpp)
    if(associated(cpatch%today_gpp_pot))       deallocate(cpatch%today_gpp_pot)
    if(associated(cpatch%today_gpp_max))       deallocate(cpatch%today_gpp_max)
    if(associated(cpatch%growth_respiration))  deallocate(cpatch%growth_respiration)
    if(associated(cpatch%storage_respiration)) deallocate(cpatch%storage_respiration)
    if(associated(cpatch%vleaf_respiration))   deallocate(cpatch%vleaf_respiration)
    if(associated(cpatch%dmean_gpp         ))  deallocate(cpatch%dmean_gpp         )
    if(associated(cpatch%dmean_leaf_resp   ))  deallocate(cpatch%dmean_leaf_resp   )
    if(associated(cpatch%dmean_root_resp   ))  deallocate(cpatch%dmean_root_resp   )
    if(associated(cpatch%mmean_gpp         ))  deallocate(cpatch%mmean_gpp         )
    if(associated(cpatch%mmean_leaf_resp   ))  deallocate(cpatch%mmean_leaf_resp   )
    if(associated(cpatch%mmean_root_resp   ))  deallocate(cpatch%mmean_root_resp   )
    if(associated(cpatch%mmean_growth_resp ))  deallocate(cpatch%mmean_growth_resp )
    if(associated(cpatch%mmean_storage_resp))  deallocate(cpatch%mmean_storage_resp)
    if(associated(cpatch%mmean_vleaf_resp  ))  deallocate(cpatch%mmean_vleaf_resp  )
    if(associated(cpatch%fsn))                 deallocate(cpatch%fsn)
    if(associated(cpatch%monthly_dndt))        deallocate(cpatch%monthly_dndt)
    if(associated(cpatch%mort_rate))           deallocate(cpatch%mort_rate)
    if(associated(cpatch%mmean_mort_rate))     deallocate(cpatch%mmean_mort_rate)

    if(associated(cpatch%old_stoma_data))         deallocate(cpatch%old_stoma_data)
    if(associated(cpatch%old_stoma_vector))       deallocate(cpatch%old_stoma_vector)
    if(associated(cpatch%Psi_open))               deallocate(cpatch%Psi_open)
    if(associated(cpatch%krdepth))                deallocate(cpatch%krdepth)
    if(associated(cpatch%first_census))           deallocate(cpatch%first_census)
    if(associated(cpatch%new_recruit_flag))       deallocate(cpatch%new_recruit_flag)
    if(associated(cpatch%light_level))            deallocate(cpatch%light_level)
    if(associated(cpatch%dmean_light_level))      deallocate(cpatch%dmean_light_level)
    if(associated(cpatch%mmean_light_level))      deallocate(cpatch%mmean_light_level)
    if(associated(cpatch%light_level_beam))       deallocate(cpatch%light_level_beam)
    if(associated(cpatch%dmean_light_level_beam)) deallocate(cpatch%dmean_light_level_beam)
    if(associated(cpatch%mmean_light_level_beam)) deallocate(cpatch%mmean_light_level_beam)
    if(associated(cpatch%light_level_diff))       deallocate(cpatch%light_level_diff)
    if(associated(cpatch%dmean_light_level_diff)) deallocate(cpatch%dmean_light_level_diff)
    if(associated(cpatch%mmean_light_level_diff)) deallocate(cpatch%mmean_light_level_diff)
    if(associated(cpatch%beamext_level))          deallocate(cpatch%beamext_level)
    if(associated(cpatch%dmean_beamext_level))    deallocate(cpatch%dmean_beamext_level)
    if(associated(cpatch%mmean_beamext_level))    deallocate(cpatch%mmean_beamext_level)
    if(associated(cpatch%diffext_level))          deallocate(cpatch%diffext_level)
    if(associated(cpatch%dmean_diffext_level))    deallocate(cpatch%dmean_diffext_level)
    if(associated(cpatch%mmean_diffext_level))    deallocate(cpatch%mmean_diffext_level)
    if(associated(cpatch%lambda_light))           deallocate(cpatch%lambda_light)
    if(associated(cpatch%dmean_lambda_light))     deallocate(cpatch%dmean_lambda_light)
    if(associated(cpatch%mmean_lambda_light))     deallocate(cpatch%mmean_lambda_light)
    if(associated(cpatch%norm_par_beam))          deallocate(cpatch%norm_par_beam)
    if(associated(cpatch%dmean_norm_par_beam))    deallocate(cpatch%dmean_norm_par_beam)
    if(associated(cpatch%mmean_norm_par_beam))    deallocate(cpatch%mmean_norm_par_beam)
    if(associated(cpatch%norm_par_diff))          deallocate(cpatch%norm_par_diff)
    if(associated(cpatch%dmean_norm_par_diff))    deallocate(cpatch%dmean_norm_par_diff)
    if(associated(cpatch%mmean_norm_par_diff))    deallocate(cpatch%mmean_norm_par_diff)
    if(associated(cpatch%par_v))                  deallocate(cpatch%par_v)
    if(associated(cpatch%par_v_beam))             deallocate(cpatch%par_v_beam)
    if(associated(cpatch%par_v_diffuse))          deallocate(cpatch%par_v_diffuse)
    if(associated(cpatch%dmean_par_v))            deallocate(cpatch%dmean_par_v)
    if(associated(cpatch%dmean_par_v_beam))       deallocate(cpatch%dmean_par_v_beam)
    if(associated(cpatch%dmean_par_v_diff))       deallocate(cpatch%dmean_par_v_diff)
    if(associated(cpatch%mmean_par_v))            deallocate(cpatch%mmean_par_v)
    if(associated(cpatch%mmean_par_v_beam))       deallocate(cpatch%mmean_par_v_beam)
    if(associated(cpatch%mmean_par_v_diff))       deallocate(cpatch%mmean_par_v_diff)
    if(associated(cpatch%rshort_v))               deallocate(cpatch%rshort_v)
    if(associated(cpatch%rshort_v_beam))          deallocate(cpatch%rshort_v_beam)
    if(associated(cpatch%rshort_v_diffuse))       deallocate(cpatch%rshort_v_diffuse)
    if(associated(cpatch%rlong_v))                deallocate(cpatch%rlong_v)
    if(associated(cpatch%rlong_v_surf))           deallocate(cpatch%rlong_v_surf)
    if(associated(cpatch%rlong_v_incid))          deallocate(cpatch%rlong_v_incid)
    if(associated(cpatch%rb))                     deallocate(cpatch%rb)
    if(associated(cpatch%A_open))                 deallocate(cpatch%A_open)    
    if(associated(cpatch%A_closed))               deallocate(cpatch%A_closed)
    if(associated(cpatch%Psi_closed))             deallocate(cpatch%Psi_closed)
    if(associated(cpatch%rsw_open))               deallocate(cpatch%rsw_open)
    if(associated(cpatch%rsw_closed))             deallocate(cpatch%rsw_closed)
    if(associated(cpatch%fsw))                    deallocate(cpatch%fsw)
    if(associated(cpatch%fs_open))                deallocate(cpatch%fs_open)
    if(associated(cpatch%dmean_fs_open))          deallocate(cpatch%dmean_fs_open)
    if(associated(cpatch%dmean_fsw))              deallocate(cpatch%dmean_fsw)
    if(associated(cpatch%dmean_fsn))              deallocate(cpatch%dmean_fsn)
    if(associated(cpatch%mmean_fs_open))          deallocate(cpatch%mmean_fs_open)
    if(associated(cpatch%mmean_fsw))              deallocate(cpatch%mmean_fsw)
    if(associated(cpatch%mmean_fsn))              deallocate(cpatch%mmean_fsn)
    if(associated(cpatch%stomatal_resistance))    deallocate(cpatch%stomatal_resistance)
    if(associated(cpatch%leaf_maintenance))       deallocate(cpatch%leaf_maintenance)
    if(associated(cpatch%root_maintenance))       deallocate(cpatch%root_maintenance)
    if(associated(cpatch%mmean_leaf_maintenance)) deallocate(cpatch%mmean_leaf_maintenance)
    if(associated(cpatch%mmean_root_maintenance)) deallocate(cpatch%mmean_root_maintenance)
    if(associated(cpatch%leaf_drop))              deallocate(cpatch%leaf_drop)
    if(associated(cpatch%mmean_leaf_drop))        deallocate(cpatch%mmean_leaf_drop)
    if(associated(cpatch%bseeds))                 deallocate(cpatch%bseeds)
    if(associated(cpatch%leaf_respiration))       deallocate(cpatch%leaf_respiration)
    if(associated(cpatch%root_respiration))       deallocate(cpatch%root_respiration)
    if(associated(cpatch%hcapveg))                deallocate(cpatch%hcapveg)
    if(associated(cpatch%gpp))                    deallocate(cpatch%gpp)
    if(associated(cpatch%paw_avg))                deallocate(cpatch%paw_avg)
    if(associated(cpatch%turnover_amp))           deallocate(cpatch%turnover_amp)
    if(associated(cpatch%llspan))                 deallocate(cpatch%llspan)
    if(associated(cpatch%vm_bar))                 deallocate(cpatch%vm_bar)
    if(associated(cpatch%sla))                    deallocate(cpatch%sla)

    return
  end subroutine deallocate_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine huge_edtype(cgrid)
  
  ! =======================================
  ! zero state variables
  ! =======================================
    
    implicit none
    
    type(edtype),target :: cgrid
    
    if(associated(cgrid%lat                     )) cgrid%lat                      = large_real
    if(associated(cgrid%lon                     )) cgrid%lon                      = large_real
    if(associated(cgrid%natm                    )) cgrid%natm                     = large_integer  ! Integer
    if(associated(cgrid%xatm                    )) cgrid%xatm                     = large_integer  ! Integer
    if(associated(cgrid%yatm                    )) cgrid%yatm                     = large_integer  ! Integer
    if(associated(cgrid%ntext_soil              )) cgrid%ntext_soil               = large_integer  ! Integer
    if(associated(cgrid%lsl                     )) cgrid%lsl                      = large_integer  ! Integer
    
    if(associated(cgrid%pysi_n                  )) cgrid%pysi_n                   = large_integer  ! Integer
    if(associated(cgrid%pysi_id                 )) cgrid%pysi_id                  = large_integer  ! Integer
    if(associated(cgrid%walltime_py             )) cgrid%walltime_py              = large_double
    if(associated(cgrid%sensflux_py             )) cgrid%sensflux_py              = large_real
    if(associated(cgrid%site_adjacency          )) cgrid%site_adjacency           = large_real
    if(associated(cgrid%wbar                    )) cgrid%wbar                     = large_real
    if(associated(cgrid%Te                      )) cgrid%Te                       = large_real
    if(associated(cgrid%zbar                    )) cgrid%zbar                     = large_real
!!  if(associated(cgrid%tau                     )) cgrid%tau                      = large_real
    if(associated(cgrid%sheat                   )) cgrid%sheat                    = large_real
    if(associated(cgrid%baseflow                )) cgrid%baseflow                 = large_real
    if(associated(cgrid%runoff                  )) cgrid%runoff                   = large_real
    if(associated(cgrid%swliq                   )) cgrid%swliq                    = large_real
    
    if(associated(cgrid%ilon                    )) cgrid%ilon                     = large_integer  ! Integer
    if(associated(cgrid%ilat                    )) cgrid%ilat                     = large_integer  ! Integer
    if(associated(cgrid%total_agb               )) cgrid%total_agb                = large_real
    if(associated(cgrid%total_basal_area        )) cgrid%total_basal_area         = large_real
    if(associated(cgrid%total_agb_growth        )) cgrid%total_agb_growth         = large_real
    if(associated(cgrid%total_agb_mort          )) cgrid%total_agb_mort           = large_real
    if(associated(cgrid%total_agb_recruit       )) cgrid%total_agb_recruit        = large_real
    if(associated(cgrid%total_basal_area_growth )) cgrid%total_basal_area_growth  = large_real
    if(associated(cgrid%total_basal_area_mort   )) cgrid%total_basal_area_mort    = large_real
    if(associated(cgrid%total_basal_area_recruit)) cgrid%total_basal_area_recruit = large_real
    if(associated(cgrid%nsites                  )) cgrid%nsites                   = large_integer  ! Integer
    if(associated(cgrid%sitenums                )) cgrid%sitenums                 = large_integer  ! Integer
    if(associated(cgrid%load_adjacency          )) cgrid%load_adjacency           = large_integer  ! Integer
    if(associated(cgrid%cosz                    )) cgrid%cosz                     = large_real
    if(associated(cgrid%mean_gpp                )) cgrid%mean_gpp                 = large_real
    if(associated(cgrid%mean_precip             )) cgrid%mean_precip              = large_real
    if(associated(cgrid%mean_qprecip            )) cgrid%mean_qprecip             = large_real
    if(associated(cgrid%mean_netrad             )) cgrid%mean_netrad              = large_real
    if(associated(cgrid%cbudget_initialstorage  )) cgrid%cbudget_initialstorage   = large_real
    if(associated(cgrid%cbudget_nep             )) cgrid%cbudget_nep              = large_real
    if(associated(cgrid%nbudget_initialstorage  )) cgrid%nbudget_initialstorage   = large_real
    if(associated(cgrid%basal_area              )) cgrid%basal_area               = large_real
    if(associated(cgrid%agb                     )) cgrid%agb                      = large_real
    if(associated(cgrid%pldens                  )) cgrid%pldens                   = large_real
    if(associated(cgrid%bseeds                  )) cgrid%bseeds                   = large_real
    
!   metinput has itself many pointers, so nothing is really allocated 
!    if(associated(cgrid%metinput                )) cgrid%metinput                 = large_real

    if(associated(cgrid%met                     )) then
       cgrid%met(:)%nir_beam       = large_real
       cgrid%met(:)%nir_diffuse    = large_real
       cgrid%met(:)%par_beam       = large_real
       cgrid%met(:)%par_diffuse    = large_real
       cgrid%met(:)%atm_tmp        = large_real
       cgrid%met(:)%atm_theta      = large_real
       cgrid%met(:)%atm_enthalpy   = large_real
       cgrid%met(:)%atm_shv        = large_real
       cgrid%met(:)%theta          = large_real
       cgrid%met(:)%rshort         = large_real
       cgrid%met(:)%rshort_diffuse = large_real
       cgrid%met(:)%rlong          = large_real
       cgrid%met(:)%pcpg           = large_real
       cgrid%met(:)%qpcpg          = large_real
       cgrid%met(:)%dpcpg          = large_real
       cgrid%met(:)%vels           = large_real
       cgrid%met(:)%vels_stab      = large_real
       cgrid%met(:)%vels_unstab    = large_real
       cgrid%met(:)%prss           = large_real
       cgrid%met(:)%exner          = large_real
       cgrid%met(:)%geoht          = large_real
       cgrid%met(:)%atm_co2        = large_real
    end if                         

    if(associated(cgrid%lapse                   )) then
       cgrid%lapse(:)%nir_beam       = large_real
       cgrid%lapse(:)%nir_diffuse    = large_real
       cgrid%lapse(:)%par_beam       = large_real
       cgrid%lapse(:)%par_diffuse    = large_real
       cgrid%lapse(:)%atm_tmp        = large_real
       cgrid%lapse(:)%atm_theta      = large_real
       cgrid%lapse(:)%atm_enthalpy   = large_real
       cgrid%lapse(:)%atm_shv        = large_real
       cgrid%lapse(:)%theta          = large_real
       cgrid%lapse(:)%rshort         = large_real
       cgrid%lapse(:)%rshort_diffuse = large_real
       cgrid%lapse(:)%rlong          = large_real
       cgrid%lapse(:)%pcpg           = large_real
       cgrid%lapse(:)%qpcpg          = large_real
       cgrid%lapse(:)%dpcpg          = large_real
       cgrid%lapse(:)%vels           = large_real
       cgrid%lapse(:)%vels_stab      = large_real
       cgrid%lapse(:)%vels_unstab    = large_real
       cgrid%lapse(:)%prss           = large_real
       cgrid%lapse(:)%exner          = large_real
       cgrid%lapse(:)%geoht          = large_real
       cgrid%lapse(:)%atm_co2        = large_real
    end if                         

    if(associated(cgrid%lai                     )) cgrid%lai                      = large_real
    if(associated(cgrid%avg_lma                 )) cgrid%avg_lma                  = large_real
    if(associated(cgrid%wpa                     )) cgrid%wpa                      = large_real
    if(associated(cgrid%wai                     )) cgrid%wai                      = large_real
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    if(associated(cgrid%avg_vapor_vc            )) cgrid%avg_vapor_vc             = large_real
    if(associated(cgrid%avg_dew_cg              )) cgrid%avg_dew_cg               = large_real
    if(associated(cgrid%avg_vapor_gc            )) cgrid%avg_vapor_gc             = large_real
    if(associated(cgrid%avg_wshed_vg            )) cgrid%avg_wshed_vg             = large_real
    if(associated(cgrid%avg_intercepted         )) cgrid%avg_intercepted          = large_real
    if(associated(cgrid%avg_vapor_ac            )) cgrid%avg_vapor_ac             = large_real
    if(associated(cgrid%avg_transp              )) cgrid%avg_transp               = large_real
    if(associated(cgrid%avg_evap                )) cgrid%avg_evap                 = large_real
    if(associated(cgrid%avg_smoist_gg           )) cgrid%avg_smoist_gg            = large_real
    if(associated(cgrid%avg_smoist_gc           )) cgrid%avg_smoist_gc            = large_real
    if(associated(cgrid%avg_runoff              )) cgrid%avg_runoff               = large_real
    if(associated(cgrid%avg_drainage            )) cgrid%avg_drainage             = large_real
    if(associated(cgrid%avg_drainage_heat       )) cgrid%avg_drainage_heat        = large_real
    if(associated(cgrid%aux                     )) cgrid%aux                      = large_real
    if(associated(cgrid%aux_s                   )) cgrid%aux_s                    = large_real
    if(associated(cgrid%avg_sensible_vc         )) cgrid%avg_sensible_vc          = large_real
    if(associated(cgrid%avg_qwshed_vg           )) cgrid%avg_qwshed_vg            = large_real
    if(associated(cgrid%avg_qintercepted        )) cgrid%avg_qintercepted         = large_real
    if(associated(cgrid%avg_sensible_gc         )) cgrid%avg_sensible_gc          = large_real
    if(associated(cgrid%avg_sensible_ac         )) cgrid%avg_sensible_ac          = large_real
    if(associated(cgrid%avg_sensible_gg         )) cgrid%avg_sensible_gg          = large_real
    if(associated(cgrid%avg_runoff_heat         )) cgrid%avg_runoff_heat          = large_real

    if(associated(cgrid%max_veg_temp          )) cgrid%max_veg_temp           = large_real
    if(associated(cgrid%min_veg_temp          )) cgrid%min_veg_temp           = large_real
    if(associated(cgrid%max_soil_temp          )) cgrid%max_soil_temp           = large_real
    if(associated(cgrid%min_soil_temp          )) cgrid%min_soil_temp           = large_real

    ! Fast time state diagnostics
    if(associated(cgrid%avg_veg_energy          )) cgrid%avg_veg_energy           = large_real
    if(associated(cgrid%avg_veg_hcap            )) cgrid%avg_veg_hcap             = large_real
    if(associated(cgrid%avg_veg_temp            )) cgrid%avg_veg_temp             = large_real
    if(associated(cgrid%avg_veg_fliq            )) cgrid%avg_veg_fliq             = large_real
    if(associated(cgrid%avg_veg_water           )) cgrid%avg_veg_water            = large_real
    if(associated(cgrid%avg_can_temp            )) cgrid%avg_can_temp             = large_real
    if(associated(cgrid%avg_can_shv             )) cgrid%avg_can_shv              = large_real
    if(associated(cgrid%avg_can_co2             )) cgrid%avg_can_co2              = large_real
    if(associated(cgrid%avg_can_rhos            )) cgrid%avg_can_rhos             = large_real
    if(associated(cgrid%avg_can_prss            )) cgrid%avg_can_prss             = large_real
    if(associated(cgrid%avg_can_theta           )) cgrid%avg_can_theta            = large_real
    if(associated(cgrid%avg_can_enthalpy        )) cgrid%avg_can_enthalpy         = large_real
    if(associated(cgrid%avg_can_depth           )) cgrid%avg_can_depth            = large_real
    if(associated(cgrid%avg_soil_energy         )) cgrid%avg_soil_energy          = large_real
    if(associated(cgrid%avg_soil_water          )) cgrid%avg_soil_water           = large_real
    if(associated(cgrid%avg_soil_temp           )) cgrid%avg_soil_temp            = large_real
    if(associated(cgrid%avg_soil_fracliq        )) cgrid%avg_soil_fracliq         = large_real
    if(associated(cgrid%avg_soil_wetness        )) cgrid%avg_soil_wetness         = large_real
    if(associated(cgrid%avg_skin_temp           )) cgrid%avg_skin_temp            = large_real
    if(associated(cgrid%avg_available_water     )) cgrid%avg_available_water      = large_real
  
    if(associated(cgrid%avg_lai_ebalvars        )) cgrid%avg_lai_ebalvars         = large_real

    !!! added for NACP intercomparison (MCD)
    if(associated(cgrid%avg_sfcw_depth          )) cgrid%avg_sfcw_depth           = large_real
    if(associated(cgrid%avg_sfcw_energy         )) cgrid%avg_sfcw_energy          = large_real
    if(associated(cgrid%avg_sfcw_mass           )) cgrid%avg_sfcw_mass            = large_real
    if(associated(cgrid%avg_sfcw_tempk          )) cgrid%avg_sfcw_tempk           = large_real
    if(associated(cgrid%avg_sfcw_fracliq        )) cgrid%avg_sfcw_fracliq         = large_real
    if(associated(cgrid%avg_bdead               )) cgrid%avg_bdead                = large_real
    if(associated(cgrid%avg_balive              )) cgrid%avg_balive               = large_real
    if(associated(cgrid%avg_bleaf               )) cgrid%avg_bleaf                = large_real
    if(associated(cgrid%avg_broot               )) cgrid%avg_broot                = large_real
    if(associated(cgrid%avg_bsapwood            )) cgrid%avg_bsapwood             = large_real
    if(associated(cgrid%avg_bstorage            )) cgrid%avg_bstorage             = large_real
    if(associated(cgrid%avg_bseeds              )) cgrid%avg_bseeds               = large_real
    if(associated(cgrid%avg_fsc                 )) cgrid%avg_fsc                  = large_real
    if(associated(cgrid%avg_ssc                 )) cgrid%avg_ssc                  = large_real
    if(associated(cgrid%avg_stsc                )) cgrid%avg_stsc                 = large_real
  if(associated(cgrid%avg_fsn                 )) cgrid%avg_fsn                  = large_real
       if(associated(cgrid%avg_msn                 )) cgrid%avg_msn                  = large_real

     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf  )) cgrid%Cleaf   = large_real
     if(associated(cgrid%Croot  )) cgrid%Croot   = large_real
     if(associated(cgrid%Cstore )) cgrid%Cstore  = large_real 
     if(associated(cgrid%Ccwd   )) cgrid%Ccwd    = large_real
     if(associated(cgrid%Nleaf  )) cgrid%Nleaf   = large_real
     if(associated(cgrid%Ndead  )) cgrid%Ndead   = large_real
     if(associated(cgrid%Nroot  )) cgrid%Nroot   = large_real
     if(associated(cgrid%Nstore )) cgrid%Nstore  = large_real
     if(associated(cgrid%Ncwd   )) cgrid%Ncwd    = large_real
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf_grow       )) cgrid%Cleaf_grow   = large_real
     if(associated(cgrid%Croot_grow       )) cgrid%Croot_grow   = large_real
     if(associated(cgrid%Cdead_grow       )) cgrid%Cdead_grow   = large_real
     if(associated(cgrid%Cstore_grow      )) cgrid%Cstore_grow  = large_real
     if(associated(cgrid%Cleaf_litter_flux)) cgrid%Cleaf_litter_flux = large_real
     if(associated(cgrid%Croot_litter_flux)) cgrid%Croot_litter_flux = large_real
     if(associated(cgrid%Ccwd_flux        )) cgrid%Ccwd_flux    = large_real
     if(associated(cgrid%Nleaf_grow       )) cgrid%Nleaf_grow   = large_real
     if(associated(cgrid%Ndead_grow       )) cgrid%Ndead_grow   = large_real
     if(associated(cgrid%Nroot_grow       )) cgrid%Nroot_grow   = large_real
     if(associated(cgrid%Nstore_grow      )) cgrid%Nstore_grow  = large_real
     if(associated(cgrid%Nleaf_litter_flux)) cgrid%Nleaf_litter_flux = large_real
     if(associated(cgrid%Nroot_litter_flux)) cgrid%Nroot_litter_flux = large_real 
     if(associated(cgrid%Ncwd_flux        )) cgrid%Ncwd_flux    = large_real
     if(associated(cgrid%Nbiomass_uptake  )) cgrid%Nbiomass_uptake = large_real
     if(associated(cgrid%Ngross_min       )) cgrid%Ngross_min   = large_real
     if(associated(cgrid%Nnet_min         )) cgrid%Nnet_min     = large_real

    ! ---------------------------------------------

    if(associated(cgrid%avg_nir_beam            )) cgrid%avg_nir_beam             = large_real
    if(associated(cgrid%avg_nir_diffuse         )) cgrid%avg_nir_diffuse          = large_real
    if(associated(cgrid%avg_par_beam            )) cgrid%avg_par_beam             = large_real
    if(associated(cgrid%avg_par_diffuse         )) cgrid%avg_par_diffuse          = large_real
    if(associated(cgrid%avg_atm_tmp             )) cgrid%avg_atm_tmp              = large_real
    if(associated(cgrid%avg_atm_shv             )) cgrid%avg_atm_shv              = large_real
    if(associated(cgrid%avg_rshort              )) cgrid%avg_rshort               = large_real
    if(associated(cgrid%avg_rshort_diffuse      )) cgrid%avg_rshort_diffuse       = large_real
    if(associated(cgrid%avg_rlong               )) cgrid%avg_rlong                = large_real
    if(associated(cgrid%avg_pcpg                )) cgrid%avg_pcpg                 = large_real
    if(associated(cgrid%avg_qpcpg               )) cgrid%avg_qpcpg                = large_real
    if(associated(cgrid%avg_dpcpg               )) cgrid%avg_dpcpg                = large_real
    if(associated(cgrid%avg_vels                )) cgrid%avg_vels                 = large_real
    if(associated(cgrid%avg_atm_prss            )) cgrid%avg_atm_prss             = large_real
    if(associated(cgrid%avg_exner               )) cgrid%avg_exner                = large_real
    if(associated(cgrid%avg_geoht               )) cgrid%avg_geoht                = large_real
    if(associated(cgrid%avg_atm_co2             )) cgrid%avg_atm_co2              = large_real
    if(associated(cgrid%avg_albedt              )) cgrid%avg_albedt               = large_real
    if(associated(cgrid%avg_rlongup             )) cgrid%avg_rlongup              = large_real

    if(associated(cgrid%runoff                  )) cgrid%runoff                   = large_real
    
    if(associated(cgrid%dmean_gpp               )) cgrid%dmean_gpp                = large_real
    if(associated(cgrid%dmean_evap              )) cgrid%dmean_evap               = large_real
    if(associated(cgrid%dmean_transp            )) cgrid%dmean_transp             = large_real
    if(associated(cgrid%dmean_pcpg              )) cgrid%dmean_pcpg               = large_real
    if(associated(cgrid%dmean_runoff            )) cgrid%dmean_runoff             = large_real
    if(associated(cgrid%dmean_drainage          )) cgrid%dmean_drainage           = large_real
    if(associated(cgrid%dmean_vapor_ac          )) cgrid%dmean_vapor_ac           = large_real
    if(associated(cgrid%dmean_vapor_gc          )) cgrid%dmean_vapor_gc           = large_real
    if(associated(cgrid%dmean_vapor_vc          )) cgrid%dmean_vapor_vc           = large_real
    if(associated(cgrid%dmean_sensible_vc       )) cgrid%dmean_sensible_vc        = large_real
    if(associated(cgrid%dmean_sensible_gc       )) cgrid%dmean_sensible_gc        = large_real
    if(associated(cgrid%dmean_sensible_ac       )) cgrid%dmean_sensible_ac        = large_real
    if(associated(cgrid%dmean_plresp            )) cgrid%dmean_plresp             = large_real
    if(associated(cgrid%dmean_rh                )) cgrid%dmean_rh                 = large_real
    if(associated(cgrid%dmean_leaf_resp         )) cgrid%dmean_leaf_resp          = large_real
    if(associated(cgrid%dmean_root_resp         )) cgrid%dmean_root_resp          = large_real
    if(associated(cgrid%dmean_growth_resp       )) cgrid%dmean_growth_resp        = large_real
    if(associated(cgrid%dmean_storage_resp      )) cgrid%dmean_storage_resp       = large_real
    if(associated(cgrid%dmean_vleaf_resp        )) cgrid%dmean_vleaf_resp         = large_real
    if(associated(cgrid%dmean_nep               )) cgrid%dmean_nep                = large_real
    if(associated(cgrid%dmean_soil_temp         )) cgrid%dmean_soil_temp          = large_real
    if(associated(cgrid%dmean_soil_water        )) cgrid%dmean_soil_water         = large_real
    if(associated(cgrid%dmean_fs_open           )) cgrid%dmean_fs_open            = large_real
    if(associated(cgrid%dmean_fsw               )) cgrid%dmean_fsw                = large_real
    if(associated(cgrid%dmean_fsn               )) cgrid%dmean_fsn                = large_real
    if(associated(cgrid%dmean_gpp_dbh           )) cgrid%dmean_gpp_dbh            = large_real
    if(associated(cgrid%dmean_can_temp          )) cgrid%dmean_can_temp           = large_real
    if(associated(cgrid%dmean_can_shv           )) cgrid%dmean_can_shv            = large_real
    if(associated(cgrid%dmean_can_co2           )) cgrid%dmean_can_co2            = large_real
    if(associated(cgrid%dmean_can_rhos          )) cgrid%dmean_can_rhos           = large_real
    if(associated(cgrid%dmean_can_prss          )) cgrid%dmean_can_prss           = large_real
    if(associated(cgrid%dmean_can_theta         )) cgrid%dmean_can_theta          = large_real
    if(associated(cgrid%dmean_veg_energy        )) cgrid%dmean_veg_energy         = large_real
    if(associated(cgrid%dmean_veg_water         )) cgrid%dmean_veg_water          = large_real
    if(associated(cgrid%dmean_veg_temp          )) cgrid%dmean_veg_temp           = large_real
    if(associated(cgrid%dmean_veg_hcap          )) cgrid%dmean_veg_hcap           = large_real
    if(associated(cgrid%dmean_atm_temp          )) cgrid%dmean_atm_temp           = large_real
    if(associated(cgrid%dmean_atm_shv           )) cgrid%dmean_atm_shv            = large_real
    if(associated(cgrid%dmean_atm_prss          )) cgrid%dmean_atm_prss           = large_real
    if(associated(cgrid%dmean_atm_vels          )) cgrid%dmean_atm_vels           = large_real
    if(associated(cgrid%dmean_rshort            )) cgrid%dmean_rshort             = large_real
    if(associated(cgrid%dmean_rlong             )) cgrid%dmean_rlong              = large_real
    if(associated(cgrid%dmean_co2_residual      )) cgrid%dmean_co2_residual       = large_real
    if(associated(cgrid%dmean_energy_residual   )) cgrid%dmean_energy_residual    = large_real
    if(associated(cgrid%dmean_water_residual    )) cgrid%dmean_water_residual     = large_real
    if(associated(cgrid%lai_pft                 )) cgrid%lai_pft                  = large_real
    if(associated(cgrid%wpa_pft                 )) cgrid%wpa_pft                  = large_real
    if(associated(cgrid%wai_pft                 )) cgrid%wai_pft                  = large_real
    if(associated(cgrid%mmean_gpp               )) cgrid%mmean_gpp                = large_real
    if(associated(cgrid%mmean_evap              )) cgrid%mmean_evap               = large_real
    if(associated(cgrid%mmean_transp            )) cgrid%mmean_transp             = large_real
    if(associated(cgrid%mmean_sensible_vc       )) cgrid%mmean_sensible_vc        = large_real
    if(associated(cgrid%mmean_sensible_gc       )) cgrid%mmean_sensible_gc        = large_real
    if(associated(cgrid%mmean_sensible_ac       )) cgrid%mmean_sensible_ac        = large_real
    if(associated(cgrid%mmean_nep               )) cgrid%mmean_nep                = large_real
    if(associated(cgrid%mmean_soil_temp         )) cgrid%mmean_soil_temp          = large_real
    if(associated(cgrid%mmean_soil_water        )) cgrid%mmean_soil_water         = large_real
    if(associated(cgrid%mmean_fs_open           )) cgrid%mmean_fs_open            = large_real
    if(associated(cgrid%mmean_fsw               )) cgrid%mmean_fsw                = large_real
    if(associated(cgrid%mmean_fsn               )) cgrid%mmean_fsn                = large_real
    if(associated(cgrid%mmean_plresp            )) cgrid%mmean_plresp             = large_real
    if(associated(cgrid%mmean_rh                )) cgrid%mmean_rh                 = large_real
    if(associated(cgrid%mmean_leaf_resp         )) cgrid%mmean_leaf_resp          = large_real
    if(associated(cgrid%mmean_root_resp         )) cgrid%mmean_root_resp          = large_real
    if(associated(cgrid%mmean_growth_resp       )) cgrid%mmean_growth_resp        = large_real
    if(associated(cgrid%mmean_storage_resp      )) cgrid%mmean_storage_resp       = large_real
    if(associated(cgrid%mmean_vleaf_resp        )) cgrid%mmean_vleaf_resp         = large_real
    if(associated(cgrid%mmean_gpp_dbh           )) cgrid%mmean_gpp_dbh            = large_real
    if(associated(cgrid%mmean_lai_pft           )) cgrid%mmean_lai_pft            = large_real
    if(associated(cgrid%mmean_wpa_pft           )) cgrid%mmean_wpa_pft            = large_real
    if(associated(cgrid%mmean_wai_pft           )) cgrid%mmean_wai_pft            = large_real
    if(associated(cgrid%dmean_can_temp          )) cgrid%dmean_can_temp           = large_real
    if(associated(cgrid%dmean_can_shv           )) cgrid%dmean_can_shv            = large_real
    if(associated(cgrid%dmean_can_co2           )) cgrid%dmean_can_co2            = large_real
    if(associated(cgrid%dmean_can_rhos          )) cgrid%dmean_can_rhos           = large_real
    if(associated(cgrid%dmean_can_prss          )) cgrid%dmean_can_prss           = large_real
    if(associated(cgrid%dmean_can_theta         )) cgrid%dmean_can_theta          = large_real
    if(associated(cgrid%dmean_veg_energy        )) cgrid%dmean_veg_energy         = large_real
    if(associated(cgrid%dmean_veg_water         )) cgrid%dmean_veg_water          = large_real
    if(associated(cgrid%dmean_veg_hcap          )) cgrid%dmean_veg_hcap           = large_real
    if(associated(cgrid%dmean_veg_temp          )) cgrid%dmean_veg_temp           = large_real
    if(associated(cgrid%dmean_atm_temp          )) cgrid%dmean_atm_temp           = large_real
    if(associated(cgrid%dmean_atm_shv           )) cgrid%dmean_atm_shv            = large_real
    if(associated(cgrid%dmean_atm_prss          )) cgrid%dmean_atm_prss           = large_real
    if(associated(cgrid%dmean_atm_vels          )) cgrid%dmean_atm_vels           = large_real
    if(associated(cgrid%dmean_rlong             )) cgrid%dmean_rlong              = large_real
    if(associated(cgrid%dmean_rshort            )) cgrid%dmean_rshort             = large_real
    if(associated(cgrid%mmean_co2_residual      )) cgrid%mmean_co2_residual       = large_real
    if(associated(cgrid%mmean_energy_residual   )) cgrid%mmean_energy_residual    = large_real
    if(associated(cgrid%mmean_water_residual    )) cgrid%mmean_water_residual     = large_real
    if(associated(cgrid%bseeds_pft              )) cgrid%bseeds_pft               = large_real
    if(associated(cgrid%agb_pft                 )) cgrid%agb_pft                  = large_real
    if(associated(cgrid%ba_pft                  )) cgrid%ba_pft                   = large_real
    if(associated(cgrid%stdev_gpp               )) cgrid%stdev_gpp                = large_real
    if(associated(cgrid%stdev_evap              )) cgrid%stdev_evap               = large_real
    if(associated(cgrid%stdev_transp            )) cgrid%stdev_transp             = large_real
    if(associated(cgrid%stdev_sensible          )) cgrid%stdev_sensible           = large_real
    if(associated(cgrid%stdev_nep               )) cgrid%stdev_nep                = large_real
    if(associated(cgrid%stdev_rh                )) cgrid%stdev_rh                 = large_real
    if(associated(cgrid%disturbance_rates       )) cgrid%disturbance_rates        = large_real
    if(associated(cgrid%workload                )) cgrid%workload                 = large_real


    return
  end subroutine huge_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine huge_polygontype(cpoly)
    
    implicit none
    type(polygontype), target :: cpoly
    
    if(associated(cpoly%sipa_id                     )) cpoly%sipa_id                     = large_integer ! Integer
    if(associated(cpoly%sipa_n                      )) cpoly%sipa_n                      = large_integer ! Integer
    if(associated(cpoly%patch_count                 )) cpoly%patch_count                 = large_integer ! Integer
    if(associated(cpoly%sitenum                     )) cpoly%sitenum                     = large_integer ! Integer

    if(associated(cpoly%area                        )) cpoly%area                        = large_real
    if(associated(cpoly%patch_area                  )) cpoly%patch_area                  = large_real
    if(associated(cpoly%elevation                   )) cpoly%elevation                   = large_real
    if(associated(cpoly%slope                       )) cpoly%slope                       = large_real
    if(associated(cpoly%aspect                      )) cpoly%aspect                      = large_real

    if(associated(cpoly%num_landuse_years           )) cpoly%num_landuse_years           = large_integer  ! Integer
    if(associated(cpoly%lai_pft                     )) cpoly%lai_pft                     = large_real
    if(associated(cpoly%wpa_pft                     )) cpoly%wpa_pft                     = large_real
    if(associated(cpoly%wai_pft                     )) cpoly%wai_pft                     = large_real

    if(associated(cpoly%TCI                         )) cpoly%TCI                         = large_real
    if(associated(cpoly%pptweight                   )) cpoly%pptweight                   = large_real
    if(associated(cpoly%lsl                         )) cpoly%lsl                         = large_integer  ! Integer
    if(associated(cpoly%hydro_next                  )) cpoly%hydro_next                  = large_integer  ! Integer
    if(associated(cpoly%hydro_prev                  )) cpoly%hydro_prev                  = large_integer  ! Integer
    if(associated(cpoly%moist_W                     )) cpoly%moist_W                     = large_real
    if(associated(cpoly%moist_f                     )) cpoly%moist_f                     = large_real
    if(associated(cpoly%moist_tau                   )) cpoly%moist_tau                   = large_real
    if(associated(cpoly%moist_zi                    )) cpoly%moist_zi                    = large_real
    if(associated(cpoly%baseflow                    )) cpoly%baseflow                    = large_real
    if(associated(cpoly%ntext_soil                  )) cpoly%ntext_soil                  = large_integer  ! Integer
    if(associated(cpoly%min_monthly_temp            )) cpoly%min_monthly_temp            = large_real
    if(associated(cpoly%plantation                  )) cpoly%plantation                  = large_integer  ! Integer
    if(associated(cpoly%agri_stocking_pft           )) cpoly%agri_stocking_pft           = large_integer  ! Integer
    if(associated(cpoly%agri_stocking_density       )) cpoly%agri_stocking_density       = large_real
    if(associated(cpoly%plantation_stocking_pft     )) cpoly%plantation_stocking_pft     = large_integer  ! Integer
    if(associated(cpoly%plantation_stocking_density )) cpoly%plantation_stocking_density = large_real
    if(associated(cpoly%primary_harvest_memory      )) cpoly%primary_harvest_memory      = large_real
    if(associated(cpoly%secondary_harvest_memory    )) cpoly%secondary_harvest_memory    = large_real
    if(associated(cpoly%fire_disturbance_rate       )) cpoly%fire_disturbance_rate       = large_real
    if(associated(cpoly%ignition_rate               )) cpoly%ignition_rate               = large_real
    if(associated(cpoly%lambda_fire                 )) cpoly%lambda_fire                 = large_real

    if(associated(cpoly%phen_pars                   )) then
       cpoly%phen_pars(:)%nyears = large_integer
       ! The other variables are also pointers...
    end if

    if(associated(cpoly%treefall_disturbance_rate   )) cpoly%treefall_disturbance_rate   = large_real
    if(associated(cpoly%nat_disturbance_rate        )) cpoly%nat_disturbance_rate        = large_real
    if(associated(cpoly%nat_dist_type               )) cpoly%nat_dist_type               = large_integer  ! Integer
    if(associated(cpoly%disturbance_memory          )) cpoly%disturbance_memory          = large_real
    if(associated(cpoly%disturbance_rates           )) cpoly%disturbance_rates           = large_real
    if(associated(cpoly%loss_fraction               )) cpoly%loss_fraction               = large_real
    if(associated(cpoly%green_leaf_factor           )) cpoly%green_leaf_factor           = large_real
    if(associated(cpoly%leaf_aging_factor           )) cpoly%leaf_aging_factor           = large_real
    if(associated(cpoly%met                         )) then
       cpoly%met(:)%nir_beam          = large_real
       cpoly%met(:)%nir_diffuse       = large_real
       cpoly%met(:)%par_beam          = large_real
       cpoly%met(:)%par_diffuse       = large_real
       cpoly%met(:)%atm_tmp           = large_real
       cpoly%met(:)%atm_theta         = large_real
       cpoly%met(:)%atm_enthalpy      = large_real
       cpoly%met(:)%atm_shv           = large_real
       cpoly%met(:)%theta             = large_real
       cpoly%met(:)%rshort            = large_real
       cpoly%met(:)%rshort_diffuse    = large_real
       cpoly%met(:)%rlong             = large_real
       cpoly%met(:)%pcpg              = large_real
       cpoly%met(:)%qpcpg             = large_real
       cpoly%met(:)%dpcpg             = large_real
       cpoly%met(:)%vels              = large_real
       cpoly%met(:)%vels_stab         = large_real
       cpoly%met(:)%vels_unstab       = large_real
       cpoly%met(:)%prss              = large_real
       cpoly%met(:)%exner             = large_real
       cpoly%met(:)%geoht             = large_real
       cpoly%met(:)%atm_co2           = large_real
    end if
    if(associated(cpoly%basal_area                  )) cpoly%basal_area                  = large_real
    if(associated(cpoly%agb                         )) cpoly%agb                         = large_real
    if(associated(cpoly%pldens                      )) cpoly%pldens                      = large_real
    if(associated(cpoly%bseeds                      )) cpoly%bseeds                      = large_real
    
    if(associated(cpoly%basal_area_growth           )) cpoly%basal_area_growth           = large_real
    if(associated(cpoly%agb_growth                  )) cpoly%agb_growth                  = large_real
    if(associated(cpoly%basal_area_mort             )) cpoly%basal_area_mort             = large_real
    if(associated(cpoly%agb_mort                    )) cpoly%agb_mort                    = large_real
    if(associated(cpoly%basal_area_cut              )) cpoly%basal_area_cut              = large_real!NOT IN REGISTRY
    if(associated(cpoly%agb_cut                     )) cpoly%agb_cut                     = large_real

    if(associated(cpoly%cosaoi                      )) cpoly%cosaoi                      = large_real
    if(associated(cpoly%albedo_beam                 )) cpoly%albedo_beam                 = large_real
    if(associated(cpoly%albedo_diffuse              )) cpoly%albedo_diffuse              = large_real
    if(associated(cpoly%rlong_albedo                )) cpoly%rlong_albedo                = large_real
    
    if(associated(cpoly%albedt                      )) cpoly%albedt                      = large_real
    if(associated(cpoly%rlongup                     )) cpoly%rlongup                     = large_real
    if(associated(cpoly%daylight                    )) cpoly%daylight                    = large_real
    
    if(associated(cpoly%lai                         )) cpoly%lai                         = large_real
    if(associated(cpoly%avg_lma                     )) cpoly%avg_lma                     = large_real
    if(associated(cpoly%wpa                         )) cpoly%wpa                         = large_real
    if(associated(cpoly%wai                         )) cpoly%wai                         = large_real

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    if(associated(cpoly%avg_vapor_vc                )) cpoly%avg_vapor_vc                = large_real
    if(associated(cpoly%avg_dew_cg                  )) cpoly%avg_dew_cg                  = large_real
    if(associated(cpoly%avg_vapor_gc                )) cpoly%avg_vapor_gc                = large_real
    if(associated(cpoly%avg_wshed_vg                )) cpoly%avg_wshed_vg                = large_real
    if(associated(cpoly%avg_intercepted             )) cpoly%avg_intercepted             = large_real
    if(associated(cpoly%avg_vapor_ac                )) cpoly%avg_vapor_ac                = large_real
    if(associated(cpoly%avg_transp                  )) cpoly%avg_transp                  = large_real
    if(associated(cpoly%avg_evap                    )) cpoly%avg_evap                    = large_real
    if(associated(cpoly%avg_smoist_gg               )) cpoly%avg_smoist_gg               = large_real
    if(associated(cpoly%avg_smoist_gc               )) cpoly%avg_smoist_gc               = large_real
    if(associated(cpoly%avg_runoff                  )) cpoly%avg_runoff                  = large_real
    if(associated(cpoly%avg_drainage                )) cpoly%avg_drainage                = large_real
    if(associated(cpoly%avg_drainage_heat           )) cpoly%avg_drainage_heat           = large_real
    if(associated(cpoly%aux                         )) cpoly%aux                         = large_real
    if(associated(cpoly%aux_s                       )) cpoly%aux_s                       = large_real
    if(associated(cpoly%avg_sensible_vc             )) cpoly%avg_sensible_vc             = large_real
    if(associated(cpoly%avg_qwshed_vg               )) cpoly%avg_qwshed_vg               = large_real
    if(associated(cpoly%avg_qintercepted            )) cpoly%avg_qintercepted            = large_real
    if(associated(cpoly%avg_sensible_gc             )) cpoly%avg_sensible_gc             = large_real
    if(associated(cpoly%avg_sensible_ac             )) cpoly%avg_sensible_ac             = large_real
    if(associated(cpoly%avg_sensible_gg             )) cpoly%avg_sensible_gg             = large_real
    if(associated(cpoly%avg_runoff_heat             )) cpoly%avg_runoff_heat             = large_real
    if(associated(cpoly%avg_veg_energy              )) cpoly%avg_veg_energy              = large_real
    if(associated(cpoly%avg_veg_hcap                )) cpoly%avg_veg_hcap                = large_real
    if(associated(cpoly%avg_veg_temp                )) cpoly%avg_veg_temp                = large_real
    if(associated(cpoly%avg_veg_fliq                )) cpoly%avg_veg_fliq                = large_real
    if(associated(cpoly%avg_veg_water               )) cpoly%avg_veg_water               = large_real
    if(associated(cpoly%avg_can_temp                )) cpoly%avg_can_temp                = large_real
    if(associated(cpoly%avg_can_shv                 )) cpoly%avg_can_shv                 = large_real
    if(associated(cpoly%avg_can_co2                 )) cpoly%avg_can_co2                 = large_real
    if(associated(cpoly%avg_can_rhos                )) cpoly%avg_can_rhos                = large_real
    if(associated(cpoly%avg_can_prss                )) cpoly%avg_can_prss                = large_real
    if(associated(cpoly%avg_can_theta               )) cpoly%avg_can_theta               = large_real
    if(associated(cpoly%avg_can_enthalpy            )) cpoly%avg_can_enthalpy            = large_real
    if(associated(cpoly%avg_can_depth               )) cpoly%avg_can_depth               = large_real
    if(associated(cpoly%avg_soil_energy             )) cpoly%avg_soil_energy             = large_real
    if(associated(cpoly%avg_soil_water              )) cpoly%avg_soil_water              = large_real
    if(associated(cpoly%avg_soil_temp               )) cpoly%avg_soil_temp               = large_real
    if(associated(cpoly%avg_soil_fracliq            )) cpoly%avg_soil_fracliq            = large_real
    if(associated(cpoly%avg_soil_wetness            )) cpoly%avg_soil_wetness            = large_real
    if(associated(cpoly%avg_skin_temp               )) cpoly%avg_skin_temp               = large_real
    if(associated(cpoly%avg_available_water         )) cpoly%avg_available_water         = large_real
    if(associated(cpoly%runoff                      )) cpoly%runoff                      = large_real
    if(associated(cpoly%avg_atm_tmp                 )) cpoly%avg_atm_tmp                 = large_real
    if(associated(cpoly%avg_atm_shv                 )) cpoly%avg_atm_shv                 = large_real
    if(associated(cpoly%avg_atm_prss                )) cpoly%avg_atm_prss                = large_real

    !NACP
    if(associated(cpoly%avg_sfcw_depth              )) cpoly%avg_sfcw_depth              = large_real
    if(associated(cpoly%avg_sfcw_energy             )) cpoly%avg_sfcw_energy             = large_real
    if(associated(cpoly%avg_sfcw_mass               )) cpoly%avg_sfcw_mass               = large_real
    if(associated(cpoly%avg_sfcw_tempk              )) cpoly%avg_sfcw_tempk              = large_real
    if(associated(cpoly%avg_sfcw_fracliq            )) cpoly%avg_sfcw_fracliq            = large_real
    if(associated(cpoly%avg_fsc                     )) cpoly%avg_fsc                     = large_real
    if(associated(cpoly%avg_stsc                    )) cpoly%avg_stsc                    = large_real
    if(associated(cpoly%avg_ssc                     )) cpoly%avg_ssc                     = large_real
    if(associated(cpoly%avg_balive                  )) cpoly%avg_balive                  = large_real
    if(associated(cpoly%avg_bleaf                   )) cpoly%avg_bleaf                   = large_real
    if(associated(cpoly%avg_bsapwood                )) cpoly%avg_bsapwood                = large_real
    if(associated(cpoly%avg_broot                   )) cpoly%avg_broot                   = large_real
    if(associated(cpoly%avg_bdead                   )) cpoly%avg_bdead                   = large_real
    if(associated(cpoly%avg_fsn                     )) cpoly%avg_fsn                     = large_real
    if(associated(cpoly%avg_msn                     )) cpoly%avg_msn                     = large_real
    if(associated(cpoly%avg_bstorage                )) cpoly%avg_bstorage                = large_real
    if(associated(cpoly%avg_bseeds                  )) cpoly%avg_bseeds                  = large_real

    if(associated(cpoly%dmean_co2_residual          )) cpoly%dmean_co2_residual       = large_real
    if(associated(cpoly%dmean_energy_residual       )) cpoly%dmean_energy_residual    = large_real
    if(associated(cpoly%dmean_water_residual        )) cpoly%dmean_water_residual     = large_real
    if(associated(cpoly%mmean_co2_residual          )) cpoly%mmean_co2_residual       = large_real
    if(associated(cpoly%mmean_energy_residual       )) cpoly%mmean_energy_residual    = large_real
    if(associated(cpoly%mmean_water_residual        )) cpoly%mmean_water_residual     = large_real

    return
  end subroutine huge_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine huge_sitetype(csite)

    implicit none
    type(sitetype),target :: csite
    
    if(associated(csite%paco_id                      )) csite%paco_id                      = large_integer ! Integer
    if(associated(csite%paco_n                       )) csite%paco_n                       = large_integer ! Integer

    if(associated(csite%dist_type                    )) csite%dist_type                    = large_integer ! Integer
    if(associated(csite%age                          )) csite%age                          = large_real
    if(associated(csite%area                         )) csite%area                         = large_real
    if(associated(csite%laiarea                      )) csite%laiarea                      = large_real
    if(associated(csite%hcapveg                      )) csite%hcapveg                      = large_real
    if(associated(csite%fast_soil_C                  )) csite%fast_soil_C                  = large_real
    if(associated(csite%slow_soil_C                  )) csite%slow_soil_C                  = large_real
    if(associated(csite%structural_soil_C            )) csite%structural_soil_C            = large_real
    if(associated(csite%structural_soil_L            )) csite%structural_soil_L            = large_real
    if(associated(csite%mineralized_soil_N           )) csite%mineralized_soil_N           = large_real
    if(associated(csite%fast_soil_N                  )) csite%fast_soil_N                  = large_real
    if(associated(csite%pname                        )) csite%pname                        = ''      ! Character
    if(associated(csite%sum_dgd                      )) csite%sum_dgd                      = large_real
    if(associated(csite%sum_chd                      )) csite%sum_chd                      = large_real
    if(associated(csite%plantation                   )) csite%plantation                   = large_integer ! Integer
    if(associated(csite%cohort_count                 )) csite%cohort_count                 = large_integer ! Integer
    if(associated(csite%can_enthalpy                 )) csite%can_enthalpy                 = large_real
    if(associated(csite%can_temp                     )) csite%can_temp                     = large_real
    if(associated(csite%can_shv                      )) csite%can_shv                      = large_real
    if(associated(csite%can_co2                      )) csite%can_co2                      = large_real
    if(associated(csite%can_rhos                     )) csite%can_rhos                     = large_real
    if(associated(csite%can_prss                     )) csite%can_prss                     = large_real
    if(associated(csite%can_theta                    )) csite%can_theta                    = large_real
    if(associated(csite%can_depth                    )) csite%can_depth                    = large_real
    if(associated(csite%lambda_light                 )) csite%lambda_light                 = large_real
    if(associated(csite%dmean_lambda_light           )) csite%dmean_lambda_light           = large_real
    if(associated(csite%mmean_lambda_light           )) csite%mmean_lambda_light           = large_real
    if(associated(csite%lai                          )) csite%lai                          = large_real
    if(associated(csite%wpa                          )) csite%wpa                          = large_real
    if(associated(csite%wai                          )) csite%wai                          = large_real

    if(associated(csite%sfcwater_mass                )) csite%sfcwater_mass                = large_real
    if(associated(csite%sfcwater_energy              )) csite%sfcwater_energy              = large_real
    if(associated(csite%sfcwater_depth               )) csite%sfcwater_depth               = large_real
    if(associated(csite%rshort_s                     )) csite%rshort_s                     = large_real
    if(associated(csite%rshort_s_beam                )) csite%rshort_s_beam                = large_real
    if(associated(csite%rshort_s_diffuse             )) csite%rshort_s_diffuse             = large_real
    if(associated(csite%sfcwater_tempk               )) csite%sfcwater_tempk               = large_real
    if(associated(csite%sfcwater_fracliq             )) csite%sfcwater_fracliq             = large_real
    if(associated(csite%nlev_sfcwater                )) csite%nlev_sfcwater                = large_integer ! Integer
    if(associated(csite%ntext_soil                   )) csite%ntext_soil                   = large_integer ! Integer
    if(associated(csite%soil_energy                  )) csite%soil_energy                  = large_real
    if(associated(csite%soil_water                   )) csite%soil_water                   = large_double
    if(associated(csite%soil_tempk                   )) csite%soil_tempk                   = large_real
    if(associated(csite%soil_fracliq                 )) csite%soil_fracliq                 = large_real
    if(associated(csite%ground_shv                   )) csite%ground_shv                   = large_real
    if(associated(csite%surface_ssh                  )) csite%surface_ssh                  = large_real
    if(associated(csite%rough                        )) csite%rough                        = large_real
    if(associated(csite%A_o_max                      )) csite%A_o_max                      = large_real
    if(associated(csite%A_c_max                      )) csite%A_c_max                      = large_real
    if(associated(csite%old_stoma_data_max)) then
       csite%old_stoma_data_max(:,:)%recalc             = 1        ! Integer
       csite%old_stoma_data_max(:,:)%T_L                = large_real
       csite%old_stoma_data_max(:,:)%e_A                = large_real
       csite%old_stoma_data_max(:,:)%PAR                = large_real
       csite%old_stoma_data_max(:,:)%rb_factor          = large_real
       csite%old_stoma_data_max(:,:)%prss               = large_real
       csite%old_stoma_data_max(:,:)%phenology_factor   = large_real
       csite%old_stoma_data_max(:,:)%gsw_open           = large_real
       csite%old_stoma_data_max(:,:)%ilimit             = large_integer  ! Integer
       csite%old_stoma_data_max(:,:)%T_L_residual       = large_real
       csite%old_stoma_data_max(:,:)%e_a_residual       = large_real
       csite%old_stoma_data_max(:,:)%par_residual       = large_real
       csite%old_stoma_data_max(:,:)%rb_residual        = large_real
       csite%old_stoma_data_max(:,:)%prss_residual      = large_real
       csite%old_stoma_data_max(:,:)%leaf_residual      = large_real
       csite%old_stoma_data_max(:,:)%gsw_residual       = large_real
    end if
    if(associated(csite%old_stoma_vector_max         )) csite%old_stoma_vector_max(:,:,:)  = large_real
 

    if(associated(csite%avg_daily_temp               )) csite%avg_daily_temp               = large_real
    if(associated(csite%mean_rh                      )) csite%mean_rh                      = large_real
    if(associated(csite%dmean_rh                     )) csite%dmean_rh                     = large_real
    if(associated(csite%mmean_rh                     )) csite%mmean_rh                     = large_real
    if(associated(csite%mean_nep                     )) csite%mean_nep                     = large_real
    if(associated(csite%wbudget_loss2atm             )) csite%wbudget_loss2atm             = large_real
    if(associated(csite%wbudget_denseffect           )) csite%wbudget_denseffect           = large_real
    if(associated(csite%wbudget_precipgain           )) csite%wbudget_precipgain           = large_real
    if(associated(csite%wbudget_loss2runoff          )) csite%wbudget_loss2runoff          = large_real
    if(associated(csite%wbudget_loss2drainage        )) csite%wbudget_loss2drainage        = large_real
    if(associated(csite%wbudget_initialstorage       )) csite%wbudget_initialstorage       = large_real
    if(associated(csite%wbudget_residual             )) csite%wbudget_residual             = large_real
    if(associated(csite%ebudget_latent               )) csite%ebudget_latent               = large_real
    if(associated(csite%ebudget_loss2atm             )) csite%ebudget_loss2atm             = large_real
    if(associated(csite%ebudget_denseffect           )) csite%ebudget_denseffect           = large_real
    if(associated(csite%ebudget_loss2runoff          )) csite%ebudget_loss2runoff          = large_real
    if(associated(csite%ebudget_loss2drainage        )) csite%ebudget_loss2drainage        = large_real
    if(associated(csite%ebudget_netrad               )) csite%ebudget_netrad               = large_real
    if(associated(csite%ebudget_precipgain           )) csite%ebudget_precipgain           = large_real
    if(associated(csite%ebudget_initialstorage       )) csite%ebudget_initialstorage       = large_real
    if(associated(csite%ebudget_residual             )) csite%ebudget_residual             = large_real
    if(associated(csite%co2budget_initialstorage     )) csite%co2budget_initialstorage     = large_real
    if(associated(csite%co2budget_residual           )) csite%co2budget_residual           = large_real
    if(associated(csite%co2budget_loss2atm           )) csite%co2budget_loss2atm           = large_real
    if(associated(csite%co2budget_denseffect         )) csite%co2budget_denseffect         = large_real
    if(associated(csite%co2budget_gpp                )) csite%co2budget_gpp                = large_real
    if(associated(csite%co2budget_gpp_dbh            )) csite%co2budget_gpp_dbh            = large_real
    if(associated(csite%co2budget_plresp             )) csite%co2budget_plresp             = large_real
    if(associated(csite%co2budget_rh                 )) csite%co2budget_rh                 = large_real
    if(associated(csite%today_A_decomp               )) csite%today_A_decomp               = large_real
    if(associated(csite%today_Af_decomp              )) csite%today_Af_decomp              = large_real
    if(associated(csite%dmean_A_decomp               )) csite%dmean_A_decomp               = large_real
    if(associated(csite%dmean_Af_decomp              )) csite%dmean_Af_decomp              = large_real
    if(associated(csite%mmean_A_decomp               )) csite%mmean_A_decomp               = large_real
    if(associated(csite%mmean_Af_decomp              )) csite%mmean_Af_decomp              = large_real
    if(associated(csite%repro                        )) csite%repro                        = large_real
    if(associated(csite%veg_rough                    )) csite%veg_rough                    = large_real
    if(associated(csite%veg_height                   )) csite%veg_height                   = large_real
    if(associated(csite%fsc_in                       )) csite%fsc_in                       = large_real
    if(associated(csite%ssc_in                       )) csite%ssc_in                       = large_real
    if(associated(csite%ssl_in                       )) csite%ssl_in                       = large_real
    if(associated(csite%fsn_in                       )) csite%fsn_in                       = large_real
    if(associated(csite%total_plant_nitrogen_uptake  )) csite%total_plant_nitrogen_uptake  = large_real
    if(associated(csite%mineralized_N_input  )) csite%mineralized_N_input  = large_real
    if(associated(csite%mineralized_N_loss  )) csite%mineralized_N_loss  = large_real
    if(associated(csite%rshort_g                     )) csite%rshort_g                     = large_real
    if(associated(csite%rshort_g_beam                )) csite%rshort_g_beam                = large_real
    if(associated(csite%rshort_g_diffuse             )) csite%rshort_g_diffuse             = large_real
    if(associated(csite%rlong_g                      )) csite%rlong_g                      = large_real
    if(associated(csite%rlong_g_surf                 )) csite%rlong_g_surf                 = large_real
    if(associated(csite%rlong_g_incid                )) csite%rlong_g_incid                = large_real
    if(associated(csite%rlong_s                      )) csite%rlong_s                      = large_real
    if(associated(csite%rlong_s_surf                 )) csite%rlong_s_surf                 = large_real
    if(associated(csite%rlong_s_incid                )) csite%rlong_s_incid                = large_real
    if(associated(csite%albedt                       )) csite%albedt                       = large_real
    if(associated(csite%albedo_beam                  )) csite%albedo_beam                  = large_real
    if(associated(csite%albedo_diffuse               )) csite%albedo_diffuse               = large_real
    if(associated(csite%rlongup                      )) csite%rlongup                      = large_real
    if(associated(csite%rlong_albedo                 )) csite%rlong_albedo                 = large_real
    if(associated(csite%total_snow_depth             )) csite%total_snow_depth             = large_real
    if(associated(csite%snowfac                      )) csite%snowfac                      = large_real
    if(associated(csite%A_decomp                     )) csite%A_decomp                     = large_real
    if(associated(csite%f_decomp                     )) csite%f_decomp                     = large_real
    if(associated(csite%rh                           )) csite%rh                           = large_real
    if(associated(csite%cwd_rh                       )) csite%cwd_rh                       = large_real
    if(associated(csite%fuse_flag                    )) csite%fuse_flag                    = large_integer !Integer
    if(associated(csite%pft_density_profile          )) csite%pft_density_profile          = large_real
    if(associated(csite%plant_ag_biomass             )) csite%plant_ag_biomass             = large_real

    if(associated(csite%mean_wflux                   )) csite%mean_wflux                   = large_real
    if(associated(csite%mean_latflux                 )) csite%mean_latflux                 = large_real
    if(associated(csite%mean_hflux                   )) csite%mean_hflux                   = large_real
    if(associated(csite%mean_runoff                  )) csite%mean_runoff                  = large_real
    if(associated(csite%mean_qrunoff                 )) csite%mean_qrunoff                 = large_real

    if(associated(csite%htry                         )) csite%htry                         = large_real
    if(associated(csite%avg_rk4step                  )) csite%avg_rk4step                  = large_real
    if(associated(csite%dmean_rk4step                )) csite%dmean_rk4step                = large_real
    if(associated(csite%mmean_rk4step                )) csite%mmean_rk4step                = large_real

    if(associated(csite%avg_available_water          )) csite%avg_available_water          = large_real

    if(associated(csite%ustar                        )) csite%ustar                        = large_real
    if(associated(csite%tstar                        )) csite%tstar                        = large_real
    if(associated(csite%qstar                        )) csite%qstar                        = large_real
    if(associated(csite%cstar                        )) csite%cstar                        = large_real

    if(associated(csite%upwp                         )) csite%upwp                         = large_real
    if(associated(csite%tpwp                         )) csite%tpwp                         = large_real
    if(associated(csite%qpwp                         )) csite%qpwp                         = large_real
    if(associated(csite%cpwp                         )) csite%cpwp                         = large_real
    if(associated(csite%wpwp                         )) csite%wpwp                         = large_real

    if(associated(csite%avg_carbon_ac                )) csite%avg_carbon_ac                = large_real

    if(associated(csite%avg_vapor_vc                 )) csite%avg_vapor_vc                 = large_real
    if(associated(csite%avg_dew_cg                   )) csite%avg_dew_cg                   = large_real
    if(associated(csite%avg_vapor_gc                 )) csite%avg_vapor_gc                 = large_real
    if(associated(csite%avg_wshed_vg                 )) csite%avg_wshed_vg                 = large_real
    if(associated(csite%avg_intercepted              )) csite%avg_intercepted              = large_real
    if(associated(csite%avg_vapor_ac                 )) csite%avg_vapor_ac                 = large_real
    if(associated(csite%avg_transp                   )) csite%avg_transp                   = large_real
    if(associated(csite%avg_evap                     )) csite%avg_evap                     = large_real
    if(associated(csite%avg_netrad                   )) csite%avg_netrad                   = large_real
    if(associated(csite%avg_smoist_gg                )) csite%avg_smoist_gg                = large_real
    if(associated(csite%avg_smoist_gc                )) csite%avg_smoist_gc                = large_real
    if(associated(csite%avg_runoff                   )) csite%avg_runoff                   = large_real
    if(associated(csite%avg_drainage                 )) csite%avg_drainage                 = large_real
    if(associated(csite%avg_drainage_heat            )) csite%avg_drainage_heat            = large_real
    if(associated(csite%aux                          )) csite%aux                          = large_real
    if(associated(csite%aux_s                        )) csite%aux_s                        = large_real
    if(associated(csite%avg_sensible_vc              )) csite%avg_sensible_vc              = large_real
    if(associated(csite%avg_qwshed_vg                )) csite%avg_qwshed_vg                = large_real
    if(associated(csite%avg_qintercepted             )) csite%avg_qintercepted             = large_real
    if(associated(csite%avg_sensible_gc              )) csite%avg_sensible_gc              = large_real
    if(associated(csite%avg_sensible_ac              )) csite%avg_sensible_ac              = large_real
    if(associated(csite%avg_sensible_gg              )) csite%avg_sensible_gg              = large_real
    if(associated(csite%avg_runoff_heat              )) csite%avg_runoff_heat              = large_real
    if(associated(csite%avg_veg_energy               )) csite%avg_veg_energy               = large_real
    if(associated(csite%avg_veg_temp                 )) csite%avg_veg_temp                 = large_real
    if(associated(csite%avg_veg_fliq                 )) csite%avg_veg_fliq                 = large_real
    if(associated(csite%avg_veg_water                )) csite%avg_veg_water                = large_real

    if(associated(csite%watertable                   )) csite%watertable                   = large_real
    if(associated(csite%moist_dz                     )) csite%moist_dz                     = large_real
    if(associated(csite%ksat                         )) csite%ksat                         = large_real
    if(associated(csite%soil_sat_energy              )) csite%soil_sat_energy              = large_real
    if(associated(csite%soil_sat_water               )) csite%soil_sat_water               = large_real
    if(associated(csite%soil_sat_heat                )) csite%soil_sat_heat                = large_real
    if(associated(csite%runoff_A                     )) csite%runoff_A                     = large_real
    if(associated(csite%runoff_rate                  )) csite%runoff_rate                  = large_real
    if(associated(csite%runoff                       )) csite%runoff                       = large_real

    if(associated(csite%dmean_co2_residual      )) csite%dmean_co2_residual       = large_real
    if(associated(csite%dmean_energy_residual   )) csite%dmean_energy_residual    = large_real
    if(associated(csite%dmean_water_residual    )) csite%dmean_water_residual     = large_real
    if(associated(csite%mmean_co2_residual      )) csite%mmean_co2_residual       = large_real
    if(associated(csite%mmean_energy_residual   )) csite%mmean_energy_residual    = large_real
    if(associated(csite%mmean_water_residual    )) csite%mmean_water_residual     = large_real

    return
  end subroutine huge_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine huge_patchtype(cpatch)

    implicit none
    type(patchtype), target :: cpatch
    
    if(associated(cpatch%pft))                  cpatch%pft                 = large_integer    !Integer
    if(associated(cpatch%nplant))               cpatch%nplant              = large_real
    if(associated(cpatch%hite))                 cpatch%hite                = large_real
    if(associated(cpatch%agb))                  cpatch%agb                 = large_real
    if(associated(cpatch%basarea))              cpatch%basarea             = large_real
    if(associated(cpatch%dagb_dt))              cpatch%dagb_dt             = large_real
    if(associated(cpatch%dba_dt))               cpatch%dba_dt              = large_real
    if(associated(cpatch%ddbh_dt))              cpatch%ddbh_dt             = large_real
    if(associated(cpatch%dbh))                  cpatch%dbh                 = large_real
    if(associated(cpatch%bdead))                cpatch%bdead               = large_real
    if(associated(cpatch%bleaf))                cpatch%bleaf               = large_real
    if(associated(cpatch%phenology_status))     cpatch%phenology_status    = large_integer !Integer
    if(associated(cpatch%balive))               cpatch%balive              = large_real
    if(associated(cpatch%broot))                cpatch%broot               = large_real
    if(associated(cpatch%bsapwood))             cpatch%bsapwood            = large_real
    if(associated(cpatch%lai))                  cpatch%lai                 = large_real
    if(associated(cpatch%wpa))                  cpatch%wpa                 = large_real
    if(associated(cpatch%wai))                  cpatch%wai                 = large_real
    if(associated(cpatch%solvable))             cpatch%solvable            = .true.
    if(associated(cpatch%bstorage))             cpatch%bstorage            = large_real
    if(associated(cpatch%cb))                   cpatch%cb                  = large_real
    if(associated(cpatch%cb_max))               cpatch%cb_max              = large_real
    if(associated(cpatch%cbr_bar))              cpatch%cbr_bar             = large_real
    if(associated(cpatch%mmean_cb))             cpatch%mmean_cb            = large_real
    if(associated(cpatch%veg_energy))           cpatch%veg_energy          = large_real
    if(associated(cpatch%veg_temp))             cpatch%veg_temp            = large_real
    if(associated(cpatch%veg_fliq))             cpatch%veg_fliq            = large_real
    if(associated(cpatch%veg_water))            cpatch%veg_water           = large_real
    if(associated(cpatch%mean_gpp))             cpatch%mean_gpp            = large_real
    if(associated(cpatch%mean_leaf_resp))       cpatch%mean_leaf_resp      = large_real
    if(associated(cpatch%mean_root_resp))       cpatch%mean_root_resp      = large_real
    if(associated(cpatch%mean_growth_resp))     cpatch%mean_growth_resp      = large_real
    if(associated(cpatch%mean_storage_resp))    cpatch%mean_storage_resp      = large_real
    if(associated(cpatch%mean_vleaf_resp))      cpatch%mean_vleaf_resp      = large_real
    if(associated(cpatch%dmean_leaf_resp))      cpatch%dmean_leaf_resp     = large_real
    if(associated(cpatch%today_leaf_resp))      cpatch%today_leaf_resp     = large_real
    if(associated(cpatch%today_root_resp))      cpatch%today_root_resp     = large_real
    if(associated(cpatch%today_gpp))            cpatch%today_gpp           = large_real
    if(associated(cpatch%today_gpp_pot))        cpatch%today_gpp_pot       = large_real
    if(associated(cpatch%today_gpp_max))        cpatch%today_gpp_max       = large_real
    if(associated(cpatch%growth_respiration))   cpatch%growth_respiration  = large_real
    if(associated(cpatch%storage_respiration))  cpatch%storage_respiration = large_real
    if(associated(cpatch%vleaf_respiration))    cpatch%vleaf_respiration   = large_real
    if(associated(cpatch%dmean_gpp         ))   cpatch%dmean_gpp           = large_real
    if(associated(cpatch%dmean_leaf_resp   ))   cpatch%dmean_leaf_resp     = large_real
    if(associated(cpatch%dmean_root_resp   ))   cpatch%dmean_root_resp     = large_real
    if(associated(cpatch%mmean_gpp         ))   cpatch%mmean_gpp           = large_real
    if(associated(cpatch%mmean_leaf_resp   ))   cpatch%mmean_leaf_resp     = large_real
    if(associated(cpatch%mmean_root_resp   ))   cpatch%mmean_root_resp     = large_real
    if(associated(cpatch%mmean_growth_resp ))   cpatch%mmean_growth_resp   = large_real
    if(associated(cpatch%mmean_storage_resp))   cpatch%mmean_storage_resp  = large_real
    if(associated(cpatch%mmean_vleaf_resp  ))   cpatch%mmean_vleaf_resp    = large_real
    if(associated(cpatch%fsn))                  cpatch%fsn                 = large_real
    if(associated(cpatch%monthly_dndt))         cpatch%monthly_dndt        = large_real
    if(associated(cpatch%mort_rate))            cpatch%mort_rate           = large_real
    if(associated(cpatch%mmean_mort_rate))      cpatch%mmean_mort_rate     = large_real

    if(associated(cpatch%old_stoma_data)) then
       cpatch%old_stoma_data(:)%recalc             = 1
       cpatch%old_stoma_data(:)%T_L                = large_real
       cpatch%old_stoma_data(:)%e_A                = large_real
       cpatch%old_stoma_data(:)%PAR                = large_real
       cpatch%old_stoma_data(:)%rb_factor          = large_real
       cpatch%old_stoma_data(:)%prss               = large_real
       cpatch%old_stoma_data(:)%phenology_factor   = large_real
       cpatch%old_stoma_data(:)%gsw_open           = large_real
       cpatch%old_stoma_data(:)%ilimit             = large_integer
       cpatch%old_stoma_data(:)%T_L_residual       = large_real
       cpatch%old_stoma_data(:)%e_a_residual       = large_real
       cpatch%old_stoma_data(:)%par_residual       = large_real
       cpatch%old_stoma_data(:)%rb_residual        = large_real
       cpatch%old_stoma_data(:)%prss_residual      = large_real
       cpatch%old_stoma_data(:)%leaf_residual      = large_real
       cpatch%old_stoma_data(:)%gsw_residual       = large_real
    end if

    if(associated(cpatch%old_stoma_vector))     cpatch%old_stoma_vector    = large_real

    if(associated(cpatch%Psi_open))             cpatch%Psi_open            = large_real
    if(associated(cpatch%krdepth))              cpatch%krdepth             = large_integer  !Integer
    if(associated(cpatch%first_census))         cpatch%first_census        = large_integer  !Integer
    if(associated(cpatch%new_recruit_flag))     cpatch%new_recruit_flag    = large_integer  !Integer
    if(associated(cpatch%light_level))          cpatch%light_level         = large_real
    if(associated(cpatch%dmean_light_level))    cpatch%dmean_light_level   = large_real
    if(associated(cpatch%mmean_light_level))    cpatch%mmean_light_level   = large_real
    if(associated(cpatch%light_level_beam))       cpatch%light_level_beam       = large_real
    if(associated(cpatch%dmean_light_level_beam)) cpatch%dmean_light_level_beam = large_real
    if(associated(cpatch%mmean_light_level_beam)) cpatch%mmean_light_level_beam = large_real
    if(associated(cpatch%light_level_diff))       cpatch%light_level_diff       = large_real
    if(associated(cpatch%dmean_light_level_diff)) cpatch%dmean_light_level_diff = large_real
    if(associated(cpatch%mmean_light_level_diff)) cpatch%mmean_light_level_diff = large_real
    if(associated(cpatch%beamext_level))          cpatch%beamext_level         = large_real
    if(associated(cpatch%dmean_beamext_level))    cpatch%dmean_beamext_level   = large_real
    if(associated(cpatch%mmean_beamext_level))    cpatch%mmean_beamext_level   = large_real
    if(associated(cpatch%diffext_level))          cpatch%diffext_level         = large_real
    if(associated(cpatch%dmean_diffext_level))    cpatch%dmean_diffext_level   = large_real
    if(associated(cpatch%mmean_diffext_level))    cpatch%mmean_diffext_level   = large_real
    if(associated(cpatch%norm_par_beam))          cpatch%norm_par_beam       = large_real
    if(associated(cpatch%dmean_norm_par_beam))    cpatch%dmean_norm_par_beam = large_real
    if(associated(cpatch%mmean_norm_par_beam))    cpatch%mmean_norm_par_beam = large_real
    if(associated(cpatch%norm_par_diff))          cpatch%norm_par_diff       = large_real
    if(associated(cpatch%dmean_norm_par_diff))    cpatch%dmean_norm_par_diff = large_real
    if(associated(cpatch%mmean_norm_par_diff))    cpatch%mmean_norm_par_diff = large_real
    if(associated(cpatch%lambda_light))          cpatch%lambda_light         = large_real
    if(associated(cpatch%dmean_lambda_light))    cpatch%dmean_lambda_light   = large_real
    if(associated(cpatch%mmean_lambda_light))    cpatch%mmean_lambda_light   = large_real
    if(associated(cpatch%par_v))                cpatch%par_v               = large_real
    if(associated(cpatch%par_v_beam))           cpatch%par_v_beam          = large_real
    if(associated(cpatch%par_v_diffuse))        cpatch%par_v_diffuse       = large_real
    if(associated(cpatch%dmean_par_v))          cpatch%dmean_par_v               = large_real
    if(associated(cpatch%dmean_par_v_beam))     cpatch%dmean_par_v_beam          = large_real
    if(associated(cpatch%dmean_par_v_diff))     cpatch%dmean_par_v_diff      = large_real
    if(associated(cpatch%mmean_par_v))          cpatch%mmean_par_v               = large_real
    if(associated(cpatch%mmean_par_v_beam))     cpatch%mmean_par_v_beam          = large_real
    if(associated(cpatch%mmean_par_v_diff))     cpatch%mmean_par_v_diff      = large_real
    if(associated(cpatch%rshort_v))             cpatch%rshort_v            = large_real
    if(associated(cpatch%rshort_v_beam))        cpatch%rshort_v_beam       = large_real
    if(associated(cpatch%rshort_v_diffuse))     cpatch%rshort_v_diffuse    = large_real
    if(associated(cpatch%rlong_v))              cpatch%rlong_v             = large_real
    if(associated(cpatch%rlong_v_surf))         cpatch%rlong_v_surf        = large_real
    if(associated(cpatch%rlong_v_incid))        cpatch%rlong_v_incid       = large_real
    if(associated(cpatch%rb))                   cpatch%rb                  = large_real
    if(associated(cpatch%A_open))               cpatch%A_open              = large_real
    if(associated(cpatch%A_closed))             cpatch%A_closed            = large_real
    if(associated(cpatch%Psi_closed))           cpatch%Psi_closed          = large_real
    if(associated(cpatch%rsw_open))             cpatch%rsw_open            = large_real
    if(associated(cpatch%rsw_closed))           cpatch%rsw_closed          = large_real
    if(associated(cpatch%fsw))                  cpatch%fsw                 = large_real
    if(associated(cpatch%fs_open))              cpatch%fs_open             = large_real
    if(associated(cpatch%dmean_fs_open))        cpatch%dmean_fs_open       = large_real
    if(associated(cpatch%dmean_fsw))            cpatch%dmean_fsw           = large_real
    if(associated(cpatch%dmean_fsn))            cpatch%dmean_fsn           = large_real
    if(associated(cpatch%mmean_fs_open))        cpatch%mmean_fs_open       = large_real
    if(associated(cpatch%mmean_fsw))            cpatch%mmean_fsw           = large_real
    if(associated(cpatch%mmean_fsn))            cpatch%mmean_fsn           = large_real
    if(associated(cpatch%stomatal_resistance))  cpatch%stomatal_resistance = large_real
    if(associated(cpatch%leaf_maintenance))     cpatch%leaf_maintenance    = large_real
    if(associated(cpatch%root_maintenance))     cpatch%root_maintenance    = large_real
    if(associated(cpatch%mmean_leaf_maintenance)) cpatch%mmean_leaf_maintenance = large_real
    if(associated(cpatch%mmean_root_maintenance)) cpatch%mmean_root_maintenance = large_real
    if(associated(cpatch%leaf_drop))            cpatch%leaf_drop           = large_real
    if(associated(cpatch%mmean_leaf_drop))      cpatch%mmean_leaf_drop     = large_real
    if(associated(cpatch%bseeds))               cpatch%bseeds              = large_real
    if(associated(cpatch%leaf_respiration))     cpatch%leaf_respiration    = large_real
    if(associated(cpatch%root_respiration))     cpatch%root_respiration    = large_real
    if(associated(cpatch%hcapveg))              cpatch%hcapveg             = large_real
    if(associated(cpatch%gpp))                  cpatch%gpp                 = large_real
    if(associated(cpatch%paw_avg))              cpatch%paw_avg             = large_real
    if(associated(cpatch%turnover_amp))         cpatch%turnover_amp        = large_real
    if(associated(cpatch%llspan))               cpatch%llspan              = large_real
    if(associated(cpatch%vm_bar))               cpatch%vm_bar              = large_real
    if(associated(cpatch%sla))                  cpatch%sla                 = large_real
    ! Depricated
!    if(associated(cpatch%co_srad_h))            cpatch%co_srad_h           = large_real
!    if(associated(cpatch%co_lrad_h))            cpatch%co_lrad_h           = large_real
!    if(associated(cpatch%co_sens_h))            cpatch%co_sens_h           = large_real
!    if(associated(cpatch%co_evap_h))            cpatch%co_evap_h           = large_real
!    if(associated(cpatch%co_liqr_h))            cpatch%co_liqr_h           = large_real 

    return
  end subroutine huge_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_edtype(edin,edout,ipin,ipout)
  ! =====================================================
  ! Copying Functions
  ! =====================================================

    implicit none
    integer :: ipin,ipout

    type(edtype) :: edin,edout

    edout%lat(ipout)          = edin%lat(ipin)
    edout%lon(ipout)          = edin%lon(ipin)
    edout%xatm(ipout)         = edin%xatm(ipin)
    edout%yatm(ipout)         = edin%yatm(ipin)
    edout%natm(ipout)         = edin%natm(ipin)
    edout%ntext_soil(1:nzg,ipout) = edin%ntext_soil(1:nzg,ipin)
    edout%lsl(ipout)          = edin%lsl(ipin)
  
    return
  end subroutine copy_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_polygon(edin,oldsize,water_site)
  ! =====================================================
  ! Copying Functions2
  ! =====================================================

    implicit none
!    integer :: ifm,ipin,ipout,newsize,oldsize

    type(edtype),target :: edin
    integer :: oldsize,newsize,i
    integer, dimension(oldsize),intent(IN) :: water_site
    integer, dimension(oldsize) :: water_site_new

    water_site_new = 0
    
    newsize=0
    do i=1,oldsize
       if (water_site(i)==0) then
          newsize=newsize+1
          water_site_new(newsize)=i
       endif
    enddo

!    edout_real(1:oldsize)     = edin%lat(1:oldsize)
!    deallocate(edin%lat)
!    allocate(edin%lat(newsize))
!    edin%lat(1:newsize)       = edout_real(water_site_new(1:newsize)))
    
!    edout_real(1:oldsize)     = edin%lon(1:oldsize)
!    deallocate(edin%lon)
!    allocate(edin%lon(newsize))
!    edin%lon(1:newsize)       = edout_real(water_site_new(1:newsize)))

!    edout_real(1:oldsize)     = edin%xatm(1:oldsize)
!    deallocate(edin%xatm)
!    allocate(edin%xatm(newsize))
!    edin%xatm(1:newsize)       = edout_real(water_site_new(1:newsize)))
    
!    edout_real(1:oldsize)     = edin%yatm(1:oldsize)
!    deallocate(yatm%lat)
!    allocate(yatm%lat(newsize))
!    edin%yatm(1:newsize)       = edout_real(water_site_new(1:newsize)))
    
!    edout_real(1:oldsize)     = edin%natm(1:oldsize)
!    deallocate(edin%natm)
!    allocate(edin%natm(newsize))
!    edin%natm(1:newsize)       = edout_real(water_site_new(1:newsize)))
   

    return
  end subroutine copy_polygon
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_polygontype(polyin,polyout,isin,isout)
  ! =====================================================
  ! Copying Functions2
  ! =====================================================


    implicit none
    integer :: isin,isout
    type(polygontype),pointer :: polyin,polyout

    polyout%patch_count(isout)               = polyin%patch_count(isin)
    polyout%ntext_soil(1:nzg,isout)          = polyin%ntext_soil(1:nzg,isin)
    polyout%lsl(isout)                       = polyin%lsl(isin)
    polyout%area(isout)                      = polyin%area(isin)
    polyout%moist_f(isout)                   = polyin%moist_f(isin)
    polyout%sitenum(isout)                   = polyin%sitenum(isin)
    polyout%elevation(isout)                 = polyin%elevation(isin)
    polyout%slope(isout)                     = polyin%slope(isin)
    polyout%aspect(isout)                    = polyin%aspect(isin)
    polyout%TCI(isout)                       = polyin%TCI(isin)
    polyout%pptweight(isout)                 = polyin%pptweight(isin)


    return
  end subroutine copy_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_sitetype(sitein,siteout,ipin,ipout)

    implicit none
    integer :: ipin,ipout
    type(sitetype),pointer :: sitein,siteout

    print*,"COPY_SITETYPE_ IS NOT POPULATED YET"
    stop

    siteout%dist_type(ipout)          = sitein%dist_type(ipin)
    siteout%age(ipout)                = sitein%age(ipin)
    siteout%area(ipout)               = sitein%area(ipin)
    siteout%hcapveg(ipout)            = sitein%hcapveg(ipin)
    siteout%fast_soil_C(ipout)        = sitein%fast_soil_C(ipin)
    siteout%slow_soil_C(ipout)        = sitein%slow_soil_C(ipin)
    siteout%structural_soil_C(ipout)  = sitein%structural_soil_C(ipin)
    siteout%structural_soil_L(ipout)  = sitein%structural_soil_L(ipin)
    siteout%mineralized_soil_N(ipout) = sitein%mineralized_soil_N(ipin)
    siteout%fast_soil_N(ipout)        = sitein%fast_soil_N(ipin)
    siteout%pname(ipout)              = sitein%pname(ipin)
    siteout%sum_dgd(ipout)            = sitein%sum_dgd(ipin)
    siteout%sum_chd(ipout)            = sitein%sum_chd(ipin)
    siteout%plantation(ipout)         = sitein%plantation(ipin)
    siteout%cohort_count(ipout)       = sitein%cohort_count(ipin)
    siteout%can_enthalpy(ipout)       = sitein%can_enthalpy(ipin)
    siteout%can_temp(ipout)           = sitein%can_temp(ipin)
    siteout%can_shv(ipout)            = sitein%can_shv(ipin)
    siteout%can_co2(ipout)            = sitein%can_co2(ipin)
    siteout%can_rhos(ipout)           = sitein%can_rhos(ipin)
    siteout%can_prss(ipout)           = sitein%can_prss(ipin)
    siteout%can_theta(ipout)          = sitein%can_theta(ipin)
    siteout%can_depth(ipout)          = sitein%can_depth(ipin)
    siteout%lambda_light(ipout)        = sitein%lambda_light(ipin)
       
    if (idoutput > 0 .or. imoutput > 0) then
       siteout%dmean_rh             (ipout) = sitein%dmean_rh             (ipin)
       siteout%dmean_co2_residual   (ipout) = sitein%dmean_co2_residual   (ipin)
       siteout%dmean_energy_residual(ipout) = sitein%dmean_energy_residual(ipin)
       siteout%dmean_water_residual (ipout) = sitein%dmean_water_residual (ipin)
       siteout%dmean_lambda_light   (ipout) = sitein%dmean_lambda_light   (ipin)
    end if
    
    if (imoutput > 0) then
       siteout%mmean_rh             (ipout) = sitein%mmean_rh             (ipin)
       siteout%mmean_co2_residual   (ipout) = sitein%mmean_co2_residual   (ipin)
       siteout%mmean_energy_residual(ipout) = sitein%mmean_energy_residual(ipin)
       siteout%mmean_water_residual (ipout) = sitein%mmean_water_residual (ipin)
       siteout%mmean_lambda_light   (ipout) = sitein%mmean_lambda_light   (ipin)
    end if
  
    return
  end subroutine copy_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_sitetype_mask(sitein,siteout,logmask,masksz,newsz)

    ! This subroutine assumes that the size of vectors in siteout
    ! are the number of true elements in mask, while the size of the
    ! vectors in sitein are of the size of the mask itself. If this
    ! is not true, you will get a segmentation violation and the
    ! code will crash.args 1 and 3 must be dimension of arg 4
    ! argument 2 must be the dimension of the sum of the 3rd argument
    ! 
    ! THIS ROUTINE CURRENTLY ASSUMES THAT THE OUTPUT SITE
    ! HAS NOT ALLOCATED IT'S PATCH'S COHORT VECTORS YET, THIS
    ! IS BECAUSE THE LENGTHS OF THESE VECTORS ARE BASED ON THE
    ! DONOR PATH'S VECTOR SIZES.  DO NOT USE PRE-ALLOCATED
    ! RECIPIENTS

    implicit none

    type(sitetype),target :: sitein,siteout
    integer :: masksz,newsz
    integer,dimension(newsz) :: incmask
    logical,dimension(masksz)           :: logmask
    integer :: i,k,m,inc,ipft
    type(stoma_data),pointer :: osdi,osdo

    inc = 0
    do i=1,masksz
       if (logmask(i)) then
          inc = inc + 1
          incmask(inc) = i
       end if
    end do

    ! First do all of the true vectors
    siteout%paco_id(1:inc)              = pack(sitein%paco_id,logmask)
    siteout%paco_n(1:inc)               = pack(sitein%paco_n,logmask)
    siteout%dist_type(1:inc)            = pack(sitein%dist_type,logmask)
    siteout%age(1:inc)                  = pack(sitein%age,logmask)
    siteout%area(1:inc)                 = pack(sitein%area,logmask)
    siteout%hcapveg(1:inc)              = pack(sitein%hcapveg,logmask)
    siteout%fast_soil_C(1:inc)          = pack(sitein%fast_soil_C,logmask)
    siteout%slow_soil_C(1:inc)          = pack(sitein%slow_soil_C,logmask)
    siteout%structural_soil_C(1:inc)    = pack(sitein%structural_soil_C,logmask)
    siteout%structural_soil_L(1:inc)    = pack(sitein%structural_soil_L,logmask)
    siteout%mineralized_soil_N(1:inc)   = pack(sitein%mineralized_soil_N,logmask)
    siteout%fast_soil_N(1:inc)          = pack(sitein%fast_soil_N,logmask)
    siteout%sum_dgd(1:inc)              = pack(sitein%sum_dgd,logmask)
    siteout%sum_chd(1:inc)              = pack(sitein%sum_chd,logmask)
    siteout%plantation(1:inc)           = pack(sitein%plantation,logmask)
    siteout%cohort_count(1:inc)         = pack(sitein%cohort_count,logmask)
    siteout%can_enthalpy(1:inc)         = pack(sitein%can_enthalpy,logmask)
    siteout%can_temp(1:inc)             = pack(sitein%can_temp,logmask)
    siteout%can_shv(1:inc)              = pack(sitein%can_shv,logmask)
    siteout%can_co2(1:inc)              = pack(sitein%can_co2,logmask)
    siteout%can_rhos(1:inc)             = pack(sitein%can_rhos,logmask)
    siteout%can_prss(1:inc)             = pack(sitein%can_prss,logmask)
    siteout%can_theta(1:inc)            = pack(sitein%can_theta,logmask)
    siteout%can_depth(1:inc)            = pack(sitein%can_depth,logmask)
    siteout%lambda_light(1:inc)          = pack(sitein%lambda_light,logmask)
    siteout%lai(1:inc)                  = pack(sitein%lai,logmask)
    siteout%wpa(1:inc)                  = pack(sitein%wpa,logmask)
    siteout%wai(1:inc)                  = pack(sitein%wai,logmask)
    siteout%avg_daily_temp(1:inc)       = pack(sitein%avg_daily_temp,logmask)
    siteout%mean_rh(1:inc)              = pack(sitein%mean_rh,logmask)
    siteout%mean_nep(1:inc)             = pack(sitein%mean_nep,logmask)
    siteout%wbudget_loss2atm(1:inc)     = pack(sitein%wbudget_loss2atm,logmask)
    siteout%wbudget_denseffect(1:inc)        = pack(sitein%wbudget_denseffect,logmask)
    siteout%wbudget_precipgain(1:inc)        = pack(sitein%wbudget_precipgain,logmask)
    siteout%wbudget_loss2runoff(1:inc)       = pack(sitein%wbudget_loss2runoff,logmask)
    siteout%wbudget_loss2drainage(1:inc)     = pack(sitein%wbudget_loss2drainage,logmask)
    siteout%wbudget_initialstorage(1:inc)    = pack(sitein%wbudget_initialstorage,logmask)
    siteout%wbudget_residual(1:inc)          = pack(sitein%wbudget_residual,logmask)
    siteout%ebudget_latent(1:inc)            = pack(sitein%ebudget_latent,logmask)
    siteout%ebudget_loss2atm(1:inc)          = pack(sitein%ebudget_loss2atm,logmask)
    siteout%ebudget_denseffect(1:inc)        = pack(sitein%ebudget_denseffect,logmask)
    siteout%ebudget_loss2runoff(1:inc)       = pack(sitein%ebudget_loss2runoff,logmask)
    siteout%ebudget_loss2drainage(1:inc)     = pack(sitein%ebudget_loss2drainage,logmask)
    siteout%ebudget_netrad(1:inc)            = pack(sitein%ebudget_netrad,logmask)
    siteout%ebudget_precipgain(1:inc)        = pack(sitein%ebudget_precipgain,logmask)
    siteout%ebudget_initialstorage(1:inc)    = pack(sitein%ebudget_initialstorage,logmask)
    siteout%ebudget_residual(1:inc)          = pack(sitein%ebudget_residual,logmask)
    siteout%co2budget_initialstorage(1:inc)  = pack(sitein%co2budget_initialstorage,logmask)
    siteout%co2budget_residual(1:inc)   = pack(sitein%co2budget_residual,logmask)
    siteout%co2budget_loss2atm(1:inc)   = pack(sitein%co2budget_loss2atm,logmask)
    siteout%co2budget_denseffect(1:inc) = pack(sitein%co2budget_denseffect,logmask)
    siteout%co2budget_gpp(1:inc)        = pack(sitein%co2budget_gpp,logmask)
    siteout%co2budget_plresp(1:inc)     = pack(sitein%co2budget_plresp,logmask)
    siteout%co2budget_rh(1:inc)         = pack(sitein%co2budget_rh,logmask)
    siteout%today_A_decomp(1:inc)       = pack(sitein%today_A_decomp,logmask)
    siteout%today_Af_decomp(1:inc)      = pack(sitein%today_Af_decomp,logmask)
    siteout%veg_rough(1:inc)            = pack(sitein%veg_rough,logmask)
    siteout%veg_height(1:inc)           = pack(sitein%veg_height,logmask)
    siteout%fsc_in(1:inc)               = pack(sitein%fsc_in,logmask)
    siteout%ssc_in(1:inc)               = pack(sitein%ssc_in,logmask)
    siteout%ssl_in(1:inc)               = pack(sitein%ssl_in,logmask)
    siteout%fsn_in(1:inc)               = pack(sitein%fsn_in,logmask)
    siteout%total_plant_nitrogen_uptake(1:inc)    = pack(sitein%total_plant_nitrogen_uptake,logmask)
    siteout%mineralized_N_loss(1:inc)    = pack(sitein%mineralized_N_loss,logmask)
    siteout%mineralized_N_input(1:inc)    = pack(sitein%mineralized_N_input,logmask)
    siteout%rshort_g(1:inc)             = pack(sitein%rshort_g,logmask)
    siteout%rshort_g_beam(1:inc)        = pack(sitein%rshort_g_beam,logmask)
    siteout%rshort_g_diffuse(1:inc)     = pack(sitein%rshort_g_diffuse,logmask)
    siteout%rlong_g(1:inc)              = pack(sitein%rlong_g,logmask)
    siteout%rlong_g_surf(1:inc)         = pack(sitein%rlong_g_surf,logmask)
    siteout%rlong_g_incid(1:inc)        = pack(sitein%rlong_g_incid,logmask)
    siteout%rlong_s(1:inc)              = pack(sitein%rlong_s,logmask)
    siteout%rlong_s_surf(1:inc)         = pack(sitein%rlong_s_surf,logmask)
    siteout%rlong_s_incid(1:inc)        = pack(sitein%rlong_s_incid,logmask)
    siteout%albedt(1:inc)               = pack(sitein%albedt,logmask)
    siteout%albedo_beam(1:inc)          = pack(sitein%albedo_beam,logmask)
    siteout%albedo_diffuse(1:inc)       = pack(sitein%albedo_diffuse,logmask)
    siteout%rlongup(1:inc)              = pack(sitein%rlongup,logmask)
    siteout%rlong_albedo(1:inc)         = pack(sitein%rlong_albedo,logmask)
    siteout%total_snow_depth(1:inc)     = pack(sitein%total_snow_depth,logmask)
    siteout%snowfac(1:inc)              = pack(sitein%snowfac,logmask)
    siteout%A_decomp(1:inc)             = pack(sitein%A_decomp,logmask)
    siteout%f_decomp(1:inc)             = pack(sitein%f_decomp,logmask)
    siteout%rh(1:inc)                   = pack(sitein%rh,logmask)
    siteout%cwd_rh(1:inc)               = pack(sitein%cwd_rh,logmask)
    siteout%fuse_flag(1:inc)            = pack(sitein%fuse_flag,logmask)
    siteout%plant_ag_biomass(1:inc)     = pack(sitein%plant_ag_biomass,logmask)
    siteout%mean_wflux(1:inc)           = pack(sitein%mean_wflux,logmask)
    siteout%mean_latflux(1:inc)         = pack(sitein%mean_latflux,logmask)
    siteout%mean_hflux(1:inc)           = pack(sitein%mean_hflux,logmask)
    siteout%mean_runoff(1:inc)          = pack(sitein%mean_runoff,logmask)
    siteout%mean_qrunoff(1:inc)         = pack(sitein%mean_qrunoff,logmask)
    siteout%htry(1:inc)                 = pack(sitein%htry,logmask)
    siteout%avg_rk4step(1:inc)          = pack(sitein%avg_rk4step,logmask)
    siteout%avg_available_water(1:inc)  = pack(sitein%avg_available_water,logmask)
    siteout%ustar(1:inc)                = pack(sitein%ustar,logmask)
    siteout%tstar(1:inc)                = pack(sitein%tstar,logmask)
    siteout%qstar(1:inc)                = pack(sitein%qstar,logmask)
    siteout%cstar(1:inc)                = pack(sitein%cstar,logmask)
    
    siteout%upwp(1:inc)                 = pack(sitein%upwp,logmask)
    siteout%tpwp(1:inc)                 = pack(sitein%tpwp,logmask)
    siteout%qpwp(1:inc)                 = pack(sitein%qpwp,logmask)
    siteout%cpwp(1:inc)                 = pack(sitein%cpwp,logmask)
    siteout%wpwp(1:inc)                 = pack(sitein%wpwp,logmask)

    siteout%avg_carbon_ac(1:inc)        = pack(sitein%avg_carbon_ac,logmask)
    siteout%nlev_sfcwater(1:inc)        = pack(sitein%nlev_sfcwater,logmask)
    siteout%ground_shv(1:inc)           = pack(sitein%ground_shv,logmask)
    siteout%surface_ssh(1:inc)          = pack(sitein%surface_ssh,logmask)
    siteout%rough(1:inc)                = pack(sitein%rough,logmask)

    siteout%avg_carbon_ac(1:inc)        = pack(sitein%avg_carbon_ac,logmask)
    siteout%avg_vapor_vc(1:inc)         = pack(sitein%avg_vapor_vc,logmask)
    siteout%avg_dew_cg(1:inc)           = pack(sitein%avg_dew_cg,logmask)
    siteout%avg_vapor_gc(1:inc)         = pack(sitein%avg_vapor_gc,logmask)
    siteout%avg_wshed_vg(1:inc)         = pack(sitein%avg_wshed_vg,logmask)
    siteout%avg_intercepted(1:inc)      = pack(sitein%avg_intercepted,logmask)
    siteout%avg_vapor_ac(1:inc)         = pack(sitein%avg_vapor_ac,logmask)
    siteout%avg_transp(1:inc)           = pack(sitein%avg_transp,logmask)
    siteout%avg_evap(1:inc)             = pack(sitein%avg_evap,logmask)
    siteout%avg_netrad(1:inc)           = pack(sitein%avg_netrad,logmask)
    siteout%avg_runoff(1:inc)           = pack(sitein%avg_runoff,logmask)
    siteout%avg_drainage(1:inc)         = pack(sitein%avg_drainage,logmask)
    siteout%avg_drainage_heat(1:inc)    = pack(sitein%avg_drainage_heat,logmask)
    siteout%aux(1:inc)                  = pack(sitein%aux,logmask)
    siteout%avg_sensible_vc(1:inc)      = pack(sitein%avg_sensible_vc,logmask)
    siteout%avg_qwshed_vg(1:inc)        = pack(sitein%avg_qwshed_vg,logmask)
    siteout%avg_qintercepted(1:inc)     = pack(sitein%avg_qintercepted,logmask)
    siteout%avg_sensible_gc(1:inc)      = pack(sitein%avg_sensible_gc,logmask)
    siteout%avg_sensible_ac(1:inc)      = pack(sitein%avg_sensible_ac,logmask)
    siteout%avg_runoff_heat(1:inc)      = pack(sitein%avg_runoff_heat,logmask)
    siteout%avg_veg_energy(1:inc)       = pack(sitein%avg_veg_energy,logmask)
    siteout%avg_veg_temp(1:inc)         = pack(sitein%avg_veg_temp,logmask)
    siteout%avg_veg_fliq(1:inc)         = pack(sitein%avg_veg_fliq,logmask)
    siteout%avg_veg_water(1:inc)        = pack(sitein%avg_veg_water,logmask)
  


    ! Water layers 1:nzs
    
    do k=1,nzs
       siteout%sfcwater_mass(k,1:inc)    = pack(sitein%sfcwater_mass(k,:),logmask)
       siteout%sfcwater_energy(k,1:inc)  = pack(sitein%sfcwater_energy(k,:),logmask)
       siteout%sfcwater_depth(k,1:inc)   = pack(sitein%sfcwater_depth(k,:),logmask)
       siteout%rshort_s(k,1:inc)         = pack(sitein%rshort_s(k,:),logmask)
       siteout%rshort_s_beam(k,1:inc)    = pack(sitein%rshort_s_beam(k,:),logmask)
       siteout%rshort_s_diffuse(k,1:inc) = pack(sitein%rshort_s_diffuse(k,:),logmask)
       siteout%sfcwater_tempk(k,1:inc)   = pack(sitein%sfcwater_tempk(k,:),logmask)
       siteout%sfcwater_fracliq(k,1:inc)  = pack(sitein%sfcwater_fracliq(k,:),logmask)
    enddo

    ! Soil layers 1:nzg

    do k=1,nzg
       siteout%ntext_soil(k,1:inc)         = pack(sitein%ntext_soil(k,:),logmask)
       siteout%soil_energy(k,1:inc)        = pack(sitein%soil_energy(k,:),logmask)
       siteout%soil_water(k,1:inc)         = pack(sitein%soil_water(k,:),logmask)
       siteout%soil_tempk(k,1:inc)         = pack(sitein%soil_tempk(k,:),logmask)
       siteout%soil_fracliq(k,1:inc)       = pack(sitein%soil_fracliq(k,:),logmask)
       siteout%avg_smoist_gg(k,1:inc)      = pack(sitein%avg_smoist_gg(k,:),logmask)
       siteout%avg_smoist_gc(k,1:inc)      = pack(sitein%avg_smoist_gc(k,:),logmask)
       siteout%avg_sensible_gg(k,1:inc)    = pack(sitein%avg_sensible_gg(k,:),logmask)
       siteout%aux_s(k,1:inc)              = pack(sitein%aux_s(k,:),logmask)
    enddo

    ! pft types

    do k=1,n_pft
       siteout%repro(k,1:inc)            = pack(sitein%repro(k,:),logmask)
       siteout%A_o_max(k,1:inc)          = pack(sitein%A_o_max(k,:),logmask)
       siteout%A_c_max(k,1:inc)          = pack(sitein%A_c_max(k,:),logmask)
       
       do m=1,n_stoma_atts
          siteout%old_stoma_vector_max(m,k,1:inc) = pack(sitein%old_stoma_vector_max(m,k,:),logmask)
       end do

       do m=1,ff_ndbh
          siteout%pft_density_profile(k,m,1:inc) = pack(sitein%pft_density_profile(k,m,:),logmask)
       enddo
    enddo
    
    !dbh types
    do k=1,n_dbh
        siteout%co2budget_gpp_dbh(k,1:inc)        = pack(sitein%co2budget_gpp_dbh(k,:),logmask)
    end do

    
    ! Old_stoma_data_max type
    ! Derived type with n_pft x n_patch, so.... the intrinsic "pack" wont work

    do m=1,newsz
       k=incmask(m)
       call allocate_patchtype(siteout%patch(m),sitein%patch(k)%ncohorts)
       
       call copy_patchtype(sitein%patch(k),siteout%patch(m),1,sitein%patch(k)%ncohorts,1,sitein%patch(k)%ncohorts)

       do ipft=1,n_pft

          osdo => siteout%old_stoma_data_max(ipft,m)
          osdi => sitein%old_stoma_data_max(ipft,k)
          
          osdo%recalc           = osdi%recalc
          osdo%T_L              = osdi%T_L
          osdo%e_A              = osdi%e_A
          osdo%PAR              = osdi%PAR
          osdo%rb_factor        = osdi%rb_factor
          osdo%prss             = osdi%prss
          osdo%phenology_factor = osdi%phenology_factor
          osdo%gsw_open         = osdi%gsw_open
          osdo%ilimit           = osdi%ilimit
          osdo%T_L_residual     = osdi%T_L_residual
          osdo%e_a_residual     = osdi%e_a_residual
          osdo%par_residual     = osdi%par_residual
          osdo%rb_residual      = osdi%rb_residual
          osdo%leaf_residual    = osdi%leaf_residual
          osdo%gsw_residual     = osdi%gsw_residual

       enddo
       
    enddo
       
    if (idoutput > 0 .or. imoutput > 0) then
       siteout%dmean_rh             (1:inc) = pack(sitein%dmean_rh             ,logmask)
       siteout%dmean_co2_residual   (1:inc) = pack(sitein%dmean_co2_residual   ,logmask)
       siteout%dmean_energy_residual(1:inc) = pack(sitein%dmean_energy_residual,logmask)
       siteout%dmean_water_residual (1:inc) = pack(sitein%dmean_water_residual ,logmask)
       siteout%dmean_lambda_light   (1:inc) = pack(sitein%dmean_lambda_light   ,logmask)
       siteout%dmean_A_decomp       (1:inc) = pack(sitein%dmean_A_decomp       ,logmask)
       siteout%dmean_Af_decomp      (1:inc) = pack(sitein%dmean_Af_decomp      ,logmask)
       siteout%dmean_rk4step        (1:inc) = pack(sitein%dmean_rk4step        ,logmask)
    end if
    
    if (imoutput > 0) then
       siteout%mmean_rh             (1:inc) = pack(sitein%mmean_rh             ,logmask)
       siteout%mmean_co2_residual   (1:inc) = pack(sitein%mmean_co2_residual   ,logmask)
       siteout%mmean_energy_residual(1:inc) = pack(sitein%mmean_energy_residual,logmask)
       siteout%mmean_water_residual (1:inc) = pack(sitein%mmean_water_residual ,logmask)
       siteout%mmean_lambda_light   (1:inc) = pack(sitein%mmean_lambda_light   ,logmask)
       siteout%mmean_A_decomp       (1:inc) = pack(sitein%mmean_A_decomp       ,logmask)
       siteout%mmean_Af_decomp      (1:inc) = pack(sitein%mmean_Af_decomp      ,logmask)
       siteout%mmean_rk4step        (1:inc) = pack(sitein%mmean_rk4step        ,logmask)
    end if

    return
  end subroutine copy_sitetype_mask
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_patchtype_mask(patchin,patchout,mask,masksz,newsz)

    ! This subroutine assumes that the size of vectors in siteout
    ! are the number of true elements in mask, while the size of the
    ! vectors in sitein are of the size of the mask itself. If this
    ! is not true, you will get a segmentation violation and the
    ! code will crash.args 1 and 3 must be dimension of arg 4
    ! argument 2 must be the dimension of the sum of the 3rd argument
    ! 
    ! THIS ROUTINE CURRENTLY ASSUMES THAT THE OUTPUT SITE
    ! HAS NOT ALLOCATED IT'S PATCH'S COHORT VECTORS YET, THIS
    ! IS BECAUSE THE LENGTHS OF THESE VECTORS ARE BASED ON THE
    ! DONOR PATH'S VECTOR SIZES.  DO NOT USE PRE-ALLOCATED
    ! RECIPIENTS
    
    implicit none

    type(patchtype),target :: patchin,patchout
    integer :: masksz,newsz
    integer,dimension(newsz) :: incmask
    integer,dimension(masksz):: imask
    logical,dimension(masksz)  :: mask
    integer :: i,k,m,inc
    type(stoma_data),pointer :: osdi,osdo

    do i=1,masksz
       imask(i) = i
    enddo
    inc=count(mask)                   ! Number of true elements

    if (inc == 0) return

    incmask=pack(imask,mask)   ! List of true elements
    
    patchout%pft(1:inc)              = pack(patchin%pft,mask)
    patchout%nplant(1:inc)           = pack(patchin%nplant,mask)
    patchout%hite(1:inc)             = pack(patchin%hite,mask)
    patchout%agb(1:inc)              = pack(patchin%agb,mask)
    patchout%basarea(1:inc)          = pack(patchin%basarea,mask)
    patchout%dagb_dt(1:inc)          = pack(patchin%dagb_dt,mask)
    patchout%dba_dt(1:inc)           = pack(patchin%dba_dt,mask)
    patchout%ddbh_dt(1:inc)          = pack(patchin%ddbh_dt,mask)
    patchout%dbh(1:inc)              = pack(patchin%dbh,mask)
    patchout%bdead(1:inc)            = pack(patchin%bdead,mask)
    patchout%bleaf(1:inc)            = pack(patchin%bleaf,mask)
    patchout%phenology_status(1:inc) = pack(patchin%phenology_status,mask)
    patchout%balive(1:inc)           = pack(patchin%balive,mask)
    patchout%broot(1:inc)            = pack(patchin%broot,mask)
    patchout%bsapwood(1:inc)         = pack(patchin%bsapwood,mask)
    patchout%lai(1:inc)              = pack(patchin%lai,mask)
    patchout%wpa(1:inc)              = pack(patchin%wpa,mask)
    patchout%wai(1:inc)              = pack(patchin%wai,mask)
    patchout%solvable(1:inc)         = pack(patchin%solvable,mask)
    patchout%bstorage(1:inc)         = pack(patchin%bstorage,mask)
    patchout%cbr_bar(1:inc)          = pack(patchin%cbr_bar,mask)
    patchout%veg_energy(1:inc)       = pack(patchin%veg_energy,mask)
    patchout%veg_temp(1:inc)         = pack(patchin%veg_temp,mask)
    patchout%veg_fliq(1:inc)         = pack(patchin%veg_fliq,mask)
    patchout%veg_water(1:inc)        = pack(patchin%veg_water,mask)
    patchout%mean_gpp(1:inc)         = pack(patchin%mean_gpp,mask)
    patchout%mean_leaf_resp(1:inc)   = pack(patchin%mean_leaf_resp,mask)
    patchout%mean_root_resp(1:inc)   = pack(patchin%mean_root_resp,mask)
    patchout%mean_growth_resp(1:inc) = pack(patchin%mean_growth_resp,mask)
    patchout%mean_storage_resp(1:inc)= pack(patchin%mean_storage_resp,mask)
    patchout%mean_vleaf_resp(1:inc)  = pack(patchin%mean_vleaf_resp,mask)
    patchout%today_leaf_resp(1:inc)  = pack(patchin%today_leaf_resp,mask)
    patchout%today_root_resp(1:inc)  = pack(patchin%today_root_resp,mask)
    patchout%today_gpp(1:inc)        = pack(patchin%today_gpp,mask)
    patchout%today_gpp_pot(1:inc)    = pack(patchin%today_gpp_pot,mask)
    patchout%today_gpp_max(1:inc)    = pack(patchin%today_gpp_max,mask)
    patchout%growth_respiration(1:inc) = pack(patchin%growth_respiration,mask)
    patchout%storage_respiration(1:inc) = pack(patchin%storage_respiration,mask)
    patchout%vleaf_respiration(1:inc) = pack(patchin%vleaf_respiration,mask)
    patchout%fsn(1:inc)              = pack(patchin%fsn,mask)
    patchout%monthly_dndt(1:inc)     = pack(patchin%monthly_dndt,mask)
    
    patchout%Psi_open(1:inc)         = pack(patchin%Psi_open,mask)
    patchout%krdepth(1:inc)          = pack(patchin%krdepth,mask)
    patchout%first_census(1:inc)     = pack(patchin%first_census,mask)
    patchout%new_recruit_flag(1:inc) = pack(patchin%new_recruit_flag,mask)
    patchout%light_level(1:inc)      = pack(patchin%light_level,mask)
    patchout%light_level_beam(1:inc) = pack(patchin%light_level_beam,mask)
    patchout%light_level_diff(1:inc) = pack(patchin%light_level_diff,mask)
    patchout%beamext_level(1:inc)      = pack(patchin%beamext_level,mask)
    patchout%diffext_level(1:inc)      = pack(patchin%diffext_level,mask)
    patchout%norm_par_beam(1:inc)    = pack(patchin%norm_par_beam,mask)
    patchout%norm_par_diff(1:inc)    = pack(patchin%norm_par_diff,mask)
    patchout%lambda_light(1:inc)     = pack(patchin%lambda_light,mask)
    patchout%par_v(1:inc)            = pack(patchin%par_v,mask)
    patchout%par_v_beam(1:inc)       = pack(patchin%par_v_beam,mask)
    patchout%par_v_diffuse(1:inc)    = pack(patchin%par_v_diffuse,mask)
    patchout%rshort_v(1:inc)         = pack(patchin%rshort_v,mask)
    patchout%rshort_v_beam(1:inc)    = pack(patchin%rshort_v_beam,mask)
    patchout%rshort_v_diffuse(1:inc) = pack(patchin%rshort_v_diffuse,mask)
    patchout%rlong_v(1:inc)          = pack(patchin%rlong_v,mask)
    patchout%rlong_v_surf(1:inc)     = pack(patchin%rlong_v_surf,mask)
    patchout%rlong_v_incid(1:inc)    = pack(patchin%rlong_v_incid,mask)
    patchout%rb(1:inc)               = pack(patchin%rb,mask)
    patchout%A_open(1:inc)           = pack(patchin%A_open,mask)
    patchout%A_closed(1:inc)         = pack(patchin%A_closed,mask)
    patchout%Psi_closed(1:inc)       = pack(patchin%Psi_closed,mask)
    patchout%rsw_open(1:inc)         = pack(patchin%rsw_open,mask)
    patchout%rsw_closed(1:inc)       = pack(patchin%rsw_closed,mask)
    patchout%fsw(1:inc)              = pack(patchin%fsw,mask)
    patchout%fs_open(1:inc)          = pack(patchin%fs_open,mask)
    patchout%stomatal_resistance(1:inc) = pack(patchin%stomatal_resistance,mask)
    patchout%leaf_maintenance(1:inc) = pack(patchin%leaf_maintenance,mask)
    patchout%root_maintenance(1:inc) = pack(patchin%root_maintenance,mask)
    patchout%leaf_drop(1:inc)        = pack(patchin%leaf_drop,mask)
    patchout%bseeds(1:inc)           = pack(patchin%bseeds,mask)
    patchout%leaf_respiration(1:inc) = pack(patchin%leaf_respiration,mask)
    patchout%root_respiration(1:inc) = pack(patchin%root_respiration,mask)
    patchout%hcapveg(1:inc)          = pack(patchin%hcapveg,mask)
    patchout%gpp(1:inc)              = pack(patchin%gpp,mask)
    patchout%paw_avg(1:inc)          = pack(patchin%paw_avg,mask)
    patchout%turnover_amp(1:inc)     = pack(patchin%turnover_amp,mask)
    patchout%llspan(1:inc)           = pack(patchin%llspan,mask)
    patchout%vm_bar(1:inc)           = pack(patchin%vm_bar,mask)
    patchout%sla(1:inc)              = pack(patchin%sla,mask)    

    do m=1,inc
       k=incmask(m)
       do i = 1,13
          patchout%cb(i,m)               = patchin%cb(i,k)
          patchout%cb_max(i,m)           = patchin%cb_max(i,k)
       end do
       do i = 1,n_stoma_atts
          patchout%old_stoma_vector(i,m) = patchin%old_stoma_vector(i,k)
       enddo
       do i = 1,n_mort
          patchout%mort_rate(i,m)       = patchin%mort_rate(i,k)
       end do
    end do

    
    
    ! Copy the stoma data
    do m=1,inc
       k=incmask(m)
       
       osdo => patchout%old_stoma_data(m)
       osdi => patchin%old_stoma_data(k)

       osdo%recalc           = osdi%recalc
       osdo%T_L              = osdi%T_L
       osdo%e_A              = osdi%e_A
       osdo%PAR              = osdi%PAR
       osdo%rb_factor        = osdi%rb_factor
       osdo%prss             = osdi%prss
       osdo%phenology_factor = osdi%phenology_factor
       osdo%gsw_open         = osdi%gsw_open
       osdo%ilimit           = osdi%ilimit
       osdo%T_L_residual     = osdi%T_L_residual
       osdo%e_a_residual     = osdi%e_a_residual
       osdo%par_residual     = osdi%par_residual
       osdo%rb_residual      = osdi%rb_residual
       osdo%leaf_residual    = osdi%leaf_residual
       osdo%gsw_residual     = osdi%gsw_residual
       
    enddo
    
    if (idoutput > 0 .or. imoutput > 0) then
       patchout%dmean_fs_open         (1:inc) = pack(patchin%dmean_fs_open         ,mask)
       patchout%dmean_fsw             (1:inc) = pack(patchin%dmean_fsw             ,mask)
       patchout%dmean_fsn             (1:inc) = pack(patchin%dmean_fsn             ,mask)
       patchout%dmean_lambda_light    (1:inc) = pack(patchin%dmean_lambda_light    ,mask)
       patchout%dmean_light_level     (1:inc) = pack(patchin%dmean_light_level     ,mask)
       patchout%dmean_light_level_beam(1:inc) = pack(patchin%dmean_light_level_beam,mask)
       patchout%dmean_light_level_diff(1:inc) = pack(patchin%dmean_light_level_diff,mask)
       patchout%dmean_beamext_level   (1:inc) = pack(patchin%dmean_beamext_level   ,mask)
       patchout%dmean_diffext_level   (1:inc) = pack(patchin%dmean_diffext_level   ,mask)
       patchout%dmean_norm_par_beam   (1:inc) = pack(patchin%dmean_norm_par_beam   ,mask)
       patchout%dmean_norm_par_diff   (1:inc) = pack(patchin%dmean_norm_par_diff   ,mask)
       patchout%dmean_gpp             (1:inc) = pack(patchin%dmean_gpp             ,mask)
       patchout%dmean_leaf_resp       (1:inc) = pack(patchin%dmean_leaf_resp       ,mask)
       patchout%dmean_root_resp       (1:inc) = pack(patchin%dmean_root_resp       ,mask)
       patchout%dmean_par_v           (1:inc) = pack(patchin%dmean_par_v           ,mask)
       patchout%dmean_par_v_beam      (1:inc) = pack(patchin%dmean_par_v_beam      ,mask)
       patchout%dmean_par_v_diff      (1:inc) = pack(patchin%dmean_par_v_diff      ,mask)
    end if
    if (imoutput > 0) then
       patchout%mmean_fs_open         (1:inc) = pack(patchin%mmean_fs_open         ,mask)
       patchout%mmean_fsw             (1:inc) = pack(patchin%mmean_fsw             ,mask)
       patchout%mmean_fsn             (1:inc) = pack(patchin%mmean_fsn             ,mask)
       patchout%mmean_leaf_maintenance(1:inc) = pack(patchin%mmean_leaf_maintenance,mask)
       patchout%mmean_root_maintenance(1:inc) = pack(patchin%mmean_root_maintenance,mask)
       patchout%mmean_leaf_drop       (1:inc) = pack(patchin%mmean_leaf_drop       ,mask)
       patchout%mmean_cb              (1:inc) = pack(patchin%mmean_cb              ,mask)
       patchout%mmean_lambda_light    (1:inc) = pack(patchin%mmean_lambda_light    ,mask)
       patchout%mmean_light_level     (1:inc) = pack(patchin%mmean_light_level     ,mask)
       patchout%mmean_light_level_beam(1:inc) = pack(patchin%mmean_light_level_beam,mask)
       patchout%mmean_light_level_diff(1:inc) = pack(patchin%mmean_light_level_diff,mask)
       patchout%mmean_beamext_level   (1:inc) = pack(patchin%mmean_beamext_level   ,mask)
       patchout%mmean_diffext_level   (1:inc) = pack(patchin%mmean_diffext_level   ,mask)
       patchout%mmean_norm_par_beam   (1:inc) = pack(patchin%mmean_norm_par_beam   ,mask)
       patchout%mmean_norm_par_diff   (1:inc) = pack(patchin%mmean_norm_par_diff   ,mask)
       patchout%mmean_gpp             (1:inc) = pack(patchin%mmean_gpp             ,mask)
       patchout%mmean_leaf_resp       (1:inc) = pack(patchin%mmean_leaf_resp       ,mask)
       patchout%mmean_root_resp       (1:inc) = pack(patchin%mmean_root_resp       ,mask)
       patchout%mmean_growth_resp     (1:inc) = pack(patchin%mmean_growth_resp     ,mask)
       patchout%mmean_storage_resp    (1:inc) = pack(patchin%mmean_storage_resp    ,mask)
       patchout%mmean_vleaf_resp      (1:inc) = pack(patchin%mmean_vleaf_resp      ,mask)
       patchout%mmean_par_v           (1:inc) = pack(patchin%mmean_par_v           ,mask)
       patchout%mmean_par_v_beam      (1:inc) = pack(patchin%mmean_par_v_beam      ,mask)
       patchout%mmean_par_v_diff      (1:inc) = pack(patchin%mmean_par_v_diff      ,mask)

       do m=1,inc
         k=incmask(m)
         do i = 1,n_mort
             patchout%mmean_mort_rate(i,m) = patchin%mmean_mort_rate(i,k)
          end do
       end do
    end if
    
    return
  end subroutine copy_patchtype_mask
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_patchtype(patchin,patchout,ipin1,ipin2,ipout1,ipout2)

    implicit none
    integer :: ipin1,ipin2,ipout1,ipout2
    type(patchtype),target :: patchin,patchout
    type(stoma_data),pointer :: osdo,osdi
    integer :: iout,iin

    if (ipout2-ipout1.ne.ipin2-ipin1) then
       print*,"In copy_patchtype:"
       print*,"You specified unequal vector lengths"
       print*,"in the input and output targets"
       print*,"This cannot be..stopping"
       call fatal_error('unequal vector lengths','copy_patchtype','ed_state_vars.f90')
    endif

    ! Copy the stoma data. Added the loop back here because sometimes ipin1 < ipin2
    ! for example, when ncohorts=0

    iin = ipin1
    do iout=ipout1,ipout2

       patchout%pft(iout)              = patchin%pft(iin)
       patchout%nplant(iout)           = patchin%nplant(iin)
       patchout%hite(iout)             = patchin%hite(iin)
       patchout%agb(iout)              = patchin%agb(iin)
       patchout%basarea(iout)          = patchin%basarea(iin)
       patchout%dagb_dt(iout)          = patchin%dagb_dt(iin)
       patchout%dba_dt(iout)           = patchin%dba_dt(iin)
       patchout%ddbh_dt(iout)          = patchin%ddbh_dt(iin)
       patchout%dbh(iout)              = patchin%dbh(iin)
       patchout%bdead(iout)            = patchin%bdead(iin)
       patchout%bleaf(iout)            = patchin%bleaf(iin)
       patchout%phenology_status(iout) = patchin%phenology_status(iin)
       patchout%balive(iout)           = patchin%balive(iin)
       patchout%broot(iout)            = patchin%broot(iin)
       patchout%bsapwood(iout)         = patchin%bsapwood(iin)
       patchout%lai(iout)              = patchin%lai(iin)
       patchout%wpa(iout)              = patchin%wpa(iin)
       patchout%wai(iout)              = patchin%wai(iin)
       patchout%solvable(iout)         = patchin%solvable(iin)
       patchout%bstorage(iout)         = patchin%bstorage(iin)
       patchout%cb(:,iout)             = patchin%cb(:,iin)
       patchout%cb_max(:,iout)         = patchin%cb_max(:,iin)
       patchout%cbr_bar(iout)          = patchin%cbr_bar(iin)
       patchout%veg_energy(iout)       = patchin%veg_energy(iin)
       patchout%veg_temp(iout)         = patchin%veg_temp(iin)
       patchout%veg_fliq(iout)         = patchin%veg_fliq(iin)
       patchout%veg_water(iout)        = patchin%veg_water(iin)
       patchout%mean_gpp(iout)         = patchin%mean_gpp(iin)
       patchout%mean_leaf_resp(iout)   = patchin%mean_leaf_resp(iin)
       patchout%mean_root_resp(iout)   = patchin%mean_root_resp(iin)
       patchout%mean_growth_resp(iout) = patchin%mean_growth_resp(iin)
       patchout%mean_storage_resp(iout)= patchin%mean_storage_resp(iin)
       patchout%mean_vleaf_resp(iout)  = patchin%mean_vleaf_resp(iin)
       patchout%today_leaf_resp(iout)  = patchin%today_leaf_resp(iin)
       patchout%today_root_resp(iout)  = patchin%today_root_resp(iin)
       patchout%today_gpp(iout)        = patchin%today_gpp(iin)
       patchout%today_gpp_pot(iout)    = patchin%today_gpp_pot(iin)
       patchout%today_gpp_max(iout)    = patchin%today_gpp_max(iin)
       patchout%growth_respiration(iout) = patchin%growth_respiration(iin)
       patchout%storage_respiration(iout) = patchin%storage_respiration(iin)
       patchout%vleaf_respiration(iout) = patchin%vleaf_respiration(iin)
       patchout%fsn(iout)               = patchin%fsn(iin)
       patchout%monthly_dndt(iout)      = patchin%monthly_dndt(iin)
       patchout%mort_rate(:,iout)       = patchin%mort_rate(:,iin)
    
       patchout%Psi_open(iout)         = patchin%Psi_open(iin)
       patchout%krdepth(iout)          = patchin%krdepth(iin)
       patchout%first_census(iout)     = patchin%first_census(iin)
       patchout%new_recruit_flag(iout) = patchin%new_recruit_flag(iin)
       patchout%light_level(iout)      = patchin%light_level(iin)
       patchout%light_level_beam(iout) = patchin%light_level_beam(iin)
       patchout%light_level_diff(iout) = patchin%light_level_diff(iin)
       patchout%beamext_level(iout)    = patchin%beamext_level(iin)
       patchout%diffext_level(iout)    = patchin%diffext_level(iin)
       patchout%norm_par_beam(iout)    = patchin%norm_par_beam(iin)
       patchout%norm_par_diff(iout)    = patchin%norm_par_diff(iin)
       patchout%lambda_light(iout)     = patchin%lambda_light(iin)
       patchout%par_v(iout)            = patchin%par_v(iin)
       patchout%par_v_beam(iout)       = patchin%par_v_beam(iin)
       patchout%par_v_diffuse(iout)    = patchin%par_v_diffuse(iin)
       patchout%rshort_v(iout)         = patchin%rshort_v(iin)
       patchout%rshort_v_beam(iout)    = patchin%rshort_v_beam(iin)
       patchout%rshort_v_diffuse(iout) = patchin%rshort_v_diffuse(iin)
       patchout%rlong_v(iout)          = patchin%rlong_v(iin)
       patchout%rlong_v_surf(iout)     = patchin%rlong_v_surf(iin)
       patchout%rlong_v_incid(iout)    = patchin%rlong_v_incid(iin)
       patchout%rb(iout)               = patchin%rb(iin)
       patchout%A_open(iout)           = patchin%A_open(iin)
       patchout%A_closed(iout)         = patchin%A_closed(iin)
       patchout%Psi_closed(iout)       = patchin%Psi_closed(iin)
       patchout%rsw_open(iout)         = patchin%rsw_open(iin)
       patchout%rsw_closed(iout)       = patchin%rsw_closed(iin)
       patchout%fsw(iout)              = patchin%fsw(iin)
       patchout%fs_open(iout)          = patchin%fs_open(iin)
       patchout%stomatal_resistance(iout) = patchin%stomatal_resistance(iin)
       patchout%leaf_maintenance(iout) = patchin%leaf_maintenance(iin)
       patchout%root_maintenance(iout) = patchin%root_maintenance(iin)
       patchout%leaf_drop(iout)        = patchin%leaf_drop(iin)
       patchout%bseeds(iout)           = patchin%bseeds(iin)
       patchout%leaf_respiration(iout) = patchin%leaf_respiration(iin)
       patchout%root_respiration(iout) = patchin%root_respiration(iin)
       patchout%hcapveg(iout)          = patchin%hcapveg(iin)
       patchout%gpp(iout)              = patchin%gpp(iin)
       patchout%paw_avg(iout)          = patchin%paw_avg(iin)
       patchout%turnover_amp(iout)     = patchin%turnover_amp(iin)
       patchout%llspan(iout)           = patchin%llspan(iin)
       patchout%vm_bar(iout)           = patchin%vm_bar(iin)
       patchout%sla(iout)              = patchin%sla(iin)

       patchout%old_stoma_vector(:,iout) = patchin%old_stoma_vector(:,iin)
       
       osdo => patchout%old_stoma_data(iout)
       osdi => patchin%old_stoma_data(iin)

       osdo%recalc           = osdi%recalc
       osdo%T_L              = osdi%T_L
       osdo%e_A              = osdi%e_A
       osdo%PAR              = osdi%PAR
       osdo%rb_factor        = osdi%rb_factor
       osdo%prss             = osdi%prss
       osdo%phenology_factor = osdi%phenology_factor
       osdo%gsw_open         = osdi%gsw_open
       osdo%ilimit           = osdi%ilimit
       osdo%T_L_residual     = osdi%T_L_residual
       osdo%e_a_residual     = osdi%e_a_residual
       osdo%par_residual     = osdi%par_residual
       osdo%rb_residual      = osdi%rb_residual
       osdo%leaf_residual    = osdi%leaf_residual
       osdo%gsw_residual     = osdi%gsw_residual

       if (imoutput > 0 .or. idoutput > 0 ) then
          patchout%dmean_fs_open           (iout) = patchin%dmean_fs_open           (iin)
          patchout%dmean_fsw               (iout) = patchin%dmean_fsw               (iin)
          patchout%dmean_fsn               (iout) = patchin%dmean_fsn               (iin)
          patchout%dmean_light_level       (iout) = patchin%dmean_light_level       (iin)
          patchout%dmean_light_level_beam  (iout) = patchin%dmean_light_level_beam  (iin)
          patchout%dmean_light_level_diff  (iout) = patchin%dmean_light_level_diff  (iin)
          patchout%dmean_beamext_level     (iout) = patchin%dmean_beamext_level     (iin)
          patchout%dmean_diffext_level     (iout) = patchin%dmean_diffext_level     (iin)
          patchout%dmean_norm_par_beam     (iout) = patchin%dmean_norm_par_beam     (iin)
          patchout%dmean_norm_par_diff     (iout) = patchin%dmean_norm_par_diff     (iin)
          patchout%dmean_lambda_light      (iout) = patchin%dmean_lambda_light      (iin)
          patchout%dmean_gpp               (iout) = patchin%dmean_gpp               (iin)
          patchout%dmean_leaf_resp         (iout) = patchin%dmean_leaf_resp         (iin)
          patchout%dmean_root_resp         (iout) = patchin%dmean_root_resp         (iin)
          patchout%dmean_par_v             (iout) = patchin%dmean_par_v             (iin)
          patchout%dmean_par_v_beam        (iout) = patchin%dmean_par_v_beam        (iin)
          patchout%dmean_par_v_diff        (iout) = patchin%dmean_par_v_diff        (iin)
       end if
       if (imoutput > 0) then
          patchout%mmean_fs_open           (iout) = patchin%mmean_fs_open           (iin)
          patchout%mmean_fsw               (iout) = patchin%mmean_fsw               (iin)
          patchout%mmean_fsn               (iout) = patchin%mmean_fsn               (iin)
          patchout%mmean_leaf_maintenance  (iout) = patchin%mmean_leaf_maintenance  (iin)
          patchout%mmean_root_maintenance  (iout) = patchin%mmean_root_maintenance  (iin)
          patchout%mmean_leaf_drop         (iout) = patchin%mmean_leaf_drop         (iin)
          patchout%mmean_cb                (iout) = patchin%mmean_cb                (iin)
          patchout%mmean_light_level       (iout) = patchin%mmean_light_level       (iin)
          patchout%mmean_light_level_beam  (iout) = patchin%mmean_light_level_beam  (iin)
          patchout%mmean_light_level_diff  (iout) = patchin%mmean_light_level_diff  (iin)
          patchout%mmean_beamext_level     (iout) = patchin%mmean_beamext_level     (iin)
          patchout%mmean_diffext_level     (iout) = patchin%mmean_diffext_level     (iin)
          patchout%mmean_norm_par_beam     (iout) = patchin%mmean_norm_par_beam     (iin)
          patchout%mmean_norm_par_diff     (iout) = patchin%mmean_norm_par_diff     (iin)
          patchout%mmean_gpp               (iout) = patchin%mmean_gpp               (iin)
          patchout%mmean_leaf_resp         (iout) = patchin%mmean_leaf_resp         (iin)
          patchout%mmean_root_resp         (iout) = patchin%mmean_root_resp         (iin)
          patchout%mmean_growth_resp       (iout) = patchin%mmean_growth_resp       (iin)
          patchout%mmean_storage_resp      (iout) = patchin%mmean_storage_resp      (iin)
          patchout%mmean_vleaf_resp        (iout) = patchin%mmean_vleaf_resp        (iin)
          patchout%mmean_mort_rate       (:,iout) = patchin%mmean_mort_rate       (:,iin)
          patchout%mmean_lambda_light      (iout) = patchin%mmean_lambda_light      (iin)
          patchout%mmean_par_v             (iout) = patchin%mmean_par_v             (iin)
          patchout%mmean_par_v_beam        (iout) = patchin%mmean_par_v_beam        (iin)
          patchout%mmean_par_v_diff        (iout) = patchin%mmean_par_v_diff        (iin)
       end if

       iin = iin + 1

    end do

    return
  end subroutine copy_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  

  ! ===============================================================
  ! Define the vtables of the state/output variables
  !
  ! The various state scalars, vectors and arrays are
  ! now populate the vtable.  The vtable indexes the array
  ! gives it a name, records its dimensions, when it is to be
  ! used as output and how (averaging and such) and most importantly
  ! saves a pointer to its starting position.  If this routine is
  ! being called as a compute node in parallel, the first position
  ! is not necessarily the first position of the whole datavector,
  ! but will only be the first position of that nodes hyperslab chunk
  ! within the given continuous dataset.
  !
  ! The first number correspond to the data level:
  ! 1. Gridtype     (polygon level)
  ! 2. Polygontype  (site level)
  ! 3. Sitetype     (patch level)
  ! 4. Patchtype    (cohort level)
  ! 9. Scalar
  ! The other numbers correspond to the kind of dimension and variable.
  ! 0. Main vector ordinate only, integer.
  ! 1. Main vector ordinate only, real.
  ! 2. Soil layer
  ! 3. Surface water layer
  ! 4. PFT
  ! 5. Disturbance Type
  ! 6. DBH class
  ! 7. FF_DBH class
  ! 8. Month
  ! 9. 13 months
  !
  !
  ! Of these seven possible dimension (2-8), they may be used concurrently
  ! to partition the data into multi-dimensional spaces, but all seven
  ! will not be used simultaneously.  The highest ranks in use are 3.
  ! Each unique combination will have a call number associated with it.
  !
  ! 10    : rank 1 : polygon, integer
  ! 11    : rank 1 : polygon
  ! 12    : rank 2 : polygon, s-layer
  ! 13    : rank 2 : polygon, w-layer
  ! 14    : rank 2 : polygon, pft
  ! 14567 : rank 5 : polygon, pft, disturbance, dbh, age
  ! 146   : rank 3 : polygon, pft, dbh
  ! 15    : rank 2 : polygon, disturbance
  ! 155   : rank 3 : polygon, disturbance, disturbance
  ! 157   : rank 3 : polygon, disturbance, age
  ! 16    : rank 2 : polygon, dbh
  ! 17    : rank 2 : polygon, age
  ! 18    : rank 2 : polygon, mort
  ! 19    : rank 2 : polygon, month+1
  !
  ! 20    : rank 1 : site, integer
  ! 21    : rank 1 : site
  ! 22    : rank 2 : site, s-layer
  ! 23    : rank 2 : site, w-layer
  ! 24    : rank 2 : site, pft
  ! 246   : rank 3 : site, pft, dbh
  ! 25    : rank 2 : site, disturbance
  ! 255   : rank 3 : site, disturbance, disturbance
  ! 26    : rank 2 : site, dbh
  ! 27    : rank 2 : site, age
  ! 28    : rank 2 : site, mort
  ! 29    : rank 2 : site, month
  !
  ! 30    : rank 1 : patch, integer
  ! 31    : rank 1 : patch
  ! 32    : rank 2 : patch, s-layer
  ! 33    : rank 2 : patch, w-layer
  ! 34    : rank 2 : patch, pft
  ! 346   : rank 3 : patch, pft, ff_dbh
  ! 35    : rank 2 : patch, disturbance
  ! 36    : rank 2 : patch, dbh
  ! 37    : rank 2 : patch, age
  ! 38    : rank 2 : patch, mort
  !
  ! 41    : rank 1 : cohort
  ! 44    : rank 2 : cohort, pft
  ! 46    : rank 2 : cohort, dbh
  ! 47    : rank 2 : cohort, age
  ! 48    : rank 2 : cohort, mort
  ! 49    : rank 2 : cohort, month+1
  !
  ! 90    : rank 0 : integer scalar 
  !===================================================================
  
  subroutine filltab_alltypes

    ! =================================================
    !
    ! This subroutine is the main driver for filling
    ! filling the var_table of ED variables.  On a
    ! serial computing environment, this routine should be
    ! called near the end of the initialization process
    ! after the hierarchical tree structure has been
    ! trimmed via fusion/fission.  Similiarly, this
    ! routine should be called after any fusion/fission
    ! process, assuming that the major vtable structures
    ! have been deallocated prior to reallocation.
    !
    ! In a paralell environment, this routine should
    ! operate in a similiar fashion on each of the compute
    ! nodes.  It is designed such that the compute nodes
    ! will write hyperslabs of data in parallel to a
    ! joing HDF5 dataset as "collective-chunked" data
    ! The modifications that must be made after running this
    ! subroutine, are that the indexes should account
    ! for the offset of the current compute node.
    !
    ! =================================================
    
    
    use ed_var_tables,only:num_var,vt_info,var_table,nullify_vt_vector_pointers
    use ed_node_coms,only:mynum,mchnum,machs,nmachs,nnodetot,sendnum,recvnum,master_num
    use ed_max_dims, only: maxgrds, maxmach
    implicit none
    
    include 'mpif.h'

    integer :: ncohorts_g,npatches_g,nsites_g
    integer :: igr,ipy,isi,ipa,nv,ierr,nm,iptr
    integer,       dimension(MPI_STATUS_SIZE) :: status
    integer :: ping,uniqueid
    logical,save :: model_start = .true.
   
    type(edtype),pointer      :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer    :: csite
    type(patchtype),pointer   :: cpatch
    logical :: verbose = .false.    

    if (mynum == 1) then
       write(*,"(a)")'--- Re-hashing the IO pointer tables and mapping arrays'
    endif


    ! The first loop through populates the info tables

    do igr = 1,ngrids
       cgrid => edgrid_g(igr)
       
       if (num_var(igr)>0) then
          do nv=1,num_var(igr)
             if (associated(vt_info(nv,igr)%vt_vector)) then
                do iptr=1,vt_info(nv,igr)%nptrs
                   call nullify_vt_vector_pointers(vt_info(nv,igr)%vt_vector(iptr))
                end do
                deallocate(vt_info(nv,igr)%vt_vector)
             endif
          enddo
       endif

       num_var(igr) = 0

       cgrid%npolygons_global = cgrid%npolygons
       cgrid%nsites_global    = get_nsites(cgrid)
       cgrid%npatches_global  = get_npatches(cgrid)
       cgrid%ncohorts_global  = get_ncohorts(cgrid)
       
       cgrid%mach_cohort_offset_index = 0
       cgrid%mach_patch_offset_index  = 0
       cgrid%mach_site_offset_index   = 0
       cgrid%mach_polygon_offset_index= 0

       if (nnodetot /= 1) then

          ! Send all them sizes to root (CHANGED, NODE 1)
          

          if (mynum == 1) then
          
             gdpy(1,igr) = cgrid%npolygons_global
             gdsi(1,igr) = cgrid%nsites_global
             gdpa(1,igr) = cgrid%npatches_global
             gdco(1,igr) = cgrid%ncohorts_global

             call MPI_Send(ping,1,MPI_INTEGER,sendnum,94,MPI_COMM_WORLD,ierr)
             
             ! Have node 1 recieve the info
             do nm=2,nnodetot
                uniqueid=((igr-1)*maxmach)+nm
                call MPI_Recv(gdpy(nm,igr),1,MPI_INTEGER,machs(nm),500000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdsi(nm,igr),1,MPI_INTEGER,machs(nm),600000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdpa(nm,igr),1,MPI_INTEGER,machs(nm),700000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdco(nm,igr),1,MPI_INTEGER,machs(nm),800000+uniqueid,MPI_COMM_WORLD,status,ierr)
             enddo

             ! Broadcast all this info to the nodes
             do nm=2,nnodetot
                uniqueid=((igr-1)*maxmach)+nm
                call MPI_Send(gdpy,maxmach*maxgrds,MPI_INTEGER,machs(nm), 900000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdsi,maxmach*maxgrds,MPI_INTEGER,machs(nm),1000000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdpa,maxmach*maxgrds,MPI_INTEGER,machs(nm),1100000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdco,maxmach*maxgrds,MPI_INTEGER,machs(nm),1200000+uniqueid,MPI_COMM_WORLD,ierr)
             enddo

         else

            ! Set the blocking recieve to allow ordering, start with machine 1
            call MPI_Recv(ping,1,MPI_INTEGER,recvnum,94,MPI_COMM_WORLD,status,ierr)

            uniqueid=((igr-1)*maxmach)+mynum
            ! Send the information to node (1)
            call MPI_Send(cgrid%npolygons_global, 1,MPI_INTEGER,machs(1),500000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%nsites_global   , 1,MPI_INTEGER,machs(1),600000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%npatches_global , 1,MPI_INTEGER,machs(1),700000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%ncohorts_global , 1,MPI_INTEGER,machs(1),800000+uniqueid,MPI_COMM_WORLD,ierr)
          
            ! When this node is finished, send the blocking MPI_Send to the next machine
            if (mynum /= nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,94,MPI_COMM_WORLD,ierr)

            uniqueid=((igr-1)*maxmach)+mynum
            call MPI_Recv(gdpy,maxmach*maxgrds,MPI_INTEGER,machs(1), 900000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdsi,maxmach*maxgrds,MPI_INTEGER,machs(1),1000000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdpa,maxmach*maxgrds,MPI_INTEGER,machs(1),1100000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdco,maxmach*maxgrds,MPI_INTEGER,machs(1),1200000+uniqueid,MPI_COMM_WORLD,status,ierr)
 
         end if


         if(mynum == 1 .and. model_start .and. verbose) then
            
            print*,"Global Polygons: ",gdpy(1:nnodetot,igr)
            print*,"Global Site: "    ,gdsi(1:nnodetot,igr)
            print*,"Global Patches: " ,gdpa(1:nnodetot,igr)
            print*,"Global Cohorts: " ,gdco(1:nnodetot,igr)

         end if

         ! Calculate the offsets that each machine has
         
         py_off(1,igr) = 0
         si_off(1,igr) = 0
         pa_off(1,igr) = 0
         co_off(1,igr) = 0
         do nm=2,nnodetot
            py_off(nm,igr) = py_off(nm-1,igr) + gdpy(nm-1,igr)
            si_off(nm,igr) = si_off(nm-1,igr) + gdsi(nm-1,igr)
            pa_off(nm,igr) = pa_off(nm-1,igr) + gdpa(nm-1,igr)
            co_off(nm,igr) = co_off(nm-1,igr) + gdco(nm-1,igr)
         end do
         
         ! Calculate the total sizes of the arrays
         
         cgrid%npolygons_global = sum(gdpy(1:nnodetot,igr))
         cgrid%nsites_global    = sum(gdsi(1:nnodetot,igr))
         cgrid%npatches_global  = sum(gdpa(1:nnodetot,igr))
         cgrid%ncohorts_global  = sum(gdco(1:nnodetot,igr))

         print*,"NPOLY GLOB: ",cgrid%npolygons_global
         

         

         ! Calculate the local offsets
         
         cgrid%mach_polygon_offset_index = py_off(mynum,igr)
         cgrid%mach_site_offset_index    = si_off(mynum,igr)
         cgrid%mach_patch_offset_index   = pa_off(mynum,igr)
         cgrid%mach_cohort_offset_index  = co_off(mynum,igr)

          
       endif

       call filltab_globtype(igr)

       call filltab_edtype(igr,0)
       
       if (gdpy(mynum,igr)>0) then
          call filltab_polygontype(igr,1,0)
          call filltab_sitetype(igr,1,1,0)
          call filltab_patchtype(igr,1,1,1,0)
       endif
       
       
    enddo


    do igr = 1,ngrids

       ! Test to see if the var_table has been initialized. If it has
       ! then deallocate its pointers and reset its first flag. These
       ! will be reallocated on the first pass of the filltab_
       ! subroutines.

       cgrid => edgrid_g(igr)

       cgrid%pyglob_id = 0 + cgrid%mach_polygon_offset_index
       
       ! Determine the total number of variables for each grid
       ! These will determine the length of the vt_vector

       call filltab_edtype(igr,1)
       
       ncohorts_g = 0 + cgrid%mach_cohort_offset_index
       npatches_g = 0 + cgrid%mach_patch_offset_index
       nsites_g   = 0 + cgrid%mach_site_offset_index
       
       do ipy = 1,cgrid%npolygons
          
          cpoly => cgrid%polygon(ipy)
          
          cpoly%siglob_id = nsites_g + 0 ! This is the offset for the vtable write
          
          cgrid%pysi_id(ipy) = nsites_g + 1 ! This is the index written in the file
                                            ! for the user to reference

          cgrid%pysi_n(ipy) = cpoly%nsites

          nsites_g = nsites_g + cpoly%nsites
          
          call filltab_polygontype(igr,ipy,1)
          
          do isi = 1,cpoly%nsites
             
             csite => cpoly%site(isi)
             
             csite%paglob_id = npatches_g + 0

             cpoly%sipa_id(isi) = npatches_g + 1

             cpoly%sipa_n(isi) = csite%npatches

             npatches_g = npatches_g + csite%npatches
             
             call filltab_sitetype(igr,ipy,isi,1)
             
             do ipa = 1,csite%npatches

                cpatch => csite%patch(ipa)
                
                cpatch%coglob_id = ncohorts_g + 0

                csite%paco_id(ipa) = ncohorts_g + 1

                csite%paco_n(ipa) = cpatch%ncohorts

                ncohorts_g = ncohorts_g + cpatch%ncohorts

                if (cpatch%ncohorts > 0 ) then
                   
                   call filltab_patchtype(igr,ipy,isi,ipa,1)
                   
                endif

             enddo
             
          enddo
          
       enddo
       if (mynum.eq.1) then
          write(*,"(a)")'--- Mapping Completed'
       endif

       if (mynum.eq.1 .and. model_start .and. verbose) then
          model_start = .false.
          do nv=1,num_var(igr)
!             write(*,"(a,i4,a,i4,a,a)")'Registering: ',nv,' of',num_var(igr),'  ',vt_info(nv,igr)%name
          enddo
       endif 

    enddo


    return
  end subroutine filltab_alltypes
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_globtype(igr)

    use ed_var_tables,only:vtable_edio_r,vtable_edio_i,vtable_edio_c,metadata_edio
    use ed_misc_coms,only:expnme
    use soil_coms,only:soil,ed_nstyp,slz


    implicit none
    
    integer, intent(in) :: igr
    integer :: var_len,max_ptrs,var_len_global
    integer :: nvar
    type(edtype),pointer :: cgrid
    
    cgrid => edgrid_g(igr)
    
    var_len = 1
    var_len_global = 1
    max_ptrs = 1

    ! Note that "90" flags to the IO that this is a scalar
    ! (switched 99 by 90 just to keep the notation)
    nvar=1
    call vtable_edio_i(cgrid%npolygons_global,nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'NPOLYGONS_GLOBAL :90:hist:anal:dail:mont:year')

    call vtable_edio_i(cgrid%npolygons_global,nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'NPOLYGONS_GLOBAL :90:hist:anal:dail:mont:year')
    
    nvar=nvar+1
    call vtable_edio_i(cgrid%nsites_global,nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'NSITES_GLOBAL :90:hist:anal:dail:mont:year')
    
    call vtable_edio_i(cgrid%nsites_global,nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'NSITES_GLOBAL :90:hist:anal:dail:mont:year')
    
    nvar=nvar+1
    call vtable_edio_i(cgrid%npatches_global,nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'NPATCHES_GLOBAL :90:hist:anal:dail:mont:year')
    
    call vtable_edio_i(cgrid%npatches_global,nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'NPATCHES_GLOBAL :90:hist:anal:dail:mont:year')

    nvar=nvar+1
    call vtable_edio_i(cgrid%ncohorts_global,nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'NCOHORTS_GLOBAL :90:hist:anal:dail:mont:year')
    
    call vtable_edio_i(cgrid%ncohorts_global,nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'NCOHORTS_GLOBAL :90:hist:anal:dail:mont:year')
    nvar=nvar+1
    call vtable_edio_i(nzg,nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'NZG :90:hist:anal:dail:mont:year')

    call vtable_edio_i(nzg,nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'NZG :90:hist:anal:dail:mont:year')
    nvar=nvar+1
    var_len        = nzg
    var_len_global = nzg
    call vtable_edio_r(slz(1),nvar,igr,0,0, &
         var_len,var_len_global,max_ptrs,'SLZ :90:hist:anal:dail:mont:year')

    call vtable_edio_r(slz(1),nvar,igr,1,0, &
         var_len,var_len_global,max_ptrs,'SLZ :90:hist:anal:dail:mont:year')

!! SOIL PARAMETERS MAY BE ADDED IN THE FUTURE - 
!!    RIGHT NOW THIS DOESNT REALLY WORK - RGK 7-19-08
!!    nvar=nvar+1
!!    var_len        = ed_nstyp
!!    var_len_global = ed_nstyp
!!    call vtable_edio_r(soil(1)%slmsts,nvar,igr,0,0, &
!!         var_len,var_len_global,max_ptrs,'SOIL_POROSITY :90:hist:anal:dail:mont')

!!    call vtable_edio_r(soil(1)%slmsts,nvar,igr,1,0, &
!!         var_len,var_len_global,max_ptrs,'SOIL_POROSITY :90:hist:anal:dail:mont')
!!    call metadata_edio(nvar,igr,'Porosity of s-ls-sl-sil-l-scl-sicl-cl-sc-sic-c-p','m3/m3','12 classes')


    !! DO NOT INCLUDE TEXT OR STRING VARIABLES, SUCH AS EXP NAME
    !! IT IS PROBLEMATIC FOR SOME APPLICATIONS SUCH AS R
    

    nioglobal=nvar
    return
  end subroutine filltab_globtype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_edtype(igr,init)

    use ed_var_tables,only:vtable_edio_r,vtable_edio_i,vtable_edio_d,metadata_edio

    implicit none

    integer :: init
    integer, intent(in) :: igr
    integer :: var_len,max_ptrs,var_len_global
    integer :: nvar

    type(edtype),pointer :: cgrid

    cgrid => edgrid_g(igr)

    var_len = cgrid%npolygons
    var_len_global = cgrid%npolygons_global
    max_ptrs = 1

    
    nvar=nioglobal
    
    if (associated(cgrid%pysi_id)) then
       nvar = nvar + 1
       call vtable_edio_i(cgrid%pysi_id(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'PYSI_ID :11:hist:anal:dail:mont:year')

       call metadata_edio(nvar,igr,'Polygons first site indices','NA','ipoly')

    endif
    
    if (associated(cgrid%pysi_n)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%pysi_n(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'PYSI_N :11:hist:anal:dail:mont:year')

       call metadata_edio(nvar,igr,'Number of sites per polygon','NA','ipoly')

    endif

    if (associated(cgrid%walltime_py)) then
       nvar=nvar+1
       call vtable_edio_d(cgrid%walltime_py(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WALLTIME_PY :11:hist:anal:dail:mont:year')
       call metadata_edio(nvar,igr,'simulation walltime of each polygon for fast integration','seconds','ipoly')
       
    endif

    
    if (associated(cgrid%lat)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%lat(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LATITUDE :11:hist:anal:dail:mont:year')

       call metadata_edio(nvar,igr,'Latitude of Polygon','decimal degrees','ipoly')

    endif
    
    if (associated(cgrid%lon)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%lon(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LONGITUDE :11:hist:anal:dail:mont:year') 
       
        call metadata_edio(nvar,igr,'Longitude of Polygon','decimal degrees','ipoly')
       
    endif
    
    if (associated(cgrid%natm)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%natm(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NATM :11:hist') 

       call metadata_edio(nvar,igr,'Number of atm cells per polygon','NA','ipoly')
       
    endif
    
    if (associated(cgrid%xatm)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%xatm(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'XATM :11:hist:anal:dail:mont:year') 

       call metadata_edio(nvar,igr,'Atm. cell x-indices of polygon','NA','ipoly')

    endif
    
    if (associated(cgrid%yatm)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%yatm(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'YATM :11:hist:anal:dail:mont:year') 

       call metadata_edio(nvar,igr,'Atm cell y-indices of polygon','NA','ipoly')

    endif

    if (associated(cgrid%ntext_soil)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%ntext_soil(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NTEXT_SOIL :12:hist:anal:dail:mont:year')   
       call metadata_edio(nvar,igr,'Polygon mode soil class','OGE2 Class','ipoly-ngz')

    endif
    
    if (associated(cgrid%lsl)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%lsl(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LSL :11:hist:anal:dail:mont:year')
       call metadata_edio(nvar,igr,'Index of lowest soil layer','NA','ipoly')
       
    endif
    
    if (associated(cgrid%wbar)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%wbar(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WBAR :11:hist') 
       call metadata_edio(nvar,igr,'Polygon average topographic moisture index','NA','ipoly')
       
    endif
    
    if (associated(cgrid%Te)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Te(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TE :11:hist') 
       call metadata_edio(nvar,igr,'NA','NA','ipoly')
    endif
    
    if (associated(cgrid%zbar)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%zbar(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'ZBAR :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon average water table depth','[m]','ipoly')
    endif
 
!!==  Variable not used    
!!    if (associated(cgrid%tau)) then
!!       nvar=nvar+1
!!       call vtable_edio_r(cgrid%tau(1),nvar,igr,init,cgrid%pyglob_id, &
!!            var_len,var_len_global,max_ptrs,'TAU :11:hist') 
!!       call metadata_edio(nvar,igr,'NA','NA','ipoly')
!!    endif
    
    if (associated(cgrid%sheat)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%sheat(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'SHEAT :11:hist') 
       call metadata_edio(nvar,igr,'soil heat pool for lateral hydrology','NA','ipoly')
    endif
    
    if (associated(cgrid%baseflow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%baseflow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'BASEFLOW :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'loss of water from site to watershed discharge','kg/m2/s','ipoly')
    endif
    
    if (associated(cgrid%runoff)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%runoff(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'RUNOFF :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'NA','NA','ipoly')
    endif
    
    if (associated(cgrid%swliq)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%swliq(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'SWLIQ :11:hist:anal') 
       call metadata_edio(nvar,igr,'NA','NA','ipoly')
    endif

    if (associated(cgrid%total_agb)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_agb(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_AGB :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon Total Above Ground Biomass','[kgC/m2]','ipoly')
    endif
    
    if (associated(cgrid%total_basal_area)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_basal_area(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon Total Basal Area','[cm2/m2]','ipoly')
       
    endif

    if (associated(cgrid%total_agb_growth)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_agb_growth(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_AGB_GROWTH:11:hist:anal:year') 
        call metadata_edio(nvar,igr,'Polygon AGB gain through growth','[kgC/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%total_agb_mort)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_agb_mort(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_AGB_MORT :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon AGB lost due to mortality','[kgC/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%total_agb_recruit)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_agb_recruit(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_AGB_RECRUIT :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon AGB used to generate recruits','[kgC/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%total_basal_area_growth)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_basal_area_growth(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_GROWTH :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon basal area gained through growth ','[kgC/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%total_basal_area_mort)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_basal_area_mort(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_MORT :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon basal area lost through growth ','[cm2/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%total_basal_area_recruit)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%total_basal_area_recruit(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_RECRUIT :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Polygon basal area gained by recruits','[cm2/m2/yr]','ipoly')
    endif
    
    if (associated(cgrid%load_adjacency)) then
       nvar=nvar+1
       call vtable_edio_i(cgrid%load_adjacency(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LOAD_ADJACENCY :11:hist') 
       call metadata_edio(nvar,igr,'Load Adjacency','[NA]','ipoly')
    endif
    
    if (associated(cgrid%cosz)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%cosz(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'COSZ :11:hist') 
       call metadata_edio(nvar,igr,'Cosine of the zenith angle','[a/h]','ipoly')
    endif
    
    if (associated(cgrid%cbudget_initialstorage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%cbudget_initialstorage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CBUDGET_INITIALSTORAGE :11:hist') 
       call metadata_edio(nvar,igr,'Vegetation and soil carbon,at start of budget-averaging','[kgC/m2]','ipoly')
    endif
        
    if (associated(cgrid%cbudget_nep)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%cbudget_nep(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CBUDGET_NEP :11:hist') 
       call metadata_edio(nvar,igr,'Polygon average net ecosystem production','[kgC/m2/day]','ipoly')
    endif
    
    if (associated(cgrid%nbudget_initialstorage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%nbudget_initialstorage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NBUDGET_INITIALSTORAGE :11:hist') 
       call metadata_edio(nvar,igr,'Veg and soil nitrogen, at start of budget-averaging','[kgN/m2]','ipoly')
    endif
    
    if (associated(cgrid%basal_area)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%basal_area(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'BASAL_AREA :146:hist:anal:dail:mont:year') 
       call metadata_edio(nvar,igr,'Polygon basal area profile','[cm2/m2]','ipoly - n_dbh - n_pft')
    endif
    
    if (associated(cgrid%agb)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%agb(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AGB :146:hist:anal:dail:mont:year') 
       call metadata_edio(nvar,igr,'Polygon above ground biomass profile','[kgC/m2]','ipoly - n_dbh - n_pft')
    endif
    
    if (associated(cgrid%pldens)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%pldens(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'PLDENS :146:hist:anal:dail:mont:year') 
       call metadata_edio(nvar,igr,'Polygon plant density profile','[#/m2]','ipoly - n_dbh - n_pft')
    endif
    
    if (associated(cgrid%bseeds)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%bseeds(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'BSEEDS :146:hist:anal:dail:mont:year') 
       call metadata_edio(nvar,igr,'Polygon seed biomass','[kgC/m2]','ipoly - n_dbh - n_pft')
    endif
    
    
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    if (associated(cgrid%avg_vapor_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_vapor_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VAPOR_VC :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'polygon vegetation to canopy air vapor flux','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_dew_cg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_dew_cg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_DEW_CG :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon averaged dew to ground flux','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_vapor_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_vapor_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VAPOR_GC :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'polygon moisture flux ground to canopy air','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_wshed_vg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_wshed_vg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_WSHED_VG :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged water shed from vegetation to ground','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_intercepted)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_intercepted(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_INTERCEPTED :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged intercepted precipitation by vegetation','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_vapor_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_vapor_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VAPOR_AC :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'polygon vapor flux atmosphere to canopy air','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_transp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_transp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_TRANSP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'polygon transpiration from stomata to canopy air space','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_evap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_evap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_EVAP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon averaged evap/dew from ground and leaves to CAS','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_smoist_gg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_smoist_gg(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SMOIST_GG :12:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged soil moisture flux,layer nzg is flux with CAS','[kg/m2/s]','ipoly-nzg') 
    endif
    
    if (associated(cgrid%avg_smoist_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_smoist_gc(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SMOIST_GC :12:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged soil moisture sink to transpiration','[kg/m2/s]','ipoly-nzg') 
    endif
    
    if (associated(cgrid%avg_runoff)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_runoff(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RUNOFF :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon average surface runoff','[kg/m2/s]','NA') 
    endif
    
    if (associated(cgrid%avg_drainage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_drainage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_DRAINAGE :11:hist:anal') 
       call metadata_edio(nvar,igr,'polygon average water flux through lower soil layer','[kg/m2/s]','ipoly') 
    endif
     
    if (associated(cgrid%avg_drainage_heat)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_drainage_heat(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_DRAINAGE_HEAT :11:hist:anal') 
       call metadata_edio(nvar,igr,'polygon average internal energy loss through lower soil layer','[W/m2]','ipoly') 
    endif
   
    if (associated(cgrid%aux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%aux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AUX :11:hist:anal') 
       call metadata_edio(nvar,igr,'Auxillary variable - user discretion,see rk4_derivs.f90','[user-defined]','ipoly') 
    endif
    
    if (associated(cgrid%aux_s)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%aux_s(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AUX_S :12:hist:anal') 
       call metadata_edio(nvar,igr,'Soil layer discretized, auxillary variable, see rk4_derivs.f90','[user-defined]','ipoly-nzg') 
    endif
    
    if (associated(cgrid%avg_sensible_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sensible_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_VC :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged vegetation to canopy air sensible heat flux','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_qwshed_vg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_qwshed_vg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_QWSHED_VG :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged internal energy flux of water shed from vegetation to ground','[W/m2]','ipoly') 
    endif
     
    if (associated(cgrid%avg_qintercepted)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_qintercepted(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_QINTERCEPTED :11:hist:anal') 
       call metadata_edio(nvar,igr,&
'Polygon averaged internal energy flux of intercepted precipitation by vegetation','[W/m2]','ipoly') 
    endif
   
    if (associated(cgrid%avg_sensible_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sensible_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_GC :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged sensible heat flux ground to canopy air','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_sensible_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sensible_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_AC :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon averaged sensible heat flux atmosphere  to canopy','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_sensible_gg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sensible_gg(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_GG :12:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon averaged soil sensible heat flux,layer nzg is flux with CAS ','[W/m2]','ipoly-nzg') 
    endif
    
    if (associated(cgrid%avg_runoff_heat)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_runoff_heat(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RUNOFF_HEAT :11:hist:anal') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    ! ---------------------------------------------

    if (associated(cgrid%avg_lai_ebalvars)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_lai_ebalvars(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_LAI_EBALVARS :199:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Energy Balance Variables','[variable]','ipoly - 4 - 3') 
    endif
    
    if (associated(cgrid%avg_gpp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_gpp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_GPP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average GPP','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%lai)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%lai(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LAI :11:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Polygon  LAI','[m2/m2]','ipoly') 
    endif

    if (associated(cgrid%wpa)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%wpa(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WPA :11:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Polygon  WPA','[m2/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_lma)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_lma(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_LMA :11:hist:dail') 
       call metadata_edio(nvar,igr,'Polygon LMA','[NA]','ipoly') 
    endif

    if (associated(cgrid%wai)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%wai(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WAI :11:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Polygon  WAI','[m2/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_leaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_leaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_LEAF_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Leaf Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_root_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_root_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ROOT_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Root Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_growth_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_growth_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_GROWTH_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Growth Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_storage_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_storage_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_STORAGE_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Storage Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_vleaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_vleaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VLEAF_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Virtual Leaf Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_plant_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_plant_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_PLANT_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Plant Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_htroph_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_htroph_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_HTROPH_RESP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Heterotrohic Respiration','[umol/m2/s]','ipoly') 
    endif

    if (associated(cgrid%avg_leaf_drop)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_leaf_drop(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_LEAF_DROP :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Leaf loss','[kgC/m2/yr]','ipoly') 
    endif

    if (associated(cgrid%avg_leaf_maintenance)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_leaf_maintenance(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_LEAF_MAINTENANCE :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Leaf maintenance cost','[kgC/m2/yr]','ipoly') 
    endif

    if (associated(cgrid%avg_root_maintenance)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_root_maintenance(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ROOT_MAINTENANCE :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average fine root maintenance cost','[kgC/m2/yr]','ipoly') 
    endif


       !!! added for NACP intercomparison (MCD)
    if (associated(cgrid%avg_sfcw_depth)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sfcw_depth(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SNOWDEPTH :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Poly Avg. Snow Depth ','[m]','ipoly') 
    endif
    if (associated(cgrid%avg_sfcw_energy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sfcw_energy(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SNOWENERGY :11:hist:anal') 
       call metadata_edio(nvar,igr,'Poly Avg. Snow Energy ','[J/kg]','ipoly') 
    endif
    if (associated(cgrid%avg_sfcw_mass)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sfcw_mass(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SNOWMASS :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Poly Avg. Snow Mass (SWE) ','[kg/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_sfcw_tempk)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sfcw_tempk(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SNOWTEMP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Poly Avg. Snow Temperature','[K]','ipoly') 
    endif
    if (associated(cgrid%avg_sfcw_fracliq)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_sfcw_fracliq(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SNOWFRACLIQ :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Poly Avg. Snow liquid fraction','[proportion]','ipoly') 
    endif
    if (associated(cgrid%avg_bdead)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_bdead(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BDEAD :11:hist:opti:dail:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass - structural','[kgC/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_bstorage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_bstorage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BSTORAGE :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass - storage','[kgC/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_bseeds)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_bseeds(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BSEEDS :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass - seeds','[kgC/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_balive)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_balive(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BALIVE :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass -- living','[kgC/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_bleaf)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_bleaf(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BLEAF :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass -- leaf','[kgC/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_broot)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_broot(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BROOT :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass -- fine roots','[kgC/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_bsapwood)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_bsapwood(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_BSAPWOOD :11:hist:anal:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Biomass -- sapwood','[kgC/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_fsc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_fsc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_FSC :11:hist:anal:opti:dail:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Fast Soil Carbon','[kg/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_ssc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_stsc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SSC :11:hist:anal:opti:dail:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Slow Soil Carbon','[kg/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_stsc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_stsc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_STSC :11:hist:anal:opti:year:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Structural Soil Carbon','[kg/m2]','ipoly') 
    endif


    if (associated(cgrid%avg_fsn)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_fsn(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_FSN :11:hist:anal:opti:dail:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Fast Soil Nitrogen','[kg/m2]','ipoly') 
    endif
    if (associated(cgrid%avg_msn)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_msn(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_MSN :11:hist:anal:opti:dail:year') 
       call metadata_edio(nvar,igr,'Poly Avg. Mineralized Soil Carbon','[kg/m2]','ipoly') 
    endif




     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
    if (associated(cgrid%Cleaf)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Cleaf(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CLEAF :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon','?','ipoly') 
    endif
    if (associated(cgrid%Croot)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Croot(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CROOT :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Fine Root Carbon','?','ipoly') 
    endif
    if (associated(cgrid%Cstore)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%cstore(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CSTORE :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Storage/TNC Carbon','?','ipoly') 
    endif
    if (associated(cgrid%Ccwd)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ccwd(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CCWD :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Coarse Woody Debris Carbon','?','ipoly') 
    endif
    if (associated(cgrid%Nleaf)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nleaf(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NLEAF :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Nitrogen','?','ipoly') 
    endif
    if (associated(cgrid%Ndead)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ndead(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NDEAD :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Wood Nitrogen','?','ipoly') 
    endif
    if (associated(cgrid%Nroot)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nroot(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NROOT :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Fine Root Nitrogen','?','ipoly') 
    endif
    if (associated(cgrid%Nstore)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nstore(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NSTORE :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Storage Nitrogen','?','ipoly') 
    endif
    if (associated(cgrid%Ncwd)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ncwd(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NCWD :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon','?','ipoly') 
    endif

     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
    if (associated(cgrid%Cleaf_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Cleaf_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CLEAF_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon Growth','?','ipoly') 
    endif
   if (associated(cgrid%Croot_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Croot_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CROOT_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Fine Root Carbon Growth','?','ipoly') 
    endif
   if (associated(cgrid%Cdead_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Cdead_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CDEAD_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Wood Carbon Growth','?','ipoly') 
    endif
   if (associated(cgrid%Cstore_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Cstore_Grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CSTORE_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Storage Carbon growth','?','ipoly') 
    endif
   if (associated(cgrid%Cleaf_litter_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Cleaf_litter_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CLEAF_LITTER_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Litter Carbon Flux','?','ipoly') 
    endif
   if (associated(cgrid%Croot_litter_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Croot_litter_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CROOT_LITTER_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Fine Root Litter Carbon Flux','?','ipoly') 
    endif
   if (associated(cgrid%Ccwd_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ccwd_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'CCWD_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Coarse Woody Debris Carbon FLUX','?','ipoly') 
    endif
   if (associated(cgrid%Nleaf_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nleaf_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NLEAF_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf Nitrogen growth','?','ipoly') 
    endif
   if (associated(cgrid%Nroot_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nroot_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NROOT_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Fine Root Nitrogen Growth','?','ipoly') 
    endif
   if (associated(cgrid%Ndead_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ndead_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NDEAD_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Wood Nitrogen growth','?','ipoly') 
    endif
   if (associated(cgrid%Nstore_grow)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nstore_grow(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NSTORE_GROW :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Storage Nitrogen Growth','?','ipoly') 
    endif
   if (associated(cgrid%Nleaf_litter_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nleaf_litter_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NLEAF_LITTER_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Leaf litter Nitrogen flux','?','ipoly') 
    endif
   if (associated(cgrid%Nroot_litter_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nroot_litter_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NROOT_LITTER_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. fine root litter Nitrogen flux','?','ipoly') 
    endif
   if (associated(cgrid%Ncwd_flux)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ncwd_flux(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NCWD_FLUX :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. CWD Nitrogen flux','?','ipoly') 
    endif
   if (associated(cgrid%Nbiomass_uptake)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nbiomass_uptake(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NBIOMASS_UPTAKE :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Nitrogen Uptake','?','ipoly') 
    endif
   if (associated(cgrid%Ngross_min)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Ngross_min(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NGROSS_MIN :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Gross Nitrogen mineralization','?','ipoly') 
    endif
   if (associated(cgrid%Nnet_min)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%Nnet_min(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'NNET_MIN :11:dail') 
       call metadata_edio(nvar,igr,'Poly Avg. Net Nitrogen mineralization','?','ipoly') 
    endif

    ! ----------------------------------------------
    
    if (associated(cgrid%avg_nir_beam)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_nir_beam(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_NIR_BEAM:11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Incident Near Infrared Beam Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_nir_diffuse)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_nir_diffuse(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_NIR_DIFFUSE :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Incident Near Infrared Diffuse Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_par_beam)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_par_beam(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_PAR_BEAM :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Incident Beam Photosynthetically Active Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_par_diffuse)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_par_diffuse(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_PAR_DIFFUSE :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Incident Diffuse Photosynthetically Active Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_atm_tmp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_atm_tmp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ATM_TMP :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Atmospheric Temperature at Reference Height','[K]','ipoly') 
    endif
    
    if (associated(cgrid%avg_atm_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_atm_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ATM_SHV :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Averaged Atmospheric Specific Humidity at Reference Height','[kg/kg]','ipoly') 
    endif
    
    if (associated(cgrid%avg_rshort)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_rshort(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RSHORT :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Total Incident Shortwave Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_rshort_diffuse)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_rshort_diffuse(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RSHORT_DIFFUSE :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Diffuse Incident Shortwave Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_rlong)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_rlong(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RLONG :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Incident Longwave Radiation','[W/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_pcpg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_pcpg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_PCPG :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Total Precipitation Rate','[kg/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_qpcpg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_qpcpg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_QPCPG:11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Precipitation Internal Energy Deposition Rate','[W/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_dpcpg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_dpcpg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_DPCPG :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Total Precipitation Depth Rate ','[mm/m2/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_vels)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_vels(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VELS :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Wind Magnitude (with instability correction)','[m/s]','ipoly') 
    endif
    
    if (associated(cgrid%avg_atm_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_atm_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ATM_PRSS :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Atmospheric Pressure at Ref. Height','[Pa]','ipoly') 
    endif
    
    if (associated(cgrid%avg_exner)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_exner(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_EXNER :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Exner Correction','[????]','ipoly') 
    endif
    
    if (associated(cgrid%avg_geoht)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_geoht(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_GEOHT :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Geopotential of Met. Forcing Refernce Height','[m]','ipoly') 
    endif
    
    if (associated(cgrid%avg_atm_co2)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_atm_co2(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ATM_CO2 :11:hist:opti:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Atmospheric CO2 Concentration at Ref. Height','[ppm]','ipoly') 
    endif
    
    if (associated(cgrid%avg_albedt)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_albedt(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_ALBEDT :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Surface Albedo','[W/W]','ipoly') 
    endif
    
    if (associated(cgrid%avg_rlongup)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_rlongup(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_RLONGUP :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Upwelling Longwave Radiation','[W/m2]','ipoly') 
    endif

    if (associated(cgrid%max_veg_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%max_veg_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MAX_VEG_TEMP :11:anal') 
       call metadata_edio(nvar,igr,'Temp of the hottest cohort in the polygon','[K]','ipoly') 
    endif     

    if (associated(cgrid%min_veg_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%min_veg_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MIN_VEG_TEMP :11:anal') 
       call metadata_edio(nvar,igr,'Temp of the coldest cohort in the polygon','[K]','ipoly') 
    endif

    if (associated(cgrid%max_soil_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%max_soil_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MAX_SOIL_TEMP :11:anal') 
       call metadata_edio(nvar,igr,'Temp of the hottest soil layer in the polygon','[K]','ipoly') 
    endif    

    if (associated(cgrid%min_soil_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%min_soil_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MIN_SOIL_TEMP :11:anal') 
       call metadata_edio(nvar,igr,'Temp of the coldest soil layer in the polygon','[K]','ipoly') 
    endif    
    
    if (associated(cgrid%avg_veg_energy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_veg_energy(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VEG_ENERGY :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Internal Energy of Vegetation','[J/m2]','ipoly') 
    endif

    if (associated(cgrid%avg_veg_hcap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_veg_hcap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VEG_HCAP :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average of Vegetation heat capacity','[J/m2/K]','ipoly') 
    endif

    if (associated(cgrid%avg_veg_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_veg_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VEG_TEMP :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Temperature of Vegetation','[K]','ipoly') 
    endif

    if (associated(cgrid%avg_veg_fliq)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_veg_fliq(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VEG_FLIQ :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Liquid fraction of Vegetation Sfc Water','[--]','ipoly') 
    endif
    
    if (associated(cgrid%avg_veg_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_veg_water(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_VEG_WATER :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Resident Leaf Surface Water','[kg/m2]','ipoly') 
    endif
    
    if (associated(cgrid%avg_can_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_TEMP :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Temperature of Canopy Air Space','[K]','ipoly') 
    endif
    
    if (associated(cgrid%avg_can_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_SHV :11:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Specific Humidity of Canopy Air','[kg/kg]','NA') 
    endif
    
    if (associated(cgrid%avg_can_co2)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_co2(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_CO2 :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average CO2 mixing ratio of Canopy Air','[umol/mol]','NA') 
    endif
    
    if (associated(cgrid%avg_can_rhos)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_rhos(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_RHOS :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Density of Canopy Air','[kg/m3]','NA') 
    endif
    
    if (associated(cgrid%avg_can_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_PRSS :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Canopy Air Pressure','[Pa]','NA') 
    endif
    
    if (associated(cgrid%avg_can_theta)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_theta(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_THETA :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Canopy Air Potential temperature','[K]','NA') 
    endif
    
    if (associated(cgrid%avg_can_enthalpy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_enthalpy(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_ENTHALPY :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Canopy Air Enthalpy','[J/kg]','NA') 
    endif
    
    if (associated(cgrid%avg_can_depth)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_can_depth(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_CAN_DEPTH :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Canopy height','[m]','NA') 
    endif
    
    if (associated(cgrid%avg_soil_energy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_soil_energy(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SOIL_ENERGY :12:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Volumetric Soil Water','[m/m]','ipoly - nzg') 
    endif
    
    if (associated(cgrid%avg_soil_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_soil_water(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SOIL_WATER :12:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Volumetric Soil Water','[m/m]','ipoly - nzg') 
    endif
    
    if (associated(cgrid%avg_soil_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_soil_temp(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SOIL_TEMP :12:hist:anal:opti') 
       call metadata_edio(nvar,igr,'Polygon Average Soil Temperature','[K]','ipoly - nzg') 
    endif

    if (associated(cgrid%avg_soil_fracliq)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_soil_fracliq(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SOIL_FRACLIQ :12:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Soil Fraction Liquid','[proportion]','ipoly - nzg') 
    endif
    
    if (associated(cgrid%avg_soil_wetness)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_soil_wetness(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SOIL_WETNESS :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Soil Wetness RELATIVE TO WILTING POINT','[m3/m3]','ipoly') !relative to wilting point
    endif
    
    if (associated(cgrid%avg_skin_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_skin_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_SKIN_TEMP :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Temperature of all surfaces','[K]','ipoly') 
    endif
    
    if (associated(cgrid%avg_available_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%avg_available_water(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AVG_AVAILABLE_WATER :11:hist:anal') 
       call metadata_edio(nvar,igr,'Polygon Average Available water','[K]','ipoly') 
    endif
    
    ! Daily and monthly variables. Note that all these variables need to be stored at the
    ! history file, because the averaging can be resumed...
    
    if(associated(cgrid%dmean_gpp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_gpp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_GPP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Polygon Average Daily Integrated Gross Primary Productivity','[kgC/m2/yr]','ipoly') 
    endif
    
    
    if(associated(cgrid%dmean_pcpg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_pcpg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_PCPG :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily precipitation','[kg/m2/day]','ipoly') 
    endif

    if(associated(cgrid%dmean_runoff)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_runoff(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_RUNOFF :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily surface runoff','[kg/m2/day]','ipoly') 
    endif

    if(associated(cgrid%dmean_drainage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_drainage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_DRAINAGE :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily water flux through lower soil layer','[kg/m2/day]','ipoly') 
    endif

    if(associated(cgrid%dmean_vapor_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_vapor_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_AC :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily vapor flux atm->canopy','[kg/m2/day]','ipoly') 
    endif


    if(associated(cgrid%dmean_vapor_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_vapor_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_GC :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily vapor flux ground->canopy','[kg/m2/day]','ipoly') 
    endif


    if(associated(cgrid%dmean_vapor_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_vapor_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_VC :11:hist:dail') 
       call metadata_edio(nvar,igr,'total daily vapor flux vegetation->canopy','[kg/m2/day]','ipoly') 
    endif

    if(associated(cgrid%dmean_evap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_evap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_EVAP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean (leaf+soil) evaporation','[kg/m2/s]','ipoly') 
    endif
    
    if(associated(cgrid%dmean_transp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_transp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_TRANSP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_sensible_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_sensible_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_VC :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_sensible_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_sensible_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_GC :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_sensible_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_sensible_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_AC :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_plresp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_plresp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_PLRESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_rh)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_rh(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_RH :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_leaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_leaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_LEAF_RESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_root_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_root_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ROOT_RESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_growth_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_growth_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_GROWTH_RESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_storage_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_storage_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_STORAGE_RESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_vleaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_vleaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VLEAF_RESP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_nep)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_nep(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_NEP :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_soil_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_soil_temp(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_SOIL_TEMP :12:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_soil_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_soil_water(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_SOIL_WATER :12:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_fs_open)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_fs_open(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_FS_OPEN :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_fsw)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_fsw(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_FSW :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_fsn)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_fsn(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_FSN :11:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%dmean_gpp_dbh)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_gpp_dbh(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_GPP_DBH :16:hist:dail') 
       call metadata_edio(nvar,igr,'Polygon Averaged by DBH, Daily Integrated Gross Primary Production' &
            ,'[kgC/m2/yr]','ipoly - ndbh') 
    endif
    
    if(associated(cgrid%dmean_can_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_TEMP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_can_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_SHV :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy specific humidity','[kg/kg]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_can_co2)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_co2(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_CO2 :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy CO2 mixing ratio','[umol/mol]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_can_rhos)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_rhos(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_RHOS :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy air density','[kg/m3]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_can_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_PRSS :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy air pressure','[Pa]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_can_theta)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_can_theta(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CAN_THETA :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean canopy air potential temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_veg_energy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_veg_energy(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VEG_ENERGY :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean vegetation internal energy','[J/m2]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_veg_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_veg_water(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VEG_WATER :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean vegetation surface water','[kg/m2]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_veg_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_veg_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VEG_TEMP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean vegetation temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_veg_hcap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_veg_hcap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_VEG_HCAP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean vegetation heat capacity','[J/m2/K]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_atm_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_atm_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ATM_TEMP :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean air temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_atm_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_atm_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ATM_SHV :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean air specific humidity','[kg/kg]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_atm_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_atm_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ATM_PRSS :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean air pressure','[ Pa]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_atm_vels)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_atm_vels(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ATM_VELS :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean wind speed','[m/s]','ipoly') 
    end if

    if(associated(cgrid%dmean_rshort)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_rshort(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_RSHORT :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean shortwave radiation','[w/m2]','ipoly') 
    end if
    
    if(associated(cgrid%dmean_rlong)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_rlong(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_RLONG :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean longwave radiation','[w/m2]','ipoly') 
    end if

    if(associated(cgrid%dmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_co2_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/s]','ipoly') 
    end if

    if(associated(cgrid%dmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_water_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/s]','ipoly') 
    end if

    if(associated(cgrid%dmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%dmean_energy_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL :11:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/s]','ipoly') 
    end if

    if(associated(cgrid%lai_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%lai_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'LAI_PFT :14:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Leaf Area Index','[m2/m2]','NA') 
    endif
    
    if(associated(cgrid%wpa_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%wpa_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WPA_PFT :14:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Wood Projected Area','[m2/m2]','NA') 
    endif
    
    if(associated(cgrid%wai_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%wai_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WAI_PFT :14:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Wood Area Index','[m2/m2]','NA') 
    endif
    
    if(associated(cgrid%mmean_gpp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_gpp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_GPP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_evap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_evap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_EVAP :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly Mean (leaf+soil) Evaporation Rate','[kg/m2/s]','ipoly') 
    endif
    
    if(associated(cgrid%mmean_transp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_transp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_TRANSP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if(associated(cgrid%mmean_sensible_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_sensible_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_AC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_sensible_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_sensible_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_GC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_sensible_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_sensible_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_VC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_vapor_ac)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_vapor_ac(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_AC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_vapor_gc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_vapor_gc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_GC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_vapor_vc)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_vapor_vc(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_VC :11:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif
    
    if(associated(cgrid%mmean_nep)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_nep(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_NEP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
   
    if(associated(cgrid%mmean_soil_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_soil_temp(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_SOIL_TEMP :12:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if(associated(cgrid%mmean_soil_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_soil_water(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_SOIL_WATER :12:hist:mont')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif
    
    if(associated(cgrid%mmean_fs_open)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_fs_open(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_FS_OPEN :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_fsw)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_fsw(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_FSW :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_fsn)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_fsn(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_FSN :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if(associated(cgrid%mmean_plresp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_plresp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_PLRESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_rh)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_rh(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_RH :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_leaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_leaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_LEAF_RESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_root_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_root_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ROOT_RESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_growth_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_growth_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_GROWTH_RESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_storage_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_storage_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_STORAGE_RESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_vleaf_resp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_vleaf_resp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VLEAF_RESP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_gpp_dbh)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_gpp_dbh(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_GPP_DBH :16:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_lai_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_lai_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_LAI_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_wpa_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_wpa_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_WPA_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_wai_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_wai_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_WAI_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%mmean_can_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_TEMP :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_can_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_SHV :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy specific humidity','[kg/kg]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_can_co2)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_co2(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_CO2 :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy CO2 mixing ratio','[umol/mol]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_can_rhos)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_rhos(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_RHOS :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy air density','[kg/m3]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_can_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_PRSS :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy air pressure','[Pa]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_can_theta)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_can_theta(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CAN_THETA :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean canopy air potential temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_veg_water)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_veg_water(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VEG_WATER :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean vegetation surface water','[kg/m2]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_veg_energy)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_veg_energy(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VEG_ENERGY :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean vegetation internal energy','[J/m2]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_veg_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_veg_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VEG_TEMP :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean vegetation temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_veg_hcap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_veg_hcap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_VEG_HCAP :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean vegetation heat capacity','[J/m2/K]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_rshort)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_rshort(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_RSHORT :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean downwelling solar radiation','[w/m2]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_rlong)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_rlong(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_RLONG :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean downwelling longwave radiation','[w/m2]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_atm_temp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_atm_temp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ATM_TEMP :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean air temperature','[K]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_atm_shv)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_atm_shv(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ATM_SHV :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean air specific humidity','[kg/kg]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_atm_prss)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_atm_prss(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ATM_PRSS :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean air pressure','[ Pa]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_atm_vels)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_atm_vels(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ATM_VELS :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean wind speed','[m/s]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_pcpg)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_pcpg(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_PCPG :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean precipitation rate','[kg/m2/s]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_runoff)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_runoff(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_RUNOFF :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean runoff','[kg/m2/s]','ipoly') 
    end if
    
    if(associated(cgrid%mmean_drainage)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_drainage(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_DRAINAGE :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean drainage','[kg/m2/day]','ipoly') 
    end if

    if(associated(cgrid%mmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_co2_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/s]','ipoly') 
    end if

    if(associated(cgrid%mmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_water_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/s]','ipoly') 
    end if

    if(associated(cgrid%mmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%mmean_energy_residual(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL :11:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/s]','ipoly') 
    end if
    
    if(associated(cgrid%bseeds_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%bseeds_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'BSEEDS_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%agb_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%agb_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'AGB_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'Above-ground biomass by PFT','[kgC/m2]','NA') 
    endif
    
    if(associated(cgrid%ba_pft)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%ba_pft(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'BA_PFT :14:hist:mont') 
       call metadata_edio(nvar,igr,'Basal area by PFT','[cm2/m2]','NA') 
    endif

    if(associated(cgrid%stdev_gpp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_gpp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_GPP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%stdev_evap)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_evap(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_EVAP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%stdev_transp)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_transp(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_TRANSP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%stdev_sensible)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_sensible(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_SENSIBLE :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%stdev_nep)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_nep(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_NEP :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%stdev_rh)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%stdev_rh(1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'STDEV_RH :11:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if(associated(cgrid%disturbance_rates)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%disturbance_rates(1,1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'DISTURBANCE_RATES :155:hist:mont') 
       call metadata_edio(nvar,igr,'Disturbance Rates','[NA]','NA') 
    end if

    if(associated(cgrid%workload)) then
       nvar=nvar+1
       call vtable_edio_r(cgrid%workload(1,1),nvar,igr,init,cgrid%pyglob_id, &
            var_len,var_len_global,max_ptrs,'WORKLOAD :19:hist:mont') 
       call metadata_edio(nvar,igr,'Disturbance Rates','[NA]','NA') 
    end if
 
    if (init == 0) niogrid=nvar-nioglobal
    
    return
    
  end subroutine filltab_edtype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_polygontype(igr,ipy,init)

    use ed_var_tables,only:vtable_edio_r,vtable_edio_i,metadata_edio

    implicit none

    integer :: init
    integer, intent(in) :: igr,ipy
    integer :: var_len,max_ptrs,var_len_global
    integer :: nvar
    
    type(polygontype),pointer :: cpoly
    
    if (.not.associated(edgrid_g(igr)%polygon(ipy)%sipa_id)) then
       print*,"RETURNING FROM FILLTAB_POLYGONTYPE",igr,ipy,init
       return
    endif

    cpoly => edgrid_g(igr)%polygon(ipy)

    var_len = cpoly%nsites
    var_len_global = edgrid_g(igr)%nsites_global
    max_ptrs = edgrid_g(igr)%npolygons

    nvar=nioglobal+niogrid
    
    if (associated(cpoly%sipa_id)) then
       nvar=nvar+1
       call vtable_edio_i(cpoly%sipa_id(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'SIPA_ID :21:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%sipa_n)) then
       nvar=nvar+1
       call vtable_edio_i(cpoly%sipa_n(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'SIPA_N :21:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%patch_count)) then
       nvar=nvar+1
       call vtable_edio_i(cpoly%patch_count(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'PATCH_COUNT :21:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%sitenum)) then
       nvar=nvar+1
       call vtable_edio_i(cpoly%sitenum(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'SITENUM :21:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    

    if (associated(cpoly%lsl)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%lsl(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'LSL_SI :21:hist:dail:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif   

    if (associated(cpoly%area)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%area(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AREA_SI:21:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%patch_area)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%patch_area(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PATCH_AREA:21:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%elevation)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%elevation(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'ELEVATION :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%slope)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%slope(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'SLOPE :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%aspect)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%aspect(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'ASPECT :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%num_landuse_years)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%num_landuse_years(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'NUM_LANDUSE_YEARS :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%TCI)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%TCI(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'TCI :21:hist') 
       call metadata_edio(nvar,igr,'Topographic convergence index','[NA]','NA') 
    endif      

    if (associated(cpoly%pptweight)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%TCI(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'pptweight :21:hist') 
       call metadata_edio(nvar,igr,'precip lapse weighting','[NA]','NA') 
    endif      

    if (associated(cpoly%hydro_next)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%hydro_next(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'HYDRO_NEXT :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%hydro_prev)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%hydro_prev(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'HYDRO_PREV :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%moist_W)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%moist_W(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'MOIST_W :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%moist_f)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%moist_f(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'MOIST_F :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif  

    if (associated(cpoly%moist_tau)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%moist_tau(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'MOIST_TAU :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%moist_zi)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%moist_zi(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'MOIST_ZI :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif 

    if (associated(cpoly%baseflow)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%baseflow  (1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BASEFLOW_SI :21:hist') 
       call metadata_edio(nvar,igr,'loss of water from site to watershed discharge','[kg/m2/s]','NA') 
    endif 

!    if (associated(cpoly%metplex_beg_month)) then
!       nvar=nvar+1
!         call vtable_edio_i(cpoly%metplex_beg_month(1),nvar,igr,init,cpoly%siglob_id, &
!         var_len,var_len_global,max_ptrs,'METPLEX_BEG_MONTH :21:hist') 
!       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
!    endif

!    if (associated(cpoly%metplex_beg_year)) then
!       nvar=nvar+1
!         call vtable_edio_i(cpoly%metplex_beg_year(1),nvar,igr,init,cpoly%siglob_id, &
!         var_len,var_len_global,max_ptrs,'METPLEX_BEG_YEAR :21:hist') 
!       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
!    endif

!    if (associated(cpoly%metplex_end_year)) then
!       nvar=nvar+1
!         call vtable_edio_i(cpoly%metplex_end_year(1),nvar,igr,init,cpoly%siglob_id, &
!         var_len,var_len_global,max_ptrs,'METPLEX_END_YEAR :21:hist') 
!       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
!    endif

    if (associated(cpoly%min_monthly_temp)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%min_monthly_temp(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'MIN_MONTHLY_TEMP :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

!    if (associated(cpoly%removed_biomass)) then
!       nvar=nvar+1
!         call vtable_edio_r(cpoly%removed_biomass(1),nvar,igr,init,cpoly%siglob_id, &
!         var_len,var_len_global,max_ptrs,'REMOVED_BIOMASS :21:hist') 
!       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
!    endif 

!    if (associated(cpoly%harvested_biomass)) then
!       nvar=nvar+1
!         call vtable_edio_r(cpoly%harvested_biomass(1),nvar,igr,init,cpoly%siglob_id, &
!         var_len,var_len_global,max_ptrs,'HARVESTED_BIOMASS :21:hist') 
!       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
!    endif 

    if (associated(cpoly%plantation)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%plantation(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PLANTATION_SI :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif 

    if (associated(cpoly%agri_stocking_pft)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%agri_stocking_pft(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGRI_STOCKING_PFT :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%agri_stocking_density)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%agri_stocking_density(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGRI_STOCKING_DENSITY :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%plantation_stocking_pft)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%plantation_stocking_pft(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PLANTATION_STOCKING_PFT :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%plantation_stocking_density)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%plantation_stocking_density(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PLANTATION_STOCKING_DENSITY :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%primary_harvest_memory)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%primary_harvest_memory(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PRIMARY_HARVEST_MEMORY :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%secondary_harvest_memory)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%secondary_harvest_memory(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'SECONDARY_HARVEST_MEMORY:21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%fire_disturbance_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%fire_disturbance_rate(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'FIRE_DISTURBANCE_RATE :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%ignition_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%ignition_rate(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'IGNITION_RATE :21:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%treefall_disturbance_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%treefall_disturbance_rate(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'TREEFALL_DISTURBANCE_RATE :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%nat_disturbance_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%nat_disturbance_rate(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'NAT_DISTURBANCE_RATE :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%nat_dist_type)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%nat_dist_type(1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'NAT_DIST_TYPE :21:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%lai_pft)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%lai_pft(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'LAI_PFT_SI :24:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%wpa_pft)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%wpa_pft(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'WPA_PFT_SI :24:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%wai_pft)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%wai_pft(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'WAI_PFT_SI :24:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
 
    if (associated(cpoly%ntext_soil)) then
       nvar=nvar+1
         call vtable_edio_i(cpoly%ntext_soil(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'NTEXT_SOIL_SI :22:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%lambda_fire)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%lambda_fire(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'LAMBDA_FIRE :29:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%disturbance_memory)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%disturbance_memory(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'DISTURBANCE_MEMORY :255:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%disturbance_rates)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%disturbance_rates(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'DISTURBANCE_RATES_SI :255:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%loss_fraction)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%loss_fraction(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'LOSS_FRACTION :25:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
     
    if (associated(cpoly%green_leaf_factor)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%green_leaf_factor(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'GREEN_LEAF_FACTOR :24:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
        
    if (associated(cpoly%leaf_aging_factor)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%leaf_aging_factor(1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'LEAF_AGING_FACTOR :24:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpoly%basal_area)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%basal_area(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BASAL_AREA_SI :246:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%agb)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%agb(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGB_SI :246:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%pldens)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%pldens(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'PLDENS_SI :246:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%bseeds)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%bseeds(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BSEEDS_SI :246:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%basal_area_growth)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%basal_area_growth(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BASAL_AREA_GROWTH :246:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%agb_growth)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%agb_growth(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGB_GROWTH :246:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%basal_area_mort)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%basal_area_mort(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BASAL_AREA_MORT :246:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%basal_area_cut)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%basal_area_cut(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'BASAL_AREA_CUT :246:year:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%agb_mort)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%agb_mort(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGB_MORT :246:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpoly%agb_cut)) then
       nvar=nvar+1
         call vtable_edio_r(cpoly%agb_cut(1,1,1),nvar,igr,init,cpoly%siglob_id, &
         var_len,var_len_global,max_ptrs,'AGB_CUT :246:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
 
    if(associated(cpoly%dmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%dmean_co2_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL_SI :21:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
    end if

    if(associated(cpoly%dmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%dmean_water_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL_SI :21:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/day]','ipoly') 
    end if

    if(associated(cpoly%dmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%dmean_energy_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL_SI :21:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/day]','ipoly') 
    end if

    if(associated(cpoly%mmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%mmean_co2_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL_SI :21:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
    end if

    if(associated(cpoly%mmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%mmean_water_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL_SI :21:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/day]','ipoly') 
    end if

    if(associated(cpoly%mmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(cpoly%mmean_energy_residual(1),nvar,igr,init,cpoly%siglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL_SI :21:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/day]','ipoly') 
    end if

    if (init == 0) niopoly=nvar-niogrid-nioglobal

    return
  end subroutine filltab_polygontype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_sitetype(igr,ipy,isi,init)

    use ed_var_tables,only:vtable_edio_r,vtable_edio_i,metadata_edio,vtable_edio_d

    implicit none
    
    integer :: init
    integer,intent(in) :: igr,ipy,isi
    integer :: var_len,max_ptrs,var_len_global
    integer :: nvar

    type(sitetype),pointer :: csite

    csite => edgrid_g(igr)%polygon(ipy)%site(isi)

    var_len = csite%npatches
    var_len_global = edgrid_g(igr)%npatches_global
    max_ptrs = get_nsites(edgrid_g(igr))

    !  Set the indexes

    nvar=nioglobal+niogrid+niopoly

    if (associated(csite%paco_id)) then
       nvar=nvar+1
       call vtable_edio_i(csite%paco_id(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'PACO_ID :31:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%paco_n)) then
       nvar=nvar+1
         call vtable_edio_i(csite%paco_n(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'PACO_N :31:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(csite%dist_type)) then
       nvar=nvar+1
         call vtable_edio_i(csite%dist_type(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DIST_TYPE :31:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%age)) then
       nvar=nvar+1
         call vtable_edio_r(csite%age(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AGE :31:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%area)) then
       nvar=nvar+1
         call vtable_edio_r(csite%area(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AREA :31:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%hcapveg)) then
       nvar=nvar+1
         call vtable_edio_r(csite%hcapveg(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'HCAPVEG_PA :31:hist:year') 
       call metadata_edio(nvar,igr,'Total patch heat capacity','[J/m2/K]','NA') 
    endif

    if (associated(csite%fast_soil_C)) then
       nvar=nvar+1
         call vtable_edio_r(csite%fast_soil_C(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'FAST_SOIL_C :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%slow_soil_C)) then
       nvar=nvar+1
         call vtable_edio_r(csite%slow_soil_C(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SLOW_SOIL_C :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%structural_soil_C)) then
       nvar=nvar+1
         call vtable_edio_r(csite%structural_soil_C(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCTURAL_SOIL_C :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%structural_soil_L)) then
       nvar=nvar+1
         call vtable_edio_r(csite%structural_soil_L(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCTURAL_SOIL_L :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mineralized_soil_N)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mineralized_soil_N(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MINERALIZED_SOIL_N :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%fast_soil_N)) then
       nvar=nvar+1
         call vtable_edio_r(csite%fast_soil_N(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'FAST_SOIL_N :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(csite%sum_dgd)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sum_dgd(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SUM_DGD :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%sum_chd)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sum_chd(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SUM_CHD :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%plantation)) then
       nvar=nvar+1
         call vtable_edio_i(csite%plantation(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'PLANTATION :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%can_enthalpy)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_enthalpy(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_ENTHALPY :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air enthalpy','[J/kg]','NA') 
    endif

    if (associated(csite%can_temp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_temp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_TEMP :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air temperature','[K]','NA') 
    endif

    if (associated(csite%can_shv)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_shv(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_SHV :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air specific humidity','[kg_H2O/kg_Air]','NA') 
    endif

    if (associated(csite%can_co2)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_co2(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_CO2 :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air CO2 mixing ratio','[umol/mol]','NA') 
    endif
 
    if (associated(csite%can_rhos)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_rhos(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_RHOS :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air Density','[kg/m3]','NA') 
    endif
 
    if (associated(csite%can_prss)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_prss(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_PRSS :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air Pressure','[Pa]','NA') 
    endif
 
    if (associated(csite%can_theta)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_theta(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_THETA :31:hist') 
       call metadata_edio(nvar,igr,'Canopy air Potential temperature','[K]','NA') 
    endif
   
    if (associated(csite%can_depth)) then
       nvar=nvar+1
         call vtable_edio_r(csite%can_depth(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CAN_DEPTH :31:hist') 
       call metadata_edio(nvar,igr,'Canopy depth','[m]','NA') 
    endif
   
    if (associated(csite%lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(csite%lambda_light(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'LAMBDA_LIGHT :31:hist') 
       call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
    endif
   
    if (associated(csite%dmean_lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(csite%dmean_lambda_light(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LAMBDA_LIGHT :31:hist:dail') 
       call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
    endif
   
    if (associated(csite%mmean_lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mmean_lambda_light(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LAMBDA_LIGHT :31:hist:mont') 
       call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
    endif

    if (associated(csite%lai)) then
       nvar=nvar+1
         call vtable_edio_r(csite%lai(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'LAI_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wpa)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wpa(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WPA_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wai)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wai(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WAI_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
   
    if (associated(csite%sfcwater_mass)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sfcwater_mass(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SFCWATER_MASS :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%sfcwater_energy)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sfcwater_energy(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SFCWATER_ENERGY :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%sfcwater_depth)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sfcwater_depth(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SFCWATER_DEPTH :33:hist:opti') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','m') 
    endif

    if (associated(csite%rshort_s)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_s(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_S :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rshort_s_beam)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_s_beam(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_S_BEAM :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rshort_s_diffuse)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_s_diffuse(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_S_DIFFUSE :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%sfcwater_tempk)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sfcwater_tempk(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SFCWATER_TEMPK :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%sfcwater_fracliq)) then
       nvar=nvar+1
         call vtable_edio_r(csite%sfcwater_fracliq(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SFCWATER_FRACLIQ :33:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%nlev_sfcwater)) then
       nvar=nvar+1
         call vtable_edio_i(csite%nlev_sfcwater(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'NLEV_SFCWATER :31:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ntext_soil)) then
       nvar=nvar+1
         call vtable_edio_i(csite%ntext_soil(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'NTEXT_SOIL_PA :32:hist:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%soil_energy)) then
       nvar=nvar+1
         call vtable_edio_r(csite%soil_energy(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SOIL_ENERGY_PA :32:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%soil_water)) then
       nvar=nvar+1
         call vtable_edio_r(csite%soil_water(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SOIL_WATER_PA :32:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%soil_tempk)) then
       nvar=nvar+1
         call vtable_edio_r(csite%soil_tempk(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SOIL_TEMPK_PA :32:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%soil_fracliq)) then
       nvar=nvar+1
         call vtable_edio_r(csite%soil_fracliq(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SOIL_FRACLIQ_PA :32:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ground_shv)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ground_shv(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'GROUND_SHV :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%surface_ssh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%surface_ssh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SURFACE_SSH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rough)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rough(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'ROUGH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%A_o_max)) then
       nvar=nvar+1
         call vtable_edio_r(csite%A_o_max(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'A_O_MAX :34:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%A_c_max)) then
       nvar=nvar+1
         call vtable_edio_r(csite%A_c_max(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'A_C_MAX :34:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%avg_daily_temp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%avg_daily_temp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AVG_DAILY_TEMP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif  

    if (associated(csite%mean_rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_RH :31:hist') 
       call metadata_edio(nvar,igr,'Heterotrophic respiration','[umol/m2/s]','ipatch') 
    endif

    if (associated(csite%dmean_rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%dmean_rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_RH_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of heterotrophic respiration','[umol/m2/s]','ipatch') 
    endif

    if (associated(csite%mmean_rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mmean_rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_RH_PA :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of heterotrophic respiration','[umol/m2/s]','ipatch') 
    endif

    if (associated(csite%mean_nep)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_nep(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_NEP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_loss2atm)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_loss2atm(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2ATM :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_denseffect)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_denseffect(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_DENSEFFECT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_precipgain)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_precipgain(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_PRECIPGAIN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_loss2runoff)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_loss2runoff(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2RUNOFF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_loss2drainage)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_loss2drainage(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2DRAINAGE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_initialstorage)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_initialstorage(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_INITIALSTORAGE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%wbudget_residual)) then
       nvar=nvar+1
         call vtable_edio_r(csite%wbudget_residual(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'WBUDGET_RESIDUAL :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_latent)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_latent(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_LATENT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_loss2atm)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_loss2atm(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2ATM :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_denseffect)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_denseffect(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_DENSEFFECT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_loss2runoff)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_loss2runoff(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2RUNOFF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_loss2drainage)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_loss2drainage(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2DRAINAGE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_netrad)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_netrad(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_NETRAD :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_precipgain)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_precipgain(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_PRECIPGAIN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_initialstorage)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_initialstorage(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_INITIALSTORAGE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ebudget_residual)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ebudget_residual(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'EBUDGET_RESIDUAL :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_initialstorage)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_initialstorage(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_INITIALSTORAGE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_residual)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_residual(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_RESIDUAL :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_loss2atm)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_loss2atm(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_LOSS2ATM :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_denseffect)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_denseffect(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_DENSEFFECT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_gpp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_GPP :31:hist') 
       call metadata_edio(nvar,igr,'Patch total gross primary productivity per timestep','[umol/m2/dtlsm]','NA') 
    endif

    if (associated(csite%co2budget_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_gpp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_GPP_DBH :36:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_plresp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_plresp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_PLRESP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%co2budget_rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%co2budget_rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CO2BUDGET_RH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%today_A_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%today_A_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_A_DECOMP :31:hist') 
       call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO AT END OF DAY','[NA]','NA') 
    endif

    if (associated(csite%today_Af_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%today_Af_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_AF_DECOMP :31:hist') 
       call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO AT END OF DAY','[NA]','NA') 
    endif

    if (associated(csite%dmean_A_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%dmean_A_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_A_DECOMP :31:hist:dail') 
       call metadata_edio(nvar,igr,'A factor for decomposition','[--]','ipatch') 
    endif

    if (associated(csite%dmean_Af_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%dmean_Af_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_AF_DECOMP :31:hist:dail') 
       call metadata_edio(nvar,igr,'A factor for decomposition, including N','[--]','ipatch') 
    endif

    if (associated(csite%mmean_A_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mmean_A_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_A_DECOMP :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of A factor for decomposition','[--]','ipatch') 
    endif

    if (associated(csite%mmean_Af_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mmean_Af_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_AF_DECOMP :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean A factor for decomposition, including N','[--]','ipatch') 
    endif

    if (associated(csite%repro)) then
       nvar=nvar+1
         call vtable_edio_r(csite%repro(1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'REPRO_PA :34:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%veg_rough)) then
       nvar=nvar+1
         call vtable_edio_r(csite%veg_rough(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_ROUGH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%veg_height)) then
       nvar=nvar+1
         call vtable_edio_r(csite%veg_height (1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_HEIGHT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%fsc_in)) then
       nvar=nvar+1
         call vtable_edio_r(csite%fsc_in(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'FSC_IN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ssc_in)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ssc_in(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SSC_IN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%ssl_in)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ssl_in(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SSL_IN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%fsn_in)) then
       nvar=nvar+1
         call vtable_edio_r(csite%fsn_in(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'FSN_IN :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%total_plant_nitrogen_uptake)) then
       nvar=nvar+1
         call vtable_edio_r(csite%total_plant_nitrogen_uptake(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'TOTAL_PLANT_NITROGEN_UPTAKE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mineralized_N_loss)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mineralized_N_loss(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'NMIN_LOSS :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif


    if (associated(csite%mineralized_N_input)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mineralized_N_input(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'NMIN_INPUT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rshort_g)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_g(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_G :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rshort_g_beam)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_g_beam(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_G_BEAM :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rshort_g_diffuse)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rshort_g_diffuse(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_G_DIFFUSE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_g)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_g(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_G :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_g_surf)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_g_surf(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_G_SURF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_g_incid)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_g_incid(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_G_INCID :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_s)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_s(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_S :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_s_surf)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_s_surf(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_S_SURF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_s_incid)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_s_incid(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_S_INCID :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%albedt)) then
       nvar=nvar+1
         call vtable_edio_r(csite%albedt(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'ALBEDT :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%albedo_beam)) then
       nvar=nvar+1
         call vtable_edio_r(csite%albedo_beam(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'ALBEDO_BEAM :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%albedo_diffuse)) then
       nvar=nvar+1
         call vtable_edio_r(csite%albedo_diffuse(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'ALBEDO_DIFFUSE :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlongup)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlongup(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONGUP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rlong_albedo)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rlong_albedo(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_ALBEDO :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%total_snow_depth)) then
       nvar=nvar+1
         call vtable_edio_r(csite%total_snow_depth(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'TOTAL_SNOW_DEPTH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%snowfac)) then
       nvar=nvar+1
         call vtable_edio_r(csite%snowfac(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'SNOWFAC :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%A_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%A_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'A_DECOMP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%f_decomp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%f_decomp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'F_DECOMP :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'RH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%cwd_rh)) then
       nvar=nvar+1
         call vtable_edio_r(csite%cwd_rh(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CWD_RH :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%fuse_flag)) then
       nvar=nvar+1
         call vtable_edio_i(csite%fuse_flag(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'FUSE_FLAG :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%pft_density_profile)) then
       nvar=nvar+1
         call vtable_edio_r(csite%pft_density_profile(1,1,1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'PFT_DENSITY_PROFILE :346:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%plant_ag_biomass)) then
       nvar=nvar+1
         call vtable_edio_r(csite%plant_ag_biomass(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'PLANT_AG_BIOMASS :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(csite%mean_wflux)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_wflux(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_WFLUX :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mean_latflux)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_latflux(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_LATFLUX :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mean_hflux)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_hflux(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_HFLUX :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mean_runoff)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_runoff(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_RUNOFF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%mean_qrunoff)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mean_qrunoff(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_QRUNOFF :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(csite%htry)) then
       nvar=nvar+1
         call vtable_edio_r(csite%htry(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'HTRY :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(csite%avg_rk4step)) then
       nvar=nvar+1
         call vtable_edio_r(csite%avg_rk4step(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AVG_RK4STEP :31:hist:anal') 
       call metadata_edio(nvar,igr,'Average time step used by Runge-Kutta','[s]','ipatch') 
    endif
    
    if (associated(csite%avg_available_water)) then
       nvar=nvar+1
         call vtable_edio_r(csite%avg_available_water(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AVG_AVAILABLE_WATER_PA :31:hist:anal') 
       call metadata_edio(nvar,igr,'Average available water for transpiration','[kg/m2]','ipatch') 
    endif
    
    if (associated(csite%dmean_rk4step)) then
       nvar=nvar+1
         call vtable_edio_r(csite%dmean_rk4step(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_RK4STEP :31:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean time step used by Runge-Kutta','[s]','ipatch') 
    endif
    
    if (associated(csite%mmean_rk4step)) then
       nvar=nvar+1
         call vtable_edio_r(csite%mmean_rk4step(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_RK4STEP :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean time step used by Runge-Kutta','[s]','ipatch') 
    endif
    
    if (associated(csite%ustar)) then
       nvar=nvar+1
         call vtable_edio_r(csite%ustar(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'USTAR :31:hist') 
         call metadata_edio(nvar,igr,'Patch level friction velocity','[m/s]','ipatch') 
    endif

    if (associated(csite%tstar)) then
       nvar=nvar+1
         call vtable_edio_r(csite%tstar(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'TSTAR :31:hist') 
       call metadata_edio(nvar,igr,'patch level heat transfer atm->canopy','[K]','ipatch') 
    endif

    if (associated(csite%qstar)) then
       nvar=nvar+1
         call vtable_edio_r(csite%qstar(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'QSTAR :31:hist') 
       call metadata_edio(nvar,igr,'patch level vapor transfer atm->canopy','[kg/kg]','ipatch') 
    endif

    if (associated(csite%cstar)) then
       nvar=nvar+1
         call vtable_edio_r(csite%cstar(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'CSTAR :31:hist') 
       call metadata_edio(nvar,igr,'patch level co2 transfer atm->canopy','[ppm?]','ipatch') 
    endif
    
    if (associated(csite%upwp)) then
       nvar=nvar+1
         call vtable_edio_r(csite%upwp(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'UPWP :31:hist') 
       call metadata_edio(nvar,igr,'Vertical flux of U-direction momentum','[kg m^-1 s^-2]','ipatch') 
    endif

    if (associated(csite%tpwp)) then
       nvar=nvar+1
       call vtable_edio_r(csite%tpwp(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'TPWP :31:hist') 
       call metadata_edio(nvar,igr,'Vertical flux of Heat','[kg K m^-2 s^-1]','ipatch')
    endif

    if (associated(csite%qpwp)) then
       nvar=nvar+1
       call vtable_edio_r(csite%qpwp(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'QPWP :31:hist') 
       call metadata_edio(nvar,igr,'Vertical flux of Moisture','[kg m^-2 s% -1]','ipatch')
    endif

    if (associated(csite%cpwp)) then
       nvar=nvar+1
       call vtable_edio_r(csite%cpwp(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'CPWP :31:hist') 
       call metadata_edio(nvar,igr,'Vertical flux of Carbon','[kg umolC mol^-1 m^-2 s-1]','ipatch')
    endif
    
    if (associated(csite%wpwp)) then
       nvar=nvar+1
       call vtable_edio_r(csite%wpwp(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'WPWP :31:hist') 
       call metadata_edio(nvar,igr,'Vertical flux of W-direction momentum','[kg m^-1 s^-2]','ipatch')
    endif
    
    if (associated(csite%avg_carbon_ac)) then
       nvar=nvar+1
         call vtable_edio_r(csite%avg_carbon_ac(1),nvar,igr,init,csite%paglob_id, &
         var_len,var_len_global,max_ptrs,'AVG_CARBON_AC :31:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(csite%old_stoma_vector_max)) then
       nvar=nvar+1
       call vtable_edio_r(csite%old_stoma_vector_max(1,1,1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'OLD_STOMA_VECTOR_MAX :316:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
 
    if(associated(csite%dmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%dmean_co2_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
    end if

    if(associated(csite%dmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%dmean_water_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/day]','ipoly') 
    end if

    if(associated(csite%dmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%dmean_energy_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL_PA :31:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/day]','ipoly') 
    end if

    if(associated(csite%mmean_co2_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%mmean_co2_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL_PA :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
    end if

    if(associated(csite%mmean_water_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%mmean_water_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL_PA :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/day]','ipoly') 
    end if

    if(associated(csite%mmean_energy_residual)) then
       nvar=nvar+1
       call vtable_edio_r(csite%mmean_energy_residual(1),nvar,igr,init,csite%paglob_id, &
            var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL_PA :31:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/day]','ipoly') 
    end if


    if (init == 0) niosite=nvar-niopoly-niogrid-nioglobal

    return
    
  end subroutine filltab_sitetype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_patchtype(igr,ipy,isi,ipa,init)

    use ed_var_tables,only:vtable_edio_r,vtable_edio_i,metadata_edio

    implicit none

    integer :: init
    integer :: ip,np,is,ns,iy,ny
    integer, intent(in) :: igr,ipy,isi,ipa
    integer :: var_len,max_ptrs,var_len_global
    integer :: nvar
    type(patchtype),pointer :: cpatch
    type(sitetype),pointer :: csite
    type(polygontype),pointer :: cpoly
    
    ! -------------------------------------------------------------------------
    ! The condition exists where some patches may have no cohorts.  Therefore
    ! it is possible that the first patch in the site has no cohorts, but
    ! the others might have some cohorts.  Therefore, we want to make sure
    ! we create mapping tables if "any" of the patches have cohorts, so we
    ! will provide for this conditions when we make the mapping tables.
    ! The condition has to arrise where there are no cohorts at all
    ! in any of the patches in any of the polygons on the current machine
    ! to not initialize these vectors. Unlikely, but possible.
    ! -------------------------------------------------------------------------

    if(init==0) then
       ny = edgrid_g(igr)%npolygons

       polyloop: do iy=1,ny
          cpoly => edgrid_g(igr)%polygon(iy)
          ns = cpoly%nsites
          siteloop: do is=1,ns
             csite => cpoly%site(is)
             np = csite%npatches
             patchloop: do ip=1,np
                cpatch => edgrid_g(igr)%polygon(iy)%site(is)%patch(ip)
                if (cpatch%ncohorts > 0) exit polyloop

             end do patchloop
          end do siteloop
       end do polyloop

    else
       cpatch => edgrid_g(igr)%polygon(ipy)%site(isi)%patch(ipa)
    endif

    var_len = cpatch%ncohorts
    var_len_global = edgrid_g(igr)%ncohorts_global
    max_ptrs = edgrid_g(igr)%npatches_global


    if (cpatch%ncohorts == 0) return
    
    nvar=nioglobal+niogrid+niopoly+niosite

    if (associated(cpatch%pft)) then
       nvar=nvar+1
         call vtable_edio_i(cpatch%pft(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PFT :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Plant Functional Type','[-]','NA') 
    endif

    if (associated(cpatch%nplant)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%nplant(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'NPLANT :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Plant density','[plant/m2]','NA') 
    endif

    if (associated(cpatch%hite)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%hite(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'HITE :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%agb)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%agb(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'AGB_CO :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Above-ground biomass','[kgC/plant]','icohort') 
    endif

    if (associated(cpatch%basarea)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%basarea(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BA_CO :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Basal-area','[cm2]','icohort') 
    endif

    if (associated(cpatch%dagb_dt)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dagb_dt(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DAGB_DT :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Above-ground biomass growth','[kgC/plant/yr]','icohort') 
    endif

    if (associated(cpatch%dba_dt)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dba_dt(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DBA_DT :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'Basal-area growth','[cm2/plant/yr]','icohort') 
    endif

    if (associated(cpatch%ddbh_dt)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%ddbh_dt(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DDBH_DT :41:hist:anal:mont:year') 
       call metadata_edio(nvar,igr,'DBH growth','[cm/plant/yr]','icohort') 
    endif

    if (associated(cpatch%dbh)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dbh(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DBH :41:hist:anal:year:mont') 
       call metadata_edio(nvar,igr,'Diameter at breast height','[cm]','icohort') 
    endif

    if (associated(cpatch%bdead)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%bdead(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BDEAD :41:hist:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%bleaf)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%bleaf(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BLEAF :41:hist:year:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%phenology_status)) then
       nvar=nvar+1
         call vtable_edio_i(cpatch%phenology_status(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PHENOLOGY_STATUS :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%balive)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%balive(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BALIVE :41:hist:year:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%broot)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%broot(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BROOT :41:hist:year:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%bsapwood)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%bsapwood(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BSAPWOOD :41:hist:year:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%lai)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%lai(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LAI_CO :41:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%wpa)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%wpa(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'WPA_CO :41:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%wai)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%wai(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'WAI_CO :41:hist:dail:mont:year') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%bstorage)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%bstorage(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BSTORAGE :41:hist:year:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%cb)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%cb(1,1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'CB :49:hist:mont:year') 
       call metadata_edio(nvar,igr,'carbon balance previous 12 months+current','[kgC/plant]','13 - icohort') 
    endif

    if (associated(cpatch%cb_max)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%cb_max(1,1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'CB_MAX :49:hist:mont:year') 
       call metadata_edio(nvar,igr,'TOC carbon balance previous 12 months+current','[kgC/plant]','13 - icohort') 
    endif

    if (associated(cpatch%cbr_bar)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%cbr_bar(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'CBR_BAR :41:hist:mont:year') 
       call metadata_edio(nvar,igr,'Annual average ratio of cb/cb_max','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_cb)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_cb(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_CB :41:hist:mont:year') 
       call metadata_edio(nvar,igr,'Monthly mean of carbon balance','[kgC/plant/yr]','NA') 
    endif

    if (associated(cpatch%veg_energy)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%veg_energy(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_ENERGY :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%veg_temp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%veg_temp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_TEMP :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%veg_fliq)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%veg_fliq(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_FLIQ :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%veg_water)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%veg_water(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VEG_WATER :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mean_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mean_gpp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_GPP :41:hist:anal') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mean_leaf_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mean_leaf_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_LEAF_RESP :41:hist:anal') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mean_root_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mean_root_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MEAN_ROOT_RESP :41:hist:anal') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%today_leaf_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%today_leaf_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_LEAF_RESP :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%today_root_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%today_root_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_ROOT_RESP :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%today_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%today_gpp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_GPP :41:hist') 
       call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO PRIOR TO DAILY WRITE OUT','[umol/m2/s]','icohort') 
    endif

    if (associated(cpatch%today_gpp_pot)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%today_gpp_pot(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_GPP_POT :41:hist') 
       call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO PRIOR TO DAILY WRITE OUT','[NA]','NA') 
    endif

    if (associated(cpatch%today_gpp_max)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%today_gpp_max(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TODAY_GPP_MAX :41:hist') 
       call metadata_edio(nvar,igr,'NOT A DIANOSTIC-WILL ZERO PRIOR TO DAILY WRITE OUT','[NA]','NA') 
    endif

    if (associated(cpatch%growth_respiration)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%growth_respiration(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'GROWTH_RESPIRATION :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%storage_respiration)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%storage_respiration(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'STORAGE_RESPIRATION :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%vleaf_respiration)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%vleaf_respiration(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VLEAF_RESPIRATION :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif  

    if (associated(cpatch%dmean_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_gpp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_GPP_CO :41:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean Gross Primary Productivity','[kgC/plant/yr]','icohort') 
    endif

    if (associated(cpatch%dmean_leaf_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_leaf_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LEAF_RESP_CO :41:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean leaf respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%dmean_root_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_root_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_ROOT_RESP_CO :41:hist:dail') 
       call metadata_edio(nvar,igr,'Daily mean root respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%mmean_gpp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_gpp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_GPP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean Gross Primary Productivity','[kgC/plant/yr]','icohort') 
    endif

    if (associated(cpatch%mmean_leaf_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_leaf_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LEAF_RESP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mstructural_growthean leaf respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%mmean_root_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_root_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_ROOT_RESP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean root respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%mmean_growth_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_growth_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_GROWTH_RESP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean growth respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%mmean_storage_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_storage_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_STORAGE_RESP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean storage respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%mmean_vleaf_resp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_vleaf_resp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_VLEAF_RESP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean virtual leaf respiration','[kgC/plant/yr]','icohort') 
    end if

    if (associated(cpatch%fsn)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%fsn(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'FSN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif
    
    if (associated(cpatch%old_stoma_vector)) then
       nvar=nvar+1
       call vtable_edio_r(cpatch%old_stoma_vector(1,1),nvar,igr,init,cpatch%coglob_id, &
            var_len,var_len_global,max_ptrs,'OLD_STOMA_VECTOR :416:hist')
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
    endif

    if (associated(cpatch%monthly_dndt)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%monthly_dndt(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MONTHLY_DNDT :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mort_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mort_rate(1,1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MORT_RATE_CO :48:hist:anal:dail') 
       call metadata_edio(nvar,igr,'Mortality rates','[1/yr]','icohort') 
    endif

    if (associated(cpatch%mmean_mort_rate)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_mort_rate(1,1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_MORT_RATE_CO :48:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean mortality rate','[1/yr]','icohort') 
    endif

    if (associated(cpatch%Psi_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%Psi_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PSI_OPEN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%krdepth)) then
       nvar=nvar+1
         call vtable_edio_i(cpatch%krdepth(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'KRDEPTH :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%first_census)) then
       nvar=nvar+1
         call vtable_edio_i(cpatch%first_census(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'FIRST_CENSUS :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%new_recruit_flag)) then
       nvar=nvar+1
         call vtable_edio_i(cpatch%new_recruit_flag(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'NEW_RECRUIT_FLAG :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%light_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%light_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LIGHT_LEVEL :41:hist') 
       call metadata_edio(nvar,igr,'Relative light level','[NA]','icohort') 
    endif

    if (associated(cpatch%light_level_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%light_level_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LIGHT_LEVEL_BEAM :41:hist') 
       call metadata_edio(nvar,igr,'Relative light level, beam fraction','[NA]','icohort') 
    endif

    if (associated(cpatch%light_level_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%light_level_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LIGHT_LEVEL_DIFF :41:hist') 
       call metadata_edio(nvar,igr,'Relative light level, diffuse fraction','[NA]','icohort') 
    endif

    if (associated(cpatch%beamext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%light_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BEAMEXT_LEVEL :41:hist') 
       call metadata_edio(nvar,igr,'Beam extinction level','[NA]','icohort') 
    endif

    if (associated(cpatch%diffext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%light_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DIFFEXT_LEVEL :41:hist') 
       call metadata_edio(nvar,igr,'diff extinction level','[NA]','icohort') 
    endif

    if (associated(cpatch%norm_par_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%norm_par_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'NORM_PAR_BEAM :41:hist') 
       call metadata_edio(nvar,igr,'Relative light level, beam fraction','[NA]','icohort') 
    endif

    if (associated(cpatch%norm_par_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%norm_par_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'NORM_PAR_DIFF :41:hist') 
       call metadata_edio(nvar,igr,'Relative light level, diffuse fraction','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_light_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_light_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LIGHT_LEVEL :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of Relative light level ','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_light_level_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_light_level_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LIGHT_LEVEL_BEAM :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of Relative light level (beam)','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_light_level_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_light_level_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LIGHT_LEVEL_DIFF :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of Relative light level (diffuse)','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_beamext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_beamext_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_BEAMEXT_LEVEL :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of beam extinction level ','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_diffext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_diffext_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_DIFFEXT_LEVEL :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of diff extinction level ','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_norm_par_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_norm_par_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_NORM_PAR_BEAM :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of Relative light level (beam)','[NA]','icohort') 
    endif

    if (associated(cpatch%dmean_norm_par_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_norm_par_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_NORM_PAR_DIFF :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of Relative light level (diffuse)','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_light_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_light_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LIGHT_LEVEL :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of Relative light level ','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_light_level_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_light_level_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LIGHT_LEVEL_BEAM :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of Relative light level (beam)','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_light_level_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_light_level_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LIGHT_LEVEL_DIFF :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of Relative light level (diff)','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_beamext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_beamext_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_BEAMEXT_LEVEL :41:hist:mont') 
       call metadata_edio(nvar,igr,'Diurnal mean of beam extinction level ','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_diffext_level)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_diffext_level(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_DIFFEXT_LEVEL :41:hist:mont') 
       call metadata_edio(nvar,igr,'Diurnal mean of diff extinction level ','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_norm_par_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_norm_par_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_NORM_PAR_BEAM :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of Relative light level (beam)','[NA]','icohort') 
    endif

    if (associated(cpatch%mmean_norm_par_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_norm_par_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_NORM_PAR_DIFF :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of Relative light level (diff)','[NA]','icohort') 
    endif

    if (associated(cpatch%lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%lambda_light(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LAMBDA_LIGHT_CO :41:hist') 
       call metadata_edio(nvar,igr,'Light extinction','[m2/m2]','icohort') 
    endif

    if (associated(cpatch%dmean_lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_lambda_light(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_LAMBDA_LIGHT_CO :41:hist:dail') 
       call metadata_edio(nvar,igr,'Diurnal mean of light extinction ','[m2/m2]','icohort') 
    endif

    if (associated(cpatch%mmean_lambda_light)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_lambda_light(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LAMBDA_LIGHT_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'Monthly mean of light extinction ','[m2/m2]','icohort') 
    endif

    if (associated(cpatch%par_v)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%par_v(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PAR_V :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%par_v_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%par_v_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PAR_V_BEAM :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%par_v_diffuse)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%par_v_diffuse(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PAR_V_DIFFUSE :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_par_v)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_par_v(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_PAR_V :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_par_v_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_par_v_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_PAR_V_BEAM :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_par_v_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_par_v_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_PAR_V_DIFF :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_par_v)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_par_v(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_PAR_V :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_par_v_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_par_v_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_PAR_V_BEAM :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_par_v_diff)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_par_v_diff(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_PAR_V_DIFF :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rshort_v)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rshort_v(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_V :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rshort_v_beam)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rshort_v_beam(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_V_BEAM :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rshort_v_diffuse)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rshort_v_diffuse(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RSHORT_V_DIFFUSE :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rlong_v)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rlong_v(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_V :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rlong_v_surf)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rlong_v_surf(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_V_SURF :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rlong_v_incid)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rlong_v_incid(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RLONG_V_INCID :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rb)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rb(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RB :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%A_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%A_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'A_OPEN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%A_closed)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%A_closed(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'A_CLOSED :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%Psi_closed)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%Psi_closed(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PSI_CLOSED :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rsw_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rsw_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RSW_OPEN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%rsw_closed)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%rsw_closed(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'RSW_CLOSED :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%fsw)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%fsw(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'FSW :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%fs_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%fs_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'FS_OPEN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_fs_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_fs_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_FS_OPEN_CO :41:dail:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_fsw)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_fsw(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_FSW_CO :41:dail:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%dmean_fsn)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%dmean_fsn(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'DMEAN_FSN_CO :41:dail:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_fs_open)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_fs_open(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_FS_OPEN_CO :41:mont:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_fsw)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_fsw(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_FSW_CO :41:mont:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_fsn)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_fsn(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_FSN_CO :41:mont:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%stomatal_resistance)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%stomatal_resistance(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'STOMATAL_RESISTANCE :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%leaf_maintenance)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%leaf_maintenance(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LEAF_MAINTENANCE :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%root_maintenance)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%root_maintenance(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'ROOT_MAINTENANCE :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_leaf_maintenance)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_leaf_maintenance(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LEAF_MAINTENANCE :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_root_maintenance)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_root_maintenance(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_ROOT_MAINTENANCE :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%leaf_drop)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%leaf_drop(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LEAF_DROP :41:hist:dail') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%mmean_leaf_drop)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%mmean_leaf_drop(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'MMEAN_LEAF_DROP_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%bseeds)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%bseeds(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'BSEEDS_CO :41:hist:mont') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%leaf_respiration)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%leaf_respiration(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LEAF_RESPIRATION :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%root_respiration)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%root_respiration(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'ROOT_RESPIRATION :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%hcapveg)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%hcapveg(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'HCAPVEG :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%gpp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%gpp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'GPP :41:hist') 
       call metadata_edio(nvar,igr,'Gross Primary Production','[umol/m2/s]','icohort') 
    endif

    if (associated(cpatch%paw_avg)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%paw_avg(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'PAW_AVG :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%turnover_amp)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%turnover_amp(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'TURNOVER_AMP :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%llspan)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%llspan(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'LLSPAN :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%vm_bar)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%vm_bar(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'VM_BAR :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    if (associated(cpatch%sla)) then
       nvar=nvar+1
         call vtable_edio_r(cpatch%sla(1),nvar,igr,init,cpatch%coglob_id, &
         var_len,var_len_global,max_ptrs,'SLA :41:hist') 
       call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    endif

    return
  end subroutine filltab_patchtype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
                   ! ========= UTILITITY FUNCTIONS =========== !
!==========================================================================================!
!==========================================================================================!
  function get_nsites(cgrid)

    implicit none
    integer :: get_nsites
    integer :: ipy,isi
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly

    get_nsites = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          get_nsites = get_nsites + 1
       enddo
    enddo

    return
  end function get_nsites
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  function get_npatches(cgrid)

    implicit none
    integer :: get_npatches
    integer :: ipy,isi,ipa
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly
    type(sitetype),pointer        :: csite

    get_npatches = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          csite=>cpoly%site(isi)
          do ipa=1,csite%npatches
             get_npatches = get_npatches + 1
          enddo
       enddo
    enddo
    return
  end function get_npatches
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  function get_ncohorts(cgrid)

    implicit none
    integer :: get_ncohorts
    integer :: ipy,isi,ipa,ico
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly
    type(sitetype),pointer        :: csite
    type(patchtype),pointer       :: cpatch

    get_ncohorts = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          csite=>cpoly%site(isi)
          do ipa=1,csite%npatches
             cpatch=>csite%patch(ipa)
             get_ncohorts = get_ncohorts + cpatch%ncohorts
          enddo
       enddo
    enddo
    return
  end function get_ncohorts
!==========================================================================================!
!==========================================================================================!
end module ed_state_vars




