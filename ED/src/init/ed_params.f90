!==========================================================================================!
!==========================================================================================!
!     This is the main loader of ecosystem parameters.  Since some compilers do not under- !
! stand the assignment in the modules when the variable is not a constant (parameter),     !
! this is the safest way to guarantee it will read something (not to mention that makes    !
! compilation much faster when you want to test the sensitivity of one number).            !
!------------------------------------------------------------------------------------------!
subroutine load_ed_ecosystem_params()

   use ed_max_dims , only : n_pft               ! ! intent(in)
   use pft_coms    , only : include_these_pft   & ! intent(in)
                          , include_pft         & ! intent(out)
                          , include_pft_ag      & ! intent(out)
                          , C2B                 & ! intent(out)
                          , frost_mort          & ! intent(out)
                          , grass_pft           & ! intent(out)
                          , pft_name16          ! ! intent(out)
   use disturb_coms, only : ianth_disturb       ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer :: p
   !---------------------------------------------------------------------------------------!

   !----- Loading several parameters ------------------------------------------------------!
   call init_decomp_params()
   call init_ff_coms()
   call init_disturb_params()
   call init_physiology_params()
   call init_met_params()
   call init_lapse_params()
   call init_can_rad_params()
   call init_hydro_coms()
   call init_soil_coms()
   call init_phen_coms()
   call init_ed_misc_coms()

   !---------------------------------------------------------------------------------------!
   !      Main table of Plant functional types.  If you add some PFT, please make sure     !
   ! that you assign values for all PFT-dependent variables.  Below is a summary table of  !
   ! the main characteristics of the currently available PFTs.                             !
   !---------------------------------------------------------------------------------------!
   !  PFT | Name                                       | Grass   | Tropical | agriculture? !
   !------+--------------------------------------------+---------+----------+--------------!
   !    1 | C4 grass                                   |     yes |      yes |          yes !
   !    2 | Early tropical                             |      no |      yes |           no !
   !    3 | Mid tropical                               |      no |      yes |           no !
   !    4 | Late tropical                              |      no |      yes |           no !
   !    5 | C3 grass                                   |     yes |       no |          yes !
   !    6 | Northern pines                             |      no |       no |           no !
   !    7 | Southern pines                             |      no |       no |           no !
   !    8 | Late conifers                              |      no |       no |           no !
   !    9 | Early temperate deciduous                  |      no |       no |           no !
   !   10 | Mid temperate deciduous                    |      no |       no |           no !
   !   11 | Late temperate deciduous                   |      no |       no |           no !
   !   12 | C3 pasture                                 |     yes |       no |          yes !
   !   13 | C3 crop (e.g.,wheat, rice, soybean)        |     yes |       no |          yes !
   !   14 | C4 pasture                                 |     yes |      yes |          yes !
   !   15 | C4 crop (e.g.,corn/maize)                  |     yes |      yes |          yes !
   !   16 | Subtropical C3 grass                       |     yes |      yes |          yes !
   !   17 | Araucaria                                  |      no |      yes |           no !
   !------+--------------------------------------------+---------+----------+--------------!

   !----- Name the PFTs (no spaces, please). ----------------------------------------------!
   pft_name16( 1) = 'C4_grass        '
   pft_name16( 2) = 'Early_tropical  '
   pft_name16( 3) = 'Mid_tropical    '
   pft_name16( 4) = 'Late_tropical   '
   pft_name16( 5) = 'C3_grass        '
   pft_name16( 6) = 'North_pine      '
   pft_name16( 7) = 'South_pine      '
   pft_name16( 8) = 'Late_conifer    '
   pft_name16( 9) = 'Early_hardwood  '
   pft_name16(10) = 'Mid_hardwood    '
   pft_name16(11) = 'Late_hardwood   '
   pft_name16(12) = 'C3_pasture      '
   pft_name16(13) = 'C3_crop         '
   pft_name16(14) = 'C4_pasture      '
   pft_name16(15) = 'C4_crop         '
   pft_name16(16) = 'Subtrop_C3_grass'
   pft_name16(17) = 'Araucaria       '

   !----- Define the grass PFTs -----------------------------------------------------------!
   grass_pft=huge(1)
   grass_pft(1)=1
   grass_pft(2)=5
   grass_pft(3)=12
   grass_pft(4)=13
   grass_pft(5)=14
   grass_pft(6)=15
   grass_pft(7)=16

   !---------------------------------------------------------------------------------------!
   !    Include_pft: flag specifying to whether you want to include a plant functional     !
   !                 type (1) or whether you want it excluded (0) from the simulation.     !
   !---------------------------------------------------------------------------------------!
   include_pft = 0
   include_pft_ag = 0
   do p=1,n_pft
      if (include_these_pft(p) > 0 .and. include_these_pft(p) <= n_pft) then
         include_pft(include_these_pft(p)) = 1
      end if
   end do

   !----- Grasses can grow anywhere, including agricultural patches -----------------------!
   p=1
   do while (grass_pft(p) > 0 .and. grass_pft(p) <= n_pft)
      if (include_pft(grass_pft(p)) == 1) include_pft_ag(grass_pft(p)) = 1
      p = p+1
   end do
   if (sum(include_pft_ag) == 0 .and. ianth_disturb == 1) then
!!      call fatal_error ('No grass included in include_these_pft,'//&
      call warning ('No grass included in include_these_pft,'//&
                       &' you should have at least one kind of grass...'                   &
                       ,'load_ecosystem_params','ed_params.f90')
   end if

   !----- Assign many PFT-dependent parameters. -------------------------------------------!
   call init_pft_photo_params()
   call init_pft_resp_params()
   call init_pft_alloc_params()
   call init_pft_mort_params()
   call init_pft_nitro_params()
   call init_pft_leaf_params()
   call init_pft_repro_params()
   call init_pft_derived_params()

   !----- This must be done after defining some PFT parameters. ---------------------------!
   call init_can_air_params()

   !---------------------------------------------------------------------------------------!
   !     This should be always the last one, since it depends on variables assigned in     !
   ! the previous init_????_params.                                                        !
   !---------------------------------------------------------------------------------------!
   call init_rk4_params()
   !---------------------------------------------------------------------------------------!

   return
end subroutine load_ed_ecosystem_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns values for some variables that are in ed_misc_coms, which    !
! wouldn't fit in any of the other categories.                                             !
!------------------------------------------------------------------------------------------!
subroutine init_ed_misc_coms
   use ed_max_dims  , only : n_pft               & ! intent(in)
                           , n_dbh               & ! intent(in)
                           , n_age               ! ! intent(in)
   use consts_coms  , only : erad                & ! intent(in)
                           , pio180              ! ! intent(in)
   use ed_misc_coms , only : burnin              & ! intent(out)
                           , outputMonth         & ! intent(out)
                           , restart_target_year & ! intent(out)
                           , use_target_year     & ! intent(out)
                           , maxage              & ! intent(out)
                           , dagei               & ! intent(out)
                           , maxdbh              & ! intent(out)
                           , ddbhi               & ! intent(out)
                           , vary_elev           & ! intent(out)
                           , vary_hyd            & ! intent(out)
                           , vary_rad            & ! intent(out)
                           , max_thsums_dist     & ! intent(out)
                           , max_poihist_dist    & ! intent(out)
                           , max_poi99_dist      ! ! intent(out)
   implicit none


   !----- Flags that allow components of subgrid heterogeneity to be turned on/off --------!
   vary_elev = 1
   vary_rad = 1
   vary_hyd = 1

   !----- Number of years to ignore demography when starting a run. -----------------------!
   burnin = 0

   !----- Month to output the yearly files. -----------------------------------------------!
   outputMonth = 6

   !----- Year to read when parsing pss/css with multiple years. --------------------------!
   restart_target_year = 2000

   !----- Flag specifying whether to search for a target year in pss/css. -----------------!
   use_target_year = 0    

   !----- Maximum age [yr] to split into classes. -----------------------------------------!
   maxage = 200.

   !----- Maximum DBH [cm] to be split into classes. --------------------------------------!
   maxdbh = 100.

   !---------------------------------------------------------------------------------------!
   !     The inverse of bin classes will depend on max??? and n_???, leaving one class for !
   ! when the number exceeds the maximum.                                                  !
   !---------------------------------------------------------------------------------------!
   dagei = real(n_age-1) / maxage
   ddbhi = real(n_dbh-1) / maxdbh

   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon.  The value below is 1.25 degree at the Equator.  !
   !---------------------------------------------------------------------------------------!
   max_thsums_dist     = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Alternative method for mixing 1 grid and POI's.  Only use the grid if their is   !
   ! NOT an POI  within a user specified resolution.  Remember, this assumes there is only !
   ! 1 gridded file, and it is the first file when ied_init_mode is set to 99  (Developer  !
   ! use only).                                                                            !
   !---------------------------------------------------------------------------------------!
   max_poi99_dist     = 5.0 * erad * pio180
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      This variable is used for the history start initialisation.  This sets the       !
   ! maximum acceptable distance between the expected polygon and the polygon found in the !
   ! history file.  Units: m.                                                              !
   !---------------------------------------------------------------------------------------!
   max_poihist_dist   = 250.
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ed_misc_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines the minimum and maximum acceptable values in the meteoro-    !
! logical forcing.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine init_met_params()

   use met_driver_coms, only : rshort_min   & ! intent(out)
                             , rshort_max   & ! intent(out)
                             , rlong_min    & ! intent(out)
                             , rlong_max    & ! intent(out)
                             , atm_tmp_min  & ! intent(out)
                             , atm_tmp_max  & ! intent(out)
                             , atm_shv_min  & ! intent(out)
                             , atm_shv_max  & ! intent(out)
                             , atm_rhv_min  & ! intent(out)
                             , atm_rhv_max  & ! intent(out)
                             , atm_co2_min  & ! intent(out)
                             , atm_co2_max  & ! intent(out)
                             , prss_min     & ! intent(out)
                             , prss_max     & ! intent(out)
                             , pcpg_min     & ! intent(out)
                             , pcpg_max     & ! intent(out)
                             , vels_min     & ! intent(out)
                             , vels_max     & ! intent(out)
                             , geoht_min    & ! intent(out)
                             , geoht_max    ! ! intent(out)

   !----- Minimum and maximum acceptable shortwave radiation [W/m²]. ----------------------!
   rshort_min  = 0.
   rshort_max  = 1400.
   !----- Minimum and maximum acceptable longwave radiation [W/m²]. -----------------------!
   rlong_min   = 40.
   rlong_max   = 600.
   !----- Minimum and maximum acceptable air temperature    [   K]. -----------------------!
   atm_tmp_min = 184.     ! Lowest temperature ever measured, in Vostok Basin, Antarctica
   atm_tmp_max = 331.     ! Highest temperature ever measured, in El Azizia, Libya
   !----- Minimum and maximum acceptable air specific humidity [kg_H2O/kg_air]. -----------!
   atm_shv_min = 1.e-6    ! That corresponds to a relative humidity of 0.1% at 1000hPa
   atm_shv_max = 3.2e-2   ! That corresponds to a dew point of 32°C at 1000hPa.
   !----- Minimum and maximum acceptable CO2 mixing ratio [µmol/mol]. ---------------------!
   atm_co2_min = 100.     ! 
   atm_co2_max = 1100.    ! 
   !----- Minimum and maximum acceptable pressure [Pa]. -----------------------------------!
   prss_min =  45000. ! It may crash if you run a simulation in Mt. Everest.
   prss_max = 110000. ! It may crash if you run a simulation under water.
   !----- Minimum and maximum acceptable precipitation rates [kg/m²/s]. -------------------!
   pcpg_min     = 0.0     ! No negative precipitation is allowed
   pcpg_max     = 0.1111  ! This is a precipitation rate of 400mm/hr.
   !----- Minimum and maximum acceptable wind speed [m/s]. --------------------------------!
   vels_min     =  0.0    ! No negative wind is acceptable.
   vels_max     = 85.0    ! Maximum sustained winds recorded during Typhoon Tip (1970).
   !----- Minimum and maximum reference heights [m]. --------------------------------------!
   geoht_min    =   1.0   ! This should be above-canopy measurement, but 1.0 is okay for
                          !     grasslands...
   geoht_max    = 350.0   ! This should be not that much above the canopy, but tall towers
                          !     do exist...
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Minimum and maximum acceptable relative humidity (fraction).  This is not going   !
   ! cause the simulation to crash, instead it will just impose these numbers to the       !
   ! meteorological forcing.                                                               !
   !---------------------------------------------------------------------------------------!
   atm_rhv_min = 5.e-3 ! 0.5%
   atm_rhv_max = 1.0   ! 100.0%.  Although canopy air space can experience super-
                       !    saturation, we don't allow the air above to be super-saturated.
   return
end subroutine init_met_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines defaults for lapse rates.                                    !
!------------------------------------------------------------------------------------------!
subroutine init_lapse_params()

   use met_driver_coms, only : lapse             & ! intent(out)
                             , atm_tmp_intercept & ! intent(out)
                             , atm_tmp_slope     & ! intent(out)
                             , prec_intercept    & ! intent(out)
                             , prec_slope        & ! intent(out)
                             , humid_scenario    ! ! intent(out)

   lapse%geoht        = 0.0
   lapse%vels         = 0.0
   lapse%atm_tmp      = 0.0
   lapse%atm_theta    = 0.0
   lapse%atm_theiv    = 0.0
   lapse%atm_shv      = 0.0
   lapse%prss         = 0.0
   lapse%pcpg         = 0.0
   lapse%atm_co2      = 0.0
   lapse%rlong        = 0.0
   lapse%nir_beam     = 0.0
   lapse%nir_diffuse  = 0.0
   lapse%par_beam     = 0.0
   lapse%par_diffuse  = 0.0
   lapse%pptnorm      = 0.0

   atm_tmp_intercept = 0.0
   atm_tmp_slope     = 1.0
   prec_intercept    = 0.0
   prec_slope        = 1.0
   humid_scenario    = 0

   return
end subroutine init_lapse_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some radiation related parameters.                        !
!------------------------------------------------------------------------------------------!
subroutine init_can_rad_params()

   use canopy_radiation_coms , only : leaf_reflect_nir            & ! intent(out)
                                    , leaf_trans_nir              & ! intent(out)
                                    , leaf_scatter_nir            & ! intent(out)
                                    , leaf_reflect_vis_temperate  & ! intent(out)
                                    , leaf_trans_vis_temperate    & ! intent(out)
                                    , leaf_scatter_vis            & ! intent(out)
                                    , leaf_reflect_vis_tropics    & ! intent(out)
                                    , leaf_trans_vis_tropics      & ! intent(out)
                                    , diffuse_backscatter_vis     & ! intent(out)
                                    , diffuse_backscatter_nir     & ! intent(out)
                                    , emis_v                      & ! intent(out)
                                    , mubar                       & ! intent(out)
                                    , visible_fraction            & ! intent(out)
                                    , visible_fraction_dir        & ! intent(out)
                                    , visible_fraction_dif        & ! intent(out)
                                    , leaf_reflect_nir            & ! intent(out)
                                    , leaf_trans_nir              & ! intent(out)
                                    , rshort_twilight_min         & ! intent(out)
                                    , cosz_min                    & ! intent(out)
                                    , cosz_min8                   ! ! intent(out)
   use ed_max_dims              , only : n_pft                    ! ! intent(out)
   use pft_coms              , only : phenology                   ! ! intent(out)

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   real :: leaf_scatter_vis_temperate
   real :: leaf_scatter_vis_tropics
   real :: diffuse_bscat_vis_temp
   real :: diffuse_bscat_vis_trop
   !---------------------------------------------------------------------------------------!

   mubar                      = 1.0d0 

   visible_fraction           = 0.45
   visible_fraction_dir       = 0.43
   visible_fraction_dif       = 0.52
   leaf_reflect_nir           = 0.577
   leaf_trans_nir             = 0.248

   leaf_scatter_nir           = leaf_reflect_nir + leaf_trans_nir

   leaf_scatter_vis_temperate = leaf_reflect_vis_temperate + leaf_trans_vis_temperate

   leaf_scatter_vis_tropics   = leaf_reflect_vis_tropics   + leaf_trans_vis_tropics

   diffuse_bscat_vis_temp  = (2.0 * leaf_reflect_vis_temperate - leaf_trans_vis_temperate) &
                           / (3.0 * leaf_scatter_vis_temperate)

   diffuse_bscat_vis_trop  = (2.0 * leaf_reflect_vis_tropics   - leaf_trans_vis_tropics)   &
                           / (3.0 * leaf_scatter_vis_tropics)

   diffuse_backscatter_nir = (2.0 * leaf_reflect_nir - leaf_trans_nir)                     &
                           / (3.0 * leaf_scatter_nir)

   leaf_scatter_vis(1:4)   = leaf_scatter_vis_tropics
   leaf_scatter_vis(5:11)  = leaf_scatter_vis_temperate
   leaf_scatter_vis(12:13) = leaf_scatter_vis_temperate
   leaf_scatter_vis(14:15) = leaf_scatter_vis_tropics
   leaf_scatter_vis(16)    = leaf_scatter_vis_tropics
   leaf_scatter_vis(17)    = leaf_scatter_vis_temperate

   diffuse_backscatter_vis(1:4)   = diffuse_bscat_vis_trop
   diffuse_backscatter_vis(5:11)  = diffuse_bscat_vis_temp
   diffuse_backscatter_vis(12:13) = diffuse_bscat_vis_temp
   diffuse_backscatter_vis(14:15) = diffuse_bscat_vis_trop
   diffuse_backscatter_vis(16)    = diffuse_bscat_vis_trop
   diffuse_backscatter_vis(17)    = diffuse_bscat_vis_temp

   emis_v(1)     = 9.60d-1
   emis_v(2:4)   = 9.50d-1
   emis_v(5)     = 9.60d-1
   emis_v(6:8)   = 9.70d-1
   emis_v(9:11)  = 9.50d-1
   emis_v(12:15) = 9.60d-1
   emis_v(16)    = 9.60d-1
   emis_v(17)    = 9.70d-1

   !---------------------------------------------------------------------------------------!
   !     These variables are the thresholds for things that should be computed during the  !
   ! day time hours only.                                                                  !
   !---------------------------------------------------------------------------------------!
   rshort_twilight_min = 0.5
   cosz_min            = 0.03
   cosz_min8           = dble(cosz_min)
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_can_rad_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some canopy air related parameters.                       !
!------------------------------------------------------------------------------------------!
subroutine init_can_air_params()
   use consts_coms    , only : onethird              & ! intent(in)
                             , twothirds             & ! intent(in)
                             , vonk                  ! ! intent(in)
   use pft_coms       , only : hgt_min               & ! intent(in)
                             , hgt_max               ! ! intent(in)
   use canopy_air_coms, only : i_blyr_condct         & ! intent(in)
                             , isfclyrm              & ! intent(in)
                             , ustmin                & ! intent(in)
                             , ggfact                & ! intent(in)
                             , gamm                  & ! intent(in)
                             , gamh                  & ! intent(in)
                             , tprandtl              & ! intent(in)
                             , vkopr                 & ! intent(in)
                             , vh2vr                 & ! intent(in)
                             , vh2dh                 & ! intent(in)
                             , dry_veg_lwater        & ! intent(out)
                             , fullveg_lwater        & ! intent(out)
                             , rb_inter              & ! intent(out)
                             , rb_slope              & ! intent(out)
                             , veg_height_min        & ! intent(out)
                             , minimum_canopy_depth  & ! intent(out)
                             , minimum_canopy_depth8 & ! intent(out)
                             , exar                  & ! intent(out)
                             , covr                  & ! intent(out)
                             , ugbmin                & ! intent(out)
                             , ubmin                 & ! intent(out)
                             , exar8                 & ! intent(out)
                             , ez                    & ! intent(out)
                             , ustmin8               & ! intent(out)
                             , ugbmin8               & ! intent(out)
                             , ubmin8                & ! intent(out)
                             , ggfact8               & ! intent(out)
                             , ez8                   & ! intent(out)
                             , vh2dh8                & ! intent(out)
                             , ncanmax               & ! intent(out)
                             , dz_m97                & ! intent(out)
                             , cdrag0                & ! intent(out)
                             , pm0                   & ! intent(out)
                             , c1_m97                & ! intent(out)
                             , c2_m97                & ! intent(out)
                             , c3_m97                & ! intent(out)
                             , kvwake                & ! intent(out)
                             , alpha1_m97            & ! intent(out)
                             , alpha2_m97            & ! intent(out)
                             , psi_m97               & ! intent(out)
                             , zztop                 & ! intent(out)
                             , zzmid                 & ! intent(out)
                             , lad                   & ! intent(out)
                             , dladdz                & ! intent(out)
                             , cdrag                 & ! intent(out)
                             , pshelter              & ! intent(out)
                             , cumldrag              & ! intent(out)
                             , windm97               & ! intent(out)
                             , dz_m978               & ! intent(out)
                             , cdrag08               & ! intent(out)
                             , pm08                  & ! intent(out)
                             , c1_m978               & ! intent(out)
                             , c2_m978               & ! intent(out)
                             , c3_m978               & ! intent(out)
                             , kvwake8               & ! intent(out)
                             , alpha1_m978           & ! intent(out)
                             , alpha2_m978           & ! intent(out)
                             , psi_m978              & ! intent(out)
                             , zztop8                & ! intent(out)
                             , zzmid8                & ! intent(out)
                             , lad8                  & ! intent(out)
                             , dladdz8               & ! intent(out)
                             , cdrag8                & ! intent(out)
                             , pshelter8             & ! intent(out)
                             , cumldrag8             & ! intent(out)
                             , windm978              & ! intent(out)
                             , bl79                  & ! intent(out)
                             , csm                   & ! intent(out)
                             , csh                   & ! intent(out)
                             , dl79                  & ! intent(out)
                             , bbeta                 & ! intent(out)
                             , ribmaxod95            & ! intent(out)
                             , abh91                 & ! intent(out)
                             , bbh91                 & ! intent(out)
                             , cbh91                 & ! intent(out)
                             , dbh91                 & ! intent(out)
                             , ebh91                 & ! intent(out)
                             , fbh91                 & ! intent(out)
                             , cod                   & ! intent(out)
                             , bcod                  & ! intent(out)
                             , fm1                   & ! intent(out)
                             , ate                   & ! intent(out)
                             , atetf                 & ! intent(out)
                             , z0moz0h               & ! intent(out)
                             , z0hoz0m               & ! intent(out)
                             , ribmaxbh91            & ! intent(out)
                             , bl798                 & ! intent(out)
                             , csm8                  & ! intent(out)
                             , csh8                  & ! intent(out)
                             , dl798                 & ! intent(out)
                             , bbeta8                & ! intent(out)
                             , gamm8                 & ! intent(out)
                             , gamh8                 & ! intent(out)
                             , ribmaxod958           & ! intent(out)
                             , ribmaxbh918           & ! intent(out)
                             , tprandtl8             & ! intent(out)
                             , vkopr8                & ! intent(out)
                             , abh918                & ! intent(out)
                             , bbh918                & ! intent(out)
                             , cbh918                & ! intent(out)
                             , dbh918                & ! intent(out)
                             , ebh918                & ! intent(out)
                             , fbh918                & ! intent(out)
                             , cod8                  & ! intent(out)
                             , bcod8                 & ! intent(out)
                             , fm18                  & ! intent(out)
                             , ate8                  & ! intent(out)
                             , atetf8                & ! intent(out)
                             , z0moz0h8              & ! intent(out)
                             , z0hoz0m8              & ! intent(out)
                             , aflat_turb            & ! intent(out)
                             , aflat_lami            & ! intent(out)
                             , bflat_turb            & ! intent(out)
                             , bflat_lami            & ! intent(out)
                             , nflat_turb            & ! intent(out)
                             , nflat_lami            & ! intent(out)
                             , mflat_turb            & ! intent(out)
                             , mflat_lami            & ! intent(out)
                             , beta_r1               & ! intent(out)
                             , beta_r2               & ! intent(out)
                             , beta_re0              & ! intent(out)
                             , beta_g1               & ! intent(out)
                             , beta_g2               & ! intent(out)
                             , beta_gr0              & ! intent(out)
                             , aflat_turb8           & ! intent(out)
                             , aflat_lami8           & ! intent(out)
                             , bflat_turb8           & ! intent(out)
                             , bflat_lami8           & ! intent(out)
                             , nflat_turb8           & ! intent(out)
                             , nflat_lami8           & ! intent(out)
                             , mflat_turb8           & ! intent(out)
                             , mflat_lami8           & ! intent(out)
                             , beta_r18              & ! intent(out)
                             , beta_r28              & ! intent(out)
                             , beta_re08             & ! intent(out)
                             , beta_g18              & ! intent(out)
                             , beta_g28              & ! intent(out)
                             , beta_gr08             ! ! intent(out)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer :: ican
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Maximum leaf water that plants can hold.  Should leaf water exceed this number,    !
   ! water will be no longer intercepted by the leaves, and any value in excess of this    !
   ! will be promptly removed through shedding or throughfall.  This value is in           !
   ! kg/[m2 tree], so it will be scaled by (LAI+WAI) where needed be.                      !
   !---------------------------------------------------------------------------------------!
   fullveg_lwater = 0.11
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 tree], so it will be scaled by (LAI+WAI)    !
   ! where needed be.                                                                      !
   !---------------------------------------------------------------------------------------!
   dry_veg_lwater = 5.e-4 * fullveg_lwater
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Variables to define the vegetation aerodynamic resistance.  They are currently   !
   ! not PFT dependent.                                                                    !
   !---------------------------------------------------------------------------------------!
   select case (i_blyr_condct)
   case (-1)
      rb_slope =  25.0
      rb_inter =   0.0
   case default
      rb_slope =   0.0
      rb_inter =   1.e9
   end select

   !---------------------------------------------------------------------------------------!
   ! veg_height_min       - This is the minimum vegetation height allowed [m].  Vegetation !
   !                        height is used to calculate drag coefficients and patch        !
   !                        roughness.                                                     !
   ! minimum_canopy_depth - This is the minimum canopy depth allowed [m].  Canopy depth    !
   !                        is used to calculate the heat and moisture storage capacity in !
   !                        the canopy air space.                                          !
   !---------------------------------------------------------------------------------------!
   veg_height_min        = minval(hgt_min) ! alternative: minval(hgt_min) 
   minimum_canopy_depth  = 5.0             ! alternative: minval(hgt_min) 

   !----- This is the dimensionless exponential wind atenuation factor. -------------------!
   exar  = 2.5

   !----- This is the scaling factor of tree area index (not sure if it is used...) -------!
   covr = 2.16


   !---------------------------------------------------------------------------------------!
   !      Parameters for surface layer models.                                             !
   !---------------------------------------------------------------------------------------!
   !----- This is the minimum wind speed for boundary layer conductivity. -----------------!
   ugbmin    = 0.10
   !----- This is the minimum wind scale under stable and unstable conditions. ------------!
   ubmin     = 0.65
   !---------------------------------------------------------------------------------------!

   !----- Louis (1979) model. -------------------------------------------------------------!
   bl79        = 5.0    ! b prime parameter
   csm         = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   csh         = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   dl79        = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   bbeta       = 5.0           ! Beta 
   ribmaxod95  = 0.20          ! Maximum bulk Richardson number
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   abh91       = -1.00         ! -a from equation  (28) and (32)
   bbh91       = -twothirds    ! -b from equation  (28) and (32)
   cbh91       =  5.0          !  c from equations (28) and (32)
   dbh91       =  0.35         !  d from equations (28) and (32)
   ebh91       = -twothirds    ! - factor multiplying a*zeta in equation (32)
   fbh91       =  1.50         ! exponent in equation (32)
   cod         = cbh91/dbh91   ! c/d
   bcod        = bbh91 * cod   ! b*c/d
   fm1         = fbh91 - 1.0   ! f-1
   ate         = abh91 * ebh91 ! a * e
   atetf       = ate   * fbh91 ! a * e * f
   z0moz0h     = 1.0           ! z0(M)/z0(h)
   z0hoz0m     = 1. / z0moz0h  ! z0(M)/z0(h)
   ribmaxbh91  = 0.20          ! Maximum bulk Richardson number
   !---------------------------------------------------------------------------------------!


   
   !----- Legacy variable, we can probably remove it. -------------------------------------!
   ez  = 0.172
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !      Parameters for the aerodynamic resistance between the leaf and the canopy air    !
   ! space.  These are the A, B, n, and m parameters that define the Nusselt number for    !
   ! forced and free convection, at equations 10.7 and 10.9.  The parameters are found at  !
   ! the appendix A.5(a) and A.5(b).                                                       !
   !                                                                                       !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !---------------------------------------------------------------------------------------!
   aflat_turb = 0.600    ! A (forced convection), turbulent flow
   aflat_lami = 0.032    ! A (forced convection), laminar   flow
   bflat_turb = 0.500    ! B (free   convection), turbulent flow
   bflat_lami = 0.130    ! B (free   convection), laminar   flow
   nflat_turb = 0.500    ! n (forced convection), turbulent flow
   nflat_lami = 0.800    ! n (forced convection), laminar   flow
   mflat_turb = 0.250    ! m (free   convection), turbulent flow
   mflat_lami = onethird ! m (free   convection), laminar   flow
   !---------------------------------------------------------------------------------------!
   !     Both free and forced convection tend to underestimate the Nusselt number under    !
   ! different conditions.  Based on M08 review on the subject, I wrote the following      !
   ! functional form to expand the Nusselt number by a factor beta:                        !
   ! - beta_forced = R1 + R2 * tanh[log(Re/Re0)]                                           !
   ! - beta_free   = G1 + G2 * tanh[log(Gr/Gr0)]                                           !
   !     The values of beta change depending on the boundary layer conductance method.     !
   ! Currently only the forced convection varies, as we believe that the Reynolds number   !
   ! is the one that influences the most.                                                  ! 
   !---------------------------------------------------------------------------------------!
   select case (i_blyr_condct)
   case (-1,0)
      beta_r1  = 1.
      beta_r2  = 0.
      beta_re0 = 2000.
      beta_g1  = 1.
      beta_g2  = 0.
      beta_gr0 = 100000.
   case (1)
      beta_r1  =   7./4.
      beta_r2  =   3./4.
      beta_re0 =   2000.
      beta_g1  =   3./2.
      beta_g2  =  -1./2.
      beta_gr0 = 100000.
   case (2)
      beta_r1  =   11./2.
      beta_r2  =   9. /2.
      beta_re0 =   2000.
      beta_g1  =   3./2.
      beta_g2  =  -1./2.
      beta_gr0 = 100000.
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Define the variables that are going to be used by Massman (1997).                 !
   !---------------------------------------------------------------------------------------!
   !----- Discrete step size in canopy elevation [m]. -------------------------------------!
   dz_m97     = minval(hgt_min)
   !----- Number of canopy layers. --------------------------------------------------------!
   ncanmax    = ceiling(maxval(hgt_max)/dz_m97)
   !----- Fluid drag coefficient for turbulent flow in leaves. ----------------------------!
   cdrag0    = 0.2
   !----- Sheltering factor of fluid drag on canopies. ------------------------------------!
   pm0       = 1.0
   !----- Surface drag parameters (Massman 1997). -----------------------------------------!
   c1_m97    = 0.320 
   c2_m97    = 0.264
   c3_m97    = 15.1
   !----- Eddy diffusivity due to Von Karman Wakes in gravity flows. ----------------------!
   kvwake    = 0.0 ! 0.001
   !---------------------------------------------------------------------------------------!
   !     Alpha factors to produce the profile of sheltering factor and within canopy drag, !
   ! as suggested by Massman.                                                              !
   !---------------------------------------------------------------------------------------!
   alpha1_m97 = 0.40
   alpha2_m97 = 0.00
   !---------------------------------------------------------------------------------------!
   !      Parameter to represent the roughness sublayer effect.  According to Massman,     !
   ! assuming this to be zero means that the sublayer effects will be ignored.  Otherwise  !
   ! Raupach (1994) tried values up to 0.316.                                              !
   !---------------------------------------------------------------------------------------!
   psi_m97    = 0.316
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Allocate the scratch arrays.                                                      !
   !---------------------------------------------------------------------------------------!
   allocate(zztop    (0:ncanmax))
   allocate(zzmid    (  ncanmax))
   allocate(lad      (  ncanmax))
   allocate(dladdz   (  ncanmax))
   allocate(cdrag    (  ncanmax))
   allocate(pshelter (  ncanmax))
   allocate(cumldrag (  ncanmax))
   allocate(windm97  (  ncanmax))
   allocate(zztop8   (0:ncanmax))
   allocate(zzmid8   (  ncanmax))
   allocate(lad8     (  ncanmax))
   allocate(dladdz8  (  ncanmax))
   allocate(cdrag8   (  ncanmax))
   allocate(pshelter8(  ncanmax))
   allocate(cumldrag8(  ncanmax))
   allocate(windm978 (  ncanmax))
   !---------------------------------------------------------------------------------------!



   !----- Set the double precision variables. ---------------------------------------------!
   minimum_canopy_depth8 = dble(minimum_canopy_depth)
   exar8                 = dble(exar                )
   ubmin8                = dble(ubmin               )
   ugbmin8               = dble(ugbmin              )
   ustmin8               = dble(ustmin              )
   ggfact8               = dble(ggfact              )
   ez8                   = dble(ez                  )
   vh2dh8                = dble(vh2dh               )
   bl798                 = dble(bl79                )
   csm8                  = dble(csm                 )
   csh8                  = dble(csh                 )
   dl798                 = dble(dl79                )
   bbeta8                = dble(bbeta               )
   gamm8                 = dble(gamm                )
   gamh8                 = dble(gamh                )
   ribmaxod958           = dble(ribmaxod95          )
   ribmaxbh918           = dble(ribmaxbh91          )
   tprandtl8             = dble(tprandtl            )
   vkopr8                = dble(vkopr               )
   abh918                = dble(abh91               )
   bbh918                = dble(bbh91               )
   cbh918                = dble(cbh91               )
   dbh918                = dble(dbh91               )
   ebh918                = dble(ebh91               )
   fbh918                = dble(fbh91               )
   cod8                  = dble(cod                 )
   bcod8                 = dble(bcod                )
   fm18                  = dble(fm1                 )
   ate8                  = dble(ate                 )
   atetf8                = dble(atetf               )
   z0moz0h8              = dble(z0moz0h             )
   z0hoz0m8              = dble(z0hoz0m             )
   aflat_turb8           = dble(aflat_turb          )
   aflat_lami8           = dble(aflat_lami          )
   bflat_turb8           = dble(bflat_turb          )
   bflat_lami8           = dble(bflat_lami          )
   nflat_turb8           = dble(nflat_turb          )
   nflat_lami8           = dble(nflat_lami          )
   mflat_turb8           = dble(mflat_turb          )
   mflat_lami8           = dble(mflat_lami          )
   beta_r18              = dble(beta_r1             )
   beta_r28              = dble(beta_r2             )
   beta_re08             = dble(beta_re0            )
   beta_g18              = dble(beta_g1             )
   beta_g28              = dble(beta_g2             )
   beta_gr08             = dble(beta_gr0            )
   dz_m978               = dble(dz_m97              )
   cdrag08               = dble(cdrag0              )
   pm08                  = dble(pm0                 )
   c1_m978               = dble(c1_m97              )
   c2_m978               = dble(c2_m97              )
   c3_m978               = dble(c3_m97              )
   kvwake8               = dble(kvwake              )
   alpha1_m978           = dble(alpha1_m97          )
   alpha2_m978           = dble(alpha2_m97          )
   psi_m978              = dble(psi_m97             )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the maximum height of each layer.                                            !
   !---------------------------------------------------------------------------------------!
   do ican = 0,ncanmax
      zztop(ican)  = real(ican) * dz_m97
      zztop8(ican) = dble(ican) * dz_m978
   end do
   do ican = 1,ncanmax
      zzmid (ican) = 0.5   * (zztop (ican-1) + zztop (ican))
      zzmid8(ican) = 5.d-1 * (zztop8(ican-1) + zztop8(ican))
   end do
   !---------------------------------------------------------------------------------------!



   return
end subroutine init_can_air_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_photo_params()

   use ed_max_dims    , only : n_pft                ! ! intent(in)
   use pft_coms       , only : D0                   & ! intent(out)
                             , Vm_low_temp          & ! intent(out)
                             , Vm_high_temp         & ! intent(out)
                             , Vm0                  & ! intent(out)
                             , stomatal_slope       & ! intent(out)
                             , leaf_width           & ! intent(out)
                             , cuticular_cond       & ! intent(out)
                             , quantum_efficiency   & ! intent(out)
                             , photosyn_pathway     & ! intent(out)
                             , water_conductance    ! ! intent(out)
   use consts_coms    , only : t00                  & ! intent(in)
                             , twothirds            & ! intent(in)
                             , umol_2_mol           & ! intent(in)
                             , yr_sec               ! ! intent(in)
   use physiology_coms , only: vmfact               & ! intent(in)
                             , mfact                & ! intent(in)
                             , kfact                & ! intent(in)
                             , lwfact               & ! intent(in)
                             , thioff               ! ! intent(in)
   implicit none
   !---------------------------------------------------------------------------------------!

   D0(1:17)                  = 0.01            ! same for all PFTs

   Vm_low_temp(1)            = 5.0             ! c4 grass
   Vm_low_temp(2)            = 5.0             ! early tropical
   Vm_low_temp(3)            = 5.0             ! mid tropical
   Vm_low_temp(4)            = 5.0             ! late tropical
   Vm_low_temp(5)            = 4.7137          ! c3 grass
   Vm_low_temp(6)            = 4.7137          ! northern pines ! 5.0
   Vm_low_temp(7)            = 4.7137          ! southern pines ! 5.0
   Vm_low_temp(8)            = 4.7137          ! late conifers  ! 5.0
   Vm_low_temp(9)            = 4.7137          ! early hardwoods
   Vm_low_temp(10)           = 4.7137          ! mid hardwoods
   Vm_low_temp(11)           = 4.7137          ! late hardwoods
   Vm_low_temp(12)           = 4.7137          ! c3 pasture
   Vm_low_temp(13)           = 4.7137          ! c3 crop
   Vm_low_temp(14)           = 5.0             ! c4 pasture
   Vm_low_temp(15)           = 5.0             ! c4 crop
   Vm_low_temp(16)           = 5.0             ! subtropical C3 grass
   Vm_low_temp(17)           = 5.0             ! Araucaria

   Vm_high_temp(1)           = 100.0  + thioff ! C4
   Vm_high_temp(2)           =  45.0  + thioff ! C3
   Vm_high_temp(3)           =  45.0  + thioff ! C3
   Vm_high_temp(4)           =  45.0  + thioff ! C3
   Vm_high_temp(5)           =  45.0  + thioff ! C3
   Vm_high_temp(6)           =  45.0  + thioff ! C3
   Vm_high_temp(7)           =  45.0  + thioff ! C3
   Vm_high_temp(8)           =  45.0  + thioff ! C3
   Vm_high_temp(9)           =  45.0  + thioff ! C3
   Vm_high_temp(10)          =  45.0  + thioff ! C3
   Vm_high_temp(11)          =  45.0  + thioff ! C3
   Vm_high_temp(12)          =  45.0  + thioff ! C3
   Vm_high_temp(13)          =  45.0  + thioff ! C3
   Vm_high_temp(14)          = 100.0  + thioff ! C4
   Vm_high_temp(15)          = 100.0  + thioff ! C4
   Vm_high_temp(16)          =  45.0  + thioff ! C3
   Vm_high_temp(17)          =  45.0  + thioff ! C3

   !------ Vm0 is the maximum photosynthesis capacity in µmol/m2/s. -----------------------!
   Vm0(1)                    = 12.5            * 1.5
   Vm0(2)                    = 18.8            * vmfact
   Vm0(3)                    = 12.5            * vmfact
   Vm0(4)                    = 6.25            * vmfact
   Vm0(5)                    = 18.3            * vmfact
   Vm0(6)                    = 15.625 * 0.7264 * vmfact
   Vm0(7)                    = 15.625 * 0.7264 * vmfact
   Vm0(8)                    = 6.25   * 0.7264 * vmfact
   Vm0(9)                    = 18.25  * 1.1171 * vmfact
   Vm0(10)                   = 15.625 * 1.1171 * vmfact
   Vm0(11)                   = 6.25   * 1.1171 * vmfact
   Vm0(12:13)                = 18.3            * vmfact
   Vm0(14:15)                = 12.5            * 1.5
   Vm0(16)                   = 21.875          * vmfact
   Vm0(17)                   = 15.625 * 0.7264

   !----- Define the stomatal slope (aka the M factor). -----------------------------------!
   stomatal_slope(1)         =  6.4
   stomatal_slope(2)         =  8.0    * mfact
   stomatal_slope(3)         =  8.0    * mfact
   stomatal_slope(4)         =  8.0    * mfact
   stomatal_slope(5)         =  8.0    * mfact
   stomatal_slope(6)         =  6.3949 * mfact
   stomatal_slope(7)         =  6.3949 * mfact
   stomatal_slope(8)         =  6.3949 * mfact
   stomatal_slope(9)         =  6.3949 * mfact
   stomatal_slope(10)        =  6.3949 * mfact
   stomatal_slope(11)        =  6.3949 * mfact
   stomatal_slope(12)        =  8.0    * mfact
   stomatal_slope(13)        =  8.0    * mfact
   stomatal_slope(14)        =  6.4
   stomatal_slope(15)        =  6.4
   stomatal_slope(16)        =  8.0    * mfact
   stomatal_slope(17)        =  6.4
 
   cuticular_cond(1)         = 10000.0    ! 10000.0
   cuticular_cond(2)         = 10000.0    ! 10000.0
   cuticular_cond(3)         = 10000.0    ! 10000.0
   cuticular_cond(4)         = 10000.0    ! 10000.0
   cuticular_cond(5)         = 10000.0    ! 10000.0
   cuticular_cond(6)         = 1000.0 
   cuticular_cond(7)         = 1000.0 
   cuticular_cond(8)         = 1000.0 
   cuticular_cond(9)         = 20000.0
   cuticular_cond(10)        = 20000.0
   cuticular_cond(11)        = 20000.0
   cuticular_cond(12)        = 20000.0    ! 10000.0
   cuticular_cond(13)        = 20000.0    ! 10000.0
   cuticular_cond(14)        = 10000.0    ! 10000.0
   cuticular_cond(15)        = 10000.0    ! 10000.0
   cuticular_cond(16)        = 10000.0
   cuticular_cond(17)        = 1000.0 

   quantum_efficiency(1)     = 0.053
   quantum_efficiency(2)     = 0.08
   quantum_efficiency(3)     = 0.08
   quantum_efficiency(4)     = 0.08
   quantum_efficiency(5)     = 0.08
   quantum_efficiency(6)     = 0.08
   quantum_efficiency(7)     = 0.08
   quantum_efficiency(8)     = 0.08
   quantum_efficiency(9)     = 0.08
   quantum_efficiency(10)    = 0.08
   quantum_efficiency(11)    = 0.08
   quantum_efficiency(12)    = 0.08
   quantum_efficiency(13)    = 0.08
   quantum_efficiency(14)    = 0.053
   quantum_efficiency(15)    = 0.053
   quantum_efficiency(16)    = 0.08
   quantum_efficiency(17)    = 0.08

   !---------------------------------------------------------------------------------------!
   !     The KW parameter. Medvigy et al. (2009) and Moorcroft et al. (2001) give the      !
   ! number in m²/yr/kg_C_root.  Here we must define it in m²/s/kg_C_root.                 !
   !---------------------------------------------------------------------------------------!
   water_conductance(1:17) = 150. / yr_sec * kfact
   !---------------------------------------------------------------------------------------!


   photosyn_pathway(1)       = 4
   photosyn_pathway(2:4)     = 3
   photosyn_pathway(5)       = 3
   photosyn_pathway(6:13)    = 3
   photosyn_pathway(14:15)   = 4
   photosyn_pathway(16:17)   = 3

   !----- Leaf width [m].  This controls the boundary layer conductance. ------------------!
   leaf_width( 1)    = 0.20 * lwfact
   leaf_width( 2)    = 0.20 * lwfact
   leaf_width( 3)    = 0.20 * lwfact
   leaf_width( 4)    = 0.20 * lwfact
   leaf_width( 5)    = 0.05 * lwfact
   leaf_width( 6)    = 0.05 * lwfact
   leaf_width( 7)    = 0.05 * lwfact
   leaf_width( 8)    = 0.05 * lwfact
   leaf_width( 9)    = 0.05 * lwfact
   leaf_width(10)    = 0.05 * lwfact
   leaf_width(11)    = 0.05 * lwfact
   leaf_width(12)    = 0.05 * lwfact
   leaf_width(13)    = 0.05 * lwfact
   leaf_width(14)    = 0.20 * lwfact
   leaf_width(15)    = 0.20 * lwfact
   leaf_width(16)    = 0.20 * lwfact
   leaf_width(17)    = 0.05 * lwfact
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_photo_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_decomp_params()
   use consts_coms , only : yr_day                      ! ! intent(in)
   use decomp_coms , only : resp_opt_water              & ! intent(in)
                          , resp_water_below_opt        & ! intent(in)
                          , resp_water_above_opt        & ! intent(in)
                          , resp_temperature_increase   & ! intent(in)
                          , N_immobil_supply_scale      & ! intent(in)
                          , cwd_frac                    & ! intent(in)
                          , r_fsc                       & ! intent(in)
                          , r_stsc                      & ! intent(in)
                          , r_ssc                       & ! intent(in)
                          , K1                          & ! intent(in)
                          , K2                          & ! intent(in)
                          , K3                          ! ! intent(in)

   resp_opt_water            = 0.8938
   resp_water_below_opt      = 5.0786
   resp_water_above_opt      = 4.5139
   resp_temperature_increase = 0.0757
   N_immobil_supply_scale    = 40.0 / yr_day
   cwd_frac                  = 0.2
   r_fsc                     = 1.0
   r_stsc                    = 0.3
   r_ssc                     = 1.0
   K1                        = 4.5   / yr_day
   K2                        = 11.0  / yr_day
   K3                        = 100.2 / yr_day

   return

end subroutine init_decomp_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_resp_params()

   use physiology_coms, only : gamfact                   ! ! intent(in)
   use pft_coms       , only : growth_resp_factor        & ! intent(out)
                             , leaf_turnover_rate        & ! intent(out)
                             , root_turnover_rate        & ! intent(out)
                             , dark_respiration_factor   & ! intent(out)
                             , storage_turnover_rate     & ! intent(out)
                             , root_respiration_factor   ! ! intent(out)
   use decomp_coms    , only : f_labile                  ! ! intent(out)
   use consts_coms    , only : onesixth                  & ! intent(in)
                             , onethird                  ! ! intent(in)
   implicit none

   growth_resp_factor(1)          = onethird
   growth_resp_factor(2)          = onethird
   growth_resp_factor(3)          = onethird
   growth_resp_factor(4)          = onethird
   growth_resp_factor(5)          = onethird
   growth_resp_factor(6)          = 0.4503 ! 0.333
   growth_resp_factor(7)          = 0.4503
   growth_resp_factor(8)          = 0.4503 ! 0.333
   growth_resp_factor(9)          = 0.0
   growth_resp_factor(10)         = 0.0
   growth_resp_factor(11)         = 0.0
   growth_resp_factor(12)         = onethird
   growth_resp_factor(13)         = onethird
   growth_resp_factor(14)         = onethird
   growth_resp_factor(15)         = onethird
   growth_resp_factor(16)         = onethird
   growth_resp_factor(17)         = 0.4503

   leaf_turnover_rate(1)          = 2.0
   leaf_turnover_rate(2)          = 1.0
   leaf_turnover_rate(3)          = 0.5
   leaf_turnover_rate(4)          = onethird
   leaf_turnover_rate(5)          = 2.0
   leaf_turnover_rate(6)          = onethird
   leaf_turnover_rate(7)          = onethird
   leaf_turnover_rate(8)          = onethird
   leaf_turnover_rate(9)          = 0.0
   leaf_turnover_rate(10)         = 0.0
   leaf_turnover_rate(11)         = 0.0
   leaf_turnover_rate(12)         = 2.0
   leaf_turnover_rate(13)         = 2.0
   leaf_turnover_rate(14)         = 2.0
   leaf_turnover_rate(15)         = 2.0
   leaf_turnover_rate(16)         = 2.0
   leaf_turnover_rate(17)         = onesixth

   !----- Root turnover rate.  ------------------------------------------------------------!
   root_turnover_rate(1)          = 2.0
   root_turnover_rate(2)          = 1.0
   root_turnover_rate(3)          = 0.5
   root_turnover_rate(4)          = onethird
   root_turnover_rate(5)          = 2.0
   root_turnover_rate(6)          = 3.927218 ! 0.333
   root_turnover_rate(7)          = 4.117847 ! 0.333
   root_turnover_rate(8)          = 3.800132 ! 0.333
   root_turnover_rate(9)          = 5.772506
   root_turnover_rate(10)         = 5.083700
   root_turnover_rate(11)         = 5.070992
   root_turnover_rate(12)         = onethird
   root_turnover_rate(13)         = onethird
   root_turnover_rate(14)         = 2.0
   root_turnover_rate(15)         = 2.0
   root_turnover_rate(16)         = 2.0
   root_turnover_rate(17)         = onesixth

   dark_respiration_factor(1)     = 0.06
   dark_respiration_factor(2)     = 0.02  * gamfact
   dark_respiration_factor(3)     = 0.02  * gamfact
   dark_respiration_factor(4)     = 0.02  * gamfact
   dark_respiration_factor(5)     = 0.02
   dark_respiration_factor(6)     = 0.02
   dark_respiration_factor(7)     = 0.02
   dark_respiration_factor(8)     = 0.02
   dark_respiration_factor(9)     = 0.02
   dark_respiration_factor(10)    = 0.02
   dark_respiration_factor(11)    = 0.02
   dark_respiration_factor(12)    = 0.02
   dark_respiration_factor(13)    = 0.02
   dark_respiration_factor(14)    = 0.04
   dark_respiration_factor(15)    = 0.04
   dark_respiration_factor(16)    = 0.02  * gamfact
   dark_respiration_factor(17)    = 0.025 * gamfact

   storage_turnover_rate(1)       = 0.00 ! 0.25
   storage_turnover_rate(2)       = 0.00 ! 0.25
   storage_turnover_rate(3)       = 0.00 ! 0.25
   storage_turnover_rate(4)       = 0.00 ! 0.25
   storage_turnover_rate(5)       = 0.00 ! 0.25
   storage_turnover_rate(6)       = 0.00 ! 0.25
   storage_turnover_rate(7)       = 0.00 ! 0.25
   storage_turnover_rate(8)       = 0.00 ! 0.25
   storage_turnover_rate(9)       = 0.6243
   storage_turnover_rate(10)      = 0.6243
   storage_turnover_rate(11)      = 0.6243
   storage_turnover_rate(12)      = 0.00 ! 0.25
   storage_turnover_rate(13)      = 0.00 ! 0.25
   storage_turnover_rate(14)      = 0.00 ! 0.25
   storage_turnover_rate(15)      = 0.00 ! 0.25
   storage_turnover_rate(16)      = 0.00 ! 0.25
   storage_turnover_rate(17)      = 0.00 ! 0.25

   root_respiration_factor(1:17)  = 0.528

   f_labile(1:5)                  = 1.0
   f_labile(6:11)                 = 0.79
   f_labile(12:15)                = 1.0
   f_labile(16)                   = 1.0
   f_labile(17)                   = 0.79

   return
end subroutine init_pft_resp_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns some PFT-dependent parameters that control mortality rates.  !
!------------------------------------------------------------------------------------------!
subroutine init_pft_mort_params()

   use pft_coms    , only : mort1                      & ! intent(out)
                          , mort2                      & ! intent(out)
                          , mort3                      & ! intent(out)
                          , rho                        & ! intent(out)
                          , seedling_mortality         & ! intent(out)
                          , treefall_s_gtht            & ! intent(out)
                          , treefall_s_ltht            & ! intent(out)
                          , plant_min_temp             & ! intent(out)
                          , frost_mort                 ! ! intent(out)
   use consts_coms , only : t00                        & ! intent(in)
                          , lnexp_max                  ! ! intent(in)
   use disturb_coms, only : treefall_disturbance_rate  & ! intent(inout)
                          , time2canopy                ! ! intent(in)

   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   real     :: aquad
   real     :: bquad
   real     :: cquad
   real     :: discr
   real     :: lambda_ref
   real     :: lambda_eff
   real     :: leff_neg
   real     :: leff_pos
   !---------------------------------------------------------------------------------------!


   frost_mort(1)     = 3.0
   frost_mort(2:4)   = 3.0
   frost_mort(5)     = 3.0
   frost_mort(6:11)  = 3.0
   frost_mort(12:13) = 3.0
   frost_mort(14:15) = 3.0
   frost_mort(16:17) = 3.0


   mort1(1)  = 10.0
   mort1(2)  = 10.0
   mort1(3)  = 10.0
   mort1(4)  = 10.0
   mort1(5)  = 1.0
   mort1(6)  = 1.0
   mort1(7)  = 1.0
   mort1(8)  = 1.0
   mort1(9)  = 1.0
   mort1(10) = 1.0
   mort1(11) = 1.0
   mort1(12) = 1.0
   mort1(13) = 1.0
   mort1(14) = 10.0
   mort1(15) = 10.0
   mort1(16) = 10.0
   mort1(17) = 10.0

   mort2(1)  = 20.0
   mort2(2)  = 20.0
   mort2(3)  = 20.0
   mort2(4)  = 20.0
   mort2(5)  = 20.0
   mort2(6)  = 20.0
   mort2(7)  = 20.0
   mort2(8)  = 20.0
   mort2(9)  = 20.0
   mort2(10) = 20.0
   mort2(11) = 20.0
   mort2(12) = 20.0
   mort2(13) = 20.0
   mort2(14) = 20.0
   mort2(15) = 20.0
   mort2(16) = 20.0
   mort2(17) = 20.0

   mort3(1)  = 0.15 * (1. - rho(1) / rho(4)) 
   mort3(2)  = 0.15 * (1. - rho(2) / rho(4))  
   mort3(3)  = 0.15 * (1. - rho(3) / rho(4))
   mort3(4)  = 0.0
   mort3(5)  = 0.066
   mort3(6)  = 0.0033928
   mort3(7)  = 0.0043
   mort3(8)  = 0.0023568
   mort3(9)  = 0.006144
   mort3(10) = 0.003808
   mort3(11) = 0.00428
   mort3(12) = 0.066
   mort3(13) = 0.066
   mort3(14) = 0.037
   mort3(15) = 0.037
   mort3(16) = 0.06167
   mort3(17) = 0.0043



   !---------------------------------------------------------------------------------------!
   !     Here we check whether we need to re-calculate the treefall disturbance rate so it !
   ! is consistent with the time to reach the canopy.                                      !
   !---------------------------------------------------------------------------------------!
   if (treefall_disturbance_rate == 0.) then
      !------ No disturbance rate, set time to reach canopy to infinity. ------------------!
      time2canopy = huge(1.) 
      lambda_ref  = 0.
      lambda_eff  = 0.

   else
      lambda_ref = abs(treefall_disturbance_rate)

      if (time2canopy > 0.) then
         !---------------------------------------------------------------------------------!
         !     We are not going to knock down trees as soon as the patch is created;       !
         ! instead, we will wait until the patch age is older than time2canopy.  We want,  !
         ! however, to make the mean patch age to be 1/treefall_disturbance_rate.  The     !
         ! equation below can be retrieved by integrating the steady-state probability     !
         ! distribution function.  The equation is quadratic and the discriminant will     !
         ! never be zero and the treefall_disturbance_rate will be always positive because !
         ! the values of time2canopy and treefall_disturbance_rate have already been       !
         ! tested in ed_opspec.F90.                                                        !
         !---------------------------------------------------------------------------------!
         aquad    = time2canopy * time2canopy * lambda_ref  - 2. * time2canopy
         bquad    = 2. * time2canopy * lambda_ref - 2.
         cquad    = 2. * lambda_ref
         !------ Find the discriminant. ---------------------------------------------------!
         discr    = bquad * bquad - 4. * aquad * cquad
         leff_neg = - 0.5 * (bquad - sqrt(discr)) / aquad
         leff_pos = - 0.5 * (bquad + sqrt(discr)) / aquad
         !---------------------------------------------------------------------------------!
         !      Use the maximum value, but don't let the value to be too large otherwise   !
         ! the negative exponential will cause underflow.                                  !
         !---------------------------------------------------------------------------------!
         lambda_eff = min(lnexp_max,max(leff_neg,leff_pos))
      else
         lambda_eff = lambda_ref
      end if

   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Print out the summary.                                                             !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)')           '----------------------------------------'
   write (unit=*,fmt='(a)')           '  Treefall disturbance parameters:'
   write (unit=*,fmt='(a,1x,es12.5)') '  - LAMBDA_REF  =',lambda_ref
   write (unit=*,fmt='(a,1x,es12.5)') '  - LAMBDA_EFF  =',lambda_eff
   write (unit=*,fmt='(a,1x,es12.5)') '  - TIME2CANOPY =',time2canopy
   write (unit=*,fmt='(a)')           '----------------------------------------'
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we check whether patches should be created or the treefall should affect     !
   ! only the mortality (quasi- size-structured approximation; other disturbances may be   !
   ! turned off for a true size-structured approximation).                                 !
   !---------------------------------------------------------------------------------------!
   if (treefall_disturbance_rate < 0.) then
      !------------------------------------------------------------------------------------!
      !      We incorporate the disturbance rate into the density-independent mortality    !
      ! rate and turn off the patch-creating treefall disturbance.                         !
      !------------------------------------------------------------------------------------!
      mort3(:) = mort3(:) + lambda_eff
      treefall_disturbance_rate = 0.
   else
      treefall_disturbance_rate = lambda_eff
   end if
   !---------------------------------------------------------------------------------------!


   seedling_mortality(1)    = 0.60
   seedling_mortality(2:4)  = 0.95 
   seedling_mortality(5)    = 0.60
   seedling_mortality(6:15) = 0.95 
   seedling_mortality(16)   = 0.60 
   seedling_mortality(17)   = 0.95 

   treefall_s_gtht          = 0.0

   treefall_s_ltht(1)       = 0.25
   treefall_s_ltht(2:4)     = 0.1
   treefall_s_ltht(5)       = 0.25
   treefall_s_ltht(6:11)    = 0.1
   treefall_s_ltht(12:15)   = 0.25
   treefall_s_ltht(16)      = 0.25
   treefall_s_ltht(17)      = 0.1

   plant_min_temp(1:4)      = t00+2.5
   plant_min_temp(5:6)      = t00-80.0
   plant_min_temp(7)        = t00-10.0
   plant_min_temp(8)        = t00-60.0
   plant_min_temp(9)        = t00-80.0
   plant_min_temp(10:11)    = t00-20.0
   plant_min_temp(12:13)    = t00-80.0
   plant_min_temp(14:15)    = t00+2.5
   plant_min_temp(16)       = t00-5.0
   plant_min_temp(17)       = t00-10.0

   return
end subroutine init_pft_mort_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_alloc_params()

   use pft_coms    , only : leaf_turnover_rate    & ! intent(in)
                          , is_tropical           & ! intent(out)
                          , is_grass              & ! intent(out)
                          , rho                   & ! intent(out)
                          , SLA                   & ! intent(out)
                          , horiz_branch          & ! intent(out)
                          , q                     & ! intent(out)
                          , qsw                   & ! intent(out)
                          , init_density          & ! intent(out)
                          , agf_bs                & ! intent(out)
                          , agf_bsi               & ! intent(out)
                          , hgt_min               & ! intent(out)
                          , hgt_ref               & ! intent(out)
                          , hgt_max               & ! intent(out)
                          , min_dbh               & ! intent(out)
                          , max_dbh               & ! intent(out)
                          , b1Ht                  & ! intent(out)
                          , b2Ht                  & ! intent(out)
                          , b1Bs_small            & ! intent(out)
                          , b2Bs_small            & ! intent(out)
                          , b1Bs_big              & ! intent(out)
                          , b2Bs_big              & ! intent(out)
                          , B1Ca                  & ! intent(out)
                          , B2Ca                  & ! intent(out)
                          , bdead_crit            & ! intent(out)
                          , b1Bl                  & ! intent(out)
                          , b2Bl                  & ! intent(out)
                          , C2B                   & ! intent(out)
                          , sapwood_ratio         & ! intent(out)
                          , rbranch               & ! intent(out)
                          , rdiamet               & ! intent(out)
                          , rlength               & ! intent(out)
                          , diammin               & ! intent(out)
                          , ntrunk                & ! intent(out)
                          , conijn_a              & ! intent(out)
                          , conijn_b              & ! intent(out)
                          , conijn_c              & ! intent(out)
                          , conijn_d              ! ! intent(out)
   use allometry   , only : h2dbh                 & ! function
                          , dbh2bd                ! ! function
   use consts_coms , only : twothirds             ! ! intent(in)
   use ed_max_dims , only : n_pft                 & ! intent(in)
                          , str_len               ! ! intent(in)
   use ed_misc_coms, only : iallom                ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer                           :: ipft
   integer                           :: n
   real                              :: aux
   !----- Constants shared by both bdead and bleaf (tropical PFTs) ------------------------!
   real                  , parameter :: a1          =  -1.981
   real                  , parameter :: b1          =   1.047
   real                  , parameter :: dcrit       = 100.0
   !----- Constants used by bdead only (tropical PFTs) ------------------------------------!
   real                  , parameter :: c1d         =   0.572
   real                  , parameter :: d1d         =   0.931
   real                  , parameter :: a2d         =  -1.086
   real                  , parameter :: b2d         =   0.876
   real                  , parameter :: c2d         =   0.604
   real                  , parameter :: d2d         =   0.871
   !----- Constants used by bleaf only (tropical PFTs) ------------------------------------!
   real                  , parameter :: c1l         =  -0.584
   real                  , parameter :: d1l         =   0.550
   real                  , parameter :: a2l         =  -4.111
   real                  , parameter :: b2l         =   0.605
   real                  , parameter :: c2l         =   0.848
   real                  , parameter :: d2l         =   0.438
   !---------------------------------------------------------------------------------------!
   !     MLO.   These are the new parameters obtained by adjusting a curve that is similar !
   !            to the modified Chave's equation to include wood density effect on the     !
   !            DBH->AGB allometry as described by:                                        ! 
   !                                                                                       !
   !            Baker, T. R., and co-authors, 2004: Variation in wood density determines   !
   !               spatial patterns in Amazonian forest biomass.  Glob. Change Biol., 10,  !
   !               545-562.                                                                !
   !                                                                                       !
   !            These parameters were obtaining by splitting balive and bdead at the same  !
   !            ratio as the original ED-2.1 allometry, and optimising a function of the   !
   !            form B? = (rho / a3) * exp [a1 + a2 * ln(DBH)]                             !
   !---------------------------------------------------------------------------------------!
   real, dimension(3)    , parameter :: aleaf       = (/ -1.259299,  1.679213,  4.985562 /)
   real, dimension(3)    , parameter :: adead_small = (/ -1.494639,  2.453309,  1.597272 /)
   real, dimension(3)    , parameter :: adead_big   = (/  2.105856,  2.423031, 50.198984 /)
   !----- Other constants. ----------------------------------------------------------------!
   logical               , parameter :: write_allom = .true.
   character(len=str_len), parameter :: allom_file  = 'allom_param.txt'
   !---------------------------------------------------------------------------------------!

   !----- Carbon-to-biomass ratio of plant tissues. ---------------------------------------!
   C2B    = 2.0

   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is tropical/subtropical or    !
   ! not.                                                                                  !
   !---------------------------------------------------------------------------------------! 
   is_tropical(1:4)   = .true.
   is_tropical(5:11)  = .false.
   is_tropical(12:13) = .false.
   is_tropical(14:15) = .true.
   is_tropical(16)    = .true.
   !---------------------------------------------------------------------------------------!
   !     This uses tropical allometry for DBH->Bleaf and DBH->Bdead, but otherwise it uses !
   ! the temperate properties.                                                             !
   !---------------------------------------------------------------------------------------!
   is_tropical(17)    = .true.

   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is tree or grass              !
   !---------------------------------------------------------------------------------------! 
   is_grass(1)     = .true.
   is_grass(2:4)   = .false.
   is_grass(5)     = .true.
   is_grass(6:11)  = .false.
   is_grass(12:15) = .true.
   is_grass(16)    = .true.
   is_grass(17)    = .false.

   !---------------------------------------------------------------------------------------!
   !     Wood density.  Currently only tropical PFTs need it.  C3 grass density will be    !
   ! used only for branch area purposes.                                                   !
   !---------------------------------------------------------------------------------------!
   !---- [KIM] new tropical parameters. ---------------------------------------------------!
   rho(1)     = 0.32   ! 0.40
   rho(2)     = 0.53   ! 0.40
   rho(3)     = 0.71   ! 0.60
   rho(4)     = 0.90   ! 0.87
   rho(5)     = 0.32   ! Copied from C4 grass
   rho(6:11)  = 0.00   ! Currently not used
   rho(12:13) = 0.32
   rho(14:15) = 0.32
   rho(16)    = 0.32
   rho(17)    = 0.48
   !---------------------------------------------------------------------------------------!

   !----- Specific leaf area [m² leaf / kg C] ---------------------------------------------!
   !----- [KIM] - new tropical parameters. ------------------------------------------------!
   SLA(1:4)   = 10.0**(2.4-0.46*log10(12.0/leaf_turnover_rate(1:4))) * C2B * 0.1
   ! SLA(1:4) = 10.0**(1.6923-0.3305*log10(12.0/leaf_turnover_rate(1:4)))
   SLA(5)     = 22.0
   SLA(6)     =  6.0
   SLA(7)     =  9.0
   SLA(8)     = 10.0
   SLA(9)     = 30.0
   SLA(10)    = 24.2
   SLA(11)    = 60.0
   SLA(12:13) = 22.0
   SLA(14:15) = 10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(14:15)))) * C2B * 0.1
   ! SLA(14:15) = 10.0**(1.6923-0.3305*log10(12.0/leaf_turnover_rate(14:15)))
   SLA(16)    = 10.0**(2.4-0.46*log10(12.0/leaf_turnover_rate(16))) * C2B * 0.1
   SLA(17)    = 10.0

   !---------------------------------------------------------------------------------------!
   !    Fraction of vertical branches.  Values are from Poorter et al. (2006):             !
   !                                                                                       !
   !    Poorter, L.; Bongers, L.; Bongers, F., 2006: Architecture of 54 moist-forest tree  !
   ! species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.           !
   ! For simplicity, we assume similar numbers for temperate PFTs.                         !
   !---------------------------------------------------------------------------------------!
   horiz_branch(1)     = 0.50
   horiz_branch(2)     = 0.57
   horiz_branch(3)     = 0.39
   horiz_branch(4)     = 0.61
   horiz_branch(5)     = 0.50
   horiz_branch(6:8)   = 0.61
   horiz_branch(9)     = 0.57
   horiz_branch(10)    = 0.39
   horiz_branch(11)    = 0.61
   horiz_branch(12:15) = 0.50
   horiz_branch(16)    = 0.50
   horiz_branch(17)    = 0.61
   !---------------------------------------------------------------------------------------!


   !----- Ratio between fine roots and leaves [kg_fine_roots/kg_leaves] -------------------!
   q(1)     = 1.0
   q(2)     = 1.0
   q(3)     = 1.0
   q(4)     = 1.0
   q(5)     = 1.0
   q(6)     = 0.3463 ! 1.0
   q(7)     = 0.3463 ! 1.0
   q(8)     = 0.3463 ! 1.0
   q(9)     = 1.1274
   q(10)    = 1.1274
   q(11)    = 1.1274
   q(12:15) = 1.0
   q(16)    = 1.0
   q(17)    = 1.0

   sapwood_ratio(1:17) = 3900.0

   !---------------------------------------------------------------------------------------!
   !    Finding the ratio between sapwood and leaves [kg_sapwood/kg_leaves]                !
   !                                                                                       !
   !    KIM: ED1/ED2 codes and Moorcroft et al. had the incorrect ratio.  Since the mid-   !
   ! latitude parameters have been optimized using the wrong SLA, we keep the bug until    !
   ! it is updated...                                                                      !
   !---------------------------------------------------------------------------------------!
   qsw(1:4)    = SLA(1:4)   / sapwood_ratio(1:4)    !new is SLA(1:4)/(3900.0*2.0/1000.0)
   qsw(5:13)   = SLA(5:13)  / sapwood_ratio(5:13)
   qsw(14:15)  = SLA(14:15) / sapwood_ratio(14:15)  !new is SLA(14:15)(3900.0*2.0/1000.0)
   qsw(16)     = SLA(16)    / sapwood_ratio(16)
   qsw(17)     = SLA(17)    / sapwood_ratio(17)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Initial density of plants, for near-bare-ground simulations [# of individuals/m2]  !
   !---------------------------------------------------------------------------------------!
   init_density(1)     = 0.1
   init_density(2:4)   = 0.1
   init_density(5)     = 0.1
   init_density(6:8)   = 0.1
   init_density(9:11)  = 0.1
   init_density(12:13) = 0.1
   init_density(14:15) = 0.1
   init_density(16)    = 0.1
   init_density(17)    = 0.1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Minimum height of an individual.                                                   !
   !---------------------------------------------------------------------------------------!
   hgt_min(1)     = 0.50
   hgt_min(2:4)   = 0.50
   hgt_min(5)     = 0.15
   hgt_min(6)     = 1.50
   hgt_min(7)     = 1.50
   hgt_min(8)     = 1.50
   hgt_min(9)     = 1.50
   hgt_min(10)    = 1.50
   hgt_min(11)    = 1.50
   hgt_min(12)    = 0.15
   hgt_min(13)    = 0.15
   hgt_min(14)    = 0.50
   hgt_min(15)    = 0.50
   hgt_min(16)    = 0.50
   hgt_min(17)    = 0.50
   !---------------------------------------------------------------------------------------!



   !----- Reference height for diameter/height allometry (temperates only). ---------------!
   hgt_ref(1:5)   = 0.0
   hgt_ref(6:11)  = 1.3
   hgt_ref(12:15) = 0.0
   hgt_ref(16)    = 0.0
   hgt_ref(17)    = 0.0
   !---------------------------------------------------------------------------------------!



   !----- Fraction of structural stem that is assumed to be above ground. -----------------!
   agf_bs  = 0.7
   agf_bsi = 1. / agf_bs
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   DBH/height allometry parameters.                                                    !
   !                                                                                       !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !                                                                                       !
   !   These parameters have different meaning for tropical and temperate PFTs...          !
   !---------------------------------------------------------------------------------------!
   !----- DBH-height allometry intercept [m]. ---------------------------------------------!
   b1Ht(1:4)   = 0.37 * log(10.0)
   b1Ht(5)     = 0.4778
   b1Ht(6)     = 27.14
   b1Ht(7)     = 27.14
   b1Ht(8)     = 22.79
   b1Ht(9)     = 22.6799
   b1Ht(10)    = 25.18
   b1Ht(11)    = 23.3874
   b1Ht(12:13) = 0.4778
   b1Ht(14:15) = 0.37 * log(10.0)
   b1Ht(16)    = 0.37 * log(10.0)
   b1Ht(17)    = 0.37 * log(10.0)
   !----- DBH-height allometry slope [1/cm]. ----------------------------------------------!
   b2Ht(1:4)   = 0.64
   b2Ht(5)     = -0.75
   b2Ht(6)     = -0.03884
   b2Ht(7)     = -0.03884
   b2Ht(8)     = -0.04445 
   b2Ht(9)     = -0.06534
   b2Ht(10)    = -0.04964
   b2Ht(11)    = -0.05404
   b2Ht(12:13) = -0.75
   b2Ht(14:15) =  0.64
   b2Ht(16)    =  0.64
   b2Ht(17)    =  0.64

   !----- Maximum Height. -----------------------------------------------------------------!
   hgt_max( 1) = 1.50
   hgt_max( 2) = 35.0
   hgt_max( 3) = 35.0
   hgt_max( 4) = 35.0
   hgt_max( 5) = 0.95  * b1Ht( 5)
   hgt_max( 6) = 0.999 * b1Ht( 6)
   hgt_max( 7) = 0.999 * b1Ht( 7)
   hgt_max( 8) = 0.999 * b1Ht( 8)
   hgt_max( 9) = 0.999 * b1Ht( 9)
   hgt_max(10) = 0.999 * b1Ht(10)
   hgt_max(11) = 0.999 * b1Ht(11)
   hgt_max(12) = 0.95  * b1Ht(12)
   hgt_max(13) = 0.95  * b1Ht(13)
   hgt_max(14) = 1.50
   hgt_max(15) = 1.50
   hgt_max(16) = 1.50
   hgt_max(17) = 35.0

   !----- Maximum DBH. --------------------------------------------------------------------!
   do ipft=1,n_pft
      min_dbh(ipft) = h2dbh(hgt_min(ipft),ipft)
      max_dbh(ipft) = h2dbh(hgt_max(ipft),ipft)
   end do


   !---------------------------------------------------------------------------------------!
   !     DBH-leaf allometry.  Assign temperate PFTs outside the loop, and the tropical     !
   ! ones inside the loop.                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- DBH-leaf allometry intercept [kg leaf biomass / plant * cm^(-b2Bl)]. ------------!
   b1Bl(1:4)   = 0.0
   b1Bl(5)     = 0.08
   b1Bl(6)     = 0.024
   b1Bl(7)     = 0.024
   b1Bl(8)     = 0.0454
   b1Bl(9)     = 0.0129
   b1Bl(10)    = 0.048
   b1Bl(11)    = 0.017
   b1Bl(12:13) = 0.08
   b1Bl(14:15) = 0.0
   b1Bl(16)    = 0.0
   b1Bl(17)    = 0.0
   !-----  DBH-leaf allometry slope [dimensionless]. --------------------------------------!
   b2Bl(1:4)   = 0.0
   b2Bl(5)     = 1.0
   b2Bl(6)     = 1.899
   b2Bl(7)     = 1.899
   b2Bl(8)     = 1.6829
   b2Bl(9)     = 1.7477
   b2Bl(10)    = 1.455
   b2Bl(11)    = 1.731
   b2Bl(12:13) = 1.0
   b2Bl(14:15) = 0.0
   b2Bl(16)    = 0.0
   b2Bl(17)    = 0.0
   !------- Fill in the tropical PFTs, which are functions of wood density. ---------------!
   do ipft=1,n_pft
      if (is_tropical(ipft)) then
         select case(iallom)
         case (0)
            !---- ED-2.1 allometry. -------------------------------------------------------!
            b1Bl(ipft) = exp(a1 + c1l * b1Ht(ipft) + d1l * log(rho(ipft)))
            aux        = ( (a2l - a1) + b1Ht(ipft) * (c2l - c1l) + log(rho(ipft))          &
                         * (d2l - d1l)) * (1.0/log(dcrit))
            b2Bl(ipft) = C2B * b2l + c2l * b2Ht(ipft) + aux
         case (1)
            !---- Based on modified Chave et al. (2001) allometry. ------------------------!
            b1Bl(ipft) = C2B * exp(aleaf(1)) * rho(ipft) / aleaf(3)
            b2Bl(ipft) = aleaf(2)
         end select
      end if
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     DBH-stem allometry.  Assign temperate PFTs outside the loop, and the tropical     !
   ! ones inside the loop.                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- DBH-stem allometry intercept [kg stem biomass / plant * cm^(-b2Bs)] -------------!
   b1Bs_small(1:4)   = 0.0 
   b1Bs_small(5)     = 1.0e-5
   b1Bs_small(6)     = 0.147
   b1Bs_small(7)     = 0.147
   b1Bs_small(8)     = 0.1617
   b1Bs_small(9)     = 0.02648
   b1Bs_small(10)    = 0.1617
   b1Bs_small(11)    = 0.235
   b1Bs_small(12:13) = 1.0e-5
   b1Bs_small(14:15) = 0.0 
   b1Bs_small(16)    = 0.0 
   b1Bs_small(17)    = 0.0
   !----- DBH-stem allometry slope [dimensionless]. ---------------------------------------!
   b2Bs_small(1:4)   = 0.0
   b2Bs_small(5)     = 1.0
   b2Bs_small(6)     = 2.238
   b2Bs_small(7)     = 2.238
   b2Bs_small(8)     = 2.1536
   b2Bs_small(9)     = 2.95954
   b2Bs_small(10)    = 2.4572
   b2Bs_small(11)    = 2.2518
   b2Bs_small(12:13) = 1.0
   b2Bs_small(14:15) = 0.0
   b2Bs_small(16)    = 0.0
   b2Bs_small(17)    = 0.0
   !---------------------------------------------------------------------------------------!
   !     The temperate PFTs use the same b1Bs and b2Bs for small and big trees, copy the   !
   ! values.                                                                               !
   !---------------------------------------------------------------------------------------!
   b1Bs_big(:) = b1Bs_small(:)
   b2Bs_big(:) = b2Bs_small(:)
   !------- Fill in the tropical PFTs, which are functions of wood density. ---------------!
   do ipft = 1, n_pft
      if (is_tropical(ipft)) then
         select case (iallom)
         case (0)
            !---- ED-2.1 allometry. -------------------------------------------------------!
            b1Bs_small(ipft) = exp(a1 + c1d * b1Ht(ipft) + d1d * log(rho(ipft)))
            b1Bs_big  (ipft) = exp(a1 + c1d * log(hgt_max(ipft)) + d1d * log(rho(ipft)))

            aux              = ( (a2d - a1) + b1Ht(ipft) * (c2d - c1d) + log(rho(ipft))    &
                               * (d2d - d1d)) * (1.0/log(dcrit))
            b2Bs_small(ipft) = C2B * b2d + c2d * b2Ht(ipft) + aux

            aux              = ( (a2d - a1) + log(hgt_max(ipft)) * (c2d - c1d)             &
                               + log(rho(ipft)) * (d2d - d1d)) * (1.0/log(dcrit))
            b2Bs_big  (ipft) = C2B * b2d + aux

         case (1)
            !---- Based on modified Chave et al. (2001) allometry. ------------------------!
            b1Bs_small(ipft) = C2B * exp(adead_small(1)) * rho(ipft) / adead_small(3)
            b2Bs_small(ipft) = adead_small(2)
            b1Bs_big(ipft)   = C2B * exp(adead_big(1))   * rho(ipft) / adead_big(3)
            b2Bs_big(ipft)   = adead_big(2)

         end select
      end if

      !----- Assigned for both cases, although it is really only needed for tropical. -----!
      bdead_crit(ipft) = dbh2bd(max_dbh(ipft),ipft)
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     DBH-crown allometry.                                                              !
   !---------------------------------------------------------------------------------------!
   !----- Intercept. ----------------------------------------------------------------------!
   b1Ca(1:17) = 2.490154
   !----- Slope.  -------------------------------------------------------------------------!
   b2Ca(1:17) = 0.8068806
   !---------------------------------------------------------------------------------------!

   if (write_allom) then
      open (unit=18,file=trim(allom_file),status='replace',action='write')
      write(unit=18,fmt='(209a)') ('-',n=1,209)
      write(unit=18,fmt='(18(1x,a))') '         PFT','    Tropical','         Rho'         &
                                     ,'        b1Ht','        b2Ht','        b1Bl'         &
                                     ,'        b2Bl','  b1Bs_Small','  b2Bs_Small'         &
                                     ,'    b1Bs_Big','    b1Bs_Big','        b1Bl'         &
                                     ,'        b2Bl','     Hgt_min','     Hgt_max'         &
                                     ,'     Min_DBH','     Max_DBH','   Bdead_CRT'
      write(unit=18,fmt='(209a)') ('-',n=1,209)
      do ipft=1,n_pft
         write (unit=18,fmt='(8x,i5,12x,l1,16(1x,es12.5))')                                &
                        ipft,is_tropical(ipft),rho(ipft),b1Ht(ipft),b2Ht(ipft),b1Bl(ipft)  &
                       ,b2Bl(ipft),b1Bs_small(ipft),b2Bs_small(ipft),b1Bs_big(ipft)        &
                       ,b2Bs_big(ipft),b1Ca(ipft),b2Ca(ipft),hgt_min(ipft),hgt_max(ipft)   &
                       ,min_dbh(ipft),max_dbh(ipft),bdead_crit(ipft)
      end do
      write(unit=18,fmt='(209a)') ('-',n=1,209)
      close(unit=18,status='keep')
   end if

   !---------------------------------------------------------------------------------------!
   !    Define the branching parameters, following Järvelä (2004)                          !
   !---------------------------------------------------------------------------------------!
   !----- Branching ratio -----------------------------------------------------------------!
   rbranch(1)     = 4.24
   rbranch(2:4)   = 4.23
   rbranch(5)     = 4.24
   rbranch(6:8)   = 4.44
   rbranch(9:11)  = 4.24
   rbranch(12:15) = 4.24
   rbranch(16)    = 4.24
   rbranch(17)    = 4.44
   !----- Diameter ratio ------------------------------------------------------------------!
   rdiamet(1)     = 5.00
   rdiamet(2:4)   = 1.86
   rdiamet(5)     = 5.00
   rdiamet(6:8)   = 2.04
   rdiamet(9:11)  = 1.86
   rdiamet(12:15) = 5.00
   rdiamet(16)    = 5.00
   rdiamet(17)    = 2.04
   !----- Length ratio. Järvelä used rdiamet^2/3, so do we... -----------------------------!
   rlength(1:17)  = rdiamet(1:17)**twothirds
   !----- Minimum diameter to consider [cm]. ----------------------------------------------!
   diammin(1:17)  = 1.0
   !----- Number of trunks.  Usually this is 1. -------------------------------------------!
   ntrunk(1:17)   = 1.0
   
   !---------------------------------------------------------------------------------------!
   !     The following variables are used to fit a smooth curve in the (sparse) values     !
   ! provided by Conijn (1995). This should be definitely improved...  The fitting curve   !
   ! is a + b*erf(c*bbranch+d)
   !---------------------------------------------------------------------------------------!
   conijn_a(1)     = 1.0
   conijn_a(2:4)   = 0.96305883
   conijn_a(5)     = 1.0
   conijn_a(6:11)  = 0.96305883
   conijn_a(12:15) = 1.0
   conijn_a(16)    = 1.0
   conijn_a(17)    = 0.96305883

   conijn_b(1)     = 0.0
   conijn_b(2:4)   = -0.7178682
   conijn_b(5)     = 0.0
   conijn_b(6:11)  = -0.7178682
   conijn_b(12:15) = 0.0
   conijn_b(16)    = 0.0
   conijn_b(17)    = -0.7178682

   conijn_c(1)     = 0.0
   conijn_c(2:4)   = 0.00490734
   conijn_c(5)     = 0.0
   conijn_c(6:11)  = 0.00490734
   conijn_c(12:15) = 0.0
   conijn_c(16)    = 0.0
   conijn_c(17)    = 0.00490734

   conijn_d(1)     = 0.0
   conijn_d(2:4)   = -0.0456370
   conijn_d(5)     = 0.0
   conijn_d(6:11)  = -0.0456370
   conijn_d(12:15) = 0.0
   conijn_d(16)    = 0.0
   conijn_d(17)    = -0.0456370
   return
end subroutine init_pft_alloc_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_nitro_params()

use pft_coms, only: c2n_leaf, Vm0, SLA, &
     c2n_slow,c2n_structural,c2n_storage,c2n_stem,l2n_stem, &
     C2B,plant_N_supply_scale

implicit none

c2n_slow       = 10.0  ! Carbon to Nitrogen ratio, slow pool.
c2n_structural = 150.0 ! Carbon to Nitrogen ratio, structural pool.
c2n_storage    = 150.0 ! Carbon to Nitrogen ratio, storage pool.
c2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.
l2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.


plant_N_supply_scale = 0.5 

c2n_leaf(1:5)    = 1000.0 / ((0.11289 + 0.12947 *   Vm0(1:5)) * SLA(1:5)  )
c2n_leaf(6)      = 1000.0 / ((0.11289 + 0.12947 *     15.625) * SLA(6)    )
c2n_leaf(7)      = 1000.0 / ((0.11289 + 0.12947 *     15.625) * SLA(7)    )
c2n_leaf(8)      = 1000.0 / ((0.11289 + 0.12947 *       6.25) * SLA(8)    )
c2n_leaf(9)      = 1000.0 / ((0.11289 + 0.12947 *      18.25) * SLA(9)    )
c2n_leaf(10)     = 1000.0 / ((0.11289 + 0.12947 *     15.625) * SLA(10)   )
c2n_leaf(11)     = 1000.0 / ((0.11289 + 0.12947 *       6.25) * SLA(11)   )
c2n_leaf(12:15)  = 1000.0 / ((0.11289 + 0.12947 * Vm0(12:15)) * SLA(12:15))
c2n_leaf(16:17)  = 1000.0 / ((0.11289 + 0.12947 * Vm0(16:17)) * SLA(16:17))



return
end subroutine init_pft_nitro_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine sets up some PFT and leaf dependent properties.                        !
!------------------------------------------------------------------------------------------!
subroutine init_pft_leaf_params()
   use rk4_coms       , only : ibranch_thermo       ! ! intent(in)
   use pft_coms       , only : phenology            & ! intent(out)
                             , clumping_factor      & ! intent(out)
                             , b1Tht                & ! intent(out)
                             , b2Tht                & ! intent(out)
                             , c_grn_leaf_dry       & ! intent(out)
                             , c_ngrn_biom_dry      & ! intent(out)
                             , wat_dry_ratio_grn    & ! intent(out)
                             , wat_dry_ratio_ngrn   & ! intent(out)
                             , delta_c              ! ! intent(out)
   use consts_coms    , only : t3ple                ! ! intent(out) 
   use phenology_coms , only :iphen_scheme

   implicit none

   select case (iphen_scheme)
   case (-1)
      phenology(1:8)   = 0
      phenology(9:11)  = 2
      phenology(12:17) = 0
   case (0,1)
      phenology(1)     = 1
      phenology(2:4)   = 1
      phenology(5)     = 1
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 1
      phenology(16)    = 1
      phenology(17)    = 0
   case (2)
      phenology(1)     = 4
      phenology(2:4)   = 4
      phenology(5)     = 4
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 4
      phenology(16)    = 4
      phenology(17)    = 0
   case (3)
      phenology(1)     = 4
      phenology(2:4)   = 3
      phenology(5)     = 4
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 4
      phenology(16)    = 4
      phenology(17)    = 0
   end select

   clumping_factor(1)     = 1.000d0
   clumping_factor(2:4)   = 7.350d-1
   clumping_factor(5)     = 8.400d-1
   clumping_factor(6:8)   = 7.350d-1
   clumping_factor(9:11)  = 8.400d-1
   clumping_factor(12:13) = 8.400d-1
   clumping_factor(14:15) = 1.000d0
   clumping_factor(16)    = 8.400d-1
   clumping_factor(17)    = 7.350d-1

   !---------------------------------------------------------------------------------------!
   !      The following parameters are second sources found in Gu et al. (2007)            !
   !---------------------------------------------------------------------------------------!
   c_grn_leaf_dry(1:17)      = 3218.0    ! Jones 1992  J/(kg K)
   c_ngrn_biom_dry(1:17)     = 1256.0    ! Forest Products Laboratory 
   wat_dry_ratio_grn(1:17)   = 2.5       ! 
   !wat_dry_ratio_grn(1:17)   = 1.5       ! Ceccato et al. 2001
   wat_dry_ratio_ngrn(1:17)  = 0.7       ! Forest Products Laboratory
   !---------------------------------------------------------------------------------------!
   !     Delta-c is found using the second term of the RHS of equation 5, assuming         !
   ! T=T3ple.  This is a simplification, but the specific heat usually varies by 3J/kg/K   !
   ! between 173K and 341K, so removing the dependence on temperature is not that bad      !
   ! assumption.                                                                           !
   !---------------------------------------------------------------------------------------!
   delta_c(1:17) = 100. * wat_dry_ratio_ngrn(1:17)                                         &
                 * (-0.06191 + 2.36e-4 * t3ple - 1.33e-2 * wat_dry_ratio_ngrn(1:17))

   !---------------------------------------------------------------------------------------!
   !     These are used to compute the height of the first branch, which we assume it is   !
   ! also the relative crown depth.                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Intercept. ----------------------------------------------------------------------!
   b1Tht(1)     = 0.01 ! 0.01
   b1Tht(2:4)   = 0.70 ! 0.4359
   b1Tht(5)     = 0.01 ! 0.01
   b1Tht(6:11)  = 0.70 ! 0.4359
   b1Tht(12:16) = 0.01 ! 0.01
   b1Tht(17)    = 0.70 ! 0.4359
   !----- Slope. ----------------------------------------------------------------------!
   b2Tht(1)     = 1.0 ! 1.0
   b2Tht(2:4)   = 1.0 ! 0.878
   b2Tht(5)     = 1.0 ! 1.0
   b2Tht(6:11)  = 1.0 ! 0.878
   b2Tht(12:16) = 1.0 ! 1.0
   b2Tht(17)    = 1.0 ! 0.878
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_leaf_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine sets some reproduction-related parameters.                            !
!------------------------------------------------------------------------------------------!
subroutine init_pft_repro_params()

   use pft_coms , only : r_fract            & ! intent(out)
                       , seed_rain          & ! intent(out)
                       , nonlocal_dispersal & ! intent(out)
                       , repro_min_h        ! ! intent(out)
   implicit none

   r_fract(1)                = 1.0
   r_fract(2:4)              = 0.3
   r_fract(5)                = 1.0
   r_fract(6:11)             = 0.3
   r_fract(12:15)            = 0.3
   r_fract(16)               = 1.0
   r_fract(17)               = 0.3

   seed_rain(1:17)           = 0.01

   nonlocal_dispersal(1:5)   = 1.0
   nonlocal_dispersal(6:7)   = 0.766
   nonlocal_dispersal(8)     = 0.001
   nonlocal_dispersal(9)     = 1.0
   nonlocal_dispersal(10)    = 0.325
   nonlocal_dispersal(11)    = 0.074
   nonlocal_dispersal(16)    = 1.0
   nonlocal_dispersal(17)    = 0.766

   repro_min_h(1)            = 0.0
   repro_min_h(2:4)          = 5.0
   repro_min_h(5)            = 0.0
   repro_min_h(6:11)         = 5.0
   repro_min_h(12:15)        = 0.0
   repro_min_h(16)           = 0.0
   repro_min_h(17)           = 5.0

   return
end subroutine init_pft_repro_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some variables that depend on the definition of other     !
! PFT parameters.  As such, this should be the last init_pft subroutine to be called.      !
!------------------------------------------------------------------------------------------!
subroutine init_pft_derived_params()
   use decomp_coms          , only : f_labile             ! ! intent(in)
   use ed_max_dims          , only : n_pft                & ! intent(in)
                                   , str_len              ! ! intent(in)
   use consts_coms          , only : onesixth             & ! intent(in)
                                   , twothirds            ! ! intent(in)
   use pft_coms             , only : init_density         & ! intent(in)
                                   , c2n_leaf             & ! intent(in)
                                   , c2n_stem             & ! intent(in)
                                   , b1Ht                 & ! intent(in)
                                   , b2Ht                 & ! intent(in)
                                   , hgt_min              & ! intent(in)
                                   , hgt_ref              & ! intent(in)
                                   , q                    & ! intent(in)
                                   , qsw                  & ! intent(in)
                                   , sla                  & ! intent(in)
                                   , pft_name16           & ! intent(in)
                                   , hgt_max              & ! intent(in)
                                   , max_dbh              & ! intent(in)
                                   , min_recruit_size     & ! intent(out)
                                   , min_cohort_size      & ! intent(out)
                                   , negligible_nplant    & ! intent(out)
                                   , c2n_recruit          & ! intent(out)
                                   , lai_min              ! ! intent(out)
   use phenology_coms       , only : elongf_min           ! ! intent(in)
   use allometry            , only : h2dbh                & ! function
                                   , dbh2h                & ! function
                                   , dbh2bl               & ! function
                                   , dbh2bd               ! ! function
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer                           :: ipft
   real                              :: dbh
   real                              :: huge_dbh
   real                              :: huge_height
   real                              :: balive_min
   real                              :: bleaf_min
   real                              :: bdead_min
   real                              :: balive_max
   real                              :: bleaf_max
   real                              :: bdead_max
   real                              :: min_plant_dens
   logical               , parameter :: print_zero_table = .false.
   character(len=str_len), parameter :: zero_table_fn    = 'minimum.size.txt'
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The minimum recruitment size and the recruit carbon to nitrogen ratio.  Both      !
   ! parameters actually depend on which PFT we are solving, since grasses always have     !
   ! significantly less biomass.                                                           !
   !---------------------------------------------------------------------------------------!
   if (print_zero_table) then
      open  (unit=61,file=trim(zero_table_fn),status='replace',action='write')
      write (unit=61,fmt='(18(a,1x))')                '  PFT',        'NAME            '   &
                                              ,'     HGT_MIN','         DBH'               &
                                              ,'   BLEAF_MIN','   BDEAD_MIN'               &
                                              ,'  BALIVE_MIN','   BLEAF_MAX'               &
                                              ,'   BDEAD_MAX','  BALIVE_MAX'               &
                                              ,'   INIT_DENS','MIN_REC_SIZE'               &
                                              ,'MIN_COH_SIZE',' NEGL_NPLANT'               &
                                              ,'         SLA','     LAI_MIN'               &
                                              ,'     HGT_MAX','     MAX_DBH'
   end if
   min_plant_dens = onesixth * minval(init_density)
   do ipft = 1,n_pft

      !----- Find the DBH and carbon pools associated with a newly formed recruit. --------!
      dbh        = h2dbh(hgt_min(ipft),ipft)
      bleaf_min  = dbh2bl(dbh,ipft)
      bdead_min  = dbh2bd(dbh,ipft)
      balive_min = bleaf_min * (1.0 + q(ipft) + qsw(ipft) * hgt_min(ipft))

      !------------------------------------------------------------------------------------!
      !   Find the maximum bleaf and bdead supported.  This is to find the negligible      !
      ! nplant so we ensure that the cohort is always terminated if its mortality rate is  !
      ! very high.                                                                         !
      !------------------------------------------------------------------------------------!
      huge_dbh    = 3. * max_dbh(ipft)
      huge_height = dbh2h(ipft, max_dbh(ipft))
      bleaf_max   = dbh2bl(huge_dbh,ipft)
      bdead_max   = dbh2bd(huge_dbh,ipft)
      balive_max  = bleaf_max * (1.0 + q(ipft) + qsw(ipft) * huge_height)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    The definition of the minimum recruitment size is the minimum amount of biomass !
      ! in kgC/m² is available for new recruits.  For the time being we use the near-bare  !
      ! ground state value as the minimum recruitment size, but this may change depending  !
      ! on how well it goes.                                                               !
      !------------------------------------------------------------------------------------!
      min_recruit_size(ipft) = min_plant_dens * (bdead_min + balive_min)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Minimum size (measured as biomass of living and structural tissues) allowed in  !
      ! a cohort.  Cohorts with less biomass than this are going to be terminated.         !
      !------------------------------------------------------------------------------------! 
      min_cohort_size(ipft)  = 0.1 * min_recruit_size(ipft)
      !------------------------------------------------------------------------------------! 



      !------------------------------------------------------------------------------------! 
      !    The following variable is the absolute minimum cohort population that a cohort  !
      ! can have.  This should be used only to avoid nplant=0, but IMPORTANT: this will    !
      ! lead to a ridiculously small cohort almost guaranteed to be extinct and SHOULD BE  !
      ! USED ONLY IF THE AIM IS TO ELIMINATE THE COHORT.                                   !
      !------------------------------------------------------------------------------------! 
      negligible_nplant(ipft) = onesixth * min_cohort_size(ipft) / (bdead_max + balive_max)
      !------------------------------------------------------------------------------------! 


      !----- Find the recruit carbon to nitrogen ratio. -----------------------------------!
      c2n_recruit(ipft)      = (balive_min + bdead_min)                                    &
                             / (balive_min * ( f_labile(ipft) / c2n_leaf(ipft)             &
                             + (1.0 - f_labile(ipft)) / c2n_stem(ipft))                    &
                             + bdead_min/c2n_stem(ipft))
      !------------------------------------------------------------------------------------! 



      !------------------------------------------------------------------------------------!
      !     The minimum LAI is the LAI of a plant at the minimum cohort size that is at    !
      ! the minimum elongation factor that supports leaves.                                !
      !------------------------------------------------------------------------------------!
      lai_min(ipft) = 0.1 * min_plant_dens * sla(ipft) * bleaf_min * elongf_min
      !------------------------------------------------------------------------------------!


      if (print_zero_table) then
         write (unit=61,fmt='(i5,1x,a16,1x,16(es12.5,1x))')                                &
                                                     ipft,pft_name16(ipft),hgt_min(ipft)   &
                                                    ,dbh,bleaf_min,bdead_min,balive_min    &
                                                    ,bleaf_max,bdead_max,balive_max        &
                                                    ,init_density(ipft)                    &
                                                    ,min_recruit_size(ipft)                &
                                                    ,min_cohort_size(ipft)                 &
                                                    ,negligible_nplant(ipft)               &
                                                    ,sla(ipft),lai_min(ipft)               &
                                                    ,hgt_max(ipft),max_dbh(ipft)
      end if
      !------------------------------------------------------------------------------------!
   end do

   if (print_zero_table) then
      close (unit=61,status='keep')
   end if

   return
end subroutine init_pft_derived_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_disturb_params

   use disturb_coms , only : sm_fire                  & ! intent(in)
                           , min_new_patch_area       & ! intent(out)
                           , treefall_hite_threshold  & ! intent(out)
                           , forestry_on              & ! intent(out)
                           , agriculture_on           & ! intent(out)
                           , plantation_year          & ! intent(out)
                           , plantation_rotation      & ! intent(out)
                           , mature_harvest_age       & ! intent(out)
                           , fire_dryness_threshold   & ! intent(out)
                           , fire_smoist_threshold    & ! intent(out)
                           , fire_smoist_depth        & ! intent(out)
                           , k_fire_first             & ! intent(out)
                           , fire_parameter           & ! intent(out)
                           , min_plantation_frac      & ! intent(out)
                           , max_plantation_dist      ! ! intent(out)
   use consts_coms  , only : erad                     & ! intent(in)
                           , pio180                   ! ! intent(in)
   implicit none
   
   !----- Minimum area that a patch must have to be created. ------------------------------!
   min_new_patch_area = 0.005

   !----- Only trees above this height create a gap when they fall. -----------------------!
   treefall_hite_threshold = 10.0 

   !----- Set to 1 if to do forest harvesting. --------------------------------------------!
   forestry_on = 0

   !----- Set to 1 if to do agriculture. --------------------------------------------------!
   agriculture_on = 0

   !----- Earliest year at which plantations occur. ---------------------------------------!
   plantation_year = 1960 

   !----- Number of years that a plantation requires to reach maturity. -------------------!
   plantation_rotation = 25.0

   !----- Years that a non-plantation patch requires to reach maturity. -------------------!
   mature_harvest_age = 50.0 
   
   !---------------------------------------------------------------------------------------!
   !     If include_fire is 1, then fire may occur if total (ground + underground) water   !
   ! converted to meters falls below this threshold.                                       !
   !---------------------------------------------------------------------------------------!
   fire_dryness_threshold = 0.2

   !---------------------------------------------------------------------------------------!
   !     If include_fire is 2, then fire may occur if total (ground + underground) water   !
   ! falls below a threshold defined by the total water of a soil column with average soil !
   ! moisture equal to soilcp + (slmsts-soilcp) * fire_smoist_threshold [m3_H2O/m3_gnd]    !
   ! would have.                                                                           !
   !---------------------------------------------------------------------------------------!
   fire_smoist_threshold = sm_fire

   !----- Maximum depth that will be considered in the average soil -----------------------!
   fire_smoist_depth     = -1.0

   !----- Dimensionless parameter controlling speed of fire spread. -----------------------!
   fire_parameter = 1.0

   !----- Minimum plantation fraction to consider the site a plantation. ------------------!
   min_plantation_frac = 0.125

   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon.  The value below is 1.25 degree at the Equator.  !
   !---------------------------------------------------------------------------------------!
   max_plantation_dist = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!

   return

end subroutine init_disturb_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_physiology_params()
   use physiology_coms, only : c34smin_lint_co2  & ! intent(out)
                             , c34smax_lint_co2  & ! intent(out)
                             , c34smax_gsw       & ! intent(out)
                             , gbh_2_gbw         & ! intent(out)
                             , gbw_2_gbc         & ! intent(out)
                             , gsw_2_gsc         & ! intent(out)
                             , gsc_2_gsw         & ! intent(out)
                             , tarrh             & ! intent(out)
                             , tarrhi            & ! intent(out)
                             , compp_refkin      & ! intent(out)
                             , compp_ecoeff      & ! intent(out)
                             , vm_ecoeff         & ! intent(out)
                             , vm_tempcoeff      & ! intent(out)
                             , kco2_prefac       & ! intent(out)
                             , kco2_ecoeff       & ! intent(out)
                             , ko2_prefac        & ! intent(out)
                             , ko2_ecoeff        & ! intent(out)
                             , klowco2           & ! intent(out)
                             , o2_ref            & ! intent(out)
                             , par_twilight_min  & ! intent(out)
                             , c34smin_lint_co28 & ! intent(out)
                             , c34smax_lint_co28 & ! intent(out)
                             , c34smax_gsw8      & ! intent(out)
                             , gbh_2_gbw8        & ! intent(out)
                             , gbw_2_gbc8        & ! intent(out)
                             , gsw_2_gsc8        & ! intent(out)
                             , gsc_2_gsw8        & ! intent(out)
                             , tarrh8            & ! intent(out)
                             , tarrhi8           & ! intent(out)
                             , compp_refkin8     & ! intent(out)
                             , compp_ecoeff8     & ! intent(out)
                             , vm_ecoeff8        & ! intent(out)
                             , vm_tempcoeff8     & ! intent(out)
                             , kco2_prefac8      & ! intent(out)
                             , kco2_ecoeff8      & ! intent(out)
                             , ko2_prefac8       & ! intent(out)
                             , ko2_ecoeff8       & ! intent(out)
                             , klowco28          & ! intent(out)
                             , par_twilight_min8 & ! intent(out)
                             , o2_ref8           & ! intent(out)
                             , print_photo_debug & ! intent(out)
                             , photo_prefix      ! ! intent(out)
   use consts_coms    , only : umol_2_mol        & ! intent(in)
                             , t00               & ! intent(in)
                             , mmdoc             & ! intent(in)
                             , Watts_2_Ein       ! ! intent(in)
   implicit none



   !---------------------------------------------------------------------------------------!
   !     Bounds for internal carbon and water stomatal conductance.                        !
   !---------------------------------------------------------------------------------------!
   c34smin_lint_co2 = 0.5   * umol_2_mol ! Minimum carbon dioxide concentration [  mol/mol]
   c34smax_lint_co2 = 1200. * umol_2_mol ! Maximum carbon dioxide concentration [  mol/mol]
   c34smax_gsw      = 1.e+2              ! Max. stomatal conductance (water)    [ mol/m²/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Many parameters in the model are temperature-dependent and utilise a modified     !
   ! Arrhenius function to determine this dependence.  For that to work, reference values  !
   ! at a given temperature (tarrh, in Kelvin). 
   !---------------------------------------------------------------------------------------!
   tarrh        = 15.0+t00
   tarrhi       = 1./tarrh
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following parameter is the concentration of oxygen, assumed constant through- !
   ! out the integration.  The value is in mol/mol.                                        !
   !---------------------------------------------------------------------------------------!
   o2_ref       = 0.209
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The next two variables are the parameters for the compensation point.             !
   !---------------------------------------------------------------------------------------!
   compp_refkin =  4500.  ! Reference ratio of "kinetic parameters describing the
                          !    partioning of enzyme activity to carboxylase or
                          !    oxylase function" (F96)                           [    ----]
   compp_ecoeff = -5000.  ! "Activation energy" parameter                        [       K]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      The next two variables are the parameter for the maximum capacity of Rubisco to  !
   ! perform the carboxylase function.  The reference value is PFT-dependent and the       !
   ! values are assigned in init_pft_photo_params.                                         !
   !---------------------------------------------------------------------------------------!
   vm_ecoeff    = 3000.     ! Reference exponential coeff. for carboxylase fctn. [       K]
   vm_tempcoeff = 0.4       ! Coefficient to control the temperature dependence  [     ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      The following terms are used to find the Michaelis-Menten coefficient for CO2.   !
   !---------------------------------------------------------------------------------------!
   kco2_prefac    = 150. * umol_2_mol ! Reference CO2 concentration              [ mol/mol]
   kco2_ecoeff    = 6000.             ! Reference exponential coefficient        [       K]
   !---------------------------------------------------------------------------------------! 



   !---------------------------------------------------------------------------------------! 
   !     These terms are used to find the Michaelis-Mentencoefficient for O2.              !
   !---------------------------------------------------------------------------------------!
   ko2_prefac    = 0.250     ! Reference O2 concentration.                        [ mol/mol]
   ko2_ecoeff    =  1400.    ! Reference exponential coefficient                  [       K]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following parameter is the k coefficient in Foley et al. (1996) that is used   !
   ! to determine the CO2-limited photosynthesis for C4 grasses.  Notice that Foley et al. !
   ! (1996) didn't correct for molar mass (the mmdoc term here).                           !
   !---------------------------------------------------------------------------------------!
   klowco2      = 18000. * mmdoc ! coefficient for low CO2                       [ mol/mol]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are constants obtained in Leuning et al. (1995) and Collatz et al. (1991)   !
   ! to convert different conductivities.                                                  !
   !---------------------------------------------------------------------------------------!
   gbh_2_gbw  = 1.075           ! heat  to water  - leaf boundary layer
   gbw_2_gbc  = 1.0 / 1.4       ! water to carbon - leaf boundary layer
   gsw_2_gsc  = 1.0 / 1.6       ! water to carbon - stomata
   gsc_2_gsw  = 1./gsw_2_gsc    ! carbon to water - stomata
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is the minimum threshold for the photosynthetically active radiation, in     !
   ! µmol/m²/s to consider non-night time conditions (day time or twilight).               !
   !---------------------------------------------------------------------------------------!
   par_twilight_min = 0.5 * Watts_2_Ein ! Minimum non-nocturnal PAR.             [mol/m²/s]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Find the double precision version of the variables above.                          !
   !---------------------------------------------------------------------------------------!
   c34smin_lint_co28 = dble(c34smin_lint_co2)
   c34smax_lint_co28 = dble(c34smax_lint_co2)
   c34smax_gsw8      = dble(c34smax_gsw     )
   gbh_2_gbw8        = dble(gbh_2_gbw       )
   gbw_2_gbc8        = dble(gbw_2_gbc       )
   gsw_2_gsc8        = dble(gsw_2_gsc       )
   gsc_2_gsw8        = dble(gsc_2_gsw       )
   tarrh8            = dble(tarrh           )
   tarrhi8           = dble(tarrhi          )
   compp_refkin8     = dble(compp_refkin    )
   compp_ecoeff8     = dble(compp_ecoeff    )
   vm_ecoeff8        = dble(vm_ecoeff       )
   vm_tempcoeff8     = dble(vm_tempcoeff    )
   kco2_prefac8      = dble(kco2_prefac     )
   kco2_ecoeff8      = dble(kco2_ecoeff     )
   ko2_prefac8       = dble(ko2_prefac      )
   ko2_ecoeff8       = dble(ko2_ecoeff      )
   klowco28          = dble(klowco2         )
   o2_ref8           = dble(o2_ref          )
   par_twilight_min8 = dble(par_twilight_min)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters that control debugging output.                                         !
   !---------------------------------------------------------------------------------------!
   !----- I should print detailed debug information. --------------------------------------!
   print_photo_debug = .false.
   !----- File name prefix for the detailed information in case of debugging. -------------!
   photo_prefix      = 'photo_state_'
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_physiology_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_hydro_coms

  use hydrology_coms,only:useTOPMODEL,useRUNOFF,HydroOutputPeriod, &
       MoistRateTuning,MoistSatThresh,Moist_dWT,FracLiqRunoff, &
       GrassLAIMax,inverse_runoff_time
  use ed_misc_coms, only: ied_init_mode

  implicit none

  if(ied_init_mode == 3)then
     ! Signifies a restart from an ED2 history file
     useTOPMODEL = 1
     useRUNOFF   = 0
  else
     useTOPMODEL = 0
     useRUNOFF = 0
  endif

  HydroOutputPeriod = 96 !! multiples of dtlsm

  MoistRateTuning = 1.0    

  MoistSatThresh = 0.95     
  
  Moist_dWT = 2.0 

  FracLiqRunoff = 0.5
  
  GrassLAImax = 4.0

  inverse_runoff_time = 0.1

  return
end subroutine init_hydro_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_soil_coms
   use soil_coms      , only : ed_nstyp              & ! intent(in)
                             , isoilflg              & ! intent(in)
                             , nslcon                & ! intent(in)
                             , slxclay               & ! intent(in)
                             , slxsand               & ! intent(in)
                             , soil                  & ! intent(in)
                             , betapower             & ! intent(in)
                             , soil_class            & ! type
                             , soil8                 & ! intent(out)
                             , water_stab_thresh     & ! intent(out)
                             , snowmin               & ! intent(out)
                             , dewmax                & ! intent(out)
                             , soil_rough            & ! intent(out)
                             , snow_rough            & ! intent(out)
                             , tiny_sfcwater_mass    & ! intent(out)
                             , infiltration_method   & ! intent(out)
                             , soil_rough8           & ! intent(out)
                             , snow_rough8           & ! intent(out)
                             , betapower8            ! ! intent(out)

   use grid_coms      , only : ngrids                ! ! intent(in)
   use consts_coms    , only : grav                  & ! intent(in)
                             , hr_sec                & ! intent(in)
                             , day_sec               ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer            :: nsoil
   integer            :: ifm
   !----- Local constants. ----------------------------------------------------------------!
   real   , parameter :: fieldcp_K  = 0.1 ! hydraulic conduct. at field capacity   [mm/day]
   real   , parameter :: soilcp_MPa = 3.1 ! soil-water potential for air dry soil  [   MPa]
   real   , parameter :: soilwp_MPa = 1.5 ! soil-water potential at wilting point  [   MPa]
   !---------------------------------------------------------------------------------------!


   !----- Initialise some standard variables. ---------------------------------------------!
   water_stab_thresh   = 3.0    ! Minimum water mass to be considered stable     [   kg/m2]
   snowmin             = 3.0    ! Minimum snow mass needed to create a new layer [   kg/m2]
   dewmax              = 3.0e-5 ! Maximum dew flux rate (deprecated)             [ kg/m2/s]
   soil_rough          = 0.05   ! Soil roughness height                          [       m]
   snow_rough          = 0.001  ! Snowcover roughness height                     [       m]
   tiny_sfcwater_mass  = 1.0e-3 ! Minimum allowed mass in temporary layers       [   kg/m2]
   infiltration_method = 0      ! Infiltration method, used in rk4_derivs        [     0|1]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      UPDATED based on Cosby et al. 1984, soil clay and sand fractions from table 2,   !
   ! slpots, slmsts, slbs, and slcons calculated based on table 4 sfldcap calculated based !
   ! on Clapp and Hornberger equation 2 and a soil hydraulic conductivity of 0.1mm/day,    !
   ! soilcp (air dry soil moisture capacity) is calculated based on Cosby et al. 1984      !
   ! equation 1 and a soil-water potential of -3.1MPa.  soilwp (the wilting point soil     !
   ! moisture) is calculated based on Cosby et al 1985 equation 1 and a soil-water poten-  !
   ! tial of -1.5MPa (Equation 1 in Cosby is not correct it should be saturated moisture   !
   ! potential over moisture potential)  (NML, 2/2010)                                     !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! (1st line)          slpots        slmsts          slbs     slcpd        soilcp        !
   ! (2nd line)          soilwp        slcons       slcons0 soilcond0     soilcond1        !
   ! (3rd line)       soilcond2       sfldcap         xsand     xclay         xsilt        !
   ! (4th line)         xrobulk         slden                                              !
   !---------------------------------------------------------------------------------------!
   soil = (/                                                                               &
      !----- 1. Sand. ---------------------------------------------------------------------!
       soil_class( -0.049831046,     0.373250,     3.295000, 1421830.,  0.026183447        &
                 ,  0.032636854,  2.446421e-5,  0.000500000,   0.3000,       4.8000        &
                 ,      -2.7000,  0.132130936,        0.920,    0.030,        0.050        &
                 ,        1200.,        1600.                                      )       &
      !----- 2. Loamy sand. ---------------------------------------------------------------!
      ,soil_class( -0.067406224,     0.385630,     3.794500, 1395530.,  0.041560499        &
                 ,  0.050323046,  1.776770e-5,  0.000600000,   0.3000,       4.6600        &
                 ,      -2.6000,  0.155181959,        0.825,    0.060,        0.115        &
                 ,        1250.,        1600.                                      )       &
      !----- 3. Sandy loam. ---------------------------------------------------------------!
      ,soil_class( -0.114261521,     0.407210,     4.629000, 1350750.,  0.073495043        &
                 ,  0.085973722,  1.022660e-5,  0.000769000,   0.2900,       4.2700        &
                 ,      -2.3100,  0.194037750,        0.660,    0.110,        0.230        &
                 ,        1300.,        1600.                                      )       &
      !----- 4. Silt loam. ----------------------------------------------------------------!
      ,soil_class( -0.566500112,     0.470680,     5.552000, 1264080.,  0.150665475        &
                 ,  0.171711257,  2.501101e-6,  0.000010600,   0.2700,       3.4700        &
                 ,      -1.7400,  0.273082063,        0.200,    0.160,        0.640        &
                 ,        1400.,        1600.                                      )       &
      !----- 5. Loam. ---------------------------------------------------------------------!
      ,soil_class( -0.260075834,     0.440490,     5.646000,  1289630,  0.125192234        &
                 ,  0.142369513,  4.532431e-6,  0.002200000,   0.2800,       3.6300        &
                 ,      -1.8500,  0.246915025,        0.410,    0.170,        0.420        &
                 ,        1350.,        1600.                                      )       &
      !----- 6. Sandy clay loam. ----------------------------------------------------------!
      ,soil_class( -0.116869181,     0.411230,     7.162000, 1272490.,  0.136417267        &
                 ,  0.150969505,  6.593731e-6,  0.001500000,   0.2800,       3.7800        &
                 ,      -1.9600,  0.249629687,        0.590,    0.270,        0.140        &
                 ,        1350.,        1600.                                      )       &
      !----- 7. Silty clay loam. ----------------------------------------------------------!
      ,soil_class( -0.627769194,     0.478220,     8.408000, 1173020.,  0.228171947        &
                 ,  0.248747504,  1.435262e-6,  0.000107000,   0.2600,       2.7300        &
                 ,      -1.2000,  0.333825332,        0.100,    0.340,        0.560        &
                 ,        1500.,        1600.                                      )       &
      !----- 8. Clayey loam. --------------------------------------------------------------!
      ,soil_class( -0.281968114,     0.446980,     8.342000, 1204260.,  0.192624431        &
                 ,  0.210137962,  2.717260e-6,  0.002200000,   0.2700,       3.2300        &
                 ,      -1.5600,  0.301335491,        0.320,    0.340,        0.340        &
                 ,        1450.,        1600.                                      )       &
      !----- 9. Sandy clay. ---------------------------------------------------------------!
      ,soil_class( -0.121283019,     0.415620,     9.538000, 1198500.,  0.182198910        &
                 ,  0.196607427,  4.314507e-6,  0.000002167,   0.2700,       3.3200        &
                 ,      -1.6300,  0.286363001,        0.520,    0.420,        0.060        &
                 ,        1450.,        1600.                                      )       &
      !----- 10. Silty clay. --------------------------------------------------------------!
      ,soil_class( -0.601312179,     0.479090,    10.461000, 1111830.,  0.263228486        &
                 ,  0.282143846,  1.055191e-6,  0.000001033,   0.2500,       2.5800        &
                 ,      -1.0900,  0.360319788,        0.060,    0.470,        0.470        &
                 ,        1650.,        1600.                                      )       &
      !----- 11. Clay. --------------------------------------------------------------------!
      ,soil_class( -0.299226464,     0.454400,    12.460000, 1076200.,  0.259868987        &
                 ,  0.275459057,  1.307770e-6,  0.000001283,   0.2500,       2.4000        &
                 ,      -0.9600,  0.353255209,        0.200,    0.600,        0.200        &
                 ,        1700.,        1600.                                      )       &
      !----- 12. Peat. --------------------------------------------------------------------!
      ,soil_class( -0.534564359,     0.469200,     6.180000,  874000.,  0.167047523        &
                 ,  0.187868805,  2.357930e-6,  0.000008000,   0.0600,       0.4600        &
                 ,       0.0000,  0.285709966,       0.2000,   0.2000,       0.6000        &
                 ,         500.,         300.                                      )       &
      !----- 13. Bedrock. -----------------------------------------------------------------!
      ,soil_class(    0.0000000,     0.000000,     0.000000, 2130000.,  0.000000000        &
                 ,  0.000000000,  0.000000e+0,  0.000000000,   4.6000,       0.0000        &
                 ,       0.0000,  0.000000001,       0.0000,   0.0000,       0.0000        &
                 ,           0.,           0.                                      )       &
      !----- 14. Silt. --------------------------------------------------------------------!
      ,soil_class( -1.047128548,     0.492500,     3.862500, 1293300.,  0.112299080        &
                 ,  0.135518820,  2.046592e-6,  0.000010600,   0.2700,       3.4700        &
                 ,      -1.7400,  0.245247642,        0.075,    0.050,        0.875        &
                 ,        1400.,        1600.                                      )       &
      !----- 15. Heavy clay. --------------------------------------------------------------!
      ,soil_class( -0.322106879,     0.461200,    15.630000,  976600.,  0.296806035        &
                 ,  0.310916364,  7.286705e-7,  0.000001283,   0.2500,       2.4000        &
                 ,      -0.9600,  0.382110712,        0.100,    0.800,        0.100        &
                 ,        1700.,        1600.                                      )       &
      !----- 16. Clayey sand. -------------------------------------------------------------!
      ,soil_class( -0.176502150,     0.432325,    11.230000, 1133075.,  0.221886929        &
                 ,  0.236704039,  2.426785e-6,  0.000001283,   0.2500,       2.4000        &
                 ,      -0.9600,  0.320146708,        0.375,    0.525,        0.100        &
                 ,        1700.,        1600.                                      )       &
      !----- 17. Clayey silt. -------------------------------------------------------------!
      ,soil_class( -0.438278332,     0.467825,    11.305000, 1097575.,  0.261376708        &
                 ,  0.278711303,  1.174982e-6,  0.000001283,   0.2500,       2.4000        &
                 ,      -0.9600,  0.357014719,        0.125,    0.525,        0.350        &
                 ,        1700.,        1600.                                      )       &
   /)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !******** Correct soil_class table using sand and clay fractions (if provided)  ********!
   !      Based on Cosby et al 1984, using table 4 and equation 1 (which is incorrect it   !
   ! should be saturated moisture potential over moisture potential).  NML 2/2010          !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids
      if ( isoilflg(ifm)==2 .and. slxclay > 0. .and. slxsand > 0. .and.                    &
           (slxclay + slxsand) <= 1. ) then
         soil(nslcon)%xsand   = slxsand
         soil(nslcon)%xclay   = slxclay
         soil(nslcon)%xsilt   = 1. - slxsand - slxclay

         !----- B exponent [unitless]. ----------------------------------------------------!
         soil(nslcon)%slbs    = 3.10 + 15.7*slxclay - 0.3*slxsand

         !----- Soil moisture potential at saturation [ m ]. ------------------------------!
         soil(nslcon)%slpots  = -1. * (10.**(2.17 - 0.63*slxclay - 1.58*slxsand)) * 0.01

         !----- Hydraulic conductivity at saturation [ m/s ]. -----------------------------!
         soil(nslcon)%slcons  = (10.**(-0.60 + 1.26*slxsand - 0.64*slxclay))               &
                              * 0.0254/hr_sec
         !----- Soil moisture at saturation [ m^3/m^3 ]. ----------------------------------!
         soil(nslcon)%slmsts  = (50.5 - 14.2*slxsand - 3.7*slxclay) / 100.
         !----- Soil field capacity[ m^3/m^3 ]. -------------------------------------------!
         soil(nslcon)%sfldcap = soil(nslcon)%slmsts                                        &
                              *  ( (fieldcp_K/1000./day_sec)/soil(nslcon)%slcons)          &
                              ** (1. / (2.*soil(nslcon)%slbs+3.))
         !----- Dry soil capacity (at -3.1MPa) [ m^3/m^3 ]. -------------------------------!
         soil(nslcon)%soilcp  = soil(nslcon)%slmsts                                        &
                              *  ( -1.*soil(nslcon)%slpots / (soilcp_MPa * 1000. / grav))  &
                              ** (1. / soil(nslcon)%slbs)                                         
         !----- Wilting point capacity (at -1.5MPa) [ m^3/m^3 ]. --------------------------!
         soil(nslcon)%soilwp  = soil(nslcon)%slmsts                                        &
                              *  ( -1.*soil(nslcon)%slpots / (soilwp_MPa * 1000. / grav))  &
                              ** ( 1. / soil(nslcon)%slbs)
         !---------------------------------------------------------------------------------!
      end if
   end do
   !---------------------------------------------------------------------------------------!







   !----- Here we fill soil8, which will be used in Runge-Kutta (double precision). -------!
   do nsoil=1,ed_nstyp
      soil8(nsoil)%slpots    = dble(soil(nsoil)%slpots   )
      soil8(nsoil)%slmsts    = dble(soil(nsoil)%slmsts   )
      soil8(nsoil)%slbs      = dble(soil(nsoil)%slbs     )
      soil8(nsoil)%slcpd     = dble(soil(nsoil)%slcpd    )
      soil8(nsoil)%soilcp    = dble(soil(nsoil)%soilcp   )
      soil8(nsoil)%soilwp    = dble(soil(nsoil)%soilwp   )
      soil8(nsoil)%slcons    = dble(soil(nsoil)%slcons   )
      soil8(nsoil)%slcons0   = dble(soil(nsoil)%slcons0  )
      soil8(nsoil)%soilcond0 = dble(soil(nsoil)%soilcond0)
      soil8(nsoil)%soilcond1 = dble(soil(nsoil)%soilcond1)
      soil8(nsoil)%soilcond2 = dble(soil(nsoil)%soilcond2)
      soil8(nsoil)%sfldcap   = dble(soil(nsoil)%sfldcap  )
      soil8(nsoil)%xsand     = dble(soil(nsoil)%xsand    )
      soil8(nsoil)%xclay     = dble(soil(nsoil)%xclay    )
      soil8(nsoil)%xsilt     = dble(soil(nsoil)%xsilt    )
      soil8(nsoil)%xrobulk   = dble(soil(nsoil)%xrobulk  )
      soil8(nsoil)%slden     = dble(soil(nsoil)%slden    )
   end do
   betapower8  = dble(betapower)
   soil_rough8 = dble(soil_rough)
   snow_rough8 = dble(snow_rough)
   return
end subroutine init_soil_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_phen_coms
   use consts_coms   , only : erad                     & ! intent(in)
                            , pio180                   ! ! intent(in)
   use phenology_coms, only : retained_carbon_fraction & ! intent(out)
                            , elongf_min               & ! intent(out)
                            , theta_crit               & ! intent(out)
                            , dl_tr                    & ! intent(out)
                            , st_tr1                   & ! intent(out)
                            , st_tr2                   & ! intent(out)
                            , phen_a                   & ! intent(out)
                            , phen_b                   & ! intent(out)
                            , phen_c                   & ! intent(out)
                            , rad_turnover_int         & ! intent(out)
                            , rad_turnover_slope       & ! intent(out)
                            , vm_tran                  & ! intent(out)
                            , vm_slop                  & ! intent(out)
                            , vm_amp                   & ! intent(out)
                            , vm_min                   & ! intent(out)
                            , max_phenology_dist       & ! intent(out)
                            , radint                   & ! intent(out)
                            , radslp                   & ! intent(out)
                            , thetacrit                   ! ! intent(in)
   implicit none

 
   retained_carbon_fraction = 0.5
   elongf_min               = 0.02
   theta_crit               = thetacrit
   dl_tr                    = 655.0
   st_tr1                   = 284.3
   st_tr2                   = 275.15
  
   phen_a                   = -68.0
   phen_b                   = 638.0
   phen_c                   = -0.01

   rad_turnover_int         = dble(radint)  !-11.3868
   rad_turnover_slope       = dble(radslp)  !0.0824

   vm_tran                  = 8.5
   vm_slop                  = 7.0
   vm_amp                   = 42.0
   vm_min                   = 18.0

   !---------------------------------------------------------------------------------------!
   !     This variable is the maximum distance between the coordinates of a prescribed     !
   ! phenology file and the actual polygon that we will still consider close enough to be  !
   ! representative.  If the user wants to run with prescribed phenology and the closest   !
   ! file is farther away from the polygon than the number below, the simulation will      !
   ! stop.                                                                                 !
   !---------------------------------------------------------------------------------------!
   max_phenology_dist       = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_phen_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns the fusion and splitting parameters.                         !
!------------------------------------------------------------------------------------------!
subroutine init_ff_coms
   use fusion_fission_coms, only : niter_patfus       & ! intent(out)
                                 , hgt_class          & ! intent(out)
                                 , fusetol            & ! intent(out)
                                 , fusetol_h          & ! intent(out)
                                 , lai_fuse_tol       & ! intent(out)
                                 , lai_tol            & ! intent(out)
                                 , ff_nhgt            & ! intent(out)
                                 , coh_tolerance_max  & ! intent(out)
                                 , dark_cumlai_min    & ! intent(out)
                                 , dark_cumlai_max    & ! intent(out)
                                 , dark_cumlai_mult   & ! intent(out)
                                 , sunny_cumlai_min   & ! intent(out)
                                 , sunny_cumlai_max   & ! intent(out)
                                 , sunny_cumlai_mult  & ! intent(out)
                                 , light_toler_min    & ! intent(out)
                                 , light_toler_max    & ! intent(out)
                                 , light_toler_mult   & ! intent(out)
                                 , fuse_relax         & ! intent(out)
                                 , print_fuse_details & ! intent(out)
                                 , fuse_prefix        ! ! intent(out)
   use consts_coms        , only : onethird           & ! intent(out)
                                 , twothirds          ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real              :: exp_patfus
   real              :: exp_hgtclass
   !---------------------------------------------------------------------------------------!

   fusetol           = 0.4
   fusetol_h         = 0.5
   lai_fuse_tol      = 0.8
   lai_tol           = 1.0
   ff_nhgt           = 7
   coh_tolerance_max = 10.0    ! Original 2.0

   !----- Define the number of height classes. --------------------------------------------!
   allocate (hgt_class(ff_nhgt))
   hgt_class( 1) =  2.0
   hgt_class( 2) =  5.0
   hgt_class( 3) =  9.0
   hgt_class( 4) = 14.0
   hgt_class( 5) = 20.0
   hgt_class( 6) = 26.0
   hgt_class( 7) = 35.0

   niter_patfus       = 25
   exp_patfus         = 1. / real(niter_patfus)

   dark_cumlai_min    = 5.0
   dark_cumlai_max    = 7.0
   sunny_cumlai_min   = 0.5
   sunny_cumlai_max   = 1.0
   light_toler_min    = 0.20
   light_toler_max    = twothirds
   sunny_cumlai_mult  = (sunny_cumlai_max/sunny_cumlai_min)**exp_patfus
   dark_cumlai_mult   = (dark_cumlai_min /dark_cumlai_max )**exp_patfus
   light_toler_mult   = (light_toler_max /light_toler_min )**exp_patfus

   fuse_relax        = .false.

   !----- The following flag switches detailed debugging on. ------------------------------!
   print_fuse_details = .false.
   fuse_prefix        = 'patch_fusion_'
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ff_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assigns various parameters for the Runge-Kutta solver.  It uses many  !
! values previously assigned in other parameter initialisation, so this should be the last !
! one called.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine init_rk4_params()
   use soil_coms      , only : water_stab_thresh      & ! intent(in)
                             , snowmin                & ! intent(in)
                             , tiny_sfcwater_mass     ! ! intent(in)
   use canopy_air_coms, only : dry_veg_lwater         & ! intent(in)
                             , fullveg_lwater         ! ! intent(in)
   use met_driver_coms, only : prss_min               & ! intent(in)
                             , prss_max               ! ! intent(in)
   use consts_coms    , only : wdnsi8                 ! ! intent(in)
   use rk4_coms       , only : rk4_tolerance          & ! intent(in)
                             , ibranch_thermo         & ! intent(in)
                             , maxstp                 & ! intent(out)
                             , rk4eps                 & ! intent(out)
                             , rk4eps2                & ! intent(out)
                             , rk4epsi                & ! intent(out)
                             , hmin                   & ! intent(out)
                             , print_diags            & ! intent(out)
                             , checkbudget            & ! intent(out)
                             , debug                  & ! intent(out)
                             , toocold                & ! intent(out)
                             , toohot                 & ! intent(out)
                             , lai_to_cover           & ! intent(out)
                             , hcapveg_ref            & ! intent(out)
                             , min_height             & ! intent(out)
                             , rk4min_veg_temp        & ! intent(out)
                             , rk4water_stab_thresh   & ! intent(out)
                             , rk4tiny_sfcw_mass      & ! intent(out)
                             , rk4tiny_sfcw_depth     & ! intent(out)
                             , rk4dry_veg_lwater      & ! intent(out)
                             , rk4fullveg_lwater      & ! intent(out)
                             , rk4snowmin             & ! intent(out)
                             , rk4min_can_temp        & ! intent(out)
                             , rk4max_can_temp        & ! intent(out)
                             , rk4min_can_shv         & ! intent(out)
                             , rk4max_can_shv         & ! intent(out)
                             , rk4min_can_rvap        & ! intent(out)
                             , rk4max_can_rvap        & ! intent(out)
                             , rk4min_can_rhv         & ! intent(out)
                             , rk4max_can_rhv         & ! intent(out)
                             , rk4min_can_co2         & ! intent(out)
                             , rk4max_can_co2         & ! intent(out)
                             , rk4min_soil_temp       & ! intent(out)
                             , rk4max_soil_temp       & ! intent(out)
                             , rk4min_veg_temp        & ! intent(out)
                             , rk4max_veg_temp        & ! intent(out)
                             , rk4min_veg_lwater      & ! intent(out)
                             , rk4min_sfcw_temp       & ! intent(out)
                             , rk4max_sfcw_temp       & ! intent(out)
                             , rk4min_sfcw_moist      & ! intent(out)
                             , rk4min_virt_moist      & ! intent(out)
                             , effarea_heat           & ! intent(out)
                             , effarea_evap           & ! intent(out)
                             , effarea_transp          & ! intent(out)
                             , leaf_intercept         & ! intent(out)
                             , force_idealgas         & ! intent(out)
                             , supersat_ok            & ! intent(out)
                             , record_err             & ! intent(out)
                             , print_detailed         & ! intent(out)
                             , print_thbnd            & ! intent(out)
                             , errmax_fout            & ! intent(out)
                             , sanity_fout            & ! intent(out)
                             , thbnds_fout            & ! intent(out)
                             , detail_pref            ! ! intent(out)
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Copying some variables to the Runge-Kutta counterpart (double precision).         !
   !---------------------------------------------------------------------------------------!
   rk4water_stab_thresh  = dble(water_stab_thresh )
   rk4tiny_sfcw_mass     = dble(tiny_sfcwater_mass)
   rk4dry_veg_lwater     = dble(dry_veg_lwater    )
   rk4fullveg_lwater     = dble(fullveg_lwater    )
   rk4snowmin            = dble(snowmin           )
   rk4tiny_sfcw_depth    = rk4tiny_sfcw_mass  * wdnsi8
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    The following variables control the performance of the Runge-Kutta and Euler       !
   ! integration schemes. Think twice before changing them...                              !
   !---------------------------------------------------------------------------------------!
   maxstp      = 100000000           ! Maximum number of intermediate steps. 
   rk4eps      = dble(rk4_tolerance) ! The desired accuracy.
   rk4epsi     = 1.d0/rk4eps         ! The inverse of desired accuracy.
   rk4eps2     = rk4eps**2           ! square of the accuracy
   hmin        = 1.d-7               ! The minimum step size.
   print_diags = .false.             ! Flag to print the diagnostic check.
   checkbudget = .true.              ! Flag to check CO2, water, and energy budgets every 
                                     !     time step and stop the run in case any of these 
                                     !     budgets don't close.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Miscellaneous constants used in rk4_derivs.                                       !
   !---------------------------------------------------------------------------------------!
   debug         = .false.  ! Verbose output for debug                             [   T|F]
   toocold       = 1.5315d2 ! Minimum temperature for saturation specific hum.     [     K]
   toohot        = 3.5315d2 ! Maximum temperature for saturation specific hum.     [     K]
   lai_to_cover  = 1.5d0    ! Canopies with LAI less than this number  are assumed to be 
                            !     open, ie, some fraction of the rain-drops can reach
                            !    the soil/litter layer unimpeded.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    These two parameter will scale the cohort heat capacity inside the RK4 integrator, !
   ! to avoid having patches with heat capacity that is way too small to be computational- !
   ! ly stable and resolvable in a fast way.  If you don't want this and want to use the   !
   ! nominal heat capacity, the laziest way to turn this off is by setting hcapveg_ref to  !
   ! a small number.  Don't set it to zero, otherwise you may have FPE issues.             !
   !---------------------------------------------------------------------------------------!
   hcapveg_ref         = 3.0d3            ! Reference heat capacity value          [J/m³/K]
   min_height          = 1.5d0            ! Minimum vegetation height              [     m]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Variables used to keep track on the error.                                        !
   !---------------------------------------------------------------------------------------!
   record_err     = .false.                  ! Compute and keep track of the errors.
   print_detailed = .false.                  ! Print detailed information about the thermo-
                                             !    dynamic state.  This will create one file
                                             !    for each patch, so it is not recommended 
                                             !    for simulations that span over one month.
   print_thbnd    = .false.                  ! Make a file with thermodynamic boundaries.
   errmax_fout    = 'error_max_count.txt'    ! File with the maximum error count 
   sanity_fout    = 'sanity_check_count.txt' ! File with the sanity check count
   thbnds_fout    = 'thermo_bounds.txt'      ! File with the thermodynamic boundaries.
   detail_pref    = 'thermo_state_'          ! Prefix for the detailed thermodynamic file
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Assigning some default values for the bounds at the sanity check.  Units are      !
   ! usually the standard, but a few of them are defined differently so they can be scaled !
   ! depending on the cohort and soil grid definitions.                                    !
   !---------------------------------------------------------------------------------------!
   rk4min_can_temp   =  1.8400d2  ! Minimum canopy    temperature               [        K]
   rk4max_can_temp   =  3.4100d2  ! Maximum canopy    temperature               [        K]
   rk4min_can_shv    =  1.0000d-8 ! Minimum canopy    specific humidity         [kg/kg_air]
   rk4max_can_shv    =  4.6000d-2 ! Maximum canopy    specific humidity         [kg/kg_air]
   rk4max_can_rhv    =  1.1000d0  ! Maximum canopy    relative humidity (**)    [      ---]
   rk4min_can_co2    =  1.0000d1  ! Minimum canopy    CO2 mixing ratio          [ µmol/mol]
   rk4max_can_co2    =  2.0000d3  ! Maximum canopy    CO2 mixing ratio          [ µmol/mol]
   rk4min_soil_temp  =  1.8400d2  ! Minimum soil      temperature               [        K]
   rk4max_soil_temp  =  3.5100d2  ! Maximum soil      temperature               [        K]
   rk4min_veg_temp   =  1.8400d2  ! Minimum leaf      temperature               [        K]
   rk4max_veg_temp   =  3.4100d2  ! Maximum leaf      temperature               [        K]
   rk4min_sfcw_temp  =  1.9315d2  ! Minimum snow/pond temperature               [        K]
   rk4max_sfcw_temp  =  3.4100d2  ! Maximum snow/pond temperature               [        K]
   !.......................................................................................!
   ! (**) Please, don't be too strict here.  The model currently doesn't have radiation    !
   !      fog, so supersaturation may happen.  This is a problem we may want to address in !
   !      the future, though...                                                            !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Compute the minimum and maximum mixing ratio based on the specific humidity.      !
   !---------------------------------------------------------------------------------------!
   rk4min_can_rvap = rk4min_can_shv / (1.d0 - rk4min_can_shv)
   rk4max_can_rvap = rk4max_can_shv / (1.d0 - rk4max_can_shv)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Minimum water mass at the leaf surface.  This is given in kg/m²leaf rather than   !
   ! kg/m²ground, so we scale it with LAI.                                                 !
   !---------------------------------------------------------------------------------------!
   rk4min_veg_lwater = -rk4dry_veg_lwater         ! Minimum leaf water mass     [kg/m²leaf]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    The minimum mass of surface water and virtual layer are given in m3/m3 rather than !
   ! kg/m2.  This is because there will be exchange between the top soil layer and the     !
   ! layers above in case the mass goes below the minimum.  Since this would make the im-  !
   ! pact of such exchange dependent on the soil depth, we assign the scale a function of  !
   ! the top layer thickness.                                                              !
   !---------------------------------------------------------------------------------------!
   rk4min_sfcw_moist = -5.0000d-4                  ! Minimum water mass allowed.
   rk4min_virt_moist = -5.0000d-4                  ! Minimum water allowed at virtual pool.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These variables are assigned in ed_params.f90.  Heat area should be 2.0 for all   !
   ! PFTs (two sides of the leaves exchange heat), and the evaporation area should be 1.0  !
   ! for all PFTs (only one side of the leaf is usually covered by water).  The transpir-  !
   ! ation area should be 1.0 for hypostomatous leaves, and 2.0 for symmetrical (pines)    !
   ! and amphistomatous (araucarias) leaves.  Sometimes heat and evaporation are multi-    !
   ! plied  by 1.2 and 2.2 to account for branches and twigs.  This is not recommended,    !
   ! though, because branches and twigs do not contribute to heat storage when             !
   ! ibranch_thermo is set to zero, and they are otherwise accounted through the wood area !
   ! index.                                                                                !
   !---------------------------------------------------------------------------------------!
   effarea_heat   = 2.d0 ! Heat area: related to 2*LAI
   effarea_evap   = 1.d0 ! Evaporation area: related to LAI
   !----- Transpiration.  Adjust them to 1 or 2 according to the type of leaf. ------------!
   effarea_transp(1:5)  = 1.d0 ! Hypostomatous
   effarea_transp(6:8)  = 2.d0 ! Symmetrical
   effarea_transp(9:16) = 1.d0 ! Hypostomatous
   effarea_transp(17)   = 2.d0 ! Amphistomatous
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag is used to control evaporation and transpiration when the air is         !
   ! saturated or super-saturated.  If supersat_ok is TRUE, then evaporation and           !
   ! transpiration will continue to happen even if the air is super-saturated at the       !
   ! canopy air temperature (but not at the soil or vegetation temperature).  Otherwise,   !
   ! evaporation and transpiration will be interrupted until the air becomes sub-saturated !
   ! again.  The air can still become super-saturated because mixing with the free atmo-   !
   ! sphere will not stop, but that is unlikely.                                           !
   !---------------------------------------------------------------------------------------!
   supersat_ok = .false.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The integrator will adjust pressure every time step, including the internal ones, !
   ! to make sure the ideal gas is respected.  If set to false, it will keep pressure      !
   ! constant within on DTLSM time step, and not bother forcing the canopy air space to    !
   ! respect the ideal gas equation.                                                       !
   !---------------------------------------------------------------------------------------!
   force_idealgas = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This flag is to turn on and on the leaf interception.  Except for developer       !
   ! tests, this variable should be always true.                                           !
   !---------------------------------------------------------------------------------------!
   leaf_intercept = .true.
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_rk4_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine overwrite_with_xml_config(thisnode)
   !!! PARSE XML FILE
   use ed_max_dims,   only: n_pft
   use ed_misc_coms,  only: iedcnfgf
   use hydrology_coms,only: useTOPMODEL
   use soil_coms,     only: isoilbc
   implicit none
   integer, intent(in) :: thisnode
   integer             :: max_pft_xml
   logical             :: iamhere

   if (iedcnfgf /= '') then
      inquire (file=iedcnfgf,exist=iamhere)
      if (iamhere) then

         !! FIRST, determine number of pft's defined in xml file
         call count_pft_xml_config(trim(iedcnfgf),max_pft_xml)
         if(max_pft_xml > n_pft) then

            write(unit=*,fmt='(a)') '*********************************************'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**  Number of PFTs required by XML Config  **'
            write(unit=*,fmt='(a)') '**  exceeds the memory available           **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**  Please change n_pft in Module ed_max_dims **'
            write(unit=*,fmt='(a)') '**  and recompile                          **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '*********************************************'
            call fatal_error('Too many PFTs','overwrite_with_xml_config','ed_params.f90')
         end if

         !! SECOND, update parameter defaults from XML
         call read_ed_xml_config(trim(iedcnfgf))

         !! THIRD, recalculate any derived values based on xml
         if(useTOPMODEL == 1) then
            isoilbc = 0
            print*,"TOPMODEL enabled, setting ISOILBC to 0"
         end if


         !! FINALLY, write out copy of settings
         call write_ed_xml_config()
   !      stop
      elseif (thisnode == 1) then
            write(unit=*,fmt='(a)') '*********************************************'
            write(unit=*,fmt='(a)') '**               WARNING!                  **'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '**    XML file wasn''t found. Using default **'
            write(unit=*,fmt='(a)') '** parameters in ED.                       **'
            write(unit=*,fmt='(a)') '** (You provided '//trim(iedcnfgf)//').'
            write(unit=*,fmt='(a)') '**                                         **'
            write(unit=*,fmt='(a)') '*********************************************'
      end if
   end if  !! end XML
   return
end subroutine overwrite_with_xml_config
!==========================================================================================!
!==========================================================================================!
