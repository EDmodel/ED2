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
                          , grass_pft           ! ! intent(out)
   use disturb_coms, only : ianth_disturb       ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer :: p
   !---------------------------------------------------------------------------------------!

   !----- Loading several parameters ------------------------------------------------------!
   call init_decomp_params()
   call init_ff_coms()
   call init_disturb_params()
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
   !  PFT | Name                                | Grass   | Tropical | agriculture?        !
   !------+-------------------------------------+---------+----------+---------------------!
   !    1 | C4 grass                            |     yes |      yes |                 yes !
   !    2 | Early tropical                      |      no |      yes |                  no !
   !    3 | Mid tropical                        |      no |      yes |                  no !
   !    4 | Late tropical                       |      no |      yes |                  no !
   !    5 | C3 grass                            |     yes |       no |                 yes !
   !    6 | Northern pines                      |      no |       no |                  no !
   !    7 | Southern pines                      |      no |       no |                  no !
   !    8 | Late conifers                       |      no |       no |                  no !
   !    9 | Early temperate deciduous           |      no |       no |                  no !
   !   10 | Mid temperate deciduous             |      no |       no |                  no !
   !   11 | Late temperate deciduous            |      no |       no |                  no !
   !   12 | C3 pasture                          |     yes |       no |                 yes !
   !   13 | C3 crop (e.g.,wheat, rice, soybean) |     yes |       no |                 yes !
   !   14 | C4 pasture                          |     yes |      yes |                 yes !
   !   15 | C4 crop (e.g.,corn/maize)           |     yes |      yes |                 yes !
   !------+-------------------------------------+---------+----------+---------------------!

   !----- Defining the grass PFTs ---------------------------------------------------------!
   grass_pft=huge(1)
   grass_pft(1)=1
   grass_pft(2)=5
   grass_pft(3)=12
   grass_pft(4)=13
   grass_pft(5)=14
   grass_pft(6)=15

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
   call init_pft_mort_params()
   call init_pft_alloc_params()
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
   use ed_misc_coms , only : burnin              & ! intent(out)
                           , outputMonth         & ! intent(out)
                           , restart_target_year & ! intent(out)
                           , use_target_year     & ! intent(out)
                           , maxage              & ! intent(out)
                           , dagei               & ! intent(out)
                           , maxdbh              & ! intent(out)
                           , ddbhi               & ! intent(out)
                           , vary_elev           &
                           , vary_hyd            &
                           , vary_rad
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
   atm_shv_max = 3.e-2    ! That corresponds to a dew point of 32°C at 1000hPa.
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
   lapse%atm_enthalpy = 0.0
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
                                    , lai_min                     & ! intent(out)
                                    , tai_min                     & ! intent(out)
                                    , patch_lai_min               & ! intent(out)!!!SOLVE2
                                    , blfac_min                   ! ! intent(out)
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

   diffuse_backscatter_vis(1:4)   = diffuse_bscat_vis_trop
   diffuse_backscatter_vis(5:11)  = diffuse_bscat_vis_temp
   diffuse_backscatter_vis(12:13) = diffuse_bscat_vis_temp
   diffuse_backscatter_vis(14:15) = diffuse_bscat_vis_trop

   emis_v(1)     = 9.60d-1
   emis_v(2:4)   = 9.50d-1
   emis_v(5)     = 9.60d-1
   emis_v(6:8)   = 9.70d-1
   emis_v(9:11)  = 9.50d-1
   emis_v(12:15) = 9.60d-1

   patch_lai_min = 1.0e-4
   lai_min       = 1.0e-5
   tai_min       = 1.0e-5
   blfac_min     = 1.0e-2

   return
end subroutine init_can_rad_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some canopy air related parameters.                       !
!------------------------------------------------------------------------------------------!
subroutine init_can_air_params()
   use consts_coms    , only : twothirds             & ! intent(in)
                             , vonk                  ! ! intent(in)
   use pft_coms       , only : hgt_min               ! ! intent(in)
   use canopy_air_coms, only : icanturb              & ! intent(in)
                             , isfclyrm              & ! intent(in)
                             , dry_veg_lwater        & ! intent(out)
                             , fullveg_lwater        & ! intent(out)
                             , rb_inter              & ! intent(out)
                             , rb_slope              & ! intent(out)
                             , veg_height_min        & ! intent(out)
                             , minimum_canopy_depth  & ! intent(out)
                             , minimum_canopy_depth8 & ! intent(out)
                             , exar                  & ! intent(out)
                             , covr                  & ! intent(out)
                             , ustmin                & ! intent(out)
                             , ubmin                 & ! intent(out)
                             , exar8                 & ! intent(out)
                             , ez                    & ! intent(out)
                             , vh2vr                 & ! intent(out)
                             , vh2dh                 & ! intent(out)
                             , ustmin8               & ! intent(out)
                             , ubmin8                & ! intent(out)
                             , ez8                   & ! intent(out)
                             , vh2dh8                & ! intent(out)
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
                             , fcod                  & ! intent(out)
                             , etf                   & ! intent(out)
                             , z0moz0h               & ! intent(out)
                             , z0hoz0m               & ! intent(out)
                             , ribmaxbh91            & ! intent(out)
                             , gamm                  & ! intent(out)
                             , gamh                  & ! intent(out)
                             , tprandtl              & ! intent(out)
                             , vkopr                 & ! intent(out)
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
                             , fcod8                 & ! intent(out)
                             , etf8                  & ! intent(out)
                             , z0moz0h8              & ! intent(out)
                             , z0hoz0m8              ! ! intent(out)
   implicit none
   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 plant], so it will be always scaled by      !
   ! (LAI+WAI).                                                                            !
   !---------------------------------------------------------------------------------------!
   dry_veg_lwater = 5.e-4
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Maximum leaf water that plants can hold.  Should leaf water exceed this number,    !
   ! water will be no longer intercepted by the leaves, and any value in excess of this    !
   ! will be promptly removed through shedding.  This value is in kg/[m2 plant], so it     !
   ! will be always scaled by (LAI+WAI).                                                   !
   !---------------------------------------------------------------------------------------!
   fullveg_lwater = 0.11
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Variables to define the vegetation aerodynamic resistance.  They are currently   !
   ! not PFT dependent.                                                                    !
   !---------------------------------------------------------------------------------------!
   rb_slope = 25.0
   rb_inter = 0.0 

   !---------------------------------------------------------------------------------------!
   ! veg_height_min       - This is the minimum vegetation height allowed [m].  Vegetation !
   !                        height is used to calculate drag coefficients and patch        !
   !                        roughness.                                                     !
   ! minimum_canopy_depth - This is the minimum canopy depth allowed [m].  Canopy depth    !
   !                        is used to calculate the heat and moisture storage capacity in !
   !                        the canopy air space.                                          !
   !---------------------------------------------------------------------------------------!
   select case (icanturb)
   case (-1)
      veg_height_min        = minval(hgt_min) ! used to be 0.2
      minimum_canopy_depth  = minval(hgt_min) ! used to be 0.2

   case default
      veg_height_min        = 1.0             ! alternative: minval(hgt_min) 
      minimum_canopy_depth  = 5.0             ! alternative: minval(hgt_min) 
   end select

   !----- This is the dimensionless exponential wind atenuation factor. -------------------!
   exar  = 2.5

   !----- This is the scaling factor of tree area index (not sure if it is used...) -------!
   covr = 2.16


   !---------------------------------------------------------------------------------------!
   !      Parameters for surface layer models.                                             !
   !     For the time being we keep the old values of ustmin and ubmin when the user wants !
   ! to run using Louis (1979).  This should be tested because often the value of u* is    !
   ! forced to 0.10m/s, which means that the model is not being applied.                   !
   !---------------------------------------------------------------------------------------!
   select case (isfclyrm)
   case (1)
      !----- This is the minimum ustar under stable and unstable conditions. --------------!
      ustmin    = 0.10
      !----- This is the minimum wind scale under stable and unstable conditions. ---------!
      ubmin     = 0.65
   case default
      !----- This is the minimum ustar under stable and unstable conditions. --------------!
      ustmin    = 0.025
      !----- This is the minimum wind scale under stable and unstable conditions. ---------!
      ubmin     = 0.25
   end select
   !----- Louis (1979) model. -------------------------------------------------------------!
   bl79        = 5.0    ! b prime parameter
   csm         = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   csh         = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   dl79        = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   bbeta       = 5.0           ! Beta 
   ribmaxod95  = 0.20          ! Maximum bulk Richardson number
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   abh91       = -0.70         ! -a from equation  (28) 
   bbh91       = -0.75         ! -b from equation  (28)
   cbh91       =  5.0          !  c from equations (28) and (32)
   dbh91       =  0.35         !  d from equations (28) and (32)
   ebh91       = -1.00         ! -a from equation  (32)
   fbh91       = -twothirds    ! -b from equation  (32)
   cod         = cbh91/dbh91   ! c/d
   bcod        = bbh91 * cod   ! b*c/d
   fcod        = fbh91 * cod   ! f*c/d
   etf         = ebh91 * fbh91 ! e * f
   z0moz0h     = 1.0           ! z0(M)/z0(h)
   z0hoz0m     = 1. / z0moz0h  ! z0(M)/z0(h)
   ribmaxbh91  = 6.00          ! Maximum bulk Richardson number
   !----- Used by OD95 and BH91. ----------------------------------------------------------!
   gamm        = 16.0          ! Gamma for momentum.
   gamh        = 16.0          ! Gamma for heat.
   tprandtl    = 1.00          ! Turbulent Prandtl number.
   vkopr       = vonk/tprandtl ! Von Karman / Prandtl number
   !---------------------------------------------------------------------------------------!


   
   !----- This is the relation between displacement height and roughness when icanturb=-1. !
   ez  = 0.172

   !----- This is the conversion from veg. height to roughness when icanturb /= -1. -------!
   vh2vr = 0.13
   
   !----- This is the conversion from vegetation height to displacement height. -----------!
   vh2dh  = 0.63


   !----- Set the double precision variables. ---------------------------------------------!
   minimum_canopy_depth8 = dble(minimum_canopy_depth)
   exar8                 = dble(exar                )
   ubmin8                = dble(ubmin               )
   ustmin8               = dble(ustmin              )
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
   fcod8                 = dble(fcod                )
   etf8                  = dble(etf                 )
   z0moz0h8              = dble(z0moz0h             )
   z0hoz0m8              = dble(z0hoz0m             )

   return
end subroutine init_can_air_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_photo_params()

   use ed_max_dims, only : n_pft                ! ! intent(in)
   use pft_coms   , only : D0                   & ! intent(out)
                         , Vm_low_temp          & ! intent(out)
                         , Vm_high_temp         & ! intent(out)
                         , Vm0                  & ! intent(out)
                         , stomatal_slope       & ! intent(out)
                         , cuticular_cond       & ! intent(out)
                         , quantum_efficiency   & ! intent(out)
                         , photosyn_pathway     ! ! intent(out)

   implicit none


   D0(1:15)                  = 0.01 ! same for all PFTs

   Vm_low_temp(1)            = 5.0     ! c4 grass
   Vm_low_temp(2)            = 5.0     ! early tropical
   Vm_low_temp(3)            = 5.0     ! mid tropical
   Vm_low_temp(4)            = 5.0     ! late tropical
   Vm_low_temp(5)            = 4.7137  ! c3 grass
   Vm_low_temp(6)            = 4.7137  ! northern pines ! 5.0
   Vm_low_temp(7)            = 4.7137  ! southern pines ! 5.0
   Vm_low_temp(8)            = 4.7137  ! late conifers  ! 5.0
   Vm_low_temp(9)            = 4.7137  ! early hardwoods
   Vm_low_temp(10)           = 4.7137  ! mid hardwoods
   Vm_low_temp(11)           = 4.7137  ! late hardwoods
   Vm_low_temp(12)           = 4.7137  ! c3 pasture
   Vm_low_temp(13)           = 4.7137  ! c3 crop
   Vm_low_temp(14)           = 5.0     ! c4 pasture
   Vm_low_temp(15)           = 5.0     ! c4 crop

   Vm_high_temp(1)           = 100.0    ! C4
   Vm_high_temp(2)           =  45.0    ! C3
   Vm_high_temp(3)           =  45.0    ! C3
   Vm_high_temp(4)           =  45.0    ! C3
   Vm_high_temp(5)           =  45.0    ! C3
   Vm_high_temp(6)           =  45.0    ! C3
   Vm_high_temp(7)           =  45.0    ! C3
   Vm_high_temp(8)           =  45.0    ! C3
   Vm_high_temp(9)           =  45.0    ! C3
   Vm_high_temp(10)          =  45.0    ! C3
   Vm_high_temp(11)          =  45.0    ! C3
   Vm_high_temp(12)          =  45.0    ! C3
   Vm_high_temp(13)          =  45.0    ! C3
   Vm_high_temp(14)          = 100.0    ! C4
   Vm_high_temp(15)          = 100.0    ! C4

   Vm0(1)                    = 12.5
   Vm0(2)                    = 18.8
   Vm0(3)                    = 12.5
   Vm0(4)                    = 6.25
   Vm0(5)                    = 18.3
   Vm0(6)                    = 15.625 * 0.7264
   Vm0(7)                    = 15.625 * 0.7264
   Vm0(8)                    = 6.25   * 0.7264
   Vm0(9)                    = 18.25  * 1.1171
   Vm0(10)                   = 15.625 * 1.1171
   Vm0(11)                   = 6.25   * 1.1171
   Vm0(12:13)                = 18.3
   Vm0(14:15)                = 12.5

   stomatal_slope(1)         = 10.0
   stomatal_slope(2)         = 8.0
   stomatal_slope(3)         = 8.0
   stomatal_slope(4)         = 8.0
   stomatal_slope(5)         = 6.3949
   stomatal_slope(6)         = 6.3949 ! 8.0
   stomatal_slope(7)         = 6.3949 ! 8.0
   stomatal_slope(8)         = 6.3949 ! 8.0
   stomatal_slope(9)         = 6.3949
   stomatal_slope(10)        = 6.3949
   stomatal_slope(11)        = 6.3949
   stomatal_slope(12)        = 6.3949
   stomatal_slope(13)        = 6.3949
   stomatal_slope(14)        = 10.0 
   stomatal_slope(15)        = 10.0 

   cuticular_cond(1)         = 10000.0
   cuticular_cond(2)         = 10000.0
   cuticular_cond(3)         = 10000.0
   cuticular_cond(4)         = 10000.0
   cuticular_cond(5)         = 10000.0
   cuticular_cond(6)         = 1000.0
   cuticular_cond(7)         = 1000.0
   cuticular_cond(8)         = 1000.0
   cuticular_cond(9)         = 20000.0
   cuticular_cond(10)        = 20000.0
   cuticular_cond(11)        = 20000.0
   cuticular_cond(12)        = 10000.0
   cuticular_cond(13)        = 10000.0
   cuticular_cond(14)        = 10000.0
   cuticular_cond(15)        = 10000.0

   quantum_efficiency(1)     = 0.06
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
   quantum_efficiency(14)    = 0.06
   quantum_efficiency(15)    = 0.06

   photosyn_pathway(1)       = 4
   photosyn_pathway(2:4)     = 3
   photosyn_pathway(5)       = 3
   photosyn_pathway(6:13)    = 3
   photosyn_pathway(14:15)   = 4

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

   use pft_coms   , only : growth_resp_factor        & ! intent(out)
                         , leaf_turnover_rate        & ! intent(out)
                         , root_turnover_rate        & ! intent(out)
                         , dark_respiration_factor   & ! intent(out)
                         , storage_turnover_rate     & ! intent(out)
                         , root_respiration_factor   ! ! intent(out)
   use decomp_coms, only : f_labile                  ! ! intent(out)

   implicit none

   growth_resp_factor(1)          = 0.333
   growth_resp_factor(2)          = 0.333
   growth_resp_factor(3)          = 0.333
   growth_resp_factor(4)          = 0.333
   growth_resp_factor(5)          = 0.333
   growth_resp_factor(6)          = 0.4503 ! 0.333
   growth_resp_factor(7)          = 0.4503 ! 0.333
   growth_resp_factor(8)          = 0.4503 ! 0.333
   growth_resp_factor(9)          = 0.0
   growth_resp_factor(10)         = 0.0
   growth_resp_factor(11)         = 0.0
   growth_resp_factor(12)         = 0.333
   growth_resp_factor(13)         = 0.333
   growth_resp_factor(14)         = 0.333
   growth_resp_factor(15)         = 0.333

   leaf_turnover_rate(1)          = 2.0
   leaf_turnover_rate(2)          = 1.0
   leaf_turnover_rate(3)          = 0.5
   leaf_turnover_rate(4)          = 0.333
   leaf_turnover_rate(5)          = 2.0
   leaf_turnover_rate(6)          = 0.333
   leaf_turnover_rate(7)          = 0.333
   leaf_turnover_rate(8)          = 0.333
   leaf_turnover_rate(9)          = 0.0
   leaf_turnover_rate(10)         = 0.0
   leaf_turnover_rate(11)         = 0.0
   leaf_turnover_rate(12)         = 2.0
   leaf_turnover_rate(13)         = 2.0
   leaf_turnover_rate(14)         = 2.0
   leaf_turnover_rate(15)         = 2.0

   !----- Root turnover rate.  ------------------------------------------------------------!
   root_turnover_rate(1)          = 2.0
   root_turnover_rate(2)          = 1.0
   root_turnover_rate(3)          = 0.5
   root_turnover_rate(4)          = 0.333
   root_turnover_rate(5)          = 0.333
   root_turnover_rate(6)          = 3.927218 ! 0.333
   root_turnover_rate(7)          = 4.117847 ! 0.333
   root_turnover_rate(8)          = 3.800132 ! 0.333
   root_turnover_rate(9)          = 5.772506
   root_turnover_rate(10)         = 5.083700
   root_turnover_rate(11)         = 5.070992
   root_turnover_rate(12)         = 0.333
   root_turnover_rate(13)         = 0.333
   root_turnover_rate(14)         = 2.0
   root_turnover_rate(15)         = 2.0

   dark_respiration_factor(1)     = 0.04
   dark_respiration_factor(2)     = 0.02
   dark_respiration_factor(3)     = 0.02
   dark_respiration_factor(4)     = 0.02
   dark_respiration_factor(5)     = 0.04
   dark_respiration_factor(6)     = 0.02
   dark_respiration_factor(7)     = 0.02
   dark_respiration_factor(8)     = 0.02
   dark_respiration_factor(9)     = 0.02
   dark_respiration_factor(10)    = 0.02
   dark_respiration_factor(11)    = 0.02
   dark_respiration_factor(12)    = 0.04
   dark_respiration_factor(13)    = 0.04
   dark_respiration_factor(14)    = 0.04
   dark_respiration_factor(15)    = 0.04

   storage_turnover_rate(1)       = 0.0
   storage_turnover_rate(2)       = 0.0
   storage_turnover_rate(3)       = 0.0
   storage_turnover_rate(4)       = 0.0
   storage_turnover_rate(5)       = 0.0
   storage_turnover_rate(6)       = 0.0
   storage_turnover_rate(7)       = 0.0
   storage_turnover_rate(8)       = 0.0
   storage_turnover_rate(9)       = 0.6243
   storage_turnover_rate(10)      = 0.6243
   storage_turnover_rate(11)      = 0.6243
   storage_turnover_rate(12)      = 0.0
   storage_turnover_rate(13)      = 0.0
   storage_turnover_rate(14)      = 0.0
   storage_turnover_rate(15)      = 0.0

   root_respiration_factor        = 0.528

   f_labile(1:5)                  = 1.0
   f_labile(6:11)                 = 0.79
   f_labile(12:15)                = 1.0

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
                          , seedling_mortality         & ! intent(out)
                          , treefall_s_gtht            & ! intent(out)
                          , treefall_s_ltht            & ! intent(out)
                          , plant_min_temp             & ! intent(out)
                          , frost_mort                 ! ! intent(out)
   use consts_coms , only : t00                        ! ! intent(in)
   use disturb_coms, only : treefall_disturbance_rate  ! ! intent(inout)

   implicit none


   frost_mort(1)     = 3.0
   frost_mort(2:4)   = 3.0
   frost_mort(5)     = 3.0
   frost_mort(6:11)  = 3.0
   frost_mort(12:13) = 3.0
   frost_mort(14:15) = 3.0


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

   mort3(1)  =  0.06167    ! 0.037
   mort3(2)  =  0.06167    ! 0.037
   mort3(3)  =  0.03167    ! 0.019
   mort3(4)  =  0.0
   mort3(5)  =  0.066
   mort3(6)  =  0.0033928
   mort3(7)  =  0.0043
   mort3(8)  =  0.0023568
   mort3(9)  =  0.006144
   mort3(10) =  0.003808
   mort3(11) =  0.00428
   mort3(12) =  0.066
   mort3(13) =  0.066
   mort3(14) =  0.037
   mort3(15) =  0.037
   
   if (treefall_disturbance_rate < 0.) then
      mort3(:) = mort3(:) - treefall_disturbance_rate
      treefall_disturbance_rate = 0.
   end if
   
   seedling_mortality = 0.95

   treefall_s_gtht = 0.0

   treefall_s_ltht(1)     = 0.25
   treefall_s_ltht(2:4)   = 0.1
   treefall_s_ltht(5)     = 0.25
   treefall_s_ltht(6:11)  = 0.1
   treefall_s_ltht(12:15) = 0.25

   plant_min_temp(1:4)   = t00
   plant_min_temp(5:6)   = t00-80.0
   plant_min_temp(7)     = t00-10.0
   plant_min_temp(8)     = t00-60.0
   plant_min_temp(9)     = t00-80.0
   plant_min_temp(10:11) = t00-20.0
   plant_min_temp(12:13) = t00-80.0
   plant_min_temp(14:15) = t00

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
                          , hgt_min               & ! intent(out)
                          , b1Ht                  & ! intent(out)
                          , b2Ht                  & ! intent(out)
                          , b1Bs                  & ! intent(out)
                          , b2Bs                  & ! intent(out)
                          , b1Bl                  & ! intent(out)
                          , b2Bl                  & ! intent(out)
                          , C2B                   & ! intent(out)
                          , hgt_ref               & ! intent(out)
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
   use consts_coms , only : twothirds             ! ! intent(in)
   implicit none

   !----- Carbon-to-biomass ratio of plant tissues. ---------------------------------------!
   C2B    = 2.0

   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is tropical or not.           !
   !---------------------------------------------------------------------------------------! 
   is_tropical(1:4)   = .true.
   is_tropical(5:11)  = .false.
   is_tropical(12:13) = .false.
   is_tropical(14:15) = .true.

   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is tree or grass              !
   !---------------------------------------------------------------------------------------! 
   is_grass(1)     = .true.
   is_grass(2:4)   = .false.
   is_grass(5)     = .true.
   is_grass(6:11)  = .false.
   is_grass(12:15) = .true.

   !---------------------------------------------------------------------------------------!
   !     Wood density.  Currently only tropical PFTs need it.  C3 grass density will be    !
   ! used only for branch area purposes.                                                   !
   !---------------------------------------------------------------------------------------!
![KIM] - new tropical parameters
!   rho(1)     = 0.53
!   rho(2)     = 0.53
!   rho(3)     = 0.71
!   rho(4)     = 0.90
   rho(1)     = 0.40
   rho(2)     = 0.40
   rho(3)     = 0.60
   rho(4)     = 0.87
   rho(5)     = 0.53   ! Copied from C4 grass
   rho(6:11)  = 0.00   ! Currently not used
   rho(12:13) = 0.53
   rho(14:15) = 0.53
   !---------------------------------------------------------------------------------------!

   !----- Specific leaf area [m² leaf / kg C] ---------------------------------------------!
![KIM] - new tropical parameters
!   SLA(1:4)   = 10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(1:4)))) * C2B * 0.1
   SLA(1:4) = 10.0**(1.6923-0.3305*log10(12.0/leaf_turnover_rate(1:4)))
   SLA(5)     = 22.0
   SLA(6)     =  6.0
   SLA(7)     =  9.0
   SLA(8)     = 10.0
   SLA(9)     = 30.0
   SLA(10)    = 24.2
   SLA(11)    = 60.0
   SLA(12:13) = 22.0
!   SLA(14:15) = 10.0**((2.4-0.46*log10(12.0/leaf_turnover_rate(14:15)))) * C2B * 0.1
   SLA(14:15) = 10.0**(1.6923-0.3305*log10(12.0/leaf_turnover_rate(14:15)))

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

   sapwood_ratio(1:15) = 3900.0

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


   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Initial density of plants, for near-bare-ground simulations [# of individuals/m2]  !
   !---------------------------------------------------------------------------------------!
   init_density(1)     = 0.6
   init_density(2:4)   = 0.1
   init_density(5)     = 0.6
   init_density(6:8)   = 0.1
   init_density(9:11)  = 0.1
   init_density(12:13) = 0.1
   init_density(14:15) = 0.6 

   !---------------------------------------------------------------------------------------!
   !    Minimum height of an individual.                                                   !
   !---------------------------------------------------------------------------------------!
   hgt_min(1)     = 0.50
   hgt_min(2:4)   = 0.50
   hgt_min(5)     = 0.15
   hgt_min(6:7)   = 1.50
   hgt_min(8)     = 1.50
   hgt_min(9)     = 1.50
   hgt_min(10)    = 1.50
   hgt_min(11)    = 1.50
   hgt_min(12:13) = 0.15
   hgt_min(14:15) = 0.50

   !----- Reference height for diameter/height allometry (temperates only). ---------------!
   hgt_ref(1:5)   = 0.0
   hgt_ref(6:11)  = 1.3
   hgt_ref(12:15) = 0.0

   !----- Fraction of structural stem that is assumed to be above ground. -----------------!
   agf_bs = 0.7

   !---------------------------------------------------------------------------------------!
   !    DBH/height allometry parameters.  They are used only for temperate PFTs.           !
   !---------------------------------------------------------------------------------------!
   !----- DBH-height allometry intercept [m]. ---------------------------------------------!
   b1Ht(1:4)   = 0.0
   b1Ht(5)     = 0.4778
   b1Ht(6)     = 27.14
   b1Ht(7)     = 27.14
   b1Ht(8)     = 22.79
   b1Ht(9)     = 22.6799
   b1Ht(10)    = 25.18
   b1Ht(11)    = 23.3874
   b1Ht(12:13) = 0.4778
   b1Ht(14:15) = 0.0
   !----- DBH-height allometry slope [1/cm]. ----------------------------------------------!
   b2Ht(1:4)   = 0.0
   b2Ht(5)     = -0.75
   b2Ht(6)     = -0.03884
   b2Ht(7)     = -0.03884
   b2Ht(8)     = -0.04445 
   b2Ht(9)     = -0.06534
   b2Ht(10)    = -0.04964
   b2Ht(11)    = -0.05404
   b2Ht(12:13) = -0.75
   b2Ht(14:15) = 0.0
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
   !----- DBH-stem allometry intercept [kg stem biomass / plant * cm^(-b2Bs)] -------------!
   b1Bs(1:4)   = 0.0 
   b1Bs(5)     = 1.0e-5
   b1Bs(6)     = 0.147
   b1Bs(7)     = 0.147
   b1Bs(8)     = 0.1617
   b1Bs(9)     = 0.02648
   b1Bs(10)    = 0.1617
   b1Bs(11)    = 0.235
   b1Bs(12:13) = 1.0e-5
   b1Bs(14:15) = 0.0 
   !----- DBH-stem allometry slope [dimensionless]. ---------------------------------------!
   b2Bs(1:4)   = 0.0
   b2Bs(5)     = 1.0
   b2Bs(6)     = 2.238
   b2Bs(7)     = 2.238
   b2Bs(8)     = 2.1536
   b2Bs(9)     = 2.95954
   b2Bs(10)    = 2.4572
   b2Bs(11)    = 2.2518
   b2Bs(12:13) = 1.0
   b2Bs(14:15) = 0.0

   !---------------------------------------------------------------------------------------!
   !    Defining the branching parameters, following Järvelä (2004)                        !
   !---------------------------------------------------------------------------------------!
   !----- Branching ratio -----------------------------------------------------------------!
   rbranch(1)     = 4.24
   rbranch(2:4)   = 4.23
   rbranch(5)     = 4.24
   rbranch(6:8)   = 4.44
   rbranch(9:11)  = 4.24
   rbranch(12:15) = 4.24
   !----- Diameter ratio ------------------------------------------------------------------!
   rdiamet(1)     = 5.00
   rdiamet(2:4)   = 1.86
   rdiamet(5)     = 5.00
   rdiamet(6:8)   = 2.04
   rdiamet(9:11)  = 1.86
   rdiamet(12:15) = 5.00
   !----- Length ratio. Järvelä used rdiamet^2/3, so do we... -----------------------------!
   rlength(1:15)  = rdiamet(1:15)**twothirds
   !----- Minimum diameter to consider [cm]. ----------------------------------------------!
   diammin(1:15)  = 1.0
   !----- Number of trunks.  Usually this is 1. -------------------------------------------!
   ntrunk(1:15)   = 1.0
   
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
   conijn_b(1)     = 0.0
   conijn_b(2:4)   = -0.7178682
   conijn_b(5)     = 0.0
   conijn_b(6:11)  = -0.7178682
   conijn_b(12:15) = 0.0
   conijn_c(1)     = 0.0
   conijn_c(2:4)   = 0.00490734
   conijn_c(5)     = 0.0
   conijn_c(6:11)  = 0.00490734
   conijn_c(12:15) = 0.0
   conijn_d(1)     = 0.0
   conijn_d(2:4)   = -0.0456370
   conijn_d(5)     = 0.0
   conijn_d(6:11)  = -0.0456370
   conijn_d(12:15) = 0.0
   return
end subroutine init_pft_alloc_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_nitro_params()

use pft_coms, only: c2n_leaf, Vm0, SLA, water_conductance, &
     c2n_slow,c2n_structural,c2n_storage,c2n_stem,l2n_stem, &
     C2B,plant_N_supply_scale

implicit none

c2n_slow       = 10.0 ! Carbon to Nitrogen ratio, slow pool.
c2n_structural = 150.0 ! Carbon to Nitrogen ratio, structural pool.
c2n_storage    = 150.0 ! Carbon to Nitrogen ratio, storage pool.
c2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.
l2n_stem       = 150.0 ! Carbon to Nitrogen ratio, structural stem.


plant_N_supply_scale = 0.5 

c2n_leaf = 1000.0 / ((0.11289 + 0.12947 * Vm0) * SLA)
c2n_leaf(6) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(6))
c2n_leaf(7) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(7))
c2n_leaf(8) = 1000.0 / ((0.11289 + 0.12947 * 6.25) * SLA(8))
c2n_leaf(9) = 1000.0 / ((0.11289 + 0.12947 * 18.25) * SLA(9))
c2n_leaf(10) = 1000.0 / ((0.11289 + 0.12947 * 15.625) * SLA(10))
c2n_leaf(11) = 1000.0 / ((0.11289 + 0.12947 * 6.25) * SLA(11))



water_conductance(1:4) = 0.00476           !OLD=KW60=0.001904 KW80=0.00253505, KW100=0.00316881, KW125=0.00396101
water_conductance(5:11) = 0.00476           ! equal to KW150
water_conductance(12:13) = 0.00476
water_conductance(14:15) = 0.001904

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
                             , leaf_width           & ! intent(out)
                             , crown_depth_fraction & ! intent(out)
                             , c_grn_leaf_dry       & ! intent(out)
                             , c_ngrn_biom_dry      & ! intent(out)
                             , wat_dry_ratio_grn    & ! intent(out)
                             , wat_dry_ratio_ngrn   & ! intent(out)
                             , delta_c              ! ! intent(out)
   use consts_coms    , only : t3ple                ! ! intent(out) 
   use phenology_coms , only :iphen_scheme

   implicit none

   select case (iphen_scheme)
   case (0,1)
      phenology(1)     = 1
      phenology(2:4)   = 1
      phenology(5)     = 1
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 1
   case (2)
      phenology(1)     = 4
      phenology(2:4)   = 4
      phenology(5)     = 4
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 4
   case (3)
      phenology(1)     = 4
      phenology(2:4)   = 3
      phenology(5)     = 4
      phenology(6:8)   = 0
      phenology(9:11)  = 2
      phenology(12:15) = 4
   end select

   clumping_factor(1)     = 1.000d0
   clumping_factor(2:4)   = 7.350d-1
   clumping_factor(5)     = 8.400d-1
   clumping_factor(6:8)   = 7.350d-1
   clumping_factor(9:11)  = 8.400d-1
   clumping_factor(12:13) = 8.400d-1
   clumping_factor(14:15) = 1.000d0

   leaf_width(1:4)   = 0.20
   leaf_width(5:11)  = 0.05
   leaf_width(12:13) = 0.05
   leaf_width(14:15) = 0.20

   !---------------------------------------------------------------------------------------!
   !      The following parameters are second sources found in Gu et al. (2007)            !
   !---------------------------------------------------------------------------------------!
   c_grn_leaf_dry(1:15)      = 3218.0    ! Jones 1992  J/(kg K)
   c_ngrn_biom_dry(1:15)     = 1256.0    ! Forest Products Laboratory 
   wat_dry_ratio_grn(1:15)   = 2.5       ! 
   !wat_dry_ratio_grn(1:15)   = 1.5       ! Ceccato et al. 2001
   wat_dry_ratio_ngrn(1:15)  = 0.7       ! Forest Products Laboratory
   !---------------------------------------------------------------------------------------!
   !     Delta-c is found using the second term of the RHS of equation 5, assuming         !
   ! T=T3ple.  This is a simplification, but the specific heat usually varies by 3J/kg/K   !
   ! between 173K and 341K, so removing the dependence on temperature is not that bad      !
   ! assumption.                                                                           !
   !---------------------------------------------------------------------------------------!
   delta_c(1:15) = 100. * wat_dry_ratio_ngrn(1:15)                                         &
                 * (-0.06191 + 2.36e-4 * t3ple - 1.33e-2 * wat_dry_ratio_ngrn(1:15))

   !----- Relative height of crown. -------------------------------------------------------!
   crown_depth_fraction(1)     = 1.0    
   crown_depth_fraction(2:4)   = 0.25
   crown_depth_fraction(5)     = 1.0
   crown_depth_fraction(6:11)  = 0.35
   crown_depth_fraction(12:15) = 1.0

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
   r_fract(12:15)            = 1.0

   seed_rain(1:15)           = 0.01

   nonlocal_dispersal(1:5)   = 1.0
   nonlocal_dispersal(6:7)   = 0.766
   nonlocal_dispersal(8)     = 0.001
   nonlocal_dispersal(9)     = 1.0
   nonlocal_dispersal(10)    = 0.325
   nonlocal_dispersal(11)    = 0.074
   nonlocal_dispersal(12:15) = 1.0

   repro_min_h(1)            = 0.0
   repro_min_h(2:4)          = 5.0
   repro_min_h(5)            = 0.0
   repro_min_h(6:11)         = 5.0
   repro_min_h(12:15)        = 0.0

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
   use decomp_coms , only : f_labile             ! ! intent(in)
   use ed_max_dims , only : n_pft                ! ! intent(in)
   use consts_coms , only : onesixth             ! ! intent(in)
   use pft_coms    , only : init_density         & ! intent(in)
                          , c2n_leaf             & ! intent(in)
                          , c2n_stem             & ! intent(in)
                          , b1Ht                 & ! intent(in)
                          , b2Ht                 & ! intent(in)
                          , hgt_min              & ! intent(in)
                          , hgt_ref              & ! intent(in)
                          , q                    & ! intent(in)
                          , qsw                  & ! intent(in)
                          , max_dbh              & ! intent(out)
                          , min_recruit_size     & ! intent(out)
                          , min_cohort_size      & ! intent(out)
                          , negligible_nplant    & ! intent(out)
                          , c2n_recruit          ! ! intent(out)
   use allometry   , only : h2dbh                & ! function
                          , dbh2h                & ! function
                          , dbh2bl               & ! function
                          , dbh2bd               ! ! function
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer            :: ipft
   real               :: dbh
   real               :: huge_dbh
   real               :: huge_height
   real               :: balive_min
   real               :: bleaf_min
   real               :: bdead_min
   real               :: balive_max
   real               :: bleaf_max
   real               :: bdead_max
   real               :: min_plant_dens
   logical, parameter :: print_zero_table = .true.
   !---------------------------------------------------------------------------------------!

   !----- Maximum DBH. --------------------------------------------------------------------!
   max_dbh(1)     = 0.498
   max_dbh(2:4)   = 68.31
   max_dbh(5)     = 0.498
   max_dbh(6:11)  = log(1.0-(0.999*b1Ht(6:11)-hgt_ref(6:11))/b1Ht(6:11))/b2Ht(6:11)
   max_dbh(12:15) = 0.498


   !---------------------------------------------------------------------------------------!
   !     The minimum recruitment size and the recruit carbon to nitrogen ratio.  Both      !
   ! parameters actually depend on which PFT we are solving, since grasses always have     !
   ! significantly less biomass.                                                           !
   !---------------------------------------------------------------------------------------!
   if (print_zero_table) then
      write (unit=61,fmt='(13(a,1x))') '  PFT','     HGT_MIN','         DBH'               &
                                              ,'   BLEAF_MIN','   BDEAD_MIN'               &
                                              ,'  BALIVE_MIN','   BLEAF_MAX'               &
                                              ,'   BDEAD_MAX','  BALIVE_MAX'               &
                                              ,'   INIT_DENS','MIN_REC_SIZE'               &
                                              ,'MIN_COH_SIZE',' NEGL_NPLANT'
   end if
   min_plant_dens = onesixth * minval(init_density)
   do ipft = 1,n_pft

      !----- Find the DBH and carbon pools associated with a newly formed recruit. --------!
      dbh        = h2dbh(hgt_min(ipft),ipft)
      bleaf_min  = dbh2bl(dbh,ipft)
      bdead_min  = dbh2bd(dbh,hgt_min(ipft),ipft)
      balive_min = bleaf_min * (1.0 + q(ipft) + qsw(ipft) * hgt_min(ipft))

      !------------------------------------------------------------------------------------!
      !   Find the maximum bleaf and bdead supported.  This is to find the negligible      !
      ! nplant so we ensure that the cohort is always terminated if its mortality rate is  !
      ! very high.                                                                         !
      !------------------------------------------------------------------------------------!
      huge_dbh    = 3. * max_dbh(ipft)
      huge_height = dbh2h(ipft, max_dbh(ipft))
      bleaf_max   = dbh2bl(huge_dbh,ipft)
      bdead_max   = dbh2bd(huge_dbh,huge_height,ipft)
      balive_max  = bleaf_max * (1.0 + q(ipft) + qsw(ipft) * huge_height)

      !------------------------------------------------------------------------------------!
      !    The definition of the minimum recruitment size is the minimum amount of biomass !
      ! in kgC/m² is available for new recruits.  For the time being we use the near-bare  !
      ! ground state value as the minimum recruitment size, but this may change depending  !
      ! on how well it goes.                                                               !
      !------------------------------------------------------------------------------------!
      min_recruit_size(ipft) = min_plant_dens * (bdead_min + balive_min)
      
      !------------------------------------------------------------------------------------!
      !    Minimum size (measured as biomass of living and structural tissues) allowed in  !
      ! a cohort.  Cohorts with less biomass than this are going to be terminated.         !
      !------------------------------------------------------------------------------------! 
      min_cohort_size(ipft)  = 0.1 * min_recruit_size(ipft)

      !------------------------------------------------------------------------------------! 
      !    The following variable is the absolute minimum cohort population that a cohort  !
      ! can have.  This should be used only to avoid nplant=0, but IMPORTANT: this will    !
      ! lead to a ridiculously small cohort almost guaranteed to be extinct and SHOULD BE  !
      ! USED ONLY IF THE AIM IS TO ELIMINATE THE COHORT.                                   !
      !------------------------------------------------------------------------------------! 
      negligible_nplant(ipft) = onesixth * min_cohort_size(ipft) / (bdead_max + balive_max)

      !----- Find the recruit carbon to nitrogen ratio. -----------------------------------!
      c2n_recruit(ipft)      = (balive_min + bdead_min)                                    &
                             / (balive_min * ( f_labile(ipft) / c2n_leaf(ipft)             &
                             + (1.0 - f_labile(ipft)) / c2n_stem(ipft))                    &
                             + bdead_min/c2n_stem(ipft))
      if (print_zero_table) then
         write (unit=61,fmt='(i5,1x,12(es12.5,1x))') ipft,hgt_min(ipft),dbh,bleaf_min      &
                                                    ,bdead_min,balive_min,bleaf_max        &
                                                    ,bdead_max,balive_max                  &
                                                    ,init_density(ipft)                    &
                                                    ,min_recruit_size(ipft)                &
                                                    ,min_cohort_size(ipft)                 &
                                                    ,negligible_nplant(ipft)
      end if
   end do

end subroutine init_pft_derived_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_disturb_params

   use disturb_coms , only : patch_dynamics           & ! intent(out)
                           , min_new_patch_area       & ! intent(out)
                           , treefall_hite_threshold  & ! intent(out)
                           , treefall_age_threshold   & ! intent(out)
                           , forestry_on              & ! intent(out)
                           , agriculture_on           & ! intent(out)
                           , plantation_year          & ! intent(out)
                           , plantation_rotation      & ! intent(out)
                           , mature_harvest_age       & ! intent(out)
                           , fire_dryness_threshold   & ! intent(out)
                           , fire_smoist_threshold    & ! intent(out)
                           , fire_smoist_depth        & ! intent(out)
                           , k_fire_first             & ! intent(out)
                           , fire_parameter           ! ! intent(out)

   implicit none
   
   patch_dynamics = 1
   
   min_new_patch_area = 0.005

   !----- Only trees above this height create a gap when they fall. -----------------------!
   treefall_hite_threshold = 10.0 

   !----- Minimum patch age for treefall disturbance. -------------------------------------!
   treefall_age_threshold = 0.0

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
   fire_smoist_threshold = 0.20

   !----- Maximum depth that will be considered in the average soil -----------------------!
   fire_smoist_depth     = -0.75

   !----- Dimensionless parameter controlling speed of fire spread. -----------------------!
   fire_parameter = 1.0
   
   return

end subroutine init_disturb_params
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
                             , soil_class            & ! type
                             , soil8                 & ! intent(out)
                             , water_stab_thresh     & ! intent(out)
                             , snowmin               & ! intent(out)
                             , dewmax                & ! intent(out)
                             , soil_rough            & ! intent(out)
                             , snow_rough            & ! intent(out)
                             , min_sfcwater_mass     & ! intent(out)
                             , infiltration_method   ! ! intent(out)
                             
   use grid_coms      , only : ngrids                ! ! intent(in)
   use consts_coms    , only : grav                  

   implicit none
   !----- Local variable ------------------------------------------------------------------!
   integer :: nsoil
   integer :: ifm
!   real    :: fieldcp_MPa        ! soil-water potential at field capacity         [ MPa]
   real    :: fieldcp_K          ! hydraulic conductivity at field capacity       [mm/day]
   real    :: soilcp_MPa         ! soil-water potential for air dry soil          [ MPa]
   real    :: soilwp_MPa         ! soil-water potential at wilting point          [ MPa]
   
!   fieldcp_MPa         = 0.033 
   fieldcp_K           = 0.1
   soilcp_MPa          = 3.1   
   soilwp_MPa          = 1.5

   water_stab_thresh   = 3.0    ! Minimum water mass to be considered stable     [   kg/m2]
   snowmin             = 3.0    ! Minimum snow mass needed to create a new layer [   kg/m2]
   dewmax              = 3.0e-5 ! Maximum dew flux rate (deprecated)             [ kg/m2/s]
   soil_rough          = 0.05   ! Soil roughness height                          [       m]
   snow_rough          = 0.001  ! Snowcover roughness height                     [       m]
   min_sfcwater_mass   = 1.0e-6 ! Minimum allowed mass in temporary layers       [   kg/m2]
   infiltration_method = 0      ! Infiltration method, used in rk4_derivs        [     0|1]


!------------------------------------------------------------------------------------------!
!      UPDATED based on Cosby et al. 1984, soil clay and sand fractions from table 2,      !
! slpots, slmsts, slbs, and slcons calculated based on table 4 sfldcap calculated based on !
! Clapp and Hornberger equation 2 and a soil hydraulic conductivity of 0.1mm/day, soilcp   !
! (air dry soil moisture capacity) is calculated based on Cosby et al. 1984 equation 1 and !
! a soil-water potential of -3.1MPa.  soilwp (the wilting point soil moisture) is          !
! calculated based on Cosby et al 1985 equation 1 and a soil-water potential of -1.5MPa    !
! (Equation 1 in Cosby is not correct it should be saturated moisture potential over       !
! moisture potential)  (NML, 2/2010)                                                       !
!------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
!                         slpots     slmsts       slbs     slcpd         soilcp        soilwp        slcons       slcons0 soilcond0 soilcond1 soilcond2       sfldcap    xsand    xclay   xorgan xrobulk    slden !
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
soil= (/ soil_class(-0.049831046,  0.373250,  3.295000,  1465000.,  0.026183447,  0.032636854,  2.446420e-5,  0.000500000,   0.3000,   4.8000,  -2.7000,  0.132130936,  0.9200,  0.0300,  0.0500,   1200.,  1600.)  & ! 1. Sand
      ,  soil_class(-0.068643592,  0.386340,  3.796000,  1407000.,  0.041873866,  0.050698650,  1.751180e-5,  0.000600000,   0.3000,   4.6600,  -2.6000,  0.155720881,  0.8200,  0.0600,  0.1200,   1250.,  1600.)  & ! 2. Loamy Sand
      ,  soil_class(-0.155095787,  0.418940,  4.496000,  1344000.,  0.076932647,  0.090413463,  8.228700e-6,  0.000769000,   0.2900,   4.2700,  -2.3100,  0.199963447,  0.5800,  0.1000,  0.3200,   1300.,  1600.)  & ! 3. Sandy loam
      ,  soil_class(-0.659933234,  0.476050,  5.090000,  1273000.,  0.141599910,  0.163306007,  2.396240e-6,  0.000010600,   0.2700,   3.4700,  -1.7400,  0.266720155,  0.1700,  0.1300,  0.7000,   1400.,  1600.)  & ! 4. Silt loam
      ,  soil_class(-0.238341682,  0.437280,  5.797000,  1214000.,  0.126501190,  0.143377074,  4.732940e-6,  0.002200000,   0.2800,   3.6300,  -1.8500,  0.247334654,  0.4300,  0.1800,  0.3900,   1350.,  1600.)  & ! 5. Loam
      ,  soil_class(-0.121199269,  0.412650,  7.165000,  1177000.,  0.137648733,  0.152325872,  6.405180e-6,  0.001500000,   0.2800,   3.7800,  -1.9600,  0.250954745,  0.5800,  0.2700,  0.1500,   1350.,  1600.)  & ! 6. Sandy clay loam
      ,  soil_class(-0.627769194,  0.478220,  8.408000,  1319000.,  0.228171947,  0.248747504,  1.435260e-6,  0.000107000,   0.2600,   2.7300,  -1.2000,  0.333825332,  0.1000,  0.3400,  0.5600,   1500.,  1600.)  & ! 7. Silty clay loam
      ,  soil_class(-0.281968114,  0.446980,  8.342000,  1227000.,  0.192624431,  0.210137962,  2.717260e-6,  0.002200000,   0.2700,   3.2300,  -1.5600,  0.301335491,  0.3200,  0.3400,  0.3400,   1450.,  1600.)  & ! 8. Clay loam
      ,  soil_class(-0.121283019,  0.415620,  9.538000,  1177000.,  0.182198910,  0.196607427,  4.314510e-6,  0.000002167,   0.2700,   3.3200,  -1.6300,  0.286363001,  0.5200,  0.4200,  0.0600,   1450.,  1600.)  & ! 9. Sandy clay
      ,  soil_class(-0.601312179,  0.479090,  10.46100,  1151000.,  0.263228486,  0.282143846,  1.055190e-6,  0.000001033,   0.2500,   2.5800,  -1.0900,  0.360319788,  0.0600,  0.4700,  0.4700,   1650.,  1600.)  & !10. Silty clay
      ,  soil_class(-0.286417797,  0.452300,  12.14000,  1088000.,  0.253968998,  0.269618857,  1.427350e-6,  0.000001283,   0.2500,   2.4000,  -0.9600,  0.348432364,  0.2200,  0.5800,  0.2000,   1700.,  1600.)  & !11. Clay
      ,  soil_class(-0.534564359,  0.469200,  6.180000,  8740000.,  0.167047523,  0.187868805,  2.357930e-6,  0.000008000,   0.0600,   0.4600,   0.0000,  0.285709966,  0.2000,  0.2000,  0.6000,    500.,   300.) /) !12. Peat

!------------------------------------------------------------------------------------------!
!********* Correct soil_class table using sand and clay fractions (if provided)  **********!
! Based on Cosby et al 1984, using table 4 and equation 1 (which is incorrect it should be !
! saturated moisture potential over moisture potential).  NML 2/2010
!--------------------------------------------------------------------------------------------------!
do ifm=1,ngrids
     if (isoilflg(ifm)==2 .and. slxclay >0. .and. slxsand > 0. .and. (slxclay + slxsand) <=1. ) then
        soil(nslcon)%slbs=        3.10 + 15.7*slxclay - 0.3*slxsand                                     ! [ unitless]
        
        soil(nslcon)%slpots=      -1 * (10**(2.17 - 0.63*slxclay - 1.58*slxsand)) * 0.01                ! [ m ]
        
        soil(nslcon)%slcons=      (10**(-0.60 + 1.26*slxsand - 0.64*slxclay)) * 0.0254/(60*60)          ! [ m/s ]
        
        soil(nslcon)%slmsts=      (50.5 - 14.2*slxsand - 3.7*slxclay)/100                               ! [ m^3/m^3 ]
        
!        soil(nslcon)%sfldcap=     soil(nslcon)%slmsts * ( soil(nslcon)%slpots / (fieldcp_MPa*1000/9.80665))** & 
!                                       (1 / soil(nslcon)%slbs)                                          ! [ m^3/m^3 ]
        soil(nslcon)%sfldcap=     soil(nslcon)%slmsts * ((fieldcp_K/1000/(60*60*24))/soil(nslcon)%slcons)** & 
                                       (1 / (2*soil(nslcon)%slbs+3))                                          ! [ m^3/m^3 ]
        
        soil(nslcon)%soilcp=      soil(nslcon)%slmsts * ( -1.*soil(nslcon)%slpots / (soilcp_MPa*1000/grav))** & 
                                       (1 / soil(nslcon)%slbs)                                          ! [ m^3/m^3 ]

        soil(nslcon)%soilwp=      soil(nslcon)%slmsts * ( -1.*soil(nslcon)%slpots / (soilwp_MPa*1000/grav))** & 
                                       (1 / soil(nslcon)%slbs)                                          ! [ m^3/m^3 ]       
                                       
        write(unit=*,fmt='(a)') '*********************************************'
        write(unit=*,fmt='(a)') '**RESET SOIL PARAMS:                         '
        write(unit=*,fmt='(a,1x,f8.4)') 'new slbs ', soil(nslcon)%slbs
        write(unit=*,fmt='(a,1x,f8.4)') 'new slpots ', soil(nslcon)%slpots
        write(unit=*,fmt='(a,1x,es9.2)') 'new slcons ', soil(nslcon)%slcons
        write(unit=*,fmt='(a,1x,f8.4)') 'new slmsts ', soil(nslcon)%slmsts
        write(unit=*,fmt='(a,1x,f8.4)') 'new sfldcap ', soil(nslcon)%sfldcap
        write(unit=*,fmt='(a,1x,f8.4)') 'new soilcp ', soil(nslcon)%soilcp
        write(unit=*,fmt='(a,1x,f8.4)') 'new soilwp ', soil(nslcon)%soilwp
        write(unit=*,fmt='(a)') '*********************************************'

     endif
enddo







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
      soil8(nsoil)%xorgan    = dble(soil(nsoil)%xorgan   )
      soil8(nsoil)%xrobulk   = dble(soil(nsoil)%xrobulk  )
      soil8(nsoil)%slden     = dble(soil(nsoil)%slden    )
   end do

   return
end subroutine init_soil_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_phen_coms

   use phenology_coms,only: retained_carbon_fraction, &
        theta_crit,dl_tr,st_tr1,st_tr2,phen_a,phen_b,phen_c, &
        rad_turnover_int, rad_turnover_slope, &
        vm_tran, vm_slop, vm_amp, vm_min

   implicit none

   retained_carbon_fraction = 0.5
   theta_crit = 0.2
   dl_tr = 655.0
   st_tr1 = 284.3
   st_tr2 = 275.15
  
   phen_a = -68.0   
   phen_b = 638.0    
   phen_c = -0.01    

   rad_turnover_int   = -11.3868
   rad_turnover_slope = 0.0824

   vm_tran = 9.0
   vm_slop = 10.0
   vm_amp = 20.0
   vm_min = 15.0

   return
end subroutine init_phen_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns the fusion and splitting parameters.                         !
!------------------------------------------------------------------------------------------!
subroutine init_ff_coms
   use fusion_fission_coms, only : min_dbh_class     & ! intent(out)
                                 , maxffdbh          & ! intent(out)
                                 , dffdbhi           & ! intent(out)
                                 , min_hgt_class     & ! intent(out)
                                 , fusetol           & ! intent(out)
                                 , fusetol_h         & ! intent(out)
                                 , lai_fuse_tol      & ! intent(out)
                                 , lai_tol           & ! intent(out)
                                 , ntol              & ! intent(out)
                                 , profile_tol       & ! intent(out)
                                 , max_patch_age     & ! intent(out)
                                 , ff_ndbh           & ! intent(out)
                                 , coh_tolerance_max & ! intent(out)
                                 , pat_tolerance_max & ! intent(out)
                                 , fuse_relax        ! ! intent(out)

   implicit none

   min_dbh_class     = 0.0  
   maxffdbh          = 200.0 
   min_hgt_class     = 0.0
   fusetol           = 0.4
   fusetol_h         = 0.5
   lai_fuse_tol      = 0.8
   lai_tol           = 1.0
   ntol              = 0.001
   profile_tol       = 0.2
   max_patch_age     = 500.0
   ff_ndbh           = 20
   coh_tolerance_max = 10.0 ! Original 2.0
   pat_tolerance_max = 100.0
   fuse_relax        = .false.
   dffdbhi           = real(ff_ndbh)/maxffdbh

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
   use soil_coms            , only : water_stab_thresh     & ! intent(in)
                                   , snowmin               & ! intent(in)
                                   , min_sfcwater_mass     ! ! intent(in)
   use canopy_air_coms      , only : dry_veg_lwater        & ! intent(in)
                                   , fullveg_lwater        ! ! intent(in)
   use rk4_coms             , only : rk4_tolerance         & ! intent(in)
                                   , maxstp                & ! intent(out)
                                   , rk4eps                & ! intent(out)
                                   , rk4eps2               & ! intent(out)
                                   , rk4epsi               & ! intent(out)
                                   , hmin                  & ! intent(out)
                                   , print_diags           & ! intent(out)
                                   , checkbudget           & ! intent(out)
                                   , newsnow               & ! intent(out)
                                   , debug                 & ! intent(out)
                                   , toocold               & ! intent(out)
                                   , toohot                & ! intent(out)
                                   , lai_to_cover          & ! intent(out)
                                   , hcapveg_ref           & ! intent(out)
                                   , min_height            & ! intent(out)
                                   , rk4min_veg_temp       & ! intent(out)
                                   , rk4water_stab_thresh  & ! intent(out)
                                   , rk4min_sfcwater_mass  & ! intent(out)
                                   , rk4dry_veg_lwater     & ! intent(out)
                                   , rk4fullveg_lwater     & ! intent(out)
                                   , rk4snowmin            & ! intent(out)
                                   , rk4min_can_temp       & ! intent(out)
                                   , rk4max_can_temp       & ! intent(out)
                                   , rk4min_can_shv        & ! intent(out)
                                   , rk4max_can_shv        & ! intent(out)
                                   , rk4max_can_rhv        & ! intent(out)
                                   , rk4min_can_co2        & ! intent(out)
                                   , rk4max_can_co2        & ! intent(out)
                                   , rk4min_soil_temp      & ! intent(out)
                                   , rk4max_soil_temp      & ! intent(out)
                                   , rk4min_veg_temp       & ! intent(out)
                                   , rk4max_veg_temp       & ! intent(out)
                                   , rk4min_veg_lwater     & ! intent(out)
                                   , rk4min_sfcw_temp      & ! intent(out)
                                   , rk4max_sfcw_temp      & ! intent(out)
                                   , rk4min_sfcw_moist     & ! intent(out)
                                   , rk4min_virt_moist     & ! intent(out)
                                   , check_maxleaf         & ! intent(out)
                                   , supersat_ok           ! ! intent(out)
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Copying some variables to the Runge-Kutta counterpart (double precision).         !
   !---------------------------------------------------------------------------------------!
   rk4water_stab_thresh = dble(water_stab_thresh)
   rk4min_sfcwater_mass = dble(min_sfcwater_mass)
   rk4dry_veg_lwater    = dble(dry_veg_lwater   )
   rk4fullveg_lwater    = dble(fullveg_lwater   )
   rk4snowmin           = dble(snowmin          )
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    The following variables control the Runge-Kutta performance.  Think twice before   !
   ! changing them...                                                                      !
   !---------------------------------------------------------------------------------------!
   maxstp      = 100000000           ! Maximum number of intermediate steps. 
   rk4eps      = dble(rk4_tolerance) ! The desired accuracy.
   rk4epsi     = 1.d0/rk4eps         ! The inverse of desired accuracy.
   rk4eps2     = rk4eps**2           ! square of the accuracy
   hmin        = 1.d-7               ! The minimum step size.
   print_diags = .false.             ! Flag to print the diagnostic check.
   checkbudget = .false.             ! Flag to check CO2, water, and energy budgets every 
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
   ! ly stable and solvable in a fast way.  If you don't want this and want to use the     !
   ! nominal heat capacity, the laziest way to turn this off is by setting hcapveg_ref to  !
   ! a small number.  Don't set it to zero, otherwise you may have FPE issues.             !
   !---------------------------------------------------------------------------------------!
   hcapveg_ref         = 3.0d3  ! 3.0d-1  ! Reference heat capacity value          [J/m³/K]
   min_height          = 1.5d0            ! Minimum vegetation height              [     m]
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
   rk4min_can_co2    =  2.0000d2  ! Minimum canopy    CO2 mixing ratio          [ µmol/mol]
   rk4max_can_co2    =  1.2000d3  ! Maximum canopy    CO2 mixing ratio          [ µmol/mol]
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
   !     Minimum water mass at the leaf surface.  This is given in kg/m²leaf rather than   !
   ! kg/m²ground, so we scale it with LAI.                                                 !
   !---------------------------------------------------------------------------------------!
   rk4min_veg_lwater = -rk4dry_veg_lwater ! Minimum leaf water mass             [kg/m²leaf]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    The minimum mass of surface water and virtual layer are given in m3/m3 rather than !
   ! kg/m2.  This is because there will be exchange between the top soil layer and the     !
   ! layers above in case the mass goes below the minimum.  Since this would make the im-  !
   ! pact of such exchange dependent on the soil depth, we assign the scale a function of  !
   ! the top layer thickness.                                                              !
   !---------------------------------------------------------------------------------------!
   rk4min_sfcw_moist = -5.0000d-4 ! Minimum water mass allowed.
   rk4min_virt_moist = -5.0000d-4 ! Minimum water allowed at virtual pool.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The integrator will check whether the leaves are at their maximum holding         !
   ! capacity before letting water to be intercepted and dew/frost to remain at the leaf   !
   ! surfaces.  If FALSE, the amount of water will be adjusted outside the derivative      !
   ! calculation.  Some tests have shown that not solving this inside the integrator can   !
   ! speed up the simulation between 20% in the tropics to well above 100% in some         !
   ! temperate polygons. However, some may think that it is more appropriate to leave the  !
   ! check in the integrator, so we leave this flag so one can switch between one and the  !
   ! other and decide based on their own needs.                                            !
   !---------------------------------------------------------------------------------------!
   check_maxleaf = .false.


   !---------------------------------------------------------------------------------------!
   !    This flag is used to control evaporation and transpiration when the air is         !
   ! saturated or super-saturated.  If supersat_ok is TRUE, then evaporation and           !
   ! transpiration will continue to happen even if the air is super-saturated at the       !
   ! canopy air temperature (but not at the soil or vegetation temperature).  Otherwise,   !
   ! evaporation and transpiration will be interrupted until the air becomes sub-saturated !
   ! again.  The air can still become super-saturated because mixing with the free atmo-   !
   ! sphere will not stop.                                                                 !
   !---------------------------------------------------------------------------------------!
   supersat_ok = .false.

   !---------------------------------------------------------------------------------------!
   !     This flag is used to control the method used to define water percolation through  !
   ! snow layers.                                                                          !
   !---------------------------------------------------------------------------------------!
   newsnow = .true.

   return
end subroutine init_rk4_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine overwrite_with_xml_config(thisnode)
   !!! PARSE XML FILE
   use ed_max_dims, only: n_pft
   use ed_misc_coms, only: iedcnfgf
   
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
