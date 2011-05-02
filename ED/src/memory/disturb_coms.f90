!==========================================================================================!
!==========================================================================================!
!   module disturb_coms                                                                    !
!   This module contains variables used to control the disturbance rates.                  !
!   N.B.: Variables that are not parameters should be initialized in the subroutine        !
!         initialize_disturb_params in ed_params.f90, since some compilers don't actually  !
!         initialize variables in modules.                                                 !
!------------------------------------------------------------------------------------------!
module disturb_coms
   use ed_max_dims, only : str_len & ! intent(in)
                         , maxgrds & ! intent(in)
                         , n_pft   ! ! intent(in)
   implicit none


   !=======================================================================================!
   !=======================================================================================!
   !    General parameters.                                                                !
   !---------------------------------------------------------------------------------------!

   !------ Number of land use transitions in the input land use disturbance dataset. ------!
   integer, parameter :: num_lu_trans = 19

   !---------------------------------------------------------------------------------------!
   !     Maximum number of years in which land use can be applied.  In case the simulation !
   ! runs longer than this, the missing years will be filled with zeroes.  The first and   !
   ! last year of each is checked in landuse_init.                                         !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: max_lu_years = 1000 

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Namelist variables.                                                               !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Fire model.  Possible values are:                                                 !
   ! 0. No fires;                                                                          !
   ! 1. (deprecated) Fire will be triggered with enough biomass and integrated ground      !
   !    water depth less than a threshold.  Based on ED-1, the threshold assumes that the  !
   !    soil is 1 m, so deeper soils will need to be much drier to allow fires to happen   !
   !    and often will never allow fires because the threshold may be below the minimum    !
   !    possible soil moisture.                                                            !
   ! 2. Fire will be triggered with enough biomass and the total soil water at the top 75  !
   !    cm falls below a (relative) threshold.                                             !
   !---------------------------------------------------------------------------------------!
   integer :: include_fire

   !---------------------------------------------------------------------------------------!
   !     Anthropogenic disturbance.  1 means that anthropogenic disturbances will be       !
   ! included, whereas 0 means that it won't.                                              !
   !---------------------------------------------------------------------------------------!
   integer :: ianth_disturb


   !---------------------------------------------------------------------------------------!
   !     Treefall disturbance                                                              !
   ! > 0. Usual disturbance rate, in 1/years, at which treefall gaps form.                 !
   ! = 0. No treefall disturbance;                                                         !
   ! < 0. Treefall will be added as a mortality rate (it will kill plants, but it won't    !
   !      create a new patch).                                                             !
   !---------------------------------------------------------------------------------------!
   real :: treefall_disturbance_rate
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is the time until we start knocking down trees.  Used only when              !
   ! treefall_disturbance_rate > 0.                                                        !
   !---------------------------------------------------------------------------------------!
   real :: time2canopy
   !---------------------------------------------------------------------------------------!



   !----- The prefix for land use disturbance rates. The path and prefix must be included. !
   character(len=str_len), dimension(maxgrds) :: lu_database 
   !----- File with plantation fraction.  If no file is available, leave it blank. --------!
   character(len=str_len), dimension(maxgrds) :: plantation_file
   !----- File with initial land use area scale.  If no file is available, leave it blank. !
   character(len=str_len), dimension(maxgrds) :: lu_rescale_file
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Patch dynamics variables, to be set in ed_params.f90.                              !
   !---------------------------------------------------------------------------------------!
   !----- Minimum relative area required for a patch to be created or maintained. ---------!
   real :: min_new_patch_area 
   !----- Only trees above this height create a gap when they fall. -----------------------!
   real :: treefall_hite_threshold
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Forestry variables, to be set in ed_params.f90.                                    !
   !---------------------------------------------------------------------------------------!
   !----- Set to 1 if to do forest harvesting. --------------------------------------------!
   integer :: forestry_on
   !----- Set to 1 if to do agriculture. --------------------------------------------------!
   integer :: agriculture_on
   !----- Earliest year at which plantations occur. ---------------------------------------!
   integer :: plantation_year
   !----- Number of years that a plantation requires to reach maturity. -------------------!
   real :: plantation_rotation
   !----- Years that a non-plantation patch requires to reach maturity. -------------------!
   real :: mature_harvest_age
   !----- Minimum plantation fraction to consider the site a plantation. ------------------!
   real :: min_plantation_frac
   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon for plantation fraction.                          !
   !---------------------------------------------------------------------------------------!
   real :: max_plantation_dist
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Fire parameters.                                                                  !
   !---------------------------------------------------------------------------------------!

   !----- Dimensionless parameter controlling speed of fire spread. -----------------------!
   real :: fire_parameter

   !---------------------------------------------------------------------------------------!
   !     Fire may occur if total equivalent water depth (ground + underground) falls below !
   ! this threshold and include_fire is 1.  Units: meters.                                 !
   !---------------------------------------------------------------------------------------!
   real :: fire_dryness_threshold 

   !---------------------------------------------------------------------------------------!
   !     Fire may occur when the total water (ground + underground) converted to           !
   ! equivalent average soil moisture is below this threshold and include_fire is 2.       !
   ! Units: relative fraction.                                                             !
   !---------------------------------------------------------------------------------------!
   real :: sm_fire
   real :: fire_smoist_threshold
   !---------------------------------------------------------------------------------------!
   !     Depth to be compared with the soil average when include_fire is 2. Units: meters. !
   !---------------------------------------------------------------------------------------!
   real :: fire_smoist_depth         

   !----- k level of the deepest layer to be considered. ----------------------------------!
   integer :: k_fire_first
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Variable type that contains the land use change disturbances.                     !
   !---------------------------------------------------------------------------------------!
   type lutime
      !------ The disturbance year. -------------------------------------------------------!
      integer :: landuse_year ! the year

      !------------------------------------------------------------------------------------!
      !    The land use information:  the current columns are:                             !
      !  ====== Disturbance rates ======                                                   !
      !  1 - Cropland to pasture                                           [         1/yr] !
      !  2 - Pasture to cropland                                           [         1/yr] !
      !  3 - Pasture to primary forest                                     [         1/yr] !
      !  4 - Primary forest to pasture                                     [         1/yr] !
      !  5 - Primary forest to cropland                                    [         1/yr] !
      !  6 - Cropland to primary forest                                    [         1/yr] !
      !  7 - Secondary forest to cropland                                  [         1/yr] !
      !  8 - Cropland to secondary forest                                  [         1/yr] !
      !  9 - Secondary forest to pasture                                   [         1/yr] !
      ! 10 - Pasture to secondary forest                                   [         1/yr] !
      ! 11 - Primary forest to secondary forest                            [         1/yr] !
      !  ====== Biomass to be harvested. ======                                            !
      ! 12 - Wood harvest on primary forested land.                        [          kgC] !
      ! 13 - Wood harvest on primary forested land.                        [       kgC/m²] !
      ! 14 - Wood harvest on mature secondary forest land.                 [          kgC] !
      ! 15 - Wood harvest on mature secondary forest land.                 [       kgC/m²] !
      ! 16 - Wood harvest on primary non-forested land.                    [          kgC] !
      ! 17 - Wood harvest on primary non-forested land.                    [       kgC/m²] !
      ! 18 - Wood harvest on young secondary forest land.                  [          kgC] !
      ! 19 - Wood harvest on young secondary forest land.                  [       kgC/m²] !
      !  ====== Special flags. ======                                                      !
      ! 12 - Primary forest is harvested using the probability of harvesting when the DBH  !
      !      is above the minimum DBH.                                                     !
      ! 14 - Secondary forest is harvested using the probability of harvesting when the    !
      !      DBH is above the minimum DBH.                                                 !
      !------------------------------------------------------------------------------------!
      real, dimension(num_lu_trans) :: landuse
   end type lutime
   !=======================================================================================!
   !=======================================================================================!
end module disturb_coms
!==========================================================================================!
!==========================================================================================!
