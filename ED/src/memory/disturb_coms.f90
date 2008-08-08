Module disturb_coms
  implicit none

  ! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM
  ! See "initialize_disturb_params" for settings


  ! GENERAL PARAMETERS
  !--------------------------

  integer :: patch_dynamics    !  Set to 1 to incorporate the  
  ! effects of disturbance, and to do patch fusion.

  real :: min_new_patch_area        !  minimum fractional area 
  ! required to form a new patch.

  integer, parameter :: num_lu_trans = 19 ! number of different types 
  ! of land use transitions in George Hurtt's GLU data set.

  integer, parameter :: max_lu_years = 300 ! Used to hold the lu transition array


  integer :: include_fire ! flag specifying whether or not to include fire

  integer :: ianth_disturb ! flag specifying whether or not to include 
  ! anthropogenic disturbance.

  ! TREEFALL DISTURBANCE
  !--------------------------

  real :: treefall_disturbance_rate  ! Rate (1/years) at which treefall gaps form. Read from ED_NL.
  
  real :: treefall_hite_threshold    !  Only trees above this height create a gap when they fall.
  
  real :: treefall_age_threshold     !  Minimum patch age for treefall disturbance.
  
  ! FORESTRY
  !--------------------------
  
  integer :: forestry_on          ! Set to 1 if to do forest harvesting.
  integer :: agriculture_on       ! Set to 1 if to do agriculture.
  integer :: plantation_year      ! Earliest year at which plantations occur
  real :: plantation_rotation     ! Number of years that a plantation requires to reach maturity
  real :: mature_harvest_age      ! Years that a non-plantation patch requires to reach maturity
  
  ! FIRE
  !--------------------------
  
!!!  integer :: fire_on = 1  ! Set to 1 if to do fire  [[DEPRECATED]]  (mcd)
  real :: fire_dryness_threshold  !  (meters) Fire can occur if total soil water falls below this threshold.
  real :: fire_parameter          ! Dimensionless parameter controlling speed of fire spread.
  
  type lutime
     integer :: landuse_year ! the year
     real, dimension(num_lu_trans) :: landuse  ! the landuse information 
     !(e.g., area harvested, biomass harvested, area abandoned, 
     ! area converted to agriculture, etc.)
     type(lutime), pointer :: next_lutime ! pointer to the next year
  end type lutime


end Module disturb_coms
