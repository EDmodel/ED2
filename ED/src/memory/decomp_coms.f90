!==========================================================================================!
!==========================================================================================!
!     Module decomp_coms.   These variables control the heterotrophic respiration.         !
!                                                                                          !
!     IMPORTANT: do not initialise parameters in the module unless they are constants      !
!                (defined with the "parameter" attribute").  Not all compilers will assign !
!                the values here.  The proper location to assign the initial values is in  !
!                ed_params.f90.                                                            !
!------------------------------------------------------------------------------------------!
Module decomp_coms

   use ed_max_dims, only: n_pft

   implicit none

   !=======================================================================================!
   !=======================================================================================!
   !     Parameters that will be read in the namelist.                                     !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   ! N_DECOMP_LIM -- This controls whether decomposition can be limited by nitrogen.       !
   !                0.  No limitation                                                      !
   !                1.  ED-2.1 nitrogen limitation model.                                  !
   !---------------------------------------------------------------------------------------!
   integer :: n_decomp_lim
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! DECOMP_SCHEME -- This specifies the dependence of soil decomposition on temperature.  !
   !                  0.  ED-2.1 default, the original exponential (low temperature        !
   !                      limitation only).                                                !
   !                  1.  Lloyd and Taylor (1994) model                                    !
   !                      [[option 1 requires parameters to be set in xml]]                !
   !                  2.  Similar to ED-1.0 and CENTURY model, heterotrophic respiration   !
   !                      reaches a maximum at around 38C (using the default parameters),  !
   !                      then quickly falls to zero at around 50C.                        !
   !---------------------------------------------------------------------------------------!
   integer :: decomp_scheme
   !---------------------------------------------------------------------------------------!
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     Other variables.                                                                  !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Optimal soil porosity, as a fraction of total porosity, for heterotrophic         !
   ! respiration (dimensionless).                                                          !
   !---------------------------------------------------------------------------------------!
   real :: resp_opt_water
   !---------------------------------------------------------------------------------------!
   !     Determines rate at which heterotrophic respiration declines for relative          !
   ! porosities below resp_opt_water (dimensionless).                                      !
   !---------------------------------------------------------------------------------------!
   real :: resp_water_below_opt
   !---------------------------------------------------------------------------------------!
   !     Determines rate at which heterotrophic respiration declines for relative          !
   ! porosities above resp_opt_water (dimensionless).                                      !
   !---------------------------------------------------------------------------------------!
   real :: resp_water_above_opt
   !---------------------------------------------------------------------------------------!
   !     When DECOMP_SCHEME is set to 0 or 1, it determines how rapidly heterotrophic      !
   ! respiration increases with increasing temperature (1/K).  Notice that between 0 and 1 !
   ! there is a significant change in the functional form, so they don't necessarily have  !
   ! the same meaning...                                                                   !
   !---------------------------------------------------------------------------------------!
   real :: resp_temperature_increase
   !---------------------------------------------------------------------------------------!
   !     The following variables are used when DECOMP_SCHEME is 1. (LLoyd and Taylor 1994) !
   !---------------------------------------------------------------------------------------!
   real :: rh_lloyd_1
   real :: rh_lloyd_2
   real :: rh_lloyd_3
   !---------------------------------------------------------------------------------------!
   !     The following variables are used when DECOMP_SCHEME is 2. (based on CENTURY model !
   ! and ED-1.0).                                                                          !
   !---------------------------------------------------------------------------------------!
   real    :: rh_decay_low     !  Low temperature decay rate                        [  1/K]
   real    :: rh_decay_high    !  High temperature decay rate                       [  1/K]
   real    :: rh_low_temp      !  Low temperature reference                         [    K]
   real    :: rh_high_temp     !  High temperature reference                        [    K]
   real    :: rh_decay_dry     !  Decay rate for dry soil                           [   --]
   real    :: rh_decay_wet     !  Decay rate for wet soil                           [   --]
   real    :: rh_dry_smoist    !  Dry relative soil moisture threshold              [   --]
   real    :: rh_wet_smoist    !  Wet relative soil moisture threshold              [   --]
   real    :: rh_active_depth  !  Maximum depth for avg. temperature and moisture   [    m]
   integer :: k_rh_active      !  Index of the bottommost layer                     [   --]
   !---------------------------------------------------------------------------------------!
   !     Supply coefficient for nitrogen immobilization (1/day).                           !
   !---------------------------------------------------------------------------------------!
   real :: N_immobil_supply_scale 
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural material that goes to coarse woody debris upon mortality.  !
   ! Note that currently CWD decomposed at a rate identical to structural soil C.          !
   !---------------------------------------------------------------------------------------!
   real :: cwd_frac
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration.     !
   !---------------------------------------------------------------------------------------!
   real :: r_fsc
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration      !
   ! instead of the slow pool.                                                             !
   !---------------------------------------------------------------------------------------!
   real :: r_stsc
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration      !
   !---------------------------------------------------------------------------------------!
   real :: r_ssc
   !---------------------------------------------------------------------------------------!
   !     Intrinsic decay rate of structural pool soil carbon (1/days); this is modulated   !
   ! by Lc.                                                                                !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_stsc
   !---------------------------------------------------------------------------------------!
   !     Intrinsic decay rate of fast pool soil carbon (1/days).                           !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_fsc
   !---------------------------------------------------------------------------------------!
   !     Intrinsic decay rate of slow pool soil carbon (1/days).  This pool has already    !
   ! decayed from the structural pool.                                                     !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_ssc
   !---------------------------------------------------------------------------------------!
   !     Labile fraction of leaves, fine roots and sapwood.                                !
   !     ([[MCD]].  Moved setting of values to initialize_pft_resp_params)                 !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: f_labile 
   !---------------------------------------------------------------------------------------!

end Module decomp_coms
!==========================================================================================!
!==========================================================================================!
