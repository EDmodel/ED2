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
   !     Fraction of structural pool decomposition going to heterotrophic respiration.     !
   !---------------------------------------------------------------------------------------!
   real :: r_fsc
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration.     !
   !---------------------------------------------------------------------------------------!
   real :: r_stsc
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration.     !
   !---------------------------------------------------------------------------------------!
   real :: r_msc
   !---------------------------------------------------------------------------------------!
   !     Fraction of structural pool decomposition going to heterotrophic respiration.     !
   ! This is disabled because it must be 1.0 for the original ED-2 soil decomposition      !
   ! schemes, and it is not used by the new scheme based on RothC (Sierra et al. 2012).    !
   !---------------------------------------------------------------------------------------!
   ! real :: r_ssc
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
   !     Intrinsic decay rate of microbial soil carbon (1/days).                           !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_msc
   !---------------------------------------------------------------------------------------!
   !     Intrinsic decay rate of slow pool soil carbon (1/days).  This pool has already    !
   ! decayed from the structural pool.                                                     !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_ssc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Fraction of decay that is transferred to microbial carbon pool.  This parameter   !
   ! is used only when decomp_scheme is set to 2.                                          !
   !---------------------------------------------------------------------------------------!
   real :: fx_msc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Coefficients to obtain parameter x for the RothC model as implemented by Sierra   !
   ! et al. (2012).                                                                        !
   !---------------------------------------------------------------------------------------!
   real :: xrothc_a
   real :: xrothc_b
   real :: xrothc_c
   real :: xrothc_d
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Temporary parameters to indicate above-ground (and flammable) necromass.  In the  !
   ! future we should split the necromass pools into above- and below-ground.              !
   !---------------------------------------------------------------------------------------!
   real :: agf_fsc
   real :: agf_stsc
   !---------------------------------------------------------------------------------------!


   !----- Carbon to Nitrogen ratio, structural pool. --------------------------------------!
   real :: c2n_structural
   !----- Carbon to Nitrogen ratio, slow pool. --------------------------------------------!
   real :: c2n_slow
   !---------------------------------------------------------------------------------------!

end Module decomp_coms
!==========================================================================================!
!==========================================================================================!
