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
   ! DECOMP_SCHEME -- This specifies the soil Carbon (decomposition) model.                !
   !                  0.  ED-2.1 default, with three soil carbon pools.  Temperature the   !
   !                      original exponential (low temperature limitation only).          !
   !                  1.  Similar to option 0, except that temperature is solved using     !
   !                      the Lloyd and Taylor (1994) model.                               !
   !                      [[option 1 requires parameters to be set in xml]]                !
   !                  2.  Similar to ED-1.0 and CENTURY model, heterotrophic respiration   !
   !                      reaches a maximum at around 38C (using the default parameters),  !
   !                      then quickly falls to zero at around 50C.  It applies a similar  !
   !                      function for soil moisture, which allows higher decomposition    !
   !                      rates when it is close to the optimal, plumetting when it is     !
   !                      almost saturated.                                                !
   !                  3.  Similar to option 0. Uses empirical moisture limit equation from !
   !                      Moyano et al., 2012, Biogeosciences.                             !
   !                  4.  Similar to option 1. Uses empirical moisture limit equation from !
   !                      Moyano et al., 2012, Biogeosciences.                             !
   !                  5.  Based on Bolker et al. (1998) CENTURY model.  Five necromass     !
   !                      pools (litter aka fast, structural, microbial,                   !
   !                      humified aka slow, and passive).  Temperature and moisture       !
   !                      functions are the same as 2.                                     !
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
   !     The following variables are used when DECOMP_SCHEME is 3 or 4. (based on M12).    !
   !                                                                                       !
   !  Moyano FE, Vasilyeva N, Bouckaert L, Cook F, Craine J, Curiel Yuste J, Don A,        !
   !     Epron D, Formanek P, Franzluebbers A et al. 2012. The moisture response of soil   !
   !     heterotrophic respiration: interaction with soil properties. Biogeosciences, 9:   !
   !     1173-1182. doi:10.5194/bg-9-1173-2012 (M12).                                      !
   !---------------------------------------------------------------------------------------!
   real :: rh_moyano12_a0
   real :: rh_moyano12_a1
   real :: rh_moyano12_a2
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
   !     Fraction of microbial soil decomposition going to heterotrophic respiration.      !
   ! Two values are provided because the actual value depends on the sand content.         !
   !---------------------------------------------------------------------------------------!
   real :: r_msc_int
   real :: r_msc_slp
   !---------------------------------------------------------------------------------------!
   !     Fraction of humified soil decomposition going to heterotrophic respiration.       !
   !---------------------------------------------------------------------------------------!
   real :: r_ssc
   !---------------------------------------------------------------------------------------!
   !     Fraction of passive soil decomposition going to heterotrophic respiration.        !
   !---------------------------------------------------------------------------------------!
   real :: r_psc
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
   !     Intrinsic decay rate of slow (humified) pool soil carbon (1/days).                !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_ssc
   !---------------------------------------------------------------------------------------!
   !     Intrinsic decay rate of passive pool soil carbon (1/days).                        !
   !---------------------------------------------------------------------------------------!
   real :: decay_rate_psc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Fraction of decay that is transferred between soil pools.  These quantities are   !
   ! functions of clay content, following CENTURY, and thus both the intercept and the     !
   ! slope must be provided.  This option is used only when DECOMP_SCHEME is 2.            !
   !---------------------------------------------------------------------------------------!
   real :: fx_msc_psc_int
   real :: fx_msc_psc_slp
   real :: fx_ssc_psc_int
   real :: fx_ssc_psc_slp
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !  CENTURY parameter that controls the effect of lignin on structural decomposition.    !
   !---------------------------------------------------------------------------------------!
   real :: e_lignin
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters to indicate above-ground (and flammable) necromass.  These are used    !
   ! only during initialisation (NBG or pss/css files).                                    !
   !---------------------------------------------------------------------------------------!
   real :: agf_fsc
   real :: agf_stsc
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters to indicate the fraction of soil carbon that should be microbial and   !
   ! passive.  Humified (slow) will be 1 - f0_msc - f0_psc.  These are used only during    !
   ! initialisation (NBG or pss/css files), and only when DECOMP_SCHEME = 2.               !
   !---------------------------------------------------------------------------------------!
   real :: f0_msc
   real :: f0_psc
   real :: f0_ssc
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters to initialise soil carbon pools for near bare ground simulations when  !
   ! N limitation is turned on.                                                            !
   !---------------------------------------------------------------------------------------!
   real :: nbg_nlim_fsc
   real :: nbg_nlim_stsc
   real :: nbg_nlim_ssc
   !---------------------------------------------------------------------------------------!


   !----- Carbon to Nitrogen ratio, fast pool (for initial conditions only). --------------!
   real :: c2n_fast_0
   !----- Carbon to Nitrogen ratio, structural pool. --------------------------------------!
   real :: c2n_structural
   !----- Carbon to Nitrogen ratio, slow pool (or microbial, slow, and passive). ----------!
   real :: c2n_slow
   !---------------------------------------------------------------------------------------!

end Module decomp_coms
!==========================================================================================!
!==========================================================================================!
