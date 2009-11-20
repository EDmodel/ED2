Module decomp_coms


  ! DO NOT INITIALIZE PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM

  use ed_max_dims, only: n_pft

  implicit none

  real :: resp_opt_water !  Optimal soil porosity, as a fraction of total porosity, for heterotrophic respiration (dimensionless).

  real :: resp_water_below_opt ! Determines rate at which heterotrophic respiration declines for relative porosities below resp_opt_water (dimensionless).

  real :: resp_water_above_opt ! Determines rate at which heterotrophic respiration declines for relative porosities above resp_opt_water (dimensionless).

  real :: resp_temperature_increase ! Determines how rapidly heterotrophic respiration increases with increasing temperature (1/K).

  real :: N_immobil_supply_scale ! Supply coefficient for nitrogen immobilization (1/day)

  real :: cwd_frac ! Fraction of structural material that goes to coarse woody debris upon mortality.  Note that currently CWD decomposed at a rate identical to structural soil C.

  real :: r_fsc    ! Fraction of structural pool decomposition going to heterotrophic respiration

  real :: r_stsc   ! Fraction of structural pool decomposition going to heterotrophic respiration instead of the slow pool.

  real :: r_ssc    ! Fraction of structural pool decomposition going to heterotrophic respiration

  real :: K1       ! Intrinsic decay rate of fast pool soil carbon (1/days); this is modulated by Lc

  real :: K2       ! Intrinsic decay rate of fast pool soil carbon (1/days)

  real :: K3       ! Intrinsic decay rate of slow pool soil carbon (1/days).  This pool has already decayed from the structural pool.

  ! Labile fraction of leaves, fine roots and sapwood.
  real, dimension(n_pft) :: f_labile 
    !! moved setting of values to initialize_pft_resp_params [[MCD]]

  ! Flag specifying whether or not decomposition is to be limited by 
  ! nitrogen availability.
  integer :: n_decomp_lim

  ! Specifies whether to use Lloyd and Taylor (1994) temperature dependence
  logical :: LloydTaylor

end Module decomp_coms
