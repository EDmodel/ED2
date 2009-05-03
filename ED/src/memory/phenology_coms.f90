Module phenology_coms
  implicit none

  ! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM
  ! See "initialize_phen_coms" for initial values

  real :: retained_carbon_fraction  !  Before plants drop their leaves, they retain this fraction of their leaf carbon and nitrogen and put it into storage.

  real :: theta_crit                !  When soil porosity (relative to total soil porosity) drops below, this threshold, drought-deciduous plants drop their leaves.

  ! leaf offset parameters are from White et al. 1997 Global
  ! Biogeochemical Cycles 11(2) 217-234 
  real :: dl_tr             ! critical daylength in minutes
  real :: st_tr1            !  critical soil temp
  real :: st_tr2            ! second critical soil temp
  
  ! Phenology parameters for cold deciduous trees
  ! Botta et al. 2000, Global Change Biology, 6, 709--725
  real :: phen_a
  real :: phen_b
  real :: phen_c

  ! Flag specifying which phenology scheme to run.
  ! 0 - predictive scheme for the temperate
  ! 1 - user prescribed, remote sensing driven
  ! 2 - predictive scheme for the tropics
  integer :: iphen_scheme

  ! Flag specifying which reproduction scheme to run
  ! 0 - no reproduction
  ! 1 - original ED1 reproduction
  integer :: repro_scheme

  real :: l_fract

  !Flag specifying the first and last spring
  integer :: iphenys1,iphenysf,iphenyf1,iphenyff 

  character(len=512) :: phenpath

  ! Light-controlled
  real :: rad_turnover_int
  real :: rad_turnover_slope  

  real :: vm_tran
  real :: vm_slop
  real :: vm_amp
  real :: vm_min

  ! Derived type describing prescribed phenology.
  type prescribed_phen
     
     ! number of years for which we have prescribed phenology
     integer :: nyears
     
     ! the years for which we have prescribed phenology
     integer, dimension(:), pointer :: years
     
     ! Two parameters of the springtime logistic function describing leaf flush
     real, dimension(:), pointer :: flush_a
     real, dimension(:), pointer :: flush_b
     
     ! Two parameters of the autumn logistic function describing leaf color
     real, dimension(:), pointer :: color_a
     real, dimension(:), pointer :: color_b
     
  end type prescribed_phen
  

end Module phenology_coms
