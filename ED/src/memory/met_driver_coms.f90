!==========================================================================================!
!==========================================================================================!
!     This module contains the variables that define the meteorological state that will    !
! drive ED-2.                                                                              !
!------------------------------------------------------------------------------------------!
module met_driver_coms

   use ed_max_dims, only: str_len
   integer, parameter :: metname_len = 128
   integer, parameter :: metvars_len =  16
   logical            :: have_co2
   logical            :: have_latlon   ! Does the H5 data have lat-lon grids or does the 
                                       ! user specify  its spacing directly?
   real               :: initial_co2
   integer            :: lapse_scheme  ! which variant of lapse rate correction to use

   !---------------------------------------------------------------------------------------!
   !      Values for simple climate change scenarios                                       !
   !      X_future = intercept + slope * X_current                                         !
   !---------------------------------------------------------------------------------------!
   real               :: atm_tmp_intercept
   real               :: atm_tmp_slope
   real               :: prec_intercept
   real               :: prec_slope
   integer            :: humid_scenario    ! 0 = constant sh
                                           ! 1 = constant RH

   integer :: nformats
   character(len=metname_len), allocatable, dimension(:)   :: met_names
   character(len=metvars_len), allocatable, dimension(:,:) :: met_vars
   integer                   , allocatable, dimension(:)   :: met_nlon
   integer                   , allocatable, dimension(:)   :: met_nlat
   real                      , allocatable, dimension(:)   :: met_dx
   real                      , allocatable, dimension(:)   :: met_dy
   real                      , allocatable, dimension(:)   :: met_xmin
   real                      , allocatable, dimension(:)   :: met_ymin
   integer                   , allocatable, dimension(:)   :: met_nv
   real                      , allocatable, dimension(:,:) :: met_frq
   integer                   , allocatable, dimension(:,:) :: met_interp
   integer                   , allocatable, dimension(:)   :: metyears
   integer                                                 :: ishuffle
   integer                                                 :: metcyc1
   integer                                                 :: metcycf
   character(len=str_len)                                  :: ed_met_driver_db
   integer                                                 :: imettype
   logical                                                 :: no_ll
   real                      , allocatable, dimension(:,:) :: lat2d
   real                      , allocatable, dimension(:,:) :: lon2d


   !---------------------------------------------------------------------------------------!
   !    Declaration of the structure that will hold the input data as it is read from the  !
   ! meteorological driver.                                                                !
   !---------------------------------------------------------------------------------------!
   type met_driv_data
      real, pointer, dimension(:) :: nbdsf
      real, pointer, dimension(:) :: nddsf
      real, pointer, dimension(:) :: vbdsf
      real, pointer, dimension(:) :: vddsf
      real, pointer, dimension(:) :: prate
      real, pointer, dimension(:) :: dlwrf
      real, pointer, dimension(:) :: pres
      real, pointer, dimension(:) :: hgt
      real, pointer, dimension(:) :: ugrd
      real, pointer, dimension(:) :: vgrd
      real, pointer, dimension(:) :: sh
      real, pointer, dimension(:) :: tmp
      real, pointer, dimension(:) :: co2
   end type met_driv_data
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Separated variables that store the meteorological state as it is used by the      !
   ! model and it is updated frequently.                                                   !
   !---------------------------------------------------------------------------------------!
   type met_driv_state  
      real :: nir_beam
      real :: nir_diffuse
      real :: par_beam
      real :: par_diffuse
      real :: atm_tmp
      real :: atm_theta
      real :: atm_enthalpy
      real :: atm_shv
      real :: theta
      real :: rshort
      real :: rshort_diffuse
      real :: rlong
      real :: pcpg
      real :: qpcpg
      real :: dpcpg
      real :: vels
      real :: vels_stab
      real :: vels_unstab
      real :: prss
      real :: exner
      real :: geoht
      real :: atm_co2
      real :: pptnorm
   end type met_driv_state
   !---------------------------------------------------------------------------------------!

   type(met_driv_state) :: lapse
  
   !---------------------------------------------------------------------------------------!
   !     The following variables define the lower and upper limits accepted by ED-2.  The  !
   ! actual values can be set in ed_params.f90.                                            !
   !---------------------------------------------------------------------------------------!
   real :: rshort_min
   real :: rshort_max
   real :: rlong_min
   real :: rlong_max
   real :: atm_tmp_min
   real :: atm_tmp_max
   real :: atm_shv_min
   real :: atm_shv_max
   real :: atm_rhv_min
   real :: atm_rhv_max
   real :: atm_co2_min
   real :: atm_co2_max
   real :: prss_min
   real :: prss_max
   real :: pcpg_min
   real :: pcpg_max
   real :: vels_min
   real :: vels_max
   real :: geoht_min
   real :: geoht_max

end module met_driver_coms
