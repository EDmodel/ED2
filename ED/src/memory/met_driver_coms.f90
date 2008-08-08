Module met_driver_coms

  use max_dims, only: str_len
  integer, parameter :: metname_len = 128
  integer, parameter :: metvars_len =  16
  logical            :: have_co2
  logical            :: have_latlon   ! Does the H5 data have lat-lon
                                      ! grids or does the user specify
                                      ! its spacing directly?
  real               :: initial_co2

  integer :: nformats
  character(len=metname_len), allocatable, dimension(:) :: met_names
  integer, allocatable, dimension(:) :: met_nlon
  integer, allocatable, dimension(:) :: met_nlat
  real, allocatable, dimension(:) :: met_dx
  real, allocatable, dimension(:) :: met_dy
  real, allocatable, dimension(:) :: met_xmin
  real, allocatable, dimension(:) :: met_ymin
  integer, allocatable, dimension(:) :: met_nv
  character(len=metvars_len), allocatable, dimension(:,:) :: met_vars
  real, allocatable, dimension(:,:) :: met_frq
  integer, allocatable, dimension(:,:) :: met_interp
  integer :: metcyc1
  integer :: metcycf
  character(len=str_len) :: ed_met_driver_db
  integer :: imettype
  logical :: no_ll

  real, allocatable, dimension(:,:) :: lat2d
  real, allocatable, dimension(:,:) :: lon2d

  Type met_driv_data

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

  end Type met_driv_data

!! separated variables that store the met data (met_driv_data)
!! from those that recode it's instantaneous state (met_driv_state) [MCD]

  Type met_driv_state  
     real :: nir_beam
     real :: nir_diffuse
     real :: par_beam
     real :: par_diffuse
     real :: atm_tmp
     real :: atm_shv
     real :: rhos
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

  end Type met_driv_state

  type(met_driv_state) :: lapse

end Module met_driver_coms
