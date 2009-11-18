Module ed_max_dims

  ! This module contains very basic specification of grid dimensions and other 
  ! parameters that will be used to dimension arrays and allocate memory.
  ! Grid dimensions:
  !   MAXGRDS - Maximum number of grids
  !   NXPMAX  - Maximum number of points in x-direction
  !   NYPMAX  - Maximum number of points in y-direction
  !   NZPMAX  - Maximum number of points in z-direction
  !   NZGMAX  - Maximum number of soil levels
  !   NZSMAX  - Maximum number of snow/water levels
  !   MAXDIM  - the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX



#if defined(COUPLED)
  use grid_dims, only: brams_maxgrds => maxgrds &
                     , brams_nxpmax  => nxpmax  &
                     , brams_nypmax  => nypmax  &
                     , brams_nzpmax  => nzpmax  &
                     , brams_nzgmax  => nzgmax  &
                     , brams_maxdim  => maxdim  &
                     , brams_maxdimp => maxdimp &
                     , brams_nxyzpm  => nxyzpm  &
                     , brams_maxmach => maxmach

  implicit none

!------------------------------------------------------------------------------------------!
!    In a coupled run, these dimensions must match in ED and BRAMS, so we use preprocessor !
! here.                                                                                    !
!------------------------------------------------------------------------------------------!
  integer, parameter :: maxgrds = brams_maxgrds
  integer, parameter :: nxpmax  = brams_nxpmax 
  integer, parameter :: nypmax  = brams_nypmax 
  integer, parameter :: nzpmax  = brams_nzpmax 
  integer, parameter :: nzgmax  = brams_nzgmax 
  integer, parameter :: maxdim  = brams_maxdim 
  integer, parameter :: maxdimp = brams_maxdimp
  integer, parameter :: nxyzpm  = brams_nxyzpm 
  integer, parameter :: maxmach = brams_maxmach

#else
  integer, parameter :: maxgrds=10
  integer, parameter :: nxpmax=666
  integer, parameter :: nypmax=666
  integer, parameter :: nzpmax=132
  integer, parameter :: nzgmax=25
  integer, parameter :: maxdim=666

  ! Computed parameters (function of previous parameters)
  integer, parameter :: maxdimp=maxdim+2
  integer, parameter :: nxyzpm=nzpmax*nxpmax*nypmax

  !   MAXMACH - the maximum number of processors on a parallel run
  integer, parameter :: maxmach=3000
#endif

  ! Maximum number of temporary water layers
  integer, parameter :: nzsmax=10

  ! Suppose you want to run ED only at a few scattered test sites.  
  ! These are called sites of interest (SOIs).  
  !max_soi is the maximum allowable.
  integer, parameter :: max_soi = 10

  ! Suppose you want to run ED for a few rectangular regions.  
  ! You can run for a global region, if you want.  
  ! You can also run several sub-global regions and several SOIs.  
  !The maximum number of regions is max_ed_regions.
  integer, parameter :: max_ed_regions = 10

  ! Maximum number of atmospheric grid cells that can fit in one ED polygon.
  integer, parameter :: ed_maxatm = 625

  ! Maximum file name length
  integer, parameter :: str_len=256
  
  ! Maximum variables string length
  integer, parameter :: str_len_short=32

  ! Maximum number of sites within a polygon
  integer, parameter :: max_site=1

  integer, parameter :: n_pft = 15 ! number of plant functional types:
  ! 1   C4 grass
  ! 2   early successional broadleaf evergreen
  ! 3   mid successional broadleaf evergreen
  ! 4   late successional broadleaf evergreen
  ! 5   C3 grass
  ! 6   northern pines
  ! 7   southern pines
  ! 8   late successional conifers
  ! 9   early successional broadleaf deciduous
  ! 10  mid successional broadleaf deciduous
  ! 11  late successional broadleaf deciduous
  ! 12 - c3 pasture
  ! 13 - c3 crop (e.g.,wheat, rice, soybean) 
  ! 14 - c4 pasture
  ! 15 - c4 crop (e.g.,corn/maize)

  integer, parameter :: n_dbh = 11 ! Number of DBH bins for output quantities
                                   ! For maxdbh = 100 cm (this is defined in ed_params.f90)
                                   ! the classes will be:
                                   ! 1 - [0;10cm]
                                   ! 2 - ]10;20cm]
                                   ! 3 - ]20;30cm]
                                   ! ...
                                   ! 10 - ]90;100cm]
                                   ! 11 - ]100cm;Infinity[

  integer, parameter :: n_age = 21 ! Number of patch AGE bins for output quantities
                                   ! For maxage = 200 yr (this is defined in ed_params.f90)
                                   ! the classes will be:
                                   ! 1  - [0;10yr]
                                   ! 2  - ]10;20yr]
                                   ! 3  - ]20;30yr]
                                   ! ...
                                   ! 20 - ]190;200yr]
                                   ! 21 - ]200yr;Infinity[

  integer, parameter :: n_dist_types = 3 ! Number of disturbance types:  
                                         ! 1 - agriculture
                                         ! 2 - secondary forest
                                         ! 3 - primary forest.

  integer, parameter :: n_mort = 4 ! Number of mortality types:
                                   ! 1. Ageing, PFT-dependent but otherwise constant;
                                   ! 2. Negative carbon balance;
                                   ! 3. Treefall mortality;
                                   ! 4. Mortality due to cold weather.

  ! Maximum number of model variables
  integer, parameter :: maxvars = 250

  ! Maximum number meteorological driver xfiles
  integer, parameter :: max_met_vars = 22
  
  ! Maximum number of patches and cohorts possibly allowed initially
  integer, parameter :: huge_patch  = 2000
  integer, parameter :: huge_cohort = 20000
  integer, parameter :: max_water = 100
  
  integer,parameter :: maxpvars = 50     ! The maximum number of printable variables

  integer, parameter :: maxfiles = 6666  ! Maximum number of files

  integer, parameter :: maxlist = 3*maxfiles  ! Maximum number of files (site+patch+cohort)
  
end Module ed_max_dims
