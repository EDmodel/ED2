!==========================================================================================!
!==========================================================================================!
! Module ed_max_dims: this module contains very basic specification of grid dimensions and !
!                     other parameters that will be used to dimension arrays and allocate  !
!                     memory.                                                              !
!------------------------------------------------------------------------------------------!
module ed_max_dims

#if defined(COUPLED)
   use grid_dims, only : brams_maxgrds  => maxgrds   & ! intent(in)
                       , brams_nxpmax   => nxpmax    & ! intent(in)
                       , brams_nypmax   => nypmax    & ! intent(in)
                       , brams_nzpmax   => nzpmax    & ! intent(in)
                       , brams_nzgmax   => nzgmax    & ! intent(in)
                       , brams_maxdim   => maxdim    & ! intent(in)
                       , brams_maxdimp  => maxdimp   & ! intent(in)
                       , brams_nxyzpm   => nxyzpm    & ! intent(in)
                       , brams_maxmach  => maxmach   & ! intent(in)
                       , brams_maxfiles => maxfiles  & ! intent(in)
                       , brams_str_len  => str_len   ! ! intent(in)
#endif

   implicit none

   !---------------------------------------------------------------------------------------!
   !     Here we must check whether this is the coupled model or not.  If it is, we must   !
   ! use the very same numbers for BRAMS and ED for a few variables.  We will do this by   !
   ! assigning the BRAMS variables to ED in the coupled model.  The variables that will be !
   ! made compatible are:                                                                  !
   !                                                                                       !
   !   MAXGRDS - Maximum number of grids                                                   !
   !   NXPMAX  - Maximum number of points in x-direction                                   !
   !   NYPMAX  - Maximum number of points in y-direction                                   !
   !   NZPMAX  - Maximum number of points in z-direction                                   !
   !   NZGMAX  - Maximum number of soil levels                                             !
   !   MAXDIM  - the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX                             !
   !   MAXDIMP - MAXDIM + 2                                                                !
   !   NXYZPM  - Maximum number of volume points                                           !
   !   MAXMACH - Maximum number of cores on a parallel run.                                !
   !---------------------------------------------------------------------------------------!
#if defined(COUPLED)
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
   integer, parameter :: maxgrds = 10
   integer, parameter :: nxpmax  = 666
   integer, parameter :: nypmax  = 666
   integer, parameter :: nzpmax  = 132
   integer, parameter :: nzgmax  = 100
   integer, parameter :: maxdim  = 666
   integer, parameter :: maxdimp = maxdim + 2
   integer, parameter :: nxyzpm  = nzpmax * nxpmax * nypmax
   integer, parameter :: maxmach = 3000
#endif
   !---------------------------------------------------------------------------------------!



   !----- Maximum number of temporary water layers. ---------------------------------------!
   integer, parameter :: nzsmax=10
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Suppose you want to run ED only at a few scattered test polygons.  These are      !
   ! called polygons of interest (POIs), formerly known as sites of interest (SOI) when    !
   ! there was no distinction between sites and polygons.  max_poi is the maximum allowed. !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: max_poi = 10
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Suppose you want to run ED for a few rectangular regions.  You can run for a      !
   ! global region, if you want.  You can also run several sub-global regions and several  !
   ! POIs.  The maximum number of regions is max_ed_regions.                               !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: max_ed_regions = 10
   !---------------------------------------------------------------------------------------!


   !----- Maximum number of atmospheric grid cells that can fit in one ED polygon. --------!
   integer, parameter :: ed_maxatm = 625
   !---------------------------------------------------------------------------------------!


   !----- Maximum file name length. -------------------------------------------------------!
#if defined(COUPLED)
   integer, parameter :: str_len = brams_str_len
#else
   integer, parameter :: str_len = 300
#endif
   !---------------------------------------------------------------------------------------!


   !----- Maximum variables string length. ------------------------------------------------!
   integer, parameter :: str_len_short=32
   !---------------------------------------------------------------------------------------!


   !----- Maximum number of sites within a polygon. ---------------------------------------!
   integer, parameter :: max_site = 1
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Number of Plant Function Types (PFTs).                                             !
   !                                                                                       !
   ! 1   C4 grass                                                                          !
   ! 2   early successional broadleaf evergreen                                            !
   ! 3   mid successional broadleaf evergreen                                              !
   ! 4   late successional broadleaf evergreen                                             !
   ! 5   C3 grass                                                                          !
   ! 6   northern pines                                                                    !
   ! 7   southern pines                                                                    !
   ! 8   late successional conifers                                                        !
   ! 9   early successional broadleaf deciduous                                            !
   ! 10  mid successional broadleaf deciduous                                              !
   ! 11  late successional broadleaf deciduous                                             !
   ! 12 - c3 pasture                                                                       !
   ! 13 - c3 crop (e.g.,wheat, rice, soybean)                                              !
   ! 14 - c4 pasture                                                                       !
   ! 15 - c4 crop (e.g.,corn/maize)                                                        !
   ! 16 - Subtropical C3 grass                                                             !
   ! 17 - Araucaria                                                                        !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_pft = 17
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Number of DBH bins for output quantities                                              !
   !                                                                                       !
   ! For maxdbh = 100 cm (this is defined in ed_params.f90)                                !
   ! the classes will be:                                                                  !
   ! 1 - [0;10cm]                                                                          !
   ! 2 - ]10;20cm]                                                                         !
   ! 3 - ]20;30cm]                                                                         !
   ! ...                                                                                   !
   ! 10 - ]90;100cm]                                                                       !
   ! 11 - ]100cm;Infinity[                                                                 !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_dbh = 11 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Number of patch AGE bins for output quantities                                        !
   !                                                                                       !
   ! For maxage = 200 yr (this is defined in ed_params.f90)                                !
   ! the classes will be:                                                                  !
   ! 1  - [0;10yr]                                                                         !
   ! 2  - ]10;20yr]                                                                        !
   ! 3  - ]20;30yr]                                                                        !
   ! ...                                                                                   !
   ! 20 - ]190;200yr]                                                                      !
   ! 21 - ]200yr;Infinity[                                                                 !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_age = 21 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Number of disturbance types:                                                          !
   !                                                                                       !
   ! 1 -- Clear cut (cropland and pasture).                                                !
   ! 2 -- Forest plantation.                                                               !
   ! 3 -- Tree fall.                                                                       !
   ! 4 -- Fire.                                                                            !
   ! 5 -- Forest regrowth.                                                                 !
   ! 6 -- Logged forest.                                                                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_dist_types = 6
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Number of mortality types:                                                            !
   !                                                                                       !
   ! 1. Ageing, PFT-dependent but otherwise constant;                                      !
   ! 2. Negative carbon balance;                                                           !
   ! 3. Treefall mortality;                                                                !
   ! 4. Mortality due to cold weather.                                                     !
   ! 5. Disturbance mortality.  This is not directly applied to the cohort population,     !
   !    because this mortality is associated with the creation of a new patch, but it is   !
   !    saved here for posterior analysis.                                                 !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_mort = 5
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Number of radiation profiles.                                                       !
   !  1. PAR, Beam, Down                                                                   !
   !  2. PAR, Beam, Up (only Medvigy's radiation, zero for others).                        !
   !  3. PAR, Diff, Down                                                                   !
   !  4. PAR, Diff, Up                                                                     !
   !  5. NIR, Beam, Down                                                                   !
   !  6. NIR, Beam, Up (only Medvigy's radiation, zero for others).                        !
   !  7. NIR, Diff, Down                                                                   !
   !  8. NIR, Diff, Up                                                                     !
   !  9. TIR, Diff, Down                                                                   !
   ! 10. TIR. Diff, Up                                                                     !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: n_radprof = 10
   !---------------------------------------------------------------------------------------!



   !----- Maximum number meteorological driver xfiles. ------------------------------------!
   integer, parameter :: max_met_vars = 22
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      For restart runs, this is the maximum number of certain variables that can be    !
   ! read.                                                                                 !
   !  HUGE_POLYGON - maximum number of input polygons.                                     !
   !  HUGE_PATCH   - maximum number of input patches.                                      !
   !  HUGE_COHORT  - maximum number of input cohorts.                                      !
   !  MAX_WATER    - maximum number of soil water levels (not assigned to polygons).       !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: huge_polygon = nxpmax * nypmax
   integer, parameter :: huge_patch   = 10000
   integer, parameter :: huge_cohort  = 250000
   integer, parameter :: max_water    = 100
   !---------------------------------------------------------------------------------------!



   !----- Maximum number of land use polygons that can be read by filelist. ---------------!
   integer, parameter :: huge_lu = 99999
   !---------------------------------------------------------------------------------------!



   !----- The maximum number of printable variables. --------------------------------------!
   integer, parameter :: maxpvars = 50
   !---------------------------------------------------------------------------------------!



   !----- Maximum number of files that can be read by filelist. ---------------------------!
#if defined(COUPLED)
   integer, parameter :: maxfiles = brams_maxfiles
#else
   integer, parameter :: maxfiles = 99999
#endif
   !---------------------------------------------------------------------------------------!



   !----- Maximum number of files (site+patch+cohort). ------------------------------------!
   integer, parameter :: maxlist = 3 * maxfiles
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Although these aren't maximum dimensions, these are used to initialise the name-  !
   ! list variables.                                                                       !
   !---------------------------------------------------------------------------------------!
   real                  , parameter :: undef_real      = -999.9
   real(kind=8)          , parameter :: undef_dble      = -9.999d2
   integer               , parameter :: undef_integer   = -999
   character(len=str_len), parameter :: undef_character = 'nothing'
   character(len=str_len), parameter :: undef_path      = '/nowhere'
   logical               , parameter :: undef_logical   = .false.
   !---------------------------------------------------------------------------------------!

end module ed_max_dims
!==========================================================================================!
!==========================================================================================!
