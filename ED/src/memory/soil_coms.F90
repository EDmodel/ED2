!==========================================================================================!
!==========================================================================================!
!   This module contains soil- and surface-related properties to be used in various points !
! of the code.                                                                             !
!                                                                                          !
! IMPORTANT: DO NOT INITIALIZE PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL        !
!            ACTUALLY INITIALIZE THEM.  See "initialize_soil_coms" for some default        !
!            values.                                                                       !
!==========================================================================================!
!==========================================================================================!
module soil_coms
   use ed_max_dims  , only : str_len   & ! intent(in)
                           , maxgrds   & ! intent(in)
                           , nzgmax    ! ! intent(in)
   use grid_coms    , only : nzg       & ! intent(in)
                           , nzs       ! ! intent(in)
#if defined(COUPLED)
   use leaf_coms    , only : nstyp     & ! intent(in)
                           , nscol     & ! intent(in)
                           , nvtyp     & ! intent(in)
                           , nvtyp_teb ! ! intent(in)
#endif

   implicit none

  !----- These variables depend on whether it's a coupled or stand alone model. -----------!
#if defined(COUPLED)
   integer, parameter :: ed_nstyp = nstyp          ! total # of soil textural classes
   integer, parameter :: ed_nscol = nscol          ! total # of soil colour classes
   integer, parameter :: ed_nvtyp = nvtyp+nvtyp_teb
#else
   integer, parameter :: ed_nstyp = 17             ! total # of soil textural classes
   integer, parameter :: ed_nscol = 21             ! total # of soil colour classes
   integer, parameter :: ed_nvtyp = 21
#endif

   !=======================================================================================!
   !=======================================================================================!
   !    The following variables are assigned through the namelist.                         !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   ! ISOILBC -- This controls the soil moisture boundary condition at the bottom.  Choose  !
   !            the option according to the site characteristics.                          !
   !            0.  Flat bedrock.  Flux from the bottom of the bottommost layer is zero.   !
   !            1.  Gravitational flow (free drainage).  The flux from the bottom of the   !
   !                bottommost layer is due to gradient of height only.                    !
   !            2.  Lateral drainage.  Similar to free drainage, but the gradient is       !
   !                reduced by the slope not being completely vertical.  The reduction is  !
   !                controlled by variable SLDRAIN.  In the future options 0, 1, and 2 may !
   !                be combined into a single option.                                      !
   !            3.  Aquifer.  Soil moisture of the ficticious layer beneath the bottom is  !
   !                always at saturation.                                                  !
   !---------------------------------------------------------------------------------------!
   integer                                    :: isoilbc
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! SLDRAIN      -- This is used only when ISOILBC is set to 2.  In this case SLDRAIN is  !
   !                 the equivalent slope that will slow down drainage.  If this is set to !
   !                 zero, then lateral drainage reduces to flat bedrock, and if this is   !
   !                 set to 90, then lateral drainage becomes free drainage.  SLDRAIN must !
   !                 be between 0 and 90.                                                  !
   !                                                                                       !
   ! SLDRAIN8     -- These are auxiliary variables derived from sldrain.  SIN_SLDRAIN is   !
   ! SIN_SLDRAIN     the sine of the slope, and SLDRAIN8 and SIN_SLDRAIN8 are the double   !
   ! SIN_SLDRAIN8    precision versions of SLDRAIN and SIN_SLDRAIN.                        !
   !---------------------------------------------------------------------------------------!
   real                                       :: sldrain
   real(kind=8)                               :: sldrain8
   real                                       :: sin_sldrain
   real(kind=8)                               :: sin_sldrain8
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! ISOILFLG -- this controls which soil type input you want to use.                      !
   !             1. Read in from a dataset I will provide in the SOIL_DATABASE variable a  !
   !                few lines below.                                                       !
   !                  below.                                                               !
   !             2. No data available, I will use constant values I will provide in        !
   !                NSLCON or by prescribing the fraction of sand and clay (see SLXSAND    !
   !                and SLXCLAY).                                                          !
   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds)                :: isoilflg
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! NSLCON -- ED-2 Soil classes that the model will use when ISOILFLG is set to 2.        !
   !           Possible values are:                                                        !
   !---------------------------------------------------------------------------------------!
   !   1 -- sand                |   7 -- silty clay loam     |  13 -- bedrock              !
   !   2 -- loamy sand          |   8 -- clayey loam         |  14 -- silt                 !
   !   3 -- sandy loam          |   9 -- sandy clay          |  15 -- heavy clay           !
   !   4 -- silt loam           |  10 -- silty clay          |  16 -- clayey sand          !
   !   5 -- loam                |  11 -- clay                |  17 -- clayey silt          !
   !   6 -- sandy clay loam     |  12 -- peat                                              !
   !---------------------------------------------------------------------------------------!
   integer                                    :: nslcon
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! ISOILCOL -- LEAF-3 and ED-2 soil colour classes that the model will use when ISOILFLG !
   !             is set to 2.  Soil classes are from 1 to 20 (1 = lightest; 20 = darkest). !
   !             The values are the same as CLM-4.0.  The table is the albedo for visible  !
   !             and near infra-red.                                                       !
   !---------------------------------------------------------------------------------------!
   !                                                                                       !
   !       |-----------------------------------------------------------------------|       !
   !       |       |   Dry soil  |  Saturated  |       |   Dry soil  |  Saturated  |       !
   !       | Class |-------------+-------------| Class +-------------+-------------|       !
   !       |       |  VIS |  NIR |  VIS |  NIR |       |  VIS |  NIR |  VIS |  NIR |       !
   !       |-------+------+------+------+------+-------+------+------+------+------|       !
   !       |     1 | 0.36 | 0.61 | 0.25 | 0.50 |    11 | 0.24 | 0.37 | 0.13 | 0.26 |       !
   !       |     2 | 0.34 | 0.57 | 0.23 | 0.46 |    12 | 0.23 | 0.35 | 0.12 | 0.24 |       !
   !       |     3 | 0.32 | 0.53 | 0.21 | 0.42 |    13 | 0.22 | 0.33 | 0.11 | 0.22 |       !
   !       |     4 | 0.31 | 0.51 | 0.20 | 0.40 |    14 | 0.20 | 0.31 | 0.10 | 0.20 |       !
   !       |     5 | 0.30 | 0.49 | 0.19 | 0.38 |    15 | 0.18 | 0.29 | 0.09 | 0.18 |       !
   !       |     6 | 0.29 | 0.48 | 0.18 | 0.36 |    16 | 0.16 | 0.27 | 0.08 | 0.16 |       !
   !       |     7 | 0.28 | 0.45 | 0.17 | 0.34 |    17 | 0.14 | 0.25 | 0.07 | 0.14 |       !
   !       |     8 | 0.27 | 0.43 | 0.16 | 0.32 |    18 | 0.12 | 0.23 | 0.06 | 0.12 |       !
   !       |     9 | 0.26 | 0.41 | 0.15 | 0.30 |    19 | 0.10 | 0.21 | 0.05 | 0.10 |       !
   !       |    10 | 0.25 | 0.39 | 0.14 | 0.28 |    20 | 0.08 | 0.16 | 0.04 | 0.08 |       !
   !       |-----------------------------------------------------------------------|       !
   !                                                                                       !
   !   Soil type 21 is a special case in which we use the albedo method that used to be    !
   ! the default in ED-2.1.                                                                !
   !---------------------------------------------------------------------------------------!
   integer                                    :: isoilcol
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     These variables are used to define the soil properties when you don't want to use !
   ! the standard soil classes.                                                            !
   !                                                                                       !
   ! SLXCLAY -- Prescribed fraction of clay  [0-1]                                         !
   ! SLXSAND -- Prescribed fraction of sand  [0-1].                                        !
   !                                                                                       !
   !     They are used only when ISOILFLG is 2, both values are between 0. and 1., and     !
   ! theira sum doesn't exceed 1.  Otherwise standard ED values will be used instead.      !
   !---------------------------------------------------------------------------------------!
   real                                       :: slxclay
   real                                       :: slxsand
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !  ZROUGH -- constant roughness, in metres, if for all domain                           !
   !---------------------------------------------------------------------------------------!
   real                                       :: zrough
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Soil grid and initial conditions if no file is provided:                           !
   !                                                                                       !
   ! SLZ     - soil depth in m.  Values must be negative and go from the deepest layer to  !
   !           the top.                                                                    !
   ! SLMSTR  - this is the initial soil moisture, now given as the soil moisture index.    !
   !           Values can be fraction, in which case they will be linearly interpolated    !
   !           between the special points (e.g. 0.5 will put soil moisture half way        !
   !           between the wilting point and field capacity).                              !
   !              -1 = dry air soil moisture                                               !
   !               0 = wilting point                                                       !
   !               1 = field capacity                                                      !
   !               2 = porosity (saturation)                                               !
   ! STGOFF  - initial temperature offset (soil temperature = air temperature + offset)    !
   !---------------------------------------------------------------------------------------!
   real, dimension(nzgmax)                    :: slz
   real, dimension(nzgmax)                    :: slmstr
   real, dimension(nzgmax)                    :: stgoff
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !  Input databases                                                                      !
   !  VEG_DATABASE     -- vegetation database, used only to determine the land/water mask. !
   !                      Fill with the path and the prefix.                               !
   !  SOIL_DATABASE    -- soil database, used to determine the soil type.  Fill with the   !
   !                      path and the prefix.                                             !
   !  SOILSTATE_DB     -- Dataset in case you want to provide the initial conditions of    !
   !                      soil temperature and moisture.                                   !
   !  SOILDEPTH_DB     -- Dataset in case you want to read in soil depth information.      !
   !---------------------------------------------------------------------------------------!
   character(len=str_len), dimension(maxgrds) :: veg_database
   character(len=str_len), dimension(maxgrds) :: soil_database
   character(len=str_len)                     :: soilstate_db
   character(len=str_len)                     :: soildepth_db
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! ISOILSTATEINIT -- Variable controlling how to initialise the soil temperature and     !
   !                   moisture                                                            !
   !                   0.  Use SLMSTR and STGOFF.                                          !
   !                   1.  Read from SOILSTATE_DB.                                         !
   ! ISOILDEPTHFLG  -- Variable controlling how to initialise soil depth                   !
   !                   0.  Constant, always defined by the first SLZ layer.                !
   !                   1.  Read from SOILDEPTH_DB.                                         !
   !---------------------------------------------------------------------------------------!
   integer                                    :: isoilstateinit ! Soil state initial cond. 
   integer                                    :: isoildepthflg  ! Soil depth initial cond. 
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   ! RUNOFF_TIME -- In case a temporary surface water (TSW) is created, this is the "e-    !
   !                -folding lifetime" of the TSW in seconds due to runoff.  If you don't  !
   !                want runoff to happen, set this to 0.                                  !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   real                                       :: runoff_time
   real                                       :: runoff_time_i
   logical                                    :: simplerunoff
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!


   !---------------------------------------------------------------------------------------!
   !    These following variables will behave as parameters, but they are initialized at   !
   ! init_soil_coms (ed_params.f90), or through the XML config file.                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: soil_rough          ! soil roughness height                   [       m]
   real(kind=4) :: snow_rough          ! snowcover roughness height              [       m]
   real(kind=4) :: ny07_eq04_a         ! parameter a for snow fraction           [     ---]
   real(kind=4) :: ny07_eq04_m         ! parameter m for snow fraction           [     ---]
   real(kind=8) :: soil_rough8         ! soil roughness height                   [       m]
   real(kind=8) :: snow_rough8         ! snowcover roughness height              [       m]
   real(kind=8) :: ny07_eq04_a8        ! parameter a for snow fraction           [     ---]
   real(kind=8) :: ny07_eq04_m8        ! parameter m for snow fraction           [     ---]
   real(kind=4) :: dewmax              ! Maximum dew flux rate (deprecated)      [ kg/m2/s]
   real(kind=4) :: water_stab_thresh   ! stability threshold for RK4 integrator  [   kg/m2]
   real(kind=4) :: tiny_sfcwater_mass  ! Min. mass allowed in temporary layers   [   kg/m2]
   real(kind=4) :: snowmin             ! Min. snow mass needed to create new lyr [   kg/m2]
   integer      :: infiltration_method ! Infiltration scheme (for rk4_derivs)    [     0|1]
   real(kind=4) :: freezecoef          ! Coeff. for infiltration of frozen water [     ---]
   real(kind=8) :: freezecoef8         ! Coeff. for infiltration of frozen water [     ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Miscellaneous constants.                                                          !
   !---------------------------------------------------------------------------------------!
   integer           , parameter :: pctlcon = 1
   integer           , parameter :: nvgcon = 7 ! I don't think it is been used...
   !----- Constants from  equation E27 (Medvigy 2007) -------------------------------------!
   real(kind=8), dimension(6), parameter :: ss = (/ 1.093d-3, 2.800d-2, 3.000d-2           &
                                                  , 3.030d-4,-1.770d-7, 2.250d-9 /) 
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Variables to be initialised in sfcdata_ed (ed_init.f90).                          !
   !---------------------------------------------------------------------------------------!
   real   , dimension(ed_nstyp)        :: slden    ! dry soil density               [kg/m3]
   real   , dimension(ed_nstyp)        :: fhydraul ! vertically varying hydraulic 
                                                   !    conductivity factor
   integer, dimension(ed_nvtyp)        :: kroot    ! level in which roots are       [  ---]
   real   , dimension(nzgmax,ed_nvtyp) :: root     ! root depth                     [    m]
   real, allocatable, dimension(:,:)   :: slcons1  ! z-dep. soil sat hydraul cond   [  m/s]
   real, allocatable, dimension(:)     :: dslz     ! soil layer thickness at T pt   [    m]
   real, allocatable, dimension(:)     :: dslzo2   ! ½ soil layer thick. at T pt    [    m]
   real, allocatable, dimension(:)     :: dslzi    ! 1/dslz                         [  1/m]
   real, allocatable, dimension(:)     :: dslzidt  ! dtll / dslz                    [  s/m]
   real, allocatable, dimension(:)     :: slzt     ! soil depth at T pt             [    m]
   real, allocatable, dimension(:)     :: dslzt    ! soil layer thickness at M pt   [    m]
   real, allocatable, dimension(:)     :: dslzti   ! 1/dslzt                        [  1/m]
   real, allocatable, dimension(:)     :: dslztidt ! dtll / dslzt                   [  s/m]
   !----- The next variables are double precision version of the previous ones. -----------!
   real(kind=8), allocatable, dimension(:)   :: slz8            ! Soil levels.
   real(kind=8), allocatable, dimension(:,:) :: slcons18
   real(kind=8), allocatable, dimension(:)   :: dslz8
   real(kind=8), allocatable, dimension(:)   :: dslzo28
   real(kind=8), allocatable, dimension(:)   :: dslzi8
   real(kind=8), allocatable, dimension(:)   :: dslzidt8
   real(kind=8), allocatable, dimension(:)   :: slzt8
   real(kind=8), allocatable, dimension(:)   :: dslzt8
   real(kind=8), allocatable, dimension(:)   :: dslzti8
   real(kind=8), allocatable, dimension(:)   :: dslztidt8
   real(kind=8), dimension(20)               :: thicknet
   real(kind=8), dimension(20,20)            :: thick   
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Scratch matrix for reading the soil layer and soil depth dataset.  This will be    !
   ! deallocated before the first time step.  The number of points is set here for a       !
   ! global dataset with 1x1 degree resolution.                                            !
   !---------------------------------------------------------------------------------------!
   integer, parameter                   :: nlon_lyr=360
   integer, parameter                   :: nlat_lyr=180
   integer, allocatable, dimension(:,:) :: layer_index
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Soil Characteristics. See: Clapp & Hornberger, 1978;                                 !
   !                             McCumber & Pielke, 1981;                                  !
   !                             Pielke, 1984;                                             !
   !                             Tremback & Kessler, 1985).                                !
   !---------------------------------------------------------------------------------------!
   type soil_class
      real(kind=4) :: slpots     ! Soil moisture potential at saturation         [       m]
      real(kind=4) :: slmsts     ! Soil moisture at saturation                   [   m3/m3]
      real(kind=4) :: slbs       ! B exponent                                    [     n/d]
      real(kind=4) :: slcpd      ! Specific heat of dry soil                     [  J/m3/K]
      real(kind=4) :: soilcp     ! Dry soil capacity (at -3.1MPa)                [   m3/m3]
      real(kind=4) :: soilwp     ! Wilting point capacity (at -1.5MPa)           [   m3/m3]
      real(kind=4) :: slcons     ! hydraulic conductivity at saturation          [     m/s]
      real(kind=4) :: slcons0    ! Surface value for slcons                      [     m/s]
      real(kind=4) :: thcond0    ! First coefficient for thermal conductivity    [   W/m/K]
      real(kind=4) :: thcond1    ! Second coefficient for thermal conductivity   [   W/m/K]
      real(kind=4) :: thcond2    ! Third coefficient for thermal conductivity    [     ---]
      real(kind=4) :: thcond3    ! Fourth coefficient for thermal conductivity   [     ---]
      real(kind=4) :: sfldcap    ! Soil field capacity                           [   m3/m3]
      real(kind=4) :: xsand      ! Percentage of sand                            [     ---]
      real(kind=4) :: xclay      ! Percentage of clay                            [     ---]
      real(kind=4) :: xsilt      ! Percentage of silt                            [     ---]
      real(kind=4) :: xrobulk    ! Bulk density                                  [     ---]
      real(kind=4) :: slden      ! "Dry" soil density (porosity)                 [   kg/m3]
      real(kind=4) :: soilld     ! Soil moist. below which drought phen. happens [   m3/m3]
      real(kind=4) :: soilfr     ! Soil moist. below which fires may happen      [   m3/m3]
      real(kind=4) :: slpotwp    ! Water potential for wilting point             [       m]
      real(kind=4) :: slpotfc    ! Water potential for field capacity            [       m]
      real(kind=4) :: slpotld    ! Water pot. below which drought phen happens   [       m]
      real(kind=4) :: slpotfr    ! Water pot. below which fire happens           [       m]
   end type soil_class
   !----- Double precision version --------------------------------------------------------!
   type soil_class8
      real(kind=8) :: slpots     ! Soil moisture potential at saturation         [       m]
      real(kind=8) :: slmsts     ! Soil moisture at saturation                   [   m3/m3]
      real(kind=8) :: slbs       ! B exponent                                    [     n/d]
      real(kind=8) :: slcpd      ! Specific heat of dry soil                     [  J/m3/K]
      real(kind=8) :: soilcp     ! Dry soil capacity (at -3.1MPa)                [   m3/m3]
      real(kind=8) :: soilwp     ! Wilting point capacity (at -1.5MPa)           [   m3/m3]
      real(kind=8) :: slcons     ! hydraulic conductivity at saturation          [     m/s]
      real(kind=8) :: slcons0    ! Surface value for slcons                      [     m/s]
      real(kind=8) :: thcond0    ! First coefficient for thermal conductivity    [   W/m/K]
      real(kind=8) :: thcond1    ! Second coefficient for thermal conductivity   [   W/m/K]
      real(kind=8) :: thcond2    ! Third coefficient for thermal conductivity    [     ---]
      real(kind=8) :: thcond3    ! Fourth coefficient for thermal conductivity   [     ---]
      real(kind=8) :: sfldcap    ! Soil field capacity                           [   m3/m3]
      real(kind=8) :: xsand      ! Percentage of sand                            [     ---]
      real(kind=8) :: xclay      ! Percentage of clay                            [     ---]
      real(kind=8) :: xsilt      ! Percentage of silt                            [     ---]
      real(kind=8) :: xrobulk    ! Bulk density                                  [     ---]
      real(kind=8) :: slden      ! "Dry" soil density (porosity)                 [   kg/m3]
      real(kind=8) :: soilld     ! Soil moist. below which drought phen. happens [   m3/m3]
      real(kind=8) :: soilfr     ! Soil moist. below which fires may happen      [   m3/m3]
      real(kind=8) :: slpotwp    ! Water potential for wilting point             [       m]
      real(kind=8) :: slpotfc    ! Water potential for field capacity            [       m]
      real(kind=8) :: slpotld    ! Water pot. below which drought phen happens   [       m]
      real(kind=8) :: slpotfr    ! Water pot. below which fire happens           [       m]
   end type soil_class8
   !---------------------------------------------------------------------------------------!
   !----- To be filled in ed_params.f90. --------------------------------------------------!
   type(soil_class8), dimension(ed_nstyp)            :: soil8 
   type(soil_class) , dimension(ed_nstyp)            :: soil
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Define soil colour structure.                                                     !
   !---------------------------------------------------------------------------------------!
   type soilcol_class
      real(kind=4) :: alb_vis_dry
      real(kind=4) :: alb_nir_dry
      real(kind=4) :: alb_vis_wet
      real(kind=4) :: alb_nir_wet
      real(kind=4) :: emiss_tir
   end type soilcol_class
   !----- To be filled in ed_params.f90. --------------------------------------------------!
   type(soilcol_class), dimension(ed_nscol) :: soilcol
   !---------------------------------------------------------------------------------------!


 ! Look-up tables for vegetation and soil properties:
   type veg_class
     real :: albedv
     real :: emisv
     real :: vglai
     real :: vgdlai
     real :: vgfrac
     real :: vgdfrac
     real :: vegzo
     real :: vgdisp
     real :: rootdep
   end type veg_class

   type(veg_class), parameter, dimension(0:ed_nvtyp) ::  veget=(/  &
         !-----------------------------------------------------------------------------
         !albedo    lai       vfrac       zo      rootdep    LEAF-3 CLASS #
         !     emiss     dlai     dvfrac      zdisp          AND DESCRIPTION
         !-----------------------------------------------------------------------------
            veg_class(.14, .99, 0.0, 0.0, .00, .00,  .00,  0.1,  .0 )  & !  0  Ocean
         ,  veg_class(.14, .99, 0.0, 0.0, .00, .00,  .00,  0.1,  .0 )  & !  1  Lakes, rivers, streams (inland water)
         ,  veg_class(.40, .82, 0.0, 0.0, .00, .00,  .01,  0.1,  .0 )  & !  2  Ice cap/glacier
         ,  veg_class(.30, .86, 0.0, 0.0, .00, .00,  .05,   .1, 1.0 )  & !  3  Desert
         ,  veg_class(.10, .97, 6.0, 1.0, .80, .10, 1.00, 15.0, 1.5 )  & !  4  Evergreen needleleaf tree
         ,  veg_class(.10, .95, 6.0, 5.0, .80, .30, 1.00, 20.0, 1.5 )  & !  5  Deciduous needleleaf tree
         ,  veg_class(.20, .95, 6.0, 5.0, .80, .30,  .80, 15.0, 2.0 )  & !  6  Deciduous broadleaf tree
         ,  veg_class(.15, .95, 6.0, 1.0, .90, .50, 2.00, 20.0, 3.0 )  & !  7  Evergreen broadleaf tree
         ,  veg_class(.26, .96, 2.0, 1.5, .80, .10,  .02,   .2, 1.0 )  & !  8  Short grass
         ,  veg_class(.16, .96, 6.0, 5.5, .80, .30,  .10,  1.0, 1.0 )  & !  9  Tall grass
         ,  veg_class(.25, .96, 6.0, 5.5, .10, .10,  .10,   .5, 1.0 )  & ! 10  Semi-desert
         ,  veg_class(.20, .95, 6.0, 5.5, .60, .20,  .04,   .1, 1.0 )  & ! 11  Tundra
         ,  veg_class(.10, .97, 6.0, 1.0, .80, .20,  .10,  1.0, 1.0 )  & ! 12  Evergreen shrub
         ,  veg_class(.20, .97, 6.0, 5.0, .80, .30,  .10,  1.0, 1.0 )  & ! 13  Deciduous shrub
         ,  veg_class(.15, .96, 6.0, 3.0, .80, .20,  .80, 20.0, 2.0 )  & ! 14  Mixed woodland
         ,  veg_class(.20, .95, 6.0, 5.5, .85, .60,  .06,   .7, 1.0 )  & ! 15  Crop/mixed farming
         ,  veg_class(.18, .95, 6.0, 5.5, .80, .60,  .06,   .7, 1.0 )  & ! 16  Irrigated crop
         ,  veg_class(.12, .98, 6.0, 5.5, .80, .40,  .03,  1.0, 1.0 )  & ! 17  Bog or marsh
         ,  veg_class(.18, .96, 5.0, 4.0, .80, .20,  .51,  3.6, 1.0 )  & ! 18  Wooded grassland
         ,  veg_class(.15, .90, 4.8, 3.6, .74, .31,  .80,  1.1,  .8 )  & ! 19  Urban and built up
         ,  veg_class(.21, .95, 6.0, 0.8, .90, .50, 2.00, 20.0, 3.0 )  & ! 20  Wetland evergreen broadleaf tree
         ,  veg_class(.15, .96, 4.8, 3.6, .10, .05, 1.40,  1.1,  .8 )   /) ! 21  Very urban
 !MLO]

   !=======================================================================================!
   !=======================================================================================!

   contains
   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine allocates soil grid arrays.                                        !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_soilgrid()
      implicit none


      allocate (slcons1(0:nzg,ed_nstyp))

      allocate (dslz     (0:nzg))
      allocate (dslzo2   (0:nzg))
      allocate (dslzi    (0:nzg))
      allocate (dslzidt  (0:nzg))
      allocate (slzt     (0:nzg))
      allocate (dslzt    (  nzg))
      allocate (dslzti   (  nzg))
      allocate (dslztidt (  nzg))

      allocate (slcons18(0:nzg,ed_nstyp))

      allocate (slz8      (nzg+1))
      allocate (dslz8     (0:nzg))
      allocate (dslzo28   (0:nzg))
      allocate (dslzi8    (0:nzg))
      allocate (dslzidt8  (0:nzg))
      allocate (slzt8     (0:nzg))
      allocate (dslzt8    (  nzg))
      allocate (dslzti8   (  nzg))
      allocate (dslztidt8 (  nzg))

      return
   end subroutine alloc_soilgrid
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function determines the soil class based on the fraction of sand, clay, and  !
   ! silt separates.                                                                       !
   !---------------------------------------------------------------------------------------!
   integer function find_soil_class(sandfrac,clayfrac)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: sandfrac
      real   , intent(in) :: clayfrac
      !----- Local variables. -------------------------------------------------------------!
      integer             :: sclass
      real                :: sand
      real                :: clay
      real                :: silt
      !------------------------------------------------------------------------------------!



      !----- Define the percentage of sand, clay, and silt. -------------------------------!
      sand = 100. * sandfrac
      clay = 100. * clayfrac
      silt = 100. - sand - clay
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Here there is not much we can do other than explore where in the triangle      !
      ! space we are.                                                                      !
      !------------------------------------------------------------------------------------!
      if (silt > 100. .or. silt < 0. .or. sand > 100. .or. sand < 0. .or.                  &
          clay > 100. .or. clay < 0. ) then
         write (unit=*,fmt='(a)') '---------------------------------------------------'
         write (unit=*,fmt='(a)') ' At least one of your percentages is off-bounds...'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'SAND = ',sand,'%'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'CLAY = ',clay,'%'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'SILT = ',silt,'%'
         write (unit=*,fmt='(a)') ' This soil doesn''t fit into any category...'
         write (unit=*,fmt='(a)') '---------------------------------------------------'
         call fatal_error('Incorrect fractions of sand, clay, or silt...'                  &
                         ,'find_soil_class','soil_coms.F90')
      elseif (sand > 85.0 + 0.5 * clay) then
         sclass =  1  !----- Sand. --------------------------------------------------------!
      elseif (sand > 70.0 + clay) then
         sclass =  2  !----- Loamy sand. --------------------------------------------------!
      elseif ((clay <= 20.0 .and. sand > 52.5) .or. (clay <= 7.5 .and. silt <= 50.0)) then
         sclass =  3  !----- Sandy loam. --------------------------------------------------!
      elseif ((clay <= 27.5 .and. silt > 50.0 .and. silt <= 80.0) .or.                     &
              (silt >  80.0 .and. clay > 12.5)                    ) then
         sclass =  4  !----- Silt loam. ---------------------------------------------------!
      elseif (clay > 7.5 .and. clay <= 27.5 .and. silt > 27.5 .and. silt <= 50.0 .and.     &
              sand <= 52.5) then
         sclass =  5  !----- Loam. --------------------------------------------------------!
      elseif (clay > 20.0 .and. clay <= 35.0 .and. silt <= 27.5 .and. sand > 45.0) then
         sclass =  6  !----- Sandy clay loam. ---------------------------------------------!
      elseif (clay > 27.5 .and. clay <= 40.0 .and. sand <= 20.0) then
         sclass =  7  !----- Silty clay loam. ---------------------------------------------!
      elseif (clay > 27.5 .and. clay <= 40.0 .and. sand > 20.0 .and. sand <= 45.0) then
         sclass =  8  !----- Clayey loam. -------------------------------------------------!
      elseif (clay > 35.0 .and. sand > 45.0) then
         sclass =  9  !----- Sandy clay. --------------------------------------------------!
      elseif (clay > 40.0 .and. silt > 40.0) then
         sclass = 10  !----- Silty clay. --------------------------------------------------!
      elseif (clay <= 70.0 .and. sand <= 30.0 .and. silt <= 30.0) then
         sclass = 11  !----- Clay. --------------------------------------------------------!
      elseif ( silt > 80.0 .and. clay <= 12.5) then
         sclass = 14  !----- Silt. --------------------------------------------------------!
      elseif ( clay > 70.0) then
         sclass = 15  !----- Heavy clay. --------------------------------------------------!
      elseif ( clay > 40.0 .and. sand > 30.0 .and. sand <= 45.0) then
         sclass = 16  !----- Clayey sand. -------------------------------------------------!
      elseif ( clay > 40.0 .and. silt > 30.0 .and. silt <= 40.0) then
         sclass = 17  !----- Clayey silt. -------------------------------------------------!
      else
         write (unit=*,fmt='(a)') '---------------------------------------------------'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'SAND = ',sand,'%'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'CLAY = ',clay,'%'
         write (unit=*,fmt='(a,1x,f8.2,a)') 'SILT = ',silt,'%'
         write (unit=*,fmt='(a)') ' This soil doesn''t fit into any category...'
         write (unit=*,fmt='(a)') '---------------------------------------------------'
         call fatal_error('Failed finding the correct soil class...'                       &
                         ,'find_soil_class','soil_coms.F90')
      end if
      !------------------------------------------------------------------------------------!

      find_soil_class = sclass

      return
   end function find_soil_class
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routines converts the soil moisture index into soil moisture.           !
   !  This scale is piece-wise linear, and dependent on the soil texture.                  !
   !  -1. = dry air soil moisture (soilcp)                                                 !
   !   0. = wilting point         (soilwp)                                                 !
   !   1. = field capacity        (sfldcap)                                                !
   !   2. = porosity/saturation   (slmsts)                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function ed_soil_idx2water(soil_index,ntext)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: soil_index
      integer     , intent(in) :: ntext
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: soil_water
      !------------------------------------------------------------------------------------!


      if (soil_index < 0.0) then
         soil_water = soil(ntext)%soilcp                                                   &
                    + (soil_index + 1.0) * (soil(ntext)%soilwp  - soil(ntext)%soilcp )
      elseif (soil_index < 1.0) then
         soil_water = soil(ntext)%soilwp                                                   &
                    +  soil_index        * (soil(ntext)%sfldcap - soil(ntext)%soilwp )
      else
         soil_water = soil(ntext)%sfldcap                                                  &
                    + (soil_index - 1.0) * (soil(ntext)%slmsts  - soil(ntext)%sfldcap)
      end if

      ed_soil_idx2water = soil_water

      return
   end function ed_soil_idx2water
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil moisture to soil matric potential.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function matric_potential(nsoil,soil_water)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: nsoil      ! Soil texture                         [  idx]
      real(kind=4), intent(in) :: soil_water ! Soil moisture                        [m3/m3]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist   ! Relative soil moisture               [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist      = min(soil_water/soil(nsoil)%slmsts,1.0)
      !------------------------------------------------------------------------------------!



      !----- Find the matric potential. ---------------------------------------------------!
      matric_potential = soil(nsoil)%slpots / relmoist ** soil(nsoil)%slbs
      !------------------------------------------------------------------------------------!

      return
   end function matric_potential
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function is the double precision version of the function above, so it also  !
   ! converts soil moisture to soil matric potential.                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function matric_potential8(nsoil,soil_water)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: nsoil      ! Soil texture                         [  idx]
      real(kind=8), intent(in) :: soil_water ! Soil moisture                        [m3/m3]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=8)             :: relmoist   ! Relative soil moisture               [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist      = min(soil_water/soil8(nsoil)%slmsts,1.d0)
      !------------------------------------------------------------------------------------!



      !----- Find the matric potential. ---------------------------------------------------!
      matric_potential8 = soil8(nsoil)%slpots / relmoist ** soil8(nsoil)%slbs
      !------------------------------------------------------------------------------------!

      return
   end function matric_potential8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil moisture to hydraulic conductivity.                  !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function hydr_conduct(k,nsoil,soil_water,soil_fracliq)
      use consts_coms, only : lnexp_min ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: k            ! Layer index                        [  idx]
      integer     , intent(in) :: nsoil        ! Soil texture                       [  idx]
      real(kind=4), intent(in) :: soil_water   ! Soil moisture                      [m3/m3]
      real(kind=4), intent(in) :: soil_fracliq ! Liquid fraction                    [  ---]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist     ! Relative soil moisture             [  ---]
      real(kind=4)             :: fzcorr       ! Freezing correction                [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find correction for frozen soils. -------------------------------------------!
      fzcorr = exp( max( lnexp_min, - freezecoef * (1.0 - soil_fracliq) ) )
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist = min(soil_water/soil(nsoil)%slmsts,1.0)
      !------------------------------------------------------------------------------------!



      !----- Find the hydraulic conductivity. ---------------------------------------------!
      hydr_conduct = fzcorr * slcons1(k,nsoil) * relmoist ** (2. * soil(nsoil)%slbs + 3.)
      !------------------------------------------------------------------------------------!


      return
   end function hydr_conduct
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function is the double precision version of the function above, so it also  !
   ! converts soil moisture to hydraulic conductivity.                                     !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function hydr_conduct8(k,nsoil,soil_water,soil_fracliq)
      use consts_coms, only : lnexp_min8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: k            ! Layer index                        [  idx]
      integer     , intent(in) :: nsoil        ! Soil texture                       [  idx]
      real(kind=8), intent(in) :: soil_water   ! Soil moisture                      [m3/m3]
      real(kind=8), intent(in) :: soil_fracliq ! Liquid fraction                    [  ---]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=8)             :: relmoist     ! Relative soil moisture             [  ---]
      real(kind=8)             :: fzcorr       ! Freezing correction                [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find correction for frozen soils. -------------------------------------------!
      fzcorr = exp( max( lnexp_min8, - freezecoef8 * (1.d0 - soil_fracliq) ) )
      !------------------------------------------------------------------------------------!



      !------ Find relative soil moisture. ------------------------------------------------!
      relmoist      = min(soil_water/soil8(nsoil)%slmsts,1.d0)
      !------------------------------------------------------------------------------------!



      !----- Find the hydraulic conductivity. ---------------------------------------------!
      hydr_conduct8 = fzcorr * slcons18(k,nsoil)                                           &
                    * relmoist ** (2.d0 * soil8(nsoil)%slbs + 3.d0)
      !------------------------------------------------------------------------------------!


      return
   end function hydr_conduct8
   !=======================================================================================!
   !=======================================================================================!
end module soil_coms
!==========================================================================================!
!==========================================================================================!
