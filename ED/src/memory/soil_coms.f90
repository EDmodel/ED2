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
                           , nzgmax    & ! intent(in)
                           , ed_nstyp  & ! intent(in)
                           , ed_nscol  & ! intent(in)
                           , ed_nvtyp  ! ! intent(in)
   use grid_coms    , only : nzg       & ! intent(in)
                           , nzs       ! ! intent(in)

   implicit none


   !=======================================================================================!
   !=======================================================================================!
   !    The following variables are assigned through the namelist.                         !
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! SOIL_HYDRO_SCHEME -- This controls the pedotransfer functions in ED-2.2.              !
   !                                                                                       !
   !    0. (ED-2.2 default).  Matric potential and hydraulic conductivity are determined   !
   !       using the Campbell(-Mualem) approach.  Pedotransfer function                    !
   !       parameters are determined using the Cosby et al. (1984) equations, which        !
   !       depend only on sand, silt, and clay fractions and are fitted for US soils.      !
   !    1. (Beta).  Matric potential and hydraulic conductivity are determined             !
   !       using the Campbell(-Mualem) approach, although not assuming that                !
   !       residual moisture is zero.  Pedotransfer function parameters are determined     !
   !       using the Tomasella and Hodnett (1998) equations, which depend only on          !
   !       sand, silt, and clay fractions and are fitted for tropical soils.               !
   !    2. (Beta).  Matric potential and hydraulic conductivity are determined using the   !
   !       van Genuchten(-Mualem) approach. Pedotransfer function parameters are           !
   !       determined using the Hodnett and Tomasella (2002) equations, which depend on    !
   !       sand, silt, and clay fractions, soil organic carbon content, pH, cation         !
   !       exchange capacity, and dry bulk density.                                        !
   !                                                                                       !
   !    The reference below provides details on all the approaches above:                  !
   !                                                                                       !
   !    Marthews TR, Quesada CA, Galbraith DR, Malhi Y, Mullins CE, Hodnett MG, Dharssi I, !
   !       2014. High-resolution hydraulic parameter maps for surface soils in tropical    !
   !       South America. Geosci. Model Dev. 7 (3), 711-723, doi:10.5194/gmd-7-711-2014.   !
   !---------------------------------------------------------------------------------------!
   integer :: soil_hydro_scheme
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
   ! ISLCOLFLG -- This controls how to initialise soil colour.  This must be a list with   !
   !              N_ED_REGION+N_POI elements.  The first N_ED_REGION elements correspond   !
   !              to each gridded domain (from first to last).  Elements between           !
   !              N_ED_REGION+1 and N_ED_REGION+N_POI correspond to the polygons of        !
   !              interest (from 1 to N_POI.  Options are:                                 !
   !             1 -- Read in soil colour class from the files set in COLOR_DATABASE.      !
   !             2 -- Assign either the value set by ISOILCOL (see below).                 !
   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds)                :: islcolflg
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
   ! ISOILCOL -- LEAF-3 and ED-2 soil colour classes that the model will use when          !
   !             ISLCOLFLG is set to 2.  Soil classes are from 1 to 20 (1 = lightest;      !
   !             20 = darkest).  The values are the same as CLM-4.0.  The table is the     !
   !             albedo for visible and near infra-red.                                    !
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
   !     These variables are used to define the additional soil properties needed for      !
   ! SOIL_HYDRO_SCHEME=2.  These should not be XML parameters because they describe        !
   ! site-specific soil characteristics, and the code will be eventually modified to allow !
   ! reading these characteristics from HDF5 files (similarly to soil depth, colour, and   !
   ! texture).  Eventually SOC effects on soil will be made dynamic and consistent with    !
   ! simulated soil carbon content (and thus eliminating the need for SLSOC).              !
   !                                                                                       !
   ! SLSOC   -- Prescribed mass fraction of soil organic carbon [ kg/kg].                  !
   ! SLPH    -- Prescribed soil potential of hydrogen (pH)      [  0-14].                  !
   ! SLCEC   -- Prescribed cation exchange capacity             [mol/kg].                  !
   ! SLDBD   -- Prescribed dry bulk density                     [ kg/m3].                  !
   !                                                                                       !
   !     They are used only when ISOILFLG is 2, both values are between 0. and 1., and     !
   ! their sum doesn't exceed 1.  In case ISOILFLG is 2 but the fractions do not meet the  !
   ! criteria, ED-2 uses NSLCON instead.                                                   !
   !---------------------------------------------------------------------------------------!
   real                                       :: slsoc
   real                                       :: slph
   real                                       :: slcec
   real                                       :: sldbd
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
   !  SLCOL_DATABASE   -- If ISLCOLFLG=1, this variable specifies the path and prefix of   !
   !                      soil colour data base.                                           !
   !  SOILSTATE_DB     -- Dataset in case you want to provide the initial conditions of    !
   !                      soil temperature and moisture.                                   !
   !  SOILDEPTH_DB     -- Dataset in case you want to read in soil depth information.      !
   !---------------------------------------------------------------------------------------!
   character(len=str_len), dimension(maxgrds) :: veg_database
   character(len=str_len), dimension(maxgrds) :: soil_database
   character(len=str_len), dimension(maxgrds) :: slcol_database
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
   real(kind=4) :: hydcond_min         ! Coeff. for infiltration of frozen water [     m/s]
   real(kind=8) :: hydcond_min8        ! Coeff. for infiltration of frozen water [     m/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Miscellaneous constants.                                                          !
   !---------------------------------------------------------------------------------------!
   integer           , parameter :: pctlcon = 1
   !----- Constants from  equation E27 (Medvigy 2007) -------------------------------------!
   real(kind=8), dimension(6), parameter :: ss = (/ 1.093d-3, 2.800d-2, 3.000d-2           &
                                                  , 3.030d-4,-1.770d-7, 2.250d-9 /) 
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Variables to be initialised in sfcdata_ed (ed_init.f90).                          !
   !---------------------------------------------------------------------------------------!
   integer, dimension(ed_nvtyp)        :: kroot    ! level in which roots are       [  ---]
   real   , dimension(nzgmax,ed_nvtyp) :: root     ! root depth                     [    m]
   real, allocatable, dimension(:,:)   :: slcons1  ! z-dep. soil sat hydraul cond   [  m/s]
   real, allocatable, dimension(:)     :: dslz     ! soil layer thickness at T pt   [    m]
   real, allocatable, dimension(:)     :: dslzo2   ! Half soil layer thick. at T pt [    m]
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
   !    Soil Characteristics.  For a general overview, check (M14).  Additional references !
   ! correspond to specific parametrisations.                                              !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Brooks RH , Corey AT. 1964. Hydraulic properties of porous media. Hydrology Papers 3, !
   !    Colorado State University, Fort Collins, U.S.A (BC64).                             !
   ! Marthews TR, Quesada CA, Galbraith DR, Malhi Y, Mullins CE, Hodnett MG , Dharssi I.   !
   !    2014. High-resolution hydraulic parameter maps for surface soils in tropical South !
   !    America. Geosci. Model Dev. 7: 711-723. doi:10.5194/gmd-7-711-2014 (M14).          !
   ! Campbell GS. 1974. A simple method for determining unsaturated conductivity from      !
   !    moisture retention data. Soil Science 117: 311-314.                                !
   !    doi:10.1097/00010694-197406000-00001 (C74).                                        !
   ! Cosby BJ, Hornberger GM, Clapp RB , Ginn TR. 1984. A statistical exploration of the   !
   !    relationships of soil moisture characteristics to the physical properties of       !
   !    soils. Water Resour. Res. 20: 682-690. doi:10.1029/WR020i006p00682 (C84).          !
   ! van Genuchten MT. 1980. A closed-form equation for predicting the hydraulic           !
   !    conductivity of unsaturated soils1. Soil Sci. Soc. Am. J. 44: 892-898.             !
   !    doi:10.2136/sssaj1980.03615995004400050002x (vG80).                                !
   ! Hodnett M , Tomasella J. 2002. Marked differences between van Genuchten soil          !
   !    water-retention parameters for temperate and tropical soils: a new                 !
   !    water-retention pedo-transfer functions developed for tropical soils. Geoderma     !
   !    108: 155-180. doi:10.1016/S0016-7061(02)00105-2 (HT02).                            !
   ! Mualem Y. 1976. A new model for predicting the hydraulic conductivity of unsaturated  !
   !    porous media. Water Resour. Res., 12: 513-522. doi:10.1029/WR012i003p00513 (M76).  !
   ! Romano N , Santini A. 2002. Field. In: Methods of soil analysis: Part 4 physical      !
   !    methods (eds. Dane JH. & Topp GC.). Soil Science Society of America, Madison, WI,  !
   !    SSSA Book Series 5.4, chap. 3.3.3, pp. 721--738 (RS02).                            !
   ! Tomasella J , Hodnett MG. 1998. Estimating soil water retention characteristics from  !
   !    limited data in Brazilian Amazonia. Soil Sci. 163: 190-202.                        !
   !    doi:10.1097/00010694-199803000-00003 (TH98).                                       !
   !---------------------------------------------------------------------------------------!
   type soil_class
      character(len=4) :: key      ! Soil texture acronym                         [     ---]
      character(len=4) :: method   ! Method for pedotransfer function             [     ---]
      real(kind=4)     :: xsand    ! Volumetric fraction of sand                  [     ---]
      real(kind=4)     :: xsilt    ! Volumetric fraction of silt                  [     ---]
      real(kind=4)     :: xclay    ! Volumetric fraction of clay                  [     ---]
      real(kind=4)     :: slsoc    ! Soil organic carbon content.                 [   kg/kg]
      real(kind=4)     :: slph     ! Soil pH                                      [    0-14]
      real(kind=4)     :: slcec    ! Soil cation exchange capacity                [  mol/kg]
      real(kind=4)     :: sldbd    ! "Dry" bulk density                           [   kg/m3]
      real(kind=4)     :: soilre   ! Residual soil moisture (-Infinite potential) [   m3/m3]
      real(kind=4)     :: soilcp   ! Dry soil capacity (at -3.1MPa)               [   m3/m3]
      real(kind=4)     :: soilwp   ! Wilting point capacity (at -1.5MPa)          [   m3/m3]
      real(kind=4)     :: soilfr   ! Fire threshold                               [   m3/m3]
      real(kind=4)     :: soilld   ! Leaf drop threshold                          [   m3/m3]
      real(kind=4)     :: sfldcap  ! Soil field capacity                          [   m3/m3]
      real(kind=4)     :: soilbp   ! Soil moisture at bubble point (vG80 only)    [   m3/m3]
      real(kind=4)     :: slmsts   ! Soil moisture at saturation                  [   m3/m3]
      real(kind=4)     :: soilpo   ! Soil porosity (used only for vG80)           [   m3/m3]
      real(kind=4)     :: slpotcp  ! Water potential for dry soil                 [       m]
      real(kind=4)     :: slpotwp  ! Water potential for wilting point            [       m]
      real(kind=4)     :: slpotfr  ! Fire threshold                               [       m]
      real(kind=4)     :: slpotld  ! Leaf drop threshold                          [       m]
      real(kind=4)     :: slpotfc  ! Water potential for field capacity           [       m]
      real(kind=4)     :: slpotbp  ! Soil matric potential at bubble point        [       m]
      real(kind=4)     :: slpots   ! Soil matric potential at saturation          [       m]
      real(kind=4)     :: slpotpo  ! Soil matric potential at porosity            [       m]
      real(kind=4)     :: sltt     ! Pore tortuosity parameter                    [     n/d]
      real(kind=4)     :: slnm     ! Pore-size distr. idx: n(vG80); lambda(BC64)  [     n/d]
      real(kind=4)     :: slbs     ! 1/lambda or 1/n (aka C74's b factor)         [     n/d]
      real(kind=4)     :: slmm     ! vG80: m=1-1/n.   BC64: m=tort+2+2/lambda     [     n/d]
      real(kind=4)     :: slmu     ! mu = -1/m.                                   [     n/d]
      real(kind=4)     :: malpha   ! van Genuchten's -alpha (1/slpotbp)           [     1/m]
      real(kind=4)     :: slcons   ! hydraulic conductivity at saturation         [     m/s]
      real(kind=4)     :: fhydraul ! Surface value for slcons                     [     m/s]
      real(kind=4)     :: slcpd    ! Specific heat of dry soil                    [  J/m3/K]
      real(kind=4)     :: thcond0  ! First coefficient for thermal conductivity   [   W/m/K]
      real(kind=4)     :: thcond1  ! Second coefficient for thermal conductivity  [   W/m/K]
      real(kind=4)     :: thcond2  ! Third coefficient for thermal conductivity   [     ---]
      real(kind=4)     :: thcond3  ! Fourth coefficient for thermal conductivity  [     ---]
   end type soil_class
   !----- Double precision version --------------------------------------------------------!
   type soil_class8
      character(len=4) :: key      ! Soil texture acronym                         [     ---]
      character(len=4) :: method   ! Method for pedotransfer function             [     ---]
      real(kind=8)     :: xsand    ! Volumetric fraction of sand                  [     ---]
      real(kind=8)     :: xsilt    ! Volumetric fraction of silt                  [     ---]
      real(kind=8)     :: xclay    ! Volumetric fraction of clay                  [     ---]
      real(kind=8)     :: slsoc    ! Soil organic carbon content.                 [   kg/kg]
      real(kind=8)     :: slph     ! Soil pH                                      [    0-14]
      real(kind=8)     :: slcec    ! Soil cation exchange capacity                [  mol/kg]
      real(kind=8)     :: sldbd    ! "Dry" bulk density                           [   kg/m3]
      real(kind=8)     :: soilre   ! Residual soil moisture (-Infinite potential) [   m3/m3]
      real(kind=8)     :: soilcp   ! Dry soil capacity (at -3.1MPa)               [   m3/m3]
      real(kind=8)     :: soilwp   ! Wilting point capacity (at -1.5MPa)          [   m3/m3]
      real(kind=8)     :: soilfr   ! Fire threshold                               [   m3/m3]
      real(kind=8)     :: soilld   ! Leaf drop threshold                          [   m3/m3]
      real(kind=8)     :: sfldcap  ! Soil field capacity                          [   m3/m3]
      real(kind=8)     :: soilbp   ! Soil moisture at bubble point (vG80 only)    [   m3/m3]
      real(kind=8)     :: slmsts   ! Soil moisture at saturation                  [   m3/m3]
      real(kind=8)     :: soilpo   ! Soil porosity (used only for vG80)           [   m3/m3]
      real(kind=8)     :: slpotcp  ! Water potential for dry soil                 [       m]
      real(kind=8)     :: slpotwp  ! Water potential for wilting point            [       m]
      real(kind=8)     :: slpotfr  ! Fire threshold                               [       m]
      real(kind=8)     :: slpotld  ! Leaf drop threshold                          [       m]
      real(kind=8)     :: slpotfc  ! Water potential for field capacity           [       m]
      real(kind=8)     :: slpotbp  ! Soil matric potential at bubble point        [       m]
      real(kind=8)     :: slpots   ! Soil matric potential at saturation          [       m]
      real(kind=8)     :: slpotpo  ! Soil matric potential at porosity            [       m]
      real(kind=8)     :: sltt     ! Pore tortuosity parameter                    [     n/d]
      real(kind=8)     :: slnm     ! Pore-size distr. idx: n(vG80); lambda(BC64)  [     n/d]
      real(kind=8)     :: slbs     ! 1/lambda or 1/n (aka C74's b factor)         [     n/d]
      real(kind=8)     :: slmm     ! vG80: m=1-1/n.   BC64: m=tort+2+2/lambda     [     n/d]
      real(kind=8)     :: slmu     ! mu = -1/m.                                   [     n/d]
      real(kind=8)     :: malpha   ! van Genuchten's -alpha (1/slpotbp)           [     1/m]
      real(kind=8)     :: slcons   ! hydraulic conductivity at saturation         [     m/s]
      real(kind=8)     :: fhydraul ! Surface value for slcons                     [     m/s]
      real(kind=8)     :: slcpd    ! Specific heat of dry soil                    [  J/m3/K]
      real(kind=8)     :: thcond0  ! First coefficient for thermal conductivity   [   W/m/K]
      real(kind=8)     :: thcond1  ! Second coefficient for thermal conductivity  [   W/m/K]
      real(kind=8)     :: thcond2  ! Third coefficient for thermal conductivity   [     ---]
      real(kind=8)     :: thcond3  ! Fourth coefficient for thermal conductivity  [     ---]
   end type soil_class8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      These variables are scalars that will be assigned in ed_params.f90, but may be   !
   ! overwritten when reading multiple sites or running a history initialisation.          !
   !---------------------------------------------------------------------------------------!
   integer         , dimension(ed_nstyp) :: slhydro_ref
   character(len=4), dimension(ed_nstyp) :: slxkey_ref
   real(kind=4)    , dimension(ed_nstyp) :: slxsand_ref
   real(kind=4)    , dimension(ed_nstyp) :: slxsilt_ref
   real(kind=4)    , dimension(ed_nstyp) :: slxclay_ref
   real(kind=4)    , dimension(ed_nstyp) :: slsoc_ref
   real(kind=4)    , dimension(ed_nstyp) :: slph_ref
   real(kind=4)    , dimension(ed_nstyp) :: slcec_ref
   real(kind=4)    , dimension(ed_nstyp) :: sldbd_ref
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
      use ed_max_dims, only : ed_nstyp ! ! intent(in)
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
   !     This subroutine initialises all elements in the soil and soil8 structures, to     !
   ! avoid FPE.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine ed_init_soil()
      use ed_max_dims, only : undef_real      & ! intent(in)
                            , undef_dble      & ! intent(in)
                            , ed_nstyp        ! ! intent(in)

      implicit none
      !------ Local variables. ------------------------------------------------------------!
      integer :: s
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Initialise variables with undefined.  Every parameter should be assigned in    !
      ! ed_params.f90.                                                                     !
      !------------------------------------------------------------------------------------!
      do s=1,ed_nstyp
         !----- Single precision. ---------------------------------------------------------!
         soil(s)%key      = 'Rien'
         soil(s)%method   = 'RIEN'
         soil(s)%xsand    = undef_real
         soil(s)%xsilt    = undef_real
         soil(s)%xclay    = undef_real
         soil(s)%slsoc    = undef_real
         soil(s)%slph     = undef_real
         soil(s)%slcec    = undef_real
         soil(s)%sldbd    = undef_real
         soil(s)%soilre   = undef_real
         soil(s)%soilcp   = undef_real
         soil(s)%soilwp   = undef_real
         soil(s)%soilfr   = undef_real
         soil(s)%soilld   = undef_real
         soil(s)%sfldcap  = undef_real
         soil(s)%soilbp   = undef_real
         soil(s)%slmsts   = undef_real
         soil(s)%soilpo   = undef_real
         soil(s)%slpotcp  = undef_real
         soil(s)%slpotwp  = undef_real
         soil(s)%slpotfr  = undef_real
         soil(s)%slpotld  = undef_real
         soil(s)%slpotfc  = undef_real
         soil(s)%slpotbp  = undef_real
         soil(s)%slpots   = undef_real
         soil(s)%slpotpo  = undef_real
         soil(s)%sltt     = undef_real
         soil(s)%slnm     = undef_real
         soil(s)%slbs     = undef_real
         soil(s)%slmm     = undef_real
         soil(s)%slmu     = undef_real
         soil(s)%malpha   = undef_real
         soil(s)%slcons   = undef_real
         soil(s)%fhydraul = undef_real
         soil(s)%slcpd    = undef_real
         soil(s)%thcond0  = undef_real
         soil(s)%thcond1  = undef_real
         soil(s)%thcond2  = undef_real
         soil(s)%thcond3  = undef_real
         !---------------------------------------------------------------------------------!




         !----- Double precision. ---------------------------------------------------------!
         soil8(s)%key      = 'Rien'
         soil8(s)%method   = 'RIEN'
         soil8(s)%xsand    = undef_dble
         soil8(s)%xsilt    = undef_dble
         soil8(s)%xclay    = undef_dble
         soil8(s)%slsoc    = undef_dble
         soil8(s)%slph     = undef_dble
         soil8(s)%slcec    = undef_dble
         soil8(s)%sldbd    = undef_dble
         soil8(s)%soilre   = undef_dble
         soil8(s)%soilcp   = undef_dble
         soil8(s)%soilwp   = undef_dble
         soil8(s)%soilfr   = undef_dble
         soil8(s)%soilld   = undef_dble
         soil8(s)%sfldcap  = undef_dble
         soil8(s)%soilbp   = undef_dble
         soil8(s)%slmsts   = undef_dble
         soil8(s)%soilpo   = undef_dble
         soil8(s)%slpotcp  = undef_dble
         soil8(s)%slpotwp  = undef_dble
         soil8(s)%slpotfr  = undef_dble
         soil8(s)%slpotld  = undef_dble
         soil8(s)%slpotfc  = undef_dble
         soil8(s)%slpotbp  = undef_dble
         soil8(s)%slpots   = undef_dble
         soil8(s)%slpotpo  = undef_dble
         soil8(s)%sltt     = undef_dble
         soil8(s)%slnm     = undef_dble
         soil8(s)%slbs     = undef_dble
         soil8(s)%slmm     = undef_dble
         soil8(s)%slmu     = undef_dble
         soil8(s)%malpha   = undef_dble
         soil8(s)%slcons   = undef_dble
         soil8(s)%fhydraul = undef_dble
         soil8(s)%slcpd    = undef_dble
         soil8(s)%thcond0  = undef_dble
         soil8(s)%thcond1  = undef_dble
         soil8(s)%thcond2  = undef_dble
         soil8(s)%thcond3  = undef_dble
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_init_soil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function creates the soil table.  This used to be in ed_params.f90, but     !
   ! this is now called from sfcdata_ed, so the contents from the history file will have   !
   ! the last word in case RUNTYPE = 'HISTORY'.                                            !
   !---------------------------------------------------------------------------------------!
   subroutine ed_gen_soil_table()
      use detailed_coms  , only : idetailed   ! ! intent(in)
      use ed_max_dims    , only : str_len     & ! intent(in)
                                , ed_nstyp    ! ! intent(in)
      use phenology_coms , only : thetacrit   ! ! intent(in)
      use disturb_coms   , only : sm_fire     ! ! intent(in)
      use consts_coms    , only : grav        & ! intent(in)
                                , wdns        & ! intent(in)
                                , hr_sec      & ! intent(in)
                                , day_sec     ! ! intent(in)

      implicit none
      !----- Local variables. -------------------------------------------------------------!
      integer      :: s                 ! Soil texture flag
      logical      :: print_soil_table  ! Print parameter table?                  [    T|F]
      real(kind=4) :: soilep            ! Effective porosity (O19)                [  m3/m3]
      real(kind=4) :: slpot33           ! Potential for EP   (O19)                [      m]
      real(kind=4) :: slcons_mmhr       ! Sat. hydraulic conduct.                 [  mm/hr]
      real(kind=4) :: slcpd_mjm3k       ! Soil heat capacity                      [MJ/m3/K]
      real(kind=4) :: ksand             ! k-factor for sand (de Vries model)
      real(kind=4) :: ksilt             ! k-factor for silt (de Vries model)
      real(kind=4) :: kclay             ! k-factor for clay (de Vries model)
      real(kind=4) :: kair              ! k-factor for air  (de Vries model)
      !----- Local constants. -------------------------------------------------------------!
      real(kind=4), parameter :: fieldcp_K   =  0.1     ! hydr. cond.: field cap.  [mm/day]
      real(kind=4), parameter :: slpots_MPa  = -0.0005  ! Saturation for vG80      [   MPa]
      real(kind=4), parameter :: slpot33_MPa = -0.033   ! Potl. for soilep (O19)   [   MPa]
      real(kind=4), parameter :: slpotfc_MPa = -0.010   ! Field capacity (TH98)    [   MPa]
      real(kind=4), parameter :: slpotdg_MPa = -5.0     ! Matric pot.: airdry vG80 [   MPa]
      real(kind=4), parameter :: slpotcp_MPa = -3.1     ! Matric pot.: airdry BC64 [   MPa]
      real(kind=4), parameter :: slpotwp_MPa = -1.5     ! Matric pot.: wilting pt  [   MPa]
      real(kind=4), parameter :: sand_hcapv  =  2.128e6 ! Sand vol. heat cap.      [J/m3/K]
      real(kind=4), parameter :: clay_hcapv  =  2.385e6 ! Clay vol. heat cap.      [J/m3/K]
      real(kind=4), parameter :: silt_hcapv  =  2.256e6 ! Silt vol. heat cap. (*)  [J/m3/K]
      real(kind=4), parameter :: air_hcapv   =  1.212e3 ! Air vol. heat cap.       [J/m3/K]
      real(kind=4), parameter :: sand_thcond = 8.80     ! Sand thermal conduct.    [ W/m/K]
      real(kind=4), parameter :: clay_thcond = 2.92     ! Clay thermal conduct.    [ W/m/K]
      real(kind=4), parameter :: silt_thcond = 5.87     ! Silt therm. conduct. (*) [ W/m/K]
      real(kind=4), parameter :: air_thcond  = 0.025    ! Air thermal conduct.     [ W/m/K]
      real(kind=4), parameter :: h2o_thcond  = 0.57     ! Water thermal conduct.   [ W/m/K]
      !------ Name for the parameter table. -----------------------------------------------!
      character(len=str_len), parameter :: soil_table_fn = 'soil_properties.txt'
      !------------------------------------------------------------------------------------!
      ! (*) If anyone has the heat capacity and thermal conductivity for silt, please feel !
      !     free to add it in here, I didn't find any.  Apparently no one knows, and I've  !
      !     seen in other models that people just assume either the same as sand or the    !
      !     average.  Here I'm just using halfway.  I think the most important thing is to !
      !     take into account the soil and the air, which are the most different.          !
      !                                                                                    !
      ! Sand (quartz), clay, air, and water heat capacities and thermal conductivities     !
      ! values are from:                                                                   !
      !     Monteith and Unsworth, 2008: Environmental Physics.                            !
      !         Academic Press, Third Edition. Table 15.1, p. 292                          !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Initialise the soil and soil8 structures.                                       !
      !------------------------------------------------------------------------------------!
      call ed_init_soil()
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Assign soil properties from the reference vectors.  The reference may have     !
      ! changed from default depending on the user settings from ED2IN or based on a       !
      ! HISTORY restart.                                                                   !
      !------------------------------------------------------------------------------------!
      do s=1,ed_nstyp
         soil(s)%key   = slxkey_ref (s)
         soil(s)%xsand = slxsand_ref(s)
         soil(s)%xclay = slxclay_ref(s)
         soil(s)%xsilt = slxsilt_ref(s)
         soil(s)%slsoc = slsoc_ref  (s)
         soil(s)%slph  = slph_ref   (s)
         soil(s)%slcec = slcec_ref  (s)
         soil(s)%sldbd = sldbd_ref  (s)
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Calculate method- and texture-dependent properties.  For a general overview,   !
      ! check (M14).  Additional references correspond to specific parametrisations.       !
      !                                                                                    !
      ! References:                                                                        !
      !                                                                                    !
      ! Brooks RH , Corey AT. 1964. Hydraulic properties of porous media.                  !
      !    Hydrology Papers 3,  Colorado State University, Fort Collins, U.S.A (BC64).     !
      ! Marthews TR, Quesada CA, Galbraith DR, Malhi Y, Mullins CE, Hodnett MG ,           !
      !    Dharssi I. 2014. High-resolution hydraulic parameter maps for surface soils in  !
      !    tropical South America. Geosci. Model Dev. 7: 711-723.                          !
      !    doi:10.5194/gmd-7-711-2014 (M14).                                               !
      ! Campbell GS. 1974. A simple method for determining unsaturated conductivity from   !
      !    moisture retention data. Soil Science 117: 311-314.                             !
      !    doi:10.1097/00010694-197406000-00001 (C74).                                     !
      ! Cosby BJ, Hornberger GM, Clapp RB , Ginn TR. 1984. A statistical exploration of    !
      !    the relationships of soil moisture characteristics to the physical properties   !
      !    of soils. Water Resour. Res. 20: 682-690. doi:10.1029/WR020i006p00682 (C84).    !
      ! van Genuchten MT. 1980. A closed-form equation for predicting the hydraulic        !
      !    conductivity of unsaturated soils1. Soil Sci. Soc. Am. J. 44: 892-898.          !
      !    doi:10.2136/sssaj1980.03615995004400050002x (vG80).                             !
      ! Hodnett M , Tomasella J. 2002. Marked differences between van Genuchten soil       !
      !    water-retention parameters for temperate and tropical soils: a new              !
      !    water-retention pedo-transfer functions developed for tropical soils. Geoderma  !
      !    108: 155-180. doi:10.1016/S0016-7061(02)00105-2 (HT02).                         !
      ! Montzka C, Herbst M, Weihermuller L, Verhoef A , Vereecken H. 2017. A global data  !
      !    set of soil hydraulic properties and sub-grid variability of soil water         !
      !    retention and hydraulic conductivity curves. Earth Syst. Sci. Data, 9: 529-543. !
      !    doi:10.5194/essd-9-529-2017 (M17).                                              !
      ! Mualem Y. 1976. A new model for predicting the hydraulic conductivity of           !
      !    unsaturated porous media. Water Resour. Res., 12: 513-522.                      !
      !    doi:10.1029/WR012i003p00513 (M76).                                              !
      ! Ottoni MV, Ottoni Filho TB, Lopes-Assad MLR , Rotunno Filho OC. 2019. Pedotransfer !
      !    functions for saturated hydraulic conductivity using a database with temperate  !
      !    and tropical climate soils. J. Hydrol., 575: 1345-1358.                         !
      !    doi:10.1016/j.jhydrol.2019.05.050 (O19).                                        !
      ! Romano N , Santini A. 2002. Field. In: Methods of soil analysis: Part 4 physical   !
      !    methods (eds. Dane JH. & Topp GC.). Soil Science Society of America,            !
      !    Madison, WI, SSSA Book Series 5.4, chap. 3.3.3, pp. 721--738 (RS02).            !
      ! Schaap MG , Leij FJ. 2000. Improved prediction of unsaturated hydraulic            !
      !    conductivity with the Mualem- van Genuchten model. Soil Sci. Soc. Am. J.,       !
      !    64: 843-851. doi:10.2136/sssaj2000.643843x (SL00).                              !
      ! Tomasella J , Hodnett MG. 1998. Estimating soil water retention characteristics    !
      !    from limited data in Brazilian Amazonia. Soil Sci. 163: 190-202.                !
      !    doi:10.1097/00010694-199803000-00003 (TH98).                                    !
      !------------------------------------------------------------------------------------!
      do s=1,ed_nstyp

         !----- Check soil texture.  Peat and bedrock must be handled separately. ---------!
         select case (slhydro_ref(s))
         case (12)
            !------------------------------------------------------------------------------!
            !     Peat.  We always use BC64-M76 approach.  This modify hydraulic           !
            ! properties to account for high soil organic content.  This class becomes     !
            ! obsolete for SOIL_HYDRO_SCHEME=2 because we can account for SOC directly.    !
            !                                                                              !
            ! MLO - I noticed that most parameters do not correspond to what is            !
            !       implemented in LEAF3 (as of RAMS-6.0).  I left the LEAF3 values as     !
            !       comments next to the default values but someone running ED2 for peats  !
            !       should check. I think the LEAF3 values intuitively make more sense for !
            !       peat.                                                                  !
            !------------------------------------------------------------------------------!
            soil(s)%method = 'BC64'

            !----- Peat, use the default value from LEAF3. --------------------------------!
            soil(s)%slcons  =  8.0e-6    ! ED-2.2 2.357930e-6
            !------------------------------------------------------------------------------!

            !---- Pore tortuosity factor. Assumed 1 to be consistent with BC64. -----------!
            soil(s)%sltt = 1.0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Pore-size distribution factor (slnm, aka lambda) and its inverse (slbs,  !
            ! aka BC64's "b" factor).                                                      !
            !------------------------------------------------------------------------------!
            soil(s)%slbs    =  7.75  ! ED-2.2 6.180000
            soil(s)%slnm    = 1. / soil(s)%slbs
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Ancillary parameters used for hydraulic conductivity.                    !
            !------------------------------------------------------------------------------!
            soil(s)%slmm = 2. + soil(s)%sltt + 2. * soil(s)%slbs
            soil(s)%slmu = -1. / soil(s)%slmm
            !------------------------------------------------------------------------------!


            !----- Saturation potential [m]. ----------------------------------------------!
            soil(s)%slpots  = -0.356                ! ED-2.2 -0.534564359
            soil(s)%slpotbp = soil(s)%slpots        ! Bubbling point, assume saturation
            soil(s)%slpotpo = soil(s)%slpots        ! Porosity, assume saturation
            soil(s)%malpha  = 1. / soil(s)%slpotbp  ! Inverse of bubbling point (not used)
            !------------------------------------------------------------------------------!



            !----- Soil moisture at saturation [m3/m3]. -----------------------------------!
            soil(s)%slmsts  =  0.863          ! ED-2.2 0.469200
            soil(s)%soilbp  =  soil(s)%slmsts ! Assume the same as saturation
            soil(s)%soilpo  =  soil(s)%slmsts ! Assume the same as saturation
            !------------------------------------------------------------------------------!


            !---- Field capacity [m3/m3] and potential at field capacity [m]. -------------!
            soil(s)%sfldcap = 0.535  ! ED-2.2 0.285709966
            soil(s)%slpotfc = matric_potential(s,soil(s)%sfldcap)
            !------------------------------------------------------------------------------!


            !----- Residual moisture [m3/m3].  Ignored in C74 and the default ED2 method. -!
            soil(s)%soilre  = 0.
            !------------------------------------------------------------------------------!


            !----- Heat capacity. ---------------------------------------------------------!
            soil(s)%slcpd   =  874000.
            !------------------------------------------------------------------------------!
         case (13)
            !----- Bedrock.  Hydraulics is disabled, only heat capacity is needed. --------!
            soil(s)%method  = 'BDRK'
            soil(s)%slcons  = 0.0
            soil(s)%sltt    = 0.0
            soil(s)%slnm    = 1.0
            soil(s)%slbs    = 1.0
            soil(s)%slmm    = 1.0
            soil(s)%slmu    = 1.0
            soil(s)%malpha  = 0.0
            soil(s)%slpots  = 0.0
            soil(s)%slmsts  = 0.0
            soil(s)%slpotbp = 0.0
            soil(s)%soilbp  = 0.0
            soil(s)%slpotpo = 0.0
            soil(s)%soilpo  = 0.0
            soil(s)%soilre  = 0.0
            soil(s)%sfldcap = 0.0
            soil(s)%slpotfc = 0.0
            soil(s)%slcpd   = 2130000.
            !------------------------------------------------------------------------------!
         case (0)
            !------------------------------------------------------------------------------!
            !     Pedotransfer functions from BC64/M76.  Unless noted otherwise, the       !
            ! parameters for the functions are from C84, based on measurements in the      !
            ! United States.                                                               !
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Hydraulic conductivity at saturation [m/s].                             !
            !------------------------------------------------------------------------------!
            soil(s)%slcons  = (10.**(-0.60 + 1.26*soil(s)%xsand - 0.64*soil(s)%xclay))     &
                            * 0.0254/hr_sec
            !------------------------------------------------------------------------------!


            !---- Flag for method. --------------------------------------------------------!
            soil(s)%method = 'BC64'
            !------------------------------------------------------------------------------!


            !---- Pore tortuosity factor. Assumed 1 to be consistent with BC64. -----------!
            soil(s)%sltt = 1.0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Pore-size distribution factor (slnm, aka lambda) and its inverse         !
            ! (slbs, aka BC64's "b" factor).                                               !
            !------------------------------------------------------------------------------!
            soil(s)%slnm = 1. / (3.10 + 15.7*soil(s)%xclay - 0.3*soil(s)%xsand)
            soil(s)%slbs = 1. / soil(s)%slnm
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Ancillary parameters used for hydraulic conductivity.                    !
            !------------------------------------------------------------------------------!
            soil(s)%slmm = 2. + soil(s)%sltt + 2. * soil(s)%slbs
            soil(s)%slmu = -1. / soil(s)%slmm
            !------------------------------------------------------------------------------!


            !----- Saturation potential [m]. ----------------------------------------------!
            soil(s)%slpots  = -1.                                                          &
                            * (10.**(2.17 - 0.63*soil(s)%xclay - 1.58*soil(s)%xsand))      &
                            * 0.01
            soil(s)%slpotbp = soil(s)%slpots       ! Bubbling point, assume saturation
            soil(s)%slpotpo = soil(s)%slpots       ! Porosity, assume saturation
            soil(s)%malpha  = 1. / soil(s)%slpotbp
            !------------------------------------------------------------------------------!


            !----- Soil moisture at saturation (porosity) [m3/m3]. ------------------------!
            soil(s)%slmsts  = 0.01 * (50.5 - 14.2*soil(s)%xsand - 3.7*soil(s)%xclay)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Bubbling point soil moisture and porosity, assume saturation [m3/m3].    !
            !------------------------------------------------------------------------------!
            soil(s)%soilbp  = soil(s)%slmsts ! Bubbling point
            soil(s)%soilpo  = soil(s)%slmsts ! Porosity
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Residual moisture [m3/m3]. Ignored in C74 and the default ED2 method.    !
            !------------------------------------------------------------------------------!
            soil(s)%soilre = 0.0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Field capacity is defined based on hydraulic conductivity of 0.1        !
            ! mm/day, following RS02.                                                      !
            !------------------------------------------------------------------------------!
            soil(s)%sfldcap =  soil(s)%slmsts                                              &
                            *  ( soil(s)%slcons / ( fieldcp_K / ( wdns * day_sec ) ) )     &
                            ** soil(s)%slmu
            soil(s)%slpotfc = matric_potential(s,soil(s)%sfldcap)
            !------------------------------------------------------------------------------!

         case (1)
            !------------------------------------------------------------------------------!
            !     Pedotransfer functions from BC64/M76.  Unless noted otherwise, the       !
            ! parameters for the functions are from TH98, based on measurements in the     !
            ! Brazilian Amazon.                                                            !
            !------------------------------------------------------------------------------!



            !---- Flag for method. --------------------------------------------------------!
            soil(s)%method = 'BC64'
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Hydraulic conductivity at saturation [m/s].  Use C84 settings,          !
            ! following M14.                                                               !
            !------------------------------------------------------------------------------!
            soil(s)%slcons  = (10.**(-0.60 + 1.26*soil(s)%xsand - 0.64*soil(s)%xclay))     &
                            * 0.0254/hr_sec
            !------------------------------------------------------------------------------!


            !---- Pore tortuosity factor. Assumed 0.5, following M14. ---------------------!
            soil(s)%sltt = 0.5
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Pore-size distribution factor (slnm, aka lambda) and its inverse         !
            ! (slbs, aka BC64's "b" factor).                                               !
            !------------------------------------------------------------------------------!
            soil(s)%slnm = exp( - 1.197 - 0.417 * soil(s)%xsilt + 0.450 * soil(s)%xclay    &
                                - 8.940 * soil(s)%xsilt * soil(s)%xclay                    &
                                + 10.00 * soil(s)%xsilt * soil(s)%xsilt * soil(s)%xclay    )
            soil(s)%slbs = 1./ soil(s)%slnm
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Ancillary parameters used for hydraulic conductivity.                    !
            !------------------------------------------------------------------------------!
            soil(s)%slmm = 2. + soil(s)%sltt + 2. * soil(s)%slbs
            soil(s)%slmu = -1. / soil(s)%slmm
            !------------------------------------------------------------------------------!


            !----- Saturation potential [m]. ----------------------------------------------!
            soil(s)%slpots  = -1. / grav                                                   &
                            * ( 0.285 + 7.33 * soil(s)%xsilt * soil(s)%xsilt               &
                              - 1.30 * soil(s)%xsilt * soil(s)%xclay                       &
                              + 3.60 * soil(s)%xsilt * soil(s)%xsilt * soil(s)%xclay )
            soil(s)%slpotbp = soil(s)%slpots       ! Bubbling point, assume saturation
            soil(s)%slpotpo = soil(s)%slpots       ! Porosity, assume saturation
            soil(s)%malpha  = 1. / soil(s)%slpotbp
            !------------------------------------------------------------------------------!


            !----- Soil moisture at saturation (porosity) [m3/m3]. ------------------------!
            soil(s)%slmsts  = 0.4061  + 0.165 * soil(s)%xsilt + 0.162 * soil(s)%xclay      &
                            + 1.37e-3 * soil(s)%xsilt * soil(s)%xsilt                      &
                            + 1.80e-5 * soil(s)%xsilt * soil(s)%xsilt * soil(s)%xclay
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Bubbling point soil moisture and porosity, assume saturation [m3/m3].    !
            !------------------------------------------------------------------------------!
            soil(s)%soilbp  = soil(s)%slmsts ! Bubbling point
            soil(s)%soilpo  = soil(s)%slmsts ! Porosity
            !------------------------------------------------------------------------------!



            !----- Residual moisture [m3/m3].  --------------------------------------------!
            soil(s)%soilre = max( 0.0                                                      &
                                , - 0.02095 + 0.047 * soil(s)%xsilt                        &
                                            + 0.431 * soil(s)%xclay                        &
                                            - 0.00827 * soil(s)%xsilt * soil(s)%xclay )
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Field capacity is defined based on hydraulic conductivity of 0.1        !
            ! mm/day, following RS02.                                                      !
            !------------------------------------------------------------------------------!
            soil(s)%slpotfc = slpotfc_MPa * 1.e6 / (grav * wdns)
            soil(s)%sfldcap = soil_moisture(s,soil(s)%slpotfc)
            !------------------------------------------------------------------------------!

         case (2)
            !------------------------------------------------------------------------------!
            !     Pedotransfer functions from vG80/M76.  Unless noted otherwise,           !
            ! the parameters for the function are from HT02, based on measurements in      !
            ! the Brazilian Amazon.                                                        !
            !------------------------------------------------------------------------------!


            !---- Flag for method. --------------------------------------------------------!
            soil(s)%method = 'vG80'
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Pore tortuosity factor. M14 assumed 0.5, but there is evidence that     !
            ! this parameter should be regarded as empirical and some studies suggested    !
            ! that it should be even negative (e.g., SL00 and M17).  We follow SL00 and    !
            ! assume the parameter to be -1.0.                                             !
            !------------------------------------------------------------------------------!
            soil(s)%sltt = -1.0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Pore-size distribution factor (slnm, aka lambda) and its inverse         !
            ! (slbs, aka BC64's "b" factor).                                               !
            !------------------------------------------------------------------------------!
            soil(s)%slnm = exp( 0.62986 - 0.833 * soil(s)%xclay - 0.529 * soil(s)%slsoc    &
                              + 0.00593 * soil(s)%slph                                     &
                              + 0.700 * soil(s)%xclay * soil(s)%xclay                      &
                              - 1.400 * soil(s)%xsand * soil(s)%xsilt                      )
            soil(s)%slbs = 1./ soil(s)%slnm
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Ancillary parameters used for hydraulic conductivity.                    !
            !------------------------------------------------------------------------------!
            soil(s)%slmm = 1. - 1. / soil(s)%slnm
            soil(s)%slmu = soil(s)%slnm / (1. - soil(s)%slnm)
            !------------------------------------------------------------------------------!


            !----- Bubbling point potential [m].  Assume equivalent to saturation. --------!
            soil(s)%slpotbp = -1. / grav                                                   &
                            * exp(  0.02294 + 3.526 * soil(s)%xsilt                        &
                                 - 2.440 * soil(s)%slsoc + 0.076 * soil(s)%slcec           &
                                  + 0.11331 * soil(s)%slph                                 &
                                  - 1.90000 * soil(s)%xsilt * soil(s)%xsilt      )
            soil(s)%malpha  = 1. / soil(s)%slpotbp
            !------------------------------------------------------------------------------!



            !----- Soil moisture and potential at porosity [m3/m3]. -----------------------!
            soil(s)%soilpo  = 0.81799 + 0.099 * soil(s)%xclay                              &
                            - 3.142e-4 * soil(s)%sldbd + 0.01800 * soil(s)%slcec           &
                            + 0.00451 * soil(s)%slph                                       &
                            - 0.050 * soil(s)%xsand * soil(s)%xclay
            soil(s)%slpotpo = matric_potential(s,soil(s)%soilpo)
            !------------------------------------------------------------------------------!



            !----- Residual moisture [m3/m3].  --------------------------------------------!
            soil(s)%soilre = max( 0.0                                                      &
                                , 0.22733 - 0.164 * soil(s)%xsand                          &
                                + 0.235 * soil(s)%slcec                                    &
                                - 0.00831 * soil(s)%slph                                   &
                                + 0.18 * soil(s)%xclay * soil(s)%xclay                     &
                                + 0.26 * soil(s)%xsand * soil(s)%xclay                    )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Soil moisture at "saturation".  The vG80 approach assumes that actual    !
            ! saturation occurs when soil matric potential is zero, which causes           !
            ! singularities in many applications.  To prevent FPE errors, we assume        !
            ! that water potential zero corresponds to porosity (admittedly this is not    !
            ! entirely accurate), and impose "saturation" for ED-2.2 purposes to be        !
            ! when matric potential is -0.5 kPa, similar to the highest saturation         !
            ! potential values obtained through Cosby et al. (1984) parametrisation.       !
            !------------------------------------------------------------------------------!
            soil(s)%slpots = slpots_MPa * 1.e6 / (grav * wdns)
            soil(s)%slmsts = soil_moisture(s,soil(s)%slpots)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Soil moisture at bubbling point.  Unlike other schemes, we account       !
            ! for the differences between bubbling point and porosity.                     !
            !------------------------------------------------------------------------------!
            soil(s)%soilbp = soil_moisture(s,soil(s)%slpotbp)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Field capacity is defined based on hydraulic conductivity of 0.1        !
            ! mm/day, following RS02.                                                      !
            !------------------------------------------------------------------------------!
            soil(s)%slpotfc = slpotfc_MPa * 1.e6 / (grav * wdns)
            soil(s)%sfldcap = soil_moisture(s,soil(s)%slpotfc)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Hydraulic conductivity at saturation [m/s].  Here we follow O19,        !
            ! which depends upon the "effective porosity" (or difference between actual    !
            ! porosity and soil moisture at -0.033 MPa).                                   !
            !------------------------------------------------------------------------------!
            slpot33         = slpot33_MPa * 1.e6 / (grav * wdns)
            soilep          = max(0.,soil(s)%slmsts - soil_moisture(s,slpot33))
            soil(s)%slcons  = 19.31 / day_sec * soilep ** 1.948
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Additional derived parameters.                                              !
         !---------------------------------------------------------------------------------!
         select case (slhydro_ref(s))
         case (13)
            !----- Bedrock, do nothing. ---------------------------------------------------!
            soil(s)%slpotcp  = 0.0
            soil(s)%slpotwp  = 0.0
            soil(s)%slpotfr  = 0.0
            soil(s)%slpotld  = 0.0
            soil(s)%slpotfc  = 0.0
            soil(s)%slpotbp  = 0.0
            soil(s)%slpots   = 0.0
            soil(s)%slpotpo  = 0.0
            soil(s)%soilcp   = 0.0
            soil(s)%soilwp   = 0.0
            soil(s)%soilfr   = 0.0
            soil(s)%soilld   = 0.0
            soil(s)%sfldcap  = 0.0
            soil(s)%soilbp   = 0.0
            soil(s)%slmsts   = 0.0
            soil(s)%soilpo   = 0.0
            soil(s)%fhydraul = 0.0
            !------------------------------------------------------------------------------!
         case default
            !----- First guess, use water potential. --------------------------------------!
            soil(s)%slpotwp = slpotwp_MPa * 1.e6 / ( grav * wdns )
            select case (trim(soil(s)%method))
            case ('vG80')
               soil(s)%slpotcp = slpotdg_MPa * wdns / grav
            case default
               soil(s)%slpotcp = slpotcp_MPa * wdns / grav
            end select
            soil(s)%soilwp  = soil_moisture(s,soil(s)%slpotwp)
            soil(s)%soilcp  = soil_moisture(s,soil(s)%slpotcp)
            !----- In case soilcp is less than the residual (very unlikely), recalculate it. -!
            if (soil(s)%soilcp < soil(s)%soilre) then
               soil(s)%soilcp  = soil(s)%soilre
               soil(s)%slpotcp = matric_potential(s,soil(s)%soilcp)
            end if
            !----- Because we may have artificially increased soilcp, check soilwp. -------!
            if (soil(s)%soilwp < soil(s)%soilcp) then
               soil(s)%soilwp  = soil(s)%soilcp
               soil(s)%slpotwp = soil(s)%slpotcp
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Find two remaining properties, that depend on the user choices.          !
            !                                                                              !
            ! SOILLD/SLPOTLD.                                                              !
            !    The critical soil moisture below which drought deciduous plants start     !
            !    dropping their leaves.  The sign of input variable THETACRIT matters      !
            !    here.  If the user gave a positive number (or 0),  then the soil moisture !
            !    is a fraction above wilting point.  If it is negative, the value is the   !
            !    potential in MPa.                                                         !
            ! SOILFR/SLPOTFR.                                                              !
            !    The critical soil moisture below which fires may happen, provided that    !
            !    the user wants fires, and that there is enough biomass to burn.  The sign !
            !    of the input variable SM_FIRE matters here.  If the user gave a positive  !
            !    number (or 0), then the soil moisture is a fraction above dry air soil.   !
            !    If it is negative, the value is the potential in MPa.                     !
            !------------------------------------------------------------------------------!
            !----- Leaf drop. -------------------------------------------------------------!
            if (thetacrit >= 0.0) then
               soil(s)%soilld  = soil(s)%soilwp                                            &
                               + thetacrit * (soil(s)%slmsts-soil(s)%soilwp)
               soil(s)%slpotld = matric_potential(s,soil(s)%soilld)
            else
               soil(s)%slpotld = thetacrit * 1.e6 / (grav * wdns)
               soil(s)%soilld  = soil_moisture(s,soil(s)%slpotld)
            end if
            !----- Fire. ------------------------------------------------------------------!
            if (sm_fire >= 0.0) then
               soil(s)%soilfr  = soil(s)%soilcp + sm_fire * (soil(s)%slmsts-soil(s)%soilcp)
               soil(s)%slpotfr = matric_potential(s,soil(s)%soilfr)
            else
               soil(s)%slpotfr = sm_fire * 1.e6 / (grav * wdns)
               soil(s)%soilfr  = soil_moisture(s,soil(s)%slpotfr)
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Define hydraulic parameter decay, similar to TOPMODEL.  We currently use !
            ! the default value of 2.0, following N05's SIMTOP model.                      !
            !                                                                              !
            ! Niu GY, Yang ZL, Dickinson RE , Gulden LE. 2005. A simple TOPMODEL-based     !
            !    runoff parameterization (SIMTOP) for use in global climate models.        !
            !    J. Geophys. Res.-Atmos., 110: D21106. doi:10.1029/2005JD006111 (N05).     !
            !------------------------------------------------------------------------------!
            soil(s)%fhydraul = 2.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Heat capacity (J/m3/K).  Here we take the volume average amongst silt,      !
         ! clay, and sand, and consider the contribution of air sitting in.  In order to   !
         ! keep it simple, we assume that the air fraction won't change, although in       !
         ! reality its contribution should be a function of soil moisture.  Here we use    !
         ! the amount of air in case the soil moisture was halfway between dry air and     !
         ! saturated, so the error is not too biased.                                      !
         !---------------------------------------------------------------------------------!
         soil(s)%slcpd = (1. - soil(s)%slmsts)                                             &
                       * ( soil(s)%xsand * sand_hcapv + soil(s)%xsilt * silt_hcapv         &
                         + soil(s)%xclay * clay_hcapv )                                    &
                       + 0.5 * ( soil(s)%slmsts - soil(s)%soilcp ) * air_hcapv
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Thermal conductivity is the weighted average of thermal conductivities of  !
         ! all materials, although a further weighting factor due to thermal gradient of   !
         ! different materials.  We use the de Vries model described at:                   !
         !                                                                                 !
         ! Camillo, P., T.J. Schmugge, 1981: A computer program for the simulation of heat !
         !     and moisture flow in soils, NASA-TM-82121, Greenbelt, MD, United States.    !
         !                                                                                 !
         ! Parlange, M.B., et al., 1998: Review of heat and water movement in field soils, !
         !    Soil Till. Res., 47(1-2), 5-10.                                              !
         !                                                                                 !
         !---------------------------------------------------------------------------------!
         !---- The k-factors, assuming spherical particles. -------------------------------!
         ksand = 3. * h2o_thcond / ( 2. * h2o_thcond + sand_thcond )
         ksilt = 3. * h2o_thcond / ( 2. * h2o_thcond + silt_thcond )
         kclay = 3. * h2o_thcond / ( 2. * h2o_thcond + clay_thcond )
         kair  = 3. * h2o_thcond / ( 2. * h2o_thcond +  air_thcond )
         !---- The conductivity coefficients. ---------------------------------------------!
         soil(s)%thcond0 = (1. - soil(s)%slmsts )                                          &
                         * ( ksand * soil(s)%xsand * sand_thcond                           &
                           + ksilt * soil(s)%xsilt * silt_thcond                           &
                           + kclay * soil(s)%xclay * clay_thcond )                         &
                         + soil(s)%slmsts * kair * air_thcond
         soil(s)%thcond1 = h2o_thcond - kair * air_thcond
         soil(s)%thcond2 = (1. - soil(s)%slmsts )                                          &
                         * ( ksand * soil(s)%xsand + ksilt * soil(s)%xsilt                 &
                           + kclay * soil(s)%xclay                         )               &
                         + soil(s)%slmsts * kair
         soil(s)%thcond3 = 1. - kair
         !---------------------------------------------------------------------------------!

      end do
      !------------------------------------------------------------------------------------!



      !----- Here we fill soil8, which will be used in Runge-Kutta (double precision). ----!
      do s=1,ed_nstyp
         soil8(s)%key      = soil(s)%key
         soil8(s)%method   = soil(s)%method
         soil8(s)%xsand    = dble(soil(s)%xsand   )
         soil8(s)%xsilt    = dble(soil(s)%xsilt   )
         soil8(s)%xclay    = dble(soil(s)%xclay   )
         soil8(s)%slsoc    = dble(soil(s)%slsoc   )
         soil8(s)%slph     = dble(soil(s)%slph    )
         soil8(s)%slcec    = dble(soil(s)%slcec   )
         soil8(s)%sldbd    = dble(soil(s)%sldbd   )
         soil8(s)%soilre   = dble(soil(s)%soilre  )
         soil8(s)%soilcp   = dble(soil(s)%soilcp  )
         soil8(s)%soilwp   = dble(soil(s)%soilwp  )
         soil8(s)%soilfr   = dble(soil(s)%soilfr  )
         soil8(s)%soilld   = dble(soil(s)%soilld  )
         soil8(s)%sfldcap  = dble(soil(s)%sfldcap )
         soil8(s)%soilbp   = dble(soil(s)%soilbp  )
         soil8(s)%slmsts   = dble(soil(s)%slmsts  )
         soil8(s)%soilpo   = dble(soil(s)%soilpo  )
         soil8(s)%slpotcp  = dble(soil(s)%slpotcp )
         soil8(s)%slpotwp  = dble(soil(s)%slpotwp )
         soil8(s)%slpotfr  = dble(soil(s)%slpotfr )
         soil8(s)%slpotld  = dble(soil(s)%slpotld )
         soil8(s)%slpotfc  = dble(soil(s)%slpotfc )
         soil8(s)%slpotbp  = dble(soil(s)%slpotbp )
         soil8(s)%slpots   = dble(soil(s)%slpots  )
         soil8(s)%slpotpo  = dble(soil(s)%slpotpo )
         soil8(s)%sltt     = dble(soil(s)%sltt    )
         soil8(s)%slnm     = dble(soil(s)%slnm    )
         soil8(s)%slbs     = dble(soil(s)%slbs    )
         soil8(s)%slmm     = dble(soil(s)%slmm    )
         soil8(s)%slmu     = dble(soil(s)%slmu    )
         soil8(s)%malpha   = dble(soil(s)%malpha  )
         soil8(s)%slcons   = dble(soil(s)%slcons  )
         soil8(s)%fhydraul = dble(soil(s)%fhydraul)
         soil8(s)%slcpd    = dble(soil(s)%slcpd   )
         soil8(s)%thcond0  = dble(soil(s)%thcond0 )
         soil8(s)%thcond1  = dble(soil(s)%thcond1 )
         soil8(s)%thcond2  = dble(soil(s)%thcond2 )
         soil8(s)%thcond3  = dble(soil(s)%thcond3 )
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Decide whether to write the table with the soil properties.                    !
      !------------------------------------------------------------------------------------!
      print_soil_table = btest(idetailed,5)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Print the parameters in case the user wants it.                                !
      !------------------------------------------------------------------------------------!
      if (print_soil_table) then
         !----- Open and write header. ----------------------------------------------------!
         open (unit=26,file=trim(soil_table_fn),status='replace',action='write')
         write(unit=26,fmt='(38(a,1x))')        'ISOIL',        ' KEY',        'TYPE'      &
                                        ,'       XSAND','       XSILT','       XCLAY'      &
                                        ,'       SLSOC','        SLPH','       SLCEC'      &
                                        ,'       SLDBD','      SOILRE','      SOILCP'      &
                                        ,'      SOILWP','      SOILFR','      SOILLD'      &
                                        ,'      SOILFC','      SOILBP','      SOILPO'      &
                                        ,'     SLPOTCP','     SLPOTWP','     SLPOTFR'      &
                                        ,'     SLPOTLD','     SLPOTFC','     SLPOTBP'      &
                                        ,'     SLPOTPO','        SLTT','        SLNM'      &
                                        ,'        SLBS','        SLMM','        SLMU'      &
                                        ,'      MALPHA',' SLCONS_MMHR','    FHYDRAUL'      &
                                        ,' SLCPD_MJm3K','     THCOND0','     THCOND1'      &
                                        ,'     THCOND2','     THCOND3'
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop over soil texture types.                                               !
         !---------------------------------------------------------------------------------!
         do s=1,ed_nstyp
            !----- For some variables, we use different units to make them more legible. --!
            slcons_mmhr = soil(s)%slcons*1000.*hr_sec
            slcpd_mjm3k = soil(s)%slcpd*0.001
            !------------------------------------------------------------------------------!

            !----- Add soil characteristics. ----------------------------------------------!
            write(unit=26,fmt='(i5,1x,2(a4,1x),35(f12.5,1x))')                             &
                                      s,adjustr(soil(s)%key),adjustr(soil(s)%method)       &
                                     ,soil(s)%xsand   ,soil(s)%xsilt   ,soil(s)%xclay      &
                                     ,soil(s)%slsoc   ,soil(s)%slph    ,soil(s)%slcec      &
                                     ,soil(s)%sldbd   ,soil(s)%soilre  ,soil(s)%soilcp     &
                                     ,soil(s)%soilwp  ,soil(s)%soilfr  ,soil(s)%soilld     &
                                     ,soil(s)%sfldcap ,soil(s)%soilbp  ,soil(s)%slmsts     &
                                     ,soil(s)%slpotcp ,soil(s)%slpotwp ,soil(s)%slpotfr    &
                                     ,soil(s)%slpotld ,soil(s)%slpotfc ,soil(s)%slpotbp    &
                                     ,soil(s)%slpots  ,soil(s)%sltt    ,soil(s)%slnm       &
                                     ,soil(s)%slbs    ,soil(s)%slmm    ,soil(s)%slmu       &
                                     ,soil(s)%malpha  ,slcons_mmhr     ,soil(s)%fhydraul   &
                                     ,slcpd_mjm3k     ,soil(s)%thcond0 ,soil(s)%thcond1    &
                                     ,soil(s)%thcond2 ,soil(s)%thcond3 
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!


         !----- Close table. --------------------------------------------------------------!
         close(unit=26,status='keep')
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine ed_gen_soil_table
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
   !   2. = saturation            (slmsts)                                                 !
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
   real(kind=4) function matric_potential(s,soil_water)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: s          ! Soil texture                         [  idx]
      real(kind=4), intent(in) :: soil_water ! Soil moisture                        [m3/m3]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist   ! Relative soil moisture               [  ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Select approach based on the soil method and soil texture.                     !
      !------------------------------------------------------------------------------------!
      select case(soil(s)%method)
      case ('BDRK')
         !---- Bedrock. Return zero. ------------------------------------------------------!
         matric_potential = 0.0
         !---------------------------------------------------------------------------------!
      case ('BC64')
         !---------------------------------------------------------------------------------!
         !     Brooks and Corey (1964).                                                    !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water     - soil(s)%soilre )                                    &
                  / ( soil(s)%slmsts - soil(s)%soilre )
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!



         !----- Find the matric potential. ------------------------------------------------!
         matric_potential = soil(s)%slpots / relmoist ** soil(s)%slbs
         !---------------------------------------------------------------------------------!

      case ('vG80')
         !---------------------------------------------------------------------------------!
         !    van Genuchten (1980).                                                        !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water     - soil(s)%soilre )                                    &
                  / ( soil(s)%soilpo - soil(s)%soilre )
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!



         !----- Find the matric potential. ------------------------------------------------!
         matric_potential = soil(s)%slpotbp                                                &
                          * ( relmoist ** soil(s)%slmu - 1.0 ) ** ( 1. - soil(s)%slmm )
         !---------------------------------------------------------------------------------!
      end select
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
   real(kind=8) function matric_potential8(s,soil_water)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: s          ! Soil texture                         [  idx]
      real(kind=8), intent(in) :: soil_water ! Soil moisture                        [m3/m3]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=8)             :: relmoist   ! Relative soil moisture               [  ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Select approach based on the soil method and soil texture.                     !
      !------------------------------------------------------------------------------------!
      select case(soil8(s)%method)
      case ('BDRK')
         !---- Bedrock. Return zero. ------------------------------------------------------!
         matric_potential8 = 0.d0
         !---------------------------------------------------------------------------------!
      case ('BC64')
         !---------------------------------------------------------------------------------!
         !     Brooks and Corey (1964).                                                    !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water      - soil8(s)%soilre )                                    &
                  / ( soil8(s)%slmsts - soil8(s)%soilre )
         relmoist = max(0.d0,min(relmoist,1.d0))
         !---------------------------------------------------------------------------------!



         !----- Find the matric potential. ------------------------------------------------!
         matric_potential8 = soil8(s)%slpots / relmoist ** soil8(s)%slbs
         !---------------------------------------------------------------------------------!

      case ('vG80')
         !---------------------------------------------------------------------------------!
         !    van Genuchten (1980).                                                        !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water      - soil8(s)%soilre )                                  &
                  / ( soil8(s)%soilpo - soil8(s)%soilre )
         relmoist = max(0.d0,min(relmoist,1.d0))
         !---------------------------------------------------------------------------------!



         !----- Find the matric potential. ------------------------------------------------!
         matric_potential8 = soil8(s)%slpotbp                                              &
                           *  ( relmoist ** soil8(s)%slmu - 1.d0 )                         &
                           ** ( 1.d0 - soil8(s)%slmm )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function matric_potential8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil matric potential to soil moisture.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function soil_moisture(s,soil_mstpot)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: s           ! Soil texture                        [  idx]
      real(kind=4), intent(in) :: soil_mstpot ! Soil matric potential               [    m]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmstpot   ! Relative soil matric potential      [  ---]
      real(kind=4)             :: relmoist    ! Relative soil moisture              [  ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Select approach based on the soil method and soil texture.                     !
      !------------------------------------------------------------------------------------!
      select case(soil(s)%method)
      case ('BDRK')
         !---- Bedrock. Return zero. ------------------------------------------------------!
         soil_moisture = 0.0
         !---------------------------------------------------------------------------------!
      case ('BC64')
         !----- Brooks and Corey (1964). --------------------------------------------------!
         relmoist = (soil(s)%slpots / soil_mstpot) ** soil(s)%slnm
         !---------------------------------------------------------------------------------!


         !----- Make sure moisture is bounded. --------------------------------------------!
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!


         !------ Convert relative soil moisture into absolute moisture. -------------------!
         soil_moisture = (1.0 - relmoist) * soil(s)%soilre + relmoist * soil(s)%slmsts
         !---------------------------------------------------------------------------------!
      case ('vG80')
         !----- van Genuchten (1980). -----------------------------------------------------!
         relmstpot = soil(s)%malpha * soil_mstpot
         relmoist  = ( 1. / (1. + relmstpot ** soil(s)%slnm) ) ** soil(s)%slmm
         !---------------------------------------------------------------------------------!


         !----- Make sure moisture is bounded. --------------------------------------------!
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!


         !------ Convert relative soil moisture into absolute moisture. -------------------!
         soil_moisture = (1.0 - relmoist) * soil(s)%soilre + relmoist * soil(s)%soilpo
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function soil_moisture
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function converts soil moisture to hydraulic conductivity.                  !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function hydr_conduct(k,s,soil_water,soil_fracliq)
      use consts_coms, only : lnexp_min ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: k            ! Layer index                        [  idx]
      integer     , intent(in) :: s            ! Soil texture                       [  idx]
      real(kind=4), intent(in) :: soil_water   ! Soil moisture                      [m3/m3]
      real(kind=4), intent(in) :: soil_fracliq ! Liquid fraction                    [  ---]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=4)             :: relmoist     ! Relative soil moisture             [  ---]
      real(kind=4)             :: fzcorr       ! Freezing correction                [  ---]
      real(kind=4)             :: fibf         ! Incomplete beta function factor    [  ---]
      !------------------------------------------------------------------------------------!


      !------ Find correction for frozen soils. -------------------------------------------!
      fzcorr = exp( max( lnexp_min, - freezecoef * (1.0 - soil_fracliq) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Select approach based on the soil method and soil texture.                     !
      !------------------------------------------------------------------------------------!
      select case (soil(s)%method)
      case ('BDRK')
         !---- Bedrock. Return zero. ------------------------------------------------------!
         hydr_conduct = 0.0
         !---------------------------------------------------------------------------------!
      case ('BC64')
         !---------------------------------------------------------------------------------!
         !     Brooks and Corey (1964).                                                    !
         !---------------------------------------------------------------------------------!


         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water     - soil(s)%soilre )                                    &
                  / ( soil(s)%slmsts - soil(s)%soilre )
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!

         !----- Find the hydraulic conductivity. ------------------------------------------!
         hydr_conduct = fzcorr                                                             &
                      * max( hydcond_min, slcons1(k,s) * relmoist ** soil(s)%slmm )
         !---------------------------------------------------------------------------------!

      case ('vG80')
         !---------------------------------------------------------------------------------!
         !    van Genuchten (1980).                                                        !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water     - soil(s)%soilre )                                    &
                  / ( soil(s)%soilpo - soil(s)%soilre )
         relmoist = max(0.0,min(relmoist,1.0))
         !---------------------------------------------------------------------------------!


         !----- Find the hydraulic conductivity. ------------------------------------------!
         fibf         = 1. - (1. - relmoist**(1./soil(s)%slmm)) ** soil(s)%slmm
         hydr_conduct = slcons1(k,s) * relmoist ** soil(s)%sltt * fibf * fibf
         hydr_conduct = fzcorr * max( hydcond_min, hydr_conduct )
         !---------------------------------------------------------------------------------!
      end select
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
   real(kind=8) function hydr_conduct8(k,s,soil_water,soil_fracliq)
      use consts_coms, only : lnexp_min8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: k            ! Layer index                        [  idx]
      integer     , intent(in) :: s            ! Soil texture                       [  idx]
      real(kind=8), intent(in) :: soil_water   ! Soil moisture                      [m3/m3]
      real(kind=8), intent(in) :: soil_fracliq ! Liquid fraction                    [  ---]
      !----- Internal variables. ----------------------------------------------------------!
      real(kind=8)             :: relmoist     ! Relative soil moisture             [  ---]
      real(kind=8)             :: fzcorr       ! Freezing correction                [  ---]
      real(kind=8)             :: fibf         ! Incomplete beta function factor    [  ---]
      !------------------------------------------------------------------------------------!



      !------ Find correction for frozen soils. -------------------------------------------!
      fzcorr = exp( max( lnexp_min8, - freezecoef8 * (1.d0 - soil_fracliq) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Select approach based on the soil method and soil texture.                     !
      !------------------------------------------------------------------------------------!
      select case (soil8(s)%method)
      case ('BDRK')
         !---- Bedrock. Return zero. ------------------------------------------------------!
         hydr_conduct8 = 0.d0
         !---------------------------------------------------------------------------------!
      case ('BC64')
         !---------------------------------------------------------------------------------!
         !     Brooks and Corey (1964).                                                    !
         !---------------------------------------------------------------------------------!


         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water     - soil8(s)%soilre )                                   &
                  / ( soil(s)%slmsts - soil8(s)%soilre )
         relmoist = max(0.d0,min(relmoist,1.d0))
         !---------------------------------------------------------------------------------!

         !----- Find the hydraulic conductivity. ------------------------------------------!
         hydr_conduct8 = fzcorr                                                            &
                       * max( hydcond_min, slcons18(k,s) * relmoist ** soil8(s)%slmm )
         !---------------------------------------------------------------------------------!

      case ('vG80')
         !---------------------------------------------------------------------------------!
         !    van Genuchten (1980).                                                        !
         !---------------------------------------------------------------------------------!

         !------ Find relative soil moisture. ---------------------------------------------!
         relmoist = ( soil_water      - soil8(s)%soilre )                                  &
                  / ( soil8(s)%soilpo - soil8(s)%soilre )
         relmoist = max(0.d0,min(relmoist,1.d0))
         !---------------------------------------------------------------------------------!


         !----- Find the hydraulic conductivity. ------------------------------------------!
         fibf          = 1.d0 - (1.d0 - relmoist**(1.d0/soil8(s)%slmm)) ** soil8(s)%slmm
         hydr_conduct8 = slcons18(k,s) * relmoist ** soil8(s)%sltt * fibf * fibf
         hydr_conduct8 = fzcorr * max( hydcond_min8, hydr_conduct8 )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end function hydr_conduct8
   !=======================================================================================!
   !=======================================================================================!
end module soil_coms
!==========================================================================================!
!==========================================================================================!
