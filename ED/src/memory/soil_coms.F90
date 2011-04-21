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
                           , nvtyp     & ! intent(in)
                           , nvtyp_teb ! ! intent(in)
#endif

   implicit none

  !----- These variables depend on whether it's a coupled or stand alone model. -----------!
#if defined(COUPLED)
   integer, parameter :: ed_nstyp = nstyp          ! total # of soil textural classes
   integer, parameter :: ed_nvtyp = nvtyp+nvtyp_teb
#else
   integer, parameter :: ed_nstyp = 12             ! total # of soil textural classes
   integer, parameter :: ed_nvtyp = 21
#endif

   !---------------------------------------------------------------------------------------!
   !    The following variables are assigned through the namelist.                         !
   !---------------------------------------------------------------------------------------!
   integer                                    :: isoilbc        ! Bottom layer bnd. cond.
   integer, dimension(maxgrds)                :: isoilflg       ! Soil initialization flag.
   integer                                    :: nslcon         ! Default soil texture
   real                                       :: slxclay        ! Clay soil fraction
   real                                       :: slxsand        ! Sand soil fraction
   real                                       :: zrough         ! Default soil roughness.
   real, dimension(nzgmax)                    :: slmstr         ! Initial soil moist. frac.
   real, dimension(nzgmax)                    :: stgoff         ! Initial soil temp. offset
   real, dimension(nzgmax)                    :: slz            ! Soil levels.
   character(len=str_len), dimension(maxgrds) :: veg_database   ! Land/sea mask database
   character(len=str_len), dimension(maxgrds) :: soil_database  ! Soil texture database
   character(len=str_len)                     :: soilstate_db   ! Soil state database.
   character(len=str_len)                     :: soildepth_db   ! Soil depth database.
   integer                                    :: isoilstateinit ! Soil state initial cond. 
   integer                                    :: isoildepthflg  ! Soil depth initial cond. 
   real                                       :: runoff_time    ! Default runoff time scale.
   real                                       :: betapower      ! Power for gnd evaporation
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    These following variables will behave as parameters, but they are initialized at   !
   ! init_soil_coms (ed_params.f90), or through the XML config file.                       !
   !---------------------------------------------------------------------------------------!
   real         :: soil_rough          ! soil roughness height                   [       m]
   real         :: snow_rough          ! snowcover roughness height              [       m]
   real(kind=8) :: soil_rough8         ! soil roughness height                   [       m]
   real(kind=8) :: snow_rough8         ! snowcover roughness height              [       m]
   real         :: dewmax              ! Maximum dew flux rate (deprecated)      [ kg/m2/s]
   real         :: water_stab_thresh   ! stability threshold for RK4 integrator  [   kg/m2]
   real         :: tiny_sfcwater_mass  ! Min. mass allowed in temporary layers   [   kg/m2]
   real         :: snowmin             ! Min. snow mass needed to create new lyr [   kg/m2]
   integer      :: infiltration_method ! Infiltration scheme (for rk4_derivs)    [     0|1]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following variables are assigned in ed_params.f90 based on namelist           !
   ! variables.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) :: betapower8 ! Power for ground evaporation
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
   real   , dimension(ed_nstyp)        :: emisg    ! soil infrared emissivity
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
      real :: slpots             ! Soil moisture potential at saturation         [       m]
      real :: slmsts             ! Soil moisture at saturation                   [   m3/m3]
      real :: slbs               ! B exponent                                    [     n/d]
      real :: slcpd              ! Specific heat of dry soil                     [  J/m3/K]
      real :: soilcp             ! Dry soil capacity (at -3.1MPa)                [   m3/m3]
      real :: soilwp             ! Wilting point capacity (at -1.5MPa)           [   m3/m3]
      real :: slcons             ! hydraulic conductivity at saturation          [     m/s]
      real :: slcons0            ! Surface value for slcons                      [     m/s]
      real :: soilcond0          ! Intercept for conductivity calculation        [   N/K/s]
      real :: soilcond1          ! Linear coefficient for conductivity           [   N/K/s]
      real :: soilcond2          ! Quadratic coefficient for conductivity        [   N/K/s]
      real :: sfldcap            ! Soil field capacity                           [   m3/m3]
      real :: xsand              ! Percentage of sand                            [     ---]
      real :: xclay              ! Percentage of clay                            [     ---]
      real :: xorgan             ! Percentage of organic material                [     ---]
      real :: xrobulk            ! Bulk density                                  [     ---]
      real :: slden              ! "Dry" soil density (porosity)                 [   kg/m3]
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
      real(kind=8) :: soilcond0  ! Intercept for conductivity calculation        [   N/K/s]
      real(kind=8) :: soilcond1  ! Linear coefficient for conductivity           [   N/K/s]
      real(kind=8) :: soilcond2  ! Quadratic coefficient for conductivity        [   N/K/s]
      real(kind=8) :: sfldcap    ! Soil field capacity                           [   m3/m3]
      real(kind=8) :: xsand      ! Percentage of sand                            [     ---]
      real(kind=8) :: xclay      ! Percentage of clay                            [     ---]
      real(kind=8) :: xorgan     ! Percentage of organic material                [     ---]
      real(kind=8) :: xrobulk    ! Bulk density                                  [     ---]
      real(kind=8) :: slden      ! "Dry" soil density (porosity)                 [   kg/m3]
   end type soil_class8
   !---------------------------------------------------------------------------------------!
   !----- To be filled in ed_params.f90. --------------------------------------------------!
   type(soil_class8), dimension(ed_nstyp)            :: soil8 
   type(soil_class) , dimension(ed_nstyp)            :: soil


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


      allocate (slcons1(nzg,ed_nstyp))

      allocate (slz8     (nzg))
      allocate (dslz     (nzg))
      allocate (dslzo2   (nzg))
      allocate (dslzi    (nzg))
      allocate (dslzidt  (nzg))
      allocate (slzt     (nzg))
      allocate (dslzt    (nzg))
      allocate (dslzti   (nzg))
      allocate (dslztidt (nzg))

      allocate (slcons18(nzg,ed_nstyp))

      allocate (dslz8     (nzg))
      allocate (dslzo28   (nzg))
      allocate (dslzi8    (nzg))
      allocate (dslzidt8  (nzg))
      allocate (slzt8     (nzg))
      allocate (dslzt8    (nzg))
      allocate (dslzti8   (nzg))
      allocate (dslztidt8 (nzg))

      return
   end subroutine alloc_soilgrid
   !=======================================================================================!
   !=======================================================================================!

end Module soil_coms
!==========================================================================================!
!==========================================================================================!
