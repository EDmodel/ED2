Module soil_coms

  ! DO NOT INITIALIZE PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM
  ! see "initialize_soil_coms" for some default values

  use max_dims,   only: str_len,maxgrds,nzgmax
  use grid_coms,  only: nzg, nzs
#if defined(COUPLED)
  use leaf_coms, only : nstyp,nvtyp,nvtyp_teb
#endif

  implicit none


#if defined(COUPLED)
  integer, parameter :: ed_nstyp = nstyp          ! total # of soil textural classes
  integer, parameter :: ed_nvtyp = nvtyp+nvtyp_teb
#else
  integer, parameter :: ed_nstyp = 12    ! total # of soil textural classes
  integer, parameter :: ed_nvtyp = 21
#endif

  integer :: isoilbc ! Bottom layer boundary condition (0. bedrock or 1. free drainage)

  real, parameter :: soil_rough = 0.05 ! soil roughness height [m]
  real, parameter :: snow_rough = 0.01 ! snowcover roughness height [m]
  real :: water_stab_thresh            ! stability threshold for RK4
                                              ! integrator

  integer, dimension(maxgrds) :: isoilflg
  integer, parameter :: pctlcon = 1
  integer, parameter :: nvgcon = 7 ! I don't think it is been used...

  real, parameter :: min_sfcwater_mass = 1.0e-10 ! minimum allowed mass
  
  
  real, parameter :: snowmin = 3.0 ! Not sure what it means, looks like mass
  ! for sfcwater water.  Quantities below this threshold are ignored.

  integer :: infiltration_method=0
  real :: dewmax
  integer :: nslcon
  real :: zrough
  real, dimension(nzgmax) :: slmstr
  real, dimension(nzgmax) :: stgoff
  real, dimension(nzgmax) :: slz

  !----------------------------------------------------------------------------------!
  !   These are not really soil properties, but since there are just a few of them,  !
  ! they are included in soil_coms                                                   !
  character(len=str_len) :: veg_database
  !----------------------------------------------------------------------------------!
 
  character(len=str_len) :: soil_database

  integer :: isoilstateinit
  integer :: isoildepthflg
  character(len=str_len) :: soilstate_db
  character(len=str_len) :: soildepth_db
  real :: runoff_time

  real :: slden    (ed_nstyp) ! dry soil density [kg/m3]
!  real :: slcpd    (ed_nstyp) ! dry soil volumetric heat capacity [J/(m^3 K)]
!  real :: slbs     (ed_nstyp) ! b exponent [dimensionless]
!  real :: slcons   (ed_nstyp) ! sat soil hydraulic conductivity [m/s]
!  real :: slmsts   (ed_nstyp) ! sat volumetric moist content (porosity) [m^3_wat/m^3_tot]
!  real :: slpots   (ed_nstyp) ! sat moisture potential [m]
!  real :: soilcp   (ed_nstyp) ! minimum soil moisture [m^3_wat/m^3_tot]
  real :: emisg    (ed_nstyp) ! soil infrared emissivity
!  real :: slcons0  (ed_nstyp) ! surface value for slcons [m/s]
  real :: fhydraul (ed_nstyp) ! vertically varying hydraulic conductivity factor
!  real :: xsand    (ed_nstyp) ! soil fractional sand content
!  real :: xclay    (ed_nstyp) ! soil fractional clay content
!  real :: xorgan   (ed_nstyp) ! soil fractional organic content
!  real :: xrobulk  (ed_nstyp) ! soil bulk density [kg/m^3]
!  real :: soilcond0(ed_nstyp) ! New soil cond (8/17/00)
!  real :: soilcond1(ed_nstyp) ! New soil cond (8/17/00)
!  real :: soilcond2(ed_nstyp) ! New soil cond (8/17/00)

  real, allocatable :: slcons1(:,:) ! z-dependent soil sat hydraul cond [m/s]
  real, allocatable :: dslz     (:) ! soil layer thickness at T pt [m]
  real, allocatable :: dslzo2   (:) ! HALF soil layer thickness at T pt [m]
  real, allocatable :: dslzi   (:) ! inverse soil layer thickness at T pt [1/m]
  real, allocatable :: dslzidt  (:) ! dtll / soil layer thickness at T pt [s/m]
  real, allocatable :: slzt     (:) ! soil depth at T pt [m]
  real, allocatable :: dslzt    (:) ! soil layer thickness at M pt [m]
  real, allocatable :: dslzti  (:) ! inverse soil layer thickness at M pt [1/m]
  real, allocatable :: dslztidt (:) ! dtll / soil layer thickness at M pt [s/m]



! Some constants related to soil conductiviness
  real, dimension(6),parameter :: ss=(/1.093e-3,0.028,0.030,3.03e-4,-1.77e-7,2.25e-9/) ! Eq. (E27 from Medvigy 2007

!------------------------------------------------------------------------------------------!
!    Scratch matrix for reading the soil layer and soil depth dataset. This will be        !
! deallocated before the first time step. The number of points is set here for a global    !
! dataset with 1x1 degree resolution.  
!------------------------------------------------------------------------------------------!
  integer, parameter :: nlon_lyr=360, nlat_lyr=180
  integer, allocatable, dimension(:,:) :: layer_index

  type soil_class
!  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke,
!                        1981; Pielke, 1984; Tremback & Kessler, 1985)
     real :: slpots     !  slpots  - saturation moisture potential (m)                                   !
     real :: slmsts     !  slmsts  - saturation volumetric moisture content (m3/m3)                      !
     real :: slbs       !  slbs    - b exponent (dimensionless)                                          !
     real :: slcpd      !  slcons  - saturation soil hydraulic conductivity (m/s)                        !
     real :: soilcp     !  slcons0 - surface value for slcons (m/s)                                      !
     real :: slcons     !  slcpd   - dry soil volumetric heat capacity (J/m3/K)                          !
     real :: slcons0    !  slden   - dry soil density (kg/m3) (also total soil porosity)                 !
     real :: soilcond0
     real :: soilcond1
     real :: soilcond2
     real :: sfldcap
     real :: xsand
     real :: xclay
     real :: xorgan
     real :: xrobulk
     real :: slden
  end type soil_class
  type(soil_class), dimension(ed_nstyp), parameter :: soil =     &

!------------------------------------------------------------------------------------------------------------------------------------------------!
!  slpots  slmsts   slbs    slcpd   soilcp    slcons   slcons0  soilcond0  soilcond1  soilcond2    sfldcap   xsand   xclay  xorgan xrobulk slden !
!------------------------------------------------------------------------------------------------------------------------------------------------!
 (/ soil_class(-0.121,  0.395,  4.05, 1.465e6,  0.0321, 1.760e-4, 5.000e-4, &
 0.30, 4.80, -2.70,  0.135,  .97, .03, .00, 1200., 1600.)  & ! 1. Sand
 ,  soil_class(-0.090,  0.410,  4.38, 1.407e6,  0.0356, 1.563e-4, 6.000e-4, &
 0.30, 4.66, -2.60,  0.150,  .92, .07, .01, 1250., 1600.)  & ! 2. Loamy sand
 ,  soil_class(-0.218,  0.435,  4.90, 1.344e6,  0.0440, 3.467e-5, 7.690e-4, &
 0.29, 4.27, -2.31,  0.195,  .80, .18, .02, 1300., 1600.)  & ! 3. Sandy loam
 ,  soil_class(-0.786,  0.485,  5.30, 1.273e6,  0.0601, 7.200e-6, 1.060e-5, &
 0.27, 3.47, -1.74,  0.255,  .57, .40, .03, 1400., 1600.)  & ! 4. Silt loam
 ,  soil_class(-0.478,  0.451,  5.39, 1.214e6,  0.0580, 6.950e-6, 2.200e-3, &
 0.28, 3.63, -1.85,  0.240,  .60, .35, .05, 1350., 1600.)  & ! 5. Loam
 ,  soil_class(-0.299,  0.420,  7.12, 1.177e6,  0.0545, 6.300e-6, 1.500e-3, &
 0.28, 3.78, -1.96,  0.255,  .65, .31, .04, 1350., 1600.)  & ! 6. Sandy clay loam
 ,  soil_class(-0.356,  0.477,  7.75, 1.319e6,  0.0755, 1.700e-6, 1.070e-4, &
 0.26, 2.73, -1.20,  0.322,  .35, .59, .06, 1500., 1600.)  & ! 7. Silty clay loam
 ,  soil_class(-0.630,  0.476,  8.52, 1.227e6,  0.0664, 2.450e-6, 2.200e-3, &
 0.27, 3.23, -1.56,  0.325,  .48, .45, .07, 1450., 1600.)  & ! 8. Clay loam
 ,  soil_class(-0.153,  0.426, 10.40, 1.177e6,  0.0650, 2.167e-6, 2.167e-6, &
 0.27, 3.32, -1.63,  0.310,  .50, .42, .08, 1450., 1600.)  & ! 9. Sandy clay
 ,  soil_class(-0.490,  0.492, 10.40, 1.151e6,  0.0790, 1.033e-6, 1.033e-6, &
 0.25, 2.58, -1.09,  0.370,  .30, .61, .09, 1650., 1600.)  & !10. Silty clay
 ,  soil_class(-0.405,  0.482, 11.40, 1.088e6,  0.0825, 1.283e-6, 1.283e-6, &
 0.25, 2.40, -0.96,  0.367,  .25, .65, .10, 1700., 1600.)  & !11. Clay
 ,  soil_class(-0.356,  0.863,  7.75, 8.740e5,  0.0860, 8.000e-6, 8.000e-6, &
 0.06, 0.46,  0.00,  0.535,  .20, .20, .60,  500.,  300.) /) !12. Peat
! ,  soil_class(-0.000,  0.000,  0.00, 2.130e6,  0.0000, 0.000000, 0.000000, 4.60, 0.00,  0.00, 1.0e-10, .00, .00, .00,    0. ) /) !13. Bedrock

  integer,dimension(ed_nvtyp)        :: kroot
  
  real,dimension(nzgmax,ed_nvtyp) :: root


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


Contains

  subroutine alloc_soilgrid()

    implicit none
    
    ! Allocate soil grid arrays
    
    allocate (slcons1(nzg,ed_nstyp))
    
    allocate (dslz     (nzg))
    allocate (dslzo2   (nzg))
    allocate (dslzi    (nzg))
    allocate (dslzidt  (nzg))
    allocate (slzt     (nzg))
    allocate (dslzt    (nzg))
    allocate (dslzti   (nzg))
    allocate (dslztidt (nzg))
    
    return
  end subroutine alloc_soilgrid


end Module soil_coms
