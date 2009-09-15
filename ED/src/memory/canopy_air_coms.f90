!==========================================================================================!
!==========================================================================================!
!     Some constants, happening mostly at rk4_derivs.f90.  These used to be at ed_commons, !
! moved it here so they stay together. Still need some explanation of what these variables !
! represent.                                                                               !
!------------------------------------------------------------------------------------------!
module canopy_air_coms

   !=======================================================================================!
   !=======================================================================================!
   !     Parameter that will be read in the namelist.                                      !
   !---------------------------------------------------------------------------------------!
   integer :: icanturb ! Which canopy turbulence we will use?
                       ! -1. is the original ED-2.0 scheme
                       !  0. is the default scheme
                       !  1. will rescale reference wind-speed if it
                       !     the reference height is (inapropriately)
                       !     below the canopy
                       !  2. (recommended) uses the method of Massman 
                       !     1997 and the bulk Richardson number of 
                       !     instability.  This method will not work 
                       !     when zref<h
   integer :: isfclyrm ! Surface layer model (used to compute ustar, tstar,...)
                       !  1 - Louis, 1979: Boundary-Layer Meteor., 17, 187-202.
                       !      This is the ED-2.0 default, also used in (B)RAMS
                       !  2 - Oncley and Dudhia, 1995: Mon. Wea. Rev., 123, 3344-3357.
                       !      This is used in MM5 and WRF.
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!
   !     Parameters used in Euler and Runge-Kutta.                                         !
   !---------------------------------------------------------------------------------------!
   !----- Exponential wind atenuation factor [dimensionless]. -----------------------------!
   real         :: exar
   !----- Scaling factor of Tree Area Index, for computing wtveg [dimensionless]. ---------!
   real         :: covr
   !----- Minimum Ustar [m/s]. ------------------------------------------------------------!
   real         :: ustmin
   !----- Minimum speed for stars [m/s]. --------------------------------------------------!
   real         :: ubmin
   !----- Some parameters that were used in ED-2.0, added here for some tests. ------------!
   real         :: ez
   real         :: vh2vr
   real         :: vh2dh
   !----- Double precision version of some of these variables (for Runge-Kutta). ----------!
   real(kind=8) :: exar8
   real(kind=8) :: ustmin8
   real(kind=8) :: ubmin8
   real(kind=8) :: ez8
   real(kind=8) :: vh2dh8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Constants for new canopy turbulence.                                              !
   !---------------------------------------------------------------------------------------!
   !----- Discrete step size in canopy elevation [m]. -------------------------------------!
   real        , parameter :: dz     = 0.5
   !----- Fluid drag coefficient for turbulent flow in leaves. ----------------------------!
   real        , parameter :: Cd0    = 0.2
   !----- Sheltering factor of fluid drag on canopies. ------------------------------------!
   real        , parameter :: Pm     = 1.0
   !----- Surface drag parameters (Massman 1997). -----------------------------------------!
   real        , parameter :: c1_m97 = 0.320 
   real        , parameter :: c2_m97 = 0.264
   real        , parameter :: c3_m97 = 15.1
   !----- Eddy diffusivity due to Von Karman Wakes in gravity flows. ----------------------!
   real        , parameter :: kvwake = 0.001
   !----- Double precision version of these variables, used in the Runge-Kutta scheme. ----!
   real(kind=8), parameter :: dz8     = dble(dz)
   real(kind=8), parameter :: Cd08    = dble(Cd0)
   real(kind=8), parameter :: Pm8     = dble(Pm)
   real(kind=8), parameter :: c1_m978 = dble(c1_m97)
   real(kind=8), parameter :: c2_m978 = dble(c2_m97)
   real(kind=8), parameter :: c3_m978 = dble(c3_m97)
   real(kind=8), parameter :: kvwake8 = dble(kvwake)
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Constants for surface layer models.                                              !
   !---------------------------------------------------------------------------------------!
   !----- Louis (1979) model. -------------------------------------------------------------!
   real, parameter   :: bl79     = 5.0    ! b prime parameter
   real, parameter   :: csm      = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   real, parameter   :: csh      = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   real, parameter   :: dl79     = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   real, parameter   :: beta     = 5.0    ! Beta used by Businger et al. (1971)
   real, parameter   :: gamm     = 15.0   ! Gamma used by Businger et al. (1971) - momentum.
   real, parameter   :: gamh     = 9.0    ! Gamma used by Businger et al. (1971) - heat.
   real, parameter   :: rri      = 1./.74 ! 1/R value, used by Businger and Louis.
   real, parameter   :: ribmax   = 0.20   ! Maximum bulk Richardson number
   real, parameter   :: tprandtl = 1.00   ! Turbulent Prandtl number.
   !----- Double precision of all these variables. ----------------------------------------!
   real, parameter   :: bl798     = dble(bl79    )
   real, parameter   :: csm8      = dble(csm     )
   real, parameter   :: csh8      = dble(csh     )
   real, parameter   :: dl798     = dble(dl79    )
   real, parameter   :: beta8     = dble(beta    )
   real, parameter   :: gamm8     = dble(gamm    )
   real, parameter   :: gamh8     = dble(gamh    )
   real, parameter   :: rri8      = dble(rri     )
   real, parameter   :: ribmax8   = dble(ribmax  )
   real, parameter   :: tprandtl8 = dble(tprandtl)
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   Tuneable parameters that will be set in ed_params.f90.                              !
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 leaf], so it will be always scaled by LAI.  !
   !---------------------------------------------------------------------------------------!
   real :: dry_veg_lwater
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Maximum leaf water that plants can hold.  Should leaf water exceed this number,    !
   ! water will be no longer intercepted by the leaves, and any value in excess of this    !
   ! will be promptly removed through shedding.  This value is in kg/[m2 leaf], so it will !
   ! be always scaled by LAI.                                                              !
   !---------------------------------------------------------------------------------------!
   real :: fullveg_lwater
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Parameters to set the maximum vegetation-level aerodynamic resistance.             !
   ! rb_max = rb_inter + rb_slope * TAI.                                                   !
   !---------------------------------------------------------------------------------------!
   real :: rb_slope
   real :: rb_inter

   !---------------------------------------------------------------------------------------!
   !      This is the minimum vegetation height.  [m]                                      !
   !---------------------------------------------------------------------------------------!
   real         :: veg_height_min

   !---------------------------------------------------------------------------------------!
   !      This is the minimum canopy depth that is used to calculate the heat and moisture !
   ! storage capacity in the canopy air [m].                                               !
   !---------------------------------------------------------------------------------------!
   real         :: minimum_canopy_depth
   real(kind=8) :: minimum_canopy_depth8

   !=======================================================================================!
   !=======================================================================================!

end module canopy_air_coms
