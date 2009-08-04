!==========================================================================================!
!==========================================================================================!
!     Some constants, happening mostly at rk4_derivs.f90.  These used to be at ed_commons, !
! moved it here so they stay together. Still need some explanation of what these variables !
! represent.                                                                               !
!------------------------------------------------------------------------------------------!
module canopy_air_coms

   !=======================================================================================!
   !=======================================================================================!
   !     Parameters used in Euler and Runge-Kutta.                                         !
   !---------------------------------------------------------------------------------------!
   !----- Exponential wind atenuation factor [dimensionless]. -----------------------------!
   real, parameter ::   exar =   2.5
   !----- Scaling factor of Tree Area Index, for computing wtveg [dimensionless]. ---------!
   real, parameter ::   covr =   2.16
   !----- Minimum Ustar [m/s]. ------------------------------------------------------------!
   real, parameter :: ustmin =    .1
   !----- Minimum speed for stars [m/s]. --------------------------------------------------!
   real, parameter ::  ubmin =    .25 
   !----- Double precision version of these variables (for Runge-Kutta). ------------------!
   real(kind=8), parameter ::   exar8 =dble(  exar)
   real(kind=8), parameter :: ustmin8 =dble(ustmin)
   real(kind=8), parameter ::  ubmin8 =dble( ubmin)
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
   !      This is the minimum canopy depth that is used to calculate the heat and moisture !
   ! storage capacity in the canopy air [m].                                               !
   !---------------------------------------------------------------------------------------!
   real         :: minimum_canopy_depth
   real(kind=8) :: minimum_canopy_depth8

   !=======================================================================================!
   !=======================================================================================!

end module canopy_air_coms
