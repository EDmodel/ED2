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
   !     Constants for new canopy turbulence.  They are used only with Runge-Kutta, so     !
   ! make them double precision.                                                           !
   !---------------------------------------------------------------------------------------!
   !----- Discrete step size in canopy elevation [m]. -------------------------------------!
   real(kind=8),parameter :: dz = 5.d-1
   !----- Dynamic viscosity of air at 20°C [Pa*s]. ----------------------------------------!
   real(kind=8),parameter :: mu0 = 1.827d-5
   !----- Calibration parameter for mu at 20°C [dimensionless]. ---------------------------!
   real(kind=8),parameter :: mu0_c = 1.20d2
   !----- Fluid drag coefficient for turbulent flow in leaves. ----------------------------!
   real(kind=8),parameter :: Cd0 = 2.0d-1
   !----- Sheltering factor of fluid drag on canopies. ------------------------------------!
   real(kind=8),parameter :: Pm = 1.0d0
   !----- Surface drag parameters (Massman 1997). -----------------------------------------!
   real(kind=8),parameter :: c1_m97 = 3.20d-1 
   real(kind=8),parameter :: c2_m97 = 2.64d-1
   real(kind=8),parameter :: c3_m97 = 1.51d1
   !----- Eddy diffusivity due to Von Karman Wakes in gravity flows. ----------------------!
   real(kind=8),parameter :: kvwake = 1.d-3
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
   real(kind=8) :: minimum_canopy_depth

   !=======================================================================================!
   !=======================================================================================!

end module canopy_air_coms
