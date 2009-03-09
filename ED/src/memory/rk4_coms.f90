!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used by the Runge-Kutta integrator scheme.    !
! This module is intended to host the numerical method variables only, variables related   !
! to physical, hydrological or ecological properties should be stored somewhere else.      !
!                                                                                          !
!   MLO. I attempted to describe some variables, not sure about all of them, if the        !
!        definition is wrong, please correct it. Thanks!                                   !
!------------------------------------------------------------------------------------------!
module rk4_coms
   implicit none



   !---------------------------------------------------------------------------------------!
   !    The following variables are parameters that do not depend on the specific run.     !
   !---------------------------------------------------------------------------------------!
   
   !----- Maximum number of intermediate steps --------------------------------------------!
   integer, parameter :: maxstp = 100000000

   !----- Small number, to avoid singularities --------------------------------------------!
   real   , parameter :: tiny_offset = 1.0e-20

   !----- rk4eps is the desired accuracy (former eps), and rk4epsi is its reciprocal ------!
   real   , parameter :: rk4eps  = 0.01
   real   , parameter :: rk4epsi = 1./rk4eps
   
   !----- hmin is the minimum step size ---------------------------------------------------!
   real   , parameter :: hmin = 1.e-9
   !---------------------------------------------------------------------------------------!

   !----- Flag to print the diagnostic check. ---------------------------------------------!
   logical, parameter :: print_diags = .false.

   real   , parameter :: safety = 0.9
   real   , parameter :: pgrow = -0.2
   real   , parameter :: pshrnk = -0.25
   real   , parameter :: errcon = 1.89e-4

   !----- Constants used in rk4_derivs ----------------------------------------------------!
   real   , parameter :: leaf_h2o_thick = 0.11 ! mm
   logical, parameter :: debug  = .true.       ! Verbose output for debug (T|F)
   real   , parameter :: toocold = 193.15      ! Minimum temperature for sat., -80°C
   real   , parameter :: toohot  = 353.15      ! Maximum temperature for sat.,  80°C
   real   , parameter :: lai_to_cover = 1.5    ! Canopies with LAI less than this number 
                                               !    are assumed to be open, ie, some 
                                               !    fraction of the rain-drops can reach
                                               !    the soil/litter layer unimpeded.
   real   , parameter :: evap_area_one = 1.2   ! Evaporation are factor (1 side of leaves 
                                               !    + branches + stems) 
   real   , parameter :: evap_area_two = 2.2   ! Evaporation are factor (2 sides of leaves 
                                               !    + branches + stems) 
   !---------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   !    This is the minimum heat capacity we attempt to prognose leaf internal energy.     !
   ! Below  this value, the fluctuations could become too large.                           !
   !---------------------------------------------------------------------------------------!
   real, parameter :: hcapveg_min=1000.

   !----- These are all RK4 integrator factors. -------------------------------------------!
   real, parameter :: a2  = 0.2
   real, parameter :: a3  = 0.3
   real, parameter :: a4  = 0.6
   real, parameter :: a5  = 1.0
   real, parameter :: a6  = 0.875
   real, parameter :: b21 = 0.2
   real, parameter :: b31 = 3.0/40.0
   real, parameter :: b32 = 9.0/40.0
   real, parameter :: b41 = 0.3
   real, parameter :: b42 = -0.9
   real, parameter :: b43 = 1.2  
   real, parameter :: b51 = -11.0/54.0
   real, parameter :: b52 = 2.5
   real, parameter :: b53 = -70.0/27.0
   real, parameter :: b54 = 35.0/27.0
   real, parameter :: b61 = 1631.0/55296.0
   real, parameter :: b62 = 175.0/512.0
   real, parameter :: b63 = 575.0/13824.0
   real, parameter :: b64 = 44275.0/110592.0
   real, parameter :: b65 = 253.0/4096.0
   real, parameter :: c1  = 37.0/378.0
   real, parameter :: c3  = 250.0/621.0
   real, parameter :: c4  = 125.0/594.0
   real, parameter :: c6  = 512.0/1771.0
   real, parameter :: dc5 = -277.0/14336.0
   real, parameter :: dc1 = c1-2825.0/27648.0
   real, parameter :: dc3 = c3-18575.0/48384.0
   real, parameter :: dc4 = c4-13525.0/55296.0
   real, parameter :: dc6 = c6-0.25
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Integration limits for time. Since 'derivs' do not explicitly depend on time it   !
   ! doesn't really matter what this is as long as tend-tbeg makes sense.                  !
   !---------------------------------------------------------------------------------------!
   real   :: tbeg
   real   :: tend
   real   :: dtrk4
   real   :: dtrk4i
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    The following variables are the bounds of what we could consider barely reasonable !
   ! during the integration process. If values exceed these limits, reject the step right  !
   ! away.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rk4max_can_temp  = 341.00  ! ~10C hotter than record in El Azizia;
   real, parameter :: rk4min_can_temp  = 184.00  ! ~10C colder than record in Vostok;
   real, parameter :: rk4min_can_shv   = 0.0     ! Horribly dry
   real, parameter :: rk4max_can_shv   = 1.0     ! No dry air left, it's water vapour only.

   real, parameter :: rk4max_soil_temp = 351.00  ! ~10C hotter than rk4max_can_temp
   real, parameter :: rk4min_soil_temp = 184.00  ! Same as rk4min_soil_temp

   real, parameter :: rk4max_veg_temp  = 361.00  ! ~20C hotter than rk4max_can_temp. Your
                                                 !      steamed greens are served.

   real, parameter :: rk4min_sfcw_temp = 193.15  ! -80C
   real, parameter :: rk4min_sfcw_mass = -1.e-3  ! ???? Minimum water mass allowed.
   real, parameter :: rk4min_virt_water = -1.e-1 ! Minimum water allowed at virtual pool.
   !---------------------------------------------------------------------------------------!

end module rk4_coms
!==========================================================================================!
!==========================================================================================!
