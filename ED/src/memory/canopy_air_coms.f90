!==========================================================================================!
!==========================================================================================!
!     Some constants, happening mostly at rk4_derivs.f90.  These used to be at ed_commons, !
! moved it here so they stay together. Still need some explanation of what these variables !
! represent.                                                                               !
!------------------------------------------------------------------------------------------!
module canopy_air_coms


   real, parameter :: const1 = 116.6    
   real, parameter ::   exar =   2.5    
   real, parameter ::   covr =   2.16   
   real, parameter ::     bz =   0.91   
   real, parameter ::     hz =   0.0075 
   real, parameter ::     ez =   0.172  
   real, parameter :: ustmin =    .1    
   real, parameter ::  ubmin =    .25    ! should use ubmin=1.0 for convec case
  
   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 leaf], so it will be always scaled by LAI.  !
   !---------------------------------------------------------------------------------------!
   real :: min_veg_lwater = 1.e-3
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Maximum leaf water that plants can hold.  Should leaf water exceed this number,    !
   ! water will be no longer intercepted by the leaves, and any value in excess of this    !
   ! will be promptly removed through shedding.  This value is in kg/[m2 leaf], so it will !
   ! be always scaled by LAI.                                                              !
   !---------------------------------------------------------------------------------------!
   real :: max_veg_lwater = 0.11
   !---------------------------------------------------------------------------------------!

end module canopy_air_coms
