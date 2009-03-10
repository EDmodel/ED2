module canopy_air_coms

! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM
  
! Some constants, happening mostly at rk4_derivs.f90. 
! These used to be at ed_commons, moved it here so they stay altogether, but 
! it still needs some explanation of what these values are.
  real, parameter :: const1 = 116.6    
  real, parameter ::   exar =   2.5    
  real, parameter ::   covr =   2.16   
  real, parameter ::     bz =   0.91   
  real, parameter ::     hz =   0.0075 
  real, parameter ::     ez =   0.172  
  real, parameter :: ustmin =    .1    
  real, parameter ::  ubmin =    .25    ! should use ubmin=1.0 for convec case
  

end module canopy_air_coms
