! Defines constants used by hydrology routine
! Eventually should be migrated to config file
! left here for development phase to make adding/removing variables easier
module hydrology_constants

  !! Constants for surface runoff
  real, parameter :: m2f= 3.2808399 ! meters to feet
  real, parameter :: c1 = 21.5 !drag coefficient as fcn of water ht (m or ft) (C = c1 + c2*h) 
  real, parameter :: c2 = 3.54 !* m2f

  !! stored variables for output
  real, allocatable :: qw4out(:,:,:) ! flux of water out of a site (kg/s/m2?)
  real, allocatable :: qh4out(:,:,:)

end module hydrology_constants
