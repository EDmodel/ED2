!==========================================================================================!
!==========================================================================================!
!    This module contains the variables used by the LEAF-3 heterogeneous soil moisture     !
! initialisation.  It has the file prefixes and the original soil properties from which    !
! the soil moisture was computed.                                                          !
!------------------------------------------------------------------------------------------!
module mem_soil_moisture
   use leaf_coms, only : nstyp
   implicit none
   character (len=1)   :: soil_moist
   character (len=1)   :: soil_moist_fail
   character (len=256) :: usdata_in
   character (len=256) :: usmodel_in

   real, dimension(nstyp), parameter :: oxsand  = (/ .970, .920, .800, .570, .600, .650    &
                                                   , .350, .480, .500, .300, .250, .200 /)
   real, dimension(nstyp), parameter :: oslmsts = (/ .395, .410, .435, .485, .451, .420    &
                                                   , .477, .476, .426, .492, .482, .863 /)
   real, dimension(nstyp), parameter :: osoilcp = 0.1 - 0.07 * oxsand
end module mem_soil_moisture
!==========================================================================================!
!==========================================================================================!
