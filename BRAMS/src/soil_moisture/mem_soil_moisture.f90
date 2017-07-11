!==========================================================================================!
!==========================================================================================!
!    This module contains the variables used by the LEAF-3 heterogeneous soil moisture     !
! initialisation.  It has the file prefixes and the original soil properties from which    !
! the soil moisture was computed.                                                          !
!------------------------------------------------------------------------------------------!
module mem_soil_moisture
   use leaf_coms, only : nstyp
   use grid_dims, only : str_len
   implicit none
   character (len=1)       :: soil_moist
   character (len=1)       :: soil_moist_fail
   character (len=str_len) :: usdata_in
   character (len=str_len) :: usmodel_in

   real, dimension(nstyp), parameter :: oxsand  = (/ .970, .920, .800, .570, .600, .650    &
                                                   , .350, .480, .500, .300, .250, .200    &
                                                   , .000, .570, .250 ,.250, .250 /)
   real, dimension(nstyp), parameter :: oslmsts = (/ .395, .410, .435, .485, .451, .420    &
                                                   , .477, .476, .426, .492, .482, .863    &
                                                   , .000, .485, .482, .482, .482 /)
   real, dimension(nstyp), parameter :: oslpots = (/-.121,-.090,-.218,-.786,-.478,-.299    &
                                                   ,-.356,-.630,-.153,-.490,-.405,-.356    &
                                                   , .000,-.786,-.405,-.405,-.405 /)
   real, dimension(nstyp), parameter :: oslbs   = (/ 4.05, 4.38, 4.90, 5.30, 5.39, 7.12    &
                                                   , 7.75, 8.52,10.40,10.40,11.40, 7.75    &
                                                   , 0.00, 5.30,11.40,11.40,11.40 /)
   real, dimension(nstyp), parameter :: osoilcp = 0.1 - 0.07 * oxsand
end module mem_soil_moisture
!==========================================================================================!
!==========================================================================================!
