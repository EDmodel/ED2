!==========================================================================================!
!==========================================================================================!
!   module brams_data.f90.  This module contains the former rams_data common block.        !
!------------------------------------------------------------------------------------------!
module brams_data
   use rpost_dims, only : maxfiles & ! intent(in)
                        , maxgrds  & ! intent(in)
                        , nzpmax   ! ! intent(in)

   real   , dimension(maxfiles)                :: ftimes
   integer, dimension(4,maxgrds,maxfiles)      :: nfgpnts
   integer, dimension(maxfiles)                :: nfgrids
   integer, dimension(maxfiles)                :: ifdates
   integer, dimension(maxfiles)                :: iftimes
   real   , dimension(nzpmax,maxgrds,maxfiles) :: flevels
   real                                        :: startutc
   real                                        :: httop
   real   , dimension(maxgrds,maxfiles)        :: fdelx
   real   , dimension(maxgrds,maxfiles)        :: fdely

end module brams_data
!==========================================================================================!
!==========================================================================================!
