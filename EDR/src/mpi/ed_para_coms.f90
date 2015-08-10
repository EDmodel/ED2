!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ed_para_coms

  use ed_max_dims, only : maxmach

  !---------------------------------------------------------------------------
  integer                             :: mainnum,nmachs,iparallel,machsize
  integer, dimension(maxmach)         :: machnum
  integer                             :: loadmeth
  !---------------------------------------------------------------------------

end module ed_para_coms
