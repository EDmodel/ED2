!==========================================================================================!
!==========================================================================================!
!   module misc_coms.f90.  This module contains a collection of former common blocks that  !
! did not belong to any of the modules in rcommons.h                                       !
!------------------------------------------------------------------------------------------!
module misc_coms
   use rpost_dims, only : nxpmax & ! intent(in)
                        , nypmax ! ! intent(in)
   !----- Former getvar block. ------------------------------------------------------------!
   integer :: ierr_getvar
   integer :: ifound

   !----- Former mem block. ---------------------------------------------------------------!
   integer :: memsize4

   !----- Former grid2 block. -------------------------------------------------------------!
   real, dimension(nxpmax) :: glong
   real, dimension(nypmax) :: glatg

   !----- Former bin block. ---------------------------------------------------------------!
   integer  :: itypp
   integer  :: i0x
   integer  :: i1x
   integer  :: i2x
   integer  :: yoo

   !----- Former nread block. -------------------------------------------------------------!
   integer :: nrflag
   integer :: nfatal


   !----- Former pageout block. -----------------------------------------------------------!
   character(len=132), dimension(80) :: page

   !----- Former vform block. -------------------------------------------------------------!
   character(len=1), dimension(0:63) :: vc


end module misc_coms
!==========================================================================================!
!==========================================================================================!
