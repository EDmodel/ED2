!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module an_header

type head_table
   character(len=16) :: string
   integer :: npointer,idim_type,ngrid,nvalues
end type

type (head_table), allocatable,save :: anal_table(:)
integer, save:: nvbtab

end module an_header
