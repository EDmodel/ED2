!############################# Change Log ##################################
! 1.0.0.1
!
! 000829 CJT an_header ##
!            Added an_header module to replace anal_header.h ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

Module an_header

  type head_table
     character(len=16) :: string
     integer :: npointer,idim_type,ngrid,nvalues
  end type head_table

  type (head_table), allocatable,save :: anal_table(:)
  integer, save:: nvbtab

End Module an_header
