!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_all

  use mem_basic
  use mem_cuparm
  use mem_grid
  use mem_leaf
  use mem_micro
  use mem_radiate
  use mem_scalar
  use mem_scratch

  use mem_scratch1 ! For reproducibility - Saulo Barros

  use mem_tend
  use mem_turb
  use mem_varinit
  use mem_nestb
  use mem_oda

  use var_tables
  use io_params
  use micphys
  use ref_sounding


end module mem_all
