!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

Module mem_emiss
  use grid_dims, only : str_len
  !---------------------------------------------------------------------------
  integer             :: isource     !Flag for using emission module - EDF
  !This flag is set up in RAMSIN 
  !1=on - 0=off

  !Industrial emissions for NO, NO2, PM25, CO, SO2 and VOC's
  real                :: EINDNO,EINDNO2,EINDPM,EINDCO,EINDSO2,EINDVOC

  !Vehicular emissions for NO, NO2, PM25, CO, SO2 and VOC's
  real                :: EVEINO,EVEINO2,EVEIPM,EVEICO,EVEISO2,EVEIVOC
  
  character(len=3)    :: weekdayin   !Initial weekday of simulation
  !used in TEB and emission module
  !to determinate rates of anthropogenic
  !and emission sources - EDF
  
  real                :: EFSAT,EFSUN !emission fractions for saturdays and 
  ! sundays used in TEB and emission module. 
  
  integer             :: ichemi, &   !for photochemical module activation - EDF
       ichemi_in   !flag for reading a previous run as initial values
  
  character (len=str_len) :: CHEMDATA_IN !path for initial values reading
  !---------------------------------------------------------------------------

end Module
