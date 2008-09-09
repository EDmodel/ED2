!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module grid_dims

  ! This module contains very basic specification of grid dimensions and other 
  ! parameters that will be used to dimension arrays and allocate memory.

  ! Logo BRAMS version

  character(len=*), parameter :: BRAMS_version="BRAMS Version 4.0"

  ! Grid dimensions:
  !   MAXGRDS - Maximum number of grids
  !   NXPMAX  - Maximum number of points in x-direction
  !   NYPMAX  - Maximum number of points in y-direction
  !   NZPMAX  - Maximum number of points in z-direction
  !   NZGMAX  - Maximum number of soil levels
  !   MAXSCLR - Maximum number of additional scalars
  !   MAXHP   - Maximum number of u, v, OR t points in a single vertical
  !             level interpolated from opposite hemispheric grid
  !   MAXDIM  - the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX

  integer, parameter :: maxgrds=8
  integer, parameter :: nxpmax=303
  integer, parameter :: nypmax=303
  integer, parameter :: nzpmax=132
  integer, parameter :: nzgmax=20
  integer, parameter :: maxsclr=150
  integer, parameter :: maxhp=1000
  integer, parameter :: maxdim=303

  ! Computed parameters (function of previous parameters)

  integer, parameter :: maxdimp=maxdim+2
  integer, parameter :: nxyzpm=nzpmax*nxpmax*nypmax

  !   MAXMACH - the maximum number of processors on a parallel run

  integer, parameter :: maxmach=64

  ! TEB
  !  MAXSTEB - Maximum number of layers used in TEB
  !  MAXUBTP - Maximum number of urban types used in TEB
  integer, parameter :: maxsteb=5
  integer, parameter :: maxubtp=3

end module grid_dims
