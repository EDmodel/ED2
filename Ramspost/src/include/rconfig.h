!f90
!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!---------------------------------------------------------------------------
!  Set maximum values of parameters:

integer, parameter :: MAXGRDS=4     ! MAXGRDS - Maximum number of grids

integer, parameter :: NXPMAX=300    ! NXPMAX  - Maximum number of points 
                                    !             in x-direction
integer, parameter :: NYPMAX=300    ! NYPMAX  - Maximum number of points 
                                    !             in y-direction
integer, parameter :: NZPMAX=200    ! NZPMAX  - Maximum number of points 
                                    !             in z-direction
                                    ! If you change nzpmax, also change nplmax 
                                    ! at ramspost_A.f90
integer, parameter :: NZGMAX=20     ! NZGMAX  - Maximum number of soil levels

integer, parameter :: MAXSCLR=50    ! MAXSCLR - Maximum number of additional 
                                    !             scalars
integer, parameter :: MAXHP=1000    ! MAXHP   - Maximum number of u, v, OR t 
                                    !             points in a single vertical
                                    !             level interpolated from 
                                    !             opposite hemispheric grid
                                    
!---------------------------------------------------------------------------
!  Set MAXDIM to the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX

integer, parameter :: MAXDIM=300

!---------------------------------------------------------------------------
!  maxmach is the max number of processors that can be used in a parallel run

integer, parameter :: maxmach=1

!---------------------------------------------------------------------------
! maxcloud is the maximum number of clouds
integer, parameter :: maxclouds=6
