!======================================= Change Log =======================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
! Module grid_dims: this module contains very basic specification of grid dimensions and   !
!                   other parameters that will be used to dimension arrays and allocate    !
!                   memory.                                                                !
!------------------------------------------------------------------------------------------!
module grid_dims

   implicit none


   !---- BRAMS version. -------------------------------------------------------------------!
   character(len=19), parameter :: BRAMS_version='BRAMS Version 4.0.6'

   !---------------------------------------------------------------------------------------!
   ! Grid dimensions:                                                                      !
   !                                                                                       !
   !   MAXGRDS - Maximum number of grids                                                   !
   !   NXPMAX  - Maximum number of points in x-direction                                   !
   !   NYPMAX  - Maximum number of points in y-direction                                   !
   !   NZPMAX  - Maximum number of points in z-direction                                   !
   !   NZGMAX  - Maximum number of soil levels                                             !
   !   MAXSCLR - Maximum number of additional scalars                                      !
   !   MAXHP   - Maximum number of u, v, OR t points in a single vertical                  !
   !             level interpolated from opposite hemispheric grid                         !
   !   MAXDIM  - the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX                             !
   !   MAXDIMP - MAXDIM + 2                                                                !
   !   NXYZPM  - Maximum number of volume points                                           !
   !   MAXMACH - Maximum number of cores on a parallel run.                                !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxgrds =    8
   integer, parameter :: nxpmax  =  600
   integer, parameter :: nypmax  =  600
   integer, parameter :: nzpmax  =  132
   integer, parameter :: nzgmax  =  100
   integer, parameter :: maxsclr =  150
   integer, parameter :: maxhp   = 1000
   integer, parameter :: maxdim  = max(nxpmax,nypmax,nzpmax+10,nzgmax)
   integer, parameter :: maxdimp = maxdim+2
   integer, parameter :: nxyzpm  = nzpmax*nxpmax*nypmax
   integer, parameter :: maxmach = 128
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !  TEB                                                                                  !
   !                                                                                       !
   !  MAXSTEB - Maximum number of layers used in TEB                                       !
   !  MAXUBTP - Maximum number of urban types used in TEB                                  !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxsteb   = 5
   integer, parameter :: maxubtp   = 3
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Grell cumulus parametrisation.                                                       !
   !                                                                                       !
   !  MAXCLOUDS - maximum number of cloud spectral sizes allowed.                          !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxclouds = 6
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Output                                                                               !
   !                                                                                       !
   !  NDIM_TYPES - Number of different dimensions that output variables may have.          !
   !               (2-D, 3-D, surface patch-level, soil patch-level, etc...)               !
   !  NUMBER_DIMS - Number of dimensions for each of the ndim_types.                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   integer                       , parameter :: ndim_types = 9
   integer, dimension(ndim_types), parameter :: number_dims = (/ 1, 2, 3, 4, 4             &
                                                               , 3, 3, 4, 3 /)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  File management
   !
   !  STR_LEN  - String length for many file variables
   !  MAXFILES - Maximum number of files that may be loaded
   !---------------------------------------------------------------------------------------!
   integer, parameter :: str_len  = 300
   integer, parameter :: maxfiles = 36000
   !---------------------------------------------------------------------------------------!

end module grid_dims
!==========================================================================================!
!==========================================================================================!
