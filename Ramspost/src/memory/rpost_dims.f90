!==========================================================================================!
!==========================================================================================!
!   module rpost_dims.f90 (former rconfig.h).  This module contains the main dimensions    !
! sizes for Ramspost.                                                                      !
!------------------------------------------------------------------------------------------!
module rpost_dims

   !---------------------------------------------------------------------------------------!
   !  Set maximum values of parameters:                                                    !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxgrds=4     ! MAXGRDS - Maximum number of grids

   integer, parameter :: nxpmax=300    ! NXPMAX  - Maximum number of points 
                                       !             in x-direction
   integer, parameter :: nypmax=300    ! NYPMAX  - Maximum number of points 
                                       !             in y-direction
   integer, parameter :: nzpmax=200    ! NZPMAX  - Maximum number of points 
                                       !             in z-direction
                                       ! If you change nzpmax, also change nplmax 
                                       ! at ramspost_A.f90
   integer, parameter :: nzgmax=20     ! NZGMAX  - Maximum number of soil levels

   integer, parameter :: maxsclr=21    ! MAXSCLR - Maximum number of additional 
                                       !             scalars
   integer, parameter :: maxhp=1000    ! MAXHP   - Maximum number of u, v, OR t 
                                       !             points in a single vertical
                                       !             level interpolated from 
                                       !             opposite hemispheric grid
   !---------------------------------------------------------------------------------------!

    
    
   !---------------------------------------------------------------------------------------!
   !  Set MAXDIM to the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX.                         !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxdim = max(nxpmax,nypmax,nzpmax+10,nzgmax)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  maxmach is the max number of processors that can be used in a parallel run           !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxmach=1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! maxcloud is the maximum number of clouds                                              !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxclouds=6
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! maxcloud is the maximum number of wave numbers                                        !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: nwave = 50
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! maxpatch is the maximum number of patches                                             !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: maxpatch=25
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     strlen is the typical length of character variables, in particular those used for !
   ! file names.                                                                           !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: str_len = 256
   integer, parameter :: fnm_len = 600
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !        Additional parameters for some array dimensions.                              !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: nxyzpm   = nzpmax*nxpmax*nypmax 
   integer, parameter :: maxdimp  = maxdim+1
   integer, parameter :: nstyp    = 12
   integer, parameter :: nvtyp    = 30
   integer, parameter :: nkeep    = 90
   integer, parameter :: nke      = nkeep
   integer, parameter :: nintgm   = 12
   integer, parameter :: maxsndg  = 200 
   integer, parameter :: maxvarf  = 200
   integer, parameter :: maxsstf  = 200
   integer, parameter :: maxfiles = 2000
   integer, parameter :: nplmax   = nzpmax
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     GrADS maximum dimensions.  Because of the projection from polar-stereographic to  !
   ! regular longitude/latitude, GrADS dimensions must exceed the maximum grid size for a  !
   ! certain amount, given by stfac.                                                       !
   !---------------------------------------------------------------------------------------!
   real   , parameter :: stfac = 1.2
   integer, parameter :: maxgx = ceiling(stfac * nxpmax)
   integer, parameter :: maxgy = ceiling(stfac * nypmax)
   !---------------------------------------------------------------------------------------!

end module rpost_dims
!==========================================================================================!
!==========================================================================================!

