!==========================================================================================!
!==========================================================================================!
!  RAPP. module netcdf: contains variables used at netCDF I/O handling.                    !
!------------------------------------------------------------------------------------------!
module mod_netcdf
   
#if USE_NCDF
   use mod_maxdims, only : maxstr
   use mod_time, only : time_stt
   use netcdf
   implicit none


   integer :: ncid            ! NetCDF file ID #
   integer :: ndimensions     ! Number of dimensions defubed for this netCDF dataset
   integer :: nvariables      ! Number of variables defined for this netCDF dataset
   integer :: nglobals        ! Number of global attributes for this netCDF dataset
   integer :: unlimiteddimid  ! ID of the unlimited dimension, if there is one for this
                              !   netCDF dataset. If not, then it will return -1.

   !------ Dimensions for time variable and global attributes -----------------------------!
   integer                                                    :: timeid      ! Time var ID
   character(len=NF90_MAX_NAME)                               :: dummy_vname ! Scratch
   integer                                                    :: dimglobal
   integer                                                    :: globalid
   integer                                                    :: dimid
   !------ Dimensions for variables in general --------------------------------------------!
   integer                                                    :: varid  ! Variable ID
   integer                                                    :: ndims  ! # of dimensions
   integer                                                    :: xtype  ! Data type
   integer                     , dimension(NF90_MAX_VAR_DIMS) :: dimids ! Dimension ID
   integer                                                    :: natts  ! # of attributes
   !------ ID for useful dimensions, this is just to decide the variable type np-----------!
   integer                                                    :: idnntp
   integer                                                    :: idnnxp  
   integer                                                    :: idnnyp  
   integer                                                    :: idnnzp  
   integer                                                    :: idnzg   
   integer                                                    :: idnzs
   integer                                                    :: idnclouds   
   integer                                                    :: idnpatch
   integer                                                    :: idnwave
   integer                                                    :: idnnxpst
   integer                                                    :: idnnypst
   integer                                                    :: idnnzpst
   integer                                                    :: idtimelen

#endif


end module mod_netcdf
!==========================================================================================!
!==========================================================================================!
