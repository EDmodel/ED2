!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_mksfc


  Type sfcfile_vars

     !(nxp,nyp,nzg,npatch)
     real, pointer, dimension(:,:,:,:) :: soil_text

     !(nxp,nyp,npatch)
     real, pointer, dimension(:,:,:) :: patch_area,leaf_class,veg_ndvif
     
     !(nxp,nyp)
     real, pointer, dimension(:,:) :: topt, seatf, topzo
     
     ! TEB_SPM
     real, pointer, dimension(:,:) :: fuso
     
  End Type sfcfile_vars
   

  type (sfcfile_vars), allocatable :: sfcfile_p(:)
  
  !(nxpmax,nypmax)
  real, allocatable, dimension(:,:) :: scr1,scr2,vt2da,vt2db
  
  !(max(nxp*nyp)*nzg)
  real, allocatable, dimension(:) :: scrx
  
  !(np,np,nxp,nyp)
  real, allocatable, dimension(:,:,:,:) :: glatp,glonp,datp

  !(np,np,nxp,nyp)
  integer, allocatable, dimension(:,:,:,:) :: datq_patch
  
  !(np*np*nxp*nyp)
  integer, allocatable, dimension(:) :: ptable
  
  !(iblksizo,iblksizo)
  real, allocatable, dimension(:,:) :: dato
  
  !(iblksizo,iblksizo)
  character(len=1), allocatable, dimension(:,:) :: cdato
  
  !(ifile_max,jfile_max)
  integer, allocatable, dimension(:,:) :: nump,numpind,numpind1,numpind2
  
  integer :: npq
  
  integer, parameter :: maxsfcgrids=10
  
  ! SST file creation variables
  integer, parameter :: maxsstdata=100
  integer, dimension(maxsstdata,maxsfcgrids):: iyearvs,imonthvs,idatevs,ihourvs
  integer,dimension(maxsfcgrids)         :: nvsstf
  character(len=128), dimension(maxsstdata,maxsfcgrids)     :: vsstfil
  
  ! NDVI file creation variables
  integer, parameter :: maxndvidata=100
  integer, dimension(maxndvidata,maxsfcgrids):: iyearvn,imonthvn,idatevn,ihourvn
  integer,dimension(maxsfcgrids)         :: nvndvif
  character(len=128), dimension(maxndvidata,maxsfcgrids)     :: vndvifil
  
Contains

  subroutine alloc_sfcfile(sfcfile,nx,ny,nzg,npat)

    use teb_spm_start, only: TEB_SPM ! INTENT(IN)

    implicit none

    type (sfcfile_vars) :: sfcfile
    integer, intent(in) :: nx,ny,nzg,npat
     
    allocate (sfcfile%soil_text    (nzg,nx,ny,npat))
    
    allocate (sfcfile%patch_area   (nx,ny,npat))
    allocate (sfcfile%leaf_class  (nx,ny,npat))
    allocate (sfcfile%veg_ndvif    (nx,ny,npat))
    
    allocate (sfcfile%topt         (nx,ny))
    allocate (sfcfile%seatf        (nx,ny))
    allocate (sfcfile%topzo        (nx,ny))
    ! TEB_SPM
    if (TEB_SPM==1) then
       allocate (sfcfile%fuso         (nx,ny))
    endif
    
  end subroutine alloc_sfcfile
  
  ! ******************************************************************
  
  subroutine dealloc_sfcfile(sfcfile)
    
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)
    
    implicit none
    
    type (sfcfile_vars) :: sfcfile
    
    deallocate (sfcfile%soil_text)
    
    deallocate (sfcfile%patch_area)
    deallocate (sfcfile%leaf_class)
    deallocate (sfcfile%veg_ndvif)
    
    deallocate (sfcfile%topt)
    deallocate (sfcfile%seatf)
    deallocate (sfcfile%topzo)
    ! TEB_SPM
    if (TEB_SPM==1) then
       deallocate (sfcfile%fuso)
    endif
    
  end subroutine dealloc_sfcfile
  
End Module mem_mksfc
