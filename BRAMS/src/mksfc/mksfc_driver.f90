!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!===== TAKE INTO ACCOUNT ALL OPTIONS, HISTORY STARTS, ADDED GRIDS, ETC.

subroutine make_sfcfiles()

  use mem_mksfc, only: &
       sfcfile_p,      &
       scr1,           &
       scr2,           &
       vt2da,          &
       vt2db,          &
       scrx,           &
       npq,            &
       glatp,          &
       glonp,          &
       datq_patch,     &
       datp,           &
       nvndvif,        &
       nvsstf,         &
       glatp,          &
       alloc_sfcfile,  &
       dealloc_sfcfile

  use mem_grid, only:  &
       ngrids,         &
       nnxp,           &
       nnyp,           &
       nnxyp,          &
       npatch,         &
       nzg,            &
       ngridsh,        &
       nxtnest,        &
       deltaxn,        &
       deltayn,        &
       platn,          &
       plonn,          &
       xtn,            &
       ytn,            &
       runtype

  use io_params, only: &
       nxpmax,         &
       nypmax,         &
       iupdsst,        &
       isstflg,        &
       itoptflg,       &
       iupdndvi,       &
       ndvifn,         &
       ndviflg,        &
       isoilfn,        &
       isoilflg,       &
       ivegtfn,        &
       ivegtflg

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM !INTENT(IN)
  !

  implicit none

  integer :: ifm,icm,nvtime,ivtime,ng1,ng2,ng1t,ng2t,ng1s,ng2s
  integer :: isfcerr,itoperr,issterr,indvierr

  ! TEB_SPM
  integer :: ng1f, ng2f
  integer :: ifusoerr
  !

  ! This subroutine makes sure that all surface, topo, sst, and ndvi
  ! files required for the present run exist and are correct. 

  ! The logical choices made are as follows:

  ! For runtype = 'MAKESFC':
  !    Make all surface, topo, sst, and ndvi files for grids 1:NGRIDS.
  ! For runtype = 'INITIAL' or 'MAKEVFILE':
  !    Check for existence and correctness of all surface, topo, sst,
  !       and ndvi files for grids 1:NGRIDS.  Remake all grids if the files
  !       are incorrect for the set of files: topo, sfc/ndvi, sst
  ! For runtype = 'HISTORY':
  !    If NGRIDS > NGRIDSH, check for correctness of all sfcfiles
  !       for grids NGRIDSH+1:NGRIDS.  For any grid (ifm) with an 
  !       incorrect sfcfile, remake the file, requiring that
  !       itoptflg(ifm) = 0.
  !    Check for existence and correctness of all sst and ndvi 
  !       files for grids 1:NGRIDS.  For any grid (ifm) with
  !       incorrect files, make sst and ndvi files.



  isfcerr = 0
  itoperr = 0
  issterr = 0
  indvierr = 0
  ! TEB_SPM
  ifusoerr = 0
  !

  ! Allocate memory needed for initializing sfcfiles

  !srf : para Itanium II com runtype=history
  !allocate( sfcfile_p(ngrids) )
  if (allocated(sfcfile_p)) then 
     deallocate(sfcfile_p)
  end if
  allocate(sfcfile_p(ngrids))
  do ifm = 1,ngrids
     call alloc_sfcfile(sfcfile_p(ifm),nnxp(ifm),nnyp(ifm),nzg,npatch)
  end do

  allocate (scr1 (nxpmax,nypmax))
  allocate (scr2 (nxpmax,nypmax))
  allocate (vt2da(nxpmax,nypmax))
  allocate (vt2db(nxpmax,nypmax))
  allocate (scrx(maxval(nnxyp(1:ifm))*nzg) )

  if (runtype(1:7) == 'MAKESFC') then

     itoperr = 1
     isfcerr = 1
     issterr = 1
     indvierr = 1

     ng1=1 ; ng2=ngrids      ! sst,ndvi grid bounds
     ng1t=1 ; ng2t=ngrids    ! topo grid bounds
     ng1s=1 ; ng2s=ngrids    ! sfc grid bounds

     if (TEB_SPM==1) then
        ifusoerr = 1
        ng1f=1 ; ng2f=ngrids    ! fuso grid bounds
     end if

  elseif (runtype(1:7) == 'INITIAL' .or. runtype(1:9) == 'MAKEVFILE') then

     ! Check sfc files 
     do ifm = 1,ngrids
        call sfc_check(ifm,isfcerr)
        if(isfcerr == 1) exit
     end do

     ! Check topo files 
     do ifm = 1,ngrids
        call top_check(ifm,itoperr)
        if(itoperr == 1) exit
     end do

     ! Check sst files
     do ifm = 1,ngrids
        call sst_read(2,ifm,issterr)
        if(issterr == 1) exit
     end do

     ! Check ndvi files
     do ifm = 1,ngrids
        if(ndviflg(ifm) >= 0) call ndvi_read(2,ifm,indvierr)
        if(indvierr == 1) exit
     end do

     if (TEB_SPM==1) then
        ! Check fuso files 
        do ifm = 1,ngrids
           call fuso_check(ifm,ifusoerr)
           if(ifusoerr == 1) exit
        end do
     end if

     ! If we are making ndvi files, we must also make the sfc files (and vice versa)
     if(indvierr == 1) isfcerr = 1
     if(isfcerr == 1) indvierr = 1

     if (isfcerr==0 .and. issterr==0 .and.  &
          itoperr==0 .and.indvierr==0) then
        if (TEB_SPM==1) then
           if (ifusoerr==0) print*, '   TEB_SPM: fuso files ok'
        end if
        return
     else
        print*, 'Nonexistent or incorrect surface files for'
        print*, '   INITIAL or MAKEVFILE runtype (re)making'
        if(isfcerr == 1 ) print*, '   sfc files'
        if(itoperr == 1 ) print*, '   top files'
        if(issterr == 1 ) print*, '   sst files'
        if(indvierr == 1) print*, '   ndvi files'
        if (TEB_SPM==1) then
           if(ifusoerr == 1) print*, '   fuso files'
        end if
     end if

     ng1=1 ; ng2=ngrids
     ng1t=1 ; ng2t=ngrids
     ng1s=1 ; ng2s=ngrids
     if (TEB_SPM==1) then
        ng1f=1 ; ng2f=ngrids
     end if

  elseif (runtype(1:7) == 'HISTORY') then

     !   We can't do this set of checks until we have done a full history start
     !     since we may have to interpolate added grids from coarse grid
     !     fields. This is call from INITLZ on a history start.

     ! Check topo files for added grids. If these are not correct, only
     !   allow remakes if topo is to be interpolated from existing grids.

     do ifm = ngridsh+1,ngrids
        call top_check(ifm,itoperr)
        if (itoperr == 1) then
           if (itoptflg(ifm) /= 0) then
              print*, 'Nonexistent or incorrect TOPFILE for grid ',ifm
              print*, '   which is being added on a history start:'
              print*, '   ITOPTFLG must be set to zero in this case.'
              stop 'added grid - itoptflg'
           end if
        end if
     end do
     ng1t=ngridsh+1 ; ng2t=ngrids

     ! Check sfc files for added grids. Existing grid info is read from
     !   history file.
     do ifm = ngridsh+1,ngrids
        call sfc_check(ifm,isfcerr)
        if(isfcerr == 1) exit
     end do
     ng1s=ngridsh+1 ; ng2s=ngrids

     if (TEB_SPM==1) then
        ! Check fuso files for added grids. Existing grid info is read from
        !   history file.
        do ifm = ngridsh+1,ngrids
           call fuso_check(ifm,ifusoerr)
           if(ifusoerr == 1) exit
        end do
        ng1f=ngridsh+1 ; ng2f=ngrids
     end if

     ! Check sst files for all grids. This is a potentially time-dependent field,
     !   so we need to check if all files are there. If there is a set of files
     !   for a grid that is incomplete, remake all files.

     !   This is somewhat dangerous, since the namelist could have changed 
     !   since the previous run to
     !   specify different data files and we really have no way of 
     !   knowing. We will make the assumption that this didn't occur.
     do ifm = 1,ngrids
        call sst_read(2,ifm,issterr)
        if(issterr == 1) exit
     end do

     ! Do same for NDVI files
     do ifm = 1,ngrids
        if(ndviflg(ifm) >= 0) call ndvi_read(2,ifm,indvierr)
        if(indvierr == 1) exit
     end do

     ! If we are making ndvi files, we must also make the sfc files. We will
     !    remake all grids since we may be interpolating from coarser grid.
     if(indvierr == 1 .or. isfcerr == 1) ng1s=1

     ng1=1 ; ng2=ngrids

  end if


  ! If we got here, at least one set of files are bad. Re-make the bad ones.


  !------------------------------------------
  !  TOP (topo and roughness) file creation
  if(itoperr == 1) then
     ! do topography, topo roughness on all grids
     call toptnest(ng1t,ng2t)
     do ifm = 1,ngrids
        call top_write(ifm)
     end do
  end if

  ! TEB_SPM
  if (TEB_SPM==1) then
     !------------------------------------------
     !  FUSO  file creation
     if(ifusoerr == 1) then
        ! do FUSO (Local Time) on all grids
        call fusonest(ng1f,ng2f)
        do ifm = 1,ngrids
           call fuso_write(ifm)
        end do
     end if
  end if

  !------------------------------------------
  !  SFC (veg class, patch area, soil type) and NDVI file creation
  !      If we are making ndvi files, we must also make the sfc files (and vice versa)
  if(isfcerr == 1 .or. indvierr == 1) then

     ! If iupdndvi = 1, require that:
     !    (1) ndviflg = 1 for a grid that has nxtnest = 0,
     !    (2) ndviflg /= 2 for all other grids.

     if (iupdndvi == 1) then
        do ifm = 1,ngrids
           if (ndviflg(ifm) /= 1 .and. nxtnest(ifm) == 0) then
              print*, 'iupdndvi = 1 and ndviflg /= 1 for grid ', ifm
              stop 'iupdndvi'
           end if

           if (ndviflg(ifm) == 2) then
              print*, 'iupdndvi = 1 and ndviflg = 2 for grid ', ifm
              stop 'iupdndvi'
           end if
        end do
     end if


     do ifm = ng1s,ng2s

        if (ivegtflg(ifm) == 1 .or. isoilflg(ifm) == 1 .or. ndviflg(ifm) == 1) then   

           ! Find size of patch arrays
           call patch_array_size(npq,(xtn(2,ifm)-xtn(1,ifm))  &
                ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm)  &
                ,ndviflg(ifm),ndvifn(ifm) )

           ! Allocate arrays that need npq dimension
           allocate (glatp(npq,npq,nnxp(ifm),nnyp(ifm))  &
                ,glonp(npq,npq,nnxp(ifm),nnyp(ifm))  &
                ,datq_patch(npq,npq,nnxp(ifm),nnyp(ifm))  &
                ,datp(npq,npq,nnxp(ifm),nnyp(ifm)) )

           ! Fill lat-lon
           call patch_latlon(nnxp(ifm),nnyp(ifm)  &
                ,xtn(1,ifm),ytn(1,ifm),deltaxn(ifm),deltayn(ifm)  &
                ,platn(ifm),plonn(ifm) ) 
        end if


        ! do sfcfile
        call geonest_file(ifm)
        call sfc_write(ifm)


        ! do ndvifile
        if (ndviflg(ifm) == 1) then
           call ndvi_read_dataheader(ifm)
           nvtime = nvndvif(ifm)
        elseif (isstflg(ifm) == 0) then
           nvtime = nvndvif(nxtnest(ifm))
        else
           nvtime = 1
        end if

        do ivtime = 1,nvtime
           call ndvinest (ifm,ivtime)
           call ndvi_write (ifm,ivtime)
        end do

        ! Deallocate arrays that needed npq dimension
        if(allocated(glatp)) deallocate (glatp,glonp,datq_patch,datp)

     end do

  end if

  !------------------------------------------
  !  SST file creation
  if(issterr == 1) then

     ! If iupdsst = 1, require that:
     !    (1) isstflg = 1 for a grid that has nxtnest = 0,
     !    (2) isstflg /= 2 for all other grids.

     if (iupdsst == 1) then
        do ifm = 1,ngrids
           if (isstflg(ifm) /= 1 .and. nxtnest(ifm) == 0) then
              print*, 'iupdsst = 1 and isstflg /= 1 for grid ', ifm
              stop 'iupdsst'
           end if

           if (isstflg(ifm) == 2) then
              print*, 'iupdsst = 1 and isstflg = 2 for grid ', ifm
              stop 'iupdsst'
           end if
        end do
     end if

     do ifm = ng1,ng2
        ! do sstfile
        if (isstflg(ifm) == 1) then
           call sst_read_dataheader(ifm)
           nvtime = nvsstf(ifm)
        elseif (isstflg(ifm) == 0) then
           nvtime = nvsstf(nxtnest(ifm))
        else
           nvtime = 1
        end if

        do ivtime = 1,nvtime
           call sstnest (ifm,ivtime)
           call sst_write (ifm,ivtime)
        end do
     end do

  end if


  ! Deallocate memory needed for initializing sfcfiles

  do ifm = 1,ngrids
     call dealloc_sfcfile(sfcfile_p(ifm))
  end do

  deallocate (sfcfile_p)

  deallocate (scr1,scr2,vt2da,vt2db,scrx)

end subroutine make_sfcfiles
