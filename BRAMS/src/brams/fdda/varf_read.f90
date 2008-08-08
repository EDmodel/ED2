!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine varf_read(ivflag)

use mem_grid
use mem_varinit

implicit none

integer :: ivflag

character(len=14)  :: itotdate_current
integer :: iyears,imonths,idates,ihours,nf,ifm  &
          ,ivar_wait,nwaits,nw,ivwait1,irsleep,irslp,ifileok,icm

! See if we want to possibly wait for files to be available.
!    This will control some logic...

ivar_wait=0
if (vwait1 > 0.0 .and. vwaittot > 0.0) ivar_wait=1
if(ivar_wait == 0) then
   nwaits=1
elseif (ivar_wait == 1) then
   print*
   print*,' Will wait for varfiles if not present'
   print*,'  Check interval:',vwait1,' Fail time:',vwaittot
   print*
   nwaits=int(vwaittot/vwait1)+1
endif

if (ivflag == 0) then   ! Initialization of initial fields

![ED2
   call date_make_big(iyeara,imontha,idatea,itimea*100  &
                ,itotdate_current)
!ED2]
   wait: do nw = 1, nwaits

      ! Inventory all varf files.

      call varf_file_inv (varfpfx,iyeara,imontha,idatea,itimea)

      ! The initial time must have an exact time match.
      nvarffl=0
      do nf=1,nvarffiles
         if(itotdate_current == itotdate_varf(nf)) then
            nvarffl=nf
            exit wait
         endif
      enddo
      if (ivar_wait == 0 ) then
         print*
         print*,'No initial varfiles found with prefix: ',trim(varfpfx)
         print*
         stop 'no initial varfile'
      elseif(ivar_wait == 1 .and. nw < nwaits) then
         print*,'No initial varfiles found: ',trim(varfpfx)
         print*,'    Waiting:', vwait1,' seconds.   Total wait:',nw*vwait1
         ivwait1=nint(vwait1)
         irslp = irsleep(ivwait1)
      elseif(ivar_wait == 1 .and. nw == nwaits) then
         print*,'No initial varfiles found: ',trim(varfpfx)
         print*,'    Waited too long:', vwaittot,' seconds.'
         stop ' No initial varfile after wait.'
      endif

   enddo wait

   ! Now do actual initialization for the coarse grid
   !     and find 1D reference state
   call newgrid(1)
   call varf_update(0,ifileok,1)

   !  On all fine grids, initialize the 1-D reference state arrays,
   !  the 3-D reference state arrays,
   !  and the prognostic atmospheric fields by interpolation.

   call fmrefs1d(2,ngrids)

   do ifm = 2,ngrids
      icm = nxtnest(ifm)

      ! Get 3D reference state for this grid
      call fmrefs3d(ifm,0)

      ! Interpolate prognostic fields. These will be overwritten if the varfile
      !    exists
      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,1,0)

      call newgrid(ifm)

      ! See if this grid's varfile is created.
      call varf_update(0,ifileok,1)

      if (ifileok  ==  1) then
         ! Everything's cool...
         print*,'Initial varfile read of grid-',ifm
      else
         ! Using interpolated nudging arrays from parent grid.
         print*,'Initial interpolation of grid-',ifm
      endif

      call fmdn0(ifm)

   enddo

   ! ALF
   lastdate_iv = itotdate_varf(nvarffl)

   return

elseif (ivflag == 1) then   ! Fill nudging arrays and compute weights

   ! If a history start, we will assume a past time file is there.
   !  Take closest past time.

   ! This call was deleted :
   call date_add_to(iyeara,imontha,idatea,itimea*100  &
        ,time,'s',iyears,imonths,idates,ihours)
   
   call date_make_big(iyears,imonths,idates, &
        ihours,itotdate_current)
   
   call varf_file_inv (varfpfx,iyeara,imontha,idatea,itimea)
   
   nvarffl=0
   do nf=nvarffiles,1,-1
      if(itotdate_varf(nf) <= itotdate_current ) then
         nvarffl=nf
         exit
      endif
   enddo
   if (nvarffl == 0 ) then
      print*
      print*,'No past varfiles found on nudge fill:',trim(varfpfx)
      print*
      stop 'no past varfile on fill'
   endif

   ! Compute weighting factors for grid 1
   call varweight(nnzp(1),nnxp(1),nnyp(1),varinit_g(1)%varwts(1,1,1)  &
                 ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1))


   ! Read files

   do ifm = 1,ngrids
      icm = nxtnest(ifm)

      call newgrid(ifm)

      ! Interpolate weights to all other grids
      if (ifm > 1) call vfintrpf(ifm,1)

      ! See if this grid's varfile is created.
      call varf_update(0,ifileok,0)

      if (ifileok  ==  1) then
         ! Everything's cool...
         print*,'Varfile read of grid-',ifm
      else
         ! Using interpolated nudging arrays from parent grid.
         call vfintrpf(ifm,2)
         print*,'Interpolation of grid-',ifm
      endif

   enddo

   vtime2=varf_times(nvarffl)
   print*,'New varfile times:',nvarffl,vtime2

   ! ALF
   lastdate_iv = itotdate_varf(nvarffl)

elseif (ivflag == 2) then   ! Runtime file increment

   ! Find current date/time
   call date_add_to(iyeara,imontha,idatea,itimea*100  &
        ,time,'s',iyears,imonths,idates,ihours)
   
   call date_make_big(iyears,imonths,idates,ihours  &
        ,itotdate_current)
   
endif

! Find the next varfile in the list, waiting for it if necessary


wait2: do nw = 1, nwaits

   ! Redo the inventory in case new files showed up
   call varf_file_inv (varfpfx,iyeara,imontha,idatea,itimea)

   nvarffl=0
   do nf=1,nvarffiles
      if(itotdate_varf(nf) > itotdate_current) then
         nvarffl=nf
         exit wait2
      endif
   enddo
   if (ivar_wait == 0 ) then
      print*
      print*,'No future varfiles found with prefix: ',trim(varfpfx)
      print*
      stop 'no future varfile'
   elseif(ivar_wait == 1 .and. nw < nwaits) then
      print*,'No future varfiles found: ',trim(varfpfx)
      print*,'    Waiting:', vwait1,' seconds.   Total wait:',nw*vwait1
      ivwait1=nint(vwait1)
      irslp = irsleep(ivwait1)
   elseif(ivar_wait == 1 .and. nw == nwaits) then
      print*,'No future varfiles found: ',trim(varfpfx)
      print*,'    Waited too long:', vwaittot,' seconds.'
      stop ' No future varfile after wait.'
   endif

enddo wait2

! Read future files

do ifm = 1,ngrids
   icm = nxtnest(ifm)

   call newgrid(ifm)

   ! See if this grid's varfile is created.
   call varf_update(1,ifileok,0)

   if (ifileok  ==  1) then
      ! Everything's cool...
      print*,'Future varfile read of grid-',ifm
   else
      ! Using interpolated nudging arrays from parent grid.
      call vfintrpf(ifm,2)
      print*,'Future interpolation of grid-',ifm
   endif

enddo


vtime1=vtime2
vtime2=varf_times(nvarffl)

print*,'New varfile times:',nvarffl,vtime1,vtime2

! ALF
lastdate_iv = itotdate_varf(nvarffl)
return
end

! **********************************************************************


subroutine varf_file_inv (varpref,iyear1,imonth1,idate1,itime1)

use mem_varinit

implicit none

character(len=*) :: varpref
integer :: iyear1,imonth1,idate1,itime1

integer :: nc,nf,lnf,nvftot
integer :: inyear,inmonth,indate,inhour


character(len=128), dimension(maxnudfiles) :: fnames
character(len=128) :: vpref
character(len=14)  :: itotdate
real(kind=8) :: secs_init,secs_varf

! Get abs seconds of run start

call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through history files and make inventory

nc=len_trim(varpref)
nvftot=-1
vpref=varpref

call RAMS_filelist(fnames,trim(vpref)  &
         //'*.tag',nvftot)

if(nvftot > maxnudfiles) then
   print*,'too many varf files',nvftot,maxnudfiles
   stop 'lots_of_varf'
endif

nvarffiles=0
do nf=1,nvftot
   lnf=len_trim(fnames(nf))
   read(fnames(nf)(lnf-20:lnf-4),20) inyear,inmonth,indate,inhour
   20 format(i4,1x,i2,1x,i2,1x,i6)

   call date_make_big(inyear,inmonth,indate,inhour,itotdate)

   nvarffiles=nvarffiles+1
   fnames_varf(nvarffiles)=fnames(nf)
   itotdate_varf(nvarffiles)=itotdate

   call date_abs_secs2(inyear,inmonth,indate,inhour,secs_varf)
   varf_times(nvarffiles)=secs_varf - secs_init

enddo

call RAMS_dintsort(nvarffiles,itotdate_varf,fnames_varf)

!  start printing section
!--------------------------------------------------------------

print*,' '
print*,' '
print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  Varfile Input Inventory -------------'
print*,'-------------------------------------------------------------'
do nf=1,nvarffiles
   print 8,  nf, itotdate_varf(nf),varf_times(nf) ,trim(fnames_varf(nf))
enddo
8 format(i4,1x,a16,1x,f10.0,2x,a)
print*,'------------------------------------------------------'

!--------------------------------------------------------------

return
end
