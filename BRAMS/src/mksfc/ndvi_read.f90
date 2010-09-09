!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine ndvi_read(runflag,ifm,ierr)

use mem_grid
use mem_leaf
use io_params

implicit none

integer :: runflag,ifm,ierr

character(len=14), save  :: totdate_start,totdate_init
character(len=14)  :: totdate,totdatem
integer :: iyears,imonths,idates,ihours,nf,ng,i,j,ip
real :: timefac_ndvi
real(kind=8) :: secs_init,secs1,secs2,tsf

tsf=365.0

ierr = 0

if (runflag == 1 .or. runflag == 2) then   ! Initialization(1) or file check(2)

   ! Inventory all ndvi surface files. 
   call ndvi_file_inv (ndvifpfx,ierr)

   if(ierr == 1) then
      if(runflag == 2) return
      if(runflag == 1) stop 'ndvi_read: error on init'
   endif

   ! Find init and start date
   call date_make_big(iyeara,imontha,idatea,itimea*100  &
                ,totdate_init)
   if (runtype == 'HISTORY') then

      call date_add_to_big(totdate_init,time,'s',totdate_start)

   elseif (runtype == 'INITIAL') then
      totdate_start=totdate_init
   elseif (runtype == 'MAKEVFILE') then
      totdate_start=totdate_init
   endif

   ! Do some checks on times
   do ng = 1, ngrids
      if(itotdate_ndvi(1,ng) > totdate_start .and. indvicycdata == 0) then
         print*, 'initial ndvi file time for grid ',ng  &
               ,' later than beginning of run'
         ierr = 1 ;  return
      endif
         
      if(iupdndvi == 1 .and. nndvifiles(ng) == 1) then
         print*, 'updating ndvi values but only one ndvi file'  &
                  ,' for grid',ng
         ierr = 1 ;  return
      endif
      
      call date_add_to_big(totdate_init,timmax,'s',totdatem)

      if(iupdndvi == 1 .and. indvicycdata == 0 .and. &
               itotdate_ndvi(nndvifiles(ng),ng) < totdatem) then
         print*, 'final ndvi file time for grid ',ng  &
                  ,'earlier than end of run - making new ndvi files'
         ierr = 1 ;  return
      endif
   enddo
   
   if(ierr == 1) then
      if(runflag == 1) stop 'ndvi_read: error on init time checks'
   endif
   
   ! If we are only checking, we're done.
   if(runflag == 2) return 

   do ng = 1, ngrids
   
      ! Change the ndvi file dates to the current year. We will increment
      !   when we need to read a future file.
      if(indvicyclic == 1) then
         do nf=1,nndvifiles(ng)
            itotdate_ndvi(nf,ng)(1:4)=totdate_start(1:4)
         enddo
      endif
         
      ! Find past time file. The files are ordered in time, so we only 
      !    need to find when start time is greater than a file.
      indviflp(ng)=0
      do nf=nndvifiles(ng),1,-1
         if(totdate_start >= itotdate_ndvi(nf,ng) ) then
            indviflp(ng)=nf
            exit
         endif
      enddo
      
      indviflf(ng)=indviflp(ng)+1
   
      ! If we are cyclic, we possibly didn't find a start time later than
      !   a file time. I think it is safe to assume that we will use the 
      !   last file for the start file. Also see if future time will be on 
      !   next cycle.
      if(indvicycdata == 1 ) then
         if(indviflp(ng) == 0) indviflp(ng) = nndvifiles(ng)
         if(indviflf(ng) > nndvifiles(ng)) indviflf(ng)=1
      endif
   
      ! Read past/future time ndvi field

      call ndvi_update(0,indviflp(ng))
      
      if(iupdndvi == 1) then
         ! Read future time ndvi field if updating
         call ndvi_update(1,indviflf(ng))

         ! Compute times as number of seconds past 1 Jan 1900
         call date_abs_secs(totdate_init,secs_init)
         
         ! Get model time of past file
         totdatem=itotdate_ndvi(indviflp(ng),ng)
         if(indvicyclic == 1) then
            ! If month of past file > current month, subtract a year
            if(totdatem(5:6) > totdate_start(5:6))   &
               call date_add_to_big(totdatem,-tsf,'d',totdatem)
         endif
         call date_abs_secs(totdatem,secs1)
   
         totdatem=itotdate_ndvi(indviflf(ng),ng)
         if(indvicyclic == 1) then
            if(totdatem < totdate_start) then
               ! Future file is in next year. Update all file names for a new year
               call date_add_to_big(totdatem, tsf,'d',totdatem)
               do nf=1,nndvifiles(ng)
                  itotdate_ndvi(nf,ng)(1:4)=totdatem(1:4)
               enddo
            endif
         endif
         call date_abs_secs(totdatem,secs2)
      
         ndvitime1(ng) = secs1 - secs_init
         ndvitime2(ng) = secs2 - secs_init
      else
         do ip=1,npatch
            leaf_g(ng)%veg_ndvip(1:nnxp(ng),1:nnyp(ng),ip)=  &
            leaf_g(ng)%veg_ndvif(1:nnxp(ng),1:nnyp(ng),ip)
         enddo
         ndvitime1(ng) = 0.
         ndvitime2(ng) = 0.
      endif

      ! Fill "current" time ndvi values
      if (iupdndvi == 0) then
         timefac_ndvi = 0.
      else
         timefac_ndvi = (time - ndvitime1(ng))   &
                      / (ndvitime2(ng) - ndvitime1(ng))
      endif

      do ip=1,npatch
         do j=1,nnyp(ng)
            do i=1,nnxp(ng)
               leaf_g(ng)%veg_ndvic(i,j,ip) = leaf_g(ng)%veg_ndvip(i,j,ip)  &
                  + (leaf_g(ng)%veg_ndvif(i,j,ip)  &
                  -  leaf_g(ng)%veg_ndvip(i,j,ip)) * timefac_ndvi
            enddo
         enddo
      enddo     
      
   enddo
   
   return

elseif (runflag == 3) then   ! Runtime file increment
   
   call date_add_to_big(totdate_init,time,'s',totdate)
   
   if ( time >= ndvitime2(ifm) ) then
   
      ! Update ndvi fields
      indviflp(ifm) = indviflf(ifm)
      indviflf(ifm) = indviflp(ifm) + 1
      
      ! Compute times as number of seconds past 1 Jan 1900
      !   If cyclic, modify file date 
      call date_abs_secs(totdate_init,secs_init)
      call date_add_to_big(totdate_init,time,'s',totdate)

      totdatem=itotdate_ndvi(indviflp(ifm),ifm)
      call date_abs_secs(totdatem,secs1)
      
      ! Need to deal with cyclic. 
      if(indvicyclic == 1 .and. indviflf(ifm) > nndvifiles(ifm)) then
         indviflf(ifm)=1
         ! Update all file names for a new year
         totdatem=itotdate_ndvi(indviflf(ifm),ifm)
         call date_add_to_big(totdatem, tsf,'d',totdatem)
         do nf=1,nndvifiles(ifm)
            itotdate_ndvi(nf,ifm)(1:4)=totdatem(1:4)
         enddo
      endif
         
      totdatem=itotdate_ndvi(indviflf(ifm),ifm)
      call date_abs_secs(totdatem,secs2)
      
      ndvitime1(ifm) = secs1 - secs_init
      ndvitime2(ifm) = secs2 - secs_init
   
      ! Finally read the actual field           
      call ndvi_update(1,indviflf(ifm))

      print*,'Switched ndvi files:',ndvitime1(ifm),ndvitime2(ifm)  &
               ,secs1,secs2,secs_init

   else
   
      return
   
   endif
   
endif

return
end



subroutine ndvi_file_inv (ndvifilin,ierr)

use mem_grid
use io_params

implicit none

character(len=*) :: ndvifilin
integer :: ierr

integer :: nc,nf,lnf,nftot,ng
integer :: inyear,inmonth,indate,inhour


character(len=128), dimension(maxndvifiles) :: fnames
character(len=14)  :: itotdate
character(len=1)  :: cgrid
real(kind=8) :: secs_init,secs_file

ierr=0

! Get abs seconds of run start

call date_abs_secs2(iyeara,imontha,idatea,itimea*100,secs_init)

! Go through ndvi files and make inventory. We unfortunately have to do this
!   for all grids.

indvicyclic=0
indvicycdata=0

do ng = 1, ngrids
   
   nftot = -1
   write(cgrid,'(i1)') ng
   call RAMS_filelist(fnames,trim(ndvifilin)//'-N-*-g'//cgrid//'.vfm',nftot)
   
   if(nftot <= 0) then
      print*,'No ndvi files for grid '//cgrid
      ierr=1
      return
   endif   
   
   if(nftot > maxndvifiles) then
      print*,'too many ndvi files'
      stop 'ndvi_file_inv: lots_of_ndvi'
   endif
   

   ! We need to see if the data is to be considered "cyclic". Count how many files
   !  have the year 0000. If all of them do not, we probably have an error.
   !  indvicycdata = 1 when data is cyclic
   !  indvicyclic =1 when we are have cyclic data and will be updating in time

   nndvifiles(ng)=0
   do nf=1,nftot
      lnf=len_trim(fnames(nf))
      read(fnames(nf)(lnf-23:lnf-7),20) inyear,inmonth,indate,inhour
      20 format(i4,1x,i2,1x,i2,1x,i6)

      ! Check the file headers
      call ndvi_check_header(ng,fnames(nf),ierr)
      if(ierr == 1) return
   
      if(inyear == 0) indvicycdata=indvicycdata+1

      call date_make_big(inyear,inmonth,indate,inhour,itotdate)

      nndvifiles(ng)=nndvifiles(ng)+1
      fnames_ndvi(nndvifiles(ng),ng)=fnames(nf)
      itotdate_ndvi(nndvifiles(ng),ng)=itotdate

   enddo

   call RAMS_dintsort(nndvifiles(ng),itotdate_ndvi(1,ng),fnames_ndvi(1,ng))


   !  start printing section
   !--------------------------------------------------------------

   print*,' '
   print*,' '
   print*,' '
   print*,'-------------------------------------------------------------'
   print*,'-----------  NDVI Input File Inventory: Grid '//cgrid
   print*,'-------------------------------------------------------------'
   do nf=1,nndvifiles(ng)
      print*,  itotdate_ndvi(nf,ng),'   ',trim(fnames_ndvi(nf,ng))
   enddo
   print*,'------------------------------------------------------'

enddo

! Check the cyclic condition. Only relevant if we are updating in time.
!   WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS
if(indvicycdata > 0) then
   if(indvicycdata /= sum(nndvifiles(1:ngrids))) then
      print*, 'All ndvi surface files do not have year 0000'
      print*, 'This confuses the gods and can not occur.'
      stop 'ndvi_inv'
   endif
   indvicycdata=1
else
   indvicycdata=0
endif


! Set the main cyclic flag. Only relevant if we are updating in time.
if(iupdndvi == 1 .and. indvicycdata == 1) indvicyclic=1

return
end


!--------------------------------------------------------------

subroutine ndvi_update(iswap,nfile)

use mem_leaf
use mem_grid
use io_params

implicit none

integer :: iswap,nfile

integer,save :: iun=25

integer :: ng,nc,ip
character(len=1) :: cgrid
character(len=128) :: flnm
character(len=1) :: dummy


! Put new fields into future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then
   do ng=1,ngrids
      do ip = 1,npatch
         leaf_g(ng)%veg_ndvip(1:nnxp(ng),1:nnyp(ng),ip)=  &
            leaf_g(ng)%veg_ndvif(1:nnxp(ng),1:nnyp(ng),ip)
      enddo
   enddo
endif


! Open the input file for each grid and read field.
do ng=1,ngrids
   ! Contruct file name for this grid
   write(cgrid, '(i1)') ng
   flnm=fnames_ndvi(nfile,ng)
   nc=len_trim(flnm)-4
   flnm(nc:nc)=cgrid

   call rams_f_open(iun,flnm,'FORMATTED','OLD','READ',0)
   ! Skip header lines
   read(iun,*) dummy
   read(iun,*) dummy
   read(iun,*) dummy
   read(iun,*) dummy
   do ip = 1,npatch
      call vfirec(iun,leaf_g(ng)%veg_ndvif(:,:,ip),nnxp(ng)*nnyp(ng),'LIN')
   enddo

   ! Close the input file
   close(iun)
enddo


return
end

!****************************************************************************

subroutine ndvi_check_header(ifm,flnm,ierr)

use mem_grid

! This subroutine checks for the existence of an ndvi file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ierr is returned with a value of 1.  If the file 
! exists and is ok, ierr is returned with a value of 0.

implicit none
integer :: ifm,ierr
character*(*) flnm

integer :: iiyear,iimonth,iidate,iihour  &
   ,isfc_marker,isfc_ver,nsfx,nsfy,nsfpat
real :: sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr

ierr = 0

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)

read (25,*) isfc_marker,isfc_ver
read (25,*) iiyear,iimonth,iidate,iihour
read (25,*) nsfx,nsfy,nsfpat
read (25,*) sfdx,sfdy,sfplat,sfplon,sflat,sflon
close (25)

if (nsfx                   .ne. nnxp(ifm) .or.  &
    nsfy                   .ne. nnyp(ifm) .or.  &
    nsfpat                 .ne. npatch    .or.  &
    abs(sfdx-deltaxn(ifm)) .gt. .001      .or.  &
    abs(sfdy-deltayn(ifm)) .gt. .001      .or.  &
    abs(sfplat-platn(ifm)) .gt. .001      .or.  &
    abs(sfplon-plonn(ifm)) .gt. .001      .or.  &
    abs(sflat-glatr)       .gt. .001      .or.  &
    abs(sflon-glonr)       .gt. .001) then

   ierr = 1

   print*,'ndvifile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'npatch:',npatch,nsfpat
   print*,'deltax:',deltaxn(ifm),sfdx
   print*,'deltay:',deltayn(ifm),sfdy
   print*,'platn:',platn(ifm),sfplat
   print*,'plonn:',plonn(ifm),sfplon
   print*,'SW lat:',glatr,sflat
   print*,'SW lon:',glonr,sflon
   print*,'-------------------'

else

   ierr = 0

endif

return
end

