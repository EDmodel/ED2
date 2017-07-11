!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine sst_read(runflag,ifm,ierr)

use mem_grid
use mem_leaf
use io_params

implicit none

integer :: runflag,ifm,ierr

character(len=14), save  :: totdate_start,totdate_init
character(len=14)  :: totdate,totdatem
integer :: iyears,imonths,idates,ihours,nf,ng
real(kind=8) :: secs_init,secs1,secs2

ierr = 0

if (runflag == 1 .or. runflag == 2) then   ! Initialization(1) or file check(2)

   ! Inventory all sst surface files.
   call sst_file_inv (sstfpfx,ierr)
   if(ierr == 1) then
      if(runflag == 2) return
      if(runflag == 1) stop 'sst_read: error on init'
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
      if(itotdate_sst(1,ng) > totdate_start .and. isstcycdata == 0) then
         print*, 'initial sst file time for grid ',ng  &
               ,' later than beginning of run'
         ierr = 1 ;  return
      endif

      if(iupdsst == 1 .and. nsstfiles(ng) == 1) then
         print*, 'updating SST values but only one SST file'  &
                  ,' for grid',ng
         ierr = 1 ;  return
      endif

      call date_add_to_big(totdate_init,timmax,'s',totdatem)

      if(iupdsst == 1 .and. isstcycdata == 0 .and. &
               itotdate_sst(nsstfiles(ng),ng) < totdatem) then
         print*, 'final sst file time for grid ',ng  &
                  ,'earlier than end of run - making new sst files'
         ierr = 1 ;  return
      endif
   enddo

   if(ierr == 1) then
      if(runflag == 1) stop 'sst_read: error on init time checks'
   endif

   ! If we are only checking, we're done.
   if(runflag == 2) return

   do ng = 1, ngrids

      ! Change the sst file dates to the current year. We will increment
      !   when we need to read a future file.
      if(isstcycdata == 1) then
         do nf=1,nsstfiles(ng)
            itotdate_sst(nf,ng)(1:4)=totdate_start(1:4)
         enddo
      endif

      ! Find past time file. The files are ordered in time, so we only
      !    need to find when start time is greater than a file.
      isstflp(ng)=0
      do nf=nsstfiles(ng),1,-1
         if(totdate_start >= itotdate_sst(nf,ng) ) then
            isstflp(ng)=nf
            exit
         endif
      enddo

      isstflf(ng)=isstflp(ng)+1

      ! If we are cyclic, we possibly didn't find a start time later than
      !   a file time. I think it is safe to assume that we will use the
      !   last file for the start file. Also see if future time will be on
      !   next cycle.
      if(isstcycdata == 1 ) then
         if(isstflp(ng) == 0) isstflp(ng) = nsstfiles(ng)
         if(isstflf(ng) > nsstfiles(ng)) isstflf(ng)=1
      endif

      ! Read past time sst field

      call sst_update(0,isstflp(ng))

      if(iupdsst == 1 ) then
         ! Read future time sst field if updating
         call sst_update(1,isstflf(ng))

         ! Compute times as number of seconds past 1 Jan 1900
         call date_abs_secs(totdate_init,secs_init)

         ! Get model time of past file
         totdatem=itotdate_sst(isstflp(ng),ng)
         if(isstcyclic == 1) then
            ! If month of past file > current month, subtract a year
            if(totdatem(5:6) > totdate_start(5:6))   &
               call date_add_to_big(totdatem,dble(-365.),'d',totdatem)
         endif
         call date_abs_secs(totdatem,secs1)

         totdatem=itotdate_sst(isstflf(ng),ng)
         if(isstcyclic == 1) then
            if(totdatem < totdate_start) then
               ! Future file is in next year. Update all file names for a new year
               call date_add_to_big(totdatem,dble(365.),'d',totdatem)
               do nf=1,nsstfiles(ng)
                  itotdate_sst(nf,ng)(1:4)=totdatem(1:4)
               enddo
            endif
         endif
         call date_abs_secs(totdatem,secs2)

         ssttime1(ng) = secs1 - secs_init
         ssttime2(ng) = secs2 - secs_init
         print*,'=======3',ssttime1(ng),ssttime2(ng),secs1,secs2,secs_init
      else
         leaf_g(ng)%seatp(1:nnxp(ng),1:nnyp(ng))=  &
         leaf_g(ng)%seatf(1:nnxp(ng),1:nnyp(ng))
         ssttime1(ng) = 0.
         ssttime2(ng) = 0.
      endif


   enddo

   return

elseif (runflag == 3) then   ! Runtime file increment

   if ( time >= ssttime2(ifm) ) then

      ! Update sst fields
      isstflp(ifm) = isstflf(ifm)
      isstflf(ifm) = isstflp(ifm) + 1

      ! Compute times as number of seconds past 1 Jan 1900
      !   If cyclic, modify file date
      call date_abs_secs(totdate_init,secs_init)

      call date_add_to_big(totdate_init,time,'s',totdate)

      totdatem=itotdate_sst(isstflp(ifm),ifm)
      call date_abs_secs(totdatem,secs1)

      ! Need to deal with cyclic.
      if(isstcyclic == 1 .and. isstflf(ifm) > nsstfiles(ifm)) then
         isstflf(ifm)=1
         ! Update all file names for a new year
         totdatem=itotdate_sst(isstflf(ifm),ifm)
         call date_add_to_big(totdatem,dble(365.),'d',totdatem)
         do nf=1,nsstfiles(ifm)
            itotdate_sst(nf,ifm)(1:4)=totdatem(1:4)
         enddo
      endif

      totdatem=itotdate_sst(isstflf(ifm),ifm)
      call date_abs_secs(totdatem,secs2)

      ssttime1(ifm) = secs1 - secs_init
      ssttime2(ifm) = secs2 - secs_init

      ! Finally read the actual field
      call sst_update(1,isstflf(ifm))

      print*,'Switched sst files:',ssttime1(ifm),ssttime2(ifm)  &
               ,secs1,secs2,secs_init


   else

      return

   endif

endif

return
end subroutine sst_read

!**********************************************************

subroutine sst_read_new(runflag,ifm,ierr)

use mem_grid
use mem_leaf
use io_params

implicit none

integer :: runflag,ifm,ierr

character(len=14), save  :: totdate_start,totdate_init
character(len=14)  :: totdate,totdatem
integer :: iyears,imonths,idates,ihours,nf,ng
real(kind=8) :: secs_init,secs1,secs2



! ALF
character(len=14)  :: itotdate_current
real(kind=8) :: lastdate_sec, initialdate_sec

ierr = 0

if (runflag == 1 .or. runflag == 2) then   ! Initialization(1) or file check(2)

   ! Inventory all sst surface files.
   call sst_file_inv (sstfpfx,ierr)
   if(ierr == 1) then
      if(runflag == 2) return
      if(runflag == 1) stop 'sst_read: error on init'
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
      if(itotdate_sst(1,ng) > totdate_start .and. isstcycdata == 0) then
         print*, 'initial sst file time for grid ',ng  &
               ,' later than beginning of run'
         ierr = 1 ;  return
      endif

      if(iupdsst == 1 .and. nsstfiles(ng) == 1) then
         print*, 'updating SST values but only one SST file'  &
                  ,' for grid',ng
         ierr = 1 ;  return
      endif

      call date_add_to_big(totdate_init,timmax,'s',totdatem)
      if(iupdsst == 1 .and. isstcycdata == 0 .and. &
               itotdate_sst(nsstfiles(ng),ng) < totdatem) then
         print*, 'final sst file time for grid ',ng  &
                  ,'earlier than end of run - making new sst files'
         ierr = 1 ;  return
      endif
   enddo

   if(ierr == 1) then
      if(runflag == 1) stop 'sst_read: error on init time checks'
   endif

   ! If we are only checking, we're done.
   if(runflag == 2) return

   do ng = 1, ngrids

      ! Change the sst file dates to the current year. We will increment
      !   when we need to read a future file.
      if(isstcycdata == 1) then
         do nf=1,nsstfiles(ng)
            itotdate_sst(nf,ng)(1:4)=totdate_start(1:4)
         enddo
      endif

      ! Find past time file. The files are ordered in time, so we only
      !    need to find when start time is greater than a file.
      isstflp(ng)=0
      do nf=nsstfiles(ng),1,-1
         if(totdate_start >= itotdate_sst(nf,ng) ) then
            isstflp(ng)=nf
            exit
         endif
      enddo

      isstflf(ng)=isstflp(ng)+1

      ! If we are cyclic, we possibly didn't find a start time later than
      !   a file time. I think it is safe to assume that we will use the
      !   last file for the start file. Also see if future time will be on
      !   next cycle.
      if(isstcycdata == 1 ) then
         if(isstflp(ng) == 0) isstflp(ng) = nsstfiles(ng)
         if(isstflf(ng) > nsstfiles(ng)) isstflf(ng)=1
      endif

      ! Read past time sst field

      call sst_update(0,isstflp(ng))

      if(iupdsst == 1 ) then
         ! Read future time sst field if updating
         call sst_update(1,isstflf(ng))

         ! Compute times as number of seconds past 1 Jan 1900
         call date_abs_secs(totdate_init,secs_init)

         ! Get model time of past file
         totdatem=itotdate_sst(isstflp(ng),ng)
         if(isstcyclic == 1) then
            ! If month of past file > current month, subtract a year
            if(totdatem(5:6) > totdate_start(5:6))   &
               call date_add_to_big(totdatem,dble(-365.),'d',totdatem)
         endif
         call date_abs_secs(totdatem,secs1)

         totdatem=itotdate_sst(isstflf(ng),ng)
         if(isstcyclic == 1) then
            if(totdatem < totdate_start) then
               ! Future file is in next year. Update all file names for a new year
               call date_add_to_big(totdatem,dble(365.),'d',totdatem)
               do nf=1,nsstfiles(ng)
                  itotdate_sst(nf,ng)(1:4)=totdatem(1:4)
               enddo
            endif
         endif
         call date_abs_secs(totdatem,secs2)

         ssttime1(ng) = secs1 - secs_init
         ssttime2(ng) = secs2 - secs_init
         print*,'=======3',ssttime1(ng),ssttime2(ng),secs1,secs2,secs_init
      else
         leaf_g(ng)%seatp(1:nnxp(ng),1:nnyp(ng))=  &
         leaf_g(ng)%seatf(1:nnxp(ng),1:nnyp(ng))
         ssttime1(ng) = 0.
         ssttime2(ng) = 0.
      endif

      ! ALF
      lastdate_sst(ng) = itotdate_sst(isstflf(ng),ng)

   enddo

   return

elseif (runflag == 3) then   ! Runtime file increment

   if (time>=ssttime1(ifm) .and. (time+.0001)<timmax) then

      !Checking if all values have been read
      call date_make_big(iyeara,imontha,idatea,itimea,itotdate_current)
      call date_abs_secs(itotdate_current, initialdate_sec)
      call date_abs_secs(lastdate_sst(ifm), lastdate_sec)
      if (((lastdate_sec-initialdate_sec)+.0001)>=timmax) return

      ! Update sst fields
      isstflp(ifm) = isstflf(ifm)
      isstflf(ifm) = isstflp(ifm) + 1

      ! Compute times as number of seconds past 1 Jan 1900
      !   If cyclic, modify file date
      call date_abs_secs(totdate_init,secs_init)
      call date_add_to_big(totdate_init,time,'s',totdate)

      totdatem=itotdate_sst(isstflp(ifm),ifm)
      call date_abs_secs(totdatem,secs1)

      ! Need to deal with cyclic.
      if(isstcyclic == 1 .and. isstflf(ifm) > nsstfiles(ifm)) then
         isstflf(ifm)=1
         ! Update all file names for a new year
         totdatem=itotdate_sst(isstflf(ifm),ifm)
         call date_add_to_big(totdatem,dble(365.),'d',totdatem)
         do nf=1,nsstfiles(ifm)
            itotdate_sst(nf,ifm)(1:4)=totdatem(1:4)
         enddo
      endif

      totdatem=itotdate_sst(isstflf(ifm),ifm)
      call date_abs_secs(totdatem,secs2)

      ssttime1(ifm) = secs1 - secs_init
      ssttime2(ifm) = secs2 - secs_init

      ! Finally read the actual field
      call sst_update(1,isstflf(ifm))

      print*,'Switched sst files:',ssttime1(ifm),ssttime2(ifm)  &
               ,secs1,secs2,secs_init

      ! ALF
      lastdate_sst(ifm) = itotdate_sst(isstflf(ifm),ifm)

   else

      return

   endif

endif

return
end subroutine sst_read_new



subroutine sst_file_inv (sstfilin,ierr)

use grid_dims, only : maxfiles, str_len
use mem_grid
use io_params

implicit none

character(len=str_len) :: sstfilin
integer :: ierr

integer :: nc,nf,lnf,nftot,ng
integer :: inyear,inmonth,indate,inhour

character(len=str_len), dimension(maxfiles) :: fnames
character(len=str_len)                      :: thisfile
character(len=14)  :: itotdate
character(len=1)  :: cgrid
real(kind=8) :: secs_init,secs_file

ierr=0

! Get abs seconds of run start

 call date_abs_secs2(iyeara,imontha,idatea,itimea*100,secs_init)


! Go through sst files and make inventory. We unfortunately have to do this
!   for all grids.

isstcyclic=0
isstcycdata=0

do ng = 1, ngrids

   nftot = -1
   write(cgrid,'(i1)') ng
   thisfile = trim(sstfilin)//'-W-*-g'//cgrid//'.vfm'
   call RAMS_filelist(maxfiles,fnames,thisfile,nftot)

   if(nftot <= 0) then
      print*,'No sst files for grid '//cgrid
      ierr=1
      return
   endif

   if(nftot > maxsstfiles) then
      print*,'too many sst files'
      stop 'sst_file_inv: lots_of_sst'
   endif

   ! We need to see if the data is to be considered "cyclic". Count how many files
   !  have the year 0000. If all of them do not, we probably have an error.
   !  isstcycdata = 1 when data is cyclic
   !  isstcyclic =1 when we are have cyclic data and will be updating in time

   nsstfiles(ng)=0
   do nf=1,nftot
      lnf=len_trim(fnames(nf))
      read(fnames(nf)(lnf-23:lnf-7),fmt='(i4,1x,i2,1x,i2,1x,i6)') &
          inyear,inmonth,indate,inhour


      ! Check the file headers
      call sst_check_header(ng,fnames(nf),ierr)
      if(ierr == 1) return

      if(inyear == 0) isstcycdata=isstcycdata+1

      call date_make_big(inyear,inmonth,indate,inhour,itotdate)

      nsstfiles(ng)=nsstfiles(ng)+1
      fnames_sst(nsstfiles(ng),ng)=fnames(nf)
      itotdate_sst(nsstfiles(ng),ng)=itotdate

   enddo

   call RAMS_dintsort(nsstfiles(ng),itotdate_sst(1,ng),fnames_sst(1,ng))


   !  start printing section
   !--------------------------------------------------------------

   write(unit=*,fmt='(a)'      ) ' '
   write(unit=*,fmt='(a)'      ) ' '
   write(unit=*,fmt='(a)'      ) ' '
   write(unit=*,fmt='(a)'      )  '-------------------------------------------------------'
   write(unit=*,fmt='(a,1x,i4)')  '  SST Input File Inventory: Grid ',ng
   write(unit=*,fmt='(a)'      )  '-------------------------------------------------------'
   do nf=1,nsstfiles(ng)
      write(unit=*,fmt='(i6,2(1x,a))') nf,itotdate_sst(nf,ng),trim(fnames_sst(nf,ng))
   end do
   write(unit=*,fmt='(a)'      )  '-------------------------------------------------------'
end do

! Check the cyclic data condition.
!   WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS
if(isstcycdata > 0) then
   if(isstcycdata /= sum(nsstfiles(1:ngrids))) then
      print*, 'All sst surface files do not have year 0000'
      print*, 'This confuses the gods and can not occur.'
      stop 'sst_inv'
   endif
   isstcycdata=1
else
   isstcycdata=0
endif


! Set the main cyclic flag. Only relevant if we are updating in time.
if(iupdsst == 1 .and. isstcycdata == 1) isstcyclic=1

return
end


!--------------------------------------------------------------

subroutine sst_update(iswap,nfile)

use mem_grid
use mem_leaf
use io_params

implicit none

integer :: iswap,nfile

integer,save :: iun=25

integer :: ng,nc
character(len=1) :: cgrid
character(len=str_len) :: flnm
character(len=1) :: dummy


! Put new fields into future arrays. If iswap == 1,
!     swap future into past first

if (iswap == 1) then
   do ng=1,ngrids
      leaf_g(ng)%seatp(1:nnxp(ng),1:nnyp(ng))=  &
         leaf_g(ng)%seatf(1:nnxp(ng),1:nnyp(ng))
   enddo
endif


! Open the input file for each grid and read field.
do ng=1,ngrids
   ! Contruct file name for this grid
   write(cgrid, '(i1)') ng
   flnm=fnames_sst(nfile,ng)
   nc=len_trim(flnm)-4
   flnm(nc:nc)=cgrid

   call rams_f_open(iun,flnm,'FORMATTED','OLD','READ',0)
   ! Skip header lines
   read(25,*) dummy
   read(25,*) dummy
   read(25,*) dummy
   read(25,*) dummy
   call vfirec(25,leaf_g(ng)%seatf,nnxp(ng)*nnyp(ng),'LIN')
   !call ezcntr(leaf_g(ng)%seatf(1,1),nnxp(ng),nnyp(ng))

   ! Close the input file
   close(iun)
enddo


return
end subroutine sst_update

!****************************************************************************

subroutine sst_check_header(ifm,flnm,ierr)

use mem_grid

! This subroutine checks for the existence of an sst file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ierr is returned with a value of 1.  If the file
! exists and is ok, ierr is returned with a value of 0.

implicit none
integer :: ifm,ierr
character*(*) flnm

integer :: iiyear,iimonth,iidate,iihour  &
   ,isfc_marker,isfc_ver,nsfx,nsfy
real :: sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr

ierr = 0

!!$print*,'------------------------------------------------'
!!$print*,'---> Check grid:',ifm,' sst file. '
!!$print*,'--->   Filename:',trim(flnm)


call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)

read (25,*) isfc_marker,isfc_ver
read (25,*) iiyear,iimonth,iidate,iihour
read (25,*) nsfx,nsfy
read (25,*) sfdx,sfdy,sfplat,sfplon,sflat,sflon
close (25)

if (nsfx                   .ne. nnxp(ifm) .or.  &
    nsfy                   .ne. nnyp(ifm) .or.  &
    abs(sfdx-deltaxn(ifm)) .gt. .001      .or.  &
    abs(sfdy-deltayn(ifm)) .gt. .001      .or.  &
    abs(sfplat-platn(ifm)) .gt. .001      .or.  &
    abs(sfplon-plonn(ifm)) .gt. .001      .or.  &
    abs(sflat-glatr)       .gt. .001      .or.  &
    abs(sflon-glonr)       .gt. .001) then

   ierr = 1

   print*,'SSTfile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'deltax:',deltaxn(ifm),sfdx
   print*,'deltay:',deltayn(ifm),sfdy
   print*,'platn:',platn(ifm),sfplat
   print*,'plonn:',plonn(ifm),sfplon
   print*,'SW lat:',glatr,sflat
   print*,'SW lon:',glonr,sflon
   print*,'-------------------'

else

   ierr = 0
!!$   print*,'---> Grid:',ifm,' sstfile data okay. '
!!$   print*,'------------------------------------------------'

endif

return
end

