!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine cu_read(initflag)

!------------------------------------------------------
!  Read cumulus inversion tendencies
!------------------------------------------------------

use mem_cuparm
use mem_grid

implicit none

integer :: initflag

character(len=14)  :: itotdate_start
integer :: iyears,imonths,idates,ihours,nf,ifm

if (initflag == 1) then   ! Initialization

   ! Inventory all cu inversion files. 
   call cu_file_inv (iyeara,imontha,idatea,itimea)

   ! Find past time file

   if (runtype == 'HISTORY') then
      call date_add_to(iyeara,imontha,idatea,itimea*100  &
                   ,max(time,tcu_beg),'s',iyears,imonths,idates,ihours)
      call date_make_big(iyears,imonths,idates,ihours  &
                   ,itotdate_start)
   
   elseif (runtype == 'INITIAL') then
      call date_add_to(iyeara,imontha,idatea,itimea*100  &
                   ,tcu_beg,'s',iyears,imonths,idates,ihours)
      call date_make_big(iyears,imonths,idates,ihours  &
                   ,itotdate_start)
   endif

   do nf=1,ncufiles
      if(itotdate_start >= itotdate_cu(nf) .and.  &
            itotdate_start <  itotdate_cu(nf+1) ) then
         ncufl=nf
         exit
      endif
   enddo

   print*,'nud starting at history file:',ncufl

   ! Read initial files.
   
   call cu_update(0,ncufl)
   

elseif (initflag == 2) then   ! Runtime file increment

   if ( time >= cu_times(ncufl+1) ) ncufl = ncufl + 1

endif


! Read new files.

call cu_update(1,ncufl+1)


cutime1=cu_times(ncufl)
cutime2=cu_times(ncufl+1)

return
end


subroutine cu_file_inv (iyear1,imonth1,idate1,itime1)

use mem_cuparm
use grid_dims, only : str_len

implicit none

integer :: iyear1,imonth1,idate1,itime1

integer :: nc,nf,lnf,nhftot
integer :: inyear,inmonth,indate,inhour


character(len=str_len), dimension(maxcufiles) :: fnames
character(len=str_len) :: thisfile
character(len=14)  :: itotdate
real(kind=8) :: secs_init,secs_cu

! Get abs seconds of run start

call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through and make inventory

nhftot=-1
thisfile = cu_prefix(1:len_trim(cu_prefix))//'*.vfm'
call RAMS_filelist(maxcufiles,fnames,thisfile,nhftot)

if(nhftot > maxcufiles) then
   print*,'too many cu files'
   stop 'lots_of_cu_files'
endif

ncufiles=0
do nf=1,nhftot

   ! only save grid 1 files names and times. 

   if (index(fnames(nf),'-g1.') /= 0) then
      lnf=len_trim(fnames(nf))
      read(fnames(nf)(lnf-23:lnf-7),20) inyear,inmonth,indate,inhour
      20 format(i4,1x,i2,1x,i2,1x,i6)

      call date_make_big(inyear,inmonth,indate,inhour,itotdate)

      ncufiles=ncufiles+1
      fnames_cu(ncufiles)=fnames(nf)
      itotdate_cu(ncufiles)=itotdate
   
      call date_abs_secs2(inyear,inmonth,indate,inhour,secs_cu)
      cu_times(ncufiles)=secs_cu - secs_init
   endif

enddo

call RAMS_dintsort(ncufiles,itotdate_cu,fnames_cu)

!  start printing section
!--------------------------------------------------------------

print*,' '
print*,' '
print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  Cumulus Tendency Input File Inventory -------------'
print*,'-------------------------------------------------------------'
do nf=1,ncufiles
   print*,  itotdate_cu(nf),'   ',cu_times(nf)  &
           ,fnames_cu(nf)(1:len_trim(fnames_cu(nf)))
enddo
print*,'------------------------------------------------------'

!--------------------------------------------------------------


return
end

!******************************************************************************

subroutine cu_update(iswap,ncu)

use mem_cuparm
use mem_basic
use mem_grid

implicit none

integer :: iswap,ncu

character (len=str_len) :: cunamein
character (len=1) :: cng
integer :: ngr,nc,npts,ifm,icm
integer,save :: iun=10
logical :: there

! Put new fields into cu future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then
   do ngr=1,ngrids
      cuparm_g(ngr)%thsrcp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
         cuparm_g(ngr)%thsrcf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
      cuparm_g(ngr)%rtsrcp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
         cuparm_g(ngr)%rtsrcf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
      cuparm_g(ngr)%conprrp(1:nnxp(ngr),1:nnyp(ngr))=  &
         cuparm_g(ngr)%conprrf(1:nnxp(ngr),1:nnyp(ngr))
   enddo
endif


! Open the input file and read fields


do ngr=1,ngrids

   ifm=ngr
   icm=nxtnest(ifm)

   print*,'ncu:',ncu,fnames_cu(ncu)
   nc=len_trim(fnames_cu(ncu))
   write(cng,'(i1)') ngr
   cunamein=fnames_cu(ncu)(1:nc-5)//cng//'.vfm'
   print*,'ncu:',ncu,cunamein

   inquire (file=cunamein, exist=there)
   
   if (wt_cu_grid(ngr) > 0.) then
      if (there) then

         call rams_f_open(iun,cunamein,'FORMATTED','OLD','READ',0)

         npts=nnzp(ngr)*nnxp(ngr)*nnyp(ngr)
         call vfirec(iun,cuparm_g(ngr)%thsrcf,npts,'LIN')
         call vfirec(iun,cuparm_g(ngr)%rtsrcf,npts,'LIN')

         npts=nnxp(ngr)*nnyp(ngr)
         call vfirec(iun,cuparm_g(ngr)%conprrf,npts,'LIN')
      else
         call fmint4(cuparm_g(icm)%thsrcf  &
                    ,cuparm_g(ifm)%thsrcf  &
                    ,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
                    ,grid_g(icm)%topt,ifm,icm,'t',1)
         call fmint4(cuparm_g(icm)%rtsrcf  &
                    ,cuparm_g(ifm)%thsrcf  &
                    ,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
                    ,grid_g(icm)%topt,ifm,icm,'t',1)
         call fmint2d(icm,ifm,'t'  &
                     ,cuparm_g(icm)%conprrf  &
                     ,cuparm_g(ifm)%conprrf )
      endif
   endif

enddo

! Close the input file

close(iun)

return
end


