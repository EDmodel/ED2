!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine cond_read(initflag)

use mem_grid
use mem_varinit

implicit none

integer initflag


character(len=14)  :: itotdate_start
integer :: iyears,imonths,idates,ihours,nf,ifm


if (initflag == 1) then   ! Initialization

   ! Inventory all history files. 

   call cond_file_inv (cond_hfile,iyeara,imontha,idatea,itimea)


   ! Find past time nudging file

   if (runtype == 'HISTORY') then
      call date_add_to(iyeara,imontha,idatea,itimea*100  &
                   ,time,'s',iyears,imonths,idates,ihours)
      call date_make_big(iyears,imonths,idates,ihours  &
                   ,itotdate_start)
   elseif (runtype == 'INITIAL') then
      call date_make_big(iyeara,imontha,idatea,itimea*100  &
                   ,itotdate_start)
   endif

   do nf=1,ncondfiles
      if(itotdate_start >= itotdate_cond(nf) .and.  &
            itotdate_start <  itotdate_cond(nf+1) ) then
         ncondfl=nf
         exit
      endif
   enddo

   print*,'cond starting at history file:',ncondfl

   ! Calculate varweights just like var init
   call varweight(nnzp(1),nnxp(1),nnyp(1),varinit_g(1)%varwts_uv,varinit_g(1)%varwts_th    &
                 ,varinit_g(1)%varwts_pi,varinit_g(1)%varwts_rt,varinit_g(1)%varwts_co2    &
                 ,varinit_g(1)%varwts_rc,grid_g(1)%topt,grid_g(1)%rtgt)                   
   
   ! Read and interpolate files to new grid 1. Put stuff in varinit arrays.
   
   call cond_update(0,ncondfl)
   
   do ifm = 2,ngrids     
      call newgrid(ifm)
      call vfintrpf(ifm,1)
   enddo

elseif (initflag == 2) then   ! Runtime file increment

   if ( time >= cond_times(ncondfl+1) ) ncondfl = ncondfl + 1

endif


! Read and interpolate files to new grid 1.

call cond_update(1,ncondfl+1)

! Fill nested grid nudging arrays for grids we didn't fill directly

!do ifm = 2,ngrids
!   if(igrid_match(ifm) == -1) then
!      call newgrid(ifm)
!      call vfintrpf(ifm,2)
!   endif
!enddo

condtime1=cond_times(ncondfl)
condtime2=cond_times(ncondfl+1)

return
end



subroutine cond_file_inv (hfilin,iyear1,imonth1,idate1,itime1)

use mem_varinit
use grid_dims, only : str_len
implicit none

character(len=*) :: hfilin
integer :: iyear1,imonth1,idate1,itime1

integer :: nc,nf,lnf,nhftot
integer :: inyear,inmonth,indate,inhour


character(len=str_len), dimension(maxnudfiles) :: fnames
character(len=str_len) :: hpref,thisfile
character(len=14)  :: itotdate
real(kind=8) :: secs_init,secs_nud

! Get abs seconds of run start

call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through history files and make inventory

nc=len_trim(hfilin)
nhftot=-1
hpref=hfilin(1:nc-26)
thisfile = hpref(1:len_trim(hpref))//'????-??-??-??????-head.txt'
call RAMS_filelist(maxnudfiles,fnames,thisfile,nhftot)

if(nhftot > maxnudfiles) then
   print*,'too many cond history files'
   stop 'lots_of_cond_history'
endif

ncondfiles=0
do nf=1,nhftot
   lnf=len_trim(fnames(nf))
   read(fnames(nf)(lnf-25:lnf-9),20) inyear,inmonth,indate,inhour
   20 format(i4,1x,i2,1x,i2,1x,i6)

   call date_make_big(inyear,inmonth,indate,inhour,itotdate)

   ncondfiles=ncondfiles+1
   fnames_cond(ncondfiles)=fnames(nf)
   itotdate_cond(ncondfiles)=itotdate
   
   call date_abs_secs2(inyear,inmonth,indate,inhour,secs_nud)
   cond_times(ncondfiles)=secs_nud - secs_init

enddo

call RAMS_dintsort(ncondfiles,itotdate_cond,fnames_cond)

!  start printing section
!--------------------------------------------------------------

print*,' '
print*,' '
print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  Condensate Nudging Input File Inventory --------'
print*,'-------------------------------------------------------------'
do nf=1,ncondfiles
   print*,  itotdate_cond(nf),'   ',cond_times(nf)  &
           ,fnames_cond(nf)(1:len_trim(fnames_cond(nf)))
enddo
print*,'------------------------------------------------------'

!--------------------------------------------------------------


return
end


