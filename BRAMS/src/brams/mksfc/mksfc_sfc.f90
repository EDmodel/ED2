!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine sfc_read(ifm)

use mem_grid
use mem_leaf
use io_params

implicit none

integer :: ifm,ifileok,ipat,k,i,j
logical :: there

character(len=128) :: flnm
character(len=2) :: cgrid
character(len=1) :: dummy

! read the "sfc" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'SFC_read: file for grid ',ifm,' not there.'
   print*,'SFC_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'sfc_read: no file'
endif

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! Skip header records (already been checked)
read(25,*) dummy
read(25,*) dummy
read(25,*) dummy

do ipat = 1,npatch
   call vfirec(25,leaf_g(ifm)%patch_area(:,:,ipat),nnxp(ifm)*nnyp(ifm),'LIN')
enddo

do ipat = 1,npatch
   call vfirec(25,leaf_g(ifm)%leaf_class(:,:,ipat),nnxp(ifm)*nnyp(ifm),'LIN')
enddo

do ipat = 1,npatch
   call vfirec(25,leaf_g(ifm)%soil_text(:,:,:,ipat),nzg*nnxp(ifm)*nnyp(ifm),'LIN')
enddo

close (25)

return
end

!******************************************************************

subroutine sfc_check(ifm,ierr)

use mem_grid
use io_params

! This subroutine checks for the existence of a surface file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ifileok is returned with a value of 0.  If the file 
! exists and is ok, ifileok is returned with a value of 1.

implicit none

integer :: ifm,ierr

integer :: lc,isfc_marker,isfc_ver,nsfx,nsfy,nsfzg  &
   ,nsivegtflg,nsisoilflg,nsnofilflg,nspatch
real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr

character(len=128) :: flnm
character(len=2) :: cgrid
logical there

lc=len_trim(sfcfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'SFCfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)

read (25,*) isfc_marker,isfc_ver
read (25,100) nsfx,nsfy,nsfzg,nspatch  &
        ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
read (25,101) nsivegtflg,nsisoilflg,nsnofilflg
100  format(4i5,2f15.5,4f11.5)
101  format(5i5,2f11.5,i5,2f11.5)
close (25)


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    nsfzg                      .ne. nzg           .or.  &
    nspatch                    .ne. npatch        .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
    abs(sfplat-platn(ifm))     .gt. .001          .or.  &
    abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsivegtflg                 .ne. ivegtflg(ifm) .or.  &
    nsisoilflg                 .ne. isoilflg(ifm) .or.  &
    nsnofilflg                 .ne. nofilflg(ifm) ) then

   ierr = 1

   print*,'SFCfile mismatch on grid:',ifm
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
   print*,'ivegtflg:',ivegtflg(ifm),nsivegtflg
   print*,'isoilflg:',isoilflg(ifm),nsisoilflg
   print*,'nofilflg:',nofilflg(ifm),nsnofilflg
   print*,'-------------------'

else

   ierr = 0

endif

return
end

!*****************************************************************************

subroutine sfc_write(ifm)

use mem_mksfc
use mem_grid
use io_params

implicit none

integer :: ifm,ip,k,i,j
real :: glatr,glonr
character(len=128) :: flnm
character(len=2) :: cgrid

!     write surface characteristics, one file for each grid


write(cgrid,'(a1,i1)') 'g',ifm

flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open (25,flnm,'FORMATTED','REPLACE','WRITE',1)
rewind 25

write(25,99) 999999,3
99   format(2i8)

write(25,100) nnxp(ifm),nnyp(ifm),nzg,npatch  &
   ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
   ,glatr,glonr
100  format(4i5,2f15.5,4f11.5)

write(25,101) ivegtflg(ifm),isoilflg(ifm),nofilflg(ifm)
101  format(3i5)


do ip = 1,npatch
   call vforec(25,sfcfile_p(ifm)%patch_area(:,:,ip),nnxyp(ifm),24,scrx,'LIN')
enddo

do ip = 1,npatch
   call vforec(25,sfcfile_p(ifm)%leaf_class(:,:,ip),nnxyp(ifm),24,scrx,'LIN')
enddo

do ip = 1,npatch
   call vforec(25,sfcfile_p(ifm)%soil_text(:,:,:,ip),nzg*nnxyp(ifm),24,scrx,'LIN')
enddo

close(25)

return
end
