!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine top_read(ifm)

use mem_grid
use io_params

implicit none

integer :: ifm,i,j

character(len=128) :: flnm
character(len=2) :: cgrid
character(len=1) :: dummy
logical :: there

! read the "top" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'TOP_read: file for grid ',ifm,' not there.'
   print*,'TOP_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'top_read: no file'
endif

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! Skip header records (already been checked)
read(25,*) dummy
read(25,*) dummy
read(25,*) dummy

call vfirec(25,grid_g(ifm)%topta(1,1),nnxp(ifm)*nnyp(ifm),'LIN')
call vfirec(25,grid_g(ifm)%topzo(1,1),nnxp(ifm)*nnyp(ifm),'LIN')

close (25)

return
end

!******************************************************************

subroutine top_check(ifm,ierr)

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

integer :: lc,isfc_marker,isfc_ver,nsfx,nsfy  &
   ,nsitoptflg,nsitopsflg,nsiz0flg
real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,stoptenh,stoptwvl  &
   ,sz0max,sz0fact,glatr,glonr

character(len=128) :: flnm
character(len=2) :: cgrid
logical there

lc=len_trim(topfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'TOPfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)

read (25,*) isfc_marker,isfc_ver
read (25,*) nsfx,nsfy  &
        ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
read (25,*) nsitoptflg  &
        ,nsitopsflg,stoptenh,stoptwvl,nsiz0flg,sz0max,sz0fact
close (25)


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
    abs(sfplat-platn(ifm))     .gt. .001          .or.  &
    abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsitoptflg                 .ne. itoptflg(ifm) .or.  &
    nsitopsflg                 .ne. itopsflg(ifm) .or.  &
    abs(stoptenh-toptenh(ifm)) .gt. .001          .or.  &
    abs(stoptwvl-toptwvl(ifm)) .gt. .001          .or.  &
    nsiz0flg                   .ne. iz0flg(ifm)   .or.  &
    abs(sz0max-z0max(ifm))     .gt. .001          .or.  &
    abs(sz0fact-z0fact)        .gt. .00001) then

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
   print*,'itoptflg:',itoptflg(ifm),nsitoptflg
   print*,'itopsflg:',itopsflg(ifm),nsitopsflg
   print*,'toptenh:',toptenh(ifm),stoptenh
   print*,'toptwvl:',toptwvl(ifm),stoptwvl
   print*,'iz0flg:',iz0flg(ifm),nsiz0flg
   print*,'z0max:',z0max(ifm),sz0max
   print*,'z0fact:',z0fact,sz0fact
   print*,'-------------------'

else

   ierr = 0

endif

return
end

! ****************************************************************************

subroutine toptinit(n2,n3,ifm,topt,topzo)
implicit none
integer :: n2,n3,ifm,i,j
real, dimension(n2,n3) :: topt,topzo

! Fill the TOPT array with a default value of 0.  This default is used only
! when a standard RAMS topography dataset is not used and when no overrides
! to topography heights are defined in subroutine toptinit_user in the
! file ruser.f.

do j = 1,n3
   do i = 1,n2
      topt(i,j) = 0.
      topzo(i,j) = .0001
   enddo
enddo
return
end


!*****************************************************************************

subroutine top_write(ifm)

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

flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open (25,flnm,'FORMATTED','REPLACE','WRITE',1)
rewind 25

write(25,99) 999999,3
99   format(2i8)

write(25,100) nnxp(ifm),nnyp(ifm)  &
   ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
   ,glatr,glonr
100  format(2i5,2f15.5,4f11.5)

write(25,101) itoptflg(ifm),itopsflg(ifm),toptenh(ifm),toptwvl(ifm)  &
   ,iz0flg(ifm),z0max(ifm),z0fact
101  format(2i5,2f11.5,i5,2f11.5)


call vforec(25,sfcfile_p(ifm)%topt(1,1),nnxyp(ifm),24,scrx,'LIN')
call vforec(25,sfcfile_p(ifm)%topzo(1,1),nnxyp(ifm),24,scrx,'LIN')

close(25)

return
end
