!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!TEB
subroutine fuso_read(ifm)

use mem_grid
use mem_teb
use mem_gaspart
use io_params
  USE teb_vars_const, only: iteb
  USE mem_emiss, only : ichemi,isource !for gas emission 

implicit none

integer :: ifm,i,j

character(len=128) :: flnm
character(len=2) :: cgrid
character(len=1) :: dummy
logical :: there


! read the "fuso" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(fusfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'FUSO_read: file for grid ',ifm,' not there.'
   print*,'FUSO_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'fuso_read: no file'
endif

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! Skip header records (already been checked)
read(25,*) dummy
read(25,*) dummy
read(25,*) dummy

if(iteb==1) &
call vfirec(25,teb_g(ifm)%fuso,nnxp(ifm)*nnyp(ifm),'LIN')

if(isource==1) &
gaspart_g(ifm)%fusog(1,1)=teb_g(ifm)%fuso(1,1)

close (25)

return
end

!******************************************************************

!TEB
subroutine fuso_check(ifm,ierr)

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
   ,nsifusoflg
real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr

character(len=128) :: flnm
character(len=2) :: cgrid
logical there

lc=len_trim(fusfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(fusfiles)//'-S-'//cgrid//'.vfm'

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'FUSOfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)

read (25,*) isfc_marker,isfc_ver
read (25,*) nsfx,nsfy  &
        ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
read (25,*) nsifusoflg  
close (25)


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
    abs(sfplat-platn(ifm))     .gt. .001          .or.  &
    abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsifusoflg                 .ne. ifusflg(ifm)) then

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
   print*,'ifusflg:',ifusflg(ifm),nsifusoflg
   print*,'-------------------'

else

   ierr = 0

endif

return
end

! ****************************************************************************

! TEB
subroutine fusoinit(n2,n3,ifm,fuso,glon)
implicit none
integer :: n2,n3,ifm,i,j
real, dimension(n2,n3) :: fuso,glon

! Fill the FUSO array with a default value equivalent to Longitude difference
! from Grenwich meridian .  This default is used only
! when a standard RAMS fuso dataset is not used.

do j = 1,n3
   do i = 1,n2
      fuso(i,j) = float(int((glon(i,j)-7.5)/15))
   enddo
enddo
return
end

! ****************************************************************************

! TEB
subroutine fuso_write(ifm)

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

flnm=trim(fusfiles)//'-S-'//cgrid//'.vfm'

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open (25,flnm,'FORMATTED','REPLACE','WRITE',1)
rewind 25

write(25,99) 999999,3
99   format(2i8)

write(25,100) nnxp(ifm),nnyp(ifm)  &
   ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
   ,glatr,glonr
100  format(2i5,2f15.5,4f11.5)

write(25,101) ifusflg(ifm)
101  format(i5)

call vforec(25,sfcfile_p(ifm)%fuso,nnxyp(ifm),24,scrx,'LIN')

close(25)

return
end
