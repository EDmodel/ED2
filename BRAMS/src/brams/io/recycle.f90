!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine recycle()

use mem_grid
use mem_leaf
use mem_scratch
use var_tables
use io_params
!srf -for carma AOT recycle
  use mem_aerad, only: nwave !INTENT(IN)

use mem_cuparm, only: nclouds ! intent(in)
implicit none

character(len=128) :: flnm
!integer :: idtype,ng,nvars,lenf
integer :: ng,nvars,lenf,ierr
integer, external :: RAMS_getvar

flnm=pastfn(1:len_trim(pastfn)-9)

print*,'Reading assimilation fields from analysis file ',pastfn

call rams_read_header(flnm(1:len_trim(flnm)))

! read the requested analysis file variables

do ng=1,ngrids


   DO nvars=1,num_var(ng)
   
      if(vtab_r(nvars,ng)%irecycle == 1) then

         print*,'Reading assimilation field:',vtab_r(nvars,ng)%name  &
                                ,' for grid:',ng                     &
				,' dim:',vtab_r(nvars,ng)%idim_type  &
				,' npts:', vtab_r(nvars,ng)%npts

         ierr=RAMS_getvar(vtab_r(nvars,ng)%name,ng  &
                         ,scratch%scr1(1),scratch%scr2(1),flnm)


         if(vtab_r(nvars,ng)%idim_type == 4) then
            call unarrange_p(nnxp(ng),nnyp(ng),nzg,npatch  &
                         ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)
         elseif(vtab_r(nvars,ng)%idim_type == 5) then
            call unarrange_p(nnxp(ng),nnyp(ng),nzs,npatch  &
                         ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)

!srf
!use this for 3d atmospheric fields
         elseif(vtab_r(nvars,ng)%idim_type == 3) then
            call unarrange(nnzp(ng),nnxp(ng),nnyp(ng)  &
                         ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)
!srf
!use this for 3d (NX,NY,NWAVE) CARMA AOT fields
	 elseif(vtab_r(nvars,ng)%idim_type == 7) then
	    call rearrange_aot(nwave,nnxp(ng),nnyp(ng)  &
			 ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)

! Cumulus parameterization variables with profile and spectrum
         elseif(vtab_r(nvars,ng)%idim_type == 8) then
            call unarrange_p(nnxp(ng),nnyp(ng),nnzp(ng),nclouds  &
                         ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)

!use this for leaf 3, horizontal cumulus variables, or 2 dim:
         else
            call atob(vtab_r(nvars,ng)%npts  &
                     ,scratch%scr1(1),vtab_r(nvars,ng)%var_p)

         endif
             
      endif

   enddo

enddo
return
end

subroutine rearrange_aot(nwave, nx, ny, array_i, array_o)
   implicit none
   integer, intent(in) :: nx, ny, nwave
   real, intent(in)    :: array_i(nx,ny,nwave)
   real, intent(out)   :: array_o(nx,ny,nwave)
   
   integer :: i, j, w
   
   do w=1, nwave
      do j=1, ny
         do i=1, nx
	    array_o(i,j,w) = array_i(i,j,w)
	    !if(w==12 .and. array_o(i,j,w) > 0.01) print*,i,j,array_o(i,j,w)
	 enddo
      enddo
   enddo
!       open(19,file='aot.gra',         &
!            form='unformatted',access='direct',status='unknown',  &
!            recl=4*nx*ny)	  
!       nrec=1
!       write(19,rec=nrec) array_o(:,:,12)
!       close (19)
!       stop 333

end subroutine rearrange_aot

!*******************************************************************************

subroutine unarrange_p(n2,n3,n4,n5,a,b)
implicit none

integer :: n2,n3,n4,n5
real :: a(n2,n3,n4,n5),b(n4,n2,n3,n5)

integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(k,i,j,ip) = a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo
return
end
