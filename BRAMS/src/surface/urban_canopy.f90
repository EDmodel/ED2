!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine urban_canopy()

use mem_grid
use mem_turb
use mem_basic
use mem_scratch
use mem_tend

use node_mod

implicit none

! ut = -c abs(spd) u
! vt = -c abs(spd) v


call urb_tend(mzp,mxp,myp,ia,iz,ja,jz,jdim         &
             ,basic_g(ngrid)%up,basic_g(ngrid)%vp  &
             ,turb_g(ngrid)%cdrag,scratch%scr2     &
             ,tend_g(ngrid)%ut,tend_g(ngrid)%vt,dtlv )
             
return
end

!************************************************************************

subroutine urb_tend(m1,m2,m3,ia,iz,ja,jz,jdim,up,vp,cdrag,spd,ut,vt,dtlv)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jdim 
real, dimension(m1,m2,m3) :: up,vp,cdrag,spd,ut,vt
real :: dtlv

integer :: i,j,k

! Speed is computed at T points, then multiplied by cdrag. This is 
!     averaged to u,v points.
!   Parallel considertion: need proper u at i+1,i-1, v at j+1,j-1

do j=2,m3
   do i=2,m2
      do k=2,m1-1
         spd(k,i,j)= sqrt(  (.5*(up(k,i,j)+up(k,i-1,j)))**2  &
                           +(.5*(vp(k,i,j)+vp(k,i,j-jdim)))**2 )  &
                     *cdrag(k,i,j)
      enddo
   enddo
enddo



do j=ja,jz
   do i=ia,iz
      do k=2,m1-1
         ut(k,i,j) = ut(k,i,j) - &
            min( (spd(k,i,j)+spd(k,i+1,j))*.5*up(k,i,j)  &
                 , .3 * up(k,i,j)/dtlv   )
         vt(k,i,j) = vt(k,i,j) - &
            min( (spd(k,i,j)+spd(k,i,j+jdim))*.5*vp(k,i,j) &
                 , .3 * vp(k,i,j)/dtlv   )
         if(i==41.and.j==42.and.k==2) &
         print*,ut(k,i,j),up(k,i,j),spd(k,i,j)/cdrag(k,i,j),cdrag(k,i,j)
      enddo
   enddo
enddo

return
end

!************************************************************************

subroutine urb_drag_init()

use mem_turb
use mem_grid

implicit none

integer :: ng,k
real, allocatable :: plt(:,:)

! Compute drag coefficient from bfrac (building fraction), or whatever
!    else we put in it. Assume this will be zero if not in urban grid cell.

!  c = b*sqrt(A) / (2*(A-a))

! cdrag to be computed at T points, set to zero for non-urban cells


do ng = 1,ngrids

   call getdrag(nnzp(ng),nnxp(ng),nnyp(ng),xmn(:,ng),ymn(:,ng),zmn(:,ng)  &
               ,turb_g(ng)%cdrag)
   
allocate (plt(nnxp(ng),nnyp(ng)))
   do k=1,10
      plt(1:nnxp(ng),1:nnyp(ng))=turb_g(ng)%cdrag(k,1:nnxp(ng),1:nnyp(ng))
      call ezcntr(plt(1,1),nnxp(ng),nnyp(ng))
   enddo
deallocate(plt)
   
enddo

!call clsgks
!stop

return
end

!************************************************************************

subroutine getdrag(nzp,nxp,nyp,xm,ym,zm,cdrag)
use grid_dims, only : str_len
implicit none

integer :: nxp,nyp,nzp
real, dimension(nxp) :: xm
real, dimension(nyp) :: ym
real, dimension(nzp) :: zm
real, dimension(nzp,nxp,nyp) :: cdrag
real, allocatable, dimension(:,:,:) :: xxdrag

integer :: i,j,k,ic,jc,kc,ic1,ic2,jc1,jc2,kc1,kc2
real :: rim1_cd,rim2_cd,rjm1_cd,rjm2_cd,rkm1_cd,rkm2_cd,ric,rjc,rkc
logical :: there

character(len=str_len) :: fname='./cdrag_data'

! Read in the prepared drag coeff file

inquire(file=fname, exist=there)
if(.not. there) then
   print*,'Urban canopy drag coeff file not found.'
   print*,'File name:',trim(fname)
   stop 'no urban drag file'
endif


allocate (xxdrag(66,84,11))

call rams_c_open(trim(fname)//char(0),'r')
call rams_c_read(4,66*84*11*4,xxdrag(1,1,1))
call rams_c_close()

!do k=1,11
!call ezcntr(xxdrag(1,1,k),66,84)
!enddo
! Use the following polelat,polelon coordinates!!!

! polelat = 38.889
! polelon = -75.

cdrag(1:nzp,1:nxp,1:nyp) = 0.

do j = 2,nyp-1
   do i = 2,nxp-1
      do k = 2,nzp-1
         
         rim1_cd = (xm(i-1) + 500000. - 304527.) / 500. + 1.
         rim2_cd = (xm(i)   + 500000. - 304527.) / 500. + 1.

         rjm1_cd = ym(j-1) / 500. + 32.
         rjm2_cd = ym(j)   / 500. + 32.

         rkm1_cd = zm(k-1) / 5. + 1.
         rkm2_cd = zm(k)   / 5. + 1.
         
         if(i==41.and.j==43.and.k==2) print*,rim1_cd,rim2_cd 
         if(i==41.and.j==43.and.k==2) print*,rjm1_cd,rjm2_cd 
         if(i==41.and.j==43.and.k==2) print*,rkm1_cd,rkm2_cd 
         if(i==41.and.j==43.and.k==2) print*,xm(i),xm(i-1)
         
!         if ((rim1_cd > 1. .and. rim1_cd < 67.) .and.  &
!             (rim2_cd > 1. .and. rim2_cd < 67.) .and.  &
!             (rjm1_cd > 1. .and. rjm1_cd < 85.) .and.  &
!             (rjm2_cd > 1. .and. rjm2_cd < 85.) .and.  &
!             (rkm1_cd > 1. .and. rkm1_cd < 11.) .and.  &
!             (rkm2_cd > 1. .and. rkm2_cd < 11.)  ) then
   
   

            ic1 = max( 1, int(rim1_cd) )
            ic2 = min(66, int(rim2_cd) )
            jc1 = max( 1, int(rjm1_cd) )
            jc2 = min(83, int(rjm2_cd) )
            kc1 = max( 1, int(rkm1_cd) )
            kc2 = min(11, int(rkm2_cd) )
           
            do ic = ic1,ic2 
               do jc = jc1,jc2 
                  do kc = kc1,kc2
                     ric = float(ic)
                     rjc = float(jc)
                     rkc = float(kc)
                     cdrag(k,i,j) = cdrag(k,i,j)  &
                        +  (min(ric+1.,rim2_cd) - max(ric,rim1_cd))  &
                         * (min(rjc+1.,rjm2_cd) - max(rjc,rjm1_cd))  &
                         * (min(rkc+1.,rkm2_cd) - max(rkc,rkm1_cd))  &
                         /  ((rim2_cd - rim1_cd)  &
                           * (rjm2_cd - rjm1_cd)  &
                           * (rkm2_cd - rkm1_cd))  * xxdrag(ic,jc,kc)
                  enddo
               enddo
            enddo
            
!         endif
         
      enddo
   enddo
enddo

deallocate (xxdrag)

return
end
