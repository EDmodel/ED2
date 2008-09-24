!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
                      
subroutine first_RAMS(np1,np2,np3,ui2,vi2,pi2,ti2,ri2)

use an_header
use isan_coms
use mem_grid

implicit none

integer :: np1,np2,np3
real, dimension(np1,np2,np3) ::  ui2,vi2,pi2,ti2,ri2

integer, dimension(maxgrds) :: nnxpx,nnypx,nnzpx,nstratxx,nstratyx  &
                              ,ngbegunx,nxtnestx,nnsttopx,nnstbotx  &
                              ,ninestx,njnestx,nknestx
integer :: lenf,nv,irw,iun,ierr,ng,ngridsx,icm,ifm,k

character(len=80) :: flnma
character(len=1) :: cgrid
character(len=2) :: cng

integer, external :: cio_i,cio_f,cio_i_sca,RAMS_getvar

real, allocatable :: plt(:,:)
real, allocatable :: pltc(:,:)

! Open analysis file and read in commons

lenf=len_trim(innpr)
flnma=innpr(1:lenf)
open(10,file=flnma,form='formatted',status='old')

read(10,*) nvbtab
if(allocated(anal_table)) deallocate(anal_table)
allocate (anal_table(nvbtab))
do nv=1,nvbtab
   read(10,*)  anal_table(nv)%string   &
              ,anal_table(nv)%npointer  &
              ,anal_table(nv)%idim_type  &
              ,anal_table(nv)%ngrid  &
              ,anal_table(nv)%nvalues
enddo

! Reform filename to reflect actual analysis file rather than head file
lenf=lenf-11
write(cgrid,'(i1)') ngrid
flnma=innpr(1:lenf)//'g'//cgrid
lenf=lenf+2

! get stuff from the 1st guess run
irw=1
iun=10
ierr=cio_i_sca(iun,irw,'ngrids',ngridsx,1)
ierr=cio_i(iun,irw,'nnxp',nnxpx,ngridsx)
ierr=cio_i(iun,irw,'nnyp',nnypx,ngridsx)
ierr=cio_i(iun,irw,'nnzp',nnzpx,ngridsx)
ierr=cio_i(iun,irw,'nstratx',nstratxx,ngridsx)
ierr=cio_i(iun,irw,'nstraty',nstratyx,ngridsx)
ierr=cio_i(iun,irw,'nxtnest',nxtnestx,ngridsx)
ierr=cio_i(iun,irw,'nnsttop',nnsttopx,ngridsx)
ierr=cio_i(iun,irw,'nnstbot',nnstbotx,ngridsx)
ierr=cio_i(iun,irw,'ninest',ninestx,ngridsx)
ierr=cio_i(iun,irw,'njnest',njnestx,ngridsx)
ierr=cio_i(iun,irw,'nknest',nknestx,ngridsx)

do ng=1,ngridsx
   write(cng,1) ng
   1 format(i2.2)

   !! 1d reference state filled in the old arrays
   !!ierr=cio_f(iun,irw,'u01dn'//cng,u01dn(1,ng),nnzp(ng))
   !ierr=cio_f(iun,irw,'v01dn'//cng,v01dn(1,ng),nnzp(ng))
   !ierr=cio_f(iun,irw,'pi01dn'//cng,pi01dn(1,ng),nnzp(ng))
   !ierr=cio_f(iun,irw,'th01dn'//cng,th01dn(1,ng),nnzp(ng))
   !ierr=cio_f(iun,irw,'dn01dn'//cng,dn01dn(1,ng),nnzp(ng))
   !ierr=cio_f(iun,irw,'rt01dn'//cng,rt01dn(1,ng),nnzp(ng))

   ! 1d reference state filled in the new arrays
   ierr=cio_f(iun,irw,'pi01dn'//cng,piref(1,ng),nnzp(ng))
   ierr=cio_f(iun,irw,'th01dn'//cng,thref(1,ng),nnzp(ng))
   ierr=cio_f(iun,irw,'dn01dn'//cng,dnref(1,ng),nnzp(ng))
   ierr=cio_f(iun,irw,'rt01dn'//cng,rtref(1,ng),nnzp(ng))
enddo

close(10)

! check that grids are compatible
! would also like to check pole lat lon and topo (and more?)
do ng=1,min( nigrids, ngridsx )
   if(nnxp(ng).ne.nnxpx(ng).or.  &
      nnyp(ng).ne.nnypx(ng).or.  &
      nnzp(ng).ne.nnzpx(ng).or.  &
      nstratx(ng).ne.nstratxx(ng).or.  &
      nstraty(ng).ne.nstratyx(ng).or.  &
      nxtnest(ng).ne.nxtnestx(ng).or.  &
      nnsttop(ng).ne.nnsttopx(ng).or.  &
      nnstbot(ng).ne.nnstbotx(ng).or.  &
      ninest(ng).ne.ninestx(ng).or.  &
      njnest(ng).ne.njnestx(ng).or.  &
      nknest(ng).ne.nknestx(ng)) then
      print*,'New run and 1st guess grids incompatible',ng
      print*,'            new      1st guess'
      print*,'nnxp   ',nnxp(ng),nnxpx(ng)
      print*,'nnyp   ',nnyp(ng),nnypx(ng)
      print*,'nnzp   ',nnzp(ng),nnzpx(ng)
      print*,'nstratx',nstratx(ng),nstratxx(ng)
      print*,'nstraty',nstraty(ng),nstratyx(ng)
      print*,'nxtnest',nxtnest(ng),nxtnestx(ng)
      print*,'nnsttop',nnsttop(ng),nnsttopx(ng)
      print*,'nnstbot',nnstbot(ng),nnstbotx(ng)
      print*,'ninest ',ninest(ng),ninestx(ng)
      print*,'njnest ',njnest(ng),njnestx(ng)
      print*,'nknest ',nknest(ng),nknestx(ng)
      stop 'stop in first_RAMS'
   endif
enddo

! Read analysis files if they exist, else interpolate from the
!   parent grid of the new grid.
! Store everything in the 'A' array as nzp,nxp,nyp and things in the
!   analysis arrays as nxp,nyp,nzp.


if(ngrid <= ngridsx)then

   print*,'READING from file current grid=',ngrid,nxp,nyp,nzp
   print*,'READING from file current grid=',trim(innpr)
   ! UE_AVG and VE_AVG
   ierr=RAMS_getvar('UP',ngrid,ui2(1,1,1),rr_scr1(1),innpr(1:lenf))
   ierr=RAMS_getvar('VP',ngrid,vi2(1,1,1),rr_scr1(1),innpr(1:lenf))
   call comp_avgu(nnxp(ngrid),nnyp(ngrid),nnzp(ngrid),ui2(1,1,1))
   call comp_avgv(nnxp(ngrid),nnyp(ngrid),nnzp(ngrid),vi2(1,1,1))
print*,'222'
   ! RELHUM and THETA
   ierr=RAMS_getvar('RV',ngrid,ri2(1,1,1),rr_scr1(1),innpr(1:lenf))
   ierr=RAMS_getvar('PI',ngrid,pi2(1,1,1),rr_scr1(1),innpr(1:lenf))
   ierr=RAMS_getvar('THETA',ngrid,ti2(1,1,1),rr_scr1(1),innpr(1:lenf))
   call comp_rhfrac(nnxp(ngrid),nnyp(ngrid),nnzp(ngrid)  &
                   ,ri2(1,1,1),pi2(1,1,1),ti2(1,1,1))
print*,'333'

   ! PRESS
   call comp_press(nnxp(ngrid),nnyp(ngrid),nnzp(ngrid),pi2(1,1,1))
   call ae1t0(nnxp(ngrid)*nnyp(ngrid)*nnzp(ngrid),pi2(1,1,1),pi2(1,1,1),100.)
print*,'444'

   ! Write analysis fields to model vars in case need to interpolate from them
   ! rearrange(nzp,nxp,nyp,a,b) is b(i,j,k)=a(k,i,j)
   ! unarrange(nzp,nxp,nyp,a,b) is b(k,i,j)=a(i,j,k)
   call unarrange(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)  &
                  ,ui2(1,1,1),is_grids(ngrid)%rr_ug(1,1,1))
   call unarrange(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)  &
                  ,vi2(1,1,1),is_grids(ngrid)%rr_vg(1,1,1))
   call unarrange(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)  &
                  ,ri2(1,1,1),is_grids(ngrid)%rr_rg(1,1,1))
   call unarrange(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)  &
                  ,ti2(1,1,1),is_grids(ngrid)%rr_tg(1,1,1))
   call unarrange(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)  &
                  ,pi2(1,1,1),is_grids(ngrid)%rr_pg(1,1,1))
   print*,"ngrid,theta(1,1,1)=",ngrid,is_grids(ngrid)%rr_tg(1,1,1)
print*,'555'

   ! Compute reference state density in case we need it
   
   
   if (ngrid == 1) then
      call isan_comp_dn0(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid) &
                        ,is_grids(ngrid)%rr_pi0(1,1,1) &
                        ,is_grids(ngrid)%rr_th0(1,1,1) &
                        ,is_grids(ngrid)%rr_dn0(1,1,1) &
                        ,is_grids(ngrid)%rr_dn0u(1,1,1) &
                        ,is_grids(ngrid)%rr_dn0v(1,1,1) &
                        ,grid_g(ngrid)%topta(1,1),ngrid)
   else
   
      ifm=ngrid
      icm=nxtnest(ifm)
      call nest_interpolated_topo(nnxp(icm),nnyp(icm),nnxp(ifm),nnyp(ifm) &
            ,maxix,maxiy,ifm,grid_g(icm)%topt(1,1),rr_vt2da(1)  &
            ,rr_scr1(1),rr_scr2(1))
      call fmrefs3d_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
            ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy  &
            ,nnstbot(ifm),nnsttop(ifm),jdim  &
            ,rr_scr1(1),rr_scr2(1),rr_vt2da(1)  &
            ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
            ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
            ,is_grids(icm)%rr_th0(1,1,1),is_grids(ifm)%rr_th0(1,1,1) &
            ,is_grids(ifm)%rr_pi0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
            ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )
         
      call fmdn0_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
            ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy &
            ,rr_scr1(1),rr_scr2(1)  &
            ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
            ,is_grids(ifm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
            ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )
   endif
   
else

   ifm=ngrid
   icm=nxtnest(ifm)
   print*,'INTERPOLATING to current grid=',ifm,' from parent grid=',icm
   print*,'INTERPOLATING to current grid=',trim(innpr)
   !!!!!!!!!!!!
   allocate (plt(nnxp(ifm),nnyp(ifm)))
   allocate (pltc(nnxp(icm),nnyp(icm)))
   !!!!!!!!!!!!!!

   print*,'..icm,theta(1,1,1)=',icm,is_grids(icm)%rr_tg(1,1,1)
   print*,'..ifm,theta(1,1,1)=',ifm,is_grids(ifm)%rr_tg(1,1,1)

   ! Calculate the 3D base state
   call fmrefs1d_isan(ifm,icm,maxsigz,nnzp(ifm) &
                     ,piref(1,1),thref(1,1),dnref(1,1),rtref(1,1))

   call nest_interpolated_topo(nnxp(icm),nnyp(icm),nnxp(ifm),nnyp(ifm) &
         ,maxix,maxiy,ifm,grid_g(icm)%topt(1,1),rr_vt2da(1)  &
         ,rr_scr1(1),rr_scr2(1))
   call fmrefs3d_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
         ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy  &
         ,nnstbot(ifm),nnsttop(ifm),jdim  &
         ,rr_scr1(1),rr_scr2(1),rr_vt2da(1)  &
         ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
         ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
         ,is_grids(icm)%rr_th0(1,1,1),is_grids(ifm)%rr_th0(1,1,1) &
         ,is_grids(ifm)%rr_pi0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
         ,is_grids(ifm)%rr_dn0v(1,1,1)  &
         ,ztn(1,ifm),ztop )
         

   call fmint4_isan(is_grids(icm)%rr_tg(1,1,1)  ,is_grids(ifm)%rr_tg(1,1,1)  &
                   ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
                   ,ifm,icm,'t',1)
   print*,"....icm,theta(1,1,1)=",icm,is_grids(icm)%rr_tg(1,1,1),is_grids(icm)%rr_dn0(1,1,1)
   print*,"....ifm,theta(1,1,1)=",ifm,is_grids(ifm)%rr_tg(1,1,1),is_grids(ifm)%rr_dn0(1,1,1)
   call rearrange(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                 ,is_grids(ifm)%rr_tg(1,1,1),ti2(1,1,1))
   print*,'done with th,ti2(1,1,1)=',ti2(1,1,1)

   pltc(1:nnxp(icm),1:nnyp(icm)) = is_grids(icm)%rr_tg(2,1:nnxp(icm),1:nnyp(icm))
   call ezcntr(pltc(1,1),nnxp(icm),nnyp(icm))
   call ezcntr(ti2(1,1,2),nnxp(ifm),nnyp(ifm))


   call fmint4_isan(is_grids(icm)%rr_pg(1,1,1),is_grids(ifm)%rr_pg(1,1,1)  &
                   ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
                   ,ifm,icm,'t',0)
   print*,"....icm,p(1,1,1)=",icm,is_grids(icm)%rr_pg(1,1,1),is_grids(icm)%rr_dn0(1,1,1)
   print*,"....ifm,p(1,1,1)=",ifm,is_grids(ifm)%rr_pg(1,1,1),is_grids(ifm)%rr_dn0(1,1,1)
   print*,"....icm,p(1,1,1)max=",maxval(is_grids(icm)%rr_pg(1:nnzp(icm),1:nnxp(icm),1:nnyp(icm)))
   print*,"....ifm,p(1,1,1)max=",maxval(is_grids(ifm)%rr_pg(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)))
   print*,"....icm,p(1,1,1)min=",minval(is_grids(icm)%rr_pg(1:nnzp(icm),1:nnxp(icm),1:nnyp(icm)))
   print*,"....ifm,p(1,1,1)min=",minval(is_grids(ifm)%rr_pg(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)))
   call rearrange(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                 ,is_grids(ifm)%rr_pg(1,1,1),pi2(1,1,1))


   call fmint4_isan(is_grids(icm)%rr_rg(1,1,1),  is_grids(ifm)%rr_rg(1,1,1)  &
                   ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
                   ,ifm,icm,'t',1)
   print*,"....icm,rr(1,1,1)=",icm,is_grids(icm)%rr_r(1,1,1),is_grids(icm)%rr_dn0(1,1,1)
   print*,"....ifm,rr(1,1,1)=",ifm,is_grids(ifm)%rr_r(1,1,1),is_grids(ifm)%rr_dn0(1,1,1)
   print*,"....icm,rr(1,1,1)max=",maxval(is_grids(icm)%rr_rg(1:nnzp(icm),1:nnxp(icm),1:nnyp(icm)))
   print*,"....ifm,rr(1,1,1)max=",maxval(is_grids(ifm)%rr_rg(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)))
   print*,"....icm,rr(1,1,1)min=",minval(is_grids(icm)%rr_rg(1:nnzp(icm),1:nnxp(icm),1:nnyp(icm)))
   print*,"....ifm,rr(1,1,1)min=",minval(is_grids(ifm)%rr_rg(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)))
   call rearrange(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                      ,is_grids(ifm)%rr_rg(1,1,1),ri2(1,1,1))




   call fmint4_isan(is_grids(icm)%rr_ug(1,1,1),is_grids(ifm)%rr_ug(1,1,1)  &
                   ,is_grids(icm)%rr_dn0u(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
                   ,ifm,icm,'u',1)
   call rearrange(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                 ,is_grids(ifm)%rr_ug(1,1,1),ui2(1,1,1))


   call fmint4_isan(is_grids(icm)%rr_vg(1,1,1),is_grids(ifm)%rr_vg(1,1,1)  &
                   ,is_grids(icm)%rr_dn0v(1,1,1),is_grids(ifm)%rr_dn0v(1,1,1) &
                   ,ifm,icm,'v',1)
   call rearrange(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                 ,is_grids(ngrid)%rr_vg(1,1,1),vi2(1,1,1))

   call fmdn0_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
         ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy &
         ,rr_scr1(1),rr_scr2(1)  &
         ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
         ,is_grids(ifm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
         ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )
endif


return
end

!***************************************************************************

subroutine comp_avgu(n1,n2,n3,a)
implicit none
integer :: n1,n2,n3
real :: a(n1,n2,n3)
integer :: i,j,k
do k=1,n3
   do j=1,n2
      do i=n1,2,-1
         a(i,j,k)=0.5*(a(i,j,k)+a(i-1,j,k))
      enddo
   enddo
enddo
return
end

!***************************************************************************

subroutine comp_avgv(n1,n2,n3,a)
implicit none
integer :: n1,n2,n3
real :: a(n1,n2,n3)
integer :: i,j,k
do k=1,n3
   do j=n2,2,-1
      do i=1,n1
         a(i,j,k)=0.5*(a(i,j,k)+a(i,j-1,k))
      enddo
   enddo
enddo
return
end


!***************************************************************************

subroutine comp_rhfrac(n1,n2,n3,a,b,c)

use rconstants
use therm_lib, only: rehuil
implicit none
integer :: n1,n2,n3
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
real :: xtemp,xpress
integer :: i,j,k
do k=1,n3
   do j=1,n2
      do i=1,n1
         xtemp=c(i,j,k)*b(i,j,k)/cp
         xpress=(b(i,j,k)/cp)**cpor*p00
         a(i,j,k)=min(1.,rehuil(xpress,xtemp,a(i,j,k)))
      enddo
   enddo
enddo
return
end

!***************************************************************************

subroutine comp_press(n1,n2,n3,a)

use rconstants

implicit none
integer :: n1,n2,n3
real :: a(n1,n2,n3)
integer :: i,j,k
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=(a(i,j,k)/cp)**cpor*p00*.01
         enddo
      enddo
   enddo
return
end

!***************************************************************************

subroutine isan_comp_dn0 (n1,n2,n3,pi0,th0,dn0,dn0u,dn0v,topt,ngrd) 

use isan_coms
use rconstants
use mem_grid
use mem_scratch

implicit none

integer :: n1,n2,n3,ngrd,nmax
real :: pi0(n1,n2,n3),th0(n1,n2,n3),dn0(n1,n2,n3),topt(n2,n3)  &
       ,dn0u(n1,n2,n3),dn0v(n1,n2,n3)

integer :: i,j,k
real :: c1,c2,c3

ztop = zmn(nnzp(1)-1,1)
!print*,'ngrd,topt(10,10),ztop,ztn(10,ngrd)=',ngrd,topt(10,10),ztop,ztn(10,ngrd)
do j=1,n3
   do i=1,n2
      do k=1,n1
         vctr2(k)=ztn(k,ngrd)*(1.-topt(i,j)/ztop)+topt(i,j)
        ! if(i == 10 .and. j == 10)  &
        ! print*,'thref,pi01dn,vctr2=',thref(k,ngrd),piref(k,ngrd),vctr2(k)
      enddo
      call htint(n1,piref(1,ngrd),ztn(1,ngrd),n1,vctr11(1),vctr2(1))
      call htint(n1,thref(1,ngrd),ztn(1,ngrd),n1,vctr12(1),vctr2(1))

      do k=1,n1
         th0(k,i,j)=vctr12(k)
      enddo
      pi0(n1,i,j) = vctr11(n1)

      c1=g*2.*(1.-topt(i,j)/ztop)
      c2=(1-cpor)
      c3=cp**c2
      do k=n1-1,1,-1
         pi0(k,i,j)=pi0(k+1,i,j)  &
             +c1/((th0(k,i,j)+th0(k+1,i,j))*dzmn(k,ngrd))
      enddo

      do k=1,n1
         dn0(k,i,j)=(c3*p00)/(rgas*th0(k,i,j)*pi0(k,i,j)**c2)
        ! if(i == 10 .and. j == 10)  &
        ! print*,'th0,pi0,dn0=',th0(k,i,j),pi0(k,i,j),dn0(k,i,j)
      enddo

   enddo
enddo

do k=1,n1
   do j=1,n3
      do i=1,n2-1
         dn0u(k,i,j)=0.5*(dn0(k,i,j)+dn0(k,i+1,j))
      enddo
      dn0u(k,n2,j)=dn0(k,n2,j)
   enddo
enddo

do k=1,n1
   do i=1,n2
      do j=1,n3-1
         dn0v(k,i,j)=0.5*(dn0(k,i,j)+dn0(k,i,j+1))
      enddo
      dn0v(k,i,n3)=dn0(k,i,n3)
   enddo
enddo
   
return
end

!     *****************************************************************

subroutine fmint4_isan(var1,var2,dn0xc,dn0xf,ifm,icm,vpnt,idwt)

! Similar to fmint4 except use ISAN scratch arrays and max points

use isan_coms
use mem_grid

implicit none

integer :: ifm,icm,idwt,lll(3)

real, dimension(*) :: var1,var2,dn0xc,dn0xf
character(len=*) :: vpnt

real, allocatable :: plt(:,:) ,plt3(:,:,:),plt3b(:,:,:)

if (icm == 0) return


call fillscr(maxiz,maxix,maxiy,nnzp(icm),nnxp(icm),nnyp(icm)  &
            ,1,nnzp(icm),rr_scr1(1),var1(1))

if (idwt == 1) then
   call dnswt2(maxiz,maxix,maxiy,nnzp(icm),nnxp(icm),nnyp(icm)  &
              ,rr_scr1(1),dn0xc(1),vpnt,1)
endif

call eintp(rr_scr1(1),rr_scr2(1),maxiz,maxix,maxiy  &
          ,nnzp(ifm),nnxp(ifm),nnyp(ifm),ifm,3,vpnt,0,0)

call fillvar(maxiz,maxix,maxiy,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
            ,1,nnzp(ifm),rr_scr2(1),var2(1))

if (idwt == 1) then
   call dnswt2(nnzp(ifm),nnxp(ifm),nnyp(ifm),nnzp(ifm)  &
              ,nnxp(ifm),nnyp(ifm),var2(1),dn0xf(1),vpnt,2)
endif


call ezcntr(grid_g(icm)%topta(1,1),nnxp(icm),nnyp(icm))
call ezcntr(grid_g(ifm)%topta(1,1),nnxp(ifm),nnyp(ifm))
call ezcntr(rr_vt2da(1),nnxp(ifm),nnyp(ifm))
allocate (plt(nnxp(ifm),nnyp(ifm)))
allocate (plt3(nnzp(ifm),nnxp(ifm),nnyp(ifm)))
allocate (plt3b(nnzp(ifm),nnxp(ifm),nnyp(ifm)))
call ae1m1(nnxp(ifm)*nnyp(ifm),plt(1,1),grid_g(ifm)%topta(1,1),rr_vt2da(1))
call ezcntr(plt(1,1),nnxp(ifm),nnyp(ifm))

! interp field
call ae1(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),plt3(1,1,1),var2(1))
plt(1:nnxp(ifm),1:nnyp(ifm)) =plt3(12,1:nnxp(ifm),1:nnyp(ifm))
call ezcntr(plt(1,1),nnxp(ifm),nnyp(ifm))

!------------------
call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),var2(1),rr_vt2da(1)  &
             ,grid_g(ifm)%topta(1,1),ifm,vpnt)
!------------------

! after rtgint
call ae1(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),plt3(1,1,1),var2(1))
plt(1:nnxp(ifm),1:nnyp(ifm)) =plt3(12,1:nnxp(ifm),1:nnyp(ifm))
call ezcntr(plt(1,1),nnxp(ifm),nnyp(ifm))

! diff field
!call ae1m1(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),plt3b(1,1,1),plt3(1,1,1),var2(1))
!plt(1:nnxp(ifm),1:nnyp(ifm)) =plt3b(12,1:nnxp(ifm),1:nnyp(ifm))
!call ezcntr(plt,nnxp(ifm),nnyp(ifm))


return
end

