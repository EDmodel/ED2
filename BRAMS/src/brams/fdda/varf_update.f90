!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine varf_update(iswap,ifileok,initflag)

use mem_leaf
use mem_varinit
use mem_basic
use mem_grid
use mem_scratch
use therm_lib, only: level

implicit none

integer :: ifileok,initflag,iswap

!---------------------------------------------------------------+
!    "Variable initialization"  initialization routines
!---------------------------------------------------------------+
logical there
integer :: iver_var,nc,iyearx,imonthx,idatex,ihourx  &
          ,nxpx,nypx,nzpx,imarker,iun
real :: rlatx,wlon1x,deltaxx,deltayx,deltazx,dzratx,dzmaxx
character(len=7)   :: cgrid
character(len=128) :: flnm

!      Check and see what we are doing. If it is initial time, read
!        fields into regular arrays. If not, see if nudging will be done
!        on this grid if it is a nested grid.
if (ngrid > 1 .and. tnudcent+tnudtop < .001 .and. initflag == 0) return

! Put new fields into varinit future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then
   varinit_g(ngrid)%varup(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))=  &
      varinit_g(ngrid)%varuf(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))
   varinit_g(ngrid)%varvp(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))=  &
      varinit_g(ngrid)%varvf(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))
   varinit_g(ngrid)%varpp(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))=  &
      varinit_g(ngrid)%varpf(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))
   varinit_g(ngrid)%vartp(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))=  &
      varinit_g(ngrid)%vartf(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))
   varinit_g(ngrid)%varrp(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))=  &
      varinit_g(ngrid)%varrf(1:nnzp(ngrid),1:nnxp(ngrid),1:nnyp(ngrid))
endif

! Make data file name from tag file name
write(cgrid,'(a2,i1,a4)') '-g',ngrid,'.vfm'
nc=len_trim(fnames_varf(nvarffl))
flnm=fnames_varf(nvarffl)(1:nc-4)//trim(cgrid)

! Check for existence
inquire(file=trim(flnm),exist=there)

! Gotta have grid 1...
if (.not.there .and. ngrid == 1) then
   print*
   print*,'No grid 1 varfile found: ',trim(flnm)
   print*
   stop 'no grid 1 varfile'
endif

if(there) then
   ifileok=1
else
   ifileok=0
   return
endif


! Read the varfile fields into the "future" varinit arrays. These will be 
!   swapped to the past arrays when needed.
iun=22

call rams_f_open(iun,flnm,'FORMATTED','OLD','READ',0)

! Find varfile "version"
read(iun,*) imarker
rewind(iun)

if(imarker == 999999) then
   read(iun,*) imarker,iver_var
else
   iver_var=1
endif

read(iun,*) iyearx,imonthx,idatex,ihourx  &
     ,nxpx,nypx,nzpx,rlatx,wlon1x,deltaxx,deltayx,deltazx  &
     ,dzratx,dzmaxx

if(nxp.ne.nxpx.or.  &
   nyp.ne.nypx.or.  &
   nzp.ne.nzpx.or.  &
   abs(deltax-deltaxx).gt..001.or.  &
   abs(deltay-deltayx).gt..001.or.  &
   abs(deltaz-deltazx).gt..001.or.  &
   abs(dzrat-dzratx).gt..001.or.  & 
   abs(dzmax-dzmaxx).gt..001.or.  &
   abs(platn(ngrid)-rlatx).gt..001.or.  &
   abs(plonn(ngrid)-wlon1x).gt..001) then
   
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!    GRID MISMATCH BETWEEN VARFILE AND NAMELIST !'
   print*,'!!          RUN IS STOPPED                       !'
   print*,'!!  File:',trim(flnm)
   print*,'!!  File, Namelist values for grid:',ngrid
   print*,'!!  nxp:',nxpx,nxp
   print*,'!!  nyp:',nypx,nyp
   print*,'!!  nzp:',nzpx,nzp
   print*,'!!  deltax:',deltaxx,deltax
   print*,'!!  deltay:',deltayx,deltay
   print*,'!!  deltaz:',deltazx,deltaz
   print*,'!!  dzrat:',dzratx,dzrat
   print*,'!!  dzmax:',dzmaxx,dzmax
   print*,'!!  polelat:',rlatx,platn(ngrid)
   print*,'!!  polelon:',wlon1x,plonn(ngrid)
   PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'bad-vfile'
endif

call vfirec(iun,varinit_g(ngrid)%varuf,nxyzp,'LIN')
call vfirec(iun,varinit_g(ngrid)%varvf,nxyzp,'LIN')
call vfirec(iun,varinit_g(ngrid)%varpf,nxyzp,'LIN')
call vfirec(iun,varinit_g(ngrid)%vartf,nxyzp,'LIN')
call vfirec(iun,varinit_g(ngrid)%varrf,nxyzp,'LIN')

varinit_g(ngrid)%varrf(1:nzp,1:nxp,1:nyp)=  &
           max(1.e-8,varinit_g(ngrid)%varrf(1:nzp,1:nxp,1:nyp) )
         
if(initflag == 1 .and. iver_var == 2) then
   ! Extract snow depth from the varfile. Ignore other 2D fields for now.
   call vfirec(iun,scratch%vt2da,nxyp,'LIN')
   call vfirec(iun,scratch%vt2da,nxyp,'LIN')
   call vfirec(iun,scratch%vt2da,nxyp,'LIN')
   call vfirec(iun,leaf_g(ngrid)%snow_mass,nxyp,'LIN')
   call vfirec(iun,scratch%vt2da,nxyp,'LIN')
endif

close(iun)

! If running ADAP coord, do interpolation to Cartesian levels

if (if_adap == 1) then
   call varf_adap(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)   &
         ,varinit_g(ngrid)%varuf,varinit_g(ngrid)%varvf &
         ,varinit_g(ngrid)%varpf,varinit_g(ngrid)%vartf &
         ,varinit_g(ngrid)%varrf,grid_g(ngrid)%topta    )
endif

! Find the reference state

if(initflag == 1 .and. ngrid == 1)  &
     call varref(nzp,nxp,nyp &
         ,varinit_g(ngrid)%vartf       , varinit_g(ngrid)%varpf  &
         ,basic_g(ngrid)%pi0           , basic_g(ngrid)%th0      &
         ,varinit_g(ngrid)%varrf       , basic_g(ngrid)%dn0      &
         ,basic_g(ngrid)%dn0u          , basic_g(ngrid)%dn0v     &
         ,varinit_g(ngrid)%varuf       , varinit_g(ngrid)%varvf  &
         ,grid_g(ngrid)%topt           , grid_g(ngrid)%topu      &
         ,grid_g(ngrid)%topv           , grid_g(ngrid)%rtgt      &
         ,grid_g(ngrid)%rtgu           , grid_g(ngrid)%rtgv      &
         ,grid_g(ngrid)%topta          , level                   )

varinit_g(ngrid)%varpf(1:nzp,1:nxp,1:nyp)=  &
           varinit_g(ngrid)%varpf(1:nzp,1:nxp,1:nyp)  &
           - basic_g(ngrid)%pi0(1:nzp,1:nxp,1:nyp)
           
! If this is an initialization, put data into regular arrays

if(initflag == 1 ) then
   call atob(nxyzp,varinit_g(ngrid)%varuf,basic_g(ngrid)%uc  )
   call atob(nxyzp,varinit_g(ngrid)%varvf,basic_g(ngrid)%vc  )
   call atob(nxyzp,varinit_g(ngrid)%varpf,basic_g(ngrid)%pc  )
   call atob(nxyzp,varinit_g(ngrid)%vartf,basic_g(ngrid)%thp )
   call atob(nxyzp,varinit_g(ngrid)%varrf,basic_g(ngrid)%rtp )
endif


return
end

!     **************************************************************

subroutine varf_adap(n1,n2,n3,varu,varv,varp,vart,varr,topta)

use mem_scratch
use mem_grid
use rconstants
use therm_lib, only: virtt

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: varu,varv,varp,vart,varr
real, dimension(n2,n3) :: topta

integer :: i,j,k

! Interpolate from sigma-z varfile vertical coords to ADAP grid


do j=1,n3
   do i=1,n2
   
      do k=1,n1
         vctr10(k)=topta(i,j) + (1.-topta(i,j)/ztop)*ztn(k,ngrid)
      enddo
      vctr1(1:n1)=varu(1:n1,i,j)
      vctr2(1:n1)=varv(1:n1,i,j)
      vctr3(1:n1)=vart(1:n1,i,j)
      vctr4(1:n1)=varr(1:n1,i,j)
      call htint2(n1,vctr1,vctr10,n1,vctr11,ztn(:,ngrid))
      call htint2(n1,vctr2,vctr10,n1,vctr12,ztn(:,ngrid))
      call htint2(n1,vctr3,vctr10,n1,vctr13,ztn(:,ngrid))
      call htint2(n1,vctr4,vctr10,n1,vctr14,ztn(:,ngrid))

      ! Do hydrostatic balance
      do k=1,n1
         vctr15(k) = virtt(vctr13(k),vctr14(k))
      enddo

      vctr16(n1)= varp(n1,i,j) + g * (ztn(n1,ngrid) - vctr10(n1) )  &
               / vctr15(n1)
      do k = n1-1,1,-1
         vctr16(k) = vctr16(k+1) + g * (ztn(k+1,ngrid)-ztn(k,ngrid))  &
               /((vctr15(k)+vctr15(k+1))*.5)
      enddo
      
     varu(1:n1,i,j)= vctr11(1:n1)
     varv(1:n1,i,j)= vctr12(1:n1)
     vart(1:n1,i,j)= vctr13(1:n1)
     varr(1:n1,i,j)= vctr14(1:n1)
     varp(1:n1,i,j)= vctr16(1:n1)
      
   enddo      
enddo

return
end

!     **************************************************************

subroutine varref(n1,n2,n3,thp,pc,pi0,th0,rtp,dn0,dn0u,dn0v,uc  &
                 ,vc,topt,topu,topv,rtgt,rtgu,rtgv,topta,level)

use mem_grid
use ref_sounding
use mem_scratch
use rconstants
use therm_lib, only: virtt
                 
implicit none
integer :: n1,n2,n3,level          
real :: thp(n1,n2,n3),pc(n1,n2,n3),pi0(n1,n2,n3)  &
         ,rtp(n1,n2,n3),dn0(n1,n2,n3)  &
         ,dn0u(n1,n2,n3),dn0v(n1,n2,n3)  &
         ,uc(n1,n2,n3),vc(n1,n2,n3),topt(n2,n3),topu(n2,n3)  &
         ,topv(n2,n3),rtgt(n2,n3),rtgu(n2,n3),rtgv(n2,n3)  &
         ,th0(n1,n2,n3),topta(n2,n3)

integer :: i,j,k

!                Reference sounding is point with lowest topography
topref=1.e10
do j=1,nyp
   do i=1,nxp
      if(topta(i,j).lt.topref) then
         iref=i
         jref=j
         topref=topta(i,j)
      endif
   enddo
enddo

!  Set up 1-D reference state

if (if_adap == 0) then
   do k=1,nzp
      vctr2(k)=ztn(k,ngrid)*(1.-topref/ztop)+topref
   enddo
   call htint2(nzp,thp(:,iref,jref),vctr2,nzp,vctr1,zt)
   call htint2(nzp,uc (:,iref,jref),vctr2,nzp,u01dn(:,ngrid),zt)
   call htint2(nzp,vc (:,iref,jref),vctr2,nzp,v01dn(:,ngrid),zt)
   if (level >= 1) then
      call htint2(nzp,rtp(:,iref,jref),vctr2,nzp,rt01dn(:,ngrid),zt)
   else
      rt01dn(1:nzp,ngrid) = 0.
   endif
else
   vctr2(1:nzp)  =ztn(1:nzp,ngrid)
   vctr1(1:nzp)  =thp(1:nzp,iref,jref)
   u01dn(1:nzp,ngrid)=uc(1:nzp,iref,jref)
   v01dn(1:nzp,ngrid)=vc(1:nzp,iref,jref)
   rt01dn(1:nzp,ngrid) = 0.
   if (level >= 1) rt01dn(1:nzp,ngrid)=rtp(1:nzp,iref,jref)
endif


do k = 1,nzp
   th01dn(k,ngrid) = virtt(vctr1(k),rt01dn(k,ngrid))
end do
u01dn(1,ngrid) = u01dn(2,ngrid)
v01dn(1,ngrid) = v01dn(2,ngrid)
rt01dn(1,ngrid) = rt01dn(2,ngrid)
th01dn(1,ngrid) = th01dn(2,ngrid)

pi01dn(1,ngrid) = pc(1,iref,jref) + g * (vctr2(1) - zt(1))  &
   / (.5 * (th01dn(1,ngrid) + virtt(thp(1,iref,jref),rtp(1,iref,jref)) ))
do k = 2,nzp
  pi01dn(k,ngrid) = pi01dn(k-1,ngrid) - g / (dzm(k-1) * .5  &
     * (th01dn(k,ngrid) + th01dn(k-1,ngrid)))
enddo

do k = 1,nzp
  vctr4(k) = (pi01dn(k,ngrid) / cp) ** cpor * p00
  dn01dn(k,ngrid) = cp * vctr4(k)  &
     / (rgas * th01dn(k,ngrid) * pi01dn(k,ngrid))
enddo

!        Compute 3-D reference state from 1-D reference state

call refs3d(nzp,nxp,nyp,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

return

end

