!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine nest_interpolated_topo(nxc,nyc,nxf,nyf,maxx,maxy,ifm  &
                                 ,topo_c,topo_i,scr1,scr2)

implicit none

integer :: nxc,nyc,nxf,nyf,maxx,maxy,ifm
real, dimension(*) :: topo_c,topo_i,scr1,scr2

!    Fill topo_i with interpolated topography from coarser grid

call fillscr(1,maxx,maxy,1,nxc,nyc,1,1,scr1(1),topo_c(1))
call eintp(scr1(1),scr2(1),1,maxx,maxy,1,nxf,nyf,ifm,2,'t',0,0)
call fillvar(1,maxx,maxy,1,nxf,nyf,1,1,scr2(1),topo_i(1))

return
end


subroutine fmrefs1d_isan(ifm,icm,n0,n1  &
                   ,piref,thref,dnref,rtref)

use rconstants

implicit none

integer :: ifm,icm,n1,n0
real, dimension(n0,*) :: piref,thref,dnref,rtref

integer :: k
real :: c1,c2
real, allocatable :: vctr1(:),vctr2(:),vctr3(:),vctr4(:)

!     Interpolate the fine mesh 1-d reference state variables.

allocate(vctr1(n1),vctr2(n1),vctr3(n1),vctr4(n1))

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
if (icm >= 1) then
   do k = 1,n1
      vctr1(k) = thref(k,icm) * dnref(k,icm)
      vctr2(k) = rtref(k,icm) * dnref(k,icm)
   enddo

  call eintp(dnref(1,icm),dnref(1,ifm),n1,1,1,n1  &
     ,1,1,ifm,1,'t',0,0)
  call eintp(vctr1,vctr3,n1,1,1,n1,1,1,ifm,1,'t',0,0)
  call eintp(vctr2,vctr4,n1,1,1,n1,1,1,ifm,1,'t',0,0)

   do k = 1,n1
      thref(k,ifm) = vctr3(k) / dnref(k,ifm)
      rtref(k,ifm) = vctr4(k) / dnref(k,ifm)
      piref(k,ifm) = c2 * (dnref(k,ifm) * thref(k,ifm)) ** c1
   enddo
endif

deallocate(vctr1,vctr2,vctr3,vctr4)

return
end

!     *****************************************************************

subroutine fmrefs3d_isan (ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy,nbot,ntop,jd  &
          ,scr1,scr2,vt2da,toptf,toptc,dn0c,dn0f,th0c,th0f  &
          ,pi0f,dn0uf,dn0vf,zt,ztop)

use rconstants

implicit none
integer :: ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy,nbot,ntop,jd
real, dimension(*) :: scr1,scr2,vt2da,zt
real, dimension(n2f,n3f) :: toptf  
real, dimension(n2c,n3c) :: toptc     
real, dimension(n1f,n2f,n3f) :: dn0f,th0f,pi0f,dn0uf,dn0vf
real, dimension(n1c,n2c,n3c) :: dn0c,th0c   
real :: ztop  

real :: c1,c2,b(1)
integer :: i1,j1,i,j,k

!     Interpolate the fine mesh 3-D reference state variables.

if (icm == 0) return

call fmint3(n1c,n2c,n3c,n1f,n2f,n3f,maxiz,maxix,maxiy  &
     ,ifm,icm,nbot,ntop,jd,1,0,0,'t'  &
     ,dn0c,dn0f,dn0c,dn0f,scr1,scr2,toptf,vt2da,b(1),b(1),b(1),0)

call fmint3(n1c,n2c,n3c,n1f,n2f,n3f,maxiz,maxix,maxiy  &
     ,ifm,icm,nbot,ntop,jd,1,0,0,'t'  &
     ,th0c,th0f,dn0c,dn0f,scr1,scr2,toptf,vt2da,b(1),b(1),b(1),0)

c1 = rgas / (cp - rgas)
c2 = cp * (rgas / p00) ** c1
pi0f(1:n1f,1:n2f,1:n3f) = c2 * (dn0f(1:n1f,1:n2f,1:n3f)  &
                            *   th0f(1:n1f,1:n2f,1:n3f) ) ** c1


call fillscr(1,maxix,maxiy,1,n2c,n3c,1,1,scr1,toptc)
call eintp(scr1,scr2,1,maxix,maxiy,1,n2f,n3f,ifm,2,'t',0,0)
call fillvar(1,maxix,maxiy,1,n2f,n3f,1,1,scr2,scr1)

call rtgintrp_isan(n1f,n2f,n3f,th0f,scr1,toptf,zt,ztop)
call rtgintrp_isan(n1f,n2f,n3f,pi0f,scr1,toptf,zt,ztop)

! Define dn0u and dn0v

do j = 1,n3f
   j1 = min(j+1,n3f)
   do i = 1,n2f
      i1 = min(i+1,n2f)
      do k = 1,n1f
         dn0uf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i1,j))
         dn0vf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i,j1))
      enddo
   enddo
enddo

return
end

!*****************************************************************************

subroutine fmdn0_isan(ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy,scr1,scr2,toptf,toptc,dn0f,dn0uf,dn0vf,zt,ztop)

implicit none

integer :: ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c  &
          ,maxiz,maxix,maxiy
real, dimension(*) :: scr1,scr2,zt
real, dimension(n2f,n3f) :: toptf 
real, dimension(n2c,n3c) :: toptc     
real, dimension(n1f,n2f,n3f) :: dn0f,dn0uf,dn0vf 
real :: ztop

integer :: i,j,i1,j1,k
!     Special vertical interpolation of DN0 must be done after all other
!     3-D reference state and prognostic variables are interpolated.

if (icm == 0) return

call fillscr(1,maxix,maxiy,1,n2c,n3c,1,1,scr1,toptc)
call eintp(scr1,scr2,1,maxix,maxiy,1,n2f,n3f,ifm,2,'t',0,0)
call fillvar(1,maxix,maxiy,1,n2f,n3f,1,1,scr2,scr1)

call rtgintrp_isan(n1f,n2f,n3f,dn0f,scr1,toptf,zt,ztop)

! Define dn0u and dn0v

do j = 1,n3f
   j1 = min(j+1,n3f)
   do i = 1,n2f
      i1 = min(i+1,n2f)
      do k = 1,n1f
         dn0uf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i1,j))
         dn0vf(k,i,j) = .5 * (dn0f(k,i,j) + dn0f(k,i,j1))
      enddo
   enddo
enddo

return
end

!******************************************************************************

subroutine rtgintrp_isan(n1,n2,n3,fld,vt2da,topt,zt,ztop)

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: fld
real, dimension(n2,n3) :: vt2da,topt
real, dimension(n1) :: zt
real :: ztop

integer :: i,j,k
real, allocatable, dimension(:) :: vctr1,vctr2,vctr3

!     Do special vertical interpolation in case terrain on this grid
!     (topt) is different from what would be interpolated from the
!     coarser grid (vt2da).

allocate(vctr1(n1),vctr2(n1),vctr3(n1))

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
         vctr1(k) = zt(k) * (1. - vt2da(i,j) / ztop) + vt2da(i,j)
         vctr2(k) = zt(k) * (1. - topt(i,j) / ztop) + topt(i,j)
         vctr3(k) = fld(k,i,j)
      enddo
      call htint(n1,vctr3,vctr1,n1,fld(1,i,j),vctr2)
   enddo
enddo

deallocate(vctr1,vctr2,vctr3)

return
end
!     ****************************************************************

subroutine varfile_refstate(n1,n2,n3,thp,pc,pi0,th0,rtp,dn0  &
                 ,dn0u,dn0v,topt,rtgt,zt,ztop,piref,thref,dnref,rtref)

use rconstants
                 
implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: thp,pc,pi0,rtp,dn0,th0,dn0u,dn0v
real, dimension(n2,n3)    :: topt,topu,topv,rtgt,rtgu,rtgv
real, dimension(n1)       :: zt
real, dimension(n1)        :: piref,thref,rtref,dnref
real :: ztop

real :: topref,c1,c2,c3
real, allocatable :: vctr1(:),vctr2(:)
integer :: ir,jr,i,j,k,i1,j1

!                Reference sounding is point with lowest topography

allocate(vctr1(n1),vctr2(n1))

topref=1.e10
do j=1,n3
   do i=1,n2
      if(topt(i,j).lt.topref) then
         ir=i
         jr=j
         topref=topt(i,j)
      endif
   enddo
enddo

do k=1,n1
   vctr2(k)=zt(k)*(1.-topref/ztop)+topref
enddo

call htint2(n1,thp(1,ir,jr),vctr2,n1,vctr1,zt)
call htint2(n1,rtp(1,ir,jr),vctr2,n1,rtref(1),zt)

do k = 1,n1
   thref(k) = vctr1(k) * (1. + .61 * rtref(k))
enddo
rtref(1) = rtref(2)
thref(1) = thref(2)

piref(1) = pc(1,ir,jr) + g * (vctr2(1) - zt(1))  &
                / (.5 * (thref(1)  &
                + thp(1,ir,jr) * (1. + .61 * rtp(1,ir,jr))))
do k = 2,n1
   piref(k) = piref(k-1) - g * (zt(k)-zt(k-1))  &
             / ( .5 * (thref(k) + thref(k-1)))
enddo

do k = 1,n1
  vctr1(k) = (piref(k) / cp) ** cpor * p00
  dnref(k) = cp * vctr1(k)  &
     / (rgas * thref(k) * piref(k))
enddo

!        Compute 3-D reference state from 1-D reference state

do j=1,n3
  do i=1,n2

    do k=1,n1
      vctr2(k)=zt(k)*rtgt(i,j)+topt(i,j)
    enddo
    call htint(n1,piref(1),zt,n1,pi0(1,i,j),vctr2)
    call htint(n1,thref(1),zt,n1,th0(1,i,j),vctr2)

    c1=g*2.*(1.-topt(i,j)/ztop)
    c2=(1-cpor)
    c3=cp**c2
    do k=n1-1,1,-1
      pi0(k,i,j)=pi0(k+1,i,j)  &
                +c1*(zt(k+1)-zt(k))/(th0(k,i,j)+th0(k+1,i,j))
    enddo

    do k=1,n1
      dn0(k,i,j)=(c3*p00)/(rgas*th0(k,i,j)*pi0(k,i,j)**c2)
    enddo

  enddo
enddo

do j = 1,n3
   j1 = min(j+1,n3)
   do i = 1,n2
      i1 = min(i+1,n2)
      do k = 1,n1
         dn0u(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i1,j))
         dn0v(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i,j1))
      enddo
   enddo
enddo

deallocate(vctr1,vctr2)

return
end
