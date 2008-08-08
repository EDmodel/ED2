!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
subroutine krig(n1,n2,n3,x,y,z,ndata,xd,yd,zd,ed,dvar,ngrid  &
  ,nnanzp,topt,var,varkrg,cazmod)
  
use mem_oda

implicit none

integer :: n1,n2,n3,ndata,ngrid,nnanzp
real :: z(n1),x(n2),y(n3),var(n1,n2,n3),varkrg(n1,n2,n3)
real, dimension(ndata) :: zd,xd,yd,ed,dvar
real, dimension(n2,n3) :: topt
real :: cazmod

integer, parameter :: maxobs=7000,maxobs2=7000
real, dimension(maxobs) :: xkp,ykp,zkp,dkp,c20,szg 
real(kind=8), dimension(maxobs) :: alph,c2 
real(kind=8), dimension(maxobs2) :: c1
character(len=80) :: label
integer :: p0,p1,p2,p3

real :: sx,sy,sz,szgm,r,dmean
real, external :: gam
integer :: i,j,k,nt,ntt,nxp,nyp,ii,n,m,nkp,ni,jk,ik,ktilt,l,ip,jp
!
!     Kriging optimal interpolation subroutine
!
!     Values at grid points are based on known data using Kriging.
!     II, IND(II), NT, and NTT are indices for storing the Kriging
!     configuration and weights to be used in Kriging the simulated
!     data.


!print*,'ndata:',ndata,nnanzp

dmean=0.0
ii=1

nt=0
ntt=0
nxp=n2
nyp=n3
do j=1,nyp
   do i=1,nxp
      do k=1,nnanzp
        ! ip=11; jp=8
        ! if(i==ip.and.j==jp)print*,'------------------'
        ! print*,'------------------',i,j,k
         var(k,i,j)=0.0
         varkrg(k,i,j)=0.0
         ii=ii+1
!
!     Calculate separation vector between grid point I,J,K and data
!     points, retain those with magnitude less than RMAXKRG
!
         n=0
         do m=1,ndata
           ! if(i==ip.and.j==jp.and.k==2) then
           !    print*,'###k1:',m,dvar(m),zd(m)
           ! endif
            if( dvar(m) < -998. ) cycle
            sx=x(i)-xd(m)
            sy=y(j)-yd(m)
            sz=z(k)-zd(m)
            szg(m)=abs(topt(i,j)-ed(m))
            szg(m)=0.
            szgm=szg(m)
            r=(sx**2+sy**2+(sz*cazkrg(3,ngrid)*cazmod)**2)**0.5
            szgm=0.0
            if(r.gt.rmaxkrg(k,ngrid) ) cycle
            n=n+1
            nt=nt+1
          ! if(i==ip.and.j==jp.and.k==2) then
           ! if(i==ip.and.j==jp) then
            !  print*,'rad:',k,n,m,r,rmaxkrg(k,ngrid)
            !   print*,'rad:',k,sx,sy
            !   print*,'rad:',k,sz,(sz*cazkrg(3,ngrid)*cazmod)**2
            !   print*,'rad:',n,m,dvar(m),sx,sy
            !endif
!
!     Calculate covariance between data point and grid point
!
            c2(n)=1.0-gam(sx,sy,sz,ngrid,k,szgm,cazmod)
            xkp(n)=xd(m)
            ykp(n)=yd(m)
            zkp(n)=zd(m)
            dkp(n)=dvar(m)
         enddo

         nkp=n+1
         if(nkp.eq.1) then
            var(k,i,j)=dmean
            varkrg(k,i,j)=2.0
            go to 901
         endif
!
!     Non-bias condition
!
         c2(n+1)=1.0
         do ni=1,nkp
            c20(ni)=c2(ni)
         enddo
!
!     Calculate data-data covariance matrix (N+1 * N+1) for data
!     configuration about grid point I,J,K
!
         p0=0
         do jk=1,n
            do ik=1,jk
               sx=xkp(jk)-xkp(ik)
               sy=ykp(jk)-ykp(ik)
               sz=zkp(jk)-zkp(ik)
               p0=p0+1
               szgm=abs(szg(jk)-szg(ik))
               szgm=0.0
               c1(p0)=1.0-gam(sx,sy,sz,ngrid,k,szgm,cazmod)
!if(i==ip.and.j==jp)print*,'c1def:',jk,ik,k,p0,c1(p0)
            enddo
         enddo
         p1=1+(n*(n+1)/2)
         p2=(n+1)*(n+2)/2-1
!if(i==ip.and.j==jp)print*,'p12:',p1,p2
!
!     Last column in matrix is all 1'S; non-bias condition
!
         do p0=p1,p2
            c1(p0)=1.0
         enddo
         p3=p2+1
         c1(p3)=0.0
!
!      Call matrix solver
!
         call relms(alph,c1,c2,nkp,1,ktilt,i,j)
         if(ktilt /= 0) then
         !!!print*,'tilt:',k,i,j,ktilt
            var(k,i,j)=dmean
            varkrg(k,i,j)=2.0
            cycle
         endif
            
!
!      Calculate Kriged grid point
!
         do l=1,n
            var(k,i,j)=var(k,i,j)+alph(l)*dkp(l)
    !!     if(i==ip.and.j==jp)print*,'vars:',l,var(k,i,j),alph(l),dkp(l),ktilt
            varkrg(k,i,j)=varkrg(k,i,j)+alph(l)*c20(l)
            ntt=ntt+1
         enddo
         varkrg(k,i,j)=1.0-alph(nkp)-varkrg(k,i,j)
     !    if(i==ip.and.j==jp)print*,'vars:',var(k,i,j),varkrg(k,i,j)
         901   continue
      enddo
   enddo
enddo

return
end
!
!     ****************************************************************
!
real function gam(hx,hy,hz,ngrid,k,hzg,cazmod)

use mem_oda

implicit none

real :: hx,hy,hz,hzg,cazmod
integer :: ngrid, k


real :: h,dxmod,dymod,dx,dy,dz
integer :: is,ijs,ijs1,ijs2

!
!     Variogram calculation
!
gam=0.0
h=(hx*hx+hy*hy+hz*hz)**0.5
if(h .lt. 1.e-03) return
do 1 is=1,nstkrg(ngrid)
  ijs=is
  ijs1=ijs+nstkrg(ngrid)
  ijs2=ijs1+nstkrg(ngrid)
  dxmod=hzg/750.0*(-akrg(1,ngrid))
  dymod=dxmod
  dx=hx*caxkrg(ijs,ngrid)+hy*caxkrg(ijs1,ngrid)+  &
    hz*caxkrg(ijs2,ngrid)
  dx=abs(dx)+dxmod
  dy=hx*caykrg(ijs,ngrid)+hy*caykrg(ijs1,ngrid)+  &
    hz*caykrg(ijs2,ngrid)
  dy=abs(dy)+dymod
  dz=hx*cazkrg(ijs,ngrid)+hy*cazkrg(ijs1,ngrid)+  &
    hz*cazkrg(ijs2,ngrid)*cazmod
  h=(dx*dx+dy*dy+dz*dz)**0.5
  if(akrg(k,ngrid))10,12,11
10   gam=gam+ckrg(is,ngrid)*(1.-exp(h/akrg(k,ngrid)))
  go to 1
11   if(h .ge. akrg(k,ngrid))go to 12
  gam=gam+ckrg(is,ngrid)*(1.5*h/akrg(k,ngrid)-  &
  0.5*h*h*h/(akrg(k,ngrid)*akrg(k,ngrid)*akrg(k,ngrid)))
  go to 1
12   gam=gam+ckrg(is,ngrid)
1 continue
!
return
end
!
!     ****************************************************************
!
subroutine relms(x,a,b,m,n,ktilt,ip,jp)
implicit none
integer :: m,n,ktilt,ip,jp
!
!     Matrix solver for Kriging equations
!
real(kind=8) :: x(*),a(*),b(*)
real(kind=8) :: r,piv,tol
real :: ak
integer :: i,nm,m1,kk,k,lp,km1,ll,ij,ii,j,llb,in
!
!if(ip==4.and.jp==11) then
!print*,'relms0:',m,n
!do k=1,m
!print*,'a:',k,a(k),b(1)
!enddo
!endif


if(m .le. 0) go to 14
if(m .gt. 1) go to 1
if(a(1) .eq. 0.) go to 20
do 19 i=1,n
  x(i)=b(i)/a(1)
19 continue
ktilt=0
return
20 ktilt=1
return
!
!     Initialize
!
1 tol=0.0
ktilt=0
nm=n*m
m1=m-1
kk=0
!
!     Start triangulation
!
do 7 k=1,m1
  kk=kk+k
  ak=a(kk)
!if(ip==4.and.jp==11) print*,'relms1:',k,kk,ak,tol
  if(ak-tol)3,2,3
2   ktilt=k
!if(ip==4.and.jp==11) print*,'return1:',k,kk,ak,tol
  return
3   piv=1./ak
  ii=kk
  lp=0
  km1=k-1
  do 6 i=k,m1
    ll=ii
    ii=ii+i
    r=a(ii)*piv
    lp=lp+1
    ij=ii-km1
    do 4 j=i,m1
      ij=ij+j
      ll=ll+j
!if(ip==4.and.jp==11) print*,'def_a:',j,ij,ll,a(ij),r*a(ll)
      a(ij)=a(ij)-r*a(ll)
4     continue
    do 5 llb=k,nm,m
      in=llb+lp
      b(in)=b(in)-r*b(llb)
5     continue
6   continue
7 continue
r=a(ij)
if(r-tol)9,8,9
8 ktilt=m
return
!
!     End triangulation
!
!     Start back solution
!
9 piv=1./r
do 10 llb=m,nm,m
  x(llb)=b(llb)*piv
10 continue
i=m
kk=ij
do 13 ii=1,m1
  kk=kk-i
  piv=1./a(kk)
  i=i-1
  do 12 llb=i,nm,m
    in=llb
    r=b(in)
    ij=kk
    do 11 j=i,m1
      ij=ij+j
      in=in+1
      r=r-a(ij)*x(in)
11     continue
    x(llb)=r*piv
12   continue
13 continue
!
!     End solution
!
return
!
!     Error return
!
14 ktilt=-1
!
return
end
