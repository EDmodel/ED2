!############################# change log ##################################
! 4.3.0.2
!
! 010328 jhc gdtost2 ##
!            added second horizontal interpolation scheme which is
!            a simple bilinear one to eliminate ringing errors. ##
! 001012 mjb htint ##
!            declared passed args with length * instead of 1 to pass array
!            bounds compiler option for revu. ##
!
!###########################################################################
!  copyright (c)  1990, 1995, 1999, 2000 - all rights reserved
!  regional atmospheric modeling system - rams
!  mission research corporation / *aster division
!###########################################################################

subroutine trncl1(vctra,zlev,zsurf,htop,vctrb,vctrc,vctrd,nn)
dimension vctra(1),vctrb(1),vctrc(1),vctrd(1),zlev(1)

!     transform to z coordinates form zstar coordinates.

if(zsurf.eq.0) return
rtg=1.-zsurf/htop
kabv=2
do k=2,nn
vctrc(k)=1./(rtg*(zlev(k)-zlev(k-1)))
vctrd(k)=zlev(k)*rtg+zsurf
enddo
vctrd(1)=zsurf
do k=1,nn
kk1=kabv
do kk=kk1,nn
kabv=kk
if(vctrd(kabv).ge.zlev(k)-1)go to 10
enddo
write(6,5) id,jd
5     format(' stop in trncl1',2i5)
stop
10    continue
w1=min(1.,(vctrd(kabv)-zlev(k))*vctrc(kabv))
vctrb(k)=vctra(kabv-1)*w1+vctra(kabv)*(1.-w1)
enddo
do k=1,nn
vctra(k)=vctrb(k)
enddo
return
end

!     ******************************************************************

subroutine trncl2(vctra,zlev,zsurf,htop,vctrb,vctrc,vctrd,nn)
dimension vctra(1),vctrb(1),vctrc(1),vctrd(1),zlev(1)

!     transform to zstar coordinates from z coordinates

nnp=nn+1
if(zsurf.eq.0) return
rtg=1.-zsurf/htop
kabv=2
do k=2,nnp
vctrc(k)=zlev(k)*rtg+zsurf
vctrd(k)=1./(zlev(k)-zlev(k-1))
enddo
vctrc(1)=zsurf
do k=1,nn
kk1=kabv
do kk=kk1,nnp
kabv=kk
if(zlev(kabv)+.01.ge.vctrc(k)) go to 10
enddo
write(6,5) id,jd
5     format(' stop in trncl2',2i5)
stop
10    continue
w1=min((zlev(kabv)-vctrc(k))*vctrd(kabv),1.)
vctrb(k)=w1*vctra(kabv-1)+(1.-w1)*vctra(kabv)
enddo
do k=1,nn
vctra(k)=vctrb(k)
enddo
return
end

!     ******************************************************************

subroutine intrp(a,b,zval,z,nzp)

!     linear interpolation with z

dimension a(1),z(1)
do k=1,nzp
kabv=k
if(zval.lt.z(k)) go to 10
enddo
10    continue
if(kabv.eq.1) kabv=2
b=(a(kabv)*(zval-z(kabv-1))+a(kabv-1)*(z(kabv)-zval))  &
/(z(kabv)-z(kabv-1))
return
end

!     ******************************************************************

subroutine intrrap(a,b,plnval,pln,nzp)

!     linear interpolation with log pressure

dimension pln(1),a(1)
nz=nzp-1
do kk=1,nz
k=nzp-kk
kabv=k+1
if(plnval.lt.pln(k)) go to 10
enddo
10    continue
b=(a(kabv)*(plnval-pln(kabv-1))+a(kabv-1)*(pln(kabv)-plnval))  &
/(pln(kabv)-pln(kabv-1))
return
end

!     ******************************************************************

subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
use misc_coms, only : itypp,i0x,i1x,i2x,yoo
 yyy=1e30
 if(x2.gt.1.e19.or.x3.gt.1.e19.or.  &
   y2.gt.1.e19.or.y3.gt.1.e19)return
wt1=(xxx-x3)/(x2-x3)
wt2=1.0-wt1
istend=0
if(y4.lt.1.e19.and.x4.lt.1.e19) go to 410
yz22=wt1
yz23=wt2
yz24=0.0
istend= 1
410   if(y1.lt.1.e19.and.x1.lt.1.e19) go to 430
yz11=0.0
yz12=wt1
yz13=wt2
if(istend.eq.1)go to 480
go to 450
430   yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
if(istend.eq.  1    ) go to 470
450   yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
470   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
 go to 490
480      yyy=wt1*y2+wt2*y3
490   yoo=yyy
return
end

!     ******************************************************************

subroutine gdtost(a,ix,iy,stax,stay,staval)

!     subroutine to return stations back-interpolated values(staval)
!     from uniform grid points using overlapping-quadratics.
!     gridded values of input array a dimensioned a(ix,iy),where
!     ix=grid points in x, iy = grid points in y .  station
!     location given in terms of grid relative station x (stax)
!     and station column.
!     values greater than 1.0e30 indicate missing data.

dimension a(ix,iy),r(4),scr(4)
iy1=int(stay)-1
iy2=iy1+3
ix1=int(stax)-1
ix2=ix1+3
staval=1e30
fiym2=float(iy1)-1
fixm2=float(ix1)-1
ii=0
do 100 i=ix1,ix2
ii=ii+1
if(i.ge.1.and.i.le.ix) go to 101
scr(ii)=1e30
go to 100
101   jj=0
do 111 j=iy1,iy2
jj=jj+1
if(j.ge.1.and.j.le.iy) go to 112
r(jj)=1e30
go to 111
112   r(jj)=a(i,j)
111   continue
yy=stay-fiym2
call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
100   continue
xx=stax-fixm2
call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
return
end

!     ******************************************************************

subroutine gdtost2(a,ix,iy,stax,stay,staval)

!     subroutine to return stations back-interpolated values (staval)
!     from uniform grid points using bi-linear interpolation.
!     gridded values of input array a dimensioned a(ix,iy), where
!     ix=grid points in x, iy = grid points in y.  station
!     location given in terms of grid relative station x,y (stax,stay).

implicit none

! passed variables

integer :: ix, iy
real :: stax,stay,staval
real, dimension(ix,iy) :: a

! internal variables

integer :: i,j
real :: wtx1,wtx2,wty1,wty2


i = int(stax)
j = int(stay)
wtx2 = stax - float(i)
wty2 = stay - float(j)
wtx1 = 1. - wtx2
wty1 = 1. - wty2

staval = wtx1 * (wty1 * a(i  ,j  )+ wty2 * a(i  ,j+1)) &
       + wtx2 * (wty1 * a(i+1,j  )+ wty2 * a(i+1,j+1))
return
end subroutine gdtost2


!     ******************************************************************

subroutine htintcp(nzz1,vctra,eleva,nzz2,vctrb,elevb,  &
 vt1c,vt2c,vt1d,vt2d,vt1e,vt2e)
dimension vctra(1),vctrb(1),eleva(1),elevb(1)
dimension vt1c(1),vt1d(1),vt1e(1),vt2c(1),vt2d(1),vt2e(1)

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 35
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
l=l+1
if(l.eq.nzz1)stop 'htintcp'
go to 30
35 continue
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
vt2c(k)=vt1c(l)+(vt1c(l+1)-vt1c(l))*wt
vt2d(k)=vt1d(l)+(vt1d(l+1)-vt1d(l))*wt
vt2e(k)=vt1e(l)+(vt1e(l+1)-vt1e(l))*wt
20 continue

return
end

!     ******************************************************************

subroutine htint(nzz1,vctra,eleva,nzz2,vctrb,elevb)
dimension vctra(*),vctrb(*),eleva(*),elevb(*)

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 35
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
if(elevb(k).gt.eleva(nzz1))go to 36
l=l+1
if(l.eq.nzz1) then
  print *,'htint:nzz1',nzz1
  do kk=1,l
    print*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
  enddo
  stop 'htint'
endif
go to 30
35 continue
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
go to 20
36 continue
wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
20 continue

return
end

!     **************************************************************

subroutine htint2(nzz1,vctra,eleva,nzz2,vctrb,elevb)
dimension vctra(1),vctrb(1),eleva(1),elevb(1)

!      htint for holding values of vctrb constant under eleva(1)

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 34
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
if(elevb(k).gt.eleva(nzz1))go to 36
l=l+1
if(l.eq.nzz1)stop 'htint2'
go to 30
34   continue
vctrb(k)=vctra(1)
go to 20
35 continue
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
go to 20
36 continue
wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
20 continue

return
end

!     **************************************************************

subroutine awtcmp(x,ii,xx,ij,ityp,llb,lrb,icon,w,iord)
dimension w(6,6),x(1),xx(1),q(6)

!     ityp=0, normal, llb,lrb-left,right boundaries
!          1, isym  , llb-symmetric point,lrb-right boundary

do l=1,6
   do ll=1,6
      w(l,ll)=0.
   enddo
enddo

q12=xx(ij)

if(ityp.eq.0.or.(ityp.eq.1.and.ii.ge.llb+2)) then
  if(ii-2.lt.llb.or.ii+3.gt.lrb.or.iord.eq.2) then
    w(1,3)=(x(ii+1)-xx(ij))/(x(ii+1)-x(ii))
    w(2,3)=.5/(x(ii+1)-x(ii))
    w(1,4)=-(x(ii)-xx(ij))/(x(ii+1)-x(ii))
    w(2,4)=-.5/(x(ii+1)-x(ii))
    return
  endif
  goto 3000
endif

do l=1,6
   q(l)=x(ii-3+l)
enddo

if(ityp.eq.1.and.ii.lt.llb+2) then
   q(3)=x(ii)
   q(4)=x(ii+1)
   q(5)=x(ii+2)
   if(ii.eq.llb) then
      q(2)=2.*x(ii)-x(ii+1)
      q(1)=2.*x(ii)-x(ii+2)
   elseif(ii.eq.llb+1) then
      q(2)=x(ii-1)
      q(1)=2.*x(ii)-x(ii+1)
   endif
endif

3000 continue
goto (1000,2000) icon+1
1000 continue

w(1,1)=  &
    ((q(6)-q12)*(q(5)-q12)*(q(4)-q12)*(q(3)-q12)*(q(2)-q12  &
 ))/((q(6)-q(1))*(q(5)-q(1))*(q(4)-q(1))*(q(3)-q(1))*(q(2)  &
 -q(1)))
w(2,1)=  &
    ((((q(3)+q(2)-2.*q12)*q(4)+(q(2)-2.*q12)*q(3)-2.*q(2)*  &
 q12+3.*q12**2)*q(5)+((q(2)-2.*q12)*q(3)-2.*q(2)*q12+3.*  &
 q12**2)*q(4)-(2.*q(2)-3.*q12)*q(3)*q12+3.*q(2)*q12**2-4.*  &
 q12**3)*q(6)+(((q(2)-2.*q12)*q(3)-2.*q(2)*q12+3.*q12**2)*  &
 q(4)-(2.*q(2)-3.*q12)*q(3)*q12+3.*q(2)*q12**2-4.*q12**3)*  &
 q(5)-((2.*q(2)-3.*q12)*q(3)-3.*q(2)*q12+4.*q12**2)*q(4)*  &
 q12+(3.*q(2)-4.*q12)*q(3)*q12**2-4.*q(2)*q12**3+5.*q12**4  &
 )/(2.*(q(6)-q(1))*(q(5)-q(1))*(q(4)-q(1))*(q(3)-q(1))*(q(  &
 2)-q(1)))
w(3,1)=  &
    (((q(4)+q(3)+q(2)-3.*q12)*q(5)+(q(3)+q(2)-3.*q12)*q(4)  &
 +(q(2)-3.*q12)*q(3)-3.*q(2)*q12+6.*q12**2)*q(6)+((q(3)+q(  &
 2)-3.*q12)*q(4)+(q(2)-3.*q12)*q(3)-3.*q(2)*q12+6.*q12**2)  &
 *q(5)+((q(2)-3.*q12)*q(3)-3.*q(2)*q12+6.*q12**2)*q(4)-3.*  &
 (q(2)-2.*q12)*q(3)*q12+6.*q(2)*q12**2-10.*q12**3)/(3.*(q(  &
 6)-q(1))*(q(5)-q(1))*(q(4)-q(1))*(q(3)-q(1))*(q(2)-q(1)))
w(4,1)=  &
    ((q(5)+q(4)+q(3)+q(2)-4.*q12)*q(6)+(q(4)+q(3)+q(2)-4.*  &
 q12)*q(5)+(q(3)+q(2)-4.*q12)*q(4)+(q(2)-4.*q12)*q(3)-4.*q  &
 (2)*q12+10.*q12**2)/(4.*(q(6)-q(1))*(q(5)-q(1))*(q(4)-q(1  &
 ))*(q(3)-q(1))*(q(2)-q(1)))
w(5,1)=  &
    (q(6)+q(5)+q(4)+q(3)+q(2)-5.*q12)/(5.*(q(6)-q(1))*(q(5  &
 )-q(1))*(q(4)-q(1))*(q(3)-q(1))*(q(2)-q(1)))
w(6,1)=  &
    1./(6.*(q(6)-q(1))*(q(5)-q(1))*(q(4)-q(1))*(q(3)-q(1))  &
 *(q(2)-q(1)))
w(1,2)=  &
    (-(q(6)-q12)*(q(5)-q12)*(q(4)-q12)*(q(3)-q12)*(q(1)-  &
 q12))/((q(6)-q(2))*(q(5)-q(2))*(q(4)-q(2))*(q(3)-q(2))*(q  &
 (2)-q(1)))
w(2,2)=  &
    (-((((q(3)+q(1)-2.*q12)*q(4)+(q(1)-2.*q12)*q(3)-2.*q(1  &
 )*q12+3.*q12**2)*q(5)+((q(1)-2.*q12)*q(3)-2.*q(1)*q12+3.*  &
 q12**2)*q(4)-(2.*q(1)-3.*q12)*q(3)*q12+3.*q(1)*q12**2-4.*  &
 q12**3)*q(6)+(((q(1)-2.*q12)*q(3)-2.*q(1)*q12+3.*q12**2)*  &
 q(4)-(2.*q(1)-3.*q12)*q(3)*q12+3.*q(1)*q12**2-4.*q12**3)*  &
 q(5)-((2.*q(1)-3.*q12)*q(3)-3.*q(1)*q12+4.*q12**2)*q(4)*  &
 q12+(3.*q(1)-4.*q12)*q(3)*q12**2-4.*q(1)*q12**3+5.*q12**4  &
 ))/(2.*(q(6)-q(2))*(q(5)-q(2))*(q(4)-q(2))*(q(3)-q(2))*(q  &
 (2)-q(1)))
w(3,2)=  &
    (-(((q(4)+q(3)+q(1)-3.*q12)*q(5)+(q(3)+q(1)-3.*q12)*q(  &
 4)+(q(1)-3.*q12)*q(3)-3.*q(1)*q12+6.*q12**2)*q(6)+((q(3)+  &
 q(1)-3.*q12)*q(4)+(q(1)-3.*q12)*q(3)-3.*q(1)*q12+6.*q12**  &
 2)*q(5)+((q(1)-3.*q12)*q(3)-3.*q(1)*q12+6.*q12**2)*q(4)-  &
 3.*(q(1)-2.*q12)*q(3)*q12+6.*q(1)*q12**2-10.*q12**3))/(3.*  &
 (q(6)-q(2))*(q(5)-q(2))*(q(4)-q(2))*(q(3)-q(2))*(q(2)-q(1  &
 )))
w(4,2)=  &
    (-((q(5)+q(4)+q(3)+q(1)-4.*q12)*q(6)+(q(4)+q(3)+q(1)-  &
 4.*q12)*q(5)+(q(3)+q(1)-4.*q12)*q(4)+(q(1)-4.*q12)*q(3)-4.  &
 *q(1)*q12+10.*q12**2))/(4.*(q(6)-q(2))*(q(5)-q(2))*(q(4)-  &
 q(2))*(q(3)-q(2))*(q(2)-q(1)))
w(5,2)=  &
    (-(q(6)+q(5)+q(4)+q(3)+q(1)-5.*q12))/(5.*(q(6)-q(2))*(  &
 q(5)-q(2))*(q(4)-q(2))*(q(3)-q(2))*(q(2)-q(1)))
w(6,2)=  &
    (-1.)/(6.*(q(6)-q(2))*(q(5)-q(2))*(q(4)-q(2))*(q(3)-q(  &
 2))*(q(2)-q(1)))
w(1,3)=  &
    ((q(6)-q12)*(q(5)-q12)*(q(4)-q12)*(q(2)-q12)*(q(1)-q12  &
 ))/((q(6)-q(3))*(q(5)-q(3))*(q(4)-q(3))*(q(3)-q(2))*(q(3)  &
 -q(1)))
w(2,3)=  &
    ((((q(2)+q(1)-2.*q12)*q(4)+(q(1)-2.*q12)*q(2)-2.*q(1)*  &
 q12+3.*q12**2)*q(5)+((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*  &
 q12**2)*q(4)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*  &
 q12**3)*q(6)+(((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*q12**2)*  &
 q(4)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*q12**3)*  &
 q(5)-((2.*q(1)-3.*q12)*q(2)-3.*q(1)*q12+4.*q12**2)*q(4)*  &
 q12+(3.*q(1)-4.*q12)*q(2)*q12**2-4.*q(1)*q12**3+5.*q12**4  &
 )/(2.*(q(6)-q(3))*(q(5)-q(3))*(q(4)-q(3))*(q(3)-q(2))*(q(  &
 3)-q(1)))
w(3,3)=  &
    (((q(4)+q(2)+q(1)-3.*q12)*q(5)+(q(2)+q(1)-3.*q12)*q(4)  &
 +(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(6)+((q(2)+q(  &
 1)-3.*q12)*q(4)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)  &
 *q(5)+((q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(4)-3.*  &
 (q(1)-2.*q12)*q(2)*q12+6.*q(1)*q12**2-10.*q12**3)/(3.*(q(  &
 6)-q(3))*(q(5)-q(3))*(q(4)-q(3))*(q(3)-q(2))*(q(3)-q(1)))
w(4,3)=  &
    ((q(5)+q(4)+q(2)+q(1)-4.*q12)*q(6)+(q(4)+q(2)+q(1)-4.*  &
 q12)*q(5)+(q(2)+q(1)-4.*q12)*q(4)+(q(1)-4.*q12)*q(2)-4.*q  &
 (1)*q12+10.*q12**2)/(4.*(q(6)-q(3))*(q(5)-q(3))*(q(4)-q(3  &
 ))*(q(3)-q(2))*(q(3)-q(1)))
w(5,3)=  &
    (q(6)+q(5)+q(4)+q(2)+q(1)-5.*q12)/(5.*(q(6)-q(3))*(q(5  &
 )-q(3))*(q(4)-q(3))*(q(3)-q(2))*(q(3)-q(1)))
w(6,3)=  &
    1./(6.*(q(6)-q(3))*(q(5)-q(3))*(q(4)-q(3))*(q(3)-q(2))  &
 *(q(3)-q(1)))
w(1,4)=  &
    (-(q(6)-q12)*(q(5)-q12)*(q(3)-q12)*(q(2)-q12)*(q(1)-  &
 q12))/((q(6)-q(4))*(q(5)-q(4))*(q(4)-q(3))*(q(4)-q(2))*(q  &
 (4)-q(1)))
w(2,4)=  &
    (-((((q(2)+q(1)-2.*q12)*q(3)+(q(1)-2.*q12)*q(2)-2.*q(1  &
 )*q12+3.*q12**2)*q(5)+((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*  &
 q12**2)*q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*  &
 q12**3)*q(6)+(((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*q12**2)*  &
 q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*q12**3)*  &
 q(5)-((2.*q(1)-3.*q12)*q(2)-3.*q(1)*q12+4.*q12**2)*q(3)*  &
 q12+(3.*q(1)-4.*q12)*q(2)*q12**2-4.*q(1)*q12**3+5.*q12**4  &
 ))/(2.*(q(6)-q(4))*(q(5)-q(4))*(q(4)-q(3))*(q(4)-q(2))*(q  &
 (4)-q(1)))
w(3,4)=  &
    (-(((q(3)+q(2)+q(1)-3.*q12)*q(5)+(q(2)+q(1)-3.*q12)*q(  &
 3)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(6)+((q(2)+  &
 q(1)-3.*q12)*q(3)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**  &
 2)*q(5)+((q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(3)-  &
 3.*(q(1)-2.*q12)*q(2)*q12+6.*q(1)*q12**2-10.*q12**3))/(3.*  &
 (q(6)-q(4))*(q(5)-q(4))*(q(4)-q(3))*(q(4)-q(2))*(q(4)-q(1  &
 )))
w(4,4)=  &
    (-((q(5)+q(3)+q(2)+q(1)-4.*q12)*q(6)+(q(3)+q(2)+q(1)-  &
 4.*q12)*q(5)+(q(2)+q(1)-4.*q12)*q(3)+(q(1)-4.*q12)*q(2)-4.  &
 *q(1)*q12+10.*q12**2))/(4.*(q(6)-q(4))*(q(5)-q(4))*(q(4)-  &
 q(3))*(q(4)-q(2))*(q(4)-q(1)))
w(5,4)=  &
    (-(q(6)+q(5)+q(3)+q(2)+q(1)-5.*q12))/(5.*(q(6)-q(4))*(  &
 q(5)-q(4))*(q(4)-q(3))*(q(4)-q(2))*(q(4)-q(1)))
w(6,4)=  &
    (-1.)/(6.*(q(6)-q(4))*(q(5)-q(4))*(q(4)-q(3))*(q(4)-q(  &
 2))*(q(4)-q(1)))
w(1,5)=  &
    ((q(6)-q12)*(q(4)-q12)*(q(3)-q12)*(q(2)-q12)*(q(1)-q12  &
 ))/((q(6)-q(5))*(q(5)-q(4))*(q(5)-q(3))*(q(5)-q(2))*(q(5)  &
 -q(1)))
w(2,5)=  &
    ((((q(2)+q(1)-2.*q12)*q(3)+(q(1)-2.*q12)*q(2)-2.*q(1)*  &
 q12+3.*q12**2)*q(4)+((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*  &
 q12**2)*q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*  &
 q12**3)*q(6)+(((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*q12**2)*  &
 q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*q12**3)*  &
 q(4)-((2.*q(1)-3.*q12)*q(2)-3.*q(1)*q12+4.*q12**2)*q(3)*  &
 q12+(3.*q(1)-4.*q12)*q(2)*q12**2-4.*q(1)*q12**3+5.*q12**4  &
 )/(2.*(q(6)-q(5))*(q(5)-q(4))*(q(5)-q(3))*(q(5)-q(2))*(q(  &
 5)-q(1)))
w(3,5)=  &
    (((q(3)+q(2)+q(1)-3.*q12)*q(4)+(q(2)+q(1)-3.*q12)*q(3)  &
 +(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(6)+((q(2)+q(  &
 1)-3.*q12)*q(3)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)  &
 *q(4)+((q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(3)-3.*  &
 (q(1)-2.*q12)*q(2)*q12+6.*q(1)*q12**2-10.*q12**3)/(3.*(q(  &
 6)-q(5))*(q(5)-q(4))*(q(5)-q(3))*(q(5)-q(2))*(q(5)-q(1)))
w(4,5)=  &
    ((q(4)+q(3)+q(2)+q(1)-4.*q12)*q(6)+(q(3)+q(2)+q(1)-4.*  &
 q12)*q(4)+(q(2)+q(1)-4.*q12)*q(3)+(q(1)-4.*q12)*q(2)-4.*q  &
 (1)*q12+10.*q12**2)/(4.*(q(6)-q(5))*(q(5)-q(4))*(q(5)-q(3  &
 ))*(q(5)-q(2))*(q(5)-q(1)))
w(5,5)=  &
    (q(6)+q(4)+q(3)+q(2)+q(1)-5.*q12)/(5.*(q(6)-q(5))*(q(5  &
 )-q(4))*(q(5)-q(3))*(q(5)-q(2))*(q(5)-q(1)))
w(6,5)=  &
    1./(6.*(q(6)-q(5))*(q(5)-q(4))*(q(5)-q(3))*(q(5)-q(2))  &
 *(q(5)-q(1)))
w(1,6)=  &
    (-(q(5)-q12)*(q(4)-q12)*(q(3)-q12)*(q(2)-q12)*(q(1)-  &
 q12))/((q(6)-q(5))*(q(6)-q(4))*(q(6)-q(3))*(q(6)-q(2))*(q  &
 (6)-q(1)))
w(2,6)=  &
    (-((((q(2)+q(1)-2.*q12)*q(3)+(q(1)-2.*q12)*q(2)-2.*q(1  &
 )*q12+3.*q12**2)*q(4)+((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*  &
 q12**2)*q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*  &
 q12**3)*q(5)+(((q(1)-2.*q12)*q(2)-2.*q(1)*q12+3.*q12**2)*  &
 q(3)-(2.*q(1)-3.*q12)*q(2)*q12+3.*q(1)*q12**2-4.*q12**3)*  &
 q(4)-((2.*q(1)-3.*q12)*q(2)-3.*q(1)*q12+4.*q12**2)*q(3)*  &
 q12+(3.*q(1)-4.*q12)*q(2)*q12**2-4.*q(1)*q12**3+5.*q12**4  &
 ))/(2.*(q(6)-q(5))*(q(6)-q(4))*(q(6)-q(3))*(q(6)-q(2))*(q  &
 (6)-q(1)))
w(3,6)=  &
    (-(((q(3)+q(2)+q(1)-3.*q12)*q(4)+(q(2)+q(1)-3.*q12)*q(  &
 3)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(5)+((q(2)+  &
 q(1)-3.*q12)*q(3)+(q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**  &
 2)*q(4)+((q(1)-3.*q12)*q(2)-3.*q(1)*q12+6.*q12**2)*q(3)-  &
 3.*(q(1)-2.*q12)*q(2)*q12+6.*q(1)*q12**2-10.*q12**3))/(3.*  &
 (q(6)-q(5))*(q(6)-q(4))*(q(6)-q(3))*(q(6)-q(2))*(q(6)-q(1  &
 )))
w(4,6)=  &
    (-((q(4)+q(3)+q(2)+q(1)-4.*q12)*q(5)+(q(3)+q(2)+q(1)-  &
 4.*q12)*q(4)+(q(2)+q(1)-4.*q12)*q(3)+(q(1)-4.*q12)*q(2)-4.  &
 *q(1)*q12+10.*q12**2))/(4.*(q(6)-q(5))*(q(6)-q(4))*(q(6)-  &
 q(3))*(q(6)-q(2))*(q(6)-q(1)))
w(5,6)=  &
    (-(q(5)+q(4)+q(3)+q(2)+q(1)-5.*q12))/(5.*(q(6)-q(5))*(  &
 q(6)-q(4))*(q(6)-q(3))*(q(6)-q(2))*(q(6)-q(1)))
w(6,6)=  &
    (-1.)/(6.*(q(6)-q(5))*(q(6)-q(4))*(q(6)-q(3))*(q(6)-q(  &
 2))*(q(6)-q(1)))

goto 240

2000 continue
w( 1, 1)= 1./60.
w( 2, 1)= 2./360.
w( 3, 1)=-1./48.
w( 4, 1)=-1./144.
w( 5, 1)= 1./240.
w( 6, 1)= 1./720.
w( 1, 2)=-8./60.
w( 2, 2)=-25./360.
w( 3, 2)= 7./48.
w( 4, 2)= 11./144.
w( 5, 2)=-3./240.
w( 6, 2)=-5./720.
w( 1, 3)=37./60.
w( 2, 3)=245./360.
w( 3, 3)=-6./48.
w( 4, 3)=-28/144.
w( 5, 3)= 2./240.
w( 6, 3)= 10./720.
w( 1, 4)=37./60.
w( 2, 4)=-245./360.
w( 3, 4)=-6./48.
w( 4, 4)= 28./144.
w( 5, 4)= 2./240.
w( 6, 4)=-10./720.
w( 1, 5)=-8./60.
w( 2, 5)= 25./360.
w( 3, 5)= 7./48.
w( 4, 5)=-11./144.
w( 5, 5)=-3./240.
w( 6, 5)= 5./720.
w( 1, 6)= 1./60.
w( 2, 6)=-2./360.
w( 3, 6)=-1./48.
w( 4, 6)= 1./144.
w( 5, 6)= 1./240.
w( 6, 6)=-1./720.
do l=1,6
   do ll=1,6
      w(ll,l)=w(ll,l)/(x(ii+1)-x(ii))**(ll-1)
   enddo
enddo

240 continue
if(ityp.eq.1.and.ii.lt.llb+2)then
if(ii.eq.llb)then
do l=1,6
   w(l,4)=w(l,4)+w(l,2)
   w(l,5)=w(l,5)+w(l,1)
   w(l,2)=0.
   w(l,1)=0.
enddo
elseif(ii.eq.llb+1)then
do l=1,6
   w(l,5)=w(l,5)+w(l,1)
   w(l,1)=0.
enddo
endif
endif

return
end
