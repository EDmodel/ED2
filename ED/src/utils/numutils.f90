!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine azerov(n1)
implicit none
integer :: n,n1
real :: a1(n1),a2(n1),a3(n1),a4(n1),a5(n1)
entry azero(n1,a1)
   do n=1,n1
      a1(n)=0.
   enddo
return
entry azero2(n1,a1,a2)
   do n=1,n1
      a1(n)=0.
      a2(n)=0.
   enddo
return
entry azero3(n1,a1,a2,a3)
   do n=1,n1
      a1(n)=0.
      a2(n)=0.
      a3(n)=0.
   enddo
return
entry azero4(n1,a1,a2,a3,a4)
   do n=1,n1
      a1(n)=0.
      a2(n)=0.
      a3(n)=0.
      a4(n)=0.
   enddo
return
entry azero5(n1,a1,a2,a3,a4,a5)
   do n=1,n1
      a1(n)=0.
      a2(n)=0.
      a3(n)=0.
      a4(n)=0.
      a5(n)=0.
   enddo
return
end

![MLO ---- Similar to azerov, but for integers.
subroutine izerov(n1)
  implicit none
  integer :: n,n1
  integer :: ijk1(n1),ijk2(n1),ijk3(n1),ijk4(n1),ijk5(n1)
  entry izero(n1,ijk1)
     do n=1,n1
        ijk1(n)=0
     enddo
  return
  entry izero2(n1,ijk1,ijk2)
     do n=1,n1
        ijk1(n)=0
        ijk2(n)=0
     enddo
  return
  entry izero3(n1,ijk1,ijk2,ijk3)
     do n=1,n1
        ijk1(n)=0
        ijk2(n)=0
        ijk3(n)=0
     enddo
  return
  entry izero4(n1,ijk1,ijk2,ijk3,ijk4)
     do n=1,n1
        ijk1(n)=0
        ijk2(n)=0
        ijk3(n)=0
        ijk4(n)=0
     enddo
  return
  entry izero5(n1,ijk1,ijk2,ijk3,ijk4,ijk5)
     do n=1,n1
        ijk1(n)=0
        ijk2(n)=0
        ijk3(n)=0
        ijk4(n)=0
        ijk5(n)=0
     enddo
  return
end subroutine izerov

![MLO - Just to generate a matrix full of ones...
subroutine aonev(n1)
implicit none
integer :: n,n1
real :: a1(n1),a2(n1),a3(n1),a4(n1),a5(n1)
entry aone(n1,a1)
   do n=1,n1
      a1(n)=1.
   enddo
return
entry aone2(n1,a1,a2)
   do n=1,n1
      a1(n)=1.
      a2(n)=1.
   enddo
return
entry aone3(n1,a1,a2,a3)
   do n=1,n1
      a1(n)=1.
      a2(n)=1.
      a3(n)=1.
   enddo
return
entry aone4(n1,a1,a2,a3,a4)
   do n=1,n1
      a1(n)=1.
      a2(n)=1.
      a3(n)=1.
      a4(n)=1.
   enddo
return
entry aone5(n1,a1,a2,a3,a4,a5)
   do n=1,n1
      a1(n)=1.
      a2(n)=1.
      a3(n)=1.
      a4(n)=1.
      a5(n)=1.
   enddo
return
end subroutine aonev
!MLO]



subroutine ae1t0(n1,a,b,c)
implicit none
integer :: n1
real, dimension(n1) :: a,b
real :: c

integer :: n

do n=1,n1
   a(n)=b(n)*c
enddo

return 
end

subroutine ae1p1(npts,a,b,c)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts)
integer :: i
do i=1,npts
  a(i)=b(i)+c(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1m1(npts,a,b,c)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts)
integer :: i
do i=1,npts
  a(i)=b(i)-c(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine aen1(npts,a,b)
implicit none
integer :: npts
real :: a(npts),b(npts)
integer :: i
do i=1,npts
  a(i)=-b(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1t1(npts,a,b,c)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts)
integer :: i
do i=1,npts
  a(i)=b(i)*c(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1(npts,a,b)
implicit none
integer :: npts
real :: a(npts),b(npts)
integer :: i
do i=1,npts
  a(i)=b(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1tn1(npts,a,b,c)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts)
integer :: i
do i=1,npts
  a(i)=-b(i)*c(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1t0p1(npts,a,b,c,d)
implicit none
integer :: npts
real :: a(npts),b(npts),c,d(npts)
integer :: i
do i=1,npts
  a(i)=b(i)*c+d(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)+c(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3m3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)-c(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3t3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)*c(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1p1p1(npts,a,b,c,f)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts),f(npts)
integer :: i
do i=1,npts
  a(i)=b(i)+c(i)+f(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae1t1p1(npts,a,b,c,f)
implicit none
integer :: npts
real :: a(npts),b(npts),c(npts),f(npts)
integer :: i
do i=1,npts
  a(i)=b(i)*c(i)+f(i)
enddo
return
end
!
!     ******************************************************************
!
subroutine ae2(n2,n3,i1,i2,j1,j2,a,b)
implicit none
integer :: n2,n3,i1,i2,j1,j2
real :: a(n2,n3),b(n2,n3)
integer :: i,j
do j=j1,j2
  do i=i1,i2
    a(i,j)=b(i,j)
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3t3p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),f(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)*c(k,i,j)+f(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3t0p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c,f(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=b(k,i,j)*c+f(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine aen3t0p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c,f(n1,n2,n3)
integer :: i,j,k
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=-b(k,i,j)*c+f(k,i,j)
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae3m3d0(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),f
integer :: i,j,k

do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      a(k,i,j)=(b(k,i,j)-c(k,i,j))/f
    enddo
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine a3e2(n1,n2,n3,i1,i2,j1,j2,k,a,b)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k
real :: a(n1,n2,n3),b(n2,n3)
integer :: i,j
do j=j1,j2
  do i=i1,i2
    a(k,i,j)=b(i,j)
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine a3e1(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b)
   implicit none
   integer                  , intent(in)    :: n1,n2,n3
   integer                  , intent(in)    :: i1,i2,j1,j2,k1,k2
   real, dimension(n1,n2,n3), intent(inout) :: a
   real, dimension(n1)      , intent(in)    :: b
   integer :: i,j,k
   do j=j1,j2
      do i=i1,i2
         do k=k1,k2
            a(k,i,j)=b(k)
         end do
      end do
   end do
   return
end subroutine a3e1
!
!     ******************************************************************
!
subroutine a3e0(n1,n2,n3,i1,i2,j1,j2,k,a,b)
implicit none
integer :: n1,n2,n3,i1,i2,j1,j2,k
real :: a(n1,n2,n3),b
integer :: i,j
do j=j1,j2
  do i=i1,i2
    a(k,i,j)=b
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine alebl(n1,n2,n3,ka,kb,a,b)
implicit none
integer :: n1,n2,n3,ka,kb
real :: a(n1,n2,n3),b(n1,n2,n3)
integer :: i,j
do j=1,n3
  do i=1,n2
    a(ka,i,j)=b(kb,i,j)
  enddo
enddo
return
end
!
!     ******************************************************************
!
subroutine ae0(npts,a,b)
implicit none
integer :: npts
real :: a(npts),b
integer :: i
do i=1,npts
  a(i)=b
enddo
return
end
!
subroutine adivb(nnn,a,b,c)
implicit none
integer :: nnn,nn
real :: a(nnn),b(nnn),c(nnn)
do 1 nn=1,nnn
c(nn)=a(nn)/b(nn)
1 continue
return
end
subroutine atimb(nnn,a,b,c)
implicit none
integer :: nnn,nn
real :: a(nnn),b(nnn),c(nnn)
do 1 nn=1,nnn
c(nn)=a(nn)*b(nn)
1 continue
return
end

!     ******************************************************************

subroutine trid(var,cim1,ci,cip1,rhs,npts)

!     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
!     THE STANDARD QUICK METHOD
!
!     VAR   - VARIABLE BEING SOLVED FOR
!     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
!     CI    -   "     "       "       "  "  I     "
!     CIP1  -   "     "       "       "  "  I+1   "
!     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
!     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
!
!     WARNING: THE ARRAYS CIM1,CI, AND RHS ARE REUSED IN THIS ROUTINE.

implicit none
integer :: npts
real :: var(npts),cim1(npts),ci(npts),cip1(npts),rhs(npts)
integer :: k

cip1(1)=cip1(1)/ci(1)
do 10 k=2,npts
ci(k)=ci(k)-cim1(k)*cip1(k-1)
cip1(k)=cip1(k)/ci(k)
10 continue

rhs(1)=rhs(1)/ci(1)
do 20 k=2,npts
rhs(k)=(rhs(k)-cim1(k)*rhs(k-1))/ci(k)
20 continue

var(npts)=rhs(npts)
do 30 k=npts-1,1,-1
var(k)=rhs(k)-cip1(k)*var(k+1)
30 continue

return
end

!     ******************************************************************

subroutine trid2(var,cim1,ci,cip1,rhs,npts,scr1,scr2,scr3)

!     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
!     THE STANDARD QUICK METHOD
!
!     VAR   - VARIABLE BEING SOLVED FOR
!     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
!     CI    -   "     "       "       "  "  I     "
!     CIP1  -   "     "       "       "  "  I+1   "
!     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
!     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
!     SCR1  - SCRATCH ARRAY AT LEAST NPTS LONG
!     SCR2  - SCRATCH ARRAY "    "    "    "
!     SCR3  - SCRATCH ARRAY "    "    "    "
!
implicit none
integer :: npts
real :: var(npts),cim1(npts),ci(npts),cip1(npts),rhs(npts)
real :: scr1(npts),scr2(npts),scr3(npts)
integer :: k

scr1(1)=cip1(1)/ci(1)
scr2(1)=ci(1)
do 10 k=2,npts
scr2(k)=ci(k)-cim1(k)*scr1(k-1)
scr1(k)=cip1(k)/scr2(k)
10 continue

scr3(1)=rhs(1)/scr2(1)
do 20 k=2,npts
scr3(k)=(rhs(k)-cim1(k)*scr3(k-1))/scr2(k)
20 continue

var(npts)=scr3(npts)
do 30 k=npts-1,1,-1
var(k)=scr3(k)-scr1(k)*var(k+1)
30 continue

return
end

!     ******************************************************************

subroutine update(n,a,fa,dt)
implicit none
integer :: n,nn
real :: a(n),fa(n),dt
do 10 nn=1,n
  a(nn)=a(nn)+fa(nn)*dt
10 continue
return
end

!     ****************************************************************

subroutine accum(nxyz,arr1,arr2)
implicit none
integer :: nxyz,n
real :: arr1(nxyz),arr2(nxyz)
do n=1,nxyz
  arr1(n)=arr1(n)+arr2(n)
enddo
return
end

!     ******************************************************************

subroutine atob(n,a,b)
   implicit none
   integer, intent(in)                :: n
   real   , intent(in) , dimension(n) :: a
   real   , intent(out), dimension(n) :: b
   integer :: i
   do i=1,n
     b(i)=a(i)
   end do
   return
end subroutine atob

!     ******************************************************************

subroutine atob_log(n,a,b)
   implicit none
   integer, intent(in)                :: n
   logical, intent(in) , dimension(n) :: a
   logical, intent(out), dimension(n) :: b
   integer :: i
   do i=1,n
     b(i)=a(i)
   end do
   return
end subroutine atob_log

!     ******************************************************************

subroutine acnst(n,a,cnst)
implicit none
integer :: n
real :: a(n),cnst
integer :: nn
do 10 nn=1,n
  a(nn)=cnst
10 continue
return
end

!     ******************************************************************

real function valugp(n1,n2,n3,k,i,j,a)
implicit none
integer :: n1,n2,n3,k,i,j
real :: a(n1,n2,n3)
  valugp=a(k,i,j)
return
end

!     ******************************************************************

integer function ivalugp(n1,n2,n3,k,i,j,ia)
implicit none
integer :: n1,n2,n3,k,i,j
real :: ia(n1,n2,n3)
  ivalugp=int(ia(k,i,j))
return
end

integer function ibindec(str)
implicit none
character(len=*) :: str
integer :: inc,ic,l
ibindec=0
inc=1
l=len(str)
do ic=l,1,-1
   if(str(ic:ic).eq.'1') ibindec=ibindec+inc
   inc=inc+inc
enddo
return
end

!     ******************************************************************

integer function ibias(y,l,n)
implicit none
integer :: l,n
real :: y(l)
integer :: jd,jpow,jpow1
real :: ymin,ymax
ymin=1.e10
ymax=-1.e10
do jd=1,l
ymin=min(ymin,abs(y(jd)))
ymax=max(ymax,abs(y(jd)))
enddo
jpow=int(log10(ymax+1.e-20)+.999999)
jpow1=int(log10(ymax-ymin+1.e-20)+.999999)
ibias=2-jpow1
if(ibias+jpow.gt.4) ibias=4-jpow
10    continue
if(ibias+jpow.lt.4.and.ibias.lt.0)then
ibias=ibias+1
go to 10
endif
return
end

!     ******************************************************************

real function heav(x)
implicit none
real :: x
if(x.gt.0.)then
heav=1.
else
heav=0.
endif
return
end

!     ******************************************************************

integer function iprim(m)
implicit none
integer :: m,n
n=m
if(n.le.0)then
  print 1,n
1       format(' n=',i5,' in iprim')
  stop
endif
10    continue
if(mod(n,2).ne.0) go to 20
  n=n/2
  go to 10
20    continue
if(mod(n,3).ne.0) go to 30
  n=n/3
  go to 20
30    continue
if(mod(n,5).ne.0) go to 40
  n=n/5
  go to 30
40    continue
if(n.eq.1.and.mod(m,2).eq.0)then
  iprim=1
else
  iprim =0
endif
return
end

!     ******************************************************************

subroutine sort3(a,b,c,n)
   implicit none
   integer :: n
   real :: a(n),b(n),c(n)
   integer :: np1,k,i
   integer, external :: ismin
   real :: at,bt,ct
   np1=n+1
   
   do k=1,n
      i=ismin(np1-k,a(k),1)+k-1
      at=a(i)
      bt=b(i)
      ct=c(i)
      a(i)=a(k)
      b(i)=b(k)
      c(i)=c(k)
      a(k)=at
      b(k)=bt
      c(k)=ct
   end do
   return
end subroutine sort3
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine sorts the elements of vector a from smallest to largest.            !
!------------------------------------------------------------------------------------------!
subroutine sort_up(a,n)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                  :: n
   integer , intent(inout), dimension(n) :: a
   !----- Local variables. ----------------------------------------------------------------!
   logical ,                dimension(n) :: unlocked
   integer                               :: atmp
   integer                               :: imin
   integer                               :: k
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.

   do k=1,n
      imin        = minloc(a,1,unlocked)
      atmp        = a(imin)
      a(imin)     = a(k)
      a(k)        = atmp
      unlocked(k) = .false.
   end do
   return
end subroutine sort_up
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine sorts the elements of vector a from largest to smallest.            !
!------------------------------------------------------------------------------------------!
subroutine sort_down(a,n)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                  :: n
   integer , intent(inout), dimension(n) :: a
   !----- Local variables. ----------------------------------------------------------------!
   logical ,                dimension(n) :: unlocked
   integer                               :: atmp
   integer                               :: imax
   integer                               :: k
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.

   do k=1,n
      imax        = maxloc(a,1,unlocked)
      atmp        = a(imax)
      a(imax)     = a(k)
      a(k)        = atmp
      unlocked(k) = .false.
   end do
   return
end subroutine sort_down
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from smallest to largest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_up(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   real    , intent(in)  , dimension(nmax) :: variable
   integer , intent(out) , dimension(nmax) :: ranking
   !----- Local variables. ----------------------------------------------------------------!
   logical ,               dimension(nmax) :: unlocked
   integer                                 :: n
   integer                                 :: locmin
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.
   ranking (:) = 0
   do n=1,nmax
     locmin           = minloc(variable,1,unlocked)
     unlocked(locmin) = .false.
     ranking (locmin) = n
   end do

   return
end subroutine rank_up
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from largest to smallest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_down(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   real    , intent(in)  , dimension(nmax) :: variable
   integer , intent(out) , dimension(nmax) :: ranking
   !----- Local variables. ----------------------------------------------------------------!
   logical ,               dimension(nmax) :: unlocked
   integer                                 :: n
   integer                                 :: locmax
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.
   ranking (:) = 0
   do n=1,nmax
     locmax           = maxloc(variable,1,unlocked)
     unlocked(locmax) = .false.
     ranking (locmax) = n
   end do

   return
end subroutine rank_down
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function finds the element of the rank array that has a given rank.             !
!------------------------------------------------------------------------------------------!
integer function find_rank(ranking,nmax,rankarray)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)                  :: ranking
   integer, intent(in)                  :: nmax
   integer, intent(in), dimension(nmax) :: rankarray
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: n
   !---------------------------------------------------------------------------------------!
   find_rank=-1
   do n=1,nmax
      if (rankarray(n) == ranking) then
         find_rank=n
         return
      end if
   end do
   if (find_rank < 0) call fatal_error('Index not found','find_rank','numutils.f90')
   return
end function find_rank
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!

!
!-----------------------------------------------------------------------
!        The following functions are the FORTRAN replacements for the
!          CRAY intrinsic vector functions that perform various tasks.
!          Since they are non-existent on
!          other machines, they need to be replaced by calls to that
!          machines's functions or simulated in standard FORTRAN.
!
!       Most or all now exist in f90 and could be replaced in the 
!        few places in the code where they are used.
!-----------------------------------------------------------------------
!
!       Return sum of vector.
!
real function ssum(nn,vctr,inc)
implicit none
integer :: nn,inc
real :: vctr(*)
integer :: n,nnn
real :: sum
sum=0.
nnn=nn*inc
do 10 n=1,nnn,inc
  sum=sum+vctr(n)
10 continue
ssum=sum
return
end
!
!       Return sum of vector - Double precision.
!
real(kind=8) function dssum(nn,vctr,inc)
  implicit none
  integer :: nn,inc
  real(kind=8) :: vctr(*)
  integer :: n,nnn
  real(kind=8) :: sum
  sum=0.d0
  nnn=nn*inc
  do n=1,nnn,inc
    sum=sum+vctr(n)
  end do
  dssum=sum
  return
end function dssum
! +------------------------------------------------------------------+
!
!       return index of maximum of vector
!
integer function ismax(nn,vctr,inc)
implicit none
integer :: nn,inc
real :: vctr(*)
integer :: ism,nnn
real :: smax
ism=0
smax=-1e10
do 10 nnn=1,nn,inc
  if(vctr(nnn).gt.smax)then
    ism=nnn
    smax=vctr(nnn)
  endif
10 continue
ismax=ism
return
end
! +------------------------------------------------------------------+
!
!       return index of minimum of vector
!
integer function ismin(nn,vctr,inc)
implicit none
integer :: nn,inc
real :: vctr(*)
integer :: ism,nnn
real :: smin
ism=0
smin=1e10
do 10 nnn=1,nn,inc
  if(vctr(nnn).lt.smin)then
    ism=nnn
    smin=vctr(nnn)
  endif
10 continue
ismin=ism
return
end
! +-----------------------------------------------------------------+
!
!      Find the minimum of two double precision reals RGK 5-29-07
!
real(kind=8) function dmin2(val1,val2)
implicit none
real(kind=8) :: val1,val2
if(val1.lt.val2) then
   dmin2=val1
else
   dmin2=val2
endif
return
end
!+-------------------------------------------------------------------+
!
!     Find the maximum of two double precision reals RGK 5-29-07
!
real(kind=8) function dmax2(val1,val2)
implicit none
real(kind=8) :: val1,val2
if(val1.gt.val2) then
   dmax2=val1
else
   dmax2=val2
endif
return
end
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 => 0., else vct2.
!
real function cvmgp(vct1,vct2,vct3)
implicit none
real :: vct1,vct2,vct3
if(vct3.ge.0.0)then
cvmgp=vct1
else
cvmgp=vct2
endif
return
end
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 <= 0., else vct2.
!
real function cvmgm(vct1,vct2,vct3)
implicit none
real :: vct1,vct2,vct3
if(vct3.lt.0.0)then
cvmgm=vct1
else
cvmgm=vct2
endif
return
end
![MLO - Double precision version of cvmgp and cvmgm
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 => 0., else vct2.
!
real(kind=8) function dcvmgp(vct1,vct2,vct3)
  implicit none
  real(kind=8) :: vct1,vct2,vct3
  if(vct3.ge.0.0d0)then
    dcvmgp=vct1
  else
    dcvmgp=vct2
  endif
  return
end function dcvmgp
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 <= 0., else vct2.
!
real(kind=8) function dcvmgm(vct1,vct2,vct3)
  implicit none
  real(kind=8) :: vct1,vct2,vct3
  if(vct3.lt.0.0d0)then
    dcvmgm=vct1
  else
    dcvmgm=vct2
  endif
  return
end function dcvmgm
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 = 0., else vct2.
!
real function cvmgz(vct1,vct2,vct3)
implicit none
real :: vct1,vct2,vct3
if(vct3.eq.0.0)then
cvmgz=vct1
else
cvmgz=vct2
endif
return
end
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 ne 0., else vct2.
!
real function cvmgn(vct1,vct2,vct3)
implicit none
real :: vct1,vct2,vct3
if(vct3.ne.0.0)then
cvmgn=vct1
else
cvmgn=vct2
endif
return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the normalised value of a Gaussian distribution associated     !
! with the given cumulative distribution function. Fortran90 does not have the inverse of  !
! erf function, so we use the zeroin method (without IQI) to find the root.                !
!------------------------------------------------------------------------------------------!
real function cdf2normal(mycdf)
   use therm_lib, only : toler, maxfpo
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real   , intent(in) :: mycdf              ! The CDF we want to invert
   !----- Local variables -----------------------------------------------------------------!
   logical             :: converged          ! Convergence handle
   logical             :: bisection          ! Bisection check
   integer             :: it                 ! Iteration counter
   real                :: delta              ! Aux. variable for 2nd guess.
   real                :: normala            ! Aux. variables with previous guess
   real                :: normalz            ! Aux. variables with previous guess
   real                :: normalp            ! Aux. variables with previous guess
   real                :: normalc            ! Aux. variables with current  guess
   real                :: funa,funz          ! Function evaluation for bisection 
   real                :: func,funp          ! Function evaluation for secant 
   real                :: funnow             ! Function we want to find the root.
   !----- External functions --------------------------------------------------------------!
   real, external      :: cdf                ! Cumulative distribution function function
   !---------------------------------------------------------------------------------------!

   !----- Sanity check first... -----------------------------------------------------------!
   if (mycdf == 1.) then
      cdf2normal = 6.
      return
   elseif (mycdf == 0.) then
      cdf2normal = -6.
      return
   !----- This is a singularity for error checking, but we know it... ---------------------!
   elseif (mycdf == 0.5) then
      cdf2normal = 0.
      return
   elseif (mycdf < 0. .or. mycdf > 1.) then
      write (unit=*,fmt='(a,1x,es14.7)') 'WEIRD CDF!!! ',mycdf
      call fatal_error('Invalid input CDF!','cdf2normal','numutils.f90')
   end if

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write(unit=84,fmt='(a,1x,es14.7)') 'INPUT: ',mycdf
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   !----- First two guesses: no idea, just try -0.5 and 0.5 -------------------------------!
   normala  = -0.5
   funa     = cdf(normala) - mycdf
   normalz  =  0.5
   funz     = cdf(normalz) - mycdf
   



   !----- This is highly unlikely, but... -------------------------------------------------!
   if (funa == 0.) then
      cdf2normal = normala
      return
   elseif (funz == 0.) then
      cdf2normal = normalz
      return
   end if

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write(unit=84,fmt='(60a1)')        ('-',i=1,60)
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   
   !---------------------------------------------------------------------------------------!
   !   Finding the second guess for bisection. Checking whether the two guesses have op-   !
   ! posite signs, if not, fix the closest guess as normala and seek a better normalz.     !
   !---------------------------------------------------------------------------------------!
   if (funa * funz > 0.) then
      !----- Switching A and Z ------------------------------------------------------------!
      if (abs(funz) < abs(funa)) then
         funp    = funa
         funa    = funz
         funz    = funp
         normalp = normala
         normala = normalz
         normalz = normalp
      end if
      if (abs(funz-funa) < toler) then
         delta = 100.*toler
      else
         delta   = max(abs(funa * (normalz-normala)/(funz-funa)),100.*toler)
      end if
      bisection = .false.
      guessloop: do it=1,maxfpo
         normalz = normala + real((-1)**it * (it+3)/2) * delta
         funz    = cdf(normalz) - mycdf
         bisection = funa*funz < 0.
         if (bisection) exit guessloop
      end do guessloop
      if (.not. bisection) then
         write (unit=*,fmt='(a)') ' No second guess for you...'
         write (unit=*,fmt='(2(a,1x,es14.7))') 'normala=',normala,'funa=',funa
         write (unit=*,fmt='(2(a,1x,es14.7))') 'normalz=',normalz,'funz=',funz
         write (unit=*,fmt='(1(a,1x,es14.7))') 'delta=',delta
         call fatal_error('Failed finding the second guess for bisection'                    &
                       ,'cdf2normal','numutils.f90')
      end if
   else 
      bisection = .true.
   end if
   
   !----- Choose the closest one to be the previous guess ---------------------------------!
   if (abs(funa) < abs(funz)) then
      funp    = funa
      normalp = normala
   else
      funp    = funz
      normalp = normalz
   end if
   
   normalc    = (funz * normala - funa * normalz) / (funz-funa)
   func       = cdf(normalc) - mycdf


   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write(unit=84,fmt='(a,1x,i5,1x,a,1x,l1,1x,8(a,1x,es14.7,1x))')                          &
   !    '1STCALL: IT =', 0,'bisection=',bisection,'normala=',normala,'normalz=',normalz     &
   !                ,'normalp=',normalp,'normalc=',normalc,'funa=',funa,'funz=',funz        &
   !                ,'funp=',funp,'func=',func
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


   !----- Looping until convergence is achieved -------------------------------------------!
   converged = .false.
   bisection = .false.
   itloop: do it=1,maxfpo
   
      !------------------------------------------------------------------------------------!
      ! e1. Deciding whether to go with bisection or not. I should go with bisection if    !
      !     the secant is dangerously small (derivative too flat, which causes diver-      !
      !     gence). Also if it didn't converge fast with secant, fall back to bisection.   !
      !------------------------------------------------------------------------------------!
      bisection = it > maxfpo / 4 .or. abs(func-funp) < toler

      !------------------------------------------------------------------------------------!
      ! e2. Setting the new guess. Still not sure with which method I should go, so estab- !
      !     lish the new guess using secant. If the guess is outside the range defined by  !
      !     the A Z pair, use bisection this time.                                         !
      !------------------------------------------------------------------------------------!
      if (.not.bisection) then
          cdf2normal = (funp*normalc - func*normalp) / (funp-func)
          bisection  = abs(cdf2normal-normala) > abs(normalz-normala) .or.                 &
                       abs(cdf2normal-normalz) > abs(normalz-normala)
      end if
      if (bisection) cdf2normal = 0.5 * (normala+normalz)

      !------------------------------------------------------------------------------------!
      ! e3. Finding the new function evaluation.                                           !
      !------------------------------------------------------------------------------------!
      funnow = cdf(cdf2normal) - mycdf



      !------------------------------------------------------------------------------------!
      ! e4. Testing for convergence, depending on the method.                              !
      !------------------------------------------------------------------------------------!
      if (funnow == 0.) then
         converged = .true.
      elseif (bisection) then 
         converged = abs(cdf2normal-normala) < toler*abs(cdf2normal)
      else
         converged = abs(cdf2normal-normalc) < toler*abs(cdf2normal)
      end if
      !----- Found a good set, leaving... -------------------------------------------------!
      if (converged) exit itloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! e5. Checking which side from the A Z pair I can update, based on funnow.           !
      !------------------------------------------------------------------------------------!
      if (funnow*funa < 0.) then
         funz    = funnow
         normalz = cdf2normal
      else
         funa    = funnow
         normala = cdf2normal
      end if

      !------------------------------------------------------------------------------------!
      ! e6. Updating the Previous-Current pair for the next secant attempt.                !
      !------------------------------------------------------------------------------------!
      normalp = normalc
      funp    = func
      normalc = cdf2normal
      func    = funnow


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=84,fmt='(a,1x,i5,1x,a,1x,l1,1x,8(a,1x,es14.7,1x))')                       &
      !    'LOOP   : IT =',it,'bisection=',bisection,'normala=',normala,'normalz=',normalz  &
      !                ,'normalp=',normalp,'normalc=',normalc,'funa=',funa,'funz=',funz     &
      !                ,'funp=',funp,'func=',func
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


   end do itloop

   if (.not. converged) then
      write (unit=*,fmt='(a)') '-----------------------------------------------------------'
      write (unit=*,fmt='(a)') ' Normalised value from CDF didn''t converge!!!'
      write (unit=*,fmt='(a,1x,i5,1x,a)')    ' I gave up, after',maxfpo,'iterations...'
      write (unit=*,fmt='(a)')               ' '
      write (unit=*,fmt='(a)')               ' Input CDF.'
      write (unit=*,fmt='(a,1x,f12.4)' )     '  # mycdf   =',mycdf
      write (unit=*,fmt='(a)')               ' '
      write (unit=*,fmt='(a)')               '  Last iteration outcome.'
      write (unit=*,fmt='(2(a,1x,es14.7))' ) '  # Normala   =',normala,' # Funa=',funa
      write (unit=*,fmt='(2(a,1x,es14.7))' ) '  # Normalz   =',normalz,' # Funz=',funz
      write (unit=*,fmt='(2(a,1x,es14.7))' ) '  # Normalc   =',normalc,' # Func=',func
      write (unit=*,fmt='(2(a,1x,es14.7))' ) '  # Normalp   =',normalp,' # Funp=',funp
      write (unit=*,fmt='(2(a,1x,es14.7))' ) '  # Normal    =',cdf2normal,' # Fun =',funnow
      write (unit=*,fmt='(a)')               '  Errors.'
      write (unit=*,fmt='(a,1x,es14.7)')     '  # Secant    ='                             &
                                             ,abs(cdf2normal-normalc)/abs(cdf2normal)
      write (unit=*,fmt='(a,1x,es14.7)')     '  # Bisection ='                             &
                                             ,abs(cdf2normal-normala)/abs(cdf2normal) 
      write (unit=*,fmt='(a)') '-----------------------------------------------------------'
      
      call fatal_error('Failed finding normalised value from CDF!!!'                         &
                                      ,'cdf2normal','numutils.f90')
   !else

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=84,fmt='(a,1x,i5,1x,a,1x,l1,1x,8(a,1x,es14.7,1x))')                       &
      !    'ANSWER : IT =',it,'bisection=',bisection,'normala=',normala,'normalz=',normalz  &
      !                ,'normalp=',normalp,'normalc=',normalc,'funa=',funa,'funz=',funz     &
      !                ,'funp=',funp,'func=',func
      !write(unit=84,fmt='(60a1)')        ('-',i=1,60)
      !write(unit=84,fmt='(a)')           ' '
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   end if
   return
end function cdf2normal
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the cumulative distribution function of a given normalised     !
! value. Fortran90 (at least gfortran, intel, and pgf) contains a built-in erf function.   !
! If your compiler doesn't have, it may be needed to use iterative methods instead.        !
!------------------------------------------------------------------------------------------!
real function cdf(normal)
   use consts_coms, only : srtwoi
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in) :: normal
   !----- External function, the error function. PGI doesn't have the intrinsic erf -------!
   real, external   :: errorfun
   !---------------------------------------------------------------------------------------!
   cdf = 0.5 * (1. + errorfun(normal * srtwoi))
   return
end function cdf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function simply computes the cubic root of all numbers, including the negative    !
! ones.                                                                                    !
!------------------------------------------------------------------------------------------!
real function cbrt(x)
   use consts_coms, only: onethird
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in) :: x
   !---------------------------------------------------------------------------------------!

   if (x > 0.0) then
     cbrt=x**onethird
   else
     cbrt=-((-x)**onethird)
   end if

   return
end function cbrt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function simply computes the cubic root of all numbers, including the negative    !
! ones, for a double precision number.                                                     !
!------------------------------------------------------------------------------------------!
real(kind=8) function cbrt8(x)
   use consts_coms, only: onethird8
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   !---------------------------------------------------------------------------------------!
   if (x > 0.d0) then
     cbrt8 = x**onethird8
   else
     cbrt8 = -((-x)**onethird8)
   end if 

   return
end function cbrt8
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the error function of a variable x. Some fortran distributions !
! have the intrinsic erf function, but not all of them come with (pgi seems not to have).  !
! This function was tested against the intrinsic function for ifort and results were       !
! within the tolerance for all values tested.                                              !
!    Here we will use the adaptive quadrature method to find it. The method implemented    !
! here was entirely based on the algorithm 4.3 from:                                       !
!     BURDEN, R.L.;l FAIRES, J.D., 2005. Numerical Analysis, 8th edition. Thomson          !
!          Brooks/Cole, Belmont, CA.                                                       !
!------------------------------------------------------------------------------------------!
real function errorfun(x)
   use consts_coms, only: onethird,onesixth,sqrtpii
   use therm_lib , only: toler,maxlev 
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in) :: x
   !----- Local variables -----------------------------------------------------------------!  
   real, dimension(8)      :: vv !   Temporary storage area
   real, dimension(maxlev) :: tol,aa,hh,fa,fb,fc,ss,ll
   real                    :: a0,b0,fd,fe,s1,s2,app,toomanylevels
   integer                 :: i
   !----- Function ------------------------------------------------------------------------!
   real, external          :: expmsq
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Return right away if x is too large. erf quickly converges to 1/-1 so no need to   !
   ! compute when |x| > 9., even because the exponential would be too large.               !
   !---------------------------------------------------------------------------------------!
   if (abs(x) >= 9.) then
     errorfun = sign(1.,x)
     return
   end if
   !---------------------------------------------------------------------------------------!


   !----- Initialise integral and integration limits --------------------------------------!
   app           = 0.0            !----- This is the integral's approximate value
   a0            = min(0.,x)    !----- Lower bound
   b0            = max(0.,x)    !----- Upper bound 
   toomanylevels = real(maxlev) !----- Just for if test.
   !---------------------------------------------------------------------------------------!



   !----- Initialise scratch arrays -------------------------------------------------------!
   i      = 1
   tol(i) = 10.*toler
   aa(i)  = a0
   hh(i)  = 0.5*(b0-a0)
   fa(i)  = expmsq(aa(i))
   fc(i)  = expmsq(aa(i)+hh(i))
   fb(i)  = expmsq(b0)
   !----- Approximation from Simpson's method for entire interval -------------------------!
   ss(i)  = onethird * hh(i) * (fa(i) + 4.*fc(i) + fb(i)) 
   ll(i)  = 1.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Entering the iterative loop. I'm not using do while because people say that it is    !
   ! deprecated in fortran90.                                                              !
   !---------------------------------------------------------------------------------------!
   iterloop: do
      if (i == 0) exit iterloop

      !----- Computing Simpson's method for halves of subintervals ------------------------!
      fd = expmsq(aa(i)+0.5*hh(i))
      fe = expmsq(aa(i)+1.5*hh(i))
      s1 = onesixth * hh(i) * (fa(i) + 4.*fd + fc(i))
      s2 = onesixth * hh(i) * (fc(i) + 4.*fe + fb(i))

      !----- Saving data at this level ----------------------------------------------------!
      vv(1) = aa(i)
      vv(2) = fa(i)
      vv(3) = fc(i)
      vv(4) = fb(i)
      vv(5) = hh(i)
      vv(6) = tol(i)
      vv(7) = ss(i)
      vv(8) = ll(i)
      
      !----- Delete the level -------------------------------------------------------------!
      i     = i -1
      
      !----- Checking convergence ---------------------------------------------------------!
      if (abs(s1+s2-vv(7)) < vv(6)) then
         app = app + s1 + s2
      elseif (vv(8) >= toomanylevels) then
         call fatal_error('Too many levels, it won''t converge!','errorfun','numutils.f90')
      !----- Add one level ----------------------------------------------------------------!
      else
         !----- Data for the right half subinterval ---------------------------------------!
         i      = i + 1
         aa(i)  = vv(1) + vv(5)
         fa(i)  = vv(3)
         fc(i)  = fe
         fb(i)  = vv(4)
         hh(i)  = 0.5 * vv(5)
         tol(i) = 0.5 * vv(6) 
         ss(i)  = s2
         ll(i)  = vv(8) 
         !----- Data for the left half subinterval ----------------------------------------!
         i      = i + 1
         aa(i)  = vv(1)
         fa(i)  = vv(2)
         fc(i)  = fd
         fb(i)  = vv(3)
         hh(i)  = hh(i-1)
         tol(i) = tol(i-1)
         ss(i)  = s1
         ll(i)  = ll(i-1)
      end if
   end do iterloop
   
   !---------------------------------------------------------------------------------------!
   !     Computing the function. erf is known to be bounded by -1.and 1 and here we will   !
   ! ensure that that will be the case. Due to numerical precision the final value of      !
   ! x) may be such that 1 <= erf(x) <= 1.+toler. Although within the tolerance, this  !
   ! can cause problems elsewhere.                                                         !
   !---------------------------------------------------------------------------------------!
   errorfun = max(-1.,min(1.,2. * app * sqrtpii * sign(1.,x)))
   
   return
end function errorfun
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function simply computes f(x) = exp (-x²). Just to make errorfun neater.         !
!------------------------------------------------------------------------------------------!
real function expmsq(x)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   real, intent(in) :: x
   !---------------------------------------------------------------------------------------!
   
   expmsq=exp(-x*x)

   return
end function expmsq
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the expect value based on the normal distribution, with mean   !
! xmean and standard deviation sigmax, but only for values above xmin.  This is actually a !
! generalisation of the concept of mean, but for a subset of the distribution function.    !
! Let p(x) be the normal PDF of x for a mean value of xmean and a standard deviation       !
! sigmax.  The expected value will solve:                                                  !
!                                                                                          !
!                 Inf                                                                      !
!              INT     x p(x) dx                                                           !
!                 xmin                                                                     !
! xexpected = -------------------                                                          !
!                 Inf                                                                      !
!              INT     p(x) dx                                                             !
!                 xmin                                                                     !
!                                                                                          !
!     When xmin tends to -Infinity, xexpected tends to xmean. If xmin tends to Infinity,   !
! xexpected will tend to xmin.                                                             !
!------------------------------------------------------------------------------------------!
real function expected(xmin,xmean,sigmax)
   use consts_coms, only : srtwo      & ! intent(in)
                         , srtwoi     & ! intent(in)
                         , sqrttwopi  & ! intent(in)
                         , sqrthalfpi ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in) :: xmin
   real, intent(in) :: xmean
   real, intent(in) :: sigmax
   !----- Local variables. ----------------------------------------------------------------!
   real             :: xnorm
   real             :: expnorm
   !----- Local constants. ----------------------------------------------------------------!
   real, parameter  :: xnormmin = -5. ! This makes a cdf close to machine epsilon
   real, parameter  :: xnormmax =  5. ! This makes a cdf close to 1 - machine epsilon
   !----- External functions. -------------------------------------------------------------!
   real, external   :: cdf
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     First we find the normalised xmin, which goes into the PDF and CDF functions.     !
   !---------------------------------------------------------------------------------------!
   xnorm = (xmin - xmean) / sigmax

   !---------------------------------------------------------------------------------------!
   !     Then we integrate the probability distribution function above xnorm. If the cdf   !
   ! is very small or very close to one, we apply what we know about the limits, and skip  !
   ! the calculation.                                                                      !
   !---------------------------------------------------------------------------------------!
   if (xnorm >= xnormmax) then
      !----- cdf is too close to 1, apply the limit of x -> infinity. ---------------------!
      expected = xmin

   elseif (xnorm <= xnormmin) then
      !----- cdf is too close to 0, apply the limit of x -> - infinity. -------------------!
      expected = xmean

   else
      !----- Nice range, let's find the results. ------------------------------------------!
      expnorm   = exp(- 0.5 * xnorm * xnorm) / (sqrttwopi * (1. - cdf(xnorm)))
      expected  = xmean + sigmax * expnorm
   end if

   return
end function expected
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function converts the double precision variable into single, in a way to prevent !
! floating point exception when they are tiny.  In case the number is too small, less than !
! off, then the output value is flushed to 0.                                              !
!------------------------------------------------------------------------------------------!
real function sngloff(x,off)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   real(kind=8), intent(in) :: off
   !---------------------------------------------------------------------------------------!
   
   if (abs(x) < off) then
      sngloff = 0.
   else
      sngloff = sngl(x)
   end if
   return
end function sngloff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine returns the accumulated sum of a given vector.                        !
!------------------------------------------------------------------------------------------!
subroutine cumsum(nsiz,vec)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: nsiz
   real   , dimension(nsiz), intent(inout) :: vec
   !----- Local variables. ----------------------------------------------------------------!
   integer                :: n
   !---------------------------------------------------------------------------------------!
   do n=2,nsiz
      vec(n) = vec(n) + vec(n-1)
   end do

   return
end subroutine cumsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Check two corresponding real arrays and see values are close enough.                     !
!------------------------------------------------------------------------------------------!
integer function check_real(xx,yy,nsiz)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                 , intent(in) :: nsiz
   real   , dimension(nsiz), intent(in) :: xx
   real   , dimension(nsiz), intent(in) :: yy
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: n
   real                                 :: tol
   !---------------------------------------------------------------------------------------!
   
   !----- We first assume that we don't find any similar element. -------------------------!
   check_real = 0

   tol = min( (maxval(xx)-minval(xx)),(maxval(yy)-minval(yy)) ) * .0001
   do n = 1, nsiz
      if (abs(xx(n)-yy(n)) > tol) then
         check_real=n
         return
      end if
   end do
   return
end function check_real
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the linear system AA . X = Y for given AA and Y, using    !
! the Gaussian elimination method with partial pivoting and back-substitution.             !
!------------------------------------------------------------------------------------------!
subroutine lisys_solver(nsiz,AA,Y,X,sing)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                      , intent(in)  :: nsiz  ! matrix and vector size
   real   , dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
   real   , dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
   real   , dimension(nsiz)     , intent(out) :: X     ! unknown vector
   logical                      , intent(out) :: sing  ! The matrix was singular      [T|F]
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
   real   , dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
   real   , dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
   real                                       :: pivot  ! The pivot
   real                                       :: multip ! Multiplier
   integer                                    :: r      ! Row index
   integer                                    :: b      ! Row below index
   integer                                    :: c      ! Column index
   integer                                    :: p      ! Pivot index
   real                                       :: dumsca ! Dummy scalar, for row swapping
   !----- Local parameters. ---------------------------------------------------------------!
   real                         , parameter   :: tinyoff=1.e-20
   !---------------------------------------------------------------------------------------!
   
   !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
   EE(:,:) = AA(:,:)
   Z (:)   = Y (:)
   dumvec  = 0.
   dumsca  = 0.
   !---------------------------------------------------------------------------------------!
   !     We initialise X with a huge, non-sense value, which will become the answer when   !
   ! the matrix is singular.                                                               !
   !---------------------------------------------------------------------------------------!
   X (:)   = -huge(1.)
   !----- We first assume that everything will be fine. -----------------------------------!
   sing    = .false.

   !---------------------------------------------------------------------------------------!
   ! 1. Main elimination loop, done row by row.                                            !
   !---------------------------------------------------------------------------------------!
   elimloop: do r = 1, nsiz-1
      !------ 1a. Finding the largest element, which will become our pivot ----------------!
      p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)
      
      pivot = maxval(abs(EE(r:nsiz,r)))
      !------------------------------------------------------------------------------------!
      ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
      !     singular or almost singular, and we cannot solve it, so we switch the flag and !
      !     return.                                                                        !
      !------------------------------------------------------------------------------------!
      if (pivot < tinyoff) then
         sing = .true.
         return
      end if
      
      !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
      if (p /= r) then
         dumvec(r:nsiz) = EE(r,r:nsiz)
         dumsca         = Z(r)
         EE(r,r:nsiz)   = EE(p,r:nsiz)
         Z(r)           = Z(p)
         EE(p,r:nsiz)   = dumvec(r:nsiz)
         Z(p)           = dumsca
      end if

      !------------------------------------------------------------------------------------!
      ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
      !      zero (we won't compute that, but they will be.).                              !
      !------------------------------------------------------------------------------------!
      belowloop: do b=r+1,nsiz
         multip = EE(b,r)/EE(r,r)
         EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
         Z(b)         = Z(b)         - multip * Z(r)
      end do belowloop
   end do elimloop

   !---------------------------------------------------------------------------------------!
   ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
   !    check the last pivot too.                                                          ! 
   !---------------------------------------------------------------------------------------!
   if (abs(EE(nsiz,nsiz)) < tinyoff) then
      sing = .true.
      return
   end if

   !---------------------------------------------------------------------------------------!
   ! 3. We now perform the back-substitution, to find the solution.                        !
   !---------------------------------------------------------------------------------------!
   X(nsiz) = Z(nsiz) / EE(nsiz,nsiz)
   backsubloop: do r=nsiz-1,1,-1
      b    = r+1
      X(r) = (Z(r) - sum(EE(r,b:nsiz)*x(b:nsiz))) / EE(r,r)
   end do backsubloop

   return
end subroutine lisys_solver
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This function checks whether a number is finite or not.  This test will return true  !
! if the number is a valid one, and false if the number is either +Infinity, -Infinity, or !
! NaN.                                                                                     !
!------------------------------------------------------------------------------------------!
logical function is_finite(number)
   implicit none
   real, intent(in) :: number
   real, parameter  :: largeneg = -huge(1.)
   real, parameter  :: largepos =  huge(1.)
   is_finite = number >= largeneg .and. number <= largepos
   return
end function is_finite
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function checks whether a number is finite or not.  This test will return true  !
! if the number is a valid one, and false if the number is either +Infinity, -Infinity, or !
! NaN.                                                                                     !
!------------------------------------------------------------------------------------------!
logical function is_finite8(number)
   implicit none
   real(kind=8), intent(in) :: number
   real(kind=8), parameter  :: largeneg = -huge(1.d0)
   real(kind=8), parameter  :: largepos =  huge(1.d0)
   is_finite8 = number >= largeneg .and. number <= largepos
   return
end function is_finite8
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine extracts the diagonal of a matrix.                                   !
!------------------------------------------------------------------------------------------!
subroutine diagon(nsiz,mat,vec)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                      , intent(in)  :: nsiz
   real   , dimension(nsiz,nsiz), intent(in)  :: mat
   real   , dimension(nsiz)     , intent(out) :: vec
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: n
   !---------------------------------------------------------------------------------------!

   do n=1,nsiz
      vec(n) = mat(n,n)
   end do

   return
end subroutine diagon
!==========================================================================================!
!==========================================================================================!
