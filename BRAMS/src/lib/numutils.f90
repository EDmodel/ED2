!==========================================================================================!
!==========================================================================================!
!  Change Log                                                                              !
!  2.0.0                                                                                   !
!                                                                                          !
!------------------------------------------------------------------------------------------!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      This sub-routine flushes all elements of this array to zero.  Legacy from the old   !
! code, when vector operations didn't exist.                                               !
!------------------------------------------------------------------------------------------!
subroutine azero(nmax,arr)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)  :: nmax
   real   , dimension(nmax), intent(out) :: arr
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: n
   !---------------------------------------------------------------------------------------!

   do n=1,nmax
      arr(n) = 0.
   end do

   return
end subroutine azero
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      This sub-routine flushes all elements of this array to zero.  Legacy from the old   !
! code, when vector operations didn't exist.  The only difference between this one and     !
! azero is that the input vector here is integer.                                          !
!------------------------------------------------------------------------------------------!
subroutine izero(nmax,arr)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)  :: nmax
   integer, dimension(nmax), intent(out) :: arr
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: n
   !---------------------------------------------------------------------------------------!

   do n=1,nmax
      arr(n) = 0
   end do

   return
end subroutine izero
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      This sub-routine flushes all elements of this array to one.  Legacy from the old    !
! code, when vector operations didn't exist.                                               !
!------------------------------------------------------------------------------------------!
subroutine aone(nmax,arr)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)  :: nmax
   real   , dimension(nmax), intent(out) :: arr
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: n
   !---------------------------------------------------------------------------------------!

   do n=1,nmax
      arr(n) = 1.
   end do

   return
end subroutine aone
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
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
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine updates the matrix A using the tendency DADT applied over a time     !
! step of size DT.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine update(nsiz,a,dadt,dt)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: nsiz
   real   , dimension(nsiz), intent(inout) :: a
   real   , dimension(nsiz), intent(in)    :: dadt
   real                    , intent(in)    :: dt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: n
   !---------------------------------------------------------------------------------------!
   do n=1,nsiz
     a(n)=a(n)+dadt(n)*dt
   end do
   return
end subroutine update
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will add contrib to the total array.                                 !
!------------------------------------------------------------------------------------------!
subroutine accum(nsiz,total,contrib)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: nsiz
   real   , dimension(nsiz), intent(inout) :: total
   real   , dimension(nsiz), intent(in)    :: contrib
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: n
   !---------------------------------------------------------------------------------------!
   do n=1,nsiz
     total(n) = total(n) + contrib(n)
   end do
   return
end subroutine accum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies source array to the destination array.                 !
!------------------------------------------------------------------------------------------!
subroutine atob(nsiz,a,b)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nsiz
   real   , dimension(nsiz) , intent(in)    :: a
   real   , dimension(nsiz) , intent(inout) :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: n
   !---------------------------------------------------------------------------------------!

   do n=1,nsiz
      b(n)=a(n)
   end do

   return
end subroutine atob
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine atob_log(n,a,b)
   implicit none
   integer, intent(in)                  :: n
   logical, intent(in)   , dimension(n) :: a
   logical, intent(inout), dimension(n) :: b
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
subroutine rank_up_r(nmax,variable,ranking)
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
end subroutine rank_up_r
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from largest to smallest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_down_r(nmax,variable,ranking)
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
end subroutine rank_down_r
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from smallest to largest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_up_i(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   integer , intent(in)  , dimension(nmax) :: variable
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
end subroutine rank_up_i
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from largest to smallest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_down_i(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   integer , intent(in)  , dimension(nmax) :: variable
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
end subroutine rank_down_i
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
   if (find_rank < 0) call abort_run('Index not found','find_rank','numutils.f90')
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
      call abort_run('Invalid input CDF!','cdf2normal','numutils.f90')
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
         call abort_run('Failed finding the second guess for bisection'                    &
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
      
      call abort_run('Failed finding normalised value from CDF!!!'                         &
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
   use rconstants, only : srtwoi
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
   use rconstants, only: onethird
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
   use rconstants, only: onethird8
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
   use rconstants, only: onethird,onesixth,sqrtpii
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
         call abort_run('Too many levels, it won''t converge!','errorfun','numutils.f90')
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
   ! x) may be such that 1 <= erf(x) <= 1.+toler. Although within the tolerance, this      !
   ! can cause problems elsewhere.                                                         !
   !---------------------------------------------------------------------------------------!
   errorfun = max(-1.,min(1.,2. * app * sqrtpii * sign(1.,x)))
   
   return
end function errorfun
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
   use rconstants, only : srtwo      & ! intent(in)
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
!    This function simply computes f(x) = exp (-x^2). Just to make errorfun neater.        !
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
!     This subroutine is based on:                                                         !
!                                                                                          !
! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery: 1992. Numerical recipes !
!    in Fortran 77.  Cambridge University Press.                                           !
!------------------------------------------------------------------------------------------!
subroutine lisys_solver(nsiz,AA,Y,X,sing)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)  :: nsiz   ! matrix and vector size
   real(kind=4), dimension(nsiz,nsiz), intent(in)  :: AA     ! matrix
   real(kind=4), dimension(nsiz)     , intent(in)  :: Y      ! right-hand side vector
   real(kind=4), dimension(nsiz)     , intent(out) :: X      ! unknown vector
   logical                           , intent(out) :: sing   ! The matrix is singular [T|F]
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
   real(kind=4), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
   real(kind=4), dimension(nsiz)                   :: dumvec ! Dummy vector (row swapping)
   real(kind=4)                                    :: pivot  ! The pivot
   real(kind=4)                                    :: multip ! Multiplier
   integer                                         :: r      ! Row index
   integer                                         :: b      ! Row below index
   integer                                         :: c      ! Column index
   integer                                         :: p      ! Pivot index
   real(kind=4)                                    :: dumsca ! Dummy scalar (row swapping)
   !----- Local parameters. ---------------------------------------------------------------!
   real(kind=4)                      , parameter   :: tinyoff=1.e-20
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
!     This subroutine is the double precision version of the linear system solver above.   !
! It will solve the linear system AA . X = Y for given AA and Y, using the Gaussian        !
! elimination method with partial pivoting and back-substitution.  This subroutine is      !
! based on:                                                                                !
!                                                                                          !
! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery: 1992. Numerical recipes !
!    in Fortran 77.  Cambridge University Press.                                           !
!------------------------------------------------------------------------------------------!
subroutine lisys_solver8(nsiz,AA,Y,X,sing)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)  :: nsiz  ! matrix and vector size
   real(kind=8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
   real(kind=8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
   real(kind=8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
   logical                           , intent(out) :: sing  ! The matrix was singular [T|F]
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
   real(kind=8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
   real(kind=8), dimension(nsiz)                   :: dumvec ! Dummy vector (row swapping)
   real(kind=8)                                    :: pivot  ! The pivot
   real(kind=8)                                    :: multip ! Multiplier
   integer                                         :: r      ! Row index
   integer                                         :: b      ! Row below index
   integer                                         :: c      ! Column index
   integer                                         :: p      ! Pivot index
   real(kind=8)                                    :: dumsca ! Dummy scalar (row swapping)
   !----- Local parameters. ---------------------------------------------------------------!
   real(kind=8)                      , parameter   :: tinyoff=1.d-20
   !---------------------------------------------------------------------------------------!
   
   !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
   EE(:,:) = AA(:,:)
   Z (:)   = Y (:)
   dumvec  = 0.d0
   dumsca  = 0.d0
   !---------------------------------------------------------------------------------------!
   !     We initialise X with a huge, non-sense value, which will become the answer when   !
   ! the matrix is singular.                                                               !
   !---------------------------------------------------------------------------------------!
   X (:)   = -huge(1.d0)
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
end subroutine lisys_solver8
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





!==========================================================================================!
!==========================================================================================!
!     EIFUN8 -- This function computes the exponential integral function, defined by       !
!                                                                                          !
!                        x_                                                                !
!                        |    exp(t)                                                       !
!              Ei(x) =   |   -------- dt                                                   !
!                       _|      t                                                          !
!                      -Inf                                                                !
!                                                                                          !
!     This function checks for two approaches: series expansion, which typically works     !
! best when x is small, and the asymptotic expansion, which typically works best when x is !
! large. The approach selects the smallest result (in absolute numbers) as the most        !
! accurate method. Both the series expansion and the asymptotic expansion are provided in  !
! AS72. This approach also checks for some other edge cases, and ignores the results when  !
! the value is very negative.                                                              !
!                                                                                          !
! Reference:                                                                               !
!                                                                                          !
! Abramowitz, M., and I. A. Stegun, Eds., 1972: Handbook of mathematical functions with    !
!    formulas, graphs, and mathematical tables. 10th ed., No. 55, Applied Mathematics      !
!    Series, National Bureau of Standards, Washington, DC, USA (AS72).                     !
!                                                                                          !
!------------------------------------------------------------------------------------------!
real(kind=8) function eifun8(x)
   use rconstants, only : euler_gam8   & ! intent(in)
                        , lnexp_min8   & ! intent(in)
                        , lnexp_max8   & ! intent(in)
                        , tiny_num8    & ! intent(in)
                        , almost_zero8 ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: usum
   real(kind=8)             :: uxk
   real(kind=8)             :: uxkm1
   real(kind=8)             :: vsum
   real(kind=8)             :: vxkm1
   real(kind=8)             :: vxk
   real(kind=8)             :: ei_series
   real(kind=8)             :: ei_asymptote
   integer                  :: k
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter  :: discard8      = 1.0d+36
   integer     , parameter  :: maxiter       = 100
   !----- Polynomial coefficients. --------------------------------------------------------!
   real(kind=8), dimension(4), parameter :: apoly = (/ 8.5733287401d+00, 1.8059015973d+01  &
                                                     , 8.6347608925d+00, 2.6777373430d-01 /)
   real(kind=8), dimension(4), parameter :: bpoly = (/ 9.5733223454d+00, 2.5632956149d+01  &
                                                     , 2.1099653083d+01, 3.9584969228d+00 /)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check what to do depending on the value of x.                                     !
   !---------------------------------------------------------------------------------------!
   if (x == 0.d0) then
      !------------------------------------------------------------------------------------!
      !     Zero.  This is a singularity and the user should never call it in this case.   !
      !------------------------------------------------------------------------------------!
      stop 'Exponential integral cannot be solved for x = 0.'
      !------------------------------------------------------------------------------------!
   elseif (x <= lnexp_min8) then
      !------------------------------------------------------------------------------------!
      !    Huge negative value, the result can be set to zero.                             !
      !------------------------------------------------------------------------------------!
      eifun8 = 0.d0
      !------------------------------------------------------------------------------------!
   elseif (x <= -1.d0) then
      !------------------------------------------------------------------------------------!
      !     For negative values less than -1.0, we use the polynomial approximation        !
      ! (Equation 5.1.56 of AS72), by taking that Ei(x) = - E1(-x).                        !
      !------------------------------------------------------------------------------------!
      eifun8 = exp(x)/x                                                                    &
             * ( x * ( x * ( x * ( x - apoly(1) ) + apoly(2) ) - apoly(3) ) + apoly(4) )   &
             / ( x * ( x * ( x * ( x - bpoly(1) ) + bpoly(2) ) - bpoly(3) ) + bpoly(4) )
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !    Find both the series expansion and the asymptotic expansion, and pick the one   !
      ! with the lowest absolute value.                                                    !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !  Series expansion: Equation 5.1.11 of AS72, by taking that Ei(x) = - E1(-x).       !
      !                                                                                    !
      !                              Inf                                                   !
      !     Ei(x) = gamma + ln(x) +  SUM u(x,k),                                           !
      !                              k=1                                                   !
      !                                                                                    !
      ! where u(x,k) = [ (-1)^k * (-x)^k / (k * k!) ]                                      !
      !                                                                                    !
      !  To efficiently compute the terms inside the summation, we use that:               !
      !                                                                                    !
      !  u(x,k) = x * (k -1) / k^2 * u(x,k-1), for k >= 2.                                 !
      !------------------------------------------------------------------------------------!
      uxk   = x
      usum  = uxk
      do_expansion: do k = 2, maxiter
         !----- Update the current summation term. ----------------------------------------!
         uxkm1 = uxk
         uxk   = x * dble( k - 1 ) / dble( k * k ) * uxkm1
         !----- Check for degenerate or very large estimate. ------------------------------!
         if ( abs(uxk) > discard8 .or. abs(usum) > discard8) then
            usum = sign(discard8,usum)
            exit do_expansion
         end if
         !----- Check for convergence. ----------------------------------------------------!
         if ( any(abs(uxk) <= [ almost_zero8 * abs(usum), tiny_num8] ) ) exit do_expansion
         !----- Update summation. ---------------------------------------------------------!
         usum = usum + uxk
         !---------------------------------------------------------------------------------!
      end do do_expansion
      !----- Find the series solution. ----------------------------------------------------!
      if ( abs(usum) == discard8) then
         ei_series = sign(discard8,usum)
      else
         ei_series = euler_gam8 + log(abs(x)) + usum
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !  Asymptote expansion: Equation 5.1.51 of AS72 by taking AS72's n=1 and             !
      ! Ei(x) = -E1(-x)).                                                                  !
      !                                                                                    !
      !                          Inf                                                       !
      !     Ei(x) = exp(x) / x * SUM v(x,k),                                               !
      !                          k=0                                                       !
      !                                                                                    !
      ! where v(x,k) = k! / x^k                                                            !
      !                                                                                    !
      !  To efficiently compute the terms inside the summation, we use that:               !
      !                                                                                    !
      !  v(x,k) = k / x * v(x,n -1), for k >= 1.                                           !
      !------------------------------------------------------------------------------------!
      vxk       = 1.d0
      vsum      = vxk
      do_asymptote: do k=1,maxiter
         !----- Update the current summation term. ----------------------------------------!
         vxkm1 = vxk
         vxk   = vxkm1 * dble(k) / x
         !---------------------------------------------------------------------------------!
         !   This method can become degenerate for low x or lead to exceedinly large       !
         ! values, in these cases, halt evaluation.                                        !
         !---------------------------------------------------------------------------------!
         if ( abs(vxkm1) < abs(vxk) .or. abs(vsum) > discard8) then
            vsum = sign(discard8,vsum)
            exit do_asymptote
         end if
         !----- Check for convergence. ----------------------------------------------------!
         if ( any(abs(vxk) <= [ almost_zero8 * abs(vsum), tiny_num8] ) ) exit do_asymptote
         !----- Update summation. ---------------------------------------------------------!
         vsum = vsum + vxk
         !---------------------------------------------------------------------------------!
      end do do_asymptote
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If the solution became degenerate, skip value.                                 !
      !------------------------------------------------------------------------------------!
      if (abs(vsum) == discard8) then
         ei_asymptote = sign(discard8,vsum)
      else
         ei_asymptote = exp(x) * vsum / x
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Pick the lowest absolute value as long as the sign is reasonable.              !
      !------------------------------------------------------------------------------------!
      if (all(abs([ei_series,ei_asymptote]) == discard8)) then
         !----- Huge value, crash because this is iminent over-flow. ----------------------!
         write(unit=*,fmt='(a,1x,es12.5)') 'Attempted X =         ',x
         stop 'Exponential integral cannot be solved for large absolute x.'
         !---------------------------------------------------------------------------------!
      elseif (x < 0.d0) then
         !---------------------------------------------------------------------------------!
         !     Exponential integral is negative when x is negative, however, for some      !
         ! values between -15 < x < -14, the solutions become numerically unstable. Check  !
         ! for the most reasonable estimate.                                               !
         !---------------------------------------------------------------------------------!
         if (ei_series > 0.d0 .and. ei_asymptote > 0.d0) then
            write(unit=*,fmt='(a,1x,es12.5)') 'Attempted X                  = ',x
            write(unit=*,fmt='(a,1x,es12.5)') 'Series expansion estimate    = ',ei_series
            write(unit=*,fmt='(a,1x,es12.5)') 'Asymptote expansion estimate = ',ei_asymptote
            stop 'Exponential integral failed solving, another method might be needed.'
         elseif (ei_series > 0.d0) then
            eifun8 = ei_asymptote
         elseif (ei_asymptote > 0.d0) then
            eifun8 = ei_series
         elseif (abs(ei_series) < abs(ei_asymptote)) then
            eifun8 = ei_series
         else
            eifun8 = ei_asymptote
         end if
         !---------------------------------------------------------------------------------!
      elseif (abs(ei_series) < abs(ei_asymptote)) then
         eifun8 = ei_series
      else
         eifun8 = ei_asymptote
      end if
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   return
end function eifun8
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!    Heapsort is a robust and efficient sorting algorithm introduced by W64. The algorithm !
! implemented here is built from the Wikipedia heapsort pseudocode, which is in turn based !
! on K97.                                                                                  !
!                                                                                          !
! Williams, JWJ (1964). Algorithm 232 - Heapsort, Commun. ACM 7, 347-348.                  !
!    doi:10.1145/512274.512284 (W64).                                                      !
!                                                                                          !
! Knuth, D (1997). The Art of Computer Programming - volume 3: sort and searching.         !
!    section 5.2.3. Sorting by selection (p. 144-155). ISBN 978-0-201-89685-5 (K97).       !
!                                                                                          !
! Wikipedia link: https://en.wikipedia.org/wiki/Heapsort                                   !
!------------------------------------------------------------------------------------------!
subroutine heapsort(nx,xi,increase,xo)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)  :: nx       ! Size of input/output vectors
   real   , dimension(nx), intent(in)  :: xi       ! Input vector
   logical               , intent(in)  :: increase ! Sort from small to large?
   real   , dimension(nx), intent(out) :: xo       ! Output vector
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: i        ! Counter     (inner loop)
   integer                             :: ilwr     ! Lower index (inner loop)
   integer                             :: iupr     ! Upper index (inner loop)
   integer                             :: olwr     ! Lower index (outer loop)
   integer                             :: oupr     ! Upper index (outer loop)
   real                                :: aux      ! Placeholder for element swapping
   !---------------------------------------------------------------------------------------!


   !----- Skip routine in case this has only one element. ---------------------------------!
   if (nx < 2) then
      xo(:) = xi(:)
      return
   else if (.not. increase) then
      !----- Cheat by making numbers negative.  We switch values back in before leaving. --!
      xo(:) = -xi(:)
      !------------------------------------------------------------------------------------!
   else
      xo(:) = xi(:)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Set initial guess of lower index ilwr to half the size of the vector and iupr to   !
   ! the size of the vector. During the heap setting stage, ilwr will be reduced until it  !
   ! becomes 0, and then we start decreasing iupr until it becomes 1, at which point the   !
   ! vector becomes sorted.                                                                !
   !---------------------------------------------------------------------------------------!
   olwr  = nx/2 + 1
   oupr  = nx+1
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Main loop.                                                                        !
   !---------------------------------------------------------------------------------------!
   outer_loop: do
      !------------------------------------------------------------------------------------!
      !     Exit outer loop if we reach the upper bound has already reached 1.             !
      !------------------------------------------------------------------------------------!
      if (oupr == 2) exit outer_loop
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Check whether we are in the heap setting phase or in the retirement-and-promotion   !
      ! phase.                                                                             !
      !------------------------------------------------------------------------------------!
      if (olwr > 1) then
         !----- Heap construction. --------------------------------------------------------!
         olwr = olwr - 1
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Heap extraction.                                                           !
         !---------------------------------------------------------------------------------!
         !----- Shift upper side down one step. -------------------------------------------!
         oupr     = oupr -1
         !----- Swap indices. -------------------------------------------------------------!
         aux      = xo(oupr)
         xo(oupr) = xo(1)
         xo(1)    = aux
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Sift down step.                                                                !
      !------------------------------------------------------------------------------------!
      i = olwr
      inner_loop: do
         !----- Find the lower and right elements. ----------------------------------------!
         ilwr = 2 * i
         iupr = ilwr + 1
         !---------------------------------------------------------------------------------!

         !----- Make sure we do not exceed the heap size. ---------------------------------!
         if (iupr > oupr) exit inner_loop
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Test whether there is an upper element that is larger, and swap the order.  !
         ! Make sure that the elements are bounded before testing vector elements, to      !
         ! avoid segmentation violation.                                                   !
         !---------------------------------------------------------------------------------!
         if (iupr < oupr) then
            if (xo(ilwr) < xo(ilwr+1)) ilwr = ilwr + 1
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Test whether or not to swap elements.                                       !
         !---------------------------------------------------------------------------------!
         if (xo(i) < xo(ilwr)) then
            aux      = xo(i)
            xo(i)    = xo(ilwr)
            xo(ilwr) = aux
            i        = ilwr
         else
            exit inner_loop
         end if
         !---------------------------------------------------------------------------------!
      end do inner_loop
      !------------------------------------------------------------------------------------!
   end do outer_loop
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Before we leave, check whether this should be a high-to-low sorting.  In case so, !
   ! switch the sign again.                                                                !
   !---------------------------------------------------------------------------------------!
   if (.not. increase) xo(:) = -xo(:)
   !---------------------------------------------------------------------------------------!

   return
end subroutine heapsort
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     Function that defines the quantile given a vector.  This is a rather simple          !
! estimator, it should work reasonably well as long as x is sufficiently large and not a   !
! crazy distribution.                                                                      !
!------------------------------------------------------------------------------------------!
real function fquant(nx,x,prob)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in) :: nx
   real   , dimension(nx), intent(in) :: x
   real                  , intent(in) :: prob
   !----- Internal variables. -------------------------------------------------------------!
   real   , dimension(nx)             :: xsort
   integer                            :: il
   integer                            :: ih
   real                               :: wl
   real                               :: wh
   real                               :: ridx
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Sanity check: prob must be between 0 and 1.  If not, crash!                       !
   !---------------------------------------------------------------------------------------!
   if (prob < 0. .or. prob > 1.) then
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) '================================================='
      write(unit=*,fmt='(a)'          ) '================================================='
      write(unit=*,fmt='(a)'          ) '    In function fquant: Invalid PROB!'
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a,1x,es12.5)') '   -> Provided PROB: ',prob
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) ' '
      call fatal_error('Invalid prob setting','fquant','numutils.f90')
   end if
   !---------------------------------------------------------------------------------------!


   !----- Sort output vector. -------------------------------------------------------------!
   call heapsort(nx,x,.true.,xsort)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Find the quantile position in terms of indices.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Position without interpolation. -------------------------------------------------!
   ridx   = 1. + prob * real(nx-1)
   !----- Index just before ridx. ---------------------------------------------------------!
   il     = max(1,floor(ridx))
   !----- Index just after ridx. ----------------------------------------------------------!
   ih     = min(nx,ceiling(ridx))
   !----- Quantile is the interpolated value. ---------------------------------------------!
   if (il == ih) then
      fquant = xsort(il)
   else
      !----- Weight factors. --------------------------------------------------------------!
      wl     = ridx - real(il)
      wh     = real(ih) - ridx
      !------------------------------------------------------------------------------------!

      !----- Quantile is the weighted average. --------------------------------------------!
      fquant = (wl * xsort(il) + wh * xsort(ih)) / (wl + wh)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   return
end function fquant
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     Wrapper for the function above, with an additional mask vector to select only some   !
! elements.                                                                                !
!------------------------------------------------------------------------------------------!
real function fquant_mask(nx,x,mask,prob)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)  :: nx
   real   , dimension(nx), intent(in)  :: x
   logical, dimension(nx), intent(in)  :: mask
   real                  , intent(in)  :: prob
   !----- Internal variables. -------------------------------------------------------------!
   real   , dimension(:) , allocatable :: xuse
   integer                             :: nuse
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: fquant
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Count elements to be used, then allocate xuse and xsort.                         !
   !---------------------------------------------------------------------------------------!
   nuse        = count(mask)
   allocate (xuse(nuse))
   xuse        = pack(x,mask)
   fquant_mask = fquant(nuse,xuse,prob)
   deallocate(xuse)
   !---------------------------------------------------------------------------------------!

   return
end function fquant_mask
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     Sub-routine that solves the quadratic equation ( a * x**2 + b * x + c = 0).          !
! We test whether or not this is a trivial case that does not require solving the full     !
! equation. For the full equation, we use the approach by H02 to avoid floating point      !
! issues when solving roots. We further check whether or not the discriminant is negative. !
!                                                                                          !
!     The subroutine also requires a "undef" flag to be passed, which will flag cases      !
! in which one or both solutions are not valid. This is an argument so the solver can be   !
! used when either the largest or the smallest root is sought.                             !
!                                                                                          !
! Higham, N. J., 2002: Accuracy and Stability of Numerical Algorithms. 2nd ed., Society    !
!    for Industrial and Applied Mathematics, Philadelphia, PA, United States,              !
!    doi:10.1137/1.9780898718027 (H02).                                                    !
!------------------------------------------------------------------------------------------!
subroutine solve_quadratic(aquad,bquad,cquad,undef,root1,root2)
   use rconstants, only : tiny_num
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in)  :: aquad
   real(kind=4), intent(in)  :: bquad
   real(kind=4), intent(in)  :: cquad
   real(kind=4), intent(in)  :: undef
   real(kind=4), intent(out) :: root1
   real(kind=4), intent(out) :: root2
   !----- Internal variables. -------------------------------------------------------------!
   real(kind=4)              :: discr
   logical                   :: a_offzero
   logical                   :: b_offzero
   logical                   :: c_offzero
   !---------------------------------------------------------------------------------------!



   !----- Save logical tests. -------------------------------------------------------------!
   a_offzero = abs(aquad) >= tiny_num
   b_offzero = abs(bquad) >= tiny_num
   c_offzero = abs(cquad) >= tiny_num
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check for cases to solve.                                                         !
   !---------------------------------------------------------------------------------------!
   if (a_offzero .and. ( b_offzero .or. c_offzero ) ) then
      !------------------------------------------------------------------------------------!
      !    Quadratic equation with two non-zero solutions. Find the discriminant to find   !
      ! out whether the solutions are real (if negative, then the roots are complex).      !
      !------------------------------------------------------------------------------------!
      discr = bquad*bquad - 4.0 * aquad * cquad
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Check discriminant sign (but allow for round-off errors).                      !
      !------------------------------------------------------------------------------------!
      if (discr >= - tiny_num) then
         !----- Coerce discriminant to non-negative. --------------------------------------!
         discr = max(0.0,discr)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Follow H02's approach to find the largest root (absolute value) from the    !
         ! traditional quadratic equation, then derive the second root from the first one. !
         ! This is safe whenever b or c are non-zero.                                      !
         !---------------------------------------------------------------------------------!
         root1  = - (bquad + sign(sqrt(discr),bquad)) / ( 2. * aquad )
         root2  = cquad / ( aquad * root1 )
         !---------------------------------------------------------------------------------!
      else
         !----- Negative discriminant, return invalid roots. ------------------------------!
         root1  = undef
         root2  = undef
         !---------------------------------------------------------------------------------!
      end if
   else if (a_offzero) then
      !------------------------------------------------------------------------------------!
      !     Both bquad and cquad are nearly zero. Double root, and both have to be zero.   !
      !------------------------------------------------------------------------------------!
      root1 = 0.0
      root2 = 0.0
      !------------------------------------------------------------------------------------!
   else if (b_offzero) then
      !------------------------------------------------------------------------------------!
      !     "aquad" is not zero, not a true quadratic equation. Single root.               !
      !------------------------------------------------------------------------------------!
      root1 = - cquad / bquad
      root2 = undef
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !     Both aquad and bquad are zero, this really doesn't make any sense and should   !
      ! never happen. If it does, issue an error and stop the run.                         !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a)')           ' Quadratic equation cannot be solved!'
      write (unit=*,fmt='(a)')           ' ''aquad'' and/or ''bquad'' must be non-zero.'
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a,1x,es12.5)') ' aquad = ',aquad
      write (unit=*,fmt='(a,1x,es12.5)') ' bquad = ',bquad
      write (unit=*,fmt='(a,1x,es12.5)') ' cquad = ',cquad
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      call fatal_error(' Invalid coefficients for quadratic equation'                      &
                      ,'solve_quadratic','numutils.f90')
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine solve_quadratic
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     Sub-routine that solves the quadratic equation ( a * x**2 + b * x + c = 0).          !
! We test whether or not this is a trivial case that does not require solving the full     !
! equation. For the full equation, we use the approach by H02 to avoid floating point      !
! issues when solving roots. We further check whether or not the discriminant is negative. !
!                                                                                          !
!     The subroutine also requires a "undef" flag to be passed, which will flag cases      !
! in which one or both solutions are not valid. This is an argument so the solver can be   !
! used when either the largest or the smallest root is sought.                             !
!                                                                                          !
! Higham, N. J., 2002: Accuracy and Stability of Numerical Algorithms. 2nd ed., Society    !
!    for Industrial and Applied Mathematics, Philadelphia, PA, United States,              !
!    doi:10.1137/1.9780898718027 (H02).                                                    !
!------------------------------------------------------------------------------------------!
subroutine solve_quadratic8(aquad,bquad,cquad,undef,root1,root2)
   use rconstants, only : tiny_num8
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in)  :: aquad
   real(kind=8), intent(in)  :: bquad
   real(kind=8), intent(in)  :: cquad
   real(kind=8), intent(in)  :: undef
   real(kind=8), intent(out) :: root1
   real(kind=8), intent(out) :: root2
   !----- Internal variables. -------------------------------------------------------------!
   real(kind=8)              :: discr
   logical                   :: a_offzero
   logical                   :: b_offzero
   logical                   :: c_offzero
   !---------------------------------------------------------------------------------------!



   !----- Save logical tests. -------------------------------------------------------------!
   a_offzero = abs(aquad) >= tiny_num8
   b_offzero = abs(bquad) >= tiny_num8
   c_offzero = abs(cquad) >= tiny_num8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check for cases to solve.                                                         !
   !---------------------------------------------------------------------------------------!
   if (a_offzero .and. ( b_offzero .or. c_offzero ) ) then
      !------------------------------------------------------------------------------------!
      !    Quadratic equation with two non-zero solutions. Find the discriminant to find   !
      ! out whether the solutions are real (if negative, then the roots are complex).      !
      !------------------------------------------------------------------------------------!
      discr = bquad*bquad - 4.d0 * aquad * cquad
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Check discriminant sign (but allow for round-off errors).                      !
      !------------------------------------------------------------------------------------!
      if (discr >= - tiny_num8) then
         !----- Coerce discriminant to non-negative. --------------------------------------!
         discr = max(0.d0,discr)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Follow H02's approach to find the largest root (absolute value) from the    !
         ! traditional quadratic equation, then derive the second root from the first one. !
         ! This is safe whenever b or c are non-zero.                                      !
         !---------------------------------------------------------------------------------!
         root1  = - (bquad + sign(sqrt(discr),bquad)) / ( 2.d0 * aquad )
         root2  = cquad / ( aquad * root1 )
         !---------------------------------------------------------------------------------!
      else
         !----- Negative discriminant, return invalid roots. ------------------------------!
         root1  = undef
         root2  = undef
         !---------------------------------------------------------------------------------!
      end if
   else if (a_offzero) then
      !------------------------------------------------------------------------------------!
      !     Both bquad and cquad are nearly zero. Double root, and both have to be zero.   !
      !------------------------------------------------------------------------------------!
      root1 = 0.d0
      root2 = 0.d0
      !------------------------------------------------------------------------------------!
   else if (b_offzero) then
      !------------------------------------------------------------------------------------!
      !     "aquad" is not zero, not a true quadratic equation. Single root.               !
      !------------------------------------------------------------------------------------!
      root1 = - cquad / bquad
      root2 = undef
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !     Both aquad and bquad are zero, this really doesn't make any sense and should   !
      ! never happen. If it does, issue an error and stop the run.                         !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a)')           ' Quadratic equation cannot be solved!'
      write (unit=*,fmt='(a)')           ' ''aquad'' and/or ''bquad'' must be non-zero.'
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a,1x,es12.5)') ' aquad = ',aquad
      write (unit=*,fmt='(a,1x,es12.5)') ' bquad = ',bquad
      write (unit=*,fmt='(a,1x,es12.5)') ' cquad = ',cquad
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      call fatal_error(' Invalid coefficients for quadratic equation'                      &
                      ,'solve_quadratic8','numutils.f90')
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine solve_quadratic8
!==========================================================================================!
!==========================================================================================!

   return
end subroutine solve_quadratic8
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine extracts a vertical (z) column given a 3-D array, and the fixed      !
! indices for the x and y dimensions.                                                      !
!------------------------------------------------------------------------------------------!
subroutine array2zcol(mz,mx,my,x,y,array,vector)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)  :: mz
   integer                          , intent(in)  :: mx
   integer                          , intent(in)  :: my
   integer                          , intent(in)  :: x
   integer                          , intent(in)  :: y
   real(kind=4), dimension(mz,mx,my), intent(in)  :: array
   real(kind=4), dimension(mz)      , intent(out) :: vector
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: z
   !---------------------------------------------------------------------------------------!

   do z=1,mz
      vector(z) = array(z,x,y)
   end do

   return
end subroutine array2zcol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine extracts a longitudinal (x) column given a 3-D array, and the fixed  !
! indices for the z and y dimensions.                                                      !
!------------------------------------------------------------------------------------------!
subroutine array2xcol(mz,mx,my,z,y,array,vector)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)  :: mz
   integer                          , intent(in)  :: mx
   integer                          , intent(in)  :: my
   integer                          , intent(in)  :: z
   integer                          , intent(in)  :: y
   real(kind=4), dimension(mz,mx,my), intent(in)  :: array
   real(kind=4), dimension(mx)      , intent(out) :: vector
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: x
   !---------------------------------------------------------------------------------------!

   do x=1,mx
      vector(x) = array(z,x,y)
   end do

   return
end subroutine array2xcol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine extracts a latitudinal (y) column given a 3-D array, and the fixed   !
! indices for the z and x dimensions.                                                      !
!------------------------------------------------------------------------------------------!
subroutine array2ycol(mz,mx,my,z,x,array,vector)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)  :: mz
   integer                          , intent(in)  :: mx
   integer                          , intent(in)  :: my
   integer                          , intent(in)  :: z
   integer                          , intent(in)  :: x
   real(kind=4), dimension(mz,mx,my), intent(in)  :: array
   real(kind=4), dimension(my)      , intent(out) :: vector
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: y
   !---------------------------------------------------------------------------------------!

   do y=1,my
      vector(y) = array(z,x,y)
   end do

   return
end subroutine array2ycol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies a vertical (z) column to a 3-D array, using fixed indices for !
! the x and y dimensions.                                                                  !
!------------------------------------------------------------------------------------------!
subroutine zcol2array(mz,mx,my,x,y,vector,array)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: mz
   integer                          , intent(in)    :: mx
   integer                          , intent(in)    :: my
   integer                          , intent(in)    :: x
   integer                          , intent(in)    :: y
   real(kind=4), dimension(mz)      , intent(in)    :: vector
   real(kind=4), dimension(mz,mx,my), intent(inout) :: array
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: z
   !---------------------------------------------------------------------------------------!

   do z=1,mz
      array(z,x,y) = vector(z)
   end do

   return
end subroutine zcol2array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies a longitudinal (x) column to a 3-D array, using fixed indices !
! for the z and y dimensions.                                                              !
!------------------------------------------------------------------------------------------!
subroutine xcol2array(mz,mx,my,z,y,vector,array)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: mz
   integer                          , intent(in)    :: mx
   integer                          , intent(in)    :: my
   integer                          , intent(in)    :: z
   integer                          , intent(in)    :: y
   real(kind=4), dimension(mx)      , intent(in)    :: vector
   real(kind=4), dimension(mz,mx,my), intent(inout) :: array
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: x
   !---------------------------------------------------------------------------------------!

   do x=1,mx
      array(z,x,y) = vector(x)
   end do

   return
end subroutine xcol2array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies a latitudinal (y) column to a 3-D array, using fixed indices  !
! for the z and x dimensions.                                                              !
!------------------------------------------------------------------------------------------!
subroutine ycol2array(mz,mx,my,z,x,vector,array)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: mz
   integer                          , intent(in)    :: mx
   integer                          , intent(in)    :: my
   integer                          , intent(in)    :: z
   integer                          , intent(in)    :: x
   real(kind=4), dimension(my)      , intent(in)    :: vector
   real(kind=4), dimension(mz,mx,my), intent(inout) :: array
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: y
   !---------------------------------------------------------------------------------------!

   do y=1,my
      array(z,x,y) = vector(y)
   end do

   return
end subroutine ycol2array
!==========================================================================================!
!==========================================================================================!

