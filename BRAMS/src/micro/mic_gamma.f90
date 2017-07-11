!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function gammp(a,x)
   implicit none
   real :: a,x,gln,gammcf

   if(x < a+1.) then
       call gser(gammp,a,x,gln)
   else
       call gcf(gammcf,a,x,gln)
       gammp = 1. - gammcf
   endif

   return
end function gammp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function gammq(a,x)
   implicit none
   real :: a,x,gamser,gln

   if(x < a+1.) then
       call gser(gamser,a,x,gln)
       gammq = 1. - gamser
   else
       call gcf(gammq,a,x,gln)
   endif

   return
end function gammq
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine gcf(gammcf,a,x,gln)
   use therm_lib, only : toler
   implicit none
   real :: a,x,gammcf,gln

   real   , parameter :: tinypow=-38. ! To avoid FPE errors
   real   , external  :: gammln
   integer, parameter :: maxit_loc = 15

   real :: gold,a0,a1,b0,b1,fac,an,ana,anf,gaccel
   integer :: n

   gln = gammln(a)
   gold = 0.
   a0 = 1.
   a1 = x
   b0 = 0.
   b1 = 1.
   fac = 1.
   itloop: do n = 1, maxit_loc
      an = real(n)
      ana = an - a
      a0 = (a1+a0*ana)*fac
      b0 = (b1+b0*ana)*fac
      anf = an*fac
      a1 = x*a0 + anf*a1
      b1 = x*b0 + anf*b1
      if (a1 /= 0.) then
          fac = 1./a1
          gaccel = b1*fac
          if(abs(gaccel-gold) < toler * gaccel) exit itloop
          gold = gaccel
      end if
   end do itloop

   fac=-x+a*alog(x)-gln
   if(fac > tinypow) then
     gammcf = exp(fac)*gaccel
   else
     gammcf = 0.
   end if
   return
end subroutine gcf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine gser(gamser,a,x,gln)
   use therm_lib , only : maxit      & ! intent(in)
                        , toler      ! ! intent(in)
   use rconstants, only : lnexp_min  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in)  :: a
   real   , intent(in)  :: x
   real   , intent(out) :: gln
   real   , intent(out) :: gamser
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: ap
   real                 :: sum
   real                 :: del
   real                 :: lnexp
   integer              :: n
   !----- Local constants. ----------------------------------------------------------------!
   real   , parameter   :: tiny_offset = 1.e-30
   !----- External functions. -------------------------------------------------------------!
   real   , external    :: gammln
   !---------------------------------------------------------------------------------------!


   gln = gammln(a)
   if (x < 0.) then
      write (unit=*,fmt='(a,1x,es13.6)') 'A   = ',a
      write (unit=*,fmt='(a,1x,es13.6)') 'X   = ',x
      write (unit=*,fmt='(a,1x,es13.6)') 'GLN = ',gln
      call abort_run('Invalid X, it should be non-negative.','gser','mic_gamma.f90')
   elseif (x < tiny_offset) then
      gamser = 0.
      return
   end if

   ap = a
   sum = 1./a
   del = sum

   itloop: do n = 1, maxit
      ap = ap + 1.
      del = del*x/ap
      sum =  sum  +  del
      if(abs(del) < tiny_offset .or. abs(del) < abs(sum)*toler) exit itloop
   end do itloop

   lnexp = -x + a*log(x) - gln

   if(lnexp > lnexp_min) then
     gamser = sum * exp(lnexp)
   else
     gamser = 0.
   endif
   return
end subroutine gser
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function gammln(xx)
  implicit none
  real :: xx

  real(kind=8), dimension(6), parameter :: &
      cof=(/76.18009173d0, -86.50532033d0, 24.01409822d0,  &
           -1.231739516d0, .120858003d-2, -.536382d-5/)
  real(kind=8), parameter :: stp=2.50662827465d0, half=0.5d0, one=1.0d0, fpf=5.5d0

  real(kind=8) :: x,tmp,ser
  integer :: j

  x=dble(xx)-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
      x=x+one
      if (x < 1.d-10) then
        x=sign(1.d-10,x)
      end if
      ser=ser+cof(j)/x
  enddo
  gammln=sngl(tmp+dlog(stp*ser))
  return
end function gammln
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine avint(x,y,n,xlo,xup,ans)
implicit none
integer :: n
real :: x(n),y(n),xlo,xup,ans

real(kind=8) ::r3,rp5,sum,syl,syl2,syl3,syu,syu2,syu3,x1,x2,x3  &
,x12,x13,x23,term1,term2,term3,a,b,c,ca,cb,cc

integer :: i,inlft,inrt,istart,istop
real :: slope,fl,fr

ans =0.0
if(xlo.lt.xup) goto 3
if(xlo.eq.xup) goto 100
if(xlo.gt.xup) goto 200
3 if(n.lt.2) goto 215
do i=2,n
   if(x(i).le.x(i-1)) goto 210
   if(x(i).gt.xup) goto 6
enddo
6 continue
if(n.ge.3) goto 9

!     special n=2 case

slope = (y(2)-y(1))/(x(2)-x(1))
fl = y(1) + slope*(xlo-x(1))
fr = y(2) + slope*(xup-x(2))
ans = 0.5*(fl+fr)*(xup-xlo)
return
9 continue
if(x(n-2).lt.xlo)  goto 205
if(x(3).gt.xup)    goto 205
i = 1
10 if(x(i).ge.xlo) goto 15
i = i+1
goto 10
15 inlft = i
i = n
20 if(x(i).le.xup) goto 25
i = i-1
goto 20
25 inrt = i
if((inrt-inlft).lt.2) goto 205
istart = inlft
if(inlft.eq.1) istart = 2
istop  = inrt
if(inrt.eq.n)  istop  = n-1

r3 = 3.0d0
rp5= 0.5d0
sum = 0.0
syl = xlo
syl2= syl*syl
syl3= syl2*syl

do i=istart,istop
   x1 = x(i-1)
   x2 = x(i)
   x3 = x(i+1)
   x12 = x1-x2
   x13 = x1-x3
   x23 = x2-x3
   term1 = dble(y(i-1))/(x12*x13)
   term2 =-dble(y(i)) /(x12*x23)
   term3 = dble(y(i+1))/(x13*x23)
   a = term1+term2+term3
   b = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3
   c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
   if(i.le.istart) goto 30
   if(i.gt.istart) goto 35
30    ca = a
   cb = b
   cc = c
   goto 40
35    ca = 0.5*(a+ca)
   cb = 0.5*(b+cb)
   cc = 0.5*(c+cc)
40    syu = x2
   syu2= syu*syu
   syu3= syu2*syu
   sum = sum + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2)  &
             + cc*(syu-syl)
   ca  = a
   cb  = b
   cc  = c
   syl = syu
   syl2= syu2
   syl3= syu3
enddo
syu = xup
ans = sum + ca*(syu**3-syl3)/r3 + cb*rp5*(syu**2-syl2)  &
          + cc*(syu-syl)
100 return
200 print*, 'Upper limit of integration not greater than lower limit.'
stop 'avint2'
205 print*, 'Less than 3 function values between integration limits.'
stop 'avint3'
210 print*, 'Abscissas not strictly increasing.'
stop 'avint4'
215 print*, 'Less than 2 function values were supplied.'
stop 'avint5'
end
