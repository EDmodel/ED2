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

!******************************************************************
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
   if (find_rank < 0) call fatal_error('Index not found','find_rank','numutils.f90')
   return
end function find_rank
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     This sub-routine orders the elements of a real vector given the index vector.        !
!------------------------------------------------------------------------------------------!
subroutine order_real(nmax,variable,idx)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nmax
   real   , dimension(nmax) , intent(inout) :: variable
   integer, dimension(nmax) , intent(in)    :: idx
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: n
   real   , dimension(nmax)                 :: vtemp
   !---------------------------------------------------------------------------------------!

   vtemp(:)  = variable(:)
   do n=1,nmax
      variable(idx(n)) = vtemp(n)
   end do

   return
end subroutine order_real
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     This sub-routine orders the elements of a integer vector given the index vector.     !
!------------------------------------------------------------------------------------------!
subroutine order_integer(nmax,variable,idx)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nmax
   integer, dimension(nmax) , intent(inout) :: variable
   integer, dimension(nmax) , intent(in)    :: idx
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: n
   integer, dimension(nmax)                 :: vtemp
   !---------------------------------------------------------------------------------------!

   vtemp(:)  = variable(:)
   do n=1,nmax
      variable(idx(n)) = vtemp(n)
   end do

   return
end subroutine order_integer
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!   This function simply computes the cube root of all numbers, including the negative     !
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
   use consts_coms, only : euler_gam8   & ! intent(in)
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
   use consts_coms, only : tiny_num
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
   use consts_coms, only : tiny_num8
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
