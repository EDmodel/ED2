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
!                        0                                                                 !
!                                                                                          !
!     This function is based on:                                                           !
!                                                                                          !
! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery: 1992. Numerical recipes !
!    in Fortran 77.  Cambridge University Press, section 6.3 p. 215-219.                   !
!                                                                                          !
! with the difference that we also solve for negative numbers.  Zero cannot be solved, so  !
! if this happens, or if the sought number would lead to infinity, we stop the model.      !
!------------------------------------------------------------------------------------------!
real(kind=8) function eifun8(x)
   use consts_coms, only : euler_gam8 & ! intent(in)
                         , lnexp_min8 & ! intent(in)
                         , lnexp_max8 & ! intent(in)
                         , tiny_num8  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: sum
   real(kind=8)             :: term
   real(kind=8)             :: fact
   real(kind=8)             :: prev
   real(kind=8)             :: diter
   integer                  :: iter
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter  :: powerlim  = 1.5d+01
   real(kind=8), parameter  :: converge  = 1.0d-7
   integer     , parameter  :: maxiter   = 100
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check what to do depending on the value of x.                                     !
   !---------------------------------------------------------------------------------------!
   if (x == 0.d0) then
      !------------------------------------------------------------------------------------!
      !     Zero.  This is a singularity and the user should never call it in this case.   !
      ! That's sad, but we ought to quit this run and tell the user why the run crashed.   !
      !------------------------------------------------------------------------------------!
      call fatal_error('Exponential integral cannot be solved for x = 0.'                  &
                      ,'eifun8','numutils.f90')
   elseif (x >= lnexp_max8) then
      !----- Huge value, crash because this is iminent over-flow. -------------------------!
      write(unit=*,fmt='(a,1x,es12.5)') 'Attempted X =         ',x
      write(unit=*,fmt='(a,1x,es12.5)') 'Maximum acceptable X =',lnexp_max8
      call fatal_error('Exponential integral cannot be solved for x = 0.'                  &
                      ,'eifun8','numutils.f90')
   elseif (abs(x) <= lnexp_min8) then
      !----- Huge negative number, the result can be rounded to zero. ---------------------!
      eifun8 = 0.d0
   elseif (abs(x) <= tiny_num8) then
      !----- The number is too close to zero, bypass iterative methods. -------------------!
      eifun8 = euler_gam8 + log(abs(x))
   elseif (abs(x) <= powerlim) then
      !------------------------------------------------------------------------------------!
      !    Input x is small, so we use the power method.                                   !
      !------------------------------------------------------------------------------------!
      fact      = 1.d0
      sum       = 0.d0
      powerloop: do iter=1,maxiter
         diter = dble(iter)
         fact  = fact * x / diter
         term  = fact / diter
         sum   = sum + term
         !----- If the term is tiny, we have reached convergence, quit the loop. ----------!
         if (abs(term) < converge * abs(sum)) exit powerloop
      end do powerloop
      eifun8   = euler_gam8 + log(abs(x)) + sum
   else
      !------------------------------------------------------------------------------------!
      !    Input x is large, so we use the asymptotic approximation.                       !
      !------------------------------------------------------------------------------------!
      sum       = 0.d0
      term      = 1.d0
      asymploop: do iter=1,maxiter
         diter = dble(iter)
         prev  = term
         term  = term * diter / x
         if (abs(term) < converge) then
            !----- The term is tiny, we have reached convergence, quit the loop. ----------!
            exit asymploop
         elseif (abs(term) >= abs(prev)) then
            !------------------------------------------------------------------------------!
            !   Series is diverging, we are probably reaching round-off errors, we better  !
            ! stop now.                                                                    !
            !------------------------------------------------------------------------------!
            sum = sum - prev
            exit asymploop
         else
            sum = sum + term
         end if
      end do asymploop
      eifun8 = exp(x) * (1.d0 + sum) / x
   end if

   return
end function eifun8
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!    Heapsort is a robust and efficient sorting algorithm.  For more details, check        !
! The Numerical Recipes Book (chapter 8).                                                  !
!------------------------------------------------------------------------------------------!
subroutine heapsort(nx,xi,increase,xo)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer               , intent(in)  :: nx       ! Size of input/output vectors
   real   , dimension(nx), intent(in)  :: xi       ! Input vector
   logical               , intent(in)  :: increase ! Sort from small to large?
   real   , dimension(nx), intent(out) :: xo       ! Output vector
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: i        ! Counter
   integer                             :: ir       ! Index of selected data
   integer                             :: j        ! Index of selected data
   integer                             :: l        ! Index of selected data
   real                                :: aux      ! Placeholder
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
   !    The index l will be decremented from its initial value down to 1 during the        !
   ! "hiring" (heap creation) phase.  Once it reaches 1, the index ir will be decremented  !
   ! from its initial value down to 1 during the "retirement-and-promotion" (heap          !
   ! selection) phase.                                                                     !
   !---------------------------------------------------------------------------------------!
   l  = nx/2 + 1
   ir = nx
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Main loop.                                                                        !
   !---------------------------------------------------------------------------------------!
   outer_loop: do
      !------------------------------------------------------------------------------------!
      !      Check whether we are in the hiring phase or in the retirement-and-promotion   !
      ! phase.                                                                             !
      !------------------------------------------------------------------------------------!
      if (l > 1) then
         !---------------------------------------------------------------------------------!
         !     Still in hiring phase.                                                      !
         !---------------------------------------------------------------------------------!
         l   = l - 1
         aux = xo(l)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    In the retirement-and-promotion phase.                                       !
         !---------------------------------------------------------------------------------!
         !----- Clear a space at end of array. --------------------------------------------!
         aux    = xo(ir)
         !----- Retire the top of the heap into it. ---------------------------------------!
         xo(ir) = xo(1)
         !----- Decrease the size of the corporation. -------------------------------------!
         ir      = ir -1
         !----- Check how we are doing with promotions. -----------------------------------!
         if (ir == 1) then
            !----- Done with the last promotion.  The least competent worker of all! ------!
            xo(1) = aux
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Time to leave the loop (and the sub-routine).                           !
            !------------------------------------------------------------------------------!
            exit outer_loop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Whether in hiring phase or promotion phase, we here set up to sift down element !
      ! aux to its proper level.                                                           !
      !------------------------------------------------------------------------------------!
      i = l
      j = l+1
      inner_loop: do
         if (j > ir) exit inner_loop

         !----- Compare to the better underling. ------------------------------------------!
         if (j < ir) then
            if(xo(j) < xo(j+1)) j = j + 1
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Check whether to demote aux or not.                                         !
         !---------------------------------------------------------------------------------!
         if (aux < xo(j)) then
            !----- Demote aux. ------------------------------------------------------------!
            xo(i) = xo(j)
            i     = j
            j     = j + j
            !------------------------------------------------------------------------------!
         else
            !----- This is aux's level.  Set j to terminate the sift-down. ----------------!
            j     = ir + 1
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do inner_loop
      !------------------------------------------------------------------------------------!

      !----- Put aux into its slot. -------------------------------------------------------!
      xo(i) = aux
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
! FUNCTION bpow01
!\brief Safe power estimate to avoid floating point exceptions
!\author Marcos Longo 3 March 2021
!\details This function to calculate power functions for numbers bounded between 0 and 1
!!        safely.  It uses that y = x ** a = exp(a * ln(x)), and use the safe log limits
!!        to avoid FPE errors.  This function "should" work for values greater than 1 too,
!!        but don't use if for negative numbers.
!------------------------------------------------------------------------------------------!
real(kind=4) function bpow01(x,a)
   use consts_coms, only : lnexp_min & ! intent(in)
                         , lnexp_max ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in) :: x
   real(kind=4), intent(in) :: a
   !----- Internal variables. -------------------------------------------------------------!
   real(kind=4)             :: lnexp
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Check if x is greater than zero (if it is exactly zero, we cannot use the log    !
   ! approach).                                                                            !
   !---------------------------------------------------------------------------------------!
   if (x > 0.) then


      !----- Find the bounded term inside the exponential. --------------------------------!
      lnexp = max(lnexp_min,min(lnexp_max,a * log(x)))
      !------------------------------------------------------------------------------------!


      !----- Report result. ---------------------------------------------------------------!
      bpow01 = exp(lnexp)
      !------------------------------------------------------------------------------------!
   else if (x == 0. .and. a > 0.) then
      !----- By definition 0^a = 0 as long as a > 0. --------------------------------------!
      bpow01 = 0.
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !    Invalid input data, stop everything.                                            !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt='(a)'          ) '-----------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' Invalid variables for bpow01!'
      write(unit=*,fmt='(a)'          ) '-----------------------------------------------'
      write(unit=*,fmt='(a,1x,es12.5)') ' x = ',x
      write(unit=*,fmt='(a,1x,es12.5)') ' a = ',a
      write(unit=*,fmt='(a)'          ) '-----------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' If x >= 0., then make sure a >= 0. '
      write(unit=*,fmt='(a)'          ) ' If x = 0., then make sure a > 0. '
      write(unit=*,fmt='(a)'          ) ' Negative x values are not allowed. '
      write(unit=*,fmt='(a)'          ) '-----------------------------------------------'
      call fatal_error('Invalid values for x and/or a.','bpow01','numutils.f90')
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   return
end function bpow01
!==========================================================================================!
!==========================================================================================!

