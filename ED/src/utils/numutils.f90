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


