!===================================== Change Log =========================================!
! 2.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!                      SYSTEM DEPENDENT ROUTINES                                           !
!                                                                                          !
!    This module contains short utility routines that are not of the FORTRAN 77 standard   !
! and may differ from system to system.  These include bit manipulation, I/O, JCL calls,   !
! and vector functions.                                                                    !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Routine to get command line argument.                                                !
!------------------------------------------------------------------------------------------!
subroutine ugetarg(i,arg)
   implicit none
   integer          :: i
   character(len=*) :: arg
       

#if defined(HP)
   call getarg(i+1,arg)
#else
   call getarg(i,arg)
#endif

   return
end subroutine ugetarg
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
!     Routine returns CPU time.  Called with ICALL=1 at beginning of timestep, ICALL=2 at  !
! end of timestep.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine timing(icall,t1)
   implicit none
   integer, intent(in)   :: icall
   real   , intent(out)  :: t1
   real   , dimension(2) :: et(2)

#if defined(VAX)
   integer :: iad0
#endif

#if defined(IBM)
   real   , external     :: mclock
#elif defined(CRAY)
   real   , external     :: cputime
#elif defined(__APPLE__)
   real                  :: etime
#elif defined(PC_GFORTRAN)
   real                  :: etime
#elif defined(__GFORTRAN__)
   real                  :: etime
#else
   real   ,external      :: etime
#endif

   select case (icall)
   !----- Start call. ---------------------------------------------------------------------!
   case (1)
#if defined(CRAY)
      call cpu_time(T1)
#elif defined(VAX)
      iad0=0
      call lib$init_timer(iad0)
#elif defined(IBM)
      T1=mclock(et)/100.
#else
      T1=ETIME(et)
#endif

   !----- End call. -----------------------------------------------------------------------!
   case (2)
#if defined(VAX)
      call LIB$SHOW_TIMER(IAD0,2)
#elif defined(CRAY)
      call cpu_time(T1)
#elif defined(IBM)
      T1=mclock(et)/100.
#else
      T1=ETIME(et)
#endif

   end select

   return
end subroutine timing
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     Function that checks whether a number is NaN.  This works with PGI, Intel, and GNU.  !
!------------------------------------------------------------------------------------------!
logical function isnan_real(x)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   real, intent(in) :: x
   !---------------------------------------------------------------------------------------!

#if defined(PGI)
   isnan_real = x /= x
#else
   isnan_real = isnan(x)
#endif

   return
end function isnan_real
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     Function that checks whether a number is NaN.  This works with PGI, Intel, and GNU.  !
!------------------------------------------------------------------------------------------!
logical function isnan_dble(x)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   !---------------------------------------------------------------------------------------!

#if defined(PGI)
   isnan_dble = x /= x
#else
   isnan_dble = isnan(x)
#endif

   return
end function isnan_dble
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     Function that waits for a few seconds before moving on.  This may be useful for OMP  !
! instructions.  Most fortran distributions should work fine with built-in sleep, but it   !
! may require libraries.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine wait_sec(x)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer, intent(in) :: x
   !---------------------------------------------------------------------------------------!

   call sleep(x)

   return
end subroutine wait_sec
!==========================================================================================!
!==========================================================================================!
