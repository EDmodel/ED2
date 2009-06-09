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
!     Returns the endian 'type' of machine.                                                !
!------------------------------------------------------------------------------------------!
subroutine endian(mach_type)
   implicit none
   character(len=*) :: mach_type

#if defined(ALPHA) || defined(PC_NT1)
   mach_type='little_endian'
#else
   mach_type='big_endian'
#endif

   return
end subroutine endian
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function returns a factor for the length of a random access record length unit. !
! For example, IBM specifies bytes while SGI uses words, so specify the open statement in  !
! WORDS, then multiply by this returned value.                                             !
!------------------------------------------------------------------------------------------!
integer function iran_recsize()
   implicit none
      
#if defined(ALPHA) || defined(CRAY) 
   iran_recsize = 1
#else
   iran_recsize = 4
#endif

   return
end function iran_recsize
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
#elif defined(MAC_OS_X)
   real                  :: etime
#else
   real   , external     :: etime
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
!    This subroutine swaps the order of two bytes in integers of kind=2.                   !
!------------------------------------------------------------------------------------------!
subroutine dcw_swap16 (a,n)
   implicit none

   integer                       , intent(in)    :: n
   integer(kind=2) , dimension(n), intent(inout) :: a

#if defined(SGI)

   integer(kind=2)                               :: itemp
   integer                                       :: i
   character(len=1), dimension(2)                :: jtemp
   character(len=1)                              :: ktemp
   equivalence  (itemp,jtemp(1))

   do i=1,n
      itemp    = a(i)
      ktemp    = jtemp(1)
      jtemp(1) = jtemp(2)
      jtemp(2) = ktemp
      a(i)     = itemp
   end do
#endif

   return
end subroutine dcw_swap16
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine swaps the order of bytes in integers or real numbers of kind=4.      !
!------------------------------------------------------------------------------------------!
subroutine dcw_swap32 (a,n)
   implicit none

   integer                       , intent(in)    :: n
   integer(kind=4) , dimension(n), intent(inout) :: a


#if defined(SGI)

   integer(kind=4)                               :: itemp
   integer                                       :: i
   character(len=1), dimension(4)                :: jtemp
   character(len=1)                              :: ktemp
   equivalence (jtemp(1),itemp)

   do i=1,n
      itemp    = a(i)
      ktemp    = jtemp(4)
      jtemp(4) = jtemp(1)
      jtemp(1) = ktemp
      ktemp    = jtemp(3)
      jtemp(3) = jtemp(2)
      jtemp(2) = ktemp
      a(i)     = itemp
   end do
#endif

   return
end subroutine dcw_swap32
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine swaps the order of bytes in real numbers of kind=8.                  !
!------------------------------------------------------------------------------------------!
subroutine dcw_swap64 (a,n)
  implicit none

   integer                       , intent(in)    :: n
   integer(kind=8) , dimension(n), intent(inout) :: a

#if defined(SGI)

   integer(kind=8)                               :: itemp
   integer                                       :: i
   character(len=1), dimension(8)                :: jtemp
   character(len=1)                              :: ktemp
   equivalence (jtemp(1),itemp)

   do i = 1,n
      itemp    = a(i)
      ktemp    = jtemp(8)
      jtemp(8) = jtemp(1)
      jtemp(1) = ktemp
      ktemp    = jtemp(7)
      jtemp(7) = jtemp(2)
      jtemp(2) = ktemp
      ktemp    = jtemp(6)
      jtemp(6) = jtemp(3)
      jtemp(3) = ktemp
      ktemp    = jtemp(5)
      jtemp(5) = jtemp(4)
      jtemp(4) = ktemp
      a(i)     = itemp
   end do
#endif
   return
end subroutine dcw_swap64
!==========================================================================================!
!==========================================================================================!
