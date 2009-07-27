!==========================================================================================!
!==========================================================================================!
!   Function walltime.  This function will define the elapsed time since the time stored   !
! at the argument, using the system clock for that.                                        !
!------------------------------------------------------------------------------------------!
subroutine walltime(wbefore,wnow)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real        , intent(in)  :: wbefore
   real        , intent(out) :: wnow
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: wbefore8
   real(kind=8)             :: wnow8
   integer                  :: cnt
   integer                  :: crate
   !----- External functions. -------------------------------------------------------------!
   real       , external    :: sngloff
   !---------------------------------------------------------------------------------------!
   wbefore8 = dble(wbefore)

   call system_clock(count=cnt,count_rate=crate)
   wnow8 = dble(cnt)/dble(crate) - wbefore8

   wnow = sngloff(wnow8,1.d-20)
   return
end subroutine walltime
