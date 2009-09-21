!==========================================================================================!
!==========================================================================================!
!   Subroutine ncep_alloc.                                                                 !
!   This subroutine will allocate all structures with the correct size for data handling.  !
!------------------------------------------------------------------------------------------!
subroutine ncep_alloc(month,year)
   use mod_ioopts , only : radfrq        & ! intent(in)
                         , inpfrq        ! ! intent(in)
   use mod_grid   , only : ssxp          & ! intent(in)
                         , ssyp          & ! intent(in)
                         , sstp          & ! intent(out)
                         , t_1st         & ! intent(out)
                         , tlast         ! ! intent(out)
   use mod_model  , only : ngrids        & ! intent(in)
                         , this_time     ! ! intent(in)
   use mod_ncep   , only : ncep_g        & ! intent(inout)
                         , flux_g        & ! intent(in)
                         , state_g       & ! intent(in)
                         , rain_g        & ! intent(in)
                         , alloc_ncep    & ! subroutine
                         , nullify_ncep  & ! subroutine
                         , init_ncep     ! ! subroutine
   use rconstants , only : day_sec       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: month         ! This month
   integer, intent(in)  :: year          ! This year
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: ng            ! Grid counter
   integer              :: ng4           ! Aux. grid counter
   integer              :: ndays         ! Number of days of this month
   integer              :: doy_1st       ! Day of year of the first day of this month
   !----- External functions. -------------------------------------------------------------!
   integer, external    :: monndays      ! Function to find ndays.
   integer, external    :: dayofyear     ! Function to find the day of year .
   !---------------------------------------------------------------------------------------!

   write (unit=*,fmt='(a)') '     - Allocating the NCEP data structure...'

   !---------------------------------------------------------------------------------------!
   !    Finding the appropriate time size.  We add an extra time that corresponds to the   !
   ! first time of the next month. This is used for the time interpolation.                !
   !---------------------------------------------------------------------------------------!
   ndays          = monndays(month,year)
   do ng=1,ngrids
      select case (ng)
      case (1,3)
         sstp(ng) = nint(day_sec/inpfrq) * ndays
      case default
         sstp(ng) = nint(day_sec/radfrq) * ndays
      end select
   end do

   !----- Finding the offset for the time considering the month we are in. ----------------!
   doy_1st         = dayofyear(month,1,year)

   do ng=1, ngrids

      select case (ng)
      case (1,3)
         t_1st(ng) = nint(day_sec/inpfrq) * (doy_1st - 1) + 1
         tlast(ng) = t_1st(ng) + sstp(ng) - 1
      case default
         t_1st(ng) = 1        ! Not really used...
         tlast(ng) = sstp(ng) ! Not really used...
      end select
   end do

   write (unit=*,fmt='(3(a,1x,i6),2(a,1x,a))')                                             &
      '         [|] Time info (Grid 1): T_1ST=',t_1st(1),'. TLAST=',tlast(1)               &
                    ,'. SSTP=',sstp(1)                                                     &
                    ,'. TIME_1ST=',trim(this_time(t_1st(1),1)%gradsstamp)                    &
                    ,'. TIMELAST=',trim(this_time(tlast(1),1)%gradsstamp)        

   !----- Now we proceed with the actual structure allocation. ----------------------------!
   allocate(ncep_g(ngrids)) 
   do ng=1,ngrids
      ng4 = min (ng,4)

      call nullify_ncep(ncep_g(ng))
      call alloc_ncep(ncep_g(ng),ssxp(ng),ssyp(ng),sstp(ng),state_g(ng4),flux_g(ng4)       &
                     ,rain_g(ng4))
      call init_ncep(ncep_g(ng),sstp(ng))
   end do

   return
end subroutine ncep_alloc
!==========================================================================================!
!==========================================================================================!






