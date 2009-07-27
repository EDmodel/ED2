!==========================================================================================!
!==========================================================================================!
!     Subroutine interp_driver                                                             !
!                                                                                          !
!     This subroutine will control the interpolation from one grid to the other, using an  !
! objetive analysis.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine interp_driver()
   use mod_interp , only : mxgauss        & ! intent(inout)
                         , mygauss        & ! intent(inout)
                         , mxlola         & ! intent(inout)
                         , mylola         & ! intent(inout)
                         , interp_buffer  & ! intent(inout)
                         , alloc_interp   & ! subroutine
                         , nullify_interp & ! subroutine
                         , init_interp    & ! subroutine
                         , dealloc_interp ! ! subroutine
   use mod_grid   , only : ssxp           & ! intent(in)
                         , ssyp           ! ! intent(in)
   use mod_ioopts , only : intype         & ! intent(in)
                         , inpfrq         ! ! intent(in)
   use rconstants , only : day_sec        ! ! intent(in)
   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   integer       :: maxtimes
   !----- Saved variables. ----------------------------------------------------------------!
   logical, save :: first_time = .true.
   !---------------------------------------------------------------------------------------!

   !------ Choosing which interpolation is needed depending on the dataset. ---------------!
   select case (trim(intype))
   case ('ncep')
      if (first_time) then
         write (unit=*,fmt='(a)') '     - Finding some objective analysis parameters...'

         !---------------------------------------------------------------------------------!
         !     This is the maximum number of times we will ever deal with, allocate the    !
         ! buffer with this size...                                                        !
         !---------------------------------------------------------------------------------!
         maxtimes = 31 * (day_sec/inpfrq) + 1 

         !---------------------------------------------------------------------------------!
         !      Assigning the alias for number of points in both the Lon/Lat and the       !
         ! Gaussian grid.                                                                  !
         !---------------------------------------------------------------------------------!
         mxgauss = ssxp(1)
         mygauss = ssyp(1)
         mxlola  = ssxp(3)
         mylola  = ssyp(3)

         !---------------------------------------------------------------------------------!
         !     If this is the first time we call, we must allocate the buffer and find the !
         ! position of the lon/lat grid with respect to the Gaussian grid.                 !
         !---------------------------------------------------------------------------------!
         call nullify_interp(interp_buffer)
         call alloc_interp(interp_buffer,maxtimes)
         call init_interp(interp_buffer,first_time)

         !----- Map the position of the lon/lat grid relative to Gaussian grid. -----------!
         call map_lolaxgauss()

         !----- Find some dimension-dependent variables. ----------------------------------!
         call assign_interp_dims()

         !----- Switching the flag to false so it will never enter here again... ----------!
         first_time = .false.
      end if

      !------------------------------------------------------------------------------------!
      !    Now we will actually interpolate the variables from the Lon/Lat to the Gaussian !
      ! grid.                                                                              !
      !------------------------------------------------------------------------------------!
      call lola_2_gauss()

      !------------------------------------------------------------------------------------!
      !    Now we will perform the time intepolation for the radiation fluxes.             !
      !------------------------------------------------------------------------------------!
      call time_interp()

   end select

end subroutine interp_driver
!==========================================================================================!
!==========================================================================================!
