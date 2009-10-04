!==========================================================================================!
!==========================================================================================!
!     Subroutine interp_driver                                                             !
!                                                                                          !
!     This subroutine will control the interpolation from one grid to the other, using an  !
! objetive analysis.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine interp_driver(month)
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
                         , ssyp           & ! intent(in)
                         , sstp           ! ! intent(in)
   use mod_ioopts , only : intype         & ! intent(in)
                         , inpfrq         ! ! intent(in)
   use mod_model  , only : ngrids         ! ! intent(in)
   use mod_ncep   , only : ncep_g         ! ! intent(inout)
   use rconstants , only : day_sec        ! ! intent(in)
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer       :: month
   !----- Local variables. ----------------------------------------------------------------!
   integer       :: maxtimes
   integer       :: x
   integer       :: y
   integer       :: ng
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
      
      !------------------------------------------------------------------------------------!
      !    Finally we loop over the realisation grids to downscale precipitation.          !
      !------------------------------------------------------------------------------------!
      do ng = 4,ngrids
         do x=1,mxgauss
            do y=1,mygauss
               call a1e3(ssxp(1),ssyp(1),sstp(1),x,y,ncep_g(1)%prate,ncep_g(1)%prate1d)

               call downscale_precip(sstp(1),sstp(ng),month,ng-3                           &
                                    ,ncep_g(1)%prate1d,ncep_g(ng)%prate1d)

               call a3e1(ssxp(ng),ssyp(ng),sstp(ng),x,y,ncep_g(ng)%prate1d,ncep_g(ng)%prate)
            end do
         end do
      end do

   end select

   return
end subroutine interp_driver
!==========================================================================================!
!==========================================================================================!
