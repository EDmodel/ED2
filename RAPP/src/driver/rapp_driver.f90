!==========================================================================================!
!==========================================================================================!
!  Subroutine rapp_driver                                                            !
!                                                                                          !
!     This is the main MEVI driver, this is actually a wrapper for the basic steps of data !
! conversion, namely the data reading, the necessary calculation, such as horizontal and/or!
! vertical interpolation and diagnostic value calculation and writing the output. The      !
! decision on how to do each step is made upon checking the namelist input.                !
!------------------------------------------------------------------------------------------!
subroutine rapp_driver()
   use mod_ioopts , only : intype  & ! intent(in)
                         , inpath  & ! intent(in)
                         , outpref & ! intent(in)
                         , iyeara  & ! intent(in)
                         , iyearz  ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: year
   !---------------------------------------------------------------------------------------!

   select case (trim(intype))
   case ('ncep')
      !------------------------------------------------------------------------------------!
      !     Here we will load the data for each year and discard everything after generat- !
      ! ing the data.  This way we save memory and use most of the structure from MEVI.    !
      !------------------------------------------------------------------------------------! 
      do year=iyeara,iyearz
         !----- 1. Loading the information table. -----------------------------------------!
         call ncep_fill_infotable(year,inpath)

         !----- 2. Filling the grid coordinate structure. ---------------------------------!
         call ncep_coordinates()
         
         !----- 3. Load the data into the matrices. ---------------------------------------!
         ! call load_variables(inpref)
         
         !----- 4. Interpolate the data from lon/lat grid to Gaussian. --------------------!
         ! call interp_gauss_driver()

         !----- 5. Write the output for ED. -----------------------------------------------!
         ! call ed_output(outpref)
      end do

   case default
      call fatal_error('Invalid intype '//trim(intype)//'!!! Can''t move on.'              &
                      ,'rapp_driver','rapp_driver.f90')
   end select

   return
end subroutine rapp_driver
!==========================================================================================!
!==========================================================================================!
