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
                         , iyeara  & ! intent(in)
                         , iyearz  ! ! intent(inout)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: month
   integer           :: year
   integer           :: nnn
   logical           :: success
   !----- External functions. -------------------------------------------------------------!
   logical, external :: ncep_loadvars
   !---------------------------------------------------------------------------------------!

   select case (trim(intype))
   case ('ncep')
      !------------------------------------------------------------------------------------!
      !     Here we will load the data for each year and discard everything after generat- !
      ! ing the data.  This way we save memory and use most of the structure from MEVI.    !
      !------------------------------------------------------------------------------------! 
      yearloop: do year=iyeara,iyearz
         !----- 1. Loading the information table. -----------------------------------------!
         call ncep_fill_infotable(year)

         !----- 2. Allocate the data structures. ------------------------------------------!
         call ncep_coordinates()

         !----- In order to save memory, we load month by month. --------------------------!
         monthloop: do month=1,12
 
            write(unit=*,fmt='(2(a,1x,i5,1x))') ' [+] Processing data.  Month=',month      &
                                                                       ,'Year=',year
            !----- 3. Filling the grid coordinate structure. ------------------------------!
            call ncep_alloc(month,year)
         
            !------------------------------------------------------------------------------! 
            !      4. Load the data into the matrices.  In case the user does not have an  !
            !         year after the given iyearz, we will need to reduce the output in    !
            !         one year.                                                            !
            !------------------------------------------------------------------------------! 
            success=ncep_loadvars(month,year)
            if (.not. success) then
               write (unit=*,fmt='(92a)') ('-',nnn=1,92)
               write (unit=*,fmt='(a,1x,i5,a)') 'I couldn''t process YEAR=',year,'...'
               write (unit=*,fmt='(a)') 'I need at least one point in the following year!'
               write (unit=*,fmt='(a,1x,i5,1x,a)')                                         &
                                     'Quitting without processing from ',year,'onwards...'
               write (unit=*,fmt='(92a)') ('-',nnn=1,92)
               
               !----- 4a. Cleaning the arrays... ------------------------------------------!
               call dealloc_driver(.false.)
               call dealloc_driver(.true.)

               exit yearloop
            end if
            !----- 5. Interpolate the data from lon/lat grid to Gaussian. -----------------!
            ! call interp_gauss_driver()

            !----- 6. Write the output for ED. --------------------------------------------!
            ! call ed_output()
         
            !----- 7. Deallocate ncep structure before next month is called. --------------!
            call dealloc_driver(.false.)
         end do monthloop

         !----- 8. Deallocate everything before next year is called. ----------------------!
         call dealloc_driver(.true.)
      end do yearloop 

   case default
      call fatal_error('Invalid intype '//trim(intype)//'!!! Can''t move on.'              &
                      ,'rapp_driver','rapp_driver.f90')
   end select

   return
end subroutine rapp_driver
!==========================================================================================!
!==========================================================================================!
