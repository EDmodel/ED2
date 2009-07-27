!==========================================================================================!
!==========================================================================================!
!  Subroutine rapp_driver                                                                  !
!                                                                                          !
!     This is the main RAPP driver, this is actually a wrapper for the basic steps of data !
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

         !---------------------------------------------------------------------------------!
         ! 1. Deallocate everything before next year is called.                            !
         !---------------------------------------------------------------------------------!
         call dealloc_driver(.true.)

         !---------------------------------------------------------------------------------!
         ! 2. Loading the information table.                                               !
         !---------------------------------------------------------------------------------!
         call ncep_fill_infotable(year)

         !---------------------------------------------------------------------------------!
         ! 3. Allocate the data structures.                                                !
         !---------------------------------------------------------------------------------!
         call ncep_coordinates()

         !----- In order to save memory, we load month by month. --------------------------!
         monthloop: do month=1,12
 
            write(unit=*,fmt='(2(a,1x,i5,1x))') ' [+] Processing data.  Month=',month      &
                                                                       ,'Year=',year

            !------------------------------------------------------------------------------!
            ! 4. Filling the grid coordinate structure.                                    !
            !------------------------------------------------------------------------------!
            call ncep_alloc(month,year)

            !------------------------------------------------------------------------------!
            ! 5. Load the data into the matrices.  In case the user does not have an year  !
            !    after the given iyearz, we will need to reduce the output in one year.    !
            !------------------------------------------------------------------------------!
            success=ncep_loadvars()
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

            !------------------------------------------------------------------------------!
            ! 6. Interpolate the state variable data from lon/lat grid to Gaussian, and    !
            !    the fluxes to the more refined time scale.                                !
            !------------------------------------------------------------------------------!
            call interp_driver()

            !------------------------------------------------------------------------------!
            ! 7. Write the output for ED.                                                  !
            !------------------------------------------------------------------------------!
            call ncep_output(month,year)
         
            !------------------------------------------------------------------------------!
            ! 8. Deallocate ncep structure before next month is called.                    !
            !------------------------------------------------------------------------------!
            call dealloc_driver(.false.)
         end do monthloop
      end do yearloop 

      !----- Before we finish, we create the header. --------------------------------------!
      call ed_metd_header()

   case default
      call fatal_error('Invalid intype '//trim(intype)//'!!! Can''t move on.'              &
                      ,'rapp_driver','rapp_driver.f90')
   end select

   return
end subroutine rapp_driver
!==========================================================================================!
!==========================================================================================!
