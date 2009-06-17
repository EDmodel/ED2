!==========================================================================================!
!==========================================================================================!
!   Subroutine dealloc_driver.                                                             !
!   This subroutine will deallocate all structures so we can proceed to next year.         !
!------------------------------------------------------------------------------------------!
subroutine dealloc_driver(new_year)
   use an_header, only : nfiles           & ! intent(in)
                       , info_table       & ! intent(inout)
                       , dealloc_anheader ! ! subroutine
   use mod_grid , only : grid_g           & ! intent(inout)
                       , dealloc_grid     ! ! subroutine
   use mod_model, only : ngrids           ! ! intent(in)
   use mod_ncep , only : ncep_g           & ! intent(inout)
                       , dealloc_ncep     ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in) :: new_year
   !----- Local variables -----------------------------------------------------------------!
   integer             :: nf
   integer             :: ng
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !   Deallocating file information table and grid structure, if necessary, it will be    !
   ! allocated again soon.                                                                 !
   !---------------------------------------------------------------------------------------!
   if (new_year) then
      if (allocated(info_table)) then
         do nf=1,nfiles
            call dealloc_anheader(info_table(nf))
         end do
         deallocate(info_table)
      end if
      
      if (allocated(grid_g)) then
         do ng=1,ngrids
            call dealloc_grid(grid_g(ng))
         end do
         deallocate(grid_g)
      end if
   !---------------------------------------------------------------------------------------!
   !   Otherwise, check whether any data structure is allocated.                           !
   !---------------------------------------------------------------------------------------!
   else if (allocated(ncep_g)) then
      do ng=1,ngrids
         call dealloc_ncep(ncep_g(ng))
      end do
      deallocate(ncep_g)
   end if

   return
end subroutine dealloc_driver
!==========================================================================================!
!==========================================================================================!
