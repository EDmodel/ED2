subroutine ncep_coordinates()

   use mod_grid       , only : grid_g       & ! structure
                             , alloc_grid   & ! subroutine
                             , nullify_grid & ! subroutine
                             , zero_grid    ! ! subroutine
   use mod_model      , only : ngrids       & ! intent(in)
                             , nnxp         & ! intent(in)
                             , nnyp         & ! intent(in)
                             , nnzp         & ! intent(in)
                             , xtn          & ! intent(in)
                             , ytn          ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ifm
   integer :: x
   integer :: y
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    First thing here is to allocate the grid structure. Then we will loop through the  !
   ! grid and allocate the structure elements. We will do it only for grids that exist in  !
   ! every file (so to convert WRF files from nested grids you will probably need to run   !
   ! MEVI for each grid separately).                                                       !
   !---------------------------------------------------------------------------------------!
   allocate (grid_g(ngrids))

   gridloop: do ifm=1,ngrids
      !------------------------------------------------------------------------------------!
      !    Here we have three steps for safe allocation. First nullify, then allocate, and !
      ! finally we initialise the structure.                                               !
      !------------------------------------------------------------------------------------!
      call nullify_grid(grid_g(ifm))
      call alloc_grid(grid_g(ifm),nnxp(ifm),nnyp(ifm))
      call zero_grid(grid_g(ifm))

      !------------------------------------------------------------------------------------!
      !   Now we fill with longitude and latitude, and these should never change (RAPP     !
      ! doesn't support moving grids).                                                     !
      !------------------------------------------------------------------------------------!
      xloop: do x=1,nnxp(ifm)
         yloop: do y=1,nnyp(ifm)
            grid_g(ifm)%lon(x,y) = xtn(x,ifm)
            grid_g(ifm)%lat(x,y) = ytn(y,ifm)
         end do yloop
      end do xloop
      !------------------------------------------------------------------------------------!

   end do gridloop

   return
end subroutine ncep_coordinates
!==========================================================================================!
!==========================================================================================!
