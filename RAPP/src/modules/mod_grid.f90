!==========================================================================================!
!==========================================================================================!
! Module grid: This module contains the input grid structure.                              !
!------------------------------------------------------------------------------------------!
module mod_grid

   use mod_maxdims, only : maxgrds

   implicit none

   type grid_struct
      real, pointer, dimension(:,:  ) :: lon
      real, pointer, dimension(:,:  ) :: lat
      real, pointer, dimension(:,:  ) :: lev
   end type grid_struct

   !---------------------------------------------------------------------------------------!
   !    Full domain grid.                                                                  !
   !---------------------------------------------------------------------------------------!
   type(grid_struct), dimension(:), allocatable :: grid_g
   
   
   !----- Output dataset information. -----------------------------------------------------!
   integer, dimension(maxgrds)        :: ssxp    ! Number of points in the X direction;
   integer, dimension(maxgrds)        :: ssyp    ! Number of points in the Y direction;
   integer, dimension(maxgrds)        :: sszp    ! Number of points in the Z direction;
   integer, dimension(maxgrds)        :: sstp    ! Number of points in time;
   
   !----- Indices defining the grid edges. ------------------------------------------------!
   integer, dimension(maxgrds)        :: x_1st   ! First point to use in X direction
   integer, dimension(maxgrds)        :: xlast   ! Last  point to use in X direction
   integer, dimension(maxgrds)        :: y_1st   ! First point to use in Y direction
   integer, dimension(maxgrds)        :: ylast   ! Last  point to use in Y direction
   integer, dimension(maxgrds)        :: z_1st   ! First point to use in Z direction
   integer, dimension(maxgrds)        :: zlast   ! Last  point to use in Z direction
   integer, dimension(maxgrds)        :: t_1st   ! First time  to use
   integer, dimension(maxgrds)        :: tlast   ! Last  time  to use

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_grid(grid,nx,ny)
      implicit none 
      !----- Arguments. -------------------------------------------------------------------!
      type(grid_struct), intent(inout) :: grid
      integer          , intent(in)    :: nx,ny
      !------------------------------------------------------------------------------------!

      allocate(grid%lon (nx,ny   ))
      allocate(grid%lat (nx,ny   ))
      allocate(grid%lev (nx,ny   ))

      return
   end subroutine alloc_grid
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_grid(grid)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(grid_struct), intent(inout) :: grid
      !------------------------------------------------------------------------------------!

      if (associated(grid%lon )) nullify(grid%lon )
      if (associated(grid%lat )) nullify(grid%lat )
      if (associated(grid%lev )) nullify(grid%lev )

      return
   end subroutine nullify_grid
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine zero_grid(grid)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(grid_struct), intent(inout) :: grid
      !------------------------------------------------------------------------------------!

      if (associated(grid%lon )) grid%lon  = 0.0
      if (associated(grid%lat )) grid%lat  = 0.0
      if (associated(grid%lat )) grid%lev  = 0.0

      return
   end subroutine zero_grid
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_grid(grid)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(grid_struct), intent(inout) :: grid
      !------------------------------------------------------------------------------------!

      if (associated(grid%lon )) deallocate(grid%lon )
      if (associated(grid%lat )) deallocate(grid%lat )
      if (associated(grid%lev )) deallocate(grid%lev )

      return
   end subroutine dealloc_grid
   !=======================================================================================!
   !=======================================================================================!
end module mod_grid
!==========================================================================================!
!==========================================================================================!

