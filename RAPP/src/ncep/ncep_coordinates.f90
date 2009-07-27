subroutine ncep_coordinates()

   use mod_grid       , only : grid_g       & ! structure
                             , alloc_grid   & ! subroutine
                             , nullify_grid & ! subroutine
                             , zero_grid    & ! subroutine
                             , ssxp         & ! intent(out)
                             , ssyp         & ! intent(out)
                             , sszp         & ! intent(out)
                             , x_1st        & ! intent(out)
                             , xlast        & ! intent(out)
                             , y_1st        & ! intent(out)
                             , ylast        & ! intent(out)
                             , z_1st        & ! intent(out)
                             , zlast        ! ! intent(out)
   use mod_model      , only : ngrids       & ! intent(in)
                             , nnxp         & ! intent(in)
                             , nnyp         & ! intent(in)
                             , nnzp         & ! intent(in)
                             , xtn          & ! intent(in)
                             , ytn          ! ! intent(in)
   use mod_ioopts     , only : lonw         & ! intent(in)
                             , lone         & ! intent(in)
                             , lats         & ! intent(in)
                             , latn         & ! intent(in)
                             , ref_hgt      ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ifm,nv
   integer :: x,xfull
   integer :: y,yfull
   !---------------------------------------------------------------------------------------!
   
   write (unit=*,fmt='(92a)') ('-',nv=1,92)
   write (unit=*,fmt='(a)')   '[+] Allocating grid structure: '

   
   !----- Assigning some default numbers for all edges except time. -----------------------!
   z_1st(1:ngrids) = 1
   zlast(1:ngrids) = nnzp(1:ngrids)
   x_1st(1:ngrids) = 1
   xlast(1:ngrids) = nnxp(1:ngrids)
   y_1st(1:ngrids) = 1
   ylast(1:ngrids) = nnyp(1:ngrids)
   
   !----- Finding the edges for the Gaussian grids, according to the user's preference. ---!
   x_1st(1) = minloc(abs(xtn(1:nnxp(1),1)-lonw),dim=1)
   xlast(1) = minloc(abs(xtn(1:nnxp(1),1)-lone),dim=1)
   y_1st(1) = minloc(abs(ytn(1:nnyp(1),1)-latn),dim=1)
   ylast(1) = minloc(abs(ytn(1:nnyp(1),1)-lats),dim=1)
   
   x_1st(2) = x_1st(1)
   xlast(2) = xlast(1)
   y_1st(2) = y_1st(1)
   ylast(2) = ylast(1)
   xtn(:,2) = xtn(:,1)
   ytn(:,2) = ytn(:,1)
   !----- Finding the sub-domain size for all grids. --------------------------------------!
   do ifm=1,ngrids
      ssxp(ifm) = xlast(ifm) - x_1st(ifm) + 1
      ssyp(ifm) = ylast(ifm) - y_1st(ifm) + 1
      sszp(ifm) = zlast(ifm) - z_1st(ifm) + 1
   end do


   !---------------------------------------------------------------------------------------!
   !    First thing here is to allocate the grid structure. Then we will loop through the  !
   ! grid and allocate the structure elements. We will do it only for grids that exist in  !
   ! every file (so to convert WRF files from nested grids you will probably need to run   !
   ! MEVI for each grid separately).                                                       !
   !---------------------------------------------------------------------------------------!
   allocate (grid_g(ngrids))

   gridloop: do ifm=1,ngrids

      write (unit=*,fmt='(a,1x,i5,a)')   '     - Grid:',ifm,'...'
      write (unit=*,fmt='(a,1x,i5,a)')   '       [-] # of points in X:',ssxp(ifm),'...'
      write (unit=*,fmt='(a,1x,i5,a)')   '       [-] # of points in Y:',ssyp(ifm),'...'
      write (unit=*,fmt='(a,1x,f9.3,a)') '       [-] Westernmost longitude:'               &
                                                                 ,xtn(x_1st(ifm),ifm),'...'
      write (unit=*,fmt='(a,1x,f9.3,a)') '       [-] Easternmost longitude:'               &
                                                                 ,xtn(xlast(ifm),ifm),'...'
      write (unit=*,fmt='(a,1x,f9.3,a)') '       [-] Northernmost latitude:'               &
                                                                 ,ytn(y_1st(ifm),ifm),'...'
      write (unit=*,fmt='(a,1x,f9.3,a)') '       [-] Southernmost latitude:'               &
                                                                 ,ytn(ylast(ifm),ifm),'...'
      write (unit=*,fmt='(a,1x,f9.3,a)') '       [-] Reference height:',ref_hgt,'...'

      !------------------------------------------------------------------------------------!
      !    Here we have three steps for safe allocation. First nullify, then allocate, and !
      ! finally we initialise the structure.                                               !
      !------------------------------------------------------------------------------------!
      call nullify_grid(grid_g(ifm))
      call alloc_grid(grid_g(ifm),ssxp(ifm),ssyp(ifm))
      call zero_grid(grid_g(ifm))

      !------------------------------------------------------------------------------------!
      !   Now we fill with longitude and latitude, and these should never change (RAPP     !
      ! doesn't support moving grids).                                                     !
      !------------------------------------------------------------------------------------!
      yloop: do y=1,ssyp(ifm)
         yfull = y_1st(ifm) + y -1
         xloop: do x=1,ssxp(ifm)
            xfull = x_1st(ifm) + x -1
            grid_g(ifm)%lon(x,y) = xtn(xfull,ifm)
            grid_g(ifm)%lat(x,y) = ytn(yfull,ifm)
            grid_g(ifm)%lev(x,y) = ref_hgt
         end do xloop
      end do yloop
      !------------------------------------------------------------------------------------!

   end do gridloop
   return
end subroutine ncep_coordinates
!==========================================================================================!
!==========================================================================================!
