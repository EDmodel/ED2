!==========================================================================================!
!==========================================================================================!
!    Function cci.lieberman.                                                               !
!    This function computes the crown closure index for each individual, following the     !
! method presented by                                                                      !
!                                                                                          !
!     Lieberman, M., D. Lieberman, R. Peralta, G. S. Hartshorn, 1995: Canopy closure and   !
! the distribution of tropical forest tree species at La Selva, Costa Rica.  J. Trop.      !
! Ecol., 11 (2), 161--177.                                                                 !
!------------------------------------------------------------------------------------------!
subroutine cci_lieberman(nxyz,xyz,radius,closure,cci)
   implicit none
   !----- Variable declaration. -----------------------------------------------------------!
   integer                        , intent(in)  :: nxyz
   real(kind=8), dimension(nxyz,3), intent(in)  :: xyz
   real(kind=8)                   , intent(in)  :: radius
   logical                        , intent(in)  :: closure
   real(kind=8), dimension(nxyz)  , intent(out) :: cci
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: m
   integer                              :: n
   real(kind=8)                         :: dx
   real(kind=8)                         :: dy
   real(kind=8)                         :: dz
   real(kind=8)                         :: dr
   !---------------------------------------------------------------------------------------!


   !----- Initialise cci. -----------------------------------------------------------------!
   cci(:) = 0.d0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Loop through each element, find CCI.                                              !
   !---------------------------------------------------------------------------------------!
   oloop: do m=1,nxyz
      !------------------------------------------------------------------------------------!
      !     Check every grid element.                                                      !
      !------------------------------------------------------------------------------------!
      iloop: do n=1,nxyz
         !----- Find distances. -----------------------------------------------------------!
         dx = xyz(n,1) - xyz(m,1)
         dy = xyz(n,2) - xyz(m,2)
         dz = xyz(n,3) - xyz(m,3)
         dr = sqrt(dx*dx + dy*dy)
         !---------------------------------------------------------------------------------!


         !----- Check whether point is within radius. -------------------------------------!
         if (dr > 0.d0 .and. dr <= radius) then
            !----- Check whether it contributes to CCI or CII. ----------------------------!
            if (closure .and. dz > 0.d0) then
               cci(m) = cci(m) + dz / sqrt(dz*dz + dr*dr)
            else if ((.not. closure) .and. dz < 0.d0) then
               cci(m) = cci(m) - dz / sqrt(dz*dz + dr*dr)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do iloop
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine cci_lieberman
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Function cci_lieberman_mat.                                                           !
!    This function is similar to the function above, except that the data are given in     !
! matrix form.  This allows the routine to run significantly faster. It also checks for    !
! undefined numbers (and ignore them).                                                     !
!                                                                                          !
!     Lieberman, M., D. Lieberman, R. Peralta, G. S. Hartshorn, 1995: Canopy closure and   !
! the distribution of tropical forest tree species at La Selva, Costa Rica.  J. Trop.      !
! Ecol., 11 (2), 161--177.                                                                 !
!------------------------------------------------------------------------------------------!
subroutine cci_lieberman_mat(nx,ny,z,dxy,radius,undef,closure,cci)
   implicit none
   !----- Variable declaration. -----------------------------------------------------------!
   integer                       , intent(in)  :: nx
   integer                       , intent(in)  :: ny
   real(kind=8), dimension(nx,ny), intent(in)  :: z
   real(kind=8)                  , intent(in)  :: dxy
   real(kind=8)                  , intent(in)  :: radius
   real(kind=8)                  , intent(in)  :: undef
   logical                       , intent(in)  :: closure
   real(kind=8), dimension(nx,ny), intent(out) :: cci
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: x0
   integer                              :: y0
   integer                              :: xt
   integer                              :: yt
   integer                              :: x
   integer                              :: y
   integer                              :: off
   real(kind=8)                         :: dx
   real(kind=8)                         :: dy
   real(kind=8)                         :: dz
   real(kind=8)                         :: dr
   !---------------------------------------------------------------------------------------!


   !----- Initialise cci. -----------------------------------------------------------------!
   where (z == undef)
      cci = undef
   elsewhere
      cci = 0.d0
   end where
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Find out the offset for the x and y window.                                      !
   !---------------------------------------------------------------------------------------!
   off = ceiling(radius/dxy)
   !---------------------------------------------------------------------------------------!

   
   !---------------------------------------------------------------------------------------!
   !     Loop through each element, find CCI.                                              !
   !---------------------------------------------------------------------------------------!
   y0loop: do y0=1,ny
      x0loop: do x0=1,nx

         !----- Skip point if z(x0,y0) is undefined. --------------------------------------!
         if (z(x0,y0) == undef) cycle x0loop
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Loop through the window size, assume cyclic boundary conditions.            !
         !---------------------------------------------------------------------------------!
         ytloop: do yt=y0-off,y0+off
            !----- Find the actual x index to use. ----------------------------------------!
            y = 1 + modulo(yt-1,ny)
            !------------------------------------------------------------------------------!
            xtloop: do xt=x0-off,x0+off
               !----- Find the actual x index to use. -------------------------------------!
               x = 1 + modulo(xt-1,nx)
               !---------------------------------------------------------------------------!



               !----- Skip point if z(x0,y0) is undefined. --------------------------------!
               if (z(x,y) == undef) cycle xtloop
               !---------------------------------------------------------------------------!



               !----- Find distances. -----------------------------------------------------!
               dx = abs(xt - x0) * dxy
               dy = abs(yt - y0) * dxy
               dz = z(x,y) - z(x0,y0)
               dr = sqrt(dx*dx + dy*dy)
               !---------------------------------------------------------------------------!
   

               !----- Check whether point is within radius. -------------------------------!
               if (dr > 0.d0 .and. dr <= radius) then
                  !----- Check whether it contributes to CCI or CII. ----------------------!
                  if (closure .and. dz > 0.d0) then
                     cci(x0,y0) = cci(x0,y0) + dz / sqrt(dz*dz + dr*dr)
                  else if ((.not. closure) .and. dz < 0.d0) then
                     cci(x0,y0) = cci(x0,y0) - dz / sqrt(dz*dz + dr*dr)
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do xtloop
            !------------------------------------------------------------------------------!
         end do ytloop
         !---------------------------------------------------------------------------------!

      end do x0loop
      !------------------------------------------------------------------------------------!
   end do y0loop
   !---------------------------------------------------------------------------------------!

   return
end subroutine cci_lieberman_mat
!------------------------------------------------------------------------------------------!
