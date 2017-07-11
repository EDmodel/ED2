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
               cci(m) = cci(m) + sin(dz / sqrt(dz*dz + dr*dr))
            else if ((.not. closure) .and. dz < 0.d0) then
               cci(m) = cci(m) - sin(dz / sqrt(dz*dz + dr*dr))
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do iloop
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

end subroutine cci_lieberman
!------------------------------------------------------------------------------------------!
