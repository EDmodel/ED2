!==========================================================================================!
!==========================================================================================!
!     Module with utilities to calculate the horizontal shading factor due to taller       !
! patches, and split patches based on the probability of being shaded by a taller          !
! neighbouring patch.                                                                      !
!------------------------------------------------------------------------------------------!
module hrzshade_utils
   implicit none

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates some derived quantities needed to run the horizontal   !
   ! shading routines.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine init_cci_variables()
      use canopy_radiation_coms, only : cci_gapmin     & ! intent(in)
                                      , cci_pixres     & ! intent(in)
                                      , cci_gapmin     & ! intent(in)
                                      , cci_gapsize    & ! intent(in)
                                      , cci_nretn      & ! intent(in)
                                      , cci_nxy        & ! intent(out)
                                      , cci_ngap       & ! intent(out)
                                      , map_ztch       & ! intent(out)
                                      , map_cci        & ! intent(out)
                                      , gap_x0         & ! intent(out)
                                      , gap_y0         & ! intent(out)
                                      , gap_ipa        ! ! intent(out)
      use disturb_cosm         , only : min_patch_area ! ! intent(in)
      !----- Local variables. -------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return
   end subroutine init_cci_variables
   !=======================================================================================!
   !=======================================================================================!


   contains


   !=======================================================================================!
   !=======================================================================================!
   !     This calculates the crown closure index and the correction factor to be applied   !
   ! to incoming direct radiation based on shading from neighbouring patches.  To find     !
   ! this, we develop a pseudo-landscape and determine the probability of a gap belonging  !
   ! to each patch to be shaded by neighbouring patches.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine split_hrzshade(csite)
      use ed_state_vars         , only : sitetype            & ! structure
                                       , patchtype           ! ! structure
                                       , allocate_sitetype   & ! subroutine
                                       , deallocate_sitetype & ! subroutine
                                       , copy_sitetype_mask  & ! subroutine
                                       , copy_sitetype       ! ! subroutine
      use canopy_radiation_coms , only : ihrzrad             & ! intent(in)
                                       , cci_radius          & ! intent(in)
                                       , cci_pixres          & ! intent(in)
                                       , cci_gapsize         & ! intent(in)
                                       , cci_gapmin          & ! intent(in)
                                       , cci_density         ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)       , pointer     :: tsite        ! Scratch site
      integer                            :: ipa          ! Counters
      real                               :: nxy          ! Size of the pseudo-raster
      real                               :: npixel       ! Number of pixels
      real                               :: absarea      ! Absolute area of pseudo-raster
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the area needed to represent all patches with at least cci_gapmin gaps.  !
      !------------------------------------------------------------------------------------!
      nxy = ceiling(sqrt(cci_gapmin / minval(csite%area)))
      !------------------------------------------------------------------------------------!



      return
   end subroutine split_hrzshade
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !    Function cci.lieberman.                                                            !
   !    This function computes the crown closure index for each individual, following the  !
   ! method presented by                                                                   !
   !                                                                                       !
   !     Lieberman, M., D. Lieberman, R. Peralta, G. S. Hartshorn, 1995: Canopy closure    !
   !        and the distribution of tropical forest tree species at La Selva, Costa Rica.  !
   !        J. Trop. Ecol., 11 (2), 161--177.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine cci_lieberman(nxyz,xyz,cci)
      use canopy_radiation_coms, only : cci_radius ! ! intent(out)
      implicit none
      !----- Variable declaration. --------------------------------------------------------!
      integer                        , intent(in)  :: nxyz
      real(kind=8), dimension(nxyz,3), intent(in)  :: xyz
      real(kind=8)                   , intent(in)  :: radius
      real(kind=8), dimension(nxyz)  , intent(out) :: cci
      !----- Local variables. -------------------------------------------------------------!
      integer                                      :: m
      integer                                      :: n
      real(kind=8)                                 :: dx
      real(kind=8)                                 :: dy
      real(kind=8)                                 :: dz
      real(kind=8)                                 :: dr
      !------------------------------------------------------------------------------------!


      !----- Initialise cci. --------------------------------------------------------------!
      cci(:) = 0.d0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through each element, find CCI.                                           !
      !------------------------------------------------------------------------------------!
      oloop: do m=1,nxyz
         !---------------------------------------------------------------------------------!
         !     Check every grid element.                                                   !
         !---------------------------------------------------------------------------------!
         iloop: do n=1,nxyz
            !----- Find distances. --------------------------------------------------------!
            dx = xyz(n,1) - xyz(m,1)
            dy = xyz(n,2) - xyz(m,2)
            dz = xyz(n,3) - xyz(m,3)
            dr = sqrt(dx*dx + dy*dy)
            !------------------------------------------------------------------------------!


            !----- Check whether point is within radius and taller than current point. ----!
            if (dr > 0.d0 .and. dr <= radius .and. dz > 0.d0) then
               cci(m) = cci(m) + sin(dz / sqrt(dz*dz + dr*dr))
            end if
            !------------------------------------------------------------------------------!
         end do iloop
         !---------------------------------------------------------------------------------!
      end do oloop
      !------------------------------------------------------------------------------------!

   end subroutine cci_lieberman
   !=======================================================================================!
   !=======================================================================================!
end module hrzshade_utils
!==========================================================================================!
!==========================================================================================!
