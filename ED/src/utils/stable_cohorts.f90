!==========================================================================================!
!==========================================================================================!
!     This subroutine will go through all cohorts in a polygon, and assign them as either  !
! numerically safe or unsafe.  The idea of to use a unique test to decide which cohorts    !
! we will solve in photosynthesis, radiation, and the energy balance (RK4 or Euler).       !
!------------------------------------------------------------------------------------------!
subroutine flag_stable_cohorts(cgrid)
   use ed_state_vars         , only : edtype          & ! structure
                                    , polygontype     & ! structure
                                    , sitetype        & ! structure
                                    , patchtype       ! ! structure
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid  ! Current grid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly  ! Current polygon 
   type(sitetype)   , pointer  :: csite  ! Current site    
   type(patchtype)  , pointer  :: cpatch ! Current patch   
   integer                     :: ipy    ! Polygon index   
   integer                     :: isi    ! Site index      
   integer                     :: ipa    ! Patch index     
   integer                     :: ico    ! Cohort index    
   integer                     :: k      ! Vertical index  
   !----- External functions. -------------------------------------------------------------!
   logical          , external :: is_solvable
   !---------------------------------------------------------------------------------------!

   do ipy=1, cgrid%npolygons

      cpoly => cgrid%polygon(ipy)
      do isi=1,cpoly%nsites

         csite => cpoly%site(isi)
         do ipa=1,csite%npatches
            
            !----- Here we integrate the total snow/surface water depth, if any. ----------!
            csite%total_snow_depth(ipa) = 0
            do k=1,csite%nlev_sfcwater(ipa)
               csite%total_snow_depth(ipa) =  csite%total_snow_depth(ipa)                  &
                                           +  csite%sfcwater_depth(k,ipa)
            end do

            !----- Here we actually run the check. ----------------------------------------!
            cpatch => csite%patch(ipa)
            do ico=1, cpatch%ncohorts
               cpatch%solvable(ico) = is_solvable(cpatch%lai(ico),cpatch%bai(ico)          &
                                                 ,cpatch%hite(ico),cpatch%hcapveg(ico)     &
                                                 ,csite%total_snow_depth(ipa))
            end do

         end do

      end do

   end do

   return
end subroutine flag_stable_cohorts
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This logical function simply checks whether a cohort can be solved in a stable       !
! way.  A cohort is  solved by the RK4 integrator only when it satisfies the following     !
! three conditions:                                                                        !
! 1. The cohort total leaf+branch area index is not too low;                               !
! 2. The cohort heat capacity is not too low;                                              !
! 3. The cohort leaves aren't completely buried in snow.                                   !
!                                                                                          !
! IMPORTANT: This condition must be the same for photosynthesis, radiation, and energy     !
!            balance.                                                                      !
!------------------------------------------------------------------------------------------!
logical function is_solvable(lai,bai,hite,hcapveg,snow_depth)
   use canopy_radiation_coms , only : tai_min         ! ! intent(in)
   use rk4_coms              , only : hcapveg_coh_min ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real    , intent(in) :: lai        ! Leaf area index                       [  m²leaf/m²]
   real    , intent(in) :: bai        ! Branchwood projected area index       [  m²wood/m²]
   real    , intent(in) :: hite       ! Plant height                          [          m]
   real    , intent(in) :: hcapveg    ! Stem(trunk) projected area index      [  m²wood/m²]
   real    , intent(in) :: snow_depth ! Effective branch area index           [  m²wood/m²]
   !---------------------------------------------------------------------------------------!
   is_solvable = (lai+bai) > tai_min         .and.                                         &
                 hcapveg   > hcapveg_coh_min .and.                                         &
                 hite      > snow_depth
   return
end function is_solvable
!==========================================================================================!
!==========================================================================================!
