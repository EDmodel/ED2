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
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid         ! Current grid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly         ! Current polygon
   type(sitetype)   , pointer  :: csite         ! Current site
   type(patchtype)  , pointer  :: cpatch        ! Current patch
   integer                     :: ipy           ! Polygon index
   integer                     :: isi           ! Site index
   integer                     :: ipa           ! Patch index
   integer                     :: ico           ! Cohort index
   integer                     :: k             ! Vertical index
   !----- External functions. -------------------------------------------------------------!
   logical          , external :: is_resolvable ! Cohort can be resolved.
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy=1, cgrid%npolygons

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi=1,cpoly%nsites

         csite => cpoly%site(isi)
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            cohortloop: do ico=1, cpatch%ncohorts

               !----- Check whether we can resolve this cohort. ---------------------------!
               cpatch%resolvable(ico) = is_resolvable(csite,ipa,ico                        &
                                                     ,cpoly%green_leaf_factor(:,isi))
               !---------------------------------------------------------------------------!

            end do cohortloop
            !------------------------------------------------------------------------------!
         end do patchloop
      end do siteloop
   end do polyloop

   return
end subroutine flag_stable_cohorts
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function finds whether a cohort can or cannot be solved.  This was taken out of !
! flag_stable_cohorts so we can call it in other places and change the cohort status as    !
! soon as the conditions change.  We will skip the cohort if any of these conditions       !
! happens:                                                                                 !
! 1. The cohort leaf biomass is much less than what it would have when all                 !
!    leaves are fully flushed;                                                             !
! 2. The cohort is buried in snow or under water                                           !
! 3. The phenology-based green leaf factor is 0 (plants are not allowed to                 !
!    have any leaves;                                                                      !
! 4. The cohort is extremely sparse.                                                       !
!------------------------------------------------------------------------------------------!
logical function is_resolvable(csite,ipa,ico,green_leaf_factor)
   use ed_state_vars , only : sitetype        & ! structure
                            , patchtype       ! ! structure
   use phenology_coms, only : elongf_min      ! ! intent(in)
   use pft_coms      , only : lai_min         ! ! intent(in)
   use ed_max_dims   , only : n_pft           ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target     :: csite             ! Current site
   integer               , intent(in) :: ipa               ! Patch index
   integer               , intent(in) :: ico               ! Cohort index
   real, dimension(n_pft), intent(in) :: green_leaf_factor ! Cold phenology scale
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype)  , pointer         :: cpatch       ! Current patch
   integer                            :: ipft         ! Cohort PFT
   logical                            :: exposed      ! Cohort is above snow    [      T|F]
   logical                            :: green        ! Cohort has some leaves. [      T|F]
   logical                            :: nottoosparse ! Cohort may have leaves. [      T|F]
   !---------------------------------------------------------------------------------------!


   cpatch => csite%patch(ipa)

   ipft      = cpatch%pft(ico)

   !---------------------------------------------------------------------------------------!
   !     Check the three conditions mentioned above.                                       !
   !---------------------------------------------------------------------------------------!
   !----- 1. Check for relative leaf biomass. ---------------------------------------------!
   green        = cpatch%elongf(ico) * green_leaf_factor(ipft) >= elongf_min
   !----- 2. Check for cohort height relative to snow/water depth. ------------------------!
   exposed      = cpatch%hite(ico)  > csite%total_sfcw_depth(ipa)
   !----- 3. Check whether this cohort is not extremely sparse. ---------------------------!
   nottoosparse = cpatch%lai(ico) > lai_min(ipft)
   !---------------------------------------------------------------------------------------!
   
   is_resolvable = green .and. exposed .and. nottoosparse
   
   return
end function is_resolvable
!==========================================================================================!
!==========================================================================================!
