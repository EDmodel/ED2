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
   use allometry             , only : dbh2bl          ! ! function
   use canopy_radiation_coms , only : blfac_min       ! ! intent(in)
   use pft_coms              , only : lai_min         ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid        ! Current grid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly        ! Current polygon
   type(sitetype)   , pointer  :: csite        ! Current site
   type(patchtype)  , pointer  :: cpatch       ! Current patch
   integer                     :: ipy          ! Polygon index
   integer                     :: isi          ! Site index
   integer                     :: ipa          ! Patch index
   integer                     :: ico          ! Cohort index
   integer                     :: ipft         ! Cohort plant functional type
   integer                     :: k            ! Vertical index
   logical                     :: exposed      ! Cohort is above the snow/water [      T|F]
   logical                     :: green        ! Cohort has some leaves.        [      T|F]
   logical                     :: notwinter    ! Cohort may have leaves.        [      T|F]
   logical                     :: nottoosparse ! Cohort may have leaves.        [      T|F]
   real                        :: bleaf_pot    ! Maximum possible leaf biomass. [kgC/plant]
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy=1, cgrid%npolygons

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi=1,cpoly%nsites

         csite => cpoly%site(isi)
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- Here we integrate the total snow/surface water depth, if any. ----------!
            csite%total_snow_depth(ipa) = 0.
            sfcwloop: do k=1,csite%nlev_sfcwater(ipa)
               csite%total_snow_depth(ipa) =  csite%total_snow_depth(ipa)                  &
                                           +  csite%sfcwater_depth(k,ipa)
            end do sfcwloop

            !------------------------------------------------------------------------------!
            !   Here we actually run the check.  We will skip the cohort if any of these   !
            ! conditions happens:                                                          !
            ! 1. The cohort leaf biomass is much less than what it would have when all     !
            !    leaves are fully flushed;                                                 !
            ! 2. The cohort is buried in snow or under water                               !
            ! 3. The phenology-based green leaf factor is 0 (plants are not allowed to     !
            !    have any leaves;                                                          !
            ! 4. The cohort is not extremely sparse.                                       !
            !------------------------------------------------------------------------------!
            cohortloop: do ico=1, cpatch%ncohorts
               ipft      = cpatch%pft(ico)
               !----- Determine the maximum leaf biomass for this plant size and PFT. -----!
               bleaf_pot = dbh2bl(cpatch%dbh(ico),ipft)

               !---------------------------------------------------------------------------!
               !     Check the three conditions mentioned above.                           !
               !---------------------------------------------------------------------------!
               !----- 1. Check for relative leaf biomass. ---------------------------------!
               green        = cpatch%bleaf(ico) >= blfac_min * bleaf_pot
               !----- 2. Check for cohort height relative to snow/water depth. ------------!
               exposed      = cpatch%hite(ico)  > csite%total_snow_depth(ipa)
               !----- 3. Check for phenology. ---------------------------------------------!
               notwinter    = cpoly%green_leaf_factor(cpatch%pft(ico),isi) > 0.0
               !----- 4. Check whether this cohort is not extremely sparse. ---------------!
               nottoosparse = cpatch%lai(ico) > lai_min(ipft)
               !---------------------------------------------------------------------------!

               !----- We solve this cohort only when the four logical variables are true. -!
               cpatch%solvable(ico) = exposed .and. green .and. notwinter .and. nottoosparse
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
