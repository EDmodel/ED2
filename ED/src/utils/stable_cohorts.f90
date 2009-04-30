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
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid  ! Current grid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly     ! Current polygon 
   type(sitetype)   , pointer  :: csite     ! Current site    
   type(patchtype)  , pointer  :: cpatch    ! Current patch   
   integer                     :: ipy       ! Polygon index   
   integer                     :: isi       ! Site index      
   integer                     :: ipa       ! Patch index     
   integer                     :: ico       ! Cohort index    
   integer                     :: k         ! Vertical index  
   real                        :: bleaf_pot ! Maximum possible leaf biomass.   [ kgC/plant]
   !----- External functions. -------------------------------------------------------------!
   ! logical          , external :: is_solvable
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
               bleaf_pot = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico))
               cpatch%solvable(ico) =                                                      &
                     cpatch%hite(ico)  > csite%total_snow_depth(ipa)                .and.  &
                     cpatch%bleaf(ico) >= blfac_min * bleaf_pot
            end do

         end do

      end do

   end do

   return
end subroutine flag_stable_cohorts
!==========================================================================================!
!==========================================================================================!
