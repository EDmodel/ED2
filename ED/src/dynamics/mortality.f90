!==========================================================================================!
!==========================================================================================!
module mortality
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the total PFT-dependent mortality rate:                   !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_rates(cpatch,ico,avg_daily_temp, patch_age)
      use ed_state_vars , only : patchtype                  ! ! Structure
      use pft_coms      , only : mort0                      & ! intent(in)
                               , mort1                      & ! intent(in)
                               , mort2                      & ! intent(in)
                               , mort3                      & ! intent(in)
                               , plant_min_temp             & ! intent(in)
                               , frost_mort                 ! ! intent(in)
      use disturb_coms  , only : treefall_disturbance_rate  & ! intent(in)
                               , treefall_hite_threshold    & ! intent(in)
                               , time2canopy                ! ! intent(in)
      use ed_misc_coms  , only : current_time               ! ! intent(in)
      use ed_max_dims   , only : n_pft                      ! ! intent(in)
      use consts_coms   , only : lnexp_min                  & ! intent(in)
                               , lnexp_max                  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target     :: cpatch          ! Current patch
      integer        , intent(in) :: ico             ! Current cohort ID
      real           , intent(in) :: avg_daily_temp  ! Mean temperature yesterday
      real           , intent(in) :: patch_age       ! Patch age
      !----- Local variables --------------------------------------------------------------!
      integer                     :: ipft            ! PFT 
      real                        :: temp_dep        ! Temp. function  (frost mortality)
      real                        :: expmort         ! Carbon-balance term
      !------------------------------------------------------------------------------------!


      !----- Assume happy end, all plants survive... --------------------------------------!
      cpatch%mort_rate(1:4,ico) = 0.0
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      ! 1.  Ageing, PFT-dependent but otherwise constant.                                  !
      !------------------------------------------------------------------------------------!
      cpatch%mort_rate(1,ico) = mort3(ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 2.  Mortality rates due to negative carbon balance.                                !
      !------------------------------------------------------------------------------------!
      expmort = max( lnexp_min, min( lnexp_max                                             &
                                   , mort2(ipft) * ( cpatch%cbr_bar(ico) - mort0(ipft) ) ) )
      cpatch%mort_rate(2,ico) = mort1(ipft) / (1. + exp(expmort))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 3.  Mortality due to treefall.                                                     !
      !------------------------------------------------------------------------------------!
      if (cpatch%hite(ico) <= treefall_hite_threshold .and. patch_age > time2canopy) then
         cpatch%mort_rate(3,ico) = treefall_disturbance_rate
      else
         cpatch%mort_rate(3,ico) = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 4.   Mortality due to cold, after:                                                 !
      !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions !
      !           of land-use change, CO2 fertilization, and climate variability to the    !
      !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              !
      !           doi: 10.1111/j.1365-2486.2006.01254.x                                    !
      !------------------------------------------------------------------------------------!
      temp_dep = max( 0.0, min( 1.0, 1.0 - (avg_daily_temp - plant_min_temp(ipft)) / 5.0) )
      cpatch%mort_rate(4,ico) = frost_mort(ipft) * temp_dep
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 5. Disturbance rate mortality.  This is not used by the cohort dynamics, instead   !
      !    this is just to account for the lost density due to the patch creation.  This   !
      !    mortality will be determined by the disturbance_mortality subroutine, not here. !
      !------------------------------------------------------------------------------------!
      !cpatch%mort_rate(5,ico) = TBD
      !------------------------------------------------------------------------------------!

      return
   end subroutine mortality_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine determines the mortality rates associated with the current        !
   ! disturbance.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine disturbance_mortality(csite,ipa,area_loss,mindbh_harvest)
      use ed_state_vars, only : sitetype      & ! structure
                              , patchtype     ! ! structure
      use ed_max_dims  , only : n_pft         & ! intent(in)
                              , n_dist_types  ! ! intent(in)
      use disturb_coms , only : min_oldgrowth ! ! intent(in)
      use consts_coms  , only : lnexp_max     & ! intent(in)
                              , tiny_num      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                         , target      :: csite
      integer                                , intent(in)  :: ipa
      real          , dimension(n_dist_types), intent(in)  :: area_loss
      real          , dimension(n_pft)       , intent(in)  :: mindbh_harvest
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                        , pointer     :: cpatch
      integer                                              :: ico
      integer                                              :: new_lu
      integer                                              :: dist_path
      real                                                 :: f_survival
      real           , dimension(:)          , allocatable :: a_factor
      !------------------------------------------------------------------------------------!


      !----- Current patch, in case it is empty, return. ----------------------------------!
      cpatch => csite%patch(ipa)
      if (cpatch%ncohorts == 0) return
      !------------------------------------------------------------------------------------!


      !----- Allocate the "a_factor", which will integrate all disturbances. --------------!
      allocate(a_factor(cpatch%ncohorts))
      a_factor(:) = 0.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over new disturbance types, add survivors from each disturbance type.     !
      !------------------------------------------------------------------------------------!
      do new_lu=1,n_dist_types
         if (area_loss(new_lu) > tiny_num) then
            do ico=1,cpatch%ncohorts
              f_survival    = survivorship(new_lu,csite%dist_type(ipa),mindbh_harvest      &
                                          ,cpatch,ico)
              a_factor(ico) = a_factor(ico)                                                &
                            + ( 1.0 - f_survival ) * area_loss(new_lu) / csite%area(ipa)
            end do
         end if
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over cohorts, and find mortality.                                         !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         if ( a_factor(ico) < (1.0 - epsilon(1.0)) ) then
            cpatch%mort_rate(5,ico) = lnexp_max
         else
            cpatch%mort_rate(5,ico) = log( 1.0 / (1.0 - a_factor(ico)) )
         end if
      end do
      !------------------------------------------------------------------------------------!

      deallocate(a_factor)
      return
   end subroutine disturbance_mortality
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the survivorship rate associated with a disturbance.       !
   !  Input variables:                                                                     !
   !  -- new_lu/old_lu: the disturbance/land use type after/before disturbance             !
   !     1. Pasture.                                                                       !
   !     2. Forest plantation.                                                             !
   !     3. Tree fall.                                                                     !
   !     4. Fire.                                                                          !
   !     5. Forest regrowth.                                                               !
   !     6. Logging (tree felling).                                                        !
   !     7. Logging (collateral damage).                                                   !
   !     8. Cropland.                                                                      !
   !  -- mindbh_harvest: minimum DBH for selective logging.  All trees above threshold     !
   !                     will be logged in the tree felling patch.                         !
   !  -- cpatch: current patch.                                                            !
   !  -- ico: index for current cohort.                                                    !
   !---------------------------------------------------------------------------------------!
   real function survivorship(new_lu,old_lu,mindbh_harvest,cpatch,ico)
      use ed_state_vars, only : patchtype                ! ! structure
      use disturb_coms , only : treefall_hite_threshold  & ! intent(in)
                              , fire_hite_threshold      & ! intent(in)
                              , min_oldgrowth            ! ! intent(in)
      use pft_coms     , only : treefall_s_ltht          & ! intent(in)
                              , treefall_s_gtht          & ! intent(in)
                              , fire_s_ltht              & ! intent(in)
                              , fire_s_gtht              & ! intent(in)
                              , felling_s_gtharv         & ! intent(in)
                              , felling_s_ltharv         & ! intent(in)
                              , skid_s_ltharv            & ! intent(in)
                              , skid_s_gtharv            ! ! intent(in)
      use ed_max_dims  , only : n_pft                    ! ! intent(in)
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)                 , target     :: cpatch
      real          , dimension(n_pft), intent(in) :: mindbh_harvest
      integer                         , intent(in) :: ico
      integer                         , intent(in) :: new_lu
      integer                         , intent(in) :: old_lu
      !----- Local variables. -------------------------------------------------------------!
      integer                                      :: ipft
      !------------------------------------------------------------------------------------!


      !----- Alias for current PFT. -------------------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Survivorship depends on the disturbance type (new_lu), and within each        !
      ! disturbance type, it may depend on the disturbance pathway (dist_path), the PFT,   !
      ! and size.                                                                          !
      !------------------------------------------------------------------------------------!
      select case(new_lu)
      case (1,2,8)
         !---------------------------------------------------------------------------------!
         !     Clear cut (cropland/pasture/forest plantation).  Nothing survives.          !
         !---------------------------------------------------------------------------------!
         survivorship = 0.0
         !---------------------------------------------------------------------------------!
      case (3)
         !---------------------------------------------------------------------------------!
         !     Tree fall.  Mortality depends on the cohort height and PFT.                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < treefall_hite_threshold) then
            survivorship = treefall_s_ltht(ipft)
         else
            survivorship = treefall_s_gtht(ipft)
         end if
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     Fire.  For now fire mortality depends on the cohort height and PFT.         !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < fire_hite_threshold) then
            survivorship = fire_s_ltht(ipft)
         else
            survivorship = fire_s_gtht(ipft)
         end if
         !---------------------------------------------------------------------------------!

      case (5)
         !---------------------------------------------------------------------------------!
         !     Abandonment (secondary regrowth).  Two paths are possible: abandonment      !
         ! occurs after one last harvest, or the field/plantation is left as is.           !
         !---------------------------------------------------------------------------------!
         select case (old_lu)
         case (8)
            !----- Cropland: final harvest. -----------------------------------------------!
            survivorship = 1.0
            !------------------------------------------------------------------------------!
         case (2)
            !------------------------------------------------------------------------------!
            !     Forest plantation.  Assume typical tree fall mortality.                  !
            !------------------------------------------------------------------------------!
            if (cpatch%hite(ico) < treefall_hite_threshold) then
               survivorship = treefall_s_ltht(ipft)
            else
               survivorship = treefall_s_gtht(ipft)
            end if
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !    The only other option is pasture.  Leaving it as a default: everything    !
            ! survives.                                                                    !
            !------------------------------------------------------------------------------!
            survivorship = 1.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      case (6)
         !---------------------------------------------------------------------------------!
         !     Tree felling.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
            survivorship = felling_s_gtharv(ipft)
         else
            survivorship = felling_s_ltharv(ipft)
         end if
         !---------------------------------------------------------------------------------!
      case (7)
         !---------------------------------------------------------------------------------!
         !     Collateral damage from logging (skid trails, roads).  The damage currently  !
         ! takes into account the minimum harvest size, because presumably loggers avoid   !
         ! damaging trees that may be potentially harvested.                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
            survivorship = skid_s_gtharv(ipft)
         else
            survivorship = skid_s_ltharv(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function survivorship
   !=======================================================================================!
   !=======================================================================================!
end module mortality
!==========================================================================================!
!==========================================================================================!
