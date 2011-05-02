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
   subroutine mortality_rates(cpatch,ipa,ico,avg_daily_temp, patch_age)
      use ed_state_vars , only : patchtype                  ! ! Structure
      use pft_coms      , only : mort1                      & ! intent(in)
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
      integer        , intent(in) :: ipa             ! Current patch  ID
      integer        , intent(in) :: ico             ! Current cohort ID
      real           , intent(in) :: avg_daily_temp  ! Mean temperature yesterday
      real           , intent(in) :: patch_age
      !----- Local variables --------------------------------------------------------------!
      integer                     :: ipft            ! PFT 
      real                        :: threshtemp      ! Cold threshold temperature
      real                        :: expmort         ! Carbon-balance term
      !------------------------------------------------------------------------------------!


      !----- Assume happy end, all plants survive... --------------------------------------!
      cpatch%mort_rate(:,ico) = 0.0
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      ! 1.  Ageing, PFT-dependent but otherwise constant.                                  !
      !------------------------------------------------------------------------------------!
      cpatch%mort_rate(1,ico) = mort3(ipft)


      !------------------------------------------------------------------------------------!
      ! 2.  Mortality rates due to negative carbon balance.                                !
      !------------------------------------------------------------------------------------!
      expmort                 = max( lnexp_min, min( lnexp_max                             &
                                                   , mort2(ipft) * cpatch%cbr_bar(ico)))
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
      ! 4.   Mortality due to cold.                                                         !
      !------------------------------------------------------------------------------------!
      threshtemp = 5.0 + plant_min_temp(ipft)
      if(avg_daily_temp < threshtemp)then
         cpatch%mort_rate(4,ico) = frost_mort(ipft)
         if (avg_daily_temp > plant_min_temp(ipft)) then
            cpatch%mort_rate(4,ico) = cpatch%mort_rate(4,ico)                              &
                                    * (1.0 - (avg_daily_temp - plant_min_temp(ipft))       &
                                             / (threshtemp - plant_min_temp(ipft)) )
         end if
      else
         cpatch%mort_rate(4,ico) = 0.
      end if

      return
   end subroutine mortality_rates
   !=======================================================================================!
   !=======================================================================================!
end module mortality
!==========================================================================================!
!==========================================================================================!
