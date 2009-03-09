!==========================================================================================!
!==========================================================================================!
!    This function computes the total PFT-dependent mortality rate due to 4 factors:       !
! 1. Low carbon;                                                                           !
! 2. Treefall;                                                                             !
! 3. Cold mortality (currently off...);                                                    !
! 4. Aging(?).                                                                             !
!------------------------------------------------------------------------------------------!
real function mortality_rates_ar(cpatch,ipa,ico,avg_daily_temp)
   use ed_state_vars , only : patchtype                  ! ! Structure
   use pft_coms      , only : mort1                      & ! intent(in)
                            , mort2                      & ! intent(in)
                            , mort3                      & ! intent(in)
                            , plant_min_temp             & ! intent(in)
                            , frost_mort                 ! ! intent(in)
   use disturb_coms  , only : treefall_disturbance_rate  & ! intent(in)
                            , treefall_hite_threshold    ! ! intent(in)
   use misc_coms     , only : current_time               ! ! intent(in)
   use max_dims      , only : n_pft                      ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(patchtype)          , target     :: cpatch           ! Current patch
   integer                  , intent(in) :: ipa              ! Current patch  ID
   integer                  , intent(in) :: ico              ! Current cohort ID
   real                     , intent(in) :: avg_daily_temp   ! Mean temperature yesterday
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: ipft              ! PFT 
   real                                 :: threshtemp        ! Cold threshold temperature
   real                                 :: carbon_mort       ! Mortality due to carbon
   real                                 :: treefall_mort     ! Mortality due to treefall
   real                                 :: cold_mort         ! Mortality due to cold
   logical, dimension(n_pft), save      :: first_time=.true.
   !---------------------------------------------------------------------------------------!


   !----- Assume happy end, all plants survive... -----------------------------------------!
   mortality_rates_ar = 0.0
   ipft = cpatch%pft(ico)

   !---------------------------------------------------------------------------------------!
   !    Check mortality due to carbon balance (Density-dependent).                         !
   !---------------------------------------------------------------------------------------!
   carbon_mort = mort1(ipft) / (1. + exp(mort2(ipft) * cpatch%cbr_bar(ico)))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Mortality due to treefall.                                                        !
   !---------------------------------------------------------------------------------------!
   if(cpatch%hite(ico) <= treefall_hite_threshold) then
      treefall_mort = treefall_disturbance_rate
   else
      treefall_mort = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Mortality due to cold.                                                            !
   !---------------------------------------------------------------------------------------!
   threshtemp = 5.0 + plant_min_temp(ipft)
   if(avg_daily_temp < threshtemp)then
      cold_mort = frost_mort
      if (avg_daily_temp > plant_min_temp(ipft)) then
         cold_mort = cold_mort * (1.0 - (avg_daily_temp - plant_min_temp(ipft))            &
                                      / (threshtemp - plant_min_temp(ipft)) )
      end if
   else
      cold_mort = 0.
   end if


   !if (first_time(ipft)) then
   !   first_time(ipft) = .false.
   !   write (unit=20+ipft,fmt='(a10,8(1x,a12))')                                           &
   !      &'      TIME','       PATCH','      COHORT','      CARBON','    TREEFALL'         &
   !      &            ,'        COLD','     CBR_BAR','      HEIGHT','   YTDY_TEMP'
   !end if
   
   !write (unit=20+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),6(1x,es12.5))')                     &
   !     current_time%month,'/',current_time%date,'/',current_time%year,ipa,ico             &
   !    ,carbon_mort,treefall_mort,cold_mort                                                &
   !    ,cpatch%cbr_bar(ico),cpatch%hite(ico),avg_daily_temp

   !---------------------------------------------------------------------------------------!
   !    Find the total, density independent, mortality rate.  I don't know why, but cold   !
   ! mortality was commented out. I left it commented, but it would be good to check that. !
   !---------------------------------------------------------------------------------------!
   mortality_rates_ar = mort3(ipft) + carbon_mort + treefall_mort + cold_mort

   return
end function mortality_rates_ar
!==========================================================================================!
!==========================================================================================!
