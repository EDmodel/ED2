real function mortality_rates_ar(cpatch,ico,avg_daily_temp)

  use ed_state_vars,only:patchtype
  use pft_coms, only: mort1, mort2, mort3, plant_min_temp, frost_mort
  use disturb_coms, only: treefall_disturbance_rate, treefall_hite_threshold
  use consts_coms, only: t00

  implicit none
  type(patchtype),target :: cpatch
  integer :: ico,ipft
  real :: mintemp,threshtemp,cold_mort
  real, intent(in) :: avg_daily_temp
  real(kind=8) :: hugelog


  hugelog = log(huge(1.d0))

  mortality_rates_ar = 0.0

  ipft = cpatch%pft(ico)

  !--------Carbon Balance (Density-dependent)
  !    The double precision and the min statement were necessary here because of the 
  ! exponential, which can cause overflow.
  mortality_rates_ar = mortality_rates_ar + sngl(dble(mort1(ipft )) /  &
       (1.d0 + dexp(min(dble(mort2(ipft )) * dble(cpatch%cbr_bar(ico)),hugelog))))
  
  !--------Treefall
  if(cpatch%hite(ico) <= treefall_hite_threshold)then
     mortality_rates_ar = mortality_rates_ar + treefall_disturbance_rate
  endif

  !--------Cold
  mintemp = avg_daily_temp - t00
  threshtemp = 5.0 + plant_min_temp(cpatch%pft(ico))
  if(mintemp < threshtemp)then
     cold_mort = frost_mort
     if(mintemp > plant_min_temp(cpatch%pft(ico)))then
        cold_mort = cold_mort * (1.0 -   &
             (mintemp - plant_min_temp(cpatch%pft(ico)))  &
             /(threshtemp - plant_min_temp(cpatch%pft(ico))))
     endif
!     mortality_rates_ar = mortality_rates_ar + cold_mort
  endif

  !-------Density independent
  mortality_rates_ar = mortality_rates_ar + mort3(cpatch%pft(ico))

  return
end function mortality_rates_ar
