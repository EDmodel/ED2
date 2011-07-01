!==========================================================================================!
!==========================================================================================!
!    This subroutine copies that variables that are integrated by the Runge-Kutta solver   !
! to a buffer structure.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine copy_patch_init(sourcesite,ipa,targetp)
   use ed_state_vars         , only : sitetype               & ! structure
                                    , patchtype              ! ! structure
   use grid_coms             , only : nzg                    & ! intent(in)
                                    , nzs                    ! ! intent(in) 
   use ed_misc_coms          , only : fast_diagnostics       ! ! intent(in)
   use consts_coms           , only : cpi8                   & ! intent(in)
                                    , ep8                    & ! intent(in)
                                    , cp8                    & ! intent(in)
                                    , epim18                 & ! intent(in)
                                    , alvl8                  & ! intent(in)
                                    , rdry8                  & ! intent(in)
                                    , rdryi8                 & ! intent(in)
                                    , p00i8                  & ! intent(in)
                                    , rocp8                  ! ! intent(in)
   use rk4_coms              , only : rk4patchtype           & ! structure
                                    , rk4site                & ! structure
                                    , rk4eps                 & ! intent(in)
                                    , any_resolvable         & ! intent(out)
                                    , zoveg                  & ! intent(out)
                                    , zveg                   & ! intent(out)
                                    , wcapcan                & ! intent(out)
                                    , wcapcani               & ! intent(out)
                                    , rk4water_stab_thresh   & ! intent(in)
                                    , rk4tiny_sfcw_mass      & ! intent(in)
                                    , checkbudget            & ! intent(in)
                                    , print_detailed         & ! intent(in)
                                    , rk4min_soil_water      & ! intent(in)
                                    , rk4max_soil_water      & ! intent(in)
                                    , find_derived_thbounds  & ! sub-routine
                                    , reset_rk4_fluxes       ! ! sub-routine
   use ed_max_dims           , only : n_pft                  ! ! intent(in)
   use therm_lib8            , only : qwtk8                  & ! subroutine
                                    , thetaeiv8              & ! function
                                    , idealdenssh8           & ! function
                                    , rehuil8                & ! function
                                    , rslif8                 & ! function
                                    , reducedpress8          ! ! function
   use soil_coms             , only : soil8                  ! ! intent(in)
   use ed_therm_lib          , only : ed_grndvap8            ! ! subroutine
   use canopy_air_coms       , only : i_blyr_condct          ! ! intent(in)
   use canopy_struct_dynamics, only : canopy_turbulence8     ! ! subroutine
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)    , target     :: targetp
   type(sitetype)        , target     :: sourcesite
   integer               , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch
   real(kind=8)                       :: rsat
   real(kind=8)                       :: sum_sfcw_mass
   integer                            :: ico
   integer                            :: ipft
   integer                            :: k
   integer                            :: ksn
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Between time steps the pressure may change because of change in atmospheric       !
   ! pressure, which means that temperature is not conserved.  Potential temperature and   !
   ! equivalent potential temperature, on the other hand, are conserved because there is   !
   ! no heat flux between time steps.  So we use these instead to start all other vari-    !
   ! ables.                                                                                !
   !---------------------------------------------------------------------------------------!
   !----- Update thermo variables that are conserved between steps. -----------------------!
   targetp%can_theta    = dble(sourcesite%can_theta(ipa))
   targetp%can_theiv    = dble(sourcesite%can_theiv(ipa))
   targetp%can_shv      = dble(sourcesite%can_shv(ipa))
   targetp%can_co2      = dble(sourcesite%can_co2(ipa))
   targetp%can_depth    = dble(sourcesite%can_depth(ipa))
   targetp%can_rvap     = targetp%can_shv / (1.d0 - targetp%can_shv)

   !----- Update the canopy pressure and Exner function. ----------------------------------!
   targetp%can_prss     = reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv &
                                       ,rk4site%geoht,targetp%can_theta,targetp%can_shv    &
                                       ,targetp%can_depth)
   targetp%can_exner    = cp8 * (targetp%can_prss * p00i8) ** rocp8
   !---------------------------------------------------------------------------------------!


   !----- Update the vegetation properties used for roughness. ----------------------------!
   targetp%veg_height   = dble(sourcesite%veg_height  (ipa))
   targetp%veg_displace = dble(sourcesite%veg_displace(ipa))
   targetp%veg_rough    = dble(sourcesite%veg_rough   (ipa))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Update the natural logarithm of theta_eiv, temperature, density, relative         !
   ! humidity, and the saturation specific humidity.                                       !
   !---------------------------------------------------------------------------------------!
   targetp%can_lntheta  = log(targetp%can_theta)
   targetp%can_temp     = cpi8 * targetp%can_theta * targetp%can_exner
   targetp%can_rhos     = idealdenssh8(targetp%can_prss,targetp%can_temp,targetp%can_shv)
   targetp%can_rhv      = rehuil8(targetp%can_prss,targetp%can_temp,targetp%can_rvap)
   rsat                 = rslif8(targetp%can_prss,targetp%can_temp)
   targetp%can_ssh      = rsat / (1.d0 + rsat)
   !---------------------------------------------------------------------------------------!

   !----- Find the lower and upper bounds for the derived properties. ---------------------!
   call find_derived_thbounds(targetp%can_rhos,targetp%can_theta,targetp%can_temp          &
                             ,targetp%can_shv,targetp%can_rvap,targetp%can_prss            &
                             ,targetp%can_depth)

   !----- Impose a non-sense number for flag_wflxgc. --------------------------------------!
   targetp%flag_wflxgc  = -1

   !---------------------------------------------------------------------------------------!
   !     Soil properties.                                                                  !
   !---------------------------------------------------------------------------------------!
   do k = rk4site%lsl, nzg
      !------------------------------------------------------------------------------------!
      !     Soil water may be slightly off-bounds.  This may happen when the soil moisture !
      ! is exactly at the bounds, but then it is copied to single precision and then back  !
      ! to double precision.  Therefore at this time only we must ensure that we bound it, !
      ! otherwise the model will crash due to the round-off error.                         !            
      !------------------------------------------------------------------------------------!
      targetp%soil_water(k)   = min( rk4max_soil_water(k)                                  &
                                   , max( rk4min_soil_water(k)                             &
                                        , dble(sourcesite%soil_water(k,ipa)) ) )
      targetp%soil_energy(k)  = dble(sourcesite%soil_energy(k,ipa))
      targetp%soil_tempk(k)   = dble(sourcesite%soil_tempk(k,ipa))
      targetp%soil_fracliq(k) = dble(sourcesite%soil_fracliq(k,ipa))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Copy the surface water information.  The only non-trivial one is the energy,      !
   ! which is saved as J/kg outside the integration, but must be converted to J/m² because !
   ! this linearises the differential equations and make the solution more stable.         !
   !---------------------------------------------------------------------------------------!
   targetp%nlev_sfcwater = sourcesite%nlev_sfcwater(ipa)
   ksn                   = targetp%nlev_sfcwater
   sum_sfcw_mass  = 0.d0
   do k = 1, nzs
      targetp%sfcwater_mass(k)    = max(0.d0,dble(sourcesite%sfcwater_mass(k,ipa)))
      targetp%sfcwater_depth(k)   = dble(sourcesite%sfcwater_depth(k,ipa))
      targetp%sfcwater_energy(k)  = dble(sourcesite%sfcwater_energy(k,ipa))                &
                                  * dble(sourcesite%sfcwater_mass(k,ipa))
      targetp%sfcwater_tempk(k)   = dble(sourcesite%sfcwater_tempk(k,ipa))
      targetp%sfcwater_fracliq(k) = dble(sourcesite%sfcwater_fracliq(k,ipa))
      sum_sfcw_mass  = sum_sfcw_mass  + targetp%sfcwater_mass (k)
   end do
   !----- Define the temporary surface water flag. ----------------------------------------!
   if (targetp%nlev_sfcwater == 0) then
      !----- No layer. --------------------------------------------------------------------!
      targetp%flag_sfcwater = 0
   elseif (sum_sfcw_mass < rk4water_stab_thresh) then
      !----- There is water, but the amount is very small. --------------------------------!
      targetp%flag_sfcwater = 1
   else
      !----- There is enough water. -------------------------------------------------------!
      targetp%flag_sfcwater = 2
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the total snow depth and fraction of canopy buried in snow.                  !
   !---------------------------------------------------------------------------------------!
   targetp%snowfac          = dble(sourcesite%snowfac(ipa))
   targetp%total_sfcw_depth = dble(sourcesite%total_sfcw_depth(ipa))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute the ground temperature and specific humidity.                             !
   !---------------------------------------------------------------------------------------!
   k = max(1,ksn)
   call ed_grndvap8(ksn,targetp%soil_water(nzg),targetp%soil_tempk(nzg)                    &
                   ,targetp%soil_fracliq(nzg),targetp%sfcwater_tempk(k)                    &
                   ,targetp%sfcwater_fracliq(k),targetp%can_prss,targetp%can_shv           &
                   ,targetp%ground_shv,targetp%ground_ssh,targetp%ground_temp              &
                   ,targetp%ground_fliq,targetp%ggsoil)
   !---------------------------------------------------------------------------------------!



   !----- Initialise some turbulence properties. ------------------------------------------!
   targetp%upwp          = dble(sourcesite%upwp  (ipa))
   targetp%wpwp          = dble(sourcesite%wpwp  (ipa))
   targetp%tpwp          = dble(sourcesite%tpwp  (ipa))
   targetp%qpwp          = dble(sourcesite%qpwp  (ipa))
   targetp%cpwp          = dble(sourcesite%cpwp  (ipa))
   !---------------------------------------------------------------------------------------!
  

   !---------------------------------------------------------------------------------------!
   !     The virtual pools should be always zero, they are temporary entities used to      !
   ! store the water that falls as throughfall or shedding in the middle of a time step.   !
   ! The temperature is liquid fraction is initialised as the soil temperature, or the     !
   ! temporary surface water temperature of the top most layer, if it exists.              !
   !---------------------------------------------------------------------------------------!
   targetp%virtual_water  = 0.0d0
   targetp%virtual_energy = 0.0d0
   targetp%virtual_depth  = 0.0d0
   if (ksn == 0) then
      targetp%virtual_tempk   = targetp%soil_tempk(nzg)
      targetp%virtual_fracliq = targetp%soil_fracliq(nzg)
   else
      targetp%virtual_tempk   = targetp%sfcwater_tempk(ksn)
      targetp%virtual_fracliq = targetp%sfcwater_tempk(ksn)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Here we find the minimum patch-level leaf heat capacity.  If the total patch leaf !
   ! heat capacity is less than this, we scale the cohorts heat capacity inside the        !
   ! integrator, so it preserves the proportional heat capacity and prevents the pool to   !
   ! be too small.                                                                         !
   !---------------------------------------------------------------------------------------!
   cpatch => sourcesite%patch(ipa)
   
   any_resolvable = .false.
   do ico=1, cpatch%ncohorts
      !----- Copying the flag that determines whether this cohort is numerically stable. --!
      targetp%leaf_resolvable(ico) = cpatch%leaf_resolvable(ico)
      targetp%wood_resolvable(ico) = cpatch%wood_resolvable(ico)
      if (targetp%leaf_resolvable(ico) .or. targetp%wood_resolvable(ico)) then
         any_resolvable = .true.
      end if
   end do


   !---------------------------------------------------------------------------------------!
   !     Loop over cohorts.  We also use this loop to compute the fraction of open canopy, !
   ! which is initialised with 1. in case this is an empty patch.                          !
   !---------------------------------------------------------------------------------------!
   targetp%opencan_frac   = dble(sourcesite%opencan_frac(ipa))
   do ico = 1,cpatch%ncohorts
      ipft=cpatch%pft(ico)

      !----- Copy the plant density. ------------------------------------------------------!
      targetp%nplant(ico)     = dble(cpatch%nplant(ico))

      !----- Copy the leaf area index and total (leaf+branch+twig) area index. ------------!
      targetp%lai(ico)        = dble(cpatch%lai(ico))
      targetp%wai(ico)        = dble(cpatch%wai(ico))
      targetp%wpa(ico)        = dble(cpatch%wpa(ico))
      targetp%tai(ico)        = targetp%lai(ico) + targetp%wai(ico)
      targetp%crown_area(ico) = dble(cpatch%crown_area(ico))
      targetp%elongf(ico)     = dble(cpatch%elongf(ico)) * rk4site%green_leaf_factor(ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether the leaves of this cohort are considered "safe" or not.  In case !
      ! it is, we copy the water, heat capacity, and energy then compute the temperature   !
      ! and the fraction of leaf water.  Otherwise, just fill with some safe values, but   !
      ! the leaves won't be really solved.                                                 !
      !------------------------------------------------------------------------------------!
      if (targetp%leaf_resolvable(ico)) then
         targetp%leaf_energy(ico) = dble(cpatch%leaf_energy(ico))
         targetp%leaf_water (ico) = max(0.d0,dble(cpatch%leaf_water (ico)))
         targetp%leaf_hcap  (ico) = dble(cpatch%leaf_hcap  (ico))

         call qwtk8(targetp%leaf_energy(ico),targetp%leaf_water(ico)                       &
                   ,targetp%leaf_hcap(ico),targetp%leaf_temp(ico),targetp%leaf_fliq(ico))
      else
         targetp%leaf_fliq  (ico) = dble(cpatch%leaf_fliq  (ico))
         targetp%leaf_temp  (ico) = dble(cpatch%leaf_temp  (ico))
         targetp%leaf_water (ico) = dble(cpatch%leaf_water (ico))
         targetp%leaf_hcap  (ico) = dble(cpatch%leaf_hcap  (ico))
         targetp%leaf_energy(ico) = targetp%leaf_hcap(ico) * targetp%leaf_temp(ico)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether the wood of this cohort is considered "safe" or not.  In case    !
      ! it is, we copy the water, heat capacity, and energy then compute the temperature   !
      ! and the fraction of leaf water.  Otherwise, just fill with some safe values, but   !
      ! the wood won't be really solved.                                                   !
      !------------------------------------------------------------------------------------!
      if (targetp%wood_resolvable(ico)) then
         targetp%wood_energy(ico) = dble(cpatch%wood_energy(ico))
         targetp%wood_water (ico) = max(0.d0,dble(cpatch%wood_water (ico)))
         targetp%wood_hcap  (ico) = dble(cpatch%wood_hcap  (ico))

         call qwtk8(targetp%wood_energy(ico),targetp%wood_water(ico)                       &
                   ,targetp%wood_hcap(ico),targetp%wood_temp(ico),targetp%wood_fliq(ico))
      else
         targetp%wood_fliq  (ico) = dble(cpatch%wood_fliq  (ico))
         targetp%wood_temp  (ico) = dble(cpatch%wood_temp  (ico))
         targetp%wood_water (ico) = dble(cpatch%wood_water (ico))
         targetp%wood_hcap  (ico) = dble(cpatch%wood_hcap  (ico))
         targetp%wood_energy(ico) = targetp%wood_hcap(ico) * targetp%wood_temp(ico)
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute the leaf intercellular specific humidity, assumed to be at saturation. !
      !------------------------------------------------------------------------------------!
      targetp%lint_shv(ico) = rslif8(targetp%can_prss,targetp%leaf_temp(ico))
      targetp%lint_shv(ico) = targetp%lint_shv(ico) / (1.d0 + targetp%lint_shv(ico))
      !------------------------------------------------------------------------------------!



      !------ Copy the stomatal conductances and the fraction of open stomata. ------------!
      targetp%fs_open   (ico) = dble(cpatch%fs_open   (ico))
      targetp%gsw_open  (ico) = dble(cpatch%gsw_open  (ico))
      targetp%gsw_closed(ico) = dble(cpatch%gsw_closed(ico))
      !------------------------------------------------------------------------------------!



      !------ Copy the net absorbed radiation (short wave and long wave). -----------------!
      targetp%rshort_l   (ico) = dble(cpatch%rshort_l (ico))
      targetp%rlong_l    (ico) = dble(cpatch%rlong_l  (ico))
      targetp%rshort_w   (ico) = dble(cpatch%rshort_w (ico))
      targetp%rlong_w    (ico) = dble(cpatch%rlong_w  (ico))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Initialise psi_open and psi_closed with zeroes, they will be averaged over the !
      ! course of one Runge-Kutta time step.                                               !
      !------------------------------------------------------------------------------------!
      targetp%psi_open(ico)   = 0.d0
      targetp%psi_closed(ico) = 0.d0
      !------------------------------------------------------------------------------------!
      
   end do
   !---------------------------------------------------------------------------------------!



   !----- Initialise the characteristic properties, and the heat capacities. --------------!
   call canopy_turbulence8(sourcesite,targetp,ipa)
   !---------------------------------------------------------------------------------------!

   !----- Diagnostics variables -----------------------------------------------------------!
   if(fast_diagnostics) then
      !------------------------------------------------------------------------------------!
      !     The "budget" variables are not copied here because they are integrated outside !
      ! RK4.  Inside RK4 we only want the contribution of those variables during the span  !
      ! of one time step.                                                                  !
      !------------------------------------------------------------------------------------!
      targetp%avg_carbon_ac      = dble(sourcesite%avg_carbon_ac(ipa)     )
      targetp%avg_vapor_lc       = dble(sourcesite%avg_vapor_lc(ipa)      )
      targetp%avg_vapor_wc       = dble(sourcesite%avg_vapor_wc(ipa)      )
      targetp%avg_dew_cg         = dble(sourcesite%avg_dew_cg(ipa)        )
      targetp%avg_vapor_gc       = dble(sourcesite%avg_vapor_gc(ipa)      )
      targetp%avg_wshed_vg       = dble(sourcesite%avg_wshed_vg(ipa)      )
      targetp%avg_intercepted    = dble(sourcesite%avg_intercepted(ipa)   )
      targetp%avg_throughfall    = dble(sourcesite%avg_throughfall(ipa)   )
      targetp%avg_vapor_ac       = dble(sourcesite%avg_vapor_ac(ipa)      )
      targetp%avg_transp         = dble(sourcesite%avg_transp(ipa)        )
      targetp%avg_evap           = dble(sourcesite%avg_evap(ipa)          )
      targetp%avg_drainage       = dble(sourcesite%avg_drainage(ipa)      )
      targetp%avg_drainage_heat  = dble(sourcesite%avg_drainage_heat(ipa) )
      targetp%avg_rshort_gnd     = dble(sourcesite%avg_rshort_gnd(ipa)    )
      targetp%avg_rlong_gnd      = dble(sourcesite%avg_rlong_gnd(ipa)     )
      targetp%avg_sensible_lc    = dble(sourcesite%avg_sensible_lc(ipa)   )
      targetp%avg_sensible_wc    = dble(sourcesite%avg_sensible_wc(ipa)   )
      targetp%avg_qwshed_vg      = dble(sourcesite%avg_qwshed_vg(ipa)     )
      targetp%avg_qintercepted   = dble(sourcesite%avg_qintercepted(ipa)  )
      targetp%avg_qthroughfall   = dble(sourcesite%avg_qthroughfall(ipa)  )
      targetp%avg_sensible_gc    = dble(sourcesite%avg_sensible_gc(ipa)   )
      targetp%avg_sensible_ac    = dble(sourcesite%avg_sensible_ac(ipa)   )

      do k = rk4site%lsl, nzg
         targetp%avg_sensible_gg(k) = dble(sourcesite%avg_sensible_gg(k,ipa))
         targetp%avg_smoist_gg(k)   = dble(sourcesite%avg_smoist_gg(k,ipa)  )
         targetp%avg_transloss(k)   = dble(sourcesite%avg_transloss(k,ipa)  )
      end do
   end if
   if (checkbudget) then
      targetp%co2budget_storage     = dble(sourcesite%co2budget_initialstorage(ipa))
      targetp%ebudget_storage       = dble(sourcesite%ebudget_initialstorage(ipa))
      targetp%wbudget_storage       = dble(sourcesite%wbudget_initialstorage(ipa))
      targetp%co2budget_loss2atm    = 0.d0
      targetp%ebudget_loss2atm      = 0.d0
      targetp%ebudget_loss2drainage = 0.d0
      targetp%ebudget_loss2runoff   = 0.d0
      targetp%wbudget_loss2atm      = 0.d0
      targetp%wbudget_loss2drainage = 0.d0
      targetp%wbudget_loss2runoff   = 0.d0
   end if


   if (print_detailed) call reset_rk4_fluxes(targetp)
   return
end subroutine copy_patch_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the carbon fluxes, which could not be copied by the time we    !
! called copy_patch_init.                                                                  !
!------------------------------------------------------------------------------------------!
subroutine copy_patch_init_carbon(sourcesite,ipa,targetp)
   use ed_state_vars        , only : sitetype              & ! structure
                                   , patchtype             ! ! structure
   use consts_coms          , only : day_sec8              & ! intent(in)
                                   , umol_2_kgC8           ! ! intent(in)
   use rk4_coms             , only : rk4patchtype          ! ! structure
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)    , target     :: targetp
   type(sitetype)        , target     :: sourcesite
   integer               , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch
   integer                            :: ico
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Here we copy the cohort level variables that are part of the carbon budget.       !
   !---------------------------------------------------------------------------------------!
   cpatch => sourcesite%patch(ipa)
   do ico = 1,cpatch%ncohorts
  
      !----- Copy the variables that are already in µmol/m²/s. ----------------------------!
      targetp%gpp         (ico) = dble(cpatch%gpp                (ico))
      targetp%leaf_resp   (ico) = dble(cpatch%leaf_respiration   (ico))
      targetp%root_resp   (ico) = dble(cpatch%root_respiration   (ico))

      !------------------------------------------------------------------------------------!
      !     The following variables are in kgC/plant/day, convert them to µmol/m²/s.       !
      !------------------------------------------------------------------------------------!
      targetp%growth_resp (ico) = dble(cpatch%growth_respiration (ico))                    &
                                * targetp%nplant(ico) / (day_sec8 * umol_2_kgC8)
      targetp%storage_resp(ico) = dble(cpatch%storage_respiration(ico))                    &
                                * targetp%nplant(ico) / (day_sec8 * umol_2_kgC8)
      targetp%vleaf_resp  (ico) = dble(cpatch%vleaf_respiration  (ico))                    &
                                * targetp%nplant(ico) / (day_sec8 * umol_2_kgC8)
   end do

   !----- Heterotrophic respiration terms. ------------------------------------------------!
   targetp%cwd_rh = dble(sourcesite%cwd_rh(ipa))
   targetp%rh     = dble(sourcesite%rh    (ipa))

   return
end subroutine copy_patch_init_carbon
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function simply checks whether the relative error is large or not.               !
!------------------------------------------------------------------------------------------!
logical function large_error(err,scal)
   use rk4_coms , only : rk4eps ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real(kind=8), intent(in) :: err  ! Absolute error
   real(kind=8), intent(in) :: scal ! Characteristic scale
   !---------------------------------------------------------------------------------------!
   if(scal > 0.d0) then
      large_error = abs(err/scal)/rk4eps > 1.d0
   else
      large_error = .false.
   end if
   return
end function large_error
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine is called before the sanity check, and updates the diagnostic vari-  !
! ables, namely the temperature and liquid fraction of leaf water, soil layers and         !
! temporary snow/pond layers.                                                                      !
!------------------------------------------------------------------------------------------!
subroutine update_diagnostic_vars(initp, csite,ipa)
   use rk4_coms              , only : rk4site               & ! intent(in)
                                    , rk4tiny_sfcw_mass     & ! intent(in)
                                    , rk4min_sfcw_mass      & ! intent(in)
                                    , rk4min_virt_water     & ! intent(in)
                                    , rk4min_can_shv        & ! intent(in)
                                    , rk4min_can_theta      & ! intent(in)
                                    , rk4max_can_theta      & ! intent(in)
                                    , rk4min_can_lntheta    & ! intent(in)
                                    , rk4max_can_lntheta    & ! intent(in)
                                    , rk4min_can_temp       & ! intent(in)
                                    , rk4max_can_shv        & ! intent(in)
                                    , rk4min_veg_lwater     & ! intent(in)
                                    , rk4min_veg_temp       & ! intent(in)
                                    , rk4max_veg_temp       & ! intent(in)
                                    , rk4min_soil_temp      & ! intent(in)
                                    , rk4max_soil_temp      & ! intent(in)
                                    , rk4min_soil_water     & ! intent(in)
                                    , rk4max_soil_water     & ! intent(in)
                                    , rk4min_sfcw_temp      & ! intent(in)
                                    , rk4max_sfcw_temp      & ! intent(in)
                                    , rk4water_stab_thresh  & ! intent(in)
                                    , tiny_offset           & ! intent(in)
                                    , force_idealgas        & ! intent(in)
                                    , rk4patchtype          ! ! structure
   use ed_state_vars         , only : sitetype              & ! structure
                                    , patchtype             ! ! structure
   use soil_coms             , only : soil8                 & ! intent(in)
                                    , dslz8                 & ! intent(in)
                                    , dslzi8                ! ! intent(in)
   use grid_coms             , only : nzg                   & ! intent(in)
                                    , nzs                   ! ! intent(in)
   use therm_lib8            , only : qwtk8                 & ! subroutine
                                    , qtk8                  & ! subroutine
                                    , thetaeiv8             & ! function
                                    , rehuil8               & ! function
                                    , rslif8                & ! function
                                    , thrhsh2temp8          ! ! function
   use consts_coms           , only : alvl8                 & ! intent(in)
                                    , wdns8                 & ! intent(in)
                                    , rdryi8                & ! intent(in)
                                    , rdry8                 & ! intent(in)
                                    , epim18                & ! intent(in)
                                    , toodry8               & ! intent(in)
                                    , cp8                   & ! intent(in)
                                    , cpi8                  & ! intent(in)
                                    , p00i8                 & ! intent(in)
                                    , rocp8                 & ! intent(in)
                                    , t3ple8                & ! intent(in)
                                    , cliq8                 & ! intent(in)
                                    , cice8                 & ! intent(in)
                                    , tsupercool8           ! ! intent(in)
   use canopy_struct_dynamics, only : canopy_turbulence8    ! ! subroutine
   use ed_therm_lib          , only : ed_grndvap8           ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer :: cpatch
   integer                          :: ico
   integer                          :: k
   integer                          :: ksn
   integer                          :: kclosest
   logical                          :: ok_shv
   logical                          :: ok_theta
   logical                          :: ok_ground
   logical                          :: ok_sfcw
   logical                          :: ok_leaf
   logical                          :: ok_wood
   real(kind=8)                     :: soilhcap
   real(kind=8)                     :: int_sfcw_energy
   real(kind=8)                     :: int_virt_energy
   real(kind=8)                     :: energy_tot
   real(kind=8)                     :: wmass_tot
   real(kind=8)                     :: hcapdry_tot
   real(kind=8)                     :: rk4min_leaf_water
   real(kind=8)                     :: rk4min_wood_water
   !---------------------------------------------------------------------------------------!

   !----- Then we define some logicals to make the code cleaner. --------------------------!
   ok_shv   = initp%can_shv     >= rk4min_can_shv     .and.                                &
              initp%can_shv     <= rk4max_can_shv
   ok_theta = initp%can_lntheta >= rk4min_can_lntheta .and.                                &
              initp%can_lntheta <= rk4max_can_lntheta

   !---------------------------------------------------------------------------------------!
   !     Here we convert theta into temperature, potential temperature, and density, and   !
   ! ice-vapour equivalent potential temperature.  The latter variable (or its natural     !
   ! log) should eventually become the prognostic variable for canopy air space entropy    !
   ! when we add condensed/frozen water in the canopy air space.                           !
   !---------------------------------------------------------------------------------------!
   if (ok_shv .and. ok_theta) then

      !----- First, we update the canopy air potential temperature. -----------------------!
      initp%can_theta = exp(initp%can_lntheta)

      initp%can_rvap  = initp%can_shv / (1.d0 - initp%can_shv)

      !------------------------------------------------------------------------------------!
      !    Here we find the temperature in different ways depending on whether we are      !
      ! keeping pressure constant during one full time step or not.  If we are forcing     !
      ! ideal gas to be always respected, then we don't know the pressure until we have    !
      ! the temperature, so we compute the temperature based on potential temperature,     !
      ! density, and specific humidity, then update pressure.  Otherwise, we compute the   !
      ! temperature using the known pressure, even though this causes the ideal gas law to !
      ! not be always satisfied.                                                           !
      !------------------------------------------------------------------------------------!
      if (force_idealgas) then
         initp%can_temp  = thrhsh2temp8(initp%can_theta,initp%can_rhos,initp%can_shv)
         initp%can_prss  = initp%can_rhos * rdry8 * initp%can_temp                         &
                         * (1.d0 + epim18 * initp%can_shv)
         initp%can_exner = cp8 * (initp%can_prss * p00i8) ** rocp8
      else
         initp%can_temp  = cpi8 * initp%can_theta * initp%can_exner
      end if
      !------------------------------------------------------------------------------------!


      initp%can_rhv   = rehuil8(initp%can_prss,initp%can_temp,initp%can_rvap)
      initp%can_ssh   = rslif8(initp%can_prss,initp%can_temp)
      initp%can_ssh   = initp%can_ssh / (initp%can_ssh + 1.d0)
      initp%can_theiv = thetaeiv8(initp%can_theta,initp%can_prss,initp%can_temp            &
                                 ,initp%can_rvap,initp%can_rvap)
   elseif (initp%can_lntheta >= rk4max_can_lntheta) then
      initp%can_theta = rk4max_can_theta + 1.d0
   elseif (initp%can_lntheta <= rk4min_can_lntheta) then
      initp%can_theta = rk4min_can_theta - 1.d0
   end if
   !---------------------------------------------------------------------------------------!





   !----- Update soil temperature and liquid water fraction. ------------------------------!
   do k = rk4site%lsl, nzg
      soilhcap = soil8(rk4site%ntext_soil(k))%slcpd
      call qwtk8(initp%soil_energy(k),initp%soil_water(k)*wdns8,soilhcap                   &
                ,initp%soil_tempk(k),initp%soil_fracliq(k))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update surface water temperature and liquid water fraction, remembering that in-   !
   ! side the RK4 integration, surface water energy is in J/m². The abs is necessary be-   !
   ! cause surface mass may indeed become too negative during the integration process and  !
   ! if it happens, we want the step to be rejected.                                       !
   !---------------------------------------------------------------------------------------!
   ok_sfcw = .true.
   ksn = initp%nlev_sfcwater
   initp%total_sfcw_depth = 0.d0
   sfcwloop: do k=1,ksn
      if (initp%sfcwater_mass(k) < rk4min_sfcw_mass) then 
         !---------------------------------------------------------------------------------!
         !     Surface water mass doesn't make sense.  Skip because the step is going to   !
         ! be rejected and we may not be able to compute the temperature.                  !
         !---------------------------------------------------------------------------------!
         ok_sfcw = .false.
         cycle sfcwloop
      elseif (initp%flag_sfcwater == 1) then
         !---------------------------------------------------------------------------------!
         !     Water layer is too thin to be computationally stable, apply the quick heat  !
         ! exchange between the soil top layer and the temporary surface water.  Find the  !
         ! total internal energy of the combined pool (top soil layer plus the thin        !
         ! temporary surface water).  The units of soil properties are J/m3 for the        !
         ! internal energy, and m3/m3 for soil water, whilst the temporary surface water   !
         ! has units of J/m2 for internal energy and kg/m2 for mass.  We use the standard  !
         ! for the temporary surface water.                                                !
         !---------------------------------------------------------------------------------!
         energy_tot  = initp%sfcwater_energy(k) + initp%soil_energy(nzg) * dslz8(nzg)
         wmass_tot   = initp%sfcwater_mass  (k)                                            &
                     + initp%soil_water (nzg) * dslz8(nzg) * wdns8
         hcapdry_tot = soil8(rk4site%ntext_soil(nzg))%slcpd * dslz8(nzg)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the equilibrium temperature and liquid/ice partition.   Because we    !
         ! are assuming thermal equilibrium, the temperature and liquid fraction of the    !
         ! attempted layer is the same as the average temperature of the augmented pool.   !
         !---------------------------------------------------------------------------------!
         call qwtk8(energy_tot,wmass_tot,hcapdry_tot                                       &
                   ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
         initp%soil_tempk(nzg)   = initp%sfcwater_tempk(k)
         initp%soil_fracliq(nzg) = initp%sfcwater_fracliq(k) 
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Re-compute the internal energy of the temporary layer, using the temperature !
         ! and fraction of liquid water distribution we have just found, keeping the mass  !
         ! constant.                                                                       !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_energy(k) = initp%sfcwater_mass(k)                                 &
                                  * ( initp%sfcwater_fracliq(k)                            &
                                    * cliq8 * (initp%sfcwater_tempk(k) - tsupercool8)      &
                                    + (1.d0 - initp%sfcwater_fracliq(k))                   &
                                      * cice8 * initp%sfcwater_tempk(k)  )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Re-calculate the top soil internal energy, by removing the attempted surface !
         ! water energy from the total energy, and converting it back to J/m3.  The total  !
         ! amount of water does not need to be re-calculated at this time.                 !
         !---------------------------------------------------------------------------------!
         initp%soil_energy(nzg)  = (energy_tot - initp%sfcwater_energy(k)) * dslzi8(nzg)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Convert surface water energy from extensive quantity (J/m2) to intensive   !
         ! quantity (J/kg), then update the temperature and liquid water fraction.  We     !
         ! only do this when there is enough mass, because otherwise the amount of energy  !
         ! is too small that the code would give some unreasonable results.  Also, this is !
         ! the minimum amount of water needed for a layer to exist, anything less than     !
         ! that would be entirely absorbed by the soil layer.                              !
         !---------------------------------------------------------------------------------!
         int_sfcw_energy = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
         call qtk8(int_sfcw_energy,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      end if
      initp%total_sfcw_depth = initp%total_sfcw_depth + initp%sfcwater_depth(k)
   end do sfcwloop
   !---------------------------------------------------------------------------------------!
   !    For non-existent layers of temporary surface water, we copy the temperature and    !
   ! liquid fraction from the layer beneath.                                               !
   !---------------------------------------------------------------------------------------!
   aboveloop: do k=ksn+1,nzs
      if (k == 1) then
         initp%sfcwater_tempk  (k) = initp%soil_tempk  (nzg)
         initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
      else
         initp%sfcwater_tempk  (k) = initp%sfcwater_tempk  (k-1)
         initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
      end if
   end do aboveloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Update the fraction of the canopy covered in snow.                                !
   !---------------------------------------------------------------------------------------!
   initp%snowfac = min(9.9d-1,initp%total_sfcw_depth/initp%veg_height)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Update the virtual pool temperature and liquid water fraction, only when the      !
   ! total mass of the virtual pool makes sense and is not zero.  When the amount of water !
   ! is too small, we impose the same temperature as the exposed surface (either the top   !
   ! soil layer or the top temporary surface water layer), because the amount of energy is !
   ! too small to get an acurate temperature anyway, and all the water and energy are      !
   ! going to be absorbed by the top soil layer.                                           !
   !---------------------------------------------------------------------------------------!
   if (initp%virtual_water < rk4min_virt_water) then
      continue
   elseif (abs(initp%virtual_water) > rk4tiny_sfcw_mass) then
      int_virt_energy = initp%virtual_energy / initp%virtual_water
      call qtk8(int_virt_energy,initp%virtual_tempk,initp%virtual_fracliq)
   elseif (ksn == 0) then
      initp%virtual_tempk   = initp%soil_tempk(nzg)
      initp%virtual_fracliq = initp%soil_fracliq(nzg)
   else
      initp%virtual_tempk   = initp%sfcwater_tempk(ksn)
      initp%virtual_fracliq = initp%sfcwater_fracliq(ksn)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Compute the ground temperature and specific humidity.                             !
   !---------------------------------------------------------------------------------------!
   k = max(1,ksn)
   if (ksn == 0) then
      k = 1
      ok_ground = initp%soil_tempk(nzg) >= rk4min_soil_temp       .and.                    &
                  initp%soil_tempk(nzg) <= rk4max_soil_temp       .and.                    &
                  initp%soil_water(nzg) >= rk4min_soil_water(nzg) .and.                    &
                  initp%soil_water(nzg) <= rk4max_soil_water(nzg)
   else
      k = ksn
      ok_ground = initp%sfcwater_tempk(ksn) >= rk4min_sfcw_temp .and.                      &
                  initp%sfcwater_tempk(ksn) <= rk4max_sfcw_temp
   end if
   if (ok_ground) then
      call ed_grndvap8(ksn,initp%soil_water(nzg),initp%soil_tempk(nzg)                     &
                      ,initp%soil_fracliq(nzg),initp%sfcwater_tempk(k)                     &
                      ,initp%sfcwater_fracliq(k),initp%can_prss,initp%can_shv              &
                      ,initp%ground_shv,initp%ground_ssh,initp%ground_temp                 &
                      ,initp%ground_fliq,initp%ggsoil)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Loop over cohorts to update the leaf temperature and liquid fraction.             !
   !---------------------------------------------------------------------------------------!
   ok_leaf = .true.
   cpatch => csite%patch(ipa)
   leafloop: do ico=1,cpatch%ncohorts
      !----- Check whether the leaves of this cohort are safe... --------------------------!
      if (initp%leaf_resolvable(ico)) then

         !----- Find the minimum leaf water for this cohort. ------------------------------!
         rk4min_leaf_water = rk4min_veg_lwater * initp%lai(ico)

         !---------------------------------------------------------------------------------!
         !     Update leaf temperature and liquid fraction only if leaf water makes sense. !
         !---------------------------------------------------------------------------------!
         if (initp%leaf_water(ico) < rk4min_leaf_water) then
            ok_leaf = .false.
            cycle leafloop
         else
            call qwtk8(initp%leaf_energy(ico),initp%leaf_water(ico),initp%leaf_hcap(ico)   &
                      ,initp%leaf_temp(ico),initp%leaf_fliq(ico))
            if (initp%leaf_temp(ico) < rk4min_veg_temp .or.                                &
                initp%leaf_temp(ico) > rk4max_veg_temp) then
               ok_leaf = .false.
               cycle leafloop
            else
               !---------------------------------------------------------------------------!
               !     Compute the leaf intercellular specific humidity, assumed to be at    !
               ! saturation.                                                               !
               !---------------------------------------------------------------------------!
               initp%lint_shv(ico) = rslif8(initp%can_prss,initp%leaf_temp(ico))
               initp%lint_shv(ico) = initp%lint_shv(ico) / (1.d0 + initp%lint_shv(ico))
               !---------------------------------------------------------------------------!
            end if
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     For plants with minimal foliage or very sparse patches, fix the leaf        !
         ! temperature to the canopy air space and force leaf_water to be zero.            !
         !---------------------------------------------------------------------------------!
         initp%leaf_temp(ico)   = initp%can_temp
         initp%leaf_water(ico)  = 0.d0
         initp%leaf_energy(ico) = initp%leaf_hcap(ico) * initp%leaf_temp(ico)
         if (initp%leaf_temp(ico) == t3ple8) then
            initp%leaf_fliq(ico) = 5.d-1
         elseif (initp%leaf_temp(ico) > t3ple8) then
            initp%leaf_fliq(ico) = 1.d0
         else
            initp%leaf_fliq(ico) = 0.d0
         end if
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Compute the leaf intercellular specific humidity, assumed to be at          !
         ! saturation.                                                                     !
         !---------------------------------------------------------------------------------!
         initp%lint_shv(ico) = rslif8(initp%can_prss,initp%leaf_temp(ico))
         initp%lint_shv(ico) = initp%lint_shv(ico) / (1.d0 + initp%lint_shv(ico))
         !---------------------------------------------------------------------------------!
      end if
   end do leafloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Loop over cohorts to update the wood temperature and liquid fraction.             !
   !---------------------------------------------------------------------------------------!
   ok_wood = .true.
   cpatch => csite%patch(ipa)
   woodloop: do ico=1,cpatch%ncohorts
      !----- Check whether the leaves of this cohort are safe... --------------------------!
      if (initp%wood_resolvable(ico)) then

         !----- Find the minimum leaf water for this cohort. ------------------------------!
         rk4min_wood_water = rk4min_veg_lwater * initp%wai(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Update wood temperature and liquid fraction only if wood water makes sense. !
         !---------------------------------------------------------------------------------!
         if (initp%wood_water(ico) < rk4min_wood_water) then
            ok_wood = .false.
            cycle woodloop
         else
            call qwtk8(initp%wood_energy(ico),initp%wood_water(ico),initp%wood_hcap(ico)   &
                      ,initp%wood_temp(ico),initp%wood_fliq(ico))
            if (initp%wood_temp(ico) < rk4min_veg_temp .or.                                &
                initp%wood_temp(ico) > rk4max_veg_temp) then
               ok_leaf = .false.
               cycle woodloop
            end if
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     In case we are not solving branches, or for cohorts that are very sparse,   !
         ! fix the wood temperature to the canopy air space and force wood_water to be     !
         ! zero.                                                                           !
         !---------------------------------------------------------------------------------!
         initp%wood_temp(ico)   = initp%can_temp
         initp%wood_water(ico)  = 0.d0
         initp%wood_energy(ico) = initp%wood_hcap(ico) * initp%wood_temp(ico)
         if (initp%wood_temp(ico) == t3ple8) then
            initp%wood_fliq(ico) = 5.d-1
         elseif (initp%wood_temp(ico) > t3ple8) then
            initp%wood_fliq(ico) = 1.d0
         else
            initp%wood_fliq(ico) = 0.d0
         end if
         !---------------------------------------------------------------------------------!
      end if
   end do woodloop
   !---------------------------------------------------------------------------------------!


   !----- Compute canopy turbulence properties. -------------------------------------------!
   if (ok_leaf .and. ok_wood .and. ok_shv .and. ok_theta .and. ok_ground) then
      call canopy_turbulence8(csite,initp,ipa)
   end if
   !---------------------------------------------------------------------------------------!


   return
end subroutine update_diagnostic_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the following tasks:                                         !
! 1. Check how many layers of temporary water or snow we have, and include the virtual     !
!    pools at the topmost if needed;                                                       !
! 2. Force thermal equilibrium between topmost soil layer and a single snow/water layer    !
!    if the layer is too thin;                                                             !
! 3. Compute the amount of mass each layer has, and redistribute them accordingly.         !
! 4. Percolates excessive liquid water if needed.                                          !
!------------------------------------------------------------------------------------------!
subroutine adjust_sfcw_properties(nzg,nzs,initp,hdid,csite,ipa)

   use rk4_coms      , only : rk4patchtype          & ! structure
                            , rk4site               & ! intent(in)
                            , rk4min_sfcw_mass      & ! intent(in)
                            , rk4min_virt_water     & ! intent(in)
                            , rk4water_stab_thresh  & ! intent(in)
                            , rk4tiny_sfcw_mass     & ! intent(in)
                            , rk4tiny_sfcw_depth    & ! intent(in)
                            , rk4min_can_shv        & ! intent(in)
                            , rk4snowmin            & ! intent(in)
                            , ipercol               & ! intent(in)
                            , rk4eps2               & ! intent(in)
                            , wcapcan               & ! intent(in)
                            , wcapcani              ! ! intent(in)
   use ed_state_vars , only : sitetype              & ! structure
                            , patchtype             ! ! structure
   use soil_coms     , only : soil8                 & ! intent(in)
                            , dslz8                 & ! intent(in)
                            , dslzi8                & ! intent(in)
                            , thick                 & ! intent(in)
                            , thicknet              ! ! intent(in)
   use consts_coms   , only : cice8                 & ! intent(in)
                            , cliq8                 & ! intent(in)
                            , t3ple8                & ! intent(in)
                            , wdns8                 & ! intent(in)
                            , wdnsi8                & ! intent(in)
                            , tsupercool8           & ! intent(in)
                            , qliqt38               & ! intent(in)
                            , wdnsi8                & ! intent(in)
                            , fdnsi8                & ! intent(in)
                            , alli8                 & ! intent(in)
                            , alvi8                 ! ! intent(in)
   use therm_lib8    , only : qtk8                  & ! subroutine
                            , qwtk8                 ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp
   type(sitetype)         , target     :: csite
   integer                , intent(in) :: ipa
   real(kind=8)           , intent(in) :: hdid
   integer                , intent(in) :: nzg
   integer                , intent(in) :: nzs
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: kold
   integer                             :: newlayers
   integer                             :: nlayers
   integer                             :: ksn
   integer                             :: ksnnew
   integer                             :: k
   !----- Control variables ---------------------------------------------------------------!
   real(kind=8)                        :: hdidi
   real(kind=8)                        :: wtold
   real(kind=8)                        :: wtnew
   real(kind=8), dimension(nzs)        :: newsfcw_mass
   real(kind=8), dimension(nzs)        :: newsfcw_energy
   real(kind=8), dimension(nzs)        :: newsfcw_depth
   real(kind=8)                        :: wdiff
   real(kind=8)                        :: sum_sfcw_mass
   real(kind=8)                        :: sum_sfcw_energy
   real(kind=8)                        :: sum_sfcw_depth
   real(kind=8)                        :: energy_free
   real(kind=8)                        :: wmass_free
   real(kind=8)                        :: depth_free
   real(kind=8)                        :: energy_available
   real(kind=8)                        :: wmass_available
   real(kind=8)                        :: depth_available
   real(kind=8)                        :: energy_needed
   real(kind=8)                        :: wmass_needed
   real(kind=8)                        :: depth_needed
   real(kind=8)                        :: wmass_perc
   real(kind=8)                        :: energy_perc
   real(kind=8)                        :: depth_perc
   real(kind=8)                        :: i_energy_try
   real(kind=8)                        :: energy_try
   real(kind=8)                        :: wmass_try
   real(kind=8)                        :: depth_try
   real(kind=8)                        :: temp_try
   real(kind=8)                        :: fliq_try
   real(kind=8)                        :: energy_tot
   real(kind=8)                        :: wmass_tot
   real(kind=8)                        :: hcapdry_tot
   real(kind=8)                        :: wmass_room
   real(kind=8)                        :: energy_room
   real(kind=8)                        :: depthloss
   real(kind=8)                        :: snden
   real(kind=8)                        :: sndenmin
   real(kind=8)                        :: sndenmax
   real(kind=8)                        :: Cr               ! snow waterholding capacity
   real(kind=8)                        :: gi               ! Partial density of ice
   integer                             :: nsoil
   !----- Variables used for the water and energy budget. ---------------------------------!
   real(kind=8)                        :: wmass_cas_beg
   real(kind=8)                        :: wmass_cas_end
   real(kind=8)                        :: energy_input
   real(kind=8)                        :: wmass_virtual_beg
   real(kind=8)                        :: energy_virtual_beg
   real(kind=8)                        :: wmass_sfcw_beg
   real(kind=8)                        :: energy_sfcw_beg
   real(kind=8)                        :: wmass_soil_beg
   real(kind=8)                        :: energy_soil_beg
   real(kind=8)                        :: wmass_total_beg
   real(kind=8)                        :: energy_total_beg
   real(kind=8)                        :: wmass_virtual_end
   real(kind=8)                        :: energy_virtual_end
   real(kind=8)                        :: wmass_sfcw_end
   real(kind=8)                        :: energy_sfcw_end
   real(kind=8)                        :: wmass_soil_end
   real(kind=8)                        :: energy_soil_end
   real(kind=8)                        :: wmass_total_end
   real(kind=8)                        :: energy_total_end
   real(kind=8)                        :: wmass_total_rch
   real(kind=8)                        :: energy_total_rch
   !----- Constants -----------------------------------------------------------------------!
   logical                , parameter  :: debug   = .false.
   real(kind=8)           , parameter  :: Crmin   = 3.d-2
   real(kind=8)           , parameter  :: Crmax   = 1.d-1
   real(kind=8)           , parameter  :: ge      = 2.d2
   !---------------------------------------------------------------------------------------!


   !----- Find the inverse of the time step. ----------------------------------------------!
   hdidi      = 1.d0 / hdid
   !---------------------------------------------------------------------------------------!


   !----- Copy the original number of temporary surface water layers to ksn. --------------!
   ksn       = initp%nlev_sfcwater
   !---------------------------------------------------------------------------------------!



   !----- Copy the soil type at the topmost level to nsoil. -------------------------------!
   nsoil     = rk4site%ntext_soil(nzg)
   !---------------------------------------------------------------------------------------!
   


   !---------------------------------------------------------------------------------------!
   !      Determine the total amount of temporary surface water available as well as       !
   ! derived properties.                                                                   !
   !---------------------------------------------------------------------------------------!
   sum_sfcw_mass   = 0.d0
   sum_sfcw_energy = 0.d0
   sum_sfcw_depth  = 0.d0
   do k=1,ksn
      sum_sfcw_mass   = sum_sfcw_mass   + initp%sfcwater_mass  (k)
      sum_sfcw_energy = sum_sfcw_energy + initp%sfcwater_energy(k)
      sum_sfcw_depth  = sum_sfcw_depth  + initp%sfcwater_depth (k)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise the budget variables.                                                 !
   !---------------------------------------------------------------------------------------!
   wmass_cas_beg      = initp%can_shv * wcapcan
   energy_input       = 0.d0
   wmass_virtual_beg  = initp%virtual_water
   energy_virtual_beg = initp%virtual_energy
   wmass_sfcw_beg     = sum_sfcw_mass
   energy_sfcw_beg    = sum_sfcw_energy
   wmass_soil_beg     = initp%soil_water(nzg)  * dslz8(nzg) * wdns8
   energy_soil_beg    = initp%soil_energy(nzg) * dslz8(nzg)
   wmass_total_beg    = wmass_virtual_beg  + wmass_sfcw_beg  + wmass_soil_beg              &
                      + wmass_cas_beg
   energy_total_beg   = energy_virtual_beg + energy_sfcw_beg + energy_soil_beg
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Check the total amount of water that has just fallen plus the amount that is al-  !
   ! ready sitting over the top soil layer.  We must do this as the first step because we  !
   ! may want to eliminate this water by binding it to the top soil layer in case there is !
   ! too little water.                                                                     !
   !---------------------------------------------------------------------------------------!
   if (initp%virtual_water < rk4min_virt_water .or. sum_sfcw_mass < rk4min_sfcw_mass ) then
      !------------------------------------------------------------------------------------!
      !     Either the virtual layer or the temporary surface water has too negative mass, !
      ! so this step doesn't make sense.  We quit the sub-routine here so the sanity check !
      ! can reject this step.                                                              !
      !------------------------------------------------------------------------------------!
      return
   elseif ((initp%virtual_water + sum_sfcw_mass) < 0.d0) then
      !------------------------------------------------------------------------------------!
      !     The mass of the potential new temporary surface water is within bounds but it  !
      ! is negative.  This can only happen if the layer evaporated more than what it       !
      ! should, so we condense some of the water back from the canopy air space.  If it is !
      ! going to deplete the canopy air space specific humidity too much, then we leave    !
      ! the remainder to be stolen from the top soil, but this is dangerous because the    !
      ! soil may be too dry too.                                                           !
      !------------------------------------------------------------------------------------!
      wmass_needed               = - (initp%virtual_water  + sum_sfcw_mass  )
      energy_needed              = - (initp%virtual_energy + sum_sfcw_energy)
      depth_needed               = - (initp%virtual_depth  + sum_sfcw_depth )

      !----- Find the amount available at the canopy air space. ---------------------------!
      wmass_available            = wcapcan * (initp%can_shv - 5.d0 * rk4min_can_shv)
      if ( wmass_available > wmass_needed) then
         !---------------------------------------------------------------------------------!
         !    There is enough water vapour. The transfer will require phase change, so the !
         ! energy transfer will be a latent heat flux.  Remove the water from the canopy   !
         ! air space.                                                                      !
         !---------------------------------------------------------------------------------!
         initp%can_shv      = initp%can_shv       - wmass_needed  * wcapcani
         initp%avg_vapor_gc = initp%avg_vapor_gc  - wmass_needed * hdidi

         energy_input       = energy_needed

         wmass_free  = 0.d0
         energy_free = 0.d0
         depth_free  = 0.d0

      elseif (wmass_available > 0.d0) then
         !---------------------------------------------------------------------------------!
         !    There is not enough water vapour. Dry down to the minimum, and hope for the  !
         ! best.                                                                           !
         !---------------------------------------------------------------------------------!
         energy_available   = wmass_available * (alvi8 - initp%soil_fracliq(nzg) * alli8)
         depth_available    = wmass_available * ( initp%soil_fracliq(nzg) * wdnsi8         &
                                                + (1.d0-initp%soil_fracliq(nzg)) * fdnsi8)

         energy_input       = energy_available

         initp%can_shv      = initp%can_shv       - wmass_available * wcapcani
         initp%avg_vapor_gc = initp%avg_vapor_gc  - wmass_available * hdidi
         
         wmass_free         = wmass_available  - wmass_needed 
         energy_free        = energy_available - energy_needed
         depth_free         = depth_available  - depth_needed
      else
         !---------------------------------------------------------------------------------!
         !    There is not any water vapour.  Hope for the best.                           !
         !---------------------------------------------------------------------------------!
         wmass_free         = wmass_needed
         energy_free        = energy_needed
         depth_free         = depth_needed
      end if

      !----- Reset both the temporary surface water and the virtual layer. ----------------!
      initp%virtual_water      = 0.d0
      initp%virtual_energy     = 0.d0
      initp%virtual_depth      = 0.d0
      initp%sfcwater_mass  (:) = 0.d0
      initp%sfcwater_energy(:) = 0.d0
      initp%sfcwater_depth (:) = 0.d0
      !----- Set ksnnew to zero to force all free water to go to the soil. ----------------!
      ksnnew                   = 0
   elseif ((initp%virtual_water + sum_sfcw_mass) < rk4tiny_sfcw_mass) then
      !------------------------------------------------------------------------------------!
      !     The mass of the potential new temporary surface water is within bounds but it  !
      ! is too small to be maintained.  We add both the virtual mass and the total surface !
      ! water and dump in the free water, but set ksnnew to zero so all the water is       !
      ! infiltrated in the top soil layer.                                                 !
      !------------------------------------------------------------------------------------!
      wmass_free               = initp%virtual_water  + sum_sfcw_mass
      energy_free              = initp%virtual_energy + sum_sfcw_energy
      depth_free               = initp%virtual_depth  + sum_sfcw_depth
      !----- Reset both the temporary surface water and the virtual layer. ----------------!
      initp%virtual_water      = 0.d0
      initp%virtual_energy     = 0.d0
      initp%virtual_depth      = 0.d0
      initp%sfcwater_mass  (:) = 0.d0
      initp%sfcwater_energy(:) = 0.d0
      initp%sfcwater_depth (:) = 0.d0
      !----- Set ksnnew to zero to force all free water to go to the soil. ----------------!
      ksnnew                   = 0
   else
      !------------------------------------------------------------------------------------!
      !     The mass of the potential new temporary surface water is within bounds and     !
      ! could create at least one layer.  If there is already a temporary surface water or !
      ! snow layer, the new amount is initially put there, otherwise, we attempt to create !
      ! the first layer.                                                                   !
      !------------------------------------------------------------------------------------!
      wmass_free               = initp%virtual_water
      energy_free              = initp%virtual_energy
      depth_free               = initp%virtual_depth
      initp%virtual_water      = 0.d0
      initp%virtual_energy     = 0.d0
      initp%virtual_depth      = 0.d0
      ksnnew                   = max(ksn,1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update the prognostic and diagnostic variables by adding the free standing water.  !
   ! Then we check the size of the temporary surface water layers, and update the          !
   ! temperature in a way that ensure the layer stability.  During this process, we        !
   ! update the total temporary surface water mass, energy, and depth, which will be used  !
   ! later in the sub-routine.                                                             !
   !---------------------------------------------------------------------------------------!
   sum_sfcw_mass   = 0.d0
   sum_sfcw_energy = 0.d0
   sum_sfcw_depth  = 0.d0
   do k = ksnnew,1,-1
      !------------------------------------------------------------------------------------!
      !    Find the potential mass, energy, and depth of the temporary layer if all the    !
      ! free water became part of this layer.                                              !
      !------------------------------------------------------------------------------------!
      energy_try = initp%sfcwater_energy(k) + energy_free
      wmass_try  = initp%sfcwater_mass(k)   + wmass_free
      depth_try  = initp%sfcwater_depth(k)  + depth_free
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    In case this is a single layer, and a very thin one, we may have a hard time    !
      ! achieving numerical stability.  We can treat this case like the leaf case, in      !
      ! the sense that the water sitting on the top of the surface is in thermal           !
      ! equilibrium with the surface.                                                      !
      !------------------------------------------------------------------------------------!
      if (ksnnew == 1 .and. wmass_try < rk4water_stab_thresh) then
         !---------------------------------------------------------------------------------!
         !     Find the total internal energy of the combined pool (top soil layer plus    !
         ! the thin temporary surface water).  The units of soil properties are J/m3 for   !
         ! the internal energy, and m3/m3 for soil water, whilst the temporary surface     !
         ! water has units of J/m2 for internal energy and kg/m2 for mass.  We use the     !
         ! standard for the temporary surface water.                                       !
         !---------------------------------------------------------------------------------!
         energy_tot  = energy_try + initp%soil_energy(nzg) * dslz8(nzg)
         wmass_tot   = wmass_try  + initp%soil_water(nzg)  * dslz8(nzg) * wdns8
         hcapdry_tot = soil8(nsoil)%slcpd * dslz8(nzg)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the equilibrium temperature and liquid/ice partition.   Because we    !
         ! are assuming thermal equilibrium, the temperature and liquid fraction of the    !
         ! attempted layer is the same as the average temperature of the augmented pool.   !
         !---------------------------------------------------------------------------------!
         call qwtk8(energy_tot,wmass_tot,hcapdry_tot,temp_try,fliq_try)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Re-compute the internal energy of the temporary layer, using the temperature !
         ! and fraction of liquid water distribution we have just found, keeping the mass  !
         ! constant.                                                                       !
         !---------------------------------------------------------------------------------!
         energy_try = wmass_try * (         fliq_try  * cliq8 * (temp_try - tsupercool8)   &
                                  + (1.d0 - fliq_try) * cice8 *  temp_try                )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Re-calculate the top soil internal energy, by removing the attempted surface !
         ! water energy from the total energy, and converting it back to J/m3.  The total  !
         ! amount of water does not need to be re-calculated at this time.                 !
         !---------------------------------------------------------------------------------!
         initp%soil_energy(nzg)  = (energy_tot - energy_try) * dslzi8(nzg)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Layer is computationally stable, find temperature and liquid fraction of   !
         ! the attempted layer.                                                            !
         !---------------------------------------------------------------------------------!
         i_energy_try = energy_try / wmass_try
         call qtk8(i_energy_try,temp_try,fliq_try)
        !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Determine a first guess for the amount of mass that can be lost from this      !
      ! layer through percolation (wmass_perc).                                            !
      !------------------------------------------------------------------------------------!
      select case (ipercol)
      case (0)
         !---------------------------------------------------------------------------------!
         !     Original method, from LEAF-3.  Shed liquid in excess of a 1:9               !
         ! liquid-to-ice ratio through percolation.                                        !
         !---------------------------------------------------------------------------------!
         wmass_perc  = max(0.d0, wmass_try * (fliq_try - 1.d-1) / 9.d-1)
         !---------------------------------------------------------------------------------!
      case (1,2)
         !---------------------------------------------------------------------------------!
         !    Alternative "free" water calculation.                                        !
         !    Anderson (1976), NOAA Tech Report NWS 19.                                    !
         !---------------------------------------------------------------------------------!
         gi          = wmass_try/max(rk4tiny_sfcw_depth,depth_try) * (1.d0 - fliq_try)
         Cr          = max(Crmin, Crmin + (Crmax - Crmin) * (ge - gi) / ge)
         wmass_perc  = max(0.d0,wmass_try * (fliq_try - Cr / (1.d0 + Cr)))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Determinte whether the layer beneath the current one is another temporary      !
      ! surface water/snow layer, or the top soil layer.  In case it is the latter, we     !
      ! must check whether there is enough room for the percolate water to infiltrate      !
      ! (i.e., the soil will not become super-saturated), in which case we must reduce the !
      ! total amount of percolation.                                                       !
      !------------------------------------------------------------------------------------!
      if (k == 1) then
         !---------------------------------------------------------------------------------!
         !     Compute the available "room" for water at the top soil layer.  We must      !
         ! multiply by density and depth to make sure that the units match.                !
         !---------------------------------------------------------------------------------!
         wmass_room = max(0.d0, soil8(nsoil)%slmsts - initp%soil_water(nzg))               &
                    * wdns8 * dslz8(nzg) 
         wmass_perc = min(wmass_perc,wmass_room)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Re-calculate the total water mass and energy of this temporary surface water.  !
      ! Here we must check whether the soil layer would be with too little mass, and if    !
      ! that is the case, we will eliminate the layer by forcing the tiny left-over to go  !
      ! to the layer beneath.                                                              !
      !------------------------------------------------------------------------------------!
      if (wmass_try - wmass_perc > rk4tiny_sfcw_mass) then
         !---------------------------------------------------------------------------------!
         !      Enough mass to keep this layer.                                            !
         !---------------------------------------------------------------------------------!
         !----- Compute the internal energy and depth associated with percolated water. ---!
         energy_perc = wmass_perc * cliq8 * (temp_try - tsupercool8)
         depth_perc  = wmass_perc * wdnsi8
         !----- Find the new water mass and energy for this layer. ------------------------!
         initp%sfcwater_mass  (k) = wmass_try  - wmass_perc
         initp%sfcwater_energy(k) = energy_try - energy_perc

         !---------------------------------------------------------------------------------!
         !      Calculate density and depth of snow.  Start with the difference of depths, !
         ! but then we adjust it because the loss through percolation changes the ratio    !
         ! between ice and liquid in this layer                                            !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_depth (k) = depth_try  - depth_perc
         snden    = initp%sfcwater_mass(k)                                                 &
                  / max(rk4tiny_sfcw_depth,initp%sfcwater_depth(k))
         sndenmax = wdns8
         sndenmin = max(3.d1, 2.d2 * (wmass_free + wmass_perc) / initp%sfcwater_mass(k) )
         snden    = min(sndenmax, max(sndenmin,snden))
         initp%sfcwater_depth (k) = initp%sfcwater_mass(k) / snden
      else
         !---------------------------------------------------------------------------------!
         !      The layer would be too small, eliminate mass from this layer and send all  !
         ! mass to the layer beneath as percolated water.                                  !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_mass  (k) = 0.d0
         initp%sfcwater_energy(k) = 0.d0
         initp%sfcwater_depth (k) = 0.d0
         wmass_perc               = wmass_try
         energy_perc              = energy_try
         depth_perc               = depth_try
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate the total temporary surface water properties.                        !
      !------------------------------------------------------------------------------------!
      sum_sfcw_mass   = sum_sfcw_mass   + initp%sfcwater_mass  (k)
      sum_sfcw_energy = sum_sfcw_energy + initp%sfcwater_energy(k)
      sum_sfcw_depth  = sum_sfcw_depth  + initp%sfcwater_depth (k)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The water available for the layer beneath is going to be the total percolated  !
      ! water.                                                                             !
      !------------------------------------------------------------------------------------!
      wmass_free  = wmass_perc
      energy_free = energy_perc
      depth_free  = depth_perc
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     There may be a tiny amount of free standing water left.  We dump what we can in   !
   ! the soil, and if there is still some water to be removed we  evaporate what is left.  !
   !---------------------------------------------------------------------------------------!
   if (wmass_free > 0.d0) then
      wmass_room  = max(0.d0, soil8(nsoil)%slmsts - initp%soil_water(nzg))                 &
                  * wdns8 * dslz8(nzg)
      energy_room = energy_free * wmass_room / wmass_free

      if (wmass_room > wmass_free) then
         !---------------------------------------------------------------------------------!
         !     There is enough space in the top soil layer for the remaining water, put    !
         ! all the free water there.                                                       !
         !---------------------------------------------------------------------------------!
         initp%soil_water(nzg)  = initp%soil_water(nzg)  + wmass_free  * dslzi8(nzg)       &
                                * wdnsi8
         initp%soil_energy(nzg) = initp%soil_energy(nzg) + energy_free * dslzi8(nzg)

         wmass_free  = 0.d0
         energy_free = 0.d0
         depth_free  = 0.d0
      else
         !----- Remove the water that can go to the soil. ---------------------------------!
         wmass_free  = wmass_free  - wmass_room
         energy_free = energy_free - energy_room

         !----- Dump what we can dump on the top soil layer. ------------------------------!
         initp%soil_water(nzg)  = initp%soil_water(nzg)  + wmass_room  * dslzi8(nzg)       &
                                * wdnsi8
         initp%soil_energy(nzg) = initp%soil_energy(nzg) + energy_room * dslzi8(nzg)

         !----- Boil the remaining. -------------------------------------------------------!
         initp%can_shv      = initp%can_shv       + wmass_free  * wcapcani
         initp%avg_vapor_gc = initp%avg_vapor_gc  + wmass_free  * hdidi

         energy_input       = -energy_free

         wmass_free  = 0.d0
         energy_free = 0.d0
         depth_free  = 0.d0
      end if
   elseif (wmass_free < 0.d0) then
      wmass_needed     = - wmass_free
      energy_needed    = - energy_free
      depth_needed     = - depth_free

      !------ Find the amount of water that the soil can provide. -------------------------!
      wmass_available  = max(0.d0,initp%soil_water(nzg) - soil8(nsoil)%soilcp)             &
                       * wdns8 * dslz8(nzg)
      energy_available = energy_free * wmass_available / wmass_free

      if (wmass_available > wmass_needed) then
         !---------------------------------------------------------------------------------!
         !     There is enough space in the top soil layer to correct remaining negative   !
         ! water, get all the water needed there.                                          !
         !---------------------------------------------------------------------------------!
         initp%soil_water(nzg)  = initp%soil_water(nzg)  - wmass_needed  * dslzi8(nzg)     &
                                * wdnsi8
         initp%soil_energy(nzg) = initp%soil_energy(nzg) - energy_needed * dslzi8(nzg)
         wmass_needed     = 0.d0
         energy_needed    = 0.d0
         depth_needed     = 0.d0
      else
         !----- Add the water that can come from the soil. --------------------------------!
         wmass_needed  = wmass_needed  - wmass_available
         energy_needed = energy_needed - energy_available

         !----- Dump what we can dump on the top soil layer. ------------------------------!
         initp%soil_water(nzg)  = initp%soil_water(nzg)  - wmass_available  * dslzi8(nzg)  &
                                * wdnsi8
         initp%soil_energy(nzg) = initp%soil_energy(nzg) - energy_available * dslzi8(nzg)

         !----- Condense the remaining, hoping for the best. ------------------------------!
         initp%can_shv      = initp%can_shv       - wmass_needed  * wcapcani
         initp%avg_vapor_gc = initp%avg_vapor_gc  - wmass_needed  * hdidi

         energy_input       = - energy_needed

         wmass_free  = 0.d0
         energy_free = 0.d0
         depth_free  = 0.d0

      end if
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check the total amount of mass in the temporary surface water/snow, and adjust    !
   ! the number of layer accordingly.                                                      !
   !---------------------------------------------------------------------------------------!
   if (sum_sfcw_mass <= rk4tiny_sfcw_mass) then
      !----- Not enough water in the temporary surface water, eliminate all layers. -------!
      initp%nlev_sfcwater = 0
      !------------------------------------------------------------------------------------!


      !----- Update the flag for temporary surface water. ---------------------------------!
      initp%flag_sfcwater = 0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      The total mass should be either zero or greater than rk4tiny_sfcw_mass,       !
      ! but, just in case, we add any remaining energy to the top soil layer.              !
      !------------------------------------------------------------------------------------!
      initp%soil_water(nzg)  = initp%soil_water(nzg)  + sum_sfcw_mass   * dslzi8(nzg)      &
                                                      * wdnsi8
      initp%soil_energy(nzg) = initp%soil_energy(nzg) + sum_sfcw_energy * dslzi8(nzg)
      !------------------------------------------------------------------------------------!

      !----- Loop all layers and re-set all extensive variables to zero. ------------------!
      do k = 1, nzs
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
      end do
      !------------------------------------------------------------------------------------!
   elseif (sum_sfcw_mass < rk4water_stab_thresh) then


      !----- Not much water in the temporary surface water, impose a single layer. --------!
      initp%nlev_sfcwater = 1
      !------------------------------------------------------------------------------------!


      !----- Update the flag for temporary surface water. ---------------------------------!
      initp%flag_sfcwater = 1
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    If the total amount of temporary surface water is not enough to make it stable, !
      ! we impose it to have a single layer with all the ponding/snow in there.            !
      !------------------------------------------------------------------------------------!
      initp%sfcwater_mass  (1) = sum_sfcw_mass
      initp%sfcwater_energy(1) = sum_sfcw_energy
      initp%sfcwater_depth (1) = sum_sfcw_depth
      do k=2,nzs
         initp%sfcwater_mass   (k) = 0.d0
         initp%sfcwater_energy (k) = 0.d0
         initp%sfcwater_depth  (k) = 0.d0
      end do
      !------------------------------------------------------------------------------------!


   else
      !----- Update the flag for temporary surface water. ---------------------------------!
      initp%flag_sfcwater = 2
      !------------------------------------------------------------------------------------!


      !---- Check whether there is enough snow for a new layer. ---------------------------!
      nlayers   = ksnnew
      newlayers = 1
      do k = 1,ksnnew
         !---------------------------------------------------------------------------------!
         !     Check whether the layer as is meet the minimum requirements to stand as a   !
         ! new layer by itself.                                                            !
         !---------------------------------------------------------------------------------!
         if ( initp%sfcwater_mass(k)   >  rk4tiny_sfcw_mass              .and.             &
              rk4snowmin * thicknet(k) <= sum_sfcw_mass                  .and.             &
              initp%sfcwater_energy(k) <  initp%sfcwater_mass(k)*qliqt38       ) then
            newlayers = newlayers + 1
         end if
         !---------------------------------------------------------------------------------!
      end do

      !----- Newlayers is the new number of temporary surface water layers. ---------------!
      newlayers = min(newlayers, nzs, nlayers + 1)
      
      if (newlayers == 1) then
         newsfcw_mass  (1) = sum_sfcw_mass
         newsfcw_energy(1) = sum_sfcw_energy
         newsfcw_depth (1) = sum_sfcw_depth
      else
         kold  = 1
         wtnew = 1.d0
         wtold = 1.d0
         do k = 1,newlayers
            newsfcw_mass(k)   = sum_sfcw_mass * thick(k,newlayers)
            newsfcw_energy(k) = 0.d0
            newsfcw_depth(k)  = 0.d0
            !----- Find the properties of this new layer. ---------------------------------!
            find_layer: do

               !----- Difference between old and new snow ---------------------------------!
               wdiff = wtnew * newsfcw_mass(k) - wtold * initp%sfcwater_mass(kold)  

               if (wdiff > 0.d0) then
                  newsfcw_energy(k) = newsfcw_energy(k)                                    &
                                    + wtold * initp%sfcwater_energy(kold)
                  newsfcw_depth(k)  = newsfcw_depth(k)                                     &
                                    + wtold * initp%sfcwater_depth(kold)
                  wtnew  = wtnew - wtold * initp%sfcwater_mass(kold) / newsfcw_mass(k)
                  kold   = kold + 1
                  wtold  = 1.0
                  if (kold > nlayers) exit find_layer
               else
                  newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)             &
                                    * initp%sfcwater_energy(kold)                             &
                                    / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                  newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)             &
                                    * initp%sfcwater_depth(kold)                              &
                                    / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                  wtold = wtold - wtnew * newsfcw_mass(k)                                     &
                                / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                  wtnew = 1.
                  exit find_layer
               end if
            end do find_layer
         end do
      end if

      !----- Update the water/snow layer prognostic properties. ---------------------------!
      initp%nlev_sfcwater = newlayers
      do k = 1,newlayers
         initp%sfcwater_mass(k)   = newsfcw_mass(k)
         initp%sfcwater_energy(k) = newsfcw_energy(k)
         initp%sfcwater_depth(k)  = newsfcw_depth(k)
      end do
      do k = newlayers + 1, nzs
         initp%sfcwater_mass   (k) = 0.d0
         initp%sfcwater_energy (k) = 0.d0
         initp%sfcwater_depth  (k) = 0.d0
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Compute the budget variables after the adjustments.                              !
   !---------------------------------------------------------------------------------------!
   wmass_cas_end      = initp%can_shv * wcapcan
   wmass_virtual_end  = initp%virtual_water
   energy_virtual_end = initp%virtual_energy
   wmass_sfcw_end     = sum_sfcw_mass
   energy_sfcw_end    = sum_sfcw_energy
   wmass_soil_end     = initp%soil_water(nzg)  * dslz8(nzg) * wdns8
   energy_soil_end    = initp%soil_energy(nzg) * dslz8(nzg)
   wmass_total_end    = wmass_virtual_end  + wmass_sfcw_end  + wmass_soil_end              &
                      + wmass_cas_end
   energy_total_end   = energy_virtual_end + energy_sfcw_end + energy_soil_end
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Check whether energy and mass are conserved.                                     !
   !---------------------------------------------------------------------------------------!
   wmass_total_rch  = 2.d0 * abs(wmass_total_end - wmass_total_beg)                        &
                    / (abs(wmass_total_end) + abs(wmass_total_beg))
   energy_total_rch = 2.d0 * abs(energy_total_end - energy_input - energy_total_beg)       &
                    / (abs(energy_total_end - energy_input) + abs(energy_total_beg))
   if (wmass_total_rch > 1.d-6 .or. energy_total_rch > 1.d-6) then
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a)')           ' Water or energy conservation was violated!!!   '
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           ' - Initial conditions: '
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total water mass    = ',wmass_total_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + CAS mass            = ',wmass_cas_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Virtual mass        = ',wmass_virtual_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Ponding/snow mass   = ',wmass_sfcw_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Soil mass           = ',wmass_soil_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total energy        = ',energy_total_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Virtual energy      = ',energy_virtual_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Ponding/snow energy = ',energy_sfcw_beg
      write (unit=*,fmt='(a,1x,es14.7)') '   + Soil energy         = ',energy_soil_beg
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           ' - Final conditions: '
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total water mass    = ',wmass_total_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + CAS mass            = ',wmass_cas_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Virtual mass        = ',wmass_virtual_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Ponding/snow mass   = ',wmass_sfcw_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Soil mass           = ',wmass_soil_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total energy        = ',energy_total_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Virtual energy      = ',energy_virtual_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Ponding/snow energy = ',energy_sfcw_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Soil energy         = ',energy_soil_end
      write (unit=*,fmt='(a,1x,es14.7)') '   + Input energy (cond) = ',energy_input
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           ' - Relative error: '
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total water mass    = ',wmass_total_rch
      write (unit=*,fmt='(a,1x,es14.7)') '   + Total energy        = ',energy_total_rch
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           '------------------------------------------------'
      call fatal_error('Energy or water is not being conserved!!!'                         &
                      ,'adjust_sfcw_properties','rk4_misc.f90')
   end if

   return
end subroutine adjust_sfcw_properties
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will ensure that top soil properties are within the allowed range.   !
! We currently test this only for the top soil layer because this is the most likely to    !
! cause problems, but if it does happen in other layers, we could easily extend this for   !
! all layers.  Depending on its derivative, soil moisture can go under the minimum soil    !
! moisture possible (soilcp) or above the saturation (slmsts).  Both are bad things, and   !
! if the value is way off-bounds, then we leave it like that so the step can be rejected.  !
! However, if the value is just slightly off(*) these limits, we make a small exchange of  !
! moisture with the neighbouring layer.  This will prevent the soil to go outside the      !
! range in those double precision => single precision => double precision conversion.      !
!                                                                                          !
! (*) slightly off is defined as outside the range but within the desired accuracy         !
!     (rk4eps).                                                                            !
!------------------------------------------------------------------------------------------!
subroutine adjust_topsoil_properties(initp,hdid,csite,ipa)
   use rk4_coms             , only : rk4patchtype         & ! structure
                                   , rk4site              & ! intent(in)
                                   , rk4eps               & ! intent(in)
                                   , rk4tiny_sfcw_mass    & ! intent(in)
                                   , rk4min_sfcw_mass     & ! intent(in)
                                   , rk4min_can_shv       & ! intent(in)
                                   , rk4min_soil_water    & ! intent(in)
                                   , rk4max_soil_water    & ! intent(in)
                                   , wcapcan              & ! intent(in)
                                   , wcapcani             & ! intent(in)
                                   , hcapcani             ! ! intent(in)
   use ed_state_vars        , only : sitetype             & ! structure
                                   , patchtype            ! ! structure
   use consts_coms          , only : cice8                & ! intent(in)
                                   , cliq8                & ! intent(in)
                                   , alvl8                & ! intent(in)
                                   , alvi8                & ! intent(in)
                                   , alli8                & ! intent(in)
                                   , t3ple8               & ! intent(in)
                                   , wdns8                & ! intent(in)
                                   , fdnsi8               & ! intent(in)
                                   , toodry8              & ! intent(in)
                                   , tsupercool8          & ! intent(in)
                                   , qliqt38              & ! intent(in)
                                   , wdnsi8               ! ! intent(in)
   use therm_lib8           , only : qwtk8                & ! subroutine
                                   , qtk8                 ! ! subroutine
   use grid_coms            , only : nzg                  ! ! intent(in)
   use soil_coms            , only : soil8                & ! intent(in)
                                   , dslzi8               & ! intent(in)
                                   , dslz8                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp  ! Integration buffer
   type(sitetype)         , target     :: csite  ! Current site
   integer                , intent(in) :: ipa    ! Current patch ID
   real(kind=8)           , intent(in) :: hdid   ! Time step 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer    :: cpatch
   integer                             :: ico
   integer                             :: ksn
   integer                             :: kt
   integer                             :: kb
   integer                             :: kw
   integer                             :: nstop
   integer                             :: nsbeneath
   real(kind=8)                        :: hdidi
   real(kind=8)                        :: shctop
   real(kind=8)                        :: shcbeneath
   real(kind=8)                        :: water_needed
   real(kind=8)                        :: energy_needed
   real(kind=8)                        :: depth_needed
   real(kind=8)                        :: energy_excess
   real(kind=8)                        :: water_excess
   real(kind=8)                        :: depth_excess
   real(kind=8)                        :: water_available
   real(kind=8)                        :: energy_available
   real(kind=8)                        :: water_room
   real(kind=8)                        :: energy_room
   logical                             :: slightlymoist
   logical                             :: slightlydry
   !---------------------------------------------------------------------------------------!

   !----- Inverse of time increment -------------------------------------------------------!
   hdidi = 1.d0 / hdid
   
   !----- Defining some aliases that will be often used during the integration. -----------!
   kt            = nzg
   nstop         = rk4site%ntext_soil(kt)

   !----- Check whether we are just slightly off. -----------------------------------------!
   slightlymoist = initp%soil_water(kt) > soil8(nstop)%slmsts 
   slightlydry   = initp%soil_water(kt) < soil8(nstop)%soilcp

   !---------------------------------------------------------------------------------------!
   !     If we reached this point, then we may be slightly off-track.  It is very likely   !
   ! that we will need top soil layer temperature and liquid fraction, so we find them     !
   ! now...                                                                                !
   !---------------------------------------------------------------------------------------!
   shctop     = soil8(nstop)%slcpd
   call qwtk8(initp%soil_energy(kt),initp%soil_water(kt)*wdns8,shctop                      &
             ,initp%soil_tempk(kt),initp%soil_fracliq(kt))


   !---------------------------------------------------------------------------------------!
   !      If we access this IF statement, then the soil is slightly dry.  Here we will     !
   ! look for water in adjacent environments and "steal" the amount of water that the top  !
   ! soil layer needs to be exactly at soilcp.                                             !
   !---------------------------------------------------------------------------------------!
   if (slightlydry) then

      !------------------------------------------------------------------------------------!
      !     Now we find how much water we need.  Since we will need to exchange with other !
      ! environments, find it in kg/m2, the standard units.                                !
      !------------------------------------------------------------------------------------!
      water_needed  = (soil8(nstop)%soilcp - initp%soil_water(kt)) * dslz8(kt) * wdns8

      !------------------------------------------------------------------------------------!
      !    First, we check whether there is a temporary surface water.  If so, then we     !
      ! "steal" water from this layer, even if it is not enough.                           !
      !------------------------------------------------------------------------------------!
      if (initp%nlev_sfcwater > 0) then
         sfcwsrc: do
            water_available = initp%sfcwater_mass(1)
            if (water_available > water_needed) then
               !---------------------------------------------------------------------------!
               !    The surface layer has enough water to make up the difference, find the !
               ! energy and depth that will be removed from this layer.                    !
               !---------------------------------------------------------------------------!
               energy_needed = water_needed * initp%sfcwater_energy(1) / water_available
               depth_needed  = water_needed * initp%sfcwater_depth(1)  / water_available

               !---------------------------------------------------------------------------!
               !     Add the water and energy into the top layer, remove it from the surf- !
               ! ace, and quit.                                                            !
               !---------------------------------------------------------------------------!
               initp%soil_water(kt)     = initp%soil_water(kt)                             &
                                        + water_needed * dslzi8(kt) * wdnsi8
               initp%soil_energy(kt)    = initp%soil_energy(kt)                            &
                                        + energy_needed * dslzi8(kt)

               initp%sfcwater_depth(1)  = initp%sfcwater_depth(1)  - depth_needed
               initp%sfcwater_mass(1)   = initp%sfcwater_mass(1)   - water_needed
               initp%sfcwater_energy(1) = initp%sfcwater_energy(1) - energy_needed
               return

            elseif (water_available > 0.d0) then
               !---------------------------------------------------------------------------!
               !     This water will not be enough, but we use all this water, eliminate   !
               ! the layer,  and seek for other sources to fill the remainder.             !
               !---------------------------------------------------------------------------!
               water_needed = water_needed - water_available
               initp%soil_water(kt)     = initp%soil_water(kt)                             &
                                        + initp%sfcwater_mass(1) * dslzi8(kt) * wdnsi8
               initp%soil_energy(kt)    = initp%soil_energy(kt)                            &
                                        + initp%sfcwater_energy(1) * dslzi8(kt)
               !----- If there were more layers, then we move them down. ------------------!
               do kw=1,initp%nlev_sfcwater-1
                  initp%sfcwater_mass  (kw) = initp%sfcwater_mass  (kw+1)
                  initp%sfcwater_energy(kw) = initp%sfcwater_energy(kw+1)
                  initp%sfcwater_depth (kw) = initp%sfcwater_depth (kw+1)
               end do
               
               !----- Remove the top layer. -----------------------------------------------!
               kw=initp%nlev_sfcwater
               initp%sfcwater_mass  (kw) = 0.d0
               initp%sfcwater_energy(kw) = 0.d0
               initp%sfcwater_depth(kw)  = 0.d0
               initp%nlev_sfcwater       = initp%nlev_sfcwater - 1
               !----- If no surface water layer is left, look for another source... -------!
               if (initp%nlev_sfcwater == 0) exit sfcwsrc
            end if
         end do sfcwsrc
      end if

      !------------------------------------------------------------------------------------!
      !     If we hit this point, we aren't done yet and we no longer have temporary surf- !
      ! ace water.  The next candidate is the virtual layer.                               !
      !------------------------------------------------------------------------------------!
      water_available = initp%virtual_water
      if (water_available > water_needed) then
         !---------------------------------------------------------------------------------!
         !    Virtual layer has enough water to solve the problem.  Find the energy        !
         ! associated with the water that will be transferred to the top layer, and move   !
         ! it to there too.                                                                !
         !---------------------------------------------------------------------------------!
         energy_needed = water_needed * initp%virtual_energy / water_available
         depth_needed  = water_needed * initp%virtual_depth  / water_available

         !---------------------------------------------------------------------------------!
         !     Add the water and energy into the top layer, remove it from the surface,    !
         ! and quit.                                                                       !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  + water_needed  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)

         initp%virtual_depth   = initp%virtual_depth  - depth_needed
         initp%virtual_water   = initp%virtual_water  - water_needed
         initp%virtual_energy  = initp%virtual_energy - energy_needed
         return
      elseif (water_available > 0.d0) then
         !---------------------------------------------------------------------------------!
         !     This water will not be enough, but we use it and seek for other sources to  !
         ! fill the remainder.                                                             !
         !---------------------------------------------------------------------------------!
         water_needed = water_needed - water_available
         initp%soil_water(kt)  = initp%soil_water(kt)                                      &
                               + initp%virtual_water * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + initp%virtual_energy * dslzi8(kt)

         !----- Remove the virtual layer. -------------------------------------------------!
         initp%virtual_water   = 0.d0
         initp%virtual_energy  = 0.d0
         initp%virtual_depth   = 0.d0
      end if

      !------------------------------------------------------------------------------------!
      !    If we hit this point the temporary layers were not enough.  We now look for     !
      ! water in the layer immediately beneath the top soil layer.  We now find the amount !
      ! of water this layer can donate.  We will also need energy associated with the , so !
      ! we compute the temperature and liquid fraction.                                    !
      !------------------------------------------------------------------------------------!
      kb              = nzg -1 
      nsbeneath       = rk4site%ntext_soil(kb)
      water_available = (initp%soil_water(kb)-soil8(nsbeneath)%soilcp) * wdns8 * dslz8(kb)
      shcbeneath      = soil8(nsbeneath)%slcpd
      call qwtk8(initp%soil_energy(kb),initp%soil_water(kb)*wdns8,shcbeneath               &
                ,initp%soil_tempk(kb),initp%soil_fracliq(kb))
      !------------------------------------------------------------------------------------!



      if (water_available > water_needed) then
         !---------------------------------------------------------------------------------!
         !    The layer beneath has enough water.  Find the energy associated with this    !
         ! water.                                                                          !
         !---------------------------------------------------------------------------------!
         energy_needed = water_needed * (initp%soil_fracliq(kb) * cliq8                    &
                                        * (initp%soil_tempk(kb) - tsupercool8)             &
                                        + (1.d0 - initp%soil_fracliq(kb)) * cice8          &
                                        * initp%soil_tempk(kb))

         !----- Update water and energy in both layers. -----------------------------------!
         initp%soil_water (kt) = initp%soil_water (kt) + water_needed  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)
         initp%soil_water (kb) = initp%soil_water (kb) - water_needed  * dslzi8(kb)*wdnsi8
         initp%soil_energy(kb) = initp%soil_energy(kb) - energy_needed * dslzi8(kb)
         
         !---------------------------------------------------------------------------------!
         !    We must also update the soil fluxes.                                         !
         !---------------------------------------------------------------------------------!
         initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) + water_needed * hdidi
         initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) - water_needed * hdidi

         !------ Leave the subroutine. ----------------------------------------------------!
         return

      elseif (water_available > 0.d0) then
         !---------------------------------------------------------------------------------!
         !    The water in that layer will not be enough, but we extract all that we can   !
         ! to reduce the amount we still need.                                             ! 
         !---------------------------------------------------------------------------------!
         water_needed = water_needed - water_available

         !---------------------------------------------------------------------------------!
         !    Find the energy associated with the water that will be transferred.          !
         !---------------------------------------------------------------------------------!
         energy_available = water_available * (initp%soil_fracliq(kb) * cliq8              &
                                              * (initp%soil_tempk(kb) - tsupercool8)       &
                                              + (1.d0 - initp%soil_fracliq(kb)) * cice8    &
                                              * initp%soil_tempk(kb))

         initp%soil_water(kt)  = initp%soil_water(kt)                                      &
                               + water_available  * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_available * dslzi8(kt)

         initp%soil_water(kb)  = initp%soil_water(kb)                                      &
                               - water_available  * dslzi8(kb) * wdnsi8
         initp%soil_energy(kb) = initp%soil_energy(kb) - energy_available * dslzi8(kb)
         
         !---------------------------------------------------------------------------------!
         !    We must also update the soil fluxes.                                         !
         !---------------------------------------------------------------------------------!
         initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) + water_available * hdidi
         initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) - water_available * hdidi
         !---------------------------------------------------------------------------------!
      end if


      !------------------------------------------------------------------------------------!
      !     If we hit this point, our final hope is to trap some water from the canopy air !
      ! space.  The water available is not everything above the minimum mixing ratio       !
      ! because we don't want to risk having the canopy air space crashing.                !
      !------------------------------------------------------------------------------------!
      water_available = wcapcan * (initp%can_shv - 5.d0 * rk4min_can_shv)
      if (water_available > water_needed) then
         !---------------------------------------------------------------------------------!
         !    There is enough water vapour. The transfer will require phase change, so the !
         ! energy transfer will be a latent heat flux.  We use the liquid fraction to      !
         ! decide whether it is going to be instant frost or dew (or both).                !
         !---------------------------------------------------------------------------------!
         energy_needed = water_needed * (alvi8 - initp%soil_fracliq(kt) * alli8)

         !---------------------------------------------------------------------------------!
         !     Add the water and energy into the top layer, remove it from the canopy air  !
         ! space and quit.                                                                 !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  + water_needed  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)

         initp%can_shv         = initp%can_shv        - water_needed  * wcapcani

         initp%avg_vapor_gc    = initp%avg_vapor_gc  - water_needed * hdidi
         return

      elseif (water_available > 0.d0) then
         !---------------------------------------------------------------------------------!
         !     This is a critically dry situation and the only reason we won't start cry-  !
         ! ing is because we don't have enough water to waste in tears... As the final     !
         ! desperate act, we will trap any water in excess of the minimum moisture from    !
         ! the canopy air space.  Even if this is a tiny amount, it may be enough to just  !
         ! avoid the model to crash at the sanity check.                                   !
         !---------------------------------------------------------------------------------!
         energy_available = water_available * (alvi8 - initp%soil_fracliq(kt) * alli8)

         !---------------------------------------------------------------------------------!
         !     Add the water and energy into the top layer, remove it from the canopy air  !
         ! space and quit.                                                                 !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)                                      &
                               + water_available  * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_available * dslzi8(kt)

         initp%can_shv         = initp%can_shv        - water_available  * wcapcani
         initp%avg_vapor_gc    = initp%avg_vapor_gc   - water_available  * hdidi

         return
      end if

   !---------------------------------------------------------------------------------------!
   !      If we access this ELSEIF part, then the soil is slightly above saturation.       !
   !---------------------------------------------------------------------------------------!
   elseif (slightlymoist) then
      !------------------------------------------------------------------------------------!
      !     Now we find how much water the top soil layer needs to withdraw.  Since we     !
      ! will need to exchange with other environments, find it in kg/m2, the standard      !
      ! units.  Also, find the total energy that must go away with the water.              !
      !------------------------------------------------------------------------------------!
      water_excess  = (initp%soil_water(kt) - soil8(nstop)%slmsts) * dslz8(kt) * wdns8
      energy_excess = water_excess * ( initp%soil_fracliq(kt) * cliq8                      &
                                     * (initp%soil_tempk(kt) - tsupercool8)                &
                                     + (1.d0 - initp%soil_fracliq(kt)) * cice8             &
                                     * initp%soil_tempk(kt))
      depth_excess  = water_excess * ( initp%soil_fracliq(kt) * wdnsi8                     &
                                     + (1.d0 - initp%soil_fracliq(kt)) * fdnsi8)
      !------------------------------------------------------------------------------------!

      if (initp%nlev_sfcwater > 0) then
         !---------------------------------------------------------------------------------!
         !    If there is already a temporary surface water layer, we simply dump the      !
         ! excess of water in the first layer.                                             !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)     = initp%soil_water(kt)                                   &
                                  - water_excess * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt)    = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
         
         initp%sfcwater_mass(1)   = initp%sfcwater_mass(1)   + water_excess
         initp%sfcwater_energy(1) = initp%sfcwater_energy(1) + energy_excess
         initp%sfcwater_depth(1)  = initp%sfcwater_depth(1)  + depth_excess

         return
      elseif (initp%virtual_water + water_excess > rk4tiny_sfcw_mass) then
         !---------------------------------------------------------------------------------!
         !     If the virtual layer will have some significant mass after adding the water !
         ! excess, we transfer the water to there.  It will likely become the new          !
         ! temporary surface water layer.                                                  !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)     = initp%soil_water(kt)                                   &
                                  - water_excess * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt)    = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
         
         initp%virtual_water      = initp%virtual_water     + water_excess
         initp%virtual_energy     = initp%virtual_energy    + energy_excess
         initp%virtual_depth      = initp%virtual_depth     + depth_excess

         return
      end if

      !------------------------------------------------------------------------------------!
      !    If we hit this point, then the amount of water is small but it can't go to the  !
      ! preferred destinations.  Adding on virtual pool wouldn't help because the water    !
      ! would be sent back to the top soil layer in adjust_sfcw_properties call.  Thus the !
      ! excess now becomes the excess of water in the layer, plus any water left in the    !
      ! virtual layer...                                                                   !
      !    Anyway, we first eliminate any water left in the virtual layer by boiling it to !
      ! the atmosphere.  This is a tiny amount and even if supersaturation occurs, it      !
      ! shouldn't be enough to cause trouble.                                              !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water > rk4eps * rk4eps * rk4tiny_sfcw_mass) then 
         initp%can_shv      = initp%can_shv      + initp%virtual_water  * wcapcani 

         initp%avg_vapor_gc = initp%avg_vapor_gc + initp%virtual_water  * hdidi
         
         !----- Say goodbye to the virtual layer... ---------------------------------------!
         initp%virtual_energy = 0.d0
         initp%virtual_water  = 0.d0
         initp%virtual_depth  = 0.d0
      elseif (initp%virtual_water > 0.d0) then
         !---------------------------------------------------------------------------------!
         !    The amount of water is so small that round-off errors are bound to become    !
         ! more important than the error of not conserving energy and water.  Simply       !
         ! extinguish the virtual layer.                                                   !
         !---------------------------------------------------------------------------------!
         initp%virtual_energy = 0.d0
         initp%virtual_water  = 0.d0
         initp%virtual_depth  = 0.d0
      endif 

      !------------------------------------------------------------------------------------!
      !     Back to the top soil layer, we still need to decide where we should send the   !
      ! excess...  First we try to send to the layer beneath.  This transfer of water      !
      ! implies that some energy is also transferred. Since soil layers are not pure       !
      !  water, we must find the actual amount of energy associated with the water         !
      ! transfer.  So first we compute the temperature and liquid fraction of the layer    !
      ! beneath the top.                                                                   !
      !------------------------------------------------------------------------------------!
      kb         = nzg-1
      nsbeneath  = rk4site%ntext_soil(kb)
      shcbeneath = soil8(nsbeneath)%slcpd
      call qwtk8(initp%soil_energy(kb),initp%soil_water(kb)*wdns8,shcbeneath               &
                ,initp%soil_tempk(kb),initp%soil_fracliq(kb))
      water_room = (soil8(nsbeneath)%slmsts-initp%soil_water(kb)) * wdns8 * dslz8(kb)

      if (water_room > water_excess) then
         !---------------------------------------------------------------------------------!
         !    The layer beneath still has some room for this water excess, send the water  !
         ! down one level.                                                                 !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  - water_excess  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
         initp%soil_water(kb)  = initp%soil_water(kb)  + water_excess  * dslzi8(kb)*wdnsi8
         initp%soil_energy(kb) = initp%soil_energy(kb) + energy_excess * dslzi8(kb)
         !----- Update the fluxes too... --------------------------------------------------!
         initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) - water_excess * hdidi
         initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) + water_excess * hdidi
         !---------------------------------------------------------------------------------!

         return
      elseif (water_room > 0.d0) then
         !---------------------------------------------------------------------------------!
         !   The layer beneath the top layer can't take all water, but we send all that we !
         ! can to that layer so we reduce the amount that will be "boiled".  The remaining !
         ! will go to the canopy air space.  Even if some supersaturation happens, the     !
         ! excess won't be too much to cause the run to crash.                             !
         !---------------------------------------------------------------------------------!
         water_excess  = water_excess  - water_room
         energy_room   = water_room * ( initp%soil_fracliq(kt) * cliq8                     &
                                      * (initp%soil_tempk(kt) - tsupercool8)               &
                                      + (1.d0 - initp%soil_fracliq(kt)) * cice8            &
                                      * initp%soil_tempk(kt))

         !---------------------------------------------------------------------------------!
         !    The layer beneath still has some room for this water excess, send the water  !
         ! down one level.                                                                 !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  - water_room  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) - energy_room * dslzi8(kt)
         initp%soil_water(kb)  = initp%soil_water(kb)  + water_room  * dslzi8(kb)*wdnsi8
         initp%soil_energy(kb) = initp%soil_energy(kb) + energy_room * dslzi8(kb)
         !----- Update the fluxes too... --------------------------------------------------!
         initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) - water_room * hdidi
         initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) + water_room * hdidi
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     The water that is about to leave will do it through "boiling".  Find the    !
         ! latent heat associated with this phase change.                                  !
         !---------------------------------------------------------------------------------!
         energy_excess = water_excess * (alvi8 - initp%soil_fracliq(kt) * alli8)

         !----- Sending the water and energy to the canopy. -------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  - water_excess  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
         initp%can_shv         = initp%can_shv         + water_excess  * wcapcani
         initp%avg_vapor_gc    = initp%avg_vapor_gc    + water_excess  * hdidi
      end if
   end if

   return
end subroutine adjust_topsoil_properties
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will ensure that leaf water is positively defined.  Depending on its  !
! derivative, it can go under zero, in which case we must correct the derivatives rather   !
! than forcing it to be zero.  This guarantees mass conservation.  Likewise, if in the end !
! of the step the leaf water is over the maximum, we remove the excess through shedding.   !
!    After this is checked, we then update the remaining leaf properties, namely the       !
! temperature and liquid water fraction.                                                   !
!------------------------------------------------------------------------------------------!
subroutine adjust_veg_properties(initp,hdid,csite,ipa)
   use rk4_coms             , only : rk4patchtype       & ! structure
                                   , rk4site            & ! intent(in)
                                   , rk4eps             & ! intent(in)
                                   , rk4min_veg_lwater  & ! intent(in)
                                   , rk4min_veg_temp    & ! intent(in)
                                   , rk4max_veg_temp    & ! intent(in)
                                   , hcapcani           & ! intent(in)
                                   , wcapcani           & ! intent(in)
                                   , rk4leaf_drywhc     & ! intent(in)
                                   , rk4leaf_maxwhc     & ! intent(in)
                                   , print_detailed     ! ! intent(in)
   use ed_state_vars        , only : sitetype           & ! structure
                                   , patchtype          ! ! structure
   use ed_misc_coms         , only : fast_diagnostics     ! ! intent(in)
   use consts_coms          , only : cice8              & ! intent(in)
                                   , cliq8              & ! intent(in)
                                   , alvl8              & ! intent(in)
                                   , alvi8              & ! intent(in)
                                   , alli8              & ! intent(in)
                                   , t3ple8             & ! intent(in)
                                   , tsupercool8        & ! intent(in)
                                   , qliqt38            & ! intent(in)
                                   , wdnsi8             & ! intent(in)
                                   , fdnsi8             ! ! intent(in)
   use therm_lib8           , only : qwtk8              ! ! subroutine
   use grid_coms            , only : nzg                ! ! intent(in)
   use soil_coms            , only : dslzi8             ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp  ! Integration buffer
   type(sitetype)         , target     :: csite  ! Current site
   integer                , intent(in) :: ipa    ! Current patch ID
   real(kind=8)           , intent(in) :: hdid   ! Time step 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer    :: cpatch
   integer                             :: ico
   integer                             :: ksn
   real(kind=8)                        :: rk4min_leaf_water
   real(kind=8)                        :: rk4min_wood_water
   real(kind=8)                        :: min_leaf_water
   real(kind=8)                        :: max_leaf_water
   real(kind=8)                        :: min_wood_water
   real(kind=8)                        :: max_wood_water
   real(kind=8)                        :: leaf_wshed
   real(kind=8)                        :: leaf_qwshed
   real(kind=8)                        :: leaf_dwshed
   real(kind=8)                        :: leaf_dew
   real(kind=8)                        :: leaf_qdew
   real(kind=8)                        :: leaf_boil
   real(kind=8)                        :: leaf_qboil
   real(kind=8)                        :: leaf_wshed_tot
   real(kind=8)                        :: leaf_qwshed_tot
   real(kind=8)                        :: leaf_dwshed_tot
   real(kind=8)                        :: leaf_dew_tot
   real(kind=8)                        :: leaf_qdew_tot
   real(kind=8)                        :: leaf_boil_tot
   real(kind=8)                        :: leaf_qboil_tot
   real(kind=8)                        :: wood_wshed
   real(kind=8)                        :: wood_qwshed
   real(kind=8)                        :: wood_dwshed
   real(kind=8)                        :: wood_dew
   real(kind=8)                        :: wood_qdew
   real(kind=8)                        :: wood_boil
   real(kind=8)                        :: wood_qboil
   real(kind=8)                        :: wood_wshed_tot
   real(kind=8)                        :: wood_qwshed_tot
   real(kind=8)                        :: wood_dwshed_tot
   real(kind=8)                        :: wood_dew_tot
   real(kind=8)                        :: wood_qdew_tot
   real(kind=8)                        :: wood_boil_tot
   real(kind=8)                        :: wood_qboil_tot
   real(kind=8)                        :: old_leaf_energy
   real(kind=8)                        :: old_leaf_water
   real(kind=8)                        :: old_leaf_temp
   real(kind=8)                        :: old_leaf_fliq
   real(kind=8)                        :: old_wood_energy
   real(kind=8)                        :: old_wood_water
   real(kind=8)                        :: old_wood_temp
   real(kind=8)                        :: old_wood_fliq
   real(kind=8)                        :: hdidi
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)
   
   !----- Inverse of time increment -------------------------------------------------------!
   hdidi = 1.d0 / hdid

   !----- Initialise the total shedding. --------------------------------------------------!
   leaf_wshed_tot  = 0.d0 
   leaf_qwshed_tot = 0.d0
   leaf_dwshed_tot = 0.d0
   leaf_dew_tot    = 0.d0 
   leaf_qdew_tot   = 0.d0
   leaf_boil_tot   = 0.d0 
   leaf_qboil_tot  = 0.d0
   wood_wshed_tot  = 0.d0 
   wood_qwshed_tot = 0.d0
   wood_dwshed_tot = 0.d0
   wood_dew_tot    = 0.d0 
   wood_qdew_tot   = 0.d0
   wood_boil_tot   = 0.d0 
   wood_qboil_tot  = 0.d0

   !----- Looping over cohorts ------------------------------------------------------------!
   cohortloop: do ico=1,cpatch%ncohorts

      !------------------------------------------------------------------------------------!
      !    Check whether we can solve leaves in this cohort...                             !
      !------------------------------------------------------------------------------------!
      if (initp%leaf_resolvable(ico)) then
         !---------------------------------------------------------------------------------!
         !   Now we find the maximum leaf water possible.                                  !
         !---------------------------------------------------------------------------------!
         rk4min_leaf_water = rk4min_veg_lwater * initp%lai(ico)
         min_leaf_water    = rk4leaf_drywhc    * initp%lai(ico)
         max_leaf_water    = rk4leaf_maxwhc    * initp%lai(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    In case water is to be removed or added, we will need to update the          !
         ! leaf internal energy.  We want to preserve the temperature, though, because     !
         ! this happens as loss of internal energy (shedding) or latent heat (fast         !
         ! dew/boiling).                                                                   !
         !---------------------------------------------------------------------------------!
         call qwtk8(initp%leaf_energy(ico),initp%leaf_water(ico),initp%leaf_hcap(ico)      &
                   ,initp%leaf_temp(ico),initp%leaf_fliq(ico))
         old_leaf_energy = initp%leaf_energy(ico)
         old_leaf_water  = initp%leaf_water (ico)
         old_leaf_temp   = initp%leaf_temp  (ico)
         old_leaf_fliq   = initp%leaf_fliq  (ico)
         !---------------------------------------------------------------------------------!

         if (initp%leaf_water(ico) > max_leaf_water) then

            !------------------------------------------------------------------------------!
            !    Too much water over these leaves, we shall shed the excess to the ground. !
            !------------------------------------------------------------------------------!
            leaf_wshed  = initp%leaf_water(ico) - max_leaf_water

            leaf_qwshed = leaf_wshed                                                       &
                        * ( initp%leaf_fliq(ico) * cliq8                                   &
                          * (initp%leaf_temp(ico) - tsupercool8)                           &
                          + (1.d0-initp%leaf_fliq(ico)) * cice8 * initp%leaf_temp(ico))

            leaf_dwshed = leaf_wshed                                                       &
                        * ( initp%leaf_fliq(ico) * wdnsi8                                  &
                          + (1.d0-initp%leaf_fliq(ico)) * fdnsi8)

            !----- Add the contribution of this cohort to the total shedding. -------------!
            leaf_wshed_tot  = leaf_wshed_tot  + leaf_wshed
            leaf_qwshed_tot = leaf_qwshed_tot + leaf_qwshed
            leaf_dwshed_tot = leaf_dwshed_tot + leaf_dwshed

            !----- Update water mass and energy. ------------------------------------------!
            initp%leaf_water(ico)  = initp%leaf_water(ico)  - leaf_wshed
            initp%leaf_energy(ico) = initp%leaf_energy(ico) - leaf_qwshed

            !----- Update fluxes if needed be. --------------------------------------------!
            if (print_detailed) then
               initp%cfx_qwshed(ico) = initp%cfx_qwshed(ico) + leaf_qwshed * hdidi
            end if
 

         elseif (initp%leaf_water(ico) < min_leaf_water) then
            !------------------------------------------------------------------------------!
            !    If leaf_water is tiny and positive, exchange moisture with the air by     !
            ! donating the total amount as "boiling" (fast evaporation or sublimation).    !
            ! In case the total is tiny but negative, exchange moisture with the air,      !
            ! "stealing" moisture as fast "dew/frost" condensation.                        !
            !------------------------------------------------------------------------------!
            leaf_boil  = max(0.d0,  initp%leaf_water(ico))
            leaf_dew   = max(0.d0,- initp%leaf_water(ico))
            leaf_qboil = leaf_boil * (alvi8 - initp%leaf_fliq(ico) * alli8)
            leaf_qdew  = leaf_dew  * (alvi8 - initp%leaf_fliq(ico) * alli8)


            !----- Add the contribution of this cohort to the total boiling. --------------!
            leaf_boil_tot  = leaf_boil_tot  + leaf_boil
            leaf_dew_tot   = leaf_dew_tot   + leaf_dew
            leaf_qboil_tot = leaf_qboil_tot + leaf_qboil
            leaf_qdew_tot  = leaf_qdew_tot  + leaf_qdew

            !----- Update cohort state variables. -----------------------------------------!
            initp%leaf_water(ico)  = 0.d0
            initp%leaf_energy(ico) = initp%leaf_energy(ico)  + leaf_qdew - leaf_qboil

            !----- Update fluxes if needed be. --------------------------------------------!
            if (print_detailed) then
               initp%cfx_qwflxlc(ico) = initp%cfx_qwflxlc(ico)                             &
                                      + (leaf_qboil - leaf_qdew) * hdidi
            end if
         end if
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Check whether we can solve wood in this cohort...                               !
      !------------------------------------------------------------------------------------!
      if (initp%wood_resolvable(ico)) then
         !---------------------------------------------------------------------------------!
         !   Now we find the maximum leaf water possible.                                  !
         !---------------------------------------------------------------------------------!
         rk4min_wood_water = rk4min_veg_lwater * initp%wai(ico)
         min_wood_water    = rk4leaf_drywhc    * initp%wai(ico)
         max_wood_water    = rk4leaf_maxwhc    * initp%wai(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    In case water is to be removed or added, we will need to update the          !
         ! wood internal energy.  We want to preserve the temperature, though, because     !
         ! this happens as loss of internal energy (shedding) or latent heat (fast         !
         ! dew/boiling).                                                                   !
         !---------------------------------------------------------------------------------!
         call qwtk8(initp%wood_energy(ico),initp%wood_water(ico),initp%wood_hcap(ico)      &
                   ,initp%wood_temp(ico),initp%wood_fliq(ico))
         old_wood_energy = initp%wood_energy(ico)
         old_wood_water  = initp%wood_water (ico)
         old_wood_temp   = initp%wood_temp  (ico)
         old_wood_fliq   = initp%wood_fliq  (ico)
         !---------------------------------------------------------------------------------!

         if (initp%wood_water(ico) > max_wood_water) then

            !------------------------------------------------------------------------------!
            !    Too much water over the wood, we shall shed the excess to the ground.     !
            !------------------------------------------------------------------------------!
            wood_wshed  = initp%wood_water(ico) - max_wood_water

            wood_qwshed = wood_wshed                                                         &
                       * ( initp%wood_fliq(ico) * cliq8                                    &
                         * (initp%wood_temp(ico) - tsupercool8)                            &
                         + (1.d0-initp%wood_fliq(ico)) * cice8 * initp%wood_temp(ico))

            wood_dwshed = wood_wshed                                                         &
                       * ( initp%wood_fliq(ico) * wdnsi8                                   &
                         + (1.d0-initp%wood_fliq(ico)) * fdnsi8)

            !----- Add the contribution of this cohort to the total shedding. -------------!
            wood_wshed_tot  = wood_wshed_tot  + wood_wshed
            wood_qwshed_tot = wood_qwshed_tot + wood_qwshed
            wood_dwshed_tot = wood_dwshed_tot + wood_dwshed

            !----- Update water mass and energy. ------------------------------------------!
            initp%wood_water(ico)  = initp%wood_water(ico)  - wood_wshed
            initp%wood_energy(ico) = initp%wood_energy(ico) - wood_qwshed

            !----- Update fluxes if needed be. --------------------------------------------!
            if (print_detailed) then
               initp%cfx_qwshed(ico) = initp%cfx_qwshed(ico) + wood_qwshed * hdidi
            end if
 

         elseif (initp%wood_water(ico) < min_wood_water) then
            !------------------------------------------------------------------------------!
            !    If wood_water is tiny and positive, exchange moisture with the air by     !
            ! donating the total amount as "boiling" (fast evaporation or sublimation).    !
            ! In case the total is tiny but negative, exchange moisture with the air,      !
            ! "stealing" moisture as fast "dew/frost" condensation.                        !
            !------------------------------------------------------------------------------!
            wood_boil  = max(0.d0,  initp%wood_water(ico))
            wood_dew   = max(0.d0,- initp%wood_water(ico))
            wood_qboil = wood_boil * (alvi8 - initp%wood_fliq(ico) * alli8)
            wood_qdew  = wood_dew  * (alvi8 - initp%wood_fliq(ico) * alli8)


            !----- Add the contribution of this cohort to the total boiling. --------------!
            wood_boil_tot  = wood_boil_tot  + wood_boil
            wood_dew_tot   = wood_dew_tot   + wood_dew
            wood_qboil_tot = wood_qboil_tot + wood_qboil
            wood_qdew_tot  = wood_qdew_tot  + wood_qdew

            !----- Update cohort state variables. -----------------------------------------!
            initp%wood_water(ico)  = 0.d0
            initp%wood_energy(ico) = initp%wood_energy(ico)  + wood_qdew - wood_qboil

            !----- Update fluxes if needed be. --------------------------------------------!
            if (print_detailed) then
               initp%cfx_qwflxwc(ico) = initp%cfx_qwflxwc(ico)                             &
                                      + (wood_qboil - wood_qdew) * hdidi
            end if
         end if
      end if
      !------------------------------------------------------------------------------------!
   end do cohortloop
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    The water that fell from the leaves and branches must go somewhere...  Here we     !
   ! decide which place is the most suitable.  In case there is already a temporary        !
   ! surface water layer we can add the water there, otherwise we dump it into the virtual !
   ! layer, which may or may not become a temporary surface water layer.                   !
   !---------------------------------------------------------------------------------------!
   ksn = initp%nlev_sfcwater
   select case(initp%flag_sfcwater)
   case (0)
      !------ No temporary water, shed the water into the virtual layer. ------------------!
      initp%virtual_water        = initp%virtual_water  + leaf_wshed_tot  + wood_wshed_tot
      initp%virtual_energy       = initp%virtual_energy + leaf_qwshed_tot + wood_qwshed_tot
      initp%virtual_depth        = initp%virtual_depth  + leaf_dwshed_tot + wood_dwshed_tot
      !------------------------------------------------------------------------------------!

   case default
      !------------------------------------------------------------------------------------!
      !     There is a temporary water, shed the excess water to the temporary surface     !
      ! water layer.                                                                       !
      !------------------------------------------------------------------------------------!
      initp%sfcwater_mass(ksn)   = initp%sfcwater_mass(ksn)                                &
                                 + leaf_wshed_tot  + wood_wshed_tot
      initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn)                              &
                                 + leaf_qwshed_tot + wood_qwshed_tot
      initp%sfcwater_depth(ksn)  = initp%sfcwater_depth(ksn)                               &
                                 + leaf_dwshed_tot + wood_dwshed_tot
      !------------------------------------------------------------------------------------!

   end select

   !----- Update the canopy air specific humidity. ----------------------------------------!
   initp%can_shv  = initp%can_shv                                                          &
                  + (leaf_boil_tot + wood_boil_tot - leaf_dew_tot - wood_dew_tot)          &
                  * wcapcani


   !----- Updating output fluxes ----------------------------------------------------------!
   if (fast_diagnostics) then
      initp%avg_wshed_vg  = initp%avg_wshed_vg                                             &
                          + (leaf_wshed_tot + wood_wshed_tot)   * hdidi
      initp%avg_qwshed_vg = initp%avg_qwshed_vg                                            &
                          + (leaf_qwshed_tot + wood_qwshed_tot) * hdidi
      initp%avg_vapor_lc  = initp%avg_vapor_lc                                             &
                          + (leaf_boil_tot - leaf_dew_tot)      * hdidi
      initp%avg_vapor_wc  = initp%avg_vapor_wc                                             &
                          + (wood_boil_tot - wood_dew_tot)      * hdidi
   end if
   if (print_detailed) then
      initp%flx_wshed_vg  = initp%flx_wshed_vg                                             &
                          + (leaf_wshed_tot + wood_wshed_tot)   * hdidi
      initp%flx_qwshed_vg = initp%flx_qwshed_vg                                            &
                          + (leaf_qwshed_tot + wood_qwshed_tot) * hdidi
      initp%flx_vapor_lc  = initp%flx_vapor_lc                                             &
                          + (leaf_boil_tot - leaf_dew_tot)      * hdidi
      initp%flx_vapor_wc  = initp%flx_vapor_wc                                             &
                          + (wood_boil_tot - wood_dew_tot)      * hdidi
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine adjust_veg_properties
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine print_errmax(errmax,yerr,yscal,cpatch,y,ytemp)
   use rk4_coms              , only : rk4patchtype       & ! Structure
                                    , rk4eps             & ! intent(in)
                                    , rk4site            & ! intent(in)
                                    , checkbudget        ! ! intent(in)
   use ed_state_vars         , only : patchtype          ! ! Structure
   use grid_coms             , only : nzg                & ! intent(in)
                                    , nzs                ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target       :: yerr,yscal,y,ytemp
   type(patchtype)    , target       :: cpatch
   real(kind=8)       , intent(out)  :: errmax
   !----- Local variables -----------------------------------------------------------------!
   integer                           :: ico
   integer                           :: k
   logical                           :: troublemaker
   !----- Constants -----------------------------------------------------------------------!
   character(len=28)  , parameter    :: onefmt = '(a16,1x,3(es12.4,1x),11x,l1)'
   character(len=34)  , parameter    :: lyrfmt = '(a16,1x,i6,1x,3(es12.4,1x),11x,l1)'
   character(len=34)  , parameter    :: cohfmt = '(a16,1x,i6,1x,7(es12.4,1x),11x,l1)'
   !----- Functions -----------------------------------------------------------------------!
   logical            , external     :: large_error
   !---------------------------------------------------------------------------------------!


   write(unit=*,fmt='(80a)'    ) ('=',k=1,80)
   write(unit=*,fmt='(a)'      ) '  ..... PRINTING MAXIMUM ERROR INFORMATION: .....'
   write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) 
   write(unit=*,fmt='(a)'      ) ' Patch level variables, single layer:'
   write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
   write(unit=*,fmt='(5(a,1x))')  'Name            ','   Max.Error','   Abs.Error'&
                                &,'       Scale','Problem(T|F)'

   errmax       = max(0.0,abs(yerr%can_lntheta/yscal%can_lntheta))
   troublemaker = large_error(yerr%can_theiv,yscal%can_theiv)
   write(unit=*,fmt=onefmt) 'CAN_LNTHETA:',errmax,yerr%can_lntheta,yscal%can_lntheta       &
                                          ,troublemaker

   errmax       = max(errmax,abs(yerr%can_shv/yscal%can_shv))
   troublemaker = large_error(yerr%can_shv,yscal%can_shv)
   write(unit=*,fmt=onefmt) 'CAN_SHV:',errmax,yerr%can_shv,yscal%can_shv,troublemaker

   errmax = max(errmax,abs(yerr%can_co2/yscal%can_co2))
   troublemaker = large_error(yerr%can_co2,yscal%can_co2)
   write(unit=*,fmt=onefmt) 'CAN_CO2:',errmax,yerr%can_co2,yscal%can_co2,troublemaker

   errmax = max(errmax,abs(yerr%can_prss/yscal%can_prss))
   troublemaker = large_error(yerr%can_prss,yscal%can_prss)
   write(unit=*,fmt=onefmt) 'CAN_PRSS:',errmax,yerr%can_prss,yscal%can_prss,troublemaker

   errmax = max(errmax,abs(yerr%virtual_energy/yscal%virtual_energy))
   troublemaker = large_error(yerr%virtual_energy,yscal%virtual_energy)
   write(unit=*,fmt=onefmt) 'VIRTUAL_ENERGY:',errmax,yerr%virtual_energy                   &
                                             ,yscal%virtual_energy,troublemaker

   errmax = max(errmax,abs(yerr%virtual_water/yscal%virtual_water))
   troublemaker = large_error(yerr%virtual_water,yscal%virtual_water)
   write(unit=*,fmt=onefmt) 'VIRTUAL_WATER:',errmax,yerr%virtual_water,yscal%virtual_water &
                                            ,troublemaker

   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) ' Patch level variables, soil layers:'
   write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'               &
                                &,'   Abs.Error','       Scale','Problem(T|F)'

   do k=rk4site%lsl,nzg
      errmax = max(errmax,abs(yerr%soil_water(k)/yscal%soil_water(k)))
      troublemaker = large_error(yerr%soil_water(k),yscal%soil_water(k))
      write(unit=*,fmt=lyrfmt) 'SOIL_WATER:',k,errmax,yerr%soil_water(k)                   &
                                            ,yscal%soil_water(k),troublemaker

      errmax       = max(errmax,abs(yerr%soil_energy(k)/yscal%soil_energy(k)))
      troublemaker = large_error(yerr%soil_energy(k),yscal%soil_energy(k))
      write(unit=*,fmt=lyrfmt) 'SOIL_ENERGY:',k,errmax,yerr%soil_energy(k)                 &
                                             ,yscal%soil_energy(k),troublemaker
   enddo

   if (yerr%nlev_sfcwater > 0) then
      write(unit=*,fmt='(80a)') ('-',k=1,80)
      write(unit=*,fmt='(a)'  ) 
      write(unit=*,fmt='(80a)') ('-',k=1,80)
      write(unit=*,fmt='(a)'      ) ' Patch level variables, water/snow layers:'
      write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'      &
                                &,'   Abs.Error','       Scale','Problem(T|F)'
      do k=1,yerr%nlev_sfcwater
         errmax       = max(errmax,abs(yerr%sfcwater_energy(k)/yscal%sfcwater_energy(k)))
         troublemaker = large_error(yerr%sfcwater_energy(k),yscal%sfcwater_energy(k))
         write(unit=*,fmt=lyrfmt) 'SFCWATER_ENERGY:',k,errmax,yerr%sfcwater_energy(k)      &
                                                    ,yscal%sfcwater_energy(k),troublemaker

         errmax       = max(errmax,abs(yerr%sfcwater_mass(k)/yscal%sfcwater_mass(k)))
         troublemaker = large_error(yerr%sfcwater_mass(k),yscal%sfcwater_mass(k))
         write(unit=*,fmt=lyrfmt) 'SFCWATER_MASS:',k,errmax,yerr%sfcwater_mass(k)          &
                                                  ,yscal%sfcwater_mass(k),troublemaker
      end do
   end if

   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) ' Leaf-level variables (only the resolvable ones):'
   write(unit=*,fmt='(10(a,1x))')        'Name            ','   PFT','         LAI'        &
                                      ,'         WAI','         WPA','         TAI'        &
                                      ,'   Max.Error','   Abs.Error','       Scale'        &
                                      ,'Problem(T|F)'
   do ico = 1,cpatch%ncohorts
      if (y%leaf_resolvable(ico)) then
         errmax       = max(errmax,abs(yerr%leaf_water(ico)/yscal%leaf_water(ico)))
         troublemaker = large_error(yerr%leaf_water(ico),yscal%leaf_water(ico))
         write(unit=*,fmt=cohfmt) 'LEAF_WATER:',cpatch%pft(ico),y%lai(ico),y%wai(ico)      &
                                               ,y%wpa(ico),y%tai(ico),errmax               &
                                               ,yerr%leaf_water(ico),yscal%leaf_water(ico) &
                                               ,troublemaker
              

         errmax       = max(errmax,abs(yerr%leaf_energy(ico)/yscal%leaf_energy(ico)))
         troublemaker = large_error(yerr%leaf_energy(ico),yscal%leaf_energy(ico))
         write(unit=*,fmt=cohfmt) 'LEAF_ENERGY:',cpatch%pft(ico),cpatch%lai(ico)            &
                                                ,y%wai(ico),y%wpa(ico),y%tai(ico),errmax    &
                                                ,yerr%leaf_energy(ico)                      &
                                                ,yscal%leaf_energy(ico)                     &
                                                ,troublemaker

      end if
   end do

   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) ' Wood-level variables (only the resolvable ones):'
   write(unit=*,fmt='(10(a,1x))')        'Name            ','   PFT','         LAI'        &
                                      ,'         WAI','         WPA','         TAI'        &
                                      ,'   Max.Error','   Abs.Error','       Scale'        &
                                      ,'Problem(T|F)'
   do ico = 1,cpatch%ncohorts
      if (y%wood_resolvable(ico)) then
         errmax       = max(errmax,abs(yerr%wood_water(ico)/yscal%wood_water(ico)))
         troublemaker = large_error(yerr%wood_water(ico),yscal%wood_water(ico))
         write(unit=*,fmt=cohfmt) 'WOOD_WATER:',cpatch%pft(ico),y%lai(ico),y%wai(ico)      &
                                               ,y%wpa(ico),y%tai(ico),errmax               &
                                               ,yerr%wood_water(ico),yscal%wood_water(ico) &
                                               ,troublemaker
              

         errmax       = max(errmax,abs(yerr%wood_energy(ico)/yscal%wood_energy(ico)))
         troublemaker = large_error(yerr%wood_energy(ico),yscal%wood_energy(ico))
         write(unit=*,fmt=cohfmt) 'WOOD_ENERGY:',cpatch%pft(ico),cpatch%lai(ico)            &
                                                ,y%wai(ico),y%wpa(ico),y%tai(ico),errmax    &
                                                ,yerr%wood_energy(ico)                      &
                                                ,yscal%wood_energy(ico)                     &
                                                ,troublemaker

      end if
   end do

   !---------------------------------------------------------------------------------------!
   !     Here we just need to make sure the user is checking mass, otherwise these         !
   ! variables will not be computed.  If this turns out to be essential, we will make this !
   ! permanent and not dependent on checkbudget.  The only one that is not checked is the  !
   ! runoff, because it is computed only after a step is accepted.                         !
   !---------------------------------------------------------------------------------------!
   if (checkbudget) then
      write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
      write(unit=*,fmt='(a)'      ) 
      write(unit=*,fmt='(a)'      ) ' Budget variables, single layer:'
      write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
      write(unit=*,fmt='(5(a,1x))')  'Name            ','   Max.Error','   Abs.Error'      &
                                   &,'       Scale','Problem(T|F)'
      errmax = max(errmax                                                                  &
                  ,abs(yerr%co2budget_loss2atm/yscal%co2budget_loss2atm))
      troublemaker = large_error(yerr%co2budget_loss2atm                                   &
                                ,yscal%co2budget_loss2atm)
      write(unit=*,fmt=onefmt) 'CO2LOSS2ATM:',errmax,yerr%co2budget_loss2atm               &
                              ,yscal%co2budget_loss2atm,troublemaker

      errmax = max(errmax                                                                  &
                  ,abs(yerr%ebudget_loss2atm/yscal%ebudget_loss2atm))
      troublemaker = large_error(yerr%ebudget_loss2atm                                     &
                                ,yscal%ebudget_loss2atm)
      write(unit=*,fmt=onefmt) 'ENLOSS2ATM:',errmax,yerr%ebudget_loss2atm                  &
                              ,yscal%ebudget_loss2atm,troublemaker

      errmax = max(errmax                                                                  &
                  ,abs(yerr%wbudget_loss2atm/yscal%wbudget_loss2atm))
      troublemaker = large_error(yerr%wbudget_loss2atm                                     &
                                ,yscal%wbudget_loss2atm)
      write(unit=*,fmt=onefmt) 'H2OLOSS2ATM:',errmax,yerr%wbudget_loss2atm                 &
                              ,yscal%wbudget_loss2atm,troublemaker

      errmax = max(errmax,abs( yerr%ebudget_loss2drainage                                  &
                             / yscal%ebudget_loss2drainage))
      troublemaker = large_error(yerr%ebudget_loss2drainage                                &
                                ,yscal%ebudget_loss2drainage)
      write(unit=*,fmt=onefmt) 'ENDRAINAGE:',errmax                                        &
                              ,yerr%ebudget_loss2drainage                                  &
                              ,yscal%ebudget_loss2drainage,troublemaker

      errmax = max(errmax,abs( yerr%wbudget_loss2drainage                                  &
                             / yscal%wbudget_loss2drainage))
      troublemaker = large_error(yerr%wbudget_loss2drainage                                &
                                ,yscal%wbudget_loss2drainage)
      write(unit=*,fmt=onefmt) 'H2ODRAINAGE:',errmax                                       &
                              ,yerr%wbudget_loss2drainage                                  &
                              ,yscal%wbudget_loss2drainage,troublemaker

      errmax = max(errmax                                                                  &
                  ,abs(yerr%co2budget_storage/yscal%co2budget_storage))
      troublemaker = large_error(yerr%co2budget_storage                                    &
                                ,yscal%co2budget_storage)
      write(unit=*,fmt=onefmt) 'CO2STORAGE:',errmax,yerr%co2budget_storage                 &
                              ,yscal%co2budget_storage,troublemaker

      errmax = max(errmax                                                                  &
                  ,abs(yerr%ebudget_storage/yscal%ebudget_storage))
      troublemaker = large_error(yerr%ebudget_storage                                      &
                                ,yscal%ebudget_storage)
      write(unit=*,fmt=onefmt) 'ENSTORAGE:',errmax,yerr%ebudget_storage                    &
                              ,yscal%ebudget_storage,troublemaker

      errmax = max(errmax                                                                  &
                  ,abs(yerr%wbudget_storage/yscal%wbudget_storage))
      troublemaker = large_error(yerr%wbudget_storage                                      &
                                ,yscal%wbudget_storage)
      write(unit=*,fmt=onefmt) 'H2OSTORAGE:',errmax,yerr%wbudget_storage                   &
                              ,yscal%wbudget_storage,troublemaker
   end if

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) 

   return
end subroutine print_errmax
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine prints the patch and cohort information when the model falls apart... !
!------------------------------------------------------------------------------------------!
subroutine print_csiteipa(csite, ipa)
   use rk4_coms              , only : rk4site       ! ! intent(in)
   use ed_state_vars         , only : sitetype      & ! structure
                                    , patchtype     ! ! structure
   use ed_misc_coms          , only : current_time  ! ! intent(in)
   use grid_coms             , only : nzs           & ! intent(in)
                                    , nzg           ! ! intent(in)
   use ed_max_dims           , only : n_pft         ! ! intent(in)
   use consts_coms           , only : day_sec       & ! intent(in)
                                    , umol_2_kgC    ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: ipa
   !----- Local variable ------------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   integer                      :: ico
   integer                      :: k
   real                         :: growth_resp
   real                         :: storage_resp
   real                         :: vleaf_resp
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (csite) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,3(i2.2,a))')                                 &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%hour,':',current_time%min,':',current_time%sec,' UTC'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
 
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Leaf information (only the resolvable ones shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH','       BDEAD'   &
                            ,'       BLEAF',' LEAF_ENERGY','   LEAF_TEMP','  LEAF_WATER'
   do ico = 1,cpatch%ncohorts
      if (cpatch%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%lai(ico),cpatch%dbh(ico),cpatch%bdead(ico)        &
              ,cpatch%bleaf(ico),cpatch%leaf_energy(ico),cpatch%leaf_temp(ico)             &
              ,cpatch%leaf_water(ico)
      end if
   end do
   write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','     FS_OPEN','         FSW','         FSN'   &
                            ,'         GPP','   LEAF_RESP'
   do ico = 1,cpatch%ncohorts
      if (cpatch%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%lai(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)         &
              ,cpatch%gpp(ico),cpatch%leaf_respiration(ico)
      end if
   end do
   write (unit=*,fmt='(2(a7,1x),5(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','   ROOT_RESP',' GROWTH_RESP','   STOR_RESP'   &
                            ,'  VLEAF_RESP'
   do ico = 1,cpatch%ncohorts
      if (cpatch%leaf_resolvable(ico)) then
         growth_resp  = cpatch%growth_respiration(ico)  * cpatch%nplant(ico)               &
                      / (day_sec * umol_2_kgC)
         storage_resp = cpatch%storage_respiration(ico) * cpatch%nplant(ico)               &
                      / (day_sec * umol_2_kgC)
         vleaf_resp   = cpatch%vleaf_respiration(ico)  * cpatch%nplant(ico)                &
                      / (day_sec * umol_2_kgC)

         write(unit=*,fmt='(2(i7,1x),5(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%lai(ico),cpatch%root_respiration(ico),growth_resp,storage_resp       &
              ,vleaf_resp
      end if
   end do
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Wood information (only the resolvable ones shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','         WAI','         DBH','       BDEAD'   &
                            ,'    BSAPWOOD',' WOOD_ENERGY','   WOOD_TEMP','  WOOD_WATER'
   do ico = 1,cpatch%ncohorts
      if (cpatch%wood_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%wai(ico),cpatch%dbh(ico),cpatch%bdead(ico)        &
              ,cpatch%bsapwood(ico),cpatch%wood_energy(ico),cpatch%wood_temp(ico)          &
              ,cpatch%wood_water(ico)
      end if
   end do
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(8(a12,1x))')  '   DIST_TYPE','         AGE','        AREA'          &
                                    ,'          RH','      CWD_RH','AVGDAILY_TMP'          &
                                    ,'     SUM_CHD','     SUM_DGD'
   write (unit=*,fmt='(i12,1x,7(es12.4,1x))')  csite%dist_type(ipa),csite%age(ipa)         &
         ,csite%area(ipa),csite%rh(ipa),csite%cwd_rh(ipa),csite%avg_daily_temp(ipa)        &
         ,csite%sum_chd(ipa),csite%sum_dgd(ipa)

   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','VEG_DISPLACE'          &
                                    ,'         LAI','        HTRY','    CAN_RHOS'          &
                                    ,'   CAN_DEPTH'
   write (unit=*,fmt='(7(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)          &
                                      ,csite%veg_displace(ipa),csite%lai(ipa)              &
                                      ,csite%htry(ipa),csite%can_rhos(ipa)                 &
                                      ,csite%can_depth(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(6(a12,1x))')  '   CAN_THEIV','    CAN_TEMP','     CAN_SHV'          &
                                    ,'    CAN_PRSS','     CAN_CO2','       GGNET'
   write (unit=*,fmt='(6(es12.4,1x))') csite%can_theiv(ipa),csite%can_temp(ipa)            &
                                      ,csite%can_shv(ipa)  ,csite%can_prss(ipa)            &
                                      ,csite%can_co2(ipa)  ,csite%ggnet   (ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(9(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                    ,'       TSTAR','        ZETA','     RI_BULK'          &
                                    ,'     RLONG_G','    RSHORT_G','     RLONG_S'
   write (unit=*,fmt='(9(es12.4,1x))') csite%ustar(ipa),csite%qstar(ipa),csite%cstar(ipa)  &
                                      ,csite%tstar(ipa),csite%zeta(ipa),csite%ribulk(ipa)  &
                                      ,csite%rlong_g(ipa),csite%rshort_g(ipa)              &
                                      ,csite%rlong_s(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,a12)') '  PFT','       REPRO'
   do k=1,n_pft
      write (unit=*,fmt='(i5,1x,es12.4)') k,csite%repro(k,ipa)
   end do

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4site%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,rk4site%ntext_soil(k)                &
            ,csite%soil_energy(k,ipa),csite%soil_tempk(k,ipa),csite%soil_water(k,ipa)      &
            ,csite%soil_fracliq(k,ipa)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,6(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'       &
                                      &,'    RSHORT_S'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,6(es12.4,1x))') k,csite%sfcwater_energy(k,ipa)          &
               ,csite%sfcwater_tempk(k,ipa),csite%sfcwater_mass(k,ipa)                     &
               ,csite%sfcwater_fracliq(k,ipa),csite%sfcwater_depth(k,ipa)                  &
               ,csite%rshort_s(k,ipa)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '
   return
end subroutine print_csiteipa
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is similar to print_csite, except that it also prints the             !
! outcome of the Runge-Kutta integrator.                                                   !
!------------------------------------------------------------------------------------------!
subroutine print_rk4patch(y,csite,ipa)
   use rk4_coms              , only : rk4patchtype          & ! structure
                                    , rk4site               & ! intent(in)
                                    , rk4tiny_sfcw_mass     ! ! intent(in)
   use ed_state_vars         , only : sitetype              & ! structure
                                    , patchtype             ! ! structure
   use grid_coms             , only : nzg                   & ! intent(in)
                                    , nzs                   ! ! intent(in)
   use ed_misc_coms          , only : current_time          ! ! intent(in)
   use consts_coms           , only : pio1808               ! ! intent(in)
   use therm_lib8            , only : qtk8                  & ! subroutine
                                    , qwtk8                 ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: y
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer    :: cpatch
   integer                         :: k
   integer                         :: ico
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (rk4patch) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,3(i2.2,a))')                                 &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%hour,':',current_time%min,':',current_time%sec,' UTC'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(80a)')         ('-',k=1,80)
   write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
   write (unit=*,fmt='(a,1x,es12.4)') ' Longitude             : ',rk4site%lon
   write (unit=*,fmt='(a,1x,es12.4)') ' Latitude              : ',rk4site%lat
   write (unit=*,fmt='(a,1x,es12.4)') ' Air temperature       : ',rk4site%atm_tmp
   write (unit=*,fmt='(a,1x,es12.4)') ' Air potential temp.   : ',rk4site%atm_theta
   write (unit=*,fmt='(a,1x,es12.4)') ' Air theta_Eiv         : ',rk4site%atm_theiv
   write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov mixing ratio     : ',rk4site%atm_shv
   write (unit=*,fmt='(a,1x,es12.4)') ' CO2  mixing ratio     : ',rk4site%atm_co2
   write (unit=*,fmt='(a,1x,es12.4)') ' Pressure              : ',rk4site%atm_prss
   write (unit=*,fmt='(a,1x,es12.4)') ' Exner function        : ',rk4site%atm_exner
   write (unit=*,fmt='(a,1x,es12.4)') ' Wind speed            : ',rk4site%vels
   write (unit=*,fmt='(a,1x,es12.4)') ' Height                : ',rk4site%geoht
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. mass  flux    : ',rk4site%pcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. heat  flux    : ',rk4site%qpcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. depth flux    : ',rk4site%dpcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Downward SW radiation : ',rk4site%rshort
   write (unit=*,fmt='(a,1x,es12.4)') ' Downward LW radiation : ',rk4site%rlong
   write (unit=*,fmt='(a,1x,es12.4)') ' Zenith angle (deg)    : ',acos(rk4site%cosz)       &
                                                                 / pio1808

   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Leaf information (only those resolvable are shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','      HEIGHT','         DBH','       BDEAD'   &
                            ,'       BLEAF','     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdead(ico)       &
              ,cpatch%bleaf (ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),7(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','         GPP','   LEAF_RESP','   ROOT_RESP'   &
                            ,' GROWTH_RESP','   STOR_RESP','  VLEAF_RESP'
   do ico = 1,cpatch%ncohorts
      if (cpatch%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),7(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,y%lai(ico),y%gpp(ico),y%leaf_resp(ico),y%root_resp(ico),y%growth_resp(ico)  &
              ,y%storage_resp(ico),y%vleaf_resp(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),9(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','         WAI','         TAI',' LEAF_ENERGY'   &
             ,'  LEAF_WATER','   LEAF_HCAP','   LEAF_TEMP','   LEAF_FLIQ','    LINT_SHV'
   do ico = 1,cpatch%ncohorts
      if (y%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
               ,y%lai(ico),y%wai(ico),y%tai(ico),y%leaf_energy(ico),y%leaf_water(ico)      &
               ,y%leaf_hcap(ico),y%leaf_temp(ico),y%leaf_fliq(ico),y%lint_shv(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
              '    PFT','KRDEPTH','         LAI','      HEIGHT','   LEAF_TEMP'             &
                  ,'    VEG_WIND','  LEAF_REYNO','LEAF_GRASHOF',' LEAF_NUFORC'             &
                  ,' LEAF_NUFREE'
   do ico = 1,cpatch%ncohorts
      if (y%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico),cpatch%krdepth(ico)   &
               ,y%lai(ico),cpatch%hite(ico),y%leaf_temp(ico),y%veg_wind(ico)               &
               ,y%leaf_reynolds(ico),y%leaf_grashof(ico),y%leaf_nussforc(ico)              &
               ,y%leaf_nussfree(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),7(a12,1x))')                                               &
              '    PFT','KRDEPTH','         LAI','      HEIGHT','    LEAF_GBH'             &
                  ,'    LEAF_GBW','  GSW_CLOSED','    GSW_OPEN','     FS_OPEN'
   do ico = 1,cpatch%ncohorts
      if (y%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),7(es12.4,1x))') cpatch%pft(ico),cpatch%krdepth(ico)   &
               ,y%lai(ico),cpatch%hite(ico),y%leaf_gbh(ico),y%leaf_gbw(ico)                &
               ,y%gsw_closed(ico),y%gsw_open(ico),cpatch%fs_open(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                               &
              '    PFT','KRDEPTH','         LAI','      HEIGHT','    RSHORT_L'             &
                  ,'     RLONG_L','  PAR_L_BEAM','  PAR_L_DIFF'
   do ico = 1,cpatch%ncohorts
      if (y%leaf_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))') cpatch%pft(ico),cpatch%krdepth(ico)   &
               ,y%lai(ico),cpatch%hite(ico),y%rshort_l(ico),y%rlong_l(ico)                 &
               ,cpatch%par_l_beam(ico),cpatch%par_l_diffuse(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Wood information (only those resolvable are shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),5(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','      HEIGHT','         DBH','       BDEAD'   &
                            ,'    BSAPWOOD'
   do ico = 1,cpatch%ncohorts
      if (cpatch%wood_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),5(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdead(ico)       &
              ,cpatch%bsapwood(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','         WAI','         TAI',' WOOD_ENERGY'   &
             ,'  WOOD_WATER','   WOOD_HCAP','   WOOD_TEMP','   WOOD_FLIQ'
   do ico = 1,cpatch%ncohorts
      if (y%wood_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
               ,y%lai(ico),y%wai(ico),y%tai(ico),y%wood_energy(ico),y%wood_water(ico)      &
               ,y%wood_hcap(ico),y%wood_temp(ico),y%wood_fliq(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
              '    PFT','KRDEPTH','         WAI','      HEIGHT','   LEAF_TEMP'             &
                  ,'    VEG_WIND','  WOOD_REYNO','WOOD_GRASHOF',' WOOD_NUFORC'             &
                  ,' WOOD_NUFREE'
   do ico = 1,cpatch%ncohorts
      if (y%wood_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico),cpatch%krdepth(ico)   &
               ,y%wai(ico),cpatch%hite(ico),y%wood_temp(ico),y%veg_wind(ico)               &
               ,y%wood_reynolds(ico),y%wood_grashof(ico),y%wood_nussforc(ico)              &
               ,y%wood_nussfree(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                               &
              '    PFT','KRDEPTH','         WAI','      HEIGHT','    WOOD_GBH'             &
                  ,'    WOOD_GBW','    RSHORT_W','     RLONG_W'
   do ico = 1,cpatch%ncohorts
      if (y%wood_resolvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))') cpatch%pft(ico),cpatch%krdepth(ico)   &
               ,y%wai(ico),cpatch%hite(ico),y%wood_gbh(ico),y%wood_gbw(ico)                &
               ,y%rshort_w(ico),y%rlong_w(ico) 
      end if
   end do
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(8(a12,1x))')   '  VEG_HEIGHT','   VEG_ROUGH','VEG_DISPLACE'         &
                                     ,'   PATCH_LAI','   CAN_DEPTH','     CAN_CO2'         &
                                     ,'    CAN_PRSS','       GGNET'
                                     
   write (unit=*,fmt='(8(es12.4,1x))') y%veg_height,y%veg_rough,y%veg_displace             &
                                      ,csite%lai(ipa),y%can_depth,y%can_co2,y%can_prss     &
                                      ,y%ggnet
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(8(a12,1x))')  '    CAN_RHOS','   CAN_THEIV','   CAN_THETA'          &
                                    ,'    CAN_TEMP','     CAN_SHV','     CAN_SSH'          &
                                    ,'    CAN_RVAP','     CAN_RHV'
                                     
                                     
   write (unit=*,fmt='(8(es12.4,1x))')   y%can_rhos , y%can_theiv, y%can_theta             &
                                       , y%can_temp , y%can_shv  , y%can_ssh               &
                                       , y%can_rvap , y%can_rhv
                                       

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                    ,'       TSTAR','       ESTAR','        ZETA'          &
                                    ,'     RI_BULK'
   write (unit=*,fmt='(7(es12.4,1x))') y%ustar,y%qstar,y%cstar,y%tstar,y%estar,y%zeta      &
                                      ,y%ribulk

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(2(a12,1x))')  '          RH','      CWD_RH'
   write (unit=*,fmt='(2(es12.4,1x))') y%rh,y%cwd_rh

   write (unit=*,fmt='(80a)') ('-',k=1,80)


   write (unit=*,fmt='(5(a12,1x))')  '   FLAG_SFCW',' VIRT_ENERGY','  VIRT_WATER'          &
                                    ,'VIRTUAL_TEMP','VIRTUAL_FLIQ'
   write (unit=*,fmt='(i12,1x,4(es12.4,1x))') y%flag_sfcwater,y%virtual_energy             &
                                             ,y%virtual_water,y%virtual_tempk              &
                                             ,y%virtual_fracliq
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(4(a12,1x))')    '  GROUND_SHV','  GROUND_SSH',' GROUND_TEMP'        &
                                      ,' GROUND_FLIQ'
   write (unit=*,fmt='(4(es12.4,1x))') y%ground_shv, y%ground_ssh, y%ground_temp           &
                                      ,y%ground_fliq

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4site%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,rk4site%ntext_soil(k)                &
            ,y%soil_energy(k),y%soil_tempk(k),y%soil_water(k),y%soil_fracliq(k)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,5(es12.4,1x))') k,y%sfcwater_energy(k)                  &
               ,y%sfcwater_tempk(k),y%sfcwater_mass(k),y%sfcwater_fracliq(k)               &
               ,y%sfcwater_depth(k)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '

   !----- Printing the corresponding patch information (with some redundancy) -------------!
   call print_csiteipa(csite, ipa)

   call fatal_error('IFLAG1 problem. The model didn''t converge!','print_rk4patch'&
                 &,'rk4_integ_utils.f90')
   return
end subroutine print_rk4patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine prints the full state of a given patch, for full debugging          !
! purposes.  This will create one file for each patch.  This sub-routine will not print    !
! the temperature of each cohort, instead it will just compute the average.                !
!------------------------------------------------------------------------------------------!
subroutine print_rk4_state(initp,fluxp,csite,ipa,elapsed,hdid)
   use consts_coms  , only : t3ple8        ! ! intent(in)
   use ed_max_dims  , only : str_len       ! ! intent(in)
   use ed_misc_coms , only : current_time  ! ! intent(in)
   use ed_state_vars, only : sitetype      & ! structure
                           , patchtype     ! ! structure
   use grid_coms    , only : nzg           & ! intent(in)
                           , nzs           ! ! intent(in)
   use rk4_coms     , only : rk4patchtype  & ! structure
                           , rk4site       & ! intent(in)
                           , detail_pref   ! ! intent(in)
   use therm_lib8   , only : qwtk8         ! ! sub-routine
   use soil_coms    , only : soil8         ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(rk4patchtype)    , target     :: initp
   type(rk4patchtype)    , target     :: fluxp
   type(sitetype)        , target     :: csite
   integer               , intent(in) :: ipa
   real(kind=8)          , intent(in) :: elapsed
   real(kind=8)          , intent(in) :: hdid
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch
   type(patchtype)       , pointer    :: jpatch
   character(len=str_len)             :: detail_fout
   integer                            :: k
   integer                            :: jpa
   integer                            :: nsoil
   integer                            :: ico
   integer                            :: jco
   integer                            :: leaf_resolve
   integer                            :: wood_resolve
   logical                            :: isthere
   real(kind=8)                       :: sum_leaf_energy
   real(kind=8)                       :: sum_leaf_water
   real(kind=8)                       :: sum_leaf_hcap
   real(kind=8)                       :: sum_wood_energy
   real(kind=8)                       :: sum_wood_water
   real(kind=8)                       :: sum_wood_hcap
   real(kind=8)                       :: sum_gpp
   real(kind=8)                       :: sum_lai
   real(kind=8)                       :: sum_wai
   real(kind=8)                       :: sum_plresp
   real(kind=8)                       :: soil_rh
   real(kind=8)                       :: qintercepted
   real(kind=8)                       :: avg_leaf_temp
   real(kind=8)                       :: avg_leaf_fliq
   real(kind=8)                       :: avg_wood_temp
   real(kind=8)                       :: avg_wood_fliq
   real(kind=8)                       :: sfc_temp
   real(kind=8)                       :: elapsec
   !----- Local constants. ----------------------------------------------------------------!
   character(len=10), parameter :: phfmt='(74(a,1x))'
   character(len=48), parameter :: pbfmt='(3(i13,1x),4(es13.6,1x),3(i13,1x),64(es13.6,1x))'
   character(len=10), parameter :: chfmt='(57(a,1x))'
   character(len=48), parameter :: cbfmt='(3(i13,1x),2(es13.6,1x),3(i13,1x),49(es13.6,1x))'
   !----- Locally saved variables. --------------------------------------------------------!
   logical          , save      :: first_time=.true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     First time here.  Delete all files.                                               !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      do jpa = 1, csite%npatches
         !---------------------------------------------------------------------------------!
         ! Patch level files.                                                              !
         !---------------------------------------------------------------------------------!
         write (detail_fout,fmt='(2a,i4.4,a)') trim(detail_pref),'prk4_patch_',jpa,'.txt'
         inquire(file=trim(detail_fout),exist=isthere)
         if (isthere) then
            !---- Open the file to delete when closing. -----------------------------------!
            open (unit=83,file=trim(detail_fout),status='old',action='write')
            close(unit=83,status='delete')
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         ! Cohort level files.                                                             !
         !---------------------------------------------------------------------------------!
         jpatch => csite%patch(jpa)
         do jco = 1, jpatch%ncohorts
            write (detail_fout,fmt='(a,2(a,i4.4),a)')                                      &
                  trim(detail_pref),'crk4_patch_',jpa,'_cohort_',jco,'.txt'
            inquire(file=trim(detail_fout),exist=isthere)
            if (isthere) then
               !---- Open the file to delete when closing. --------------------------------!
               open (unit=84,file=trim(detail_fout),status='old',action='write')
               close(unit=84,status='delete')
            end if
         end do
         !---------------------------------------------------------------------------------!
      end do
      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     First we loop over all cohorts and add the vegetation energy and water.           !
   !---------------------------------------------------------------------------------------!
   sum_leaf_energy = 0.d0
   sum_leaf_water  = 0.d0
   sum_leaf_hcap   = 0.d0
   sum_wood_energy = 0.d0
   sum_wood_water  = 0.d0
   sum_wood_hcap   = 0.d0
   sum_gpp         = 0.d0
   sum_plresp      = 0.d0
   sum_lai         = 0.d0
   sum_wai         = 0.d0
   soil_rh         = initp%rh-initp%cwd_rh
   cpatch => csite%patch(ipa)
   do ico=1,cpatch%ncohorts
      if (initp%leaf_resolvable(ico)) then
         !----- Integrate vegetation properties using m2gnd rather than plant. ------------!
         sum_leaf_energy = sum_leaf_energy + initp%leaf_energy(ico)
         sum_leaf_water  = sum_leaf_water  + initp%leaf_water(ico)
         sum_leaf_hcap   = sum_leaf_hcap   + initp%leaf_hcap(ico)
         sum_gpp         = sum_gpp         + initp%gpp(ico)
         sum_plresp      = sum_plresp      + initp%leaf_resp(ico)                          &
                                           + initp%root_resp(ico)                          &
                                           + initp%growth_resp(ico)                        &
                                           + initp%storage_resp(ico)                       &
                                           + initp%vleaf_resp(ico)
      end if
      if (initp%wood_resolvable(ico)) then
         !----- Integrate vegetation properties using m2gnd rather than plant. ------------!
         sum_wood_energy = sum_wood_energy + initp%wood_energy(ico)
         sum_wood_water  = sum_wood_water  + initp%wood_water(ico)
         sum_wood_hcap   = sum_wood_hcap   + initp%wood_hcap(ico)
      end if

      !----- TAI.  We integrate all cohorts, including those that we skip. ----------------!
      sum_lai = sum_lai + initp%lai(ico)
      sum_wai = sum_wai + initp%wai(ico)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Then we find the average leaf temperature.  If none of the cohorts were solved,   !
   ! of if there is no vegetation, we assign the canopy air temperature.                   !
   !---------------------------------------------------------------------------------------!
   if (sum_leaf_energy == 0.d0) then
      avg_leaf_temp = initp%can_temp
      if (initp%can_temp == t3ple8) then
         avg_leaf_fliq = 5.d-1
      elseif (initp%can_temp > t3ple8) then
         avg_leaf_fliq = 1.d0
      else
         avg_leaf_fliq = 0.d0
      end if
   else
      call qwtk8(sum_leaf_energy,sum_leaf_water,sum_leaf_hcap,avg_leaf_temp,avg_leaf_fliq)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Then we find the average wood temperature.  If none of the cohorts were solved,   !
   ! of if there is no vegetation, we assign the canopy air temperature.                   !
   !---------------------------------------------------------------------------------------!
   if (sum_wood_energy == 0.d0) then
      avg_wood_temp = initp%can_temp
      if (initp%can_temp == t3ple8) then
         avg_wood_fliq = 5.d-1
      elseif (initp%can_temp > t3ple8) then
         avg_wood_fliq = 1.d0
      else
         avg_wood_fliq = 0.d0
      end if
   else
      call qwtk8(sum_wood_energy,sum_wood_water,sum_wood_hcap,avg_wood_temp,avg_wood_fliq)
   end if
   !---------------------------------------------------------------------------------------!


   !----- Compute the hour as elapsed seconds since midnight. -----------------------------!
   elapsec = dble(current_time%time) + elapsed
   !---------------------------------------------------------------------------------------!


   !----- Find the soil type of the top layer. --------------------------------------------!
   nsoil   = rk4site%ntext_soil(nzg)
   !---------------------------------------------------------------------------------------!



   !----- Create the file name. -----------------------------------------------------------!
   write (detail_fout,fmt='(2a,i4.4,a)') trim(detail_pref),'prk4_patch_',ipa,'.txt'
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check whether the file exists or not.  In case it doesn't, create it and add the   !
   ! header.                                                                               !
   !---------------------------------------------------------------------------------------!
   inquire(file=trim(detail_fout),exist=isthere)
   if (.not. isthere) then
      open  (unit=83,file=trim(detail_fout),status='replace',action='write')
      write (unit=83,fmt=phfmt)  '         YEAR', '        MONTH', '          DAY'         &
                               , '         TIME', '         HDID', '          LAI'         &
                               , '          WAI', '          KSN', 'FLAG.SFCWATER'         &
                               , '  FLAG.WFLXGC', '     ATM.PRSS', '     ATM.TEMP'         &
                               , '      ATM.SHV', '      ATM.CO2', '     ATM.VELS'         &
                               , '    ATM.PRATE', '   ATM.HEIGHT', '     ATM.RHOS'         &
                               , '   ATM.RELHUM', '    ATM.THETA', '    ATM.THEIV'         &
                               , '   MET.RSHORT', '    MET.RLONG', '     CAN.PRSS'         &
                               , '     CAN.TEMP', '      CAN.SHV', '      CAN.CO2'         &
                               , '    CAN.DEPTH', '     CAN.RHOS', '   CAN.RELHUM'         &
                               , '    CAN.THETA', '    CAN.THEIV', '     SFC.TEMP'         &
                               , '      SFC.SHV', '    LEAF.TEMP', '   LEAF.WATER'         &
                               , '    WOOD.TEMP', '   WOOD.WATER', '       GGBARE'         &
                               , '        GGVEG', '        GGNET', '      OPENCAN'         &
                               , '    SOIL.TEMP', '   SOIL.WATER', '       SOILCP'         &
                               , '       SOILWP', '       SOILFC', '       SLMSTS'         &
                               , '        USTAR', '        TSTAR', '        QSTAR'         &
                               , '        CSTAR', '         ZETA', '      RI_BULK'         &
                               , '   GND.RSHORT', '    GND.RLONG', '       WFLXLC'         &
                               , '       WFLXWC', '       DEWGND', '       WFLXGC'         &
                               , '       WFLXAC', '       TRANSP', '        WSHED'         &
                               , '    INTERCEPT', '  THROUGHFALL', '       HFLXGC'         &
                               , '       HFLXLC', '       HFLXWC', '       HFLXAC'         &
                               , '       CFLXAC', '        CWDRH', '       SOILRH'         &
                               , '          GPP', '       PLRESP'
                               
                               
                               
                               
      close (unit=83,status='keep')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Re-open the file at the last line, and include the current status.                !
   !---------------------------------------------------------------------------------------!
   open (unit=83,file=trim(detail_fout),status='old',action='write',position='append')
   write(unit=83,fmt=pbfmt)                                                                &
                     current_time%year     , current_time%month    , current_time%date     &
                   , elapsec               , hdid                  , sum_lai               &
                   , sum_wai               , initp%nlev_sfcwater   , initp%flag_sfcwater   &
                   , initp%flag_wflxgc     , rk4site%atm_prss      , rk4site%atm_tmp       &
                   , rk4site%atm_shv       , rk4site%atm_co2       , rk4site%vels          &
                   , rk4site%pcpg          , rk4site%geoht         , rk4site%atm_rhos      &
                   , rk4site%atm_rhv       , rk4site%atm_theta     , rk4site%atm_theiv     &
                   , rk4site%rshort        , rk4site%rlong         , initp%can_prss        &
                   , initp%can_temp        , initp%can_shv         , initp%can_co2         &
                   , initp%can_depth       , initp%can_rhos        , initp%can_rhv         &
                   , initp%can_theta       , initp%can_theiv       , initp%ground_temp     &
                   , initp%ground_shv      , avg_leaf_temp         , sum_leaf_water        &
                   , avg_wood_temp         , sum_wood_water        , initp%ggbare          &
                   , initp%ggveg           , initp%ggnet           , initp%opencan_frac    &
                   , initp%soil_tempk(nzg) , initp%soil_water(nzg) , soil8(nsoil)%soilcp   &
                   , soil8(nsoil)%soilwp   , soil8(nsoil)%sfldcap  , soil8(nsoil)%slmsts   &
                   , initp%ustar           , initp%tstar           , initp%qstar           &
                   , initp%cstar           , initp%zeta            , initp%ribulk          &
                   , fluxp%flx_rshort_gnd  , fluxp%flx_rlong_gnd   , fluxp%flx_vapor_lc    &
                   , fluxp%flx_vapor_wc    , fluxp%flx_dew_cg      , fluxp%flx_vapor_gc    &
                   , fluxp%flx_vapor_ac    , fluxp%flx_transp      , fluxp%flx_wshed_vg    &
                   , fluxp%flx_intercepted , fluxp%flx_throughfall , fluxp%flx_sensible_gc &
                   , fluxp%flx_sensible_lc , fluxp%flx_sensible_wc , fluxp%flx_sensible_ac &
                   , fluxp%flx_carbon_ac   , initp%cwd_rh          , soil_rh               &
                   , sum_gpp               , sum_plresp
                   
                   
   close(unit=83,status='keep')
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Now it is time to check and output the cohort-level files.                        !
   !---------------------------------------------------------------------------------------!
   do ico=1, cpatch%ncohorts

      !----- Find the integer version of "resolvable". ------------------------------------!
      if (initp%leaf_resolvable(ico)) then
         leaf_resolve = 1
      else
         leaf_resolve = 0
      end if
      if (initp%wood_resolvable(ico)) then
         wood_resolve = 1
      else
         wood_resolve = 0
      end if

      !----- Copy intercepted water. ------------------------------------------------------!
      qintercepted = fluxp%cfx_qintercepted(ico)

      !----- Create the file name. --------------------------------------------------------!
      write (detail_fout,fmt='(a,2(a,i4.4),a)')                                            &
                                trim(detail_pref),'crk4_patch_',ipa,'_cohort_',ico,'.txt'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check whether the file exists or not.  In case it doesn't, create it and add    !
      ! the header.                                                                        !
      !------------------------------------------------------------------------------------!
      inquire(file=trim(detail_fout),exist=isthere)
      if (.not. isthere) then
         open  (unit=84,file=trim(detail_fout),status='replace',action='write')
         write (unit=84,fmt=chfmt)                                                         &
                                 '         YEAR', '        MONTH', '          DAY'         &
                               , '         TIME', '         HDID', '          PFT'         &
                               , ' LEAF_RESOLVE', ' WOOD_RESOLVE', '       NPLANT'         &
                               , '       HEIGHT', '          LAI', '          WAI'         &
                               , '   CROWN_AREA', '  LEAF_ENERGY', '   LEAF_WATER'         &
                               , '    LEAF_HCAP', '    LEAF_TEMP', '    LEAF_FLIQ'         &
                               , '  WOOD_ENERGY', '   WOOD_WATER', '    WOOD_HCAP'         &
                               , '    WOOD_TEMP', '    WOOD_FLIQ', '     VEG_WIND'         &
                               , '      FS_OPEN', '   LEAF_REYNO', ' LEAF_GRASHOF'         &
                               , '  LEAF_NUFREE', '  LEAF_NUFORC', '   WOOD_REYNO'         &
                               , ' WOOD_GRASHOF', '  WOOD_NUFREE', '  WOOD_NUFORC'         &
                               , '     LINT_SHV', '     LEAF_GBH', '     LEAF_GBW'         &
                               , '     WOOD_GBH', '     WOOD_GBW', '     GSW_OPEN'         &
                               , '     GSW_CLOS', '          GPP', '    LEAF_RESP'         &
                               , '    ROOT_RESP', '  GROWTH_RESP', ' STORAGE_RESP'         &
                               , '   VLEAF_RESP', '     RSHORT_L', '      RLONG_L'         &
                               , '     RSHORT_W', '      RLONG_W', '       HFLXLC'         &
                               , '       HFLXWC', '      QWFLXLC', '      QWFLXWC'         &
                               , '       QWSHED', '      QTRANSP', ' QINTERCEPTED'
                               

         close (unit=84,status='keep')
      end if
      !------------------------------------------------------------------------------------!
      



      !------------------------------------------------------------------------------------!
      !     Re-open the file at the last line, and include the current status.             !
      !------------------------------------------------------------------------------------!
      open (unit=84,file=trim(detail_fout),status='old',action='write',position='append')
      write(unit=84,fmt=cbfmt)                                                             &
              current_time%year       , current_time%month      , current_time%date        &
            , elapsec                 , hdid                    , cpatch%pft(ico)          &
            , leaf_resolve            , wood_resolve            , initp%nplant(ico)        &
            , cpatch%hite(ico)        , initp%lai(ico)          , initp%wai(ico)           &
            , initp%crown_area(ico)   , initp%leaf_energy(ico)  , initp%leaf_water(ico)    &
            , initp%leaf_hcap(ico)    , initp%leaf_temp(ico)    , initp%leaf_fliq(ico)     &
            , initp%wood_energy(ico)  , initp%wood_water(ico)   , initp%wood_hcap(ico)     &
            , initp%wood_temp(ico)    , initp%wood_fliq(ico)    , initp%veg_wind(ico)      &
            , initp%fs_open(ico)      , initp%leaf_reynolds(ico), initp%leaf_grashof(ico)  &
            , initp%leaf_nussfree(ico), initp%leaf_nussforc(ico), initp%wood_reynolds(ico) &
            , initp%wood_grashof(ico) , initp%wood_nussfree(ico), initp%wood_nussforc(ico) &
            , initp%lint_shv(ico)     , initp%leaf_gbh(ico)     , initp%leaf_gbw(ico)      &
            , initp%wood_gbh(ico)     , initp%wood_gbw(ico)     , initp%gsw_open(ico)      &
            , initp%gsw_closed(ico)   , initp%gpp(ico)          , initp%leaf_resp(ico)     &
            , initp%root_resp(ico)    , initp%growth_resp(ico)  , initp%storage_resp(ico)  &
            , initp%vleaf_resp(ico)   , initp%rshort_l(ico)     , initp%rlong_l(ico)       &
            , initp%rshort_w(ico)     , initp%rlong_w(ico)      , fluxp%cfx_hflxlc(ico)    &
            , fluxp%cfx_hflxwc(ico)   , fluxp%cfx_qwflxlc(ico)  , fluxp%cfx_qwflxwc(ico)   &
            , fluxp%cfx_qwshed(ico)   , fluxp%cfx_qtransp(ico)  , qintercepted
      close(unit=84,status='keep')
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!
   return
end subroutine print_rk4_state
!==========================================================================================!
!==========================================================================================!
