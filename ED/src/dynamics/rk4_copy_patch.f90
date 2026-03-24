!==========================================================================================!
!==========================================================================================!
! MODULE: RK4_COPY_PATCH
!
!> \brief   Routines that copy contents from, to or between rk4 structures 
!> \details These subroutines copy the content either from csite to rk4patch or between two
!!          rk4patch objects, or from rk4patch objects back to csite.  Sub-routine 
!!          copy_rk4_patch had been moved to a separate module to avoid circularities.  
!!          The module was expanded to include the other routines that copy variables, to
!!          organise the code.
!------------------------------------------------------------------------------------------!
module rk4_copy_patch
   contains

   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine copies that variables that are integrated by the Runge-Kutta       !
   ! solver to a buffer structure.                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine copy_rk4patch_init(sourcesite,ipa,ibuff,targetp,vels,old_can_enthalpy        &
                                ,old_can_rhos,old_can_dmol,ebudget_prsseffect)
      use ed_state_vars         , only : sitetype               & ! structure
                                       , patchtype              ! ! structure
      use grid_coms             , only : nzg                    & ! intent(in)
                                       , nzs                    ! ! intent(in) 
      use ed_misc_coms          , only : fast_diagnostics       ! ! intent(in)
      use consts_coms           , only : ep8                    & ! intent(in)
                                       , epim18                 & ! intent(in)
                                       , mmdryi8                & ! intent(in)
                                       , rmol8                  & ! intent(in)
                                       , rdry8                  & ! intent(in)
                                       , rdryi8                 & ! intent(in)
                                       , cpdry8                 & ! intent(in)
                                       , cph2o8                 & ! intent(in)
                                       , kgCday_2_umols8        ! ! intent(in)
      use rk4_coms              , only : rk4patchtype           & ! structure
                                       , rk4site                & ! structure
                                       , rk4aux                 & ! structure
                                       , tiny_offset            & ! intent(in)
                                       , rk4water_stab_thresh   & ! intent(in)
                                       , checkbudget            & ! intent(in)
                                       , print_detailed         & ! intent(in)
                                       , reset_rk4_fluxes       ! ! sub-routine
      use rk4_misc              , only : find_derived_thbounds  ! ! sub-routine
      use ed_max_dims           , only : n_pft                  ! ! intent(in)
      use canopy_air_coms       , only : ubmin8                 ! ! intent(in)
      use therm_lib8            , only : uextcm2tl8             & ! subroutine
                                       , cmtl2uext8             & ! function
                                       , thetaeiv8              & ! function
                                       , idealdenssh8           & ! function
                                       , idealdmolsh8           & ! function
                                       , rehuil8                & ! function
                                       , qslif8                 & ! function
                                       , reducedpress8          & ! function
                                       , press2exner8           & ! function
                                       , extheta2temp8          & ! function
                                       , tq2enthalpy8           ! ! function
      use ed_therm_lib          , only : ed_grndvap8            ! ! subroutine
      use canopy_struct_dynamics, only : canopy_turbulence8     & ! subroutine
                                       , can_whccap8            ! ! subroutine
      use budget_utils          , only : ddens_dt_effect8       & ! function
                                       , find_prss_effect8      ! ! function
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype)    , target      :: targetp
      type(sitetype)        , target      :: sourcesite
      integer               , intent(in)  :: ipa
      integer               , intent(in)  :: ibuff
      real                  , intent(in)  :: vels
      real                  , intent(in)  :: old_can_enthalpy
      real                  , intent(in)  :: old_can_rhos
      real                  , intent(in)  :: old_can_dmol
      real                  , intent(out) :: ebudget_prsseffect
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)       , pointer     :: cpatch
      real(kind=8)                        :: atm_tmp_zcan
      integer                             :: ico
      integer                             :: ipft
      integer                             :: k
      integer                             :: ksn
      real(kind=8)                        :: old_can_enthalpy8
      real(kind=8)                        :: old_can_rhos8
      real(kind=8)                        :: old_can_dmol8
      real(kind=8)                        :: ebudget_prsseffect8
      !----- External function. -----------------------------------------------------------!
      real(kind=4)          , external    :: sngloff
      !------------------------------------------------------------------------------------!

      !---- Alias for the current patch. --------------------------------------------------!
      cpatch => sourcesite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !---- Convert input variables to double precision. ----------------------------------!
      old_can_enthalpy8 = dble(old_can_enthalpy)
      old_can_rhos8     = dble(old_can_rhos    )
      old_can_dmol8     = dble(old_can_dmol    )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Between time steps the pressure may change because of change in atmospheric    !
      ! pressure, which means that temperature is not conserved.  Potential temperature    !
      ! and equivalent potential temperature, on the other hand, are conserved because     !
      ! there is no heat flux between time steps.  So we use these instead to start all    !
      ! other variables.                                                                   !
      !------------------------------------------------------------------------------------!
      !----- Update thermo variables that are conserved between steps. --------------------!
      targetp%can_theta    = dble(sourcesite%can_theta(ipa))
      targetp%can_shv      = dble(sourcesite%can_shv  (ipa))
      targetp%can_co2      = dble(sourcesite%can_co2  (ipa))
      targetp%can_depth    = dble(sourcesite%can_depth(ipa))
      !------------------------------------------------------------------------------------!

      !----- Update the vegetation properties used for roughness. -------------------------!
      targetp%veg_height   = dble(sourcesite%veg_height  (ipa))
      targetp%veg_displace = dble(sourcesite%veg_displace(ipa))
      targetp%veg_rough    = dble(sourcesite%veg_rough   (ipa))
      targetp%ustar        = dble(sourcesite%ustar       (ipa))
      targetp%tstar        = dble(sourcesite%tstar       (ipa))
      targetp%qstar        = dble(sourcesite%qstar       (ipa))
      targetp%cstar        = dble(sourcesite%cstar       (ipa))
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Update the canopy pressure and Exner function.                                !
      !------------------------------------------------------------------------------------!
      targetp%can_prss  = reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv &
                                       ,rk4site%geoht,targetp%can_theta,targetp%can_shv    &
                                       ,targetp%can_depth)
      targetp%can_exner = press2exner8 (targetp%can_prss)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the pressure and Exner functions at the canopy depth, find the temper-    !
      ! ature of the air above canopy at the canopy depth, and the specific enthalpy at    !
      ! that level.                                                                        !
      !------------------------------------------------------------------------------------!
      atm_tmp_zcan         = extheta2temp8(targetp%can_exner,rk4site%atm_theta)
      targetp%atm_enthalpy = tq2enthalpy8 (atm_tmp_zcan,rk4site%atm_shv,.true.)

      !----- Get velocity for aerodynamic resistance. -------------------------------------!
      targetp%vels  = max(ubmin8,dble(vels))

      !------------------------------------------------------------------------------------!
      !      Initialise canopy air temperature and enthalpy.  Enthalpy is the actual       !
      ! prognostic variable within one time step.                                          !
      !------------------------------------------------------------------------------------!
      targetp%can_temp     = extheta2temp8(targetp%can_exner,targetp%can_theta)
      targetp%can_enthalpy = tq2enthalpy8(targetp%can_temp,targetp%can_shv,.true.)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Set up the canopy air space heat capacity at constant pressure.               !
      !------------------------------------------------------------------------------------!
      targetp%can_cp = (1.d0 - targetp%can_shv) * cpdry8 + targetp%can_shv * cph2o8
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update density, dry-air molar density, relative humidity, and the saturation   !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      targetp%can_rhos = idealdenssh8(targetp%can_prss,targetp%can_temp,targetp%can_shv)
      targetp%can_dmol = idealdmolsh8(targetp%can_prss,targetp%can_temp,targetp%can_shv)
      targetp%can_rhv  = rehuil8(targetp%can_prss,targetp%can_temp,targetp%can_shv,.true.)
      targetp%can_ssh  = qslif8(targetp%can_prss,targetp%can_temp)
      !------------------------------------------------------------------------------------!



      !----- Find the lower and upper bounds for the derived properties. ------------------!
      call find_derived_thbounds(ibuff,cpatch,targetp%can_theta,targetp%can_temp           &
                                ,targetp%can_shv,targetp%can_prss,targetp%can_depth)
      !------------------------------------------------------------------------------------!


      !----- Impose a non-sense number for flag_wflxgc. -----------------------------------!
      targetp%flag_wflxgc  = -1
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Soil properties.                                                               !
      !------------------------------------------------------------------------------------!
      do k = rk4site%lsl, nzg
         !---------------------------------------------------------------------------------!
         !     Soil water may be slightly off-bounds.  This may happen when the soil       !
         ! moisture is exactly at the bounds, but then it is copied to single precision    !
         ! and then back to double precision.  Therefore at this time only we must ensure  !
         ! that we bound it, otherwise the model will crash due to the round-off error.    !
         !---------------------------------------------------------------------------------!
         targetp%soil_water  (k) = min( rk4aux(ibuff)%rk4max_soil_water(k)                 &
                                      , max( rk4aux(ibuff)%rk4min_soil_water(k)            &
                                      , dble(sourcesite%soil_water(k,ipa)) ) )
         targetp%soil_energy (k) = dble(sourcesite%soil_energy (k,ipa))
         targetp%soil_mstpot (k) = dble(sourcesite%soil_mstpot (k,ipa))
         targetp%soil_tempk  (k) = dble(sourcesite%soil_tempk  (k,ipa))
         targetp%soil_fracliq(k) = dble(sourcesite%soil_fracliq(k,ipa))
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy the surface water information.  The only non-trivial one is the energy,   !
      ! which is saved as J/kg outside the integration, but must be converted to J/m2      !
      ! because this linearises the differential equations and make the solution more      !
      ! stable.                                                                            !
      !------------------------------------------------------------------------------------!
      targetp%nlev_sfcwater    = sourcesite%nlev_sfcwater(ipa)
      ksn                      = targetp%nlev_sfcwater
      targetp%total_sfcw_mass  = 0.d0
      do k = 1, nzs
         targetp%sfcwater_mass(k)    = max(0.d0,dble(sourcesite%sfcwater_mass(k,ipa)))
         targetp%sfcwater_depth(k)   = dble(sourcesite%sfcwater_depth(k,ipa))
         targetp%sfcwater_energy(k)  = dble(sourcesite%sfcwater_energy(k,ipa))             &
                                     * dble(sourcesite%sfcwater_mass(k,ipa))
         targetp%sfcwater_tempk(k)   = dble(sourcesite%sfcwater_tempk(k,ipa))
         targetp%sfcwater_fracliq(k) = dble(sourcesite%sfcwater_fracliq(k,ipa))
         targetp%total_sfcw_mass     = targetp%total_sfcw_mass  + targetp%sfcwater_mass (k)
      end do
      !----- Define the temporary surface water flag. -------------------------------------!
      if (targetp%nlev_sfcwater == 0) then
         !----- No layer. -----------------------------------------------------------------!
         targetp%flag_sfcwater = 0
         !---------------------------------------------------------------------------------!
      elseif (targetp%total_sfcw_mass < rk4water_stab_thresh) then
         !----- There is water, but the amount is very small. -----------------------------!
         targetp%flag_sfcwater = 1
         !---------------------------------------------------------------------------------!
      else
         !----- There is enough water. ----------------------------------------------------!
         targetp%flag_sfcwater = 2
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the total snow depth and fraction of canopy buried in snow.               !
      !------------------------------------------------------------------------------------!
      targetp%snowfac          = dble(sourcesite%snowfac(ipa))
      targetp%total_sfcw_depth = dble(sourcesite%total_sfcw_depth(ipa))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute the ground temperature and specific humidity.                          !
      !------------------------------------------------------------------------------------!
      k = max(1,ksn)
      call ed_grndvap8(ksn,targetp%soil_water(nzg),targetp%soil_tempk(nzg)                 &
                      ,targetp%soil_fracliq(nzg),targetp%sfcwater_tempk(k)                 &
                      ,targetp%snowfac,targetp%can_prss                                    &
                      ,targetp%can_shv,targetp%ground_shv,targetp%ground_ssh               &
                      ,targetp%ground_temp,targetp%ground_fliq,targetp%ggsoil)
      !------------------------------------------------------------------------------------!



      !----- Initialise some turbulence properties. ---------------------------------------!
      targetp%upwp          = dble(sourcesite%upwp  (ipa))
      targetp%wpwp          = dble(sourcesite%wpwp  (ipa))
      targetp%tpwp          = dble(sourcesite%tpwp  (ipa))
      targetp%qpwp          = dble(sourcesite%qpwp  (ipa))
      targetp%cpwp          = dble(sourcesite%cpwp  (ipa))
      !------------------------------------------------------------------------------------!
     

      !------------------------------------------------------------------------------------!
      !     The virtual pools should be always zero, they are temporary entities used to   !
      ! store the water that falls as throughfall or shedding in the middle of a time      !
      ! step.  The temperature is liquid fraction is initialised as the soil temperature,  !
      ! or the temporary surface water temperature of the top most layer, if it exists.    !
      !------------------------------------------------------------------------------------!
      targetp%virtual_water  = 0.0d0
      targetp%virtual_energy = 0.0d0
      targetp%virtual_depth  = 0.0d0
      if (ksn == 0) then
         targetp%virtual_tempk   = targetp%soil_tempk  (nzg)
         targetp%virtual_fracliq = targetp%soil_fracliq(nzg)
      else
         targetp%virtual_tempk   = targetp%sfcwater_tempk  (ksn)
         targetp%virtual_fracliq = targetp%sfcwater_fracliq(ksn)
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Here we find the minimum patch-level leaf heat capacity.  If the total patch   !
      ! leaf heat capacity is less than this, we scale the cohorts heat capacity inside    !
      ! the integrator, so it preserves the proportional heat capacity and prevents the    !
      ! pool to be too small.                                                              !
      !------------------------------------------------------------------------------------!
      rk4aux(ibuff)%any_resolvable = .false.
      do ico=1, cpatch%ncohorts
         !----- Copy the flag that determines whether this cohort is numerically stable. --!
         targetp%leaf_resolvable(ico) = cpatch%leaf_resolvable(ico)
         targetp%wood_resolvable(ico) = cpatch%wood_resolvable(ico)
         if (targetp%leaf_resolvable(ico) .or. targetp%wood_resolvable(ico)) then
            rk4aux(ibuff)%any_resolvable = .true.
         end if
         !---------------------------------------------------------------------------------!


         !----- Small/large tree flag. ----------------------------------------------------!
         targetp%is_small(ico) = cpatch%is_small(ico)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over cohorts.  We also use this loop to compute the fraction of open      !
      ! canopy, which is initialised with 1. in case this is an empty patch.               !
      !------------------------------------------------------------------------------------!
      targetp%opencan_frac   = dble(sourcesite%opencan_frac(ipa))
      do ico = 1,cpatch%ncohorts
         ipft=cpatch%pft(ico)

         !----- Copy the plant density. ---------------------------------------------------!
         targetp%nplant(ico)     = dble(cpatch%nplant(ico))

         !----- Copy the leaf area index and total (leaf+branch+twig) area index. ---------!
         targetp%lai(ico)        = dble(cpatch%lai(ico))
         targetp%wai(ico)        = dble(cpatch%wai(ico))
         targetp%tai(ico)        = targetp%lai(ico) + targetp%wai(ico)
         targetp%crown_area(ico) = dble(cpatch%crown_area(ico))
         targetp%elongf(ico)     = dble(cpatch%elongf(ico))                                &
                                 * rk4site%green_leaf_factor(ipft)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Check whether the leaves of this cohort are considered "safe" or not.  In   !
         ! case they are, we copy the water, heat capacity, and energy then compute the    !
         ! temperature and the fraction of leaf water.  Otherwise, just fill with some     !
         ! safe values, but the leaves won't be really solved.                             !
         !---------------------------------------------------------------------------------!
         if (targetp%leaf_resolvable(ico)) then
            targetp%leaf_energy   (ico) = dble(cpatch%leaf_energy(ico))
            targetp%leaf_water    (ico) = max(0.d0,dble(cpatch%leaf_water    (ico)))
            targetp%leaf_water_im2(ico) = max(0.d0,dble(cpatch%leaf_water_im2(ico)))
            targetp%leaf_hcap     (ico) = dble(cpatch%leaf_hcap  (ico))

            call uextcm2tl8(targetp%leaf_energy(ico)                                       &
                           ,targetp%leaf_water (ico) + targetp%leaf_water_im2(ico)         &
                           ,targetp%leaf_hcap(ico),targetp%leaf_temp(ico)                  &
                           ,targetp%leaf_fliq(ico))
         else
            targetp%leaf_fliq     (ico) = dble(cpatch%leaf_fliq  (ico))
            targetp%leaf_temp     (ico) = dble(cpatch%leaf_temp  (ico))
            targetp%leaf_water    (ico) = max(0.d0,dble(cpatch%leaf_water    (ico)))
            targetp%leaf_water_im2(ico) = max(0.d0,dble(cpatch%leaf_water_im2(ico)))
            targetp%leaf_hcap     (ico) = dble(cpatch%leaf_hcap  (ico))
            targetp%leaf_energy   (ico) = cmtl2uext8( targetp%leaf_hcap     (ico)          &
                                                    , targetp%leaf_water    (ico)          &
                                                    + targetp%leaf_water_im2(ico)          &
                                                    , targetp%leaf_temp     (ico)          &
                                                    , targetp%leaf_fliq     (ico) )
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Check whether the wood of this cohort is considered "safe" or not.  In case !
         ! it is, we copy the water, heat capacity, and energy then compute the temper-    !
         ! ature and the fraction of leaf water.  Otherwise, just fill with some safe      !
         ! values, but the wood won't be really solved.                                    !
         !---------------------------------------------------------------------------------!
         if (targetp%wood_resolvable(ico)) then
            targetp%wood_energy   (ico) = dble(cpatch%wood_energy(ico))
            targetp%wood_water    (ico) = max(0.d0,dble(cpatch%wood_water    (ico)))
            targetp%wood_water_im2(ico) = max(0.d0,dble(cpatch%wood_water_im2(ico)))
            targetp%wood_hcap     (ico) = dble(cpatch%wood_hcap  (ico))

            call uextcm2tl8(targetp%wood_energy(ico)                                       &
                           ,targetp%wood_water (ico) + targetp%wood_water_im2(ico)         &
                           ,targetp%wood_hcap(ico),targetp%wood_temp(ico)                  &
                           ,targetp%wood_fliq(ico))
         else
            targetp%wood_fliq     (ico) = dble(cpatch%wood_fliq  (ico))
            targetp%wood_temp     (ico) = dble(cpatch%wood_temp  (ico))
            targetp%wood_water    (ico) = max(0.d0,dble(cpatch%wood_water    (ico)))
            targetp%wood_water_im2(ico) = max(0.d0,dble(cpatch%wood_water_im2(ico)))
            targetp%wood_hcap     (ico) = dble(cpatch%wood_hcap  (ico))
            targetp%wood_energy   (ico) = cmtl2uext8( targetp%wood_hcap     (ico)          &
                                                    , targetp%wood_water    (ico)          &
                                                    + targetp%wood_water_im2(ico)          &
                                                    , targetp%wood_temp     (ico)          &
                                                    , targetp%wood_fliq     (ico) )
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   Make the combined leaf and branchwood pool.  It will be really used only if   !
         ! ibranch_thermo is set to 1.                                                     !
         !---------------------------------------------------------------------------------!
         targetp%veg_resolvable(ico) = targetp%leaf_resolvable(ico) .or.                   &
                                       targetp%wood_resolvable(ico)
         targetp%veg_energy(ico)     = targetp%leaf_energy(ico) + targetp%wood_energy(ico)
         targetp%veg_water(ico)      = targetp%leaf_water(ico)  + targetp%wood_water(ico)
         targetp%veg_water_im2(ico)  = targetp%leaf_water_im2(ico)                         &
                                     + targetp%wood_water_im2(ico)
         targetp%veg_hcap(ico)       = targetp%leaf_hcap(ico)   + targetp%wood_hcap(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Compute the leaf intercellular specific humidity, assumed to be at          !
         ! saturation.                                                                     !
         !---------------------------------------------------------------------------------!
         targetp%lint_shv(ico) = qslif8(targetp%can_prss,targetp%leaf_temp(ico))
         !---------------------------------------------------------------------------------!



         !------ Copy the stomatal conductances and the fraction of open stomata. ---------!
         targetp%fs_open   (ico) = dble(cpatch%fs_open   (ico))
         targetp%gsw_open  (ico) = dble(cpatch%gsw_open  (ico))
         targetp%gsw_closed(ico) = dble(cpatch%gsw_closed(ico))
         !---------------------------------------------------------------------------------!



         !------ Copy the net absorbed radiation (short wave and long wave). --------------!
         targetp%rshort_l   (ico) = dble(cpatch%rshort_l (ico))
         targetp%rlong_l    (ico) = dble(cpatch%rlong_l  (ico))
         targetp%rshort_w   (ico) = dble(cpatch%rshort_w (ico))
         targetp%rlong_w    (ico) = dble(cpatch%rlong_w  (ico))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Initialise psi_open and psi_closed with zeroes, they will be averaged over  !
         ! the course of one Runge-Kutta time step.                                        !
         !---------------------------------------------------------------------------------!
         targetp%psi_open(ico)   = 0.d0
         targetp%psi_closed(ico) = 0.d0
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Here we copy the cohort level variables that are part of the carbon budget.    !
      !------------------------------------------------------------------------------------!
      cpatch => sourcesite%patch(ipa)
      do ico = 1,cpatch%ncohorts
     
         !----- Copy the variables that are already in umol/m2/s. -------------------------!
         targetp%gpp         (ico) = dble(cpatch%gpp                (ico))
         targetp%leaf_resp   (ico) = dble(cpatch%leaf_respiration   (ico))
         targetp%root_resp   (ico) = dble(cpatch%root_respiration   (ico))
         targetp%stem_resp   (ico) = dble(cpatch%stem_respiration   (ico))
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following variables are in kgC/m2/day, convert them to umol/m2/s.          !
      !------------------------------------------------------------------------------------!
      targetp%commit_storage_resp  = dble(sourcesite%commit_storage_resp (ipa))            &
                                   * kgCday_2_umols8
      targetp%commit_growth_resp   = dble(sourcesite%commit_growth_resp  (ipa))            &
                                   * kgCday_2_umols8
      !------------------------------------------------------------------------------------!

      !----- Heterotrophic respiration terms, already in umol/m2/s. -----------------------!
      targetp%fgc_rh  = dble(sourcesite%fgc_rh (ipa))
      targetp%fsc_rh  = dble(sourcesite%fsc_rh (ipa))
      targetp%stgc_rh = dble(sourcesite%stgc_rh(ipa))
      targetp%stsc_rh = dble(sourcesite%stsc_rh(ipa))
      targetp%msc_rh  = dble(sourcesite%msc_rh (ipa))
      targetp%ssc_rh  = dble(sourcesite%ssc_rh (ipa))
      targetp%psc_rh  = dble(sourcesite%psc_rh (ipa))
      !------------------------------------------------------------------------------------!



      !----- Initialise the characteristic properties, and the heat capacities. -----------!
      call canopy_turbulence8(sourcesite,targetp,ipa,ibuff)
      call can_whccap8(targetp%can_rhos,targetp%can_dmol,targetp%can_depth                 &
                      ,targetp%wcapcan ,targetp%hcapcan ,targetp%ccapcan                   &
                      ,targetp%wcapcani,targetp%hcapcani,targetp%ccapcani)
      !------------------------------------------------------------------------------------!

      !----- Diagnostics variables --------------------------------------------------------!
      if(fast_diagnostics) then
         !---------------------------------------------------------------------------------!
         !     The "budget" variables are not copied here because they are integrated out- !
         ! side RK4.  Inside RK4 we only want the contribution of those variables during   !
         ! the span  of one time step.                                                     !
         !---------------------------------------------------------------------------------!
         targetp%avg_ustar          = dble(sourcesite%fmean_ustar         (ipa))
         targetp%avg_tstar          = dble(sourcesite%fmean_tstar         (ipa))
         targetp%avg_qstar          = dble(sourcesite%fmean_qstar         (ipa))
         targetp%avg_cstar          = dble(sourcesite%fmean_cstar         (ipa))
         targetp%avg_carbon_ac      = dble(sourcesite%fmean_carbon_ac     (ipa))
         targetp%avg_carbon_st      = dble(sourcesite%fmean_carbon_st     (ipa))
         targetp%avg_vapor_gc       = dble(sourcesite%fmean_vapor_gc      (ipa))
         targetp%avg_throughfall    = dble(sourcesite%fmean_throughfall   (ipa))
         targetp%avg_vapor_ac       = dble(sourcesite%fmean_vapor_ac      (ipa))
         targetp%avg_drainage       = dble(sourcesite%fmean_drainage      (ipa))
         targetp%avg_qdrainage      = dble(sourcesite%fmean_qdrainage     (ipa))
         targetp%avg_qthroughfall   = dble(sourcesite%fmean_qthroughfall  (ipa))
         targetp%avg_sensible_gc    = dble(sourcesite%fmean_sensible_gc   (ipa))
         targetp%avg_sensible_ac    = dble(sourcesite%fmean_sensible_ac   (ipa))

         do k = rk4site%lsl, nzg
            targetp%avg_sensible_gg(k) = dble(sourcesite%fmean_sensible_gg(k,ipa))
            targetp%avg_smoist_gg(k)   = dble(sourcesite%fmean_smoist_gg(k,ipa)  )
            targetp%avg_transloss(k)   = dble(sourcesite%fmean_transloss(k,ipa)  )
         end do

         do ico=1,cpatch%ncohorts
            targetp%avg_sensible_lc   (ico) = dble(cpatch%fmean_sensible_lc   (ico))
            targetp%avg_sensible_wc   (ico) = dble(cpatch%fmean_sensible_wc   (ico))
            targetp%avg_vapor_lc      (ico) = dble(cpatch%fmean_vapor_lc      (ico))
            targetp%avg_vapor_wc      (ico) = dble(cpatch%fmean_vapor_wc      (ico))
            targetp%avg_transp        (ico) = dble(cpatch%fmean_transp        (ico))
            targetp%avg_intercepted_al(ico) = dble(cpatch%fmean_intercepted_al(ico))
            targetp%avg_intercepted_aw(ico) = dble(cpatch%fmean_intercepted_aw(ico))
            targetp%avg_wshed_lg      (ico) = dble(cpatch%fmean_wshed_lg      (ico))
            targetp%avg_wshed_wg      (ico) = dble(cpatch%fmean_wshed_wg      (ico))
            targetp%avg_wflux_wl      (ico) = dble(cpatch%fmean_wflux_wl      (ico))       &
                                            * targetp%nplant(ico)
            targetp%avg_wflux_gw      (ico) = dble(cpatch%fmean_wflux_gw      (ico))       &
                                            * targetp%nplant(ico)

            do k = rk4site%lsl, nzg
               targetp%avg_wflux_gw_layer(k,ico) =                                         &
                             dble(cpatch%fmean_wflux_gw_layer(k,ico)) * targetp%nplant(ico)
            end do

         end do

      end if
      if (checkbudget) then
         !----- Initial storage from the previous time step. ------------------------------!
         targetp%co2budget_storage     = dble(sourcesite%co2budget_initialstorage(ipa))
         targetp%ebudget_storage       = dble(sourcesite%ebudget_initialstorage  (ipa))
         targetp%wbudget_storage       = dble(sourcesite%wbudget_initialstorage  (ipa))
         !---------------------------------------------------------------------------------!


         !----- These variables will be integrated over the course of a time step. --------!
         targetp%co2budget_loss2atm    = 0.d0
         targetp%ebudget_netrad        = 0.d0
         targetp%ebudget_loss2atm      = 0.d0
         targetp%ebudget_loss2drainage = 0.d0
         targetp%ebudget_loss2runoff   = 0.d0
         targetp%wbudget_loss2atm      = 0.d0
         targetp%wbudget_loss2drainage = 0.d0
         targetp%wbudget_loss2runoff   = 0.d0
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Calculate the pressure effect -- changes in bulk enthalpy associated with   !
         ! changes in pressure.                                                            !
         !---------------------------------------------------------------------------------!
         ebudget_prsseffect8  = find_prss_effect8( targetp%can_depth                       &
                                                 , old_can_rhos8    , targetp%can_rhos     &
                                                 , old_can_enthalpy8, targetp%can_enthalpy )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Density also changed when pressure updates caused enthalpy to change.  In   !
         ! the case of enthalpy, this is already incorporated in the pressure effect, but  !
         ! this effect also changes total CO2 and total water in the CAS.  We incorporate  !
         ! these changes in the initial density effect.                                    !
         !---------------------------------------------------------------------------------!
         targetp%co2budget_denseffect = ddens_dt_effect8 ( targetp%can_depth               &
                                                         , old_can_dmol8                   &
                                                         , targetp%can_dmol                &
                                                         , targetp%can_co2      )
         targetp%ebudget_denseffect   = 0.d0
         targetp%wbudget_denseffect   = ddens_dt_effect8 ( targetp%can_depth               &
                                                         , old_can_rhos8                   &
                                                         , targetp%can_rhos                &
                                                         , targetp%can_shv      )
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !----- Water deficit, always start with zero. ---------------------------------------!
      targetp%water_deficit = 0.d0
      !------------------------------------------------------------------------------------!

      !----- If writing the detailed output, reset fluxes. --------------------------------!
      if (print_detailed) call reset_rk4_fluxes(targetp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Save the pressure effect, this is the only time it happens so it does not     !
      ! need to be copied to the RK4 structures.                                           !
      !------------------------------------------------------------------------------------!
      ebudget_prsseffect = sngloff(ebudget_prsseffect8,tiny_offset)
      !------------------------------------------------------------------------------------!
      return
   end subroutine copy_rk4patch_init
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine copies the values to different buffers inside the RK4 integration  !
   ! scheme.                                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine copy_rk4_patch(sourcep, targetp, cpatch)

      use rk4_coms      , only : rk4patchtype      & ! structure
                               , checkbudget       & ! intent(in)
                               , print_detailed    & ! intent(in)
                               , rk4site
      use ed_state_vars , only : sitetype          & ! structure
                               , patchtype         ! ! structure
      use grid_coms     , only : nzg               & ! intent(in)
                               , nzs               ! ! intent(in)
      use ed_max_dims   , only : n_pft             ! ! intent(in)
      use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: sourcep
      type(rk4patchtype) , target     :: targetp
      type(patchtype)    , target     :: cpatch
      !----- Local variable ---------------------------------------------------------------!
      integer                         :: ico
      integer                         :: k
      !------------------------------------------------------------------------------------!

      targetp%can_enthalpy     = sourcep%can_enthalpy
      targetp%can_theta        = sourcep%can_theta
      targetp%can_temp         = sourcep%can_temp
      targetp%can_shv          = sourcep%can_shv
      targetp%can_co2          = sourcep%can_co2
      targetp%can_rhos         = sourcep%can_rhos
      targetp%can_dmol         = sourcep%can_dmol
      targetp%can_prss         = sourcep%can_prss
      targetp%can_exner        = sourcep%can_exner
      targetp%can_cp           = sourcep%can_cp
      targetp%can_depth        = sourcep%can_depth
      targetp%can_rhv          = sourcep%can_rhv
      targetp%can_ssh          = sourcep%can_ssh
      targetp%veg_height       = sourcep%veg_height
      targetp%veg_displace     = sourcep%veg_displace
      targetp%veg_rough        = sourcep%veg_rough
      targetp%opencan_frac     = sourcep%opencan_frac
      targetp%total_sfcw_depth = sourcep%total_sfcw_depth
      targetp%total_sfcw_mass  = sourcep%total_sfcw_mass
      targetp%snowfac          = sourcep%snowfac

   !  These are not incremented
      targetp%vels             = sourcep%vels
      targetp%atm_enthalpy     = sourcep%atm_enthalpy

      targetp%ggbare           = sourcep%ggbare
      targetp%ggveg            = sourcep%ggveg
      targetp%ggnet            = sourcep%ggnet
      targetp%ggsoil           = sourcep%ggsoil

      targetp%flag_wflxgc      = sourcep%flag_wflxgc

      targetp%virtual_water    = sourcep%virtual_water
      targetp%virtual_energy   = sourcep%virtual_energy
      targetp%virtual_depth    = sourcep%virtual_depth
      targetp%virtual_tempk    = sourcep%virtual_tempk
      targetp%virtual_fracliq  = sourcep%virtual_fracliq

      targetp%rough            = sourcep%rough
    
      targetp%upwp             = sourcep%upwp
      targetp%wpwp             = sourcep%wpwp
      targetp%tpwp             = sourcep%tpwp
      targetp%qpwp             = sourcep%qpwp
      targetp%cpwp             = sourcep%cpwp

      targetp%ground_shv       = sourcep%ground_shv
      targetp%ground_ssh       = sourcep%ground_ssh
      targetp%ground_temp      = sourcep%ground_temp
      targetp%ground_fliq      = sourcep%ground_fliq

      targetp%wcapcan          = sourcep%wcapcan
      targetp%hcapcan          = sourcep%hcapcan
      targetp%ccapcan          = sourcep%ccapcan
      targetp%wcapcani         = sourcep%wcapcani
      targetp%hcapcani         = sourcep%hcapcani
      targetp%ccapcani         = sourcep%ccapcani

      targetp%ustar            = sourcep%ustar
      targetp%cstar            = sourcep%cstar
      targetp%tstar            = sourcep%tstar
      targetp%estar            = sourcep%estar
      targetp%qstar            = sourcep%qstar
      targetp%zeta             = sourcep%zeta
      targetp%ribulk           = sourcep%ribulk
      targetp%rasveg           = sourcep%rasveg

      targetp%fgc_rh           = sourcep%fgc_rh
      targetp%fsc_rh           = sourcep%fsc_rh
      targetp%stgc_rh          = sourcep%stgc_rh
      targetp%stsc_rh          = sourcep%stsc_rh
      targetp%msc_rh           = sourcep%msc_rh
      targetp%ssc_rh           = sourcep%ssc_rh
      targetp%psc_rh           = sourcep%psc_rh

      targetp%water_deficit    = sourcep%water_deficit

      targetp%commit_storage_resp = sourcep%commit_storage_resp
      targetp%commit_growth_resp  = sourcep%commit_growth_resp

      do k=rk4site%lsl,nzg      
         targetp%soil_water            (k) = sourcep%soil_water            (k)
         targetp%soil_energy           (k) = sourcep%soil_energy           (k)
         targetp%soil_mstpot           (k) = sourcep%soil_mstpot           (k)
         targetp%soil_tempk            (k) = sourcep%soil_tempk            (k)
         targetp%soil_fracliq          (k) = sourcep%soil_fracliq          (k)
      end do

      targetp%nlev_sfcwater   = sourcep%nlev_sfcwater
      targetp%flag_sfcwater   = sourcep%flag_sfcwater

      do k=1,nzs
         targetp%sfcwater_mass     (k) = sourcep%sfcwater_mass     (k)
         targetp%sfcwater_energy   (k) = sourcep%sfcwater_energy   (k)
         targetp%sfcwater_depth    (k) = sourcep%sfcwater_depth    (k)
         targetp%sfcwater_tempk    (k) = sourcep%sfcwater_tempk    (k)
         targetp%sfcwater_fracliq  (k) = sourcep%sfcwater_fracliq  (k)
      end do

      do ico=1,cpatch%ncohorts
         targetp%leaf_resolvable   (ico) = sourcep%leaf_resolvable   (ico)
         targetp%leaf_energy       (ico) = sourcep%leaf_energy       (ico)
         targetp%leaf_water        (ico) = sourcep%leaf_water        (ico)
         targetp%leaf_water_im2    (ico) = sourcep%leaf_water_im2    (ico)
         targetp%leaf_temp         (ico) = sourcep%leaf_temp         (ico)
         targetp%leaf_fliq         (ico) = sourcep%leaf_fliq         (ico)
         targetp%leaf_hcap         (ico) = sourcep%leaf_hcap         (ico)
         targetp%leaf_reynolds     (ico) = sourcep%leaf_reynolds     (ico)
         targetp%leaf_grashof      (ico) = sourcep%leaf_grashof      (ico)
         targetp%leaf_nussfree     (ico) = sourcep%leaf_nussfree     (ico)
         targetp%leaf_nussforc     (ico) = sourcep%leaf_nussforc     (ico)
         targetp%leaf_gbh          (ico) = sourcep%leaf_gbh          (ico)
         targetp%leaf_gbw          (ico) = sourcep%leaf_gbw          (ico)
         targetp%rshort_l          (ico) = sourcep%rshort_l          (ico)
         targetp%rlong_l           (ico) = sourcep%rlong_l           (ico)

         targetp%wood_resolvable   (ico) = sourcep%wood_resolvable   (ico)
         targetp%wood_energy       (ico) = sourcep%wood_energy       (ico)
         targetp%wood_water        (ico) = sourcep%wood_water        (ico)
         targetp%wood_water_im2    (ico) = sourcep%wood_water_im2    (ico)
         targetp%wood_temp         (ico) = sourcep%wood_temp         (ico)
         targetp%wood_fliq         (ico) = sourcep%wood_fliq         (ico)
         targetp%wood_hcap         (ico) = sourcep%wood_hcap         (ico)
         targetp%wood_reynolds     (ico) = sourcep%wood_reynolds     (ico)
         targetp%wood_grashof      (ico) = sourcep%wood_grashof      (ico)
         targetp%wood_nussfree     (ico) = sourcep%wood_nussfree     (ico)
         targetp%wood_nussforc     (ico) = sourcep%wood_nussforc     (ico)
         targetp%wood_gbh          (ico) = sourcep%wood_gbh          (ico)
         targetp%wood_gbw          (ico) = sourcep%wood_gbw          (ico)
         targetp%rshort_w          (ico) = sourcep%rshort_w          (ico)
         targetp%rlong_w           (ico) = sourcep%rlong_w           (ico)

         targetp%veg_resolvable    (ico) = sourcep%veg_resolvable    (ico)
         targetp%veg_energy        (ico) = sourcep%veg_energy        (ico)
         targetp%veg_water         (ico) = sourcep%veg_water         (ico)
         targetp%veg_water_im2     (ico) = sourcep%veg_water_im2     (ico)
         targetp%veg_hcap          (ico) = sourcep%veg_hcap          (ico)

         targetp%is_small          (ico) = sourcep%is_small          (ico)
         targetp%veg_wind          (ico) = sourcep%veg_wind          (ico)
         targetp%lint_shv          (ico) = sourcep%lint_shv          (ico)
         targetp%nplant            (ico) = sourcep%nplant            (ico)
         targetp%lai               (ico) = sourcep%lai               (ico)
         targetp%wai               (ico) = sourcep%wai               (ico)
         targetp%tai               (ico) = sourcep%tai               (ico)
         targetp%crown_area        (ico) = sourcep%crown_area        (ico)
         targetp%elongf            (ico) = sourcep%elongf            (ico)
         targetp%gsw_open          (ico) = sourcep%gsw_open          (ico)
         targetp%gsw_closed        (ico) = sourcep%gsw_closed        (ico)
         targetp%psi_open          (ico) = sourcep%psi_open          (ico)
         targetp%psi_closed        (ico) = sourcep%psi_closed        (ico)
         targetp%fs_open           (ico) = sourcep%fs_open           (ico)
         targetp%gpp               (ico) = sourcep%gpp               (ico)
         targetp%leaf_resp         (ico) = sourcep%leaf_resp         (ico)
         targetp%root_resp         (ico) = sourcep%root_resp         (ico)
         targetp%stem_resp         (ico) = sourcep%stem_resp         (ico)
      end do

      if (checkbudget) then
         targetp%co2budget_storage      = sourcep%co2budget_storage
         targetp%co2budget_loss2atm     = sourcep%co2budget_loss2atm
         targetp%co2budget_denseffect   = sourcep%co2budget_denseffect
         targetp%ebudget_netrad         = sourcep%ebudget_netrad
         targetp%ebudget_loss2atm       = sourcep%ebudget_loss2atm
         targetp%ebudget_loss2drainage  = sourcep%ebudget_loss2drainage
         targetp%ebudget_loss2runoff    = sourcep%ebudget_loss2runoff
         targetp%ebudget_denseffect     = sourcep%ebudget_denseffect
         targetp%wbudget_loss2atm       = sourcep%wbudget_loss2atm
         targetp%wbudget_loss2drainage  = sourcep%wbudget_loss2drainage
         targetp%wbudget_loss2runoff    = sourcep%wbudget_loss2runoff
         targetp%ebudget_storage        = sourcep%ebudget_storage
         targetp%wbudget_storage        = sourcep%wbudget_storage
         targetp%wbudget_denseffect     = sourcep%wbudget_denseffect
      end if
      if (fast_diagnostics) then
         targetp%avg_ustar              = sourcep%avg_ustar
         targetp%avg_tstar              = sourcep%avg_tstar
         targetp%avg_qstar              = sourcep%avg_qstar
         targetp%avg_cstar              = sourcep%avg_cstar
         targetp%avg_carbon_ac          = sourcep%avg_carbon_ac
         targetp%avg_carbon_st          = sourcep%avg_carbon_st
         targetp%avg_vapor_gc           = sourcep%avg_vapor_gc
         targetp%avg_throughfall        = sourcep%avg_throughfall
         targetp%avg_vapor_ac           = sourcep%avg_vapor_ac
         targetp%avg_qthroughfall       = sourcep%avg_qthroughfall
         targetp%avg_sensible_gc        = sourcep%avg_sensible_gc
         targetp%avg_sensible_ac        = sourcep%avg_sensible_ac
         targetp%avg_drainage           = sourcep%avg_drainage
         targetp%avg_qdrainage          = sourcep%avg_qdrainage

         do k=rk4site%lsl,nzg
            targetp%avg_sensible_gg(k) = sourcep%avg_sensible_gg(k)
            targetp%avg_smoist_gg(k)   = sourcep%avg_smoist_gg(k)
            targetp%avg_transloss(k)   = sourcep%avg_transloss(k)
         end do


         do ico=1,cpatch%ncohorts
            targetp%avg_sensible_lc    (ico) = sourcep%avg_sensible_lc    (ico)
            targetp%avg_sensible_wc    (ico) = sourcep%avg_sensible_wc    (ico)
            targetp%avg_vapor_lc       (ico) = sourcep%avg_vapor_lc       (ico)
            targetp%avg_vapor_wc       (ico) = sourcep%avg_vapor_wc       (ico)
            targetp%avg_transp         (ico) = sourcep%avg_transp         (ico)
            targetp%avg_intercepted_al (ico) = sourcep%avg_intercepted_al (ico)
            targetp%avg_intercepted_aw (ico) = sourcep%avg_intercepted_aw (ico)
            targetp%avg_wshed_lg       (ico) = sourcep%avg_wshed_lg       (ico)
            targetp%avg_wshed_wg       (ico) = sourcep%avg_wshed_wg       (ico)
            targetp%avg_wflux_wl       (ico) = sourcep%avg_wflux_wl       (ico)
            targetp%avg_wflux_gw       (ico) = sourcep%avg_wflux_gw       (ico)
         end do
      end if

      if (print_detailed) then
         targetp%flx_carbon_ac          = sourcep%flx_carbon_ac
         targetp%flx_carbon_st          = sourcep%flx_carbon_st
         targetp%flx_vapor_lc           = sourcep%flx_vapor_lc
         targetp%flx_vapor_wc           = sourcep%flx_vapor_wc
         targetp%flx_vapor_gc           = sourcep%flx_vapor_gc
         targetp%flx_wshed_vg           = sourcep%flx_wshed_vg
         targetp%flx_wflux_wl           = sourcep%flx_wflux_wl
         targetp%flx_wflux_gw           = sourcep%flx_wflux_gw
         targetp%flx_intercepted        = sourcep%flx_intercepted
         targetp%flx_throughfall        = sourcep%flx_throughfall
         targetp%flx_vapor_ac           = sourcep%flx_vapor_ac
         targetp%flx_transp             = sourcep%flx_transp
         targetp%flx_rshort_gnd         = sourcep%flx_rshort_gnd
         targetp%flx_par_gnd            = sourcep%flx_par_gnd
         targetp%flx_rlong_gnd          = sourcep%flx_rlong_gnd
         targetp%flx_sensible_lc        = sourcep%flx_sensible_lc
         targetp%flx_sensible_wc        = sourcep%flx_sensible_wc
         targetp%flx_qwshed_vg          = sourcep%flx_qwshed_vg
         targetp%flx_qintercepted       = sourcep%flx_qintercepted
         targetp%flx_qthroughfall       = sourcep%flx_qthroughfall
         targetp%flx_sensible_gc        = sourcep%flx_sensible_gc
         targetp%flx_sensible_ac        = sourcep%flx_sensible_ac
         targetp%flx_drainage           = sourcep%flx_drainage
         targetp%flx_qdrainage          = sourcep%flx_qdrainage

         do k=rk4site%lsl,nzg
            targetp%flx_sensible_gg(k) = sourcep%flx_sensible_gg(k)
            targetp%flx_smoist_gg(k)   = sourcep%flx_smoist_gg(k)  
            targetp%flx_transloss(k)   = sourcep%flx_transloss(k)  
         end do

         do ico=1,cpatch%ncohorts
            targetp%cfx_hflxlc      (ico) = sourcep%cfx_hflxlc      (ico)
            targetp%cfx_hflxwc      (ico) = sourcep%cfx_hflxwc      (ico)
            targetp%cfx_qwflxlc     (ico) = sourcep%cfx_qwflxlc     (ico)
            targetp%cfx_qwflxwc     (ico) = sourcep%cfx_qwflxwc     (ico)
            targetp%cfx_qwshed      (ico) = sourcep%cfx_qwshed      (ico)
            targetp%cfx_qtransp     (ico) = sourcep%cfx_qtransp     (ico)
            targetp%cfx_qintercepted(ico) = sourcep%cfx_qintercepted(ico)
            targetp%cfx_qwflux_wl   (ico) = sourcep%cfx_qwflux_wl   (ico)
            targetp%cfx_qwflux_gw   (ico) = sourcep%cfx_qwflux_gw   (ico)
            do k=rk4site%lsl,nzg
               targetp%cfx_qwflux_gw_layer(k,ico) = sourcep%cfx_qwflux_gw_layer(k,ico)
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      return
   end subroutine copy_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp(hdid,initp,csite,ipa,nighttime,wbudget_loss2atm                 &
                          ,ebudget_netrad,ebudget_loss2atm,co2budget_loss2atm              &
                          ,wbudget_loss2drainage,ebudget_loss2drainage,wbudget_loss2runoff &
                          ,ebudget_loss2runoff,co2budget_denseffect,ebudget_denseffect     &
                          ,wbudget_denseffect)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4site              & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          & ! intent(in)
                                      , checkbudget          & ! intent(in)
                                      , ibranch_thermo       ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            ! ! structure
      use canopy_air_coms      , only : f_bndlyr_init        ! ! intent(in)
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , cpdry                & ! intent(in)
                                      , t3ple                & ! intent(in)
                                      , t3ple8               & ! intent(in)
                                      , wdns8                ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     & ! intent(in)
                                      , writing_long         & ! intent(in)
                                      , dtlsm                & ! intent(in)
                                      , dtlsm_o_frqsum       ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , dslz8                & ! intent(in)
                                      , slz8                 & ! intent(in)
                                      , slzt8                & ! intent(in)
                                      , matric_potential8    ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : thetaeiv             & ! subroutine
                                      , vpdefil              & ! subroutine
                                      , uextcm2tl            & ! subroutine
                                      , cmtl2uext            & ! subroutine
                                      , qslif                ! ! function
      use phenology_coms       , only : spot_phen            ! ! intent(in)
      use physiology_coms      , only : plant_hydro_scheme   & ! intent(in)
                                      , gbh_2_gbw            ! ! intent(in)
      use allometry            , only : h2crownbh            ! ! function
      use disturb_coms         , only : include_fire         & ! intent(in)
                                      , k_fire_first         ! ! intent(in)
      use plant_hydro          , only : twe2twi              & ! subroutine
                                      , tw2rwc               ! ! subroutine
      use rk4_misc             , only : print_rk4patch       ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target      :: initp
      type(sitetype)    , target      :: csite
      real(kind=8)      , intent(in)  :: hdid
      integer           , intent(in)  :: ipa
      logical           , intent(in)  :: nighttime
      real              , intent(out) :: wbudget_loss2atm
      real              , intent(out) :: ebudget_netrad
      real              , intent(out) :: ebudget_loss2atm
      real              , intent(out) :: co2budget_loss2atm
      real              , intent(out) :: wbudget_loss2drainage
      real              , intent(out) :: ebudget_loss2drainage
      real              , intent(out) :: wbudget_loss2runoff
      real              , intent(out) :: ebudget_loss2runoff
      real              , intent(out) :: co2budget_denseffect
      real              , intent(out) :: ebudget_denseffect
      real              , intent(out) :: wbudget_denseffect
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      integer                         :: ico
      integer                         :: ipft
      integer                         :: k
      integer                         :: ka
      integer                         :: kroot
      integer                         :: ksn
      integer                         :: kclosest
      integer                         :: nsoil
      real(kind=8)                    :: tmp_energy
      real(kind=8)                    :: available_water
      real(kind=8)                    :: gnd_water
      real(kind=8)                    :: psiplusz
      real(kind=8)                    :: mcheight
      real(kind=4)                    :: step_waterdef
      real(kind=4)                    :: can_rvap
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter         :: tendays_sec    = 10. * day_sec
      real        , parameter         :: thirtydays_sec = 30. * day_sec
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !------------------------------------------------------------------------------------!


      !----- Alias for the cohorts. -------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      !----- Alias for temporary surface water layers. ------------------------------------!
      ksn = initp%nlev_sfcwater
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
      csite%can_theta       (ipa) = sngloff(initp%can_theta       ,tiny_offset)
      csite%can_prss        (ipa) = sngloff(initp%can_prss        ,tiny_offset)
      csite%can_temp        (ipa) = sngloff(initp%can_temp        ,tiny_offset)
      csite%can_shv         (ipa) = sngloff(initp%can_shv         ,tiny_offset)
      csite%can_co2         (ipa) = sngloff(initp%can_co2         ,tiny_offset)
      csite%can_rhos        (ipa) = sngloff(initp%can_rhos        ,tiny_offset)
      csite%can_dmol        (ipa) = sngloff(initp%can_dmol        ,tiny_offset)
      csite%can_depth       (ipa) = sngloff(initp%can_depth       ,tiny_offset)
      csite%veg_displace    (ipa) = sngloff(initp%veg_displace    ,tiny_offset)
      csite%rough           (ipa) = sngloff(initp%rough           ,tiny_offset)
      csite%snowfac         (ipa) = sngloff(initp%snowfac         ,tiny_offset)
      csite%total_sfcw_depth(ipa) = sngloff(initp%total_sfcw_depth,tiny_offset)

      !------------------------------------------------------------------------------------!
      !    Find the ice-vapour equivalent potential temperature.  This is done outside the !
      ! integrator because it is an iterative method and currently we are not using it as  !
      ! a prognostic variable.                                                             !
      !------------------------------------------------------------------------------------!
      can_rvap                    = csite%can_shv(ipa) / ( 1.0 - csite%can_shv(ipa))
      csite%can_theiv(ipa)        = thetaeiv(csite%can_theta (ipa), csite%can_prss(ipa)    &
                                            ,csite%can_temp  (ipa), can_rvap               &
                                            ,can_rvap             )
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the vapour pressure deficit, which is diagnostic only.                     !
      !------------------------------------------------------------------------------------!
      csite%can_vpdef(ipa)        = vpdefil(csite%can_prss(ipa),csite%can_temp(ipa)        &
                                           ,csite%can_shv(ipa) ,.true.)
      !------------------------------------------------------------------------------------!



      !------ Copy the ground variables to the output. ------------------------------------!
      csite%ground_shv (ipa) = sngloff(initp%ground_shv , tiny_offset)
      csite%ground_ssh (ipa) = sngloff(initp%ground_ssh , tiny_offset)
      csite%ground_temp(ipa) = sngloff(initp%ground_temp, tiny_offset)
      csite%ground_fliq(ipa) = sngloff(initp%ground_fliq, tiny_offset)
      !------------------------------------------------------------------------------------!



      csite%ggbare(ipa)           = sngloff(initp%ggbare          ,tiny_offset)
      csite%ggveg (ipa)           = sngloff(initp%ggveg           ,tiny_offset)
      csite%ggnet (ipa)           = sngloff(initp%ggnet           ,tiny_offset)

      csite%ustar (ipa)           = sngloff(initp%ustar           ,tiny_offset)
      csite%tstar (ipa)           = sngloff(initp%tstar           ,tiny_offset)
      csite%qstar (ipa)           = sngloff(initp%qstar           ,tiny_offset)
      csite%cstar (ipa)           = sngloff(initp%cstar           ,tiny_offset)

      csite%zeta  (ipa)           = sngloff(initp%zeta            ,tiny_offset)
      csite%ribulk(ipa)           = sngloff(initp%ribulk          ,tiny_offset)

      csite%upwp  (ipa)           = sngloff(initp%upwp            ,tiny_offset)
      csite%wpwp  (ipa)           = sngloff(initp%wpwp            ,tiny_offset)
      csite%tpwp  (ipa)           = sngloff(initp%tpwp            ,tiny_offset)
      csite%qpwp  (ipa)           = sngloff(initp%qpwp            ,tiny_offset)
      csite%cpwp  (ipa)           = sngloff(initp%cpwp            ,tiny_offset)

      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics) then
         csite%fmean_vapor_gc        (ipa) = sngloff(initp%avg_vapor_gc       ,tiny_offset)
         csite%fmean_throughfall     (ipa) = sngloff(initp%avg_throughfall    ,tiny_offset)
         csite%fmean_vapor_ac        (ipa) = sngloff(initp%avg_vapor_ac       ,tiny_offset)
         csite%fmean_drainage        (ipa) = sngloff(initp%avg_drainage       ,tiny_offset)
         csite%fmean_qdrainage       (ipa) = sngloff(initp%avg_qdrainage      ,tiny_offset)
         csite%fmean_qthroughfall    (ipa) = sngloff(initp%avg_qthroughfall   ,tiny_offset)
         csite%fmean_sensible_gc     (ipa) = sngloff(initp%avg_sensible_gc    ,tiny_offset)
         csite%fmean_sensible_ac     (ipa) = sngloff(initp%avg_sensible_ac    ,tiny_offset)
         csite%fmean_carbon_ac       (ipa) = sngloff(initp%avg_carbon_ac      ,tiny_offset)
         csite%fmean_carbon_st       (ipa) = sngloff(initp%avg_carbon_st      ,tiny_offset)
         csite%fmean_ustar           (ipa) = sngloff(initp%avg_ustar          ,tiny_offset)
         csite%fmean_tstar           (ipa) = sngloff(initp%avg_tstar          ,tiny_offset)
         csite%fmean_qstar           (ipa) = sngloff(initp%avg_qstar          ,tiny_offset)
         csite%fmean_cstar           (ipa) = sngloff(initp%avg_cstar          ,tiny_offset)
         do k = rk4site%lsl, nzg
            csite%fmean_sensible_gg(k,ipa) = sngloff(initp%avg_sensible_gg(k) ,tiny_offset)
            csite%fmean_smoist_gg  (k,ipa) = sngloff(initp%avg_smoist_gg  (k) ,tiny_offset)
            csite%fmean_transloss  (k,ipa) = sngloff(initp%avg_transloss  (k) ,tiny_offset)
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Cohort-level variables.                                                     !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            cpatch%fmean_sensible_lc   (ico) = sngloff( initp%avg_sensible_lc   (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_sensible_wc   (ico) = sngloff( initp%avg_sensible_wc   (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_vapor_lc      (ico) = sngloff( initp%avg_vapor_lc      (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_vapor_wc      (ico) = sngloff( initp%avg_vapor_wc      (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_transp        (ico) = sngloff( initp%avg_transp        (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_intercepted_al(ico) = sngloff( initp%avg_intercepted_al(ico)      &
                                                      , tiny_offset)
            cpatch%fmean_intercepted_aw(ico) = sngloff( initp%avg_intercepted_aw(ico)      &
                                                      , tiny_offset)
            cpatch%fmean_wshed_lg      (ico) = sngloff( initp%avg_wshed_lg      (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_wshed_wg      (ico) = sngloff( initp%avg_wshed_wg      (ico)      &
                                                      , tiny_offset)
            !------------------------------------------------------------------------------!
            !     Plant hydraulic fluxes.  Convert them to kg/pl/s.                        !
            ! MLO: I kept the original units, although I would prefer to standardise all   !
            !      the fluxes to kg/m2/s.                                                  !
            !------------------------------------------------------------------------------!
            cpatch%fmean_wflux_wl      (ico) = sngloff( initp%avg_wflux_wl      (ico)      &
                                                      / initp%nplant            (ico)      &
                                                      , tiny_offset)
            cpatch%fmean_wflux_gw      (ico) = sngloff( initp%avg_wflux_gw      (ico)      &
                                                      / initp%nplant            (ico)      &
                                                      , tiny_offset)
            do k = rk4site%lsl, nzg
               cpatch%fmean_wflux_gw_layer(k,ico) =                                        &
                  sngloff( initp%avg_wflux_gw_layer(k,ico) / initp%nplant(ico)             &
                         , tiny_offset )
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      if(checkbudget) then
         co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_offset)
         co2budget_denseffect  = sngloff(initp%co2budget_denseffect ,tiny_offset)
         ebudget_netrad        = sngloff(initp%ebudget_netrad       ,tiny_offset)
         ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_offset)
         ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_offset)
         ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_offset)
         ebudget_denseffect    = sngloff(initp%ebudget_denseffect   ,tiny_offset)
         wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_offset)
         wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_offset)
         wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_offset)
         wbudget_denseffect    = sngloff(initp%wbudget_denseffect   ,tiny_offset)
      else
         co2budget_loss2atm    = 0.
         co2budget_denseffect  = 0.
         ebudget_netrad        = 0.
         ebudget_loss2atm      = 0.
         ebudget_loss2drainage = 0.
         ebudget_loss2runoff   = 0.
         ebudget_denseffect    = 0.
         wbudget_loss2atm      = 0.
         wbudget_loss2drainage = 0.
         wbudget_loss2runoff   = 0.
         wbudget_denseffect    = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is not a pure diagnostic, it is used for phenology and mortality !
      ! functions, preserve this variable and its dependencies in all contexts.            !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the water deficit.  This is done as a 30-day running average.           !
      !------------------------------------------------------------------------------------!
      step_waterdef                   = sngloff(initp%water_deficit,tiny_offset)
      csite%avg_monthly_waterdef(ipa) = csite%avg_monthly_waterdef(ipa) + step_waterdef
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     This variable is the monthly mean ground water that will be used to control    !
      ! fire disturbance.                                                                  !
      !------------------------------------------------------------------------------------!
      gnd_water = 0.d0
      !----- Add temporary surface water. -------------------------------------------------!
      do k=1,ksn
         gnd_water = gnd_water + initp%sfcwater_mass(k)
      end do
      !----- Find the bottommost layer to consider. ---------------------------------------!
      select case(include_fire)
      case (1)
         ka = rk4site%lsl
      case default
         ka = k_fire_first
      end select
      !----- Add soil moisture. -----------------------------------------------------------!
      do k=ka,nzg
         gnd_water = gnd_water + initp%soil_water(k) * dslz8(k) * wdns8
      end do
      !----- Add to the monthly mean. -----------------------------------------------------!
      csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa)                    &
                                      + sngloff(gnd_water,tiny_offset)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg - 10-day average of relative plant available water.  The relative value    !
      !           depends on whether the user wants to define phenology based on soil      !
      !           moisture or soil potential.                                              !
      !------------------------------------------------------------------------------------!
      if (spot_phen) then
         do ico = 1,cpatch%ncohorts
            ipft  = cpatch%pft(ico)
            kroot = cpatch%krdepth(ico)

            available_water = 0.d0
            do k = kroot, nzg
               nsoil            = rk4site%ntext_soil(k)
               mcheight         = 5.d-1 * ( dble(cpatch%height(ico))                       &
                                          + dble(h2crownbh(cpatch%height(ico),ipft)) )
               psiplusz         = slzt8(k) - mcheight                                      &
                                + matric_potential8(nsoil,initp%soil_water(k))
               available_water  = available_water                                          &
                                + max(0.d0,(psiplusz - soil8(nsoil)%slpotwp)) * dslz8(k)   &
                                / (soil8(nsoil)%slpotld - soil8(nsoil)%slpotwp)
            end do
            available_water     = available_water / abs(slz8(kroot))
            cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)         &
                                + sngl(available_water)*sngl(hdid)/tendays_sec
         end do
      else
         do ico = 1,cpatch%ncohorts
            available_water = 0.d0
            kroot           = cpatch%krdepth(ico)
            do k = kroot, nzg
               nsoil            = rk4site%ntext_soil(k)
               available_water  = available_water                                          &
                                + max(0.d0,(initp%soil_water(k)   - soil8(nsoil)%soilwp))  &
                                * dslz8(k) / (soil8(nsoil)%soilld - soil8(nsoil)%soilwp)
            end do
            available_water     = available_water / abs(slz8(kroot))
            cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)         &
                                + sngl(available_water)*sngl(hdid)/tendays_sec
         end do
      end if

      do k = rk4site%lsl, nzg
         csite%soil_water  (k,ipa) = sngloff(initp%soil_water  (k),tiny_offset)
         csite%soil_mstpot (k,ipa) = sngloff(initp%soil_mstpot (k),tiny_offset)
         csite%soil_energy (k,ipa) = sngloff(initp%soil_energy (k),tiny_offset)
         csite%soil_tempk  (k,ipa) = sngloff(initp%soil_tempk  (k),tiny_offset)
         csite%soil_fracliq(k,ipa) = sngloff(initp%soil_fracliq(k),tiny_offset)
      end do


      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m� inside the integrator. Convert it back !
      ! to J/kg in the layers that surface water/snow still exists.                        !
      !------------------------------------------------------------------------------------!
      csite%nlev_sfcwater(ipa)    = ksn
      csite%total_sfcw_depth(ipa) = 0.
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_depth(k,ipa)   = sngloff(initp%sfcwater_depth(k)   ,tiny_offset)
         csite%sfcwater_mass(k,ipa)    = sngloff(initp%sfcwater_mass(k)    ,tiny_offset)
         csite%sfcwater_tempk(k,ipa)   = sngloff(initp%sfcwater_tempk(k)   ,tiny_offset)
         csite%sfcwater_fracliq(k,ipa) = sngloff(initp%sfcwater_fracliq(k) ,tiny_offset)
         tmp_energy                    = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
         csite%sfcwater_energy(k,ipa)  = sngloff(tmp_energy                ,tiny_offset)
         csite%total_sfcw_depth(ipa)   =  csite%total_sfcw_depth(ipa)                      &
                                       +  csite%sfcwater_depth(k,ipa)
      end do
      !------------------------------------------------------------------------------------!
      !    For the layers that no longer exist, assign zeroes for prognostic variables,    !
      ! and something for temperature and liquid fraction (just to avoid singularities,    !
      ! and funny numbers in the output, but these values are meaningless and should never !
      ! be used).                                                                          !
      !------------------------------------------------------------------------------------!
      do k = csite%nlev_sfcwater(ipa)+1,nzs
         csite%sfcwater_energy(k,ipa)  = 0.
         csite%sfcwater_mass(k,ipa)    = 0.
         csite%sfcwater_depth(k,ipa)   = 0.
         if (k == 1) then
            csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
            csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk  (nzg,ipa)
         else
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            csite%sfcwater_tempk  (k,ipa) = csite%sfcwater_tempk  (k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Cohort variables.  Here we must check whether the cohort was really solved or  !
      ! it was skipped after being flagged as "unsafe".  In case the cohort was skipped,   !
      ! we must check whether it was because it was too small or because it was buried in  !
      ! snow.                                                                              !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !      First, update variables related with plant hydrodynamics because it can    !
         ! change leaf/wood heat capacity, which will be used later.                       !
         !---------------------------------------------------------------------------------!
         select case (plant_hydro_scheme)
         case (0)
            continue
         case default
            !----- Need to update leaf_water_im2 and wood_water_im2. ----------------------!
            cpatch%leaf_water_im2(ico) = sngloff(initp%leaf_water_im2(ico),tiny_offset)
            cpatch%wood_water_im2(ico) = sngloff(initp%wood_water_im2(ico),tiny_offset)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Update intensive internal water content.                                !
            !------------------------------------------------------------------------------!
            call twe2twi(cpatch%leaf_water_im2(ico),cpatch%wood_water_im2(ico)             &
                        ,cpatch%nplant(ico),cpatch%leaf_water_int(ico)                     &
                        ,cpatch%wood_water_int(ico))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Update rwc since it will be used to update leaf/wood heat capacity.     !
            !------------------------------------------------------------------------------!
            call tw2rwc(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)              &
                       ,cpatch%is_small(ico),cpatch%bleaf(ico),cpatch%bsapwooda(ico)       &
                       ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico)        &
                       ,cpatch%broot(ico),cpatch%dbh(ico),cpatch%pft(ico)                  &
                       ,cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Leaf and wood psi are updated in plant_hydro_driver of file             !
            ! ED/src/dynamics/plant_hydro.f90 for consistency reasons.  See that file for  !
            ! details.                                                                     !
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         select case (ibranch_thermo)
         case (1)
            !------------------------------------------------------------------------------!
            !  VEGETATION -- Leaf and branchwood were solved together, so they must remain !
            !                in thermal equilibrium.                                       !
            !------------------------------------------------------------------------------!
            if (initp%veg_resolvable(ico)) then

               !---------------------------------------------------------------------------!
               !     Copy vegetation wind.                                                 !
               !---------------------------------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    LEAVES.  It is always safe to copy internal energy and standing water, !
               !             but we must check whether leaves were truly resolved or not   !
               !             before copying the other variables.                           !
               !---------------------------------------------------------------------------!
               cpatch%leaf_water (ico) = sngloff(initp%leaf_water (ico) , tiny_offset)
               cpatch%leaf_energy(ico) = sngloff(initp%leaf_energy(ico) , tiny_offset)
               !---------------------------------------------------------------------------!


               if (initp%leaf_resolvable(ico)) then
                  !------------------------------------------------------------------------!
                  !    Leaves were solved, find the temperature and liquid fraction from   !
                  ! internal energy.                                                       !
                  !------------------------------------------------------------------------!
                  call uextcm2tl(cpatch%leaf_energy(ico)                                   &
                                ,cpatch%leaf_water(ico) + cpatch%leaf_water_im2(ico)       &
                                ,cpatch%leaf_hcap(ico),cpatch%leaf_temp(ico)               &
                                ,cpatch%leaf_fliq(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The intercellular specific humidity is always assumed to be at     !
                  ! saturation for a given temperature.  Find the saturation mixing ratio, !
                  ! then convert it to specific humidity.                                  !
                  !------------------------------------------------------------------------!
                  cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the leaf-level vapour pressure deficit using canopy pressure  !
                  ! and humitdity, but leaf temperature.                                   !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                  &
                                                  , cpatch%leaf_temp(ico)                  &
                                                  , csite%can_shv   (ipa), .true.)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Copy the conductances.                                             !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico), tiny_offset)
                  cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico), tiny_offset)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Divide the values of water demand by the time step to obtain the   !
                  ! average value over the past hdid period.                               !
                  !------------------------------------------------------------------------!
                  cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset)      &
                                         / sngl(hdid)
                  cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset)      &
                                         / sngl(hdid)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !    We solved leaf and branchwood together, the combined pool was re-   !
                  ! solvable but leaves weren't.  We copy the leaf temperature and liquid  !
                  ! fraction from the integrator, so they remain in thermal equilibrium    !
                  ! with branchwood.                                                       !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_temp(ico) = sngloff(initp%leaf_temp(ico) , tiny_offset)
                  cpatch%leaf_fliq(ico) = sngloff(initp%leaf_fliq(ico) , tiny_offset)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The intercellular specific humidity is always assumed to be at     !
                  ! saturation for a given temperature.  Find the saturation mixing ratio, !
                  ! then convert it to specific humidity.                                  !
                  !------------------------------------------------------------------------!
                  cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the leaf-level vapour pressure deficit using canopy pressure  !
                  ! and humitdity, but leaf temperature.                                   !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                  &
                                                  , cpatch%leaf_temp(ico)                  &
                                                  , csite%can_shv   (ipa), .true.)
                  !------------------------------------------------------------------------!


                  !----- Set water demand to zero. ----------------------------------------!
                  cpatch%psi_open  (ico) = 0.0
                  cpatch%psi_closed(ico) = 0.0
                  !------------------------------------------------------------------------!


                  !----- Leaf conductances cannot be zero.  Set to non-zero defaults. -----!
                  cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
                  cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    BRANCHES.  It is always safe to copy internal energy and standing      !
               !               water,  but we must check whether branches were truly       !
               !               resolved or not before copying the other variables.         !
               !---------------------------------------------------------------------------!
               cpatch%wood_water (ico) = sngloff(initp%wood_water (ico) , tiny_offset)
               cpatch%wood_energy(ico) = sngloff(initp%wood_energy(ico) , tiny_offset)
               if (initp%wood_resolvable(ico)) then
                  !------------------------------------------------------------------------!
                  !    Branches were solved, find the temperature and liquid fraction from !
                  ! internal energy.                                                       !
                  !------------------------------------------------------------------------!
                  call uextcm2tl(cpatch%wood_energy(ico)                                   &
                                ,cpatch%wood_water(ico) + cpatch%wood_water_im2(ico)       &
                                ,cpatch%wood_hcap(ico),cpatch%wood_temp(ico)               &
                                ,cpatch%wood_fliq(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Copy the conductances.                                             !
                  !------------------------------------------------------------------------!
                  cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico), tiny_offset)
                  cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico), tiny_offset)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !    We solved leaf and branchwood together, the combined pool was re-   !
                  ! solvable but leaves weren't.  We copy the leaf temperature and liquid  !
                  ! fraction from the integrator, so they remain in thermal equilibrium    !
                  ! with branchwood.                                                       !
                  !------------------------------------------------------------------------!
                  cpatch%wood_temp(ico) = sngloff(initp%wood_temp(ico) , tiny_offset)
                  cpatch%wood_fliq(ico) = sngloff(initp%wood_fliq(ico) , tiny_offset)
                  !------------------------------------------------------------------------!


                  !----- Wood conductances cannot be zero.  Set to non-zero defaults. -----!
                  cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
                  cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            elseif (cpatch%height(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the leaf and branch temperatures to the !
               ! snow temperature of the layer that is the closest to the cohort top.      !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%height(ico)) kclosest = k
               end do
               !---------------------------------------------------------------------------!


               cpatch%leaf_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               cpatch%wood_temp(ico)   = cpatch%leaf_temp(ico)

               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0.
               cpatch%wood_water(ico)  = 0.

               !---------------------------------------------------------------------------!
               !     Find the internal energy diagnostically...                            !
               !---------------------------------------------------------------------------!
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap     (ico)             &
                                                  , cpatch%leaf_water    (ico)             &
                                                  + cpatch%leaf_water_im2(ico)             &
                                                  , cpatch%leaf_temp     (ico)             &
                                                  , cpatch%leaf_fliq     (ico)             )
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap     (ico)             &
                                                  , cpatch%wood_water    (ico)             &
                                                  + cpatch%wood_water_im2(ico)             &
                                                  , cpatch%wood_temp     (ico)             &
                                                  , cpatch%wood_fliq     (ico)             )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!

               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Set water demand to zero. -------------------------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               cpatch%wood_gbw(ico) = cpatch%leaf_gbw(ico)
               cpatch%wood_gbh(ico) = cpatch%leaf_gbh(ico)
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! and branch temperatures to the canopy air space and force leaf and branch !
               ! intercepted water to be zero.                                             !
               !---------------------------------------------------------------------------!
               cpatch%leaf_temp(ico) = csite%can_temp(ipa)
               cpatch%wood_temp(ico) = cpatch%leaf_temp(ico)

               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico) = 0.5
                  cpatch%wood_fliq(ico) = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico) = 1.0
                  cpatch%wood_fliq(ico) = 1.0
               else
                  cpatch%leaf_fliq(ico) = 0.0
                  cpatch%wood_fliq(ico) = 0.0
               end if
               cpatch%leaf_water(ico)   = 0.
               cpatch%wood_water(ico)   = 0.
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the internal energy diagnostically...                            !
               !---------------------------------------------------------------------------!
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap     (ico)             &
                                                  , cpatch%leaf_water    (ico)             &
                                                  + cpatch%leaf_water_im2(ico)             &
                                                  , cpatch%leaf_temp     (ico)             &
                                                  , cpatch%leaf_fliq     (ico)             )
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap     (ico)             &
                                                  , cpatch%wood_water    (ico)             &
                                                  + cpatch%wood_water_im2(ico)             &
                                                  , cpatch%wood_temp     (ico)             &
                                                  , cpatch%wood_fliq     (ico)             )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Set water demand to zero. -------------------------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               cpatch%wood_gbw(ico) = cpatch%leaf_gbw(ico)
               cpatch%wood_gbh(ico) = cpatch%leaf_gbh(ico)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         case (0,2)
            !------------------------------------------------------------------------------!
            !  VEGETATION -- Leaf and branchwood were solved separately, so they are       !
            !                analysed independently.                                       !
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !  LEAVES                                                                      !
            !------------------------------------------------------------------------------!
            if (initp%leaf_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !     Leaves were solved, update water and internal energy, and re-         !
               ! calculate the temperature and leaf intercellular specific humidity.  The  !
               ! vegetation dry heat capacity is constant within one time step, so it      !
               ! doesn't need to be updated.                                               !
               !---------------------------------------------------------------------------!
               cpatch%leaf_water(ico)  = sngloff(initp%leaf_water(ico) , tiny_offset)
               cpatch%leaf_energy(ico) = sngloff(initp%leaf_energy(ico), tiny_offset)
               call uextcm2tl(cpatch%leaf_energy(ico)                                      &
                             ,cpatch%leaf_water(ico) + cpatch%leaf_water_im2(ico)          &
                             ,cpatch%leaf_hcap(ico),cpatch%leaf_temp(ico)                  &
                             ,cpatch%leaf_fliq(ico))

               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!



               !----- Convert the wind. ---------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Copy the conductances.                                                !
               !---------------------------------------------------------------------------!
               cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico), tiny_offset)
               cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico), tiny_offset)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Divide the values of water demand by the time step to obtain the      !
               ! average value over the past hdid period.                                  !
               !---------------------------------------------------------------------------!
               cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset)         &
                                      / sngl(hdid)
               cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset)         &
                                      / sngl(hdid)
               !---------------------------------------------------------------------------!

            elseif (cpatch%height(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the leaf temperature to the snow        !
               ! temperature of the layer that is the closest to the leaves.               !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%height(ico)) kclosest = k
               end do
               cpatch%leaf_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0.
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap     (ico)             &
                                                  , cpatch%leaf_water    (ico)             &
                                                  + cpatch%leaf_water_im2(ico)             &
                                                  , cpatch%leaf_temp     (ico)             &
                                                  , cpatch%leaf_fliq     (ico)             )
               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Set water demand to zero. -------------------------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! temperature to the canopy air space and force leaf_water to be zero.      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_temp(ico)   = csite%can_temp(ipa)
               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0.
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap     (ico)             &
                                                  , cpatch%leaf_water    (ico)             &
                                                  + cpatch%leaf_water_im2(ico)             &
                                                  , cpatch%leaf_temp     (ico)             &
                                                  , cpatch%leaf_fliq     (ico)             )
               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv  (ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind  (ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Set water demand to zero. -------------------------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!





            !------------------------------------------------------------------------------!
            !  WOOD                                                                        !
            !------------------------------------------------------------------------------!
            if (initp%wood_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !     Wood was solved, update water and internal energy, and recalculate    !
               ! the temperature.  The wood dry heat capacity is constant within one time  !
               ! step, so it doesn't need to be updated.                                   !
               !---------------------------------------------------------------------------!
               cpatch%wood_water(ico)  = sngloff(initp%wood_water(ico) , tiny_offset)
               cpatch%wood_energy(ico) = sngloff(initp%wood_energy(ico), tiny_offset)
               call uextcm2tl(cpatch%wood_energy(ico)                                      &
                             ,cpatch%wood_water(ico) + cpatch%wood_water_im2(ico)          &
                             ,cpatch%wood_hcap(ico),cpatch%wood_temp(ico)                  &
                             ,cpatch%wood_fliq(ico))

               !----- Convert the wind. ---------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Copy the conductances.                                                !
               !---------------------------------------------------------------------------!
               cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico), tiny_offset)
               cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico), tiny_offset)
               !---------------------------------------------------------------------------!

            elseif (cpatch%height(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the wood temperature to the snow        !
               ! temperature of the layer that is the closest to the branches.             !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%height(ico)) kclosest = k
               end do
               cpatch%wood_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               if (cpatch%wood_temp(ico) == t3ple) then
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%wood_temp(ico) > t3ple) then
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%wood_water(ico)  = 0.
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap     (ico)             &
                                                  , cpatch%wood_water    (ico)             &
                                                  + cpatch%wood_water_im2(ico)             &
                                                  , cpatch%wood_temp     (ico)             &
                                                  , cpatch%wood_fliq     (ico)             )
               !---------------------------------------------------------------------------!

               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     For very sparse patches of for when wood thermodynamics is off, fix   !
               ! the wood temperature to the canopy air space and force wood_water to be   !
               ! zero.                                                                     !
               !---------------------------------------------------------------------------!
               cpatch%wood_temp(ico)   = csite%can_temp(ipa)
               if (cpatch%wood_temp(ico) == t3ple) then
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%wood_temp(ico) > t3ple) then
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%wood_water(ico)  = 0.
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap     (ico)             &
                                                  , cpatch%wood_water    (ico)             &
                                                  + cpatch%wood_water_im2(ico)             &
                                                  , cpatch%wood_temp     (ico)             &
                                                  , cpatch%wood_fliq     (ico)             )
               !---------------------------------------------------------------------------!


               !----- Conductances cannot be zero.  Set to non-zero defaults. -------------!
               cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               if (.not. cpatch%leaf_resolvable(ico)) then
                  cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%leaf_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%leaf_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL LEAF_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:   ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:    ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:       ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:      ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:   '
            write (unit=*,fmt='(a,1x,es12.4)') '  - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '  - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '  - ENERGY:     ',initp%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER:      ',initp%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER_IM2:  ',initp%leaf_water_im2(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - TEMPERATURE:',initp%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - FRACLIQ:    ',initp%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - HEAT_CAP:   ',initp%leaf_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '  - ENERGY:     ',cpatch%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER:      ',cpatch%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER_IM2:  ',cpatch%leaf_water_im2(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - TEMPERATURE:',cpatch%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - FRACLIQ:    ',cpatch%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - HEAT_CAP:   ',cpatch%leaf_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_copy_patch.f90')
         end if
         if (cpatch%wood_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%wood_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL WOOD_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:   ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:    ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:       ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:      ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:   '
            write (unit=*,fmt='(a,1x,es12.4)') '  - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '  - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '  - ENERGY:     ',initp%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER:      ',initp%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER_IM2:  ',initp%wood_water_im2(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - TEMPERATURE:',initp%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - FRACLIQ:    ',initp%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - HEAT_CAP:   ',initp%wood_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '  - ENERGY:     ',cpatch%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER:      ',cpatch%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - WATER_IM2:  ',cpatch%wood_water_im2(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - TEMPERATURE:',cpatch%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - FRACLIQ:    ',cpatch%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '  - HEAT_CAP:   ',cpatch%wood_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_copy_patch.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Integrate the average state variables.  Notice that many variables (e.g.,     !
      ! temperature, density, and soil matric potential) are NOT integrated here: instead, !
      ! we find the averaged value after we normalise the average of the prognostic        !
      ! variables.  Same thing for aggregated variables at the patch level.                !
      !------------------------------------------------------------------------------------!
      csite%fmean_veg_displace(ipa) = csite%fmean_veg_displace(ipa)                        &
                                    + csite%veg_displace      (ipa) * dtlsm_o_frqsum
      csite%fmean_rough       (ipa) = csite%fmean_rough       (ipa)                        &
                                    + csite%rough             (ipa) * dtlsm_o_frqsum
      csite%fmean_can_theiv   (ipa) = csite%fmean_can_theiv   (ipa)                        &
                                    + csite%can_theiv         (ipa) * dtlsm_o_frqsum
      csite%fmean_can_theta   (ipa) = csite%fmean_can_theta   (ipa)                        &
                                    + csite%can_theta         (ipa) * dtlsm_o_frqsum
      csite%fmean_can_vpdef   (ipa) = csite%fmean_can_vpdef   (ipa)                        &
                                    + csite%can_vpdef         (ipa) * dtlsm_o_frqsum
      csite%fmean_can_shv     (ipa) = csite%fmean_can_shv     (ipa)                        &
                                    + csite%can_shv           (ipa) * dtlsm_o_frqsum
      csite%fmean_can_co2     (ipa) = csite%fmean_can_co2     (ipa)                        &
                                    + csite%can_co2           (ipa) * dtlsm_o_frqsum
      csite%fmean_can_prss    (ipa) = csite%fmean_can_prss    (ipa)                        &
                                    + csite%can_prss          (ipa) * dtlsm_o_frqsum
      csite%fmean_gnd_temp    (ipa) = csite%fmean_gnd_temp    (ipa)                        &
                                    + csite%ground_temp       (ipa) * dtlsm_o_frqsum
      csite%fmean_gnd_shv     (ipa) = csite%fmean_gnd_shv     (ipa)                        &
                                    + csite%ground_shv        (ipa) * dtlsm_o_frqsum
      csite%fmean_can_ggnd    (ipa) = csite%fmean_can_ggnd    (ipa)                        &
                                    + csite%ggnet             (ipa) * dtlsm_o_frqsum
      csite%fmean_snowfac     (ipa) = csite%fmean_snowfac     (ipa)                        &
                                    + csite%snowfac           (ipa) * dtlsm_o_frqsum
      !------------------------------------------------------------------------------------!
      !       Snow/pounding layers.  We keep track of the total, not individual layers.    !
      ! Energy will be integrated as an extensive variable, we will convert it by the      !
      ! output time only.                                                                  !
      !------------------------------------------------------------------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%fmean_sfcw_depth (ipa) = csite%fmean_sfcw_depth   (ipa)                     &
                                      + csite%sfcwater_depth   (k,ipa) * dtlsm_o_frqsum
         csite%fmean_sfcw_energy(ipa) = csite%fmean_sfcw_energy  (ipa)                     &
                                      + csite%sfcwater_energy  (k,ipa)                     &
                                      * csite%sfcwater_mass    (k,ipa) * dtlsm_o_frqsum
         csite%fmean_sfcw_mass  (ipa) = csite%fmean_sfcw_mass    (ipa)                     &
                                      + csite%sfcwater_mass    (k,ipa) * dtlsm_o_frqsum
      end do
      !------ Cohort-level variables. -----------------------------------------------------!
      do ico=1,cpatch%ncohorts
         cpatch%fmean_leaf_energy   (ico) = cpatch%fmean_leaf_energy   (ico)               &
                                          + cpatch%leaf_energy         (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_water    (ico) = cpatch%fmean_leaf_water    (ico)               &
                                          + cpatch%leaf_water          (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_water_im2(ico) = cpatch%fmean_leaf_water_im2(ico)               &
                                          + cpatch%leaf_water_im2      (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_hcap     (ico) = cpatch%fmean_leaf_hcap     (ico)               &
                                          + cpatch%leaf_hcap           (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_vpdef    (ico) = cpatch%fmean_leaf_vpdef    (ico)               &
                                          + cpatch%leaf_vpdef          (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_energy   (ico) = cpatch%fmean_wood_energy   (ico)               &
                                          + cpatch%wood_energy         (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_water    (ico) = cpatch%fmean_wood_water    (ico)               &
                                          + cpatch%wood_water          (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_water_im2(ico) = cpatch%fmean_wood_water_im2(ico)               &
                                          + cpatch%wood_water_im2      (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_hcap     (ico) = cpatch%fmean_wood_hcap     (ico)               &
                                          + cpatch%wood_hcap           (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_gsw      (ico) = cpatch%fmean_leaf_gsw      (ico)               &
                                          + cpatch%leaf_gsw            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_gbw      (ico) = cpatch%fmean_leaf_gbw      (ico)               &
                                          + cpatch%leaf_gbw            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_gbw      (ico) = cpatch%fmean_wood_gbw      (ico)               &
                                          + cpatch%wood_gbw            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_psi_open      (ico) = cpatch%fmean_psi_open      (ico)               &
                                          + cpatch%psi_open            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_psi_closed    (ico) = cpatch%fmean_psi_closed    (ico)               &
                                          + cpatch%psi_closed          (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_fs_open       (ico) = cpatch%fmean_fs_open       (ico)               &
                                          + cpatch%fs_open             (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_fsw           (ico) = cpatch%fmean_fsw           (ico)               &
                                          + cpatch%fsw                 (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_fsn           (ico) = cpatch%fmean_fsn           (ico)               &
                                          + cpatch%fsn                 (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_open        (ico) = cpatch%fmean_A_open        (ico)               &
                                          + cpatch%A_open              (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_closed      (ico) = cpatch%fmean_A_closed      (ico)               &
                                          + cpatch%A_closed            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_net         (ico) = cpatch%fmean_A_net         (ico)               &
                                          + ( ( 1. - cpatch%fs_open    (ico) )             &
                                            * cpatch%A_closed          (ico)               &
                                            + cpatch%fs_open           (ico)               &
                                            * cpatch%A_open            (ico) )               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_light       (ico) = cpatch%fmean_A_light       (ico)               &
                                          + cpatch%A_light             (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_rubp        (ico) = cpatch%fmean_A_rubp        (ico)               &
                                          + cpatch%A_rubp              (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_A_co2         (ico) = cpatch%fmean_A_co2         (ico)               &
                                          + cpatch%A_co2               (ico)               &
                                          * dtlsm_o_frqsum
         !---------------------------------------------------------------------------------!
         !     The penalty factor for water and nitrogen are meaningful only during the    !
         ! day.  For the daily means we must add only when it is daytime, so we integrate  !
         ! them here too.                                                                  !
         !---------------------------------------------------------------------------------!
         if (.not. nighttime .and. writing_long) then
            cpatch%dmean_fs_open    (ico) = cpatch%dmean_fs_open    (ico)                  &
                                          + cpatch%fs_open          (ico)                  &
                                          * dtlsm
            cpatch%dmean_fsw        (ico) = cpatch%dmean_fsw        (ico)                  &
                                          + cpatch%fsw              (ico)                  &
                                          * dtlsm
            cpatch%dmean_fsn        (ico) = cpatch%dmean_fsn        (ico)                  &
                                          + cpatch%fsn              (ico)                  &
                                          * dtlsm
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------ Soil variables. -------------------------------------------------------------!
      do k = rk4site%lsl, nzg
         csite%fmean_soil_energy(k,ipa) = csite%fmean_soil_energy(k,ipa)                   &
                                        + csite%soil_energy      (k,ipa)                   &
                                        * dtlsm_o_frqsum
         csite%fmean_soil_mstpot(k,ipa) = csite%fmean_soil_mstpot(k,ipa)                   &
                                        + csite%soil_mstpot      (k,ipa)                   &
                                        * dtlsm_o_frqsum
         csite%fmean_soil_water (k,ipa) = csite%fmean_soil_water (k,ipa)                   &
                                        + csite%soil_water       (k,ipa)                   &
                                        * dtlsm_o_frqsum
      end do
      !------------------------------------------------------------------------------------!


     return
   end subroutine initp2modelp
   !=======================================================================================!
   !=======================================================================================!

end module rk4_copy_patch
!==========================================================================================!
!==========================================================================================!
