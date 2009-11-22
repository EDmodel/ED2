!==========================================================================================!
!==========================================================================================!
!    This subroutine copies that variables that are integrated by the Runge-Kutta solver   !
! to a buffer structure.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine copy_patch_init(sourcesite,ipa,targetp)
   use ed_state_vars        , only : sitetype              & ! structure
                                   , patchtype             ! ! structure
   use grid_coms            , only : nzg                   & ! intent(in)
                                   , nzs                   ! ! intent(in) 
   use ed_misc_coms         , only : fast_diagnostics      ! ! intent(in)
   use consts_coms          , only : cpi8                  & ! intent(in)
                                   , ep8                   & ! intent(in)
                                   , cp8                   & ! intent(in)
                                   , epim18                & ! intent(in)
                                   , alvl8                 & ! intent(in)
                                   , rdry8                 & ! intent(in)
                                   , rdryi8                & ! intent(in)
                                   , p00i8                 & ! intent(in)
                                   , rocp8                 ! ! intent(in)
   use rk4_coms             , only : rk4patchtype          & ! structure
                                   , rk4site               & ! structure
                                   , hcapveg_ref           & ! intent(in)
                                   , rk4eps                & ! intent(in)
                                   , min_height            & ! intent(in)
                                   , any_solvable          & ! intent(out)
                                   , zoveg                 & ! intent(out)
                                   , zveg                  & ! intent(out)
                                   , wcapcan               & ! intent(out)
                                   , wcapcani              & ! intent(out)
                                   , rk4water_stab_thresh  & ! intent(in)
                                   , rk4min_sfcwater_mass  & ! intent(in)
                                   , checkbudget           ! ! intent(in)
   use ed_max_dims          , only : n_pft                 ! ! intent(in)
   use canopy_radiation_coms, only : tai_min               ! ! intent(in)
   use therm_lib8           , only : qwtk8                 & ! subroutine
                                   , ptqz2enthalpy8        & ! function
                                   , idealdenssh8          & ! function
                                   , reducedpress8         ! ! function
   use allometry            , only : dbh2bl                ! ! function
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)    , target     :: targetp
   type(sitetype)        , target     :: sourcesite
   integer               , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch
   real(kind=8)                       :: hvegpat_min
   real(kind=8)                       :: hcap_scale
   real(kind=8)                       :: baux,xaux,raux,jaux
   integer                            :: ico
   integer                            :: ipft
   integer                            :: k
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Between time steps the pressure may change because of change in atmospheric       !
   ! pressure, which means that enthalpy is not conserved.  Potential temperature, on the  !
   ! other hand, is conserved because there is no heat flux between time steps.  So we use !
   ! this instead to start all other variables.                                            !
   !---------------------------------------------------------------------------------------!
   !----- 1. Update thermo variables that are conserved between steps. --------------------!
   targetp%can_theta    = dble(sourcesite%can_theta(ipa))
   targetp%can_shv      = dble(sourcesite%can_shv(ipa))
   targetp%can_co2      = dble(sourcesite%can_co2(ipa))
   targetp%can_depth    = dble(sourcesite%can_depth(ipa))
   !----- 2. Update the canopy pressure based on the new atmospheric pressure. ------------!
   targetp%can_prss     = reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv &
                                       ,rk4site%geoht,targetp%can_theta,targetp%can_shv    &
                                       ,targetp%can_depth)
   !----- 3. Update temperature, enthalpy, and density. -----------------------------------!
   targetp%can_temp     = targetp%can_theta * (p00i8 * targetp%can_prss) ** rocp8
   targetp%can_enthalpy = ptqz2enthalpy8(targetp%can_prss,targetp%can_temp                 &
                                        ,targetp%can_shv,targetp%can_depth)
   targetp%can_rhos     = idealdenssh8(targetp%can_prss,targetp%can_temp,targetp%can_shv)
   !---------------------------------------------------------------------------------------!

   do k = rk4site%lsl, nzg
      targetp%soil_water(k)   = dble(sourcesite%soil_water(k,ipa))
      targetp%soil_energy(k)  = dble(sourcesite%soil_energy(k,ipa))
      targetp%soil_tempk(k)   = dble(sourcesite%soil_tempk(k,ipa))
      targetp%soil_fracliq(k) = dble(sourcesite%soil_fracliq(k,ipa))
   end do

   do k = 1, nzs
      targetp%sfcwater_mass(k)    = dble(sourcesite%sfcwater_mass(k,ipa))
      targetp%sfcwater_depth(k)   = dble(sourcesite%sfcwater_depth(k,ipa))
      !----- Converting sfcwater_energy to J/m² inside the Runge-Kutta integrator. --------!
      targetp%sfcwater_energy(k)  = dble(sourcesite%sfcwater_energy(k,ipa))                &
                                  * dble(sourcesite%sfcwater_mass(k,ipa))
      targetp%sfcwater_tempk(k)   = dble(sourcesite%sfcwater_tempk(k,ipa))
      targetp%sfcwater_fracliq(k) = dble(sourcesite%sfcwater_fracliq(k,ipa))
   end do

   targetp%ustar = dble(sourcesite%ustar(ipa))
   targetp%cstar = dble(sourcesite%cstar(ipa))
   targetp%tstar = dble(sourcesite%tstar(ipa))
   targetp%qstar = dble(sourcesite%qstar(ipa))
   targetp%estar = 0.d0

   targetp%upwp = dble(sourcesite%upwp(ipa))
   targetp%wpwp = dble(sourcesite%wpwp(ipa))
   targetp%tpwp = dble(sourcesite%tpwp(ipa))
   targetp%qpwp = dble(sourcesite%qpwp(ipa))
   targetp%cpwp = dble(sourcesite%cpwp(ipa))

  
   targetp%nlev_sfcwater = sourcesite%nlev_sfcwater(ipa)


   !----- The virtual pools should be always zero, they are temporary entities ------------!
   targetp%virtual_water = 0.0d0
   targetp%virtual_heat  = 0.0d0
   targetp%virtual_depth = 0.0d0

   if (targetp%nlev_sfcwater == 0) then
      targetp%virtual_flag = 2
   else
      if (targetp%sfcwater_mass(1) < rk4min_sfcwater_mass) then
         targetp%virtual_flag = 2
      elseif (targetp%sfcwater_mass(1) < rk4water_stab_thresh) then
         targetp%virtual_flag = 1
      else
         targetp%virtual_flag = 0
      end if
   end if

   !---------------------------------------------------------------------------------------!
   !     Here we find the minimum patch-level leaf heat capacity.  If the total patch leaf !
   ! heat capacity is less than this, we scale the cohorts heat capacity inside the        !
   ! integrator, so it preserves the proportional heat capacity and prevents the pool to   !
   ! be too small.                                                                         !
   !---------------------------------------------------------------------------------------!
   cpatch => sourcesite%patch(ipa)
   sourcesite%hcapveg(ipa) = 0.
   sourcesite%lai(ipa)     = 0.
   sourcesite%wpa(ipa)     = 0.
   sourcesite%wai(ipa)     = 0.
   do ico=1,cpatch%ncohorts
      sourcesite%hcapveg(ipa) = sourcesite%hcapveg(ipa) + cpatch%hcapveg(ico)
      sourcesite%lai(ipa)     = sourcesite%lai(ipa)     + cpatch%lai(ico)
      sourcesite%wpa(ipa)     = sourcesite%wpa(ipa)     + cpatch%wpa(ico)
      sourcesite%wai(ipa)     = sourcesite%wai(ipa)     + cpatch%wai(ico)
   end do
   
   any_solvable = .false.
   do ico=1, cpatch%ncohorts
      !----- Copying the flag that determines whether this cohort is numerically stable. --!
      targetp%solvable(ico) = cpatch%solvable(ico)
      if (targetp%solvable(ico)) any_solvable = .true.
   end do

   if ((sourcesite%lai(ipa)+sourcesite%wai(ipa)) > tai_min) then
      hvegpat_min = hcapveg_ref * max(dble(cpatch%hite(1)),min_height)
      hcap_scale  = max(1.d0,hvegpat_min / sourcesite%hcapveg(ipa))
   else
      hcap_scale  = 1.d0
   end if

   do ico = 1,cpatch%ncohorts
      ipft=cpatch%pft(ico)
      !----- Copy the leaf area index and total (leaf+branch+twig) area index. ------------!
      targetp%lai(ico)    = dble(cpatch%lai(ico))
      targetp%wai(ico)    = dble(cpatch%wai(ico))
      targetp%wpa(ico)    = dble(cpatch%wpa(ico))
      targetp%tai(ico)    = targetp%lai(ico) + dble(cpatch%wai(ico))

      !------------------------------------------------------------------------------------!
      !    If the cohort is too small, we give some extra heat capacity, so the model can  !
      ! run in a stable range inside the integrator.  At the end this extra heat capacity  !
      ! will be removed.                                                                   !
      !------------------------------------------------------------------------------------!
      targetp%hcapveg(ico) = dble(cpatch%hcapveg(ico)) * hcap_scale

      !------------------------------------------------------------------------------------!
      !     Checking whether this is considered a "safe" one or not.  In case it is, we    !
      ! copy water, temperature, and liquid fraction, and scale energy and heat capacity   !
      ! as needed.  Otherwise, just fill with some safe values, but the cohort won't be    !
      ! really solved.                                                                     !
      !------------------------------------------------------------------------------------!
      targetp%veg_water(ico)     = dble(cpatch%veg_water(ico))

      if (targetp%solvable(ico)) then
         targetp%veg_energy(ico)    = dble(cpatch%veg_energy(ico))                         &
                                    + (targetp%hcapveg(ico)-dble(cpatch%hcapveg(ico)))     &
                                    * dble(cpatch%veg_temp(ico))
         call qwtk8(targetp%veg_energy(ico),targetp%veg_water(ico),targetp%hcapveg(ico)    &
                   ,targetp%veg_temp(ico),targetp%veg_fliq(ico))
      else
         targetp%veg_fliq(ico)   = dble(cpatch%veg_fliq(ico))
         targetp%veg_temp(ico)   = dble(cpatch%veg_temp(ico))
         targetp%veg_energy(ico) = targetp%hcapveg(ico) * targetp%veg_temp(ico)
      end if
   end do

   !----- Diagnostics variables -----------------------------------------------------------!
   if(fast_diagnostics) then
      !----------------------------------------------------------------------!
      !   N.B. The "budget" variables are not copied here because they are   !
      ! integrated outside RK4.  Inside RK4 we only want the contribution    !
      ! of those variables during the span of one time step.                 !
      !----------------------------------------------------------------------!
      targetp%avg_carbon_ac      = dble(sourcesite%avg_carbon_ac(ipa)     )
      targetp%avg_vapor_vc       = dble(sourcesite%avg_vapor_vc(ipa)      )
      targetp%avg_dew_cg         = dble(sourcesite%avg_dew_cg(ipa)        )
      targetp%avg_vapor_gc       = dble(sourcesite%avg_vapor_gc(ipa)      )
      targetp%avg_wshed_vg       = dble(sourcesite%avg_wshed_vg(ipa)      )
      targetp%avg_intercepted    = dble(sourcesite%avg_intercepted(ipa)   )
      targetp%avg_vapor_ac       = dble(sourcesite%avg_vapor_ac(ipa)      )
      targetp%avg_transp         = dble(sourcesite%avg_transp(ipa)        )
      targetp%avg_evap           = dble(sourcesite%avg_evap(ipa)          )
      targetp%avg_drainage       = dble(sourcesite%avg_drainage(ipa)      )
      targetp%avg_drainage_heat  = dble(sourcesite%avg_drainage_heat(ipa) )
      targetp%avg_netrad         = dble(sourcesite%avg_netrad(ipa)        )
      targetp%avg_sensible_vc    = dble(sourcesite%avg_sensible_vc(ipa)   )
      targetp%avg_qwshed_vg      = dble(sourcesite%avg_qwshed_vg(ipa)     )
      targetp%avg_qintercepted   = dble(sourcesite%avg_qintercepted(ipa)  )
      targetp%avg_sensible_gc    = dble(sourcesite%avg_sensible_gc(ipa)   )
      targetp%avg_sensible_ac    = dble(sourcesite%avg_sensible_ac(ipa)   )

      do k = rk4site%lsl, nzg
         targetp%avg_sensible_gg(k) = dble(sourcesite%avg_sensible_gg(k,ipa))
         targetp%avg_smoist_gg(k)   = dble(sourcesite%avg_smoist_gg(k,ipa)  )
         targetp%avg_smoist_gc(k)   = dble(sourcesite%avg_smoist_gc(k,ipa)  )
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
      targetp%ebudget_latent        = 0.d0
      targetp%wbudget_loss2atm      = 0.d0
      targetp%wbudget_loss2drainage = 0.d0
      targetp%wbudget_loss2runoff   = 0.d0
   end if

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
      !----- Copy the plant density. ------------------------------------------------------!
      targetp%nplant(ico) = dble(cpatch%nplant(ico))

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
   use rk4_coms              , only : rk4site              & ! intent(in)
                                    , rk4min_sfcwater_mass & ! intent(in)
                                    , rk4min_can_shv       & ! intent(in)
                                    , rk4min_can_temp      & ! intent(in)
                                    , rk4max_can_shv       & ! intent(in)
                                    , rk4patchtype         ! ! structure
   use ed_state_vars         , only : sitetype             & ! structure
                                    , patchtype            ! ! structure
   use soil_coms             , only : soil8                ! ! intent(in)
   use grid_coms             , only : nzg                  & ! intent(in)
                                    , nzs                  ! ! intent(in)
   use therm_lib8            , only : qwtk8                & ! subroutine
                                    , qtk8                 & ! subroutine
                                    , hpqz2temp8           & ! function
                                    , ptqz2enthalpy8       & ! function
                                    , idealdenssh8         ! ! function
   use consts_coms           , only : alvl8                & ! intent(in)
                                    , wdns8                & ! intent(in)
                                    , rdryi8               & ! intent(in)
                                    , rdry8                & ! intent(in)
                                    , epim18               & ! intent(in)
                                    , toodry8              & ! intent(in)
                                    , cp8                  & ! intent(in)
                                    , cpi8                 & ! intent(in)
                                    , p008                 & ! intent(in)
                                    , rocp8                ! ! intent(in)
   use canopy_struct_dynamics, only : can_whcap8           ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer :: cpatch
   integer                          :: ico
   integer                          :: k
   real(kind=8)                     :: soilhcap
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Here we convert enthalpy into temperature, potential temperature, and density.    !
   !---------------------------------------------------------------------------------------!
   if (initp%can_shv >= rk4min_can_shv .and. initp%can_shv <= rk4max_can_shv) then 
      initp%can_temp  = hpqz2temp8(initp%can_enthalpy,initp%can_prss,initp%can_shv         &
                                  ,initp%can_depth)
   else
      initp%can_temp  = rk4min_can_temp
   end if 
   initp%can_theta = initp%can_temp * (p008 / initp%can_prss) ** rocp8
   initp%can_rhos  = idealdenssh8(initp%can_prss,initp%can_temp,initp%can_shv)
   !---------------------------------------------------------------------------------------!

   !----- Updating soil temperature and liquid water fraction. ----------------------------!
   do k = rk4site%lsl, nzg - 1
      soilhcap = soil8(csite%ntext_soil(k,ipa))%slcpd
      call qwtk8(initp%soil_energy(k),initp%soil_water(k)*wdns8,soilhcap                   &
                ,initp%soil_tempk(k),initp%soil_fracliq(k))
   end do
   !---------------------------------------------------------------------------------------!

   call can_whcap8(csite,ipa,initp%can_rhos,initp%can_depth)

   !---------------------------------------------------------------------------------------!
   !    Updating surface water temperature and liquid water fraction, remembering that in- !
   ! side the RK4 integration, surface water energy is in J/m². The abs is necessary be-   !
   ! cause surface mass may indeed become too negative during the integration process and  !
   ! if it happens, we want the step to be rejected.                                       !
   !---------------------------------------------------------------------------------------!
   do k = 1, nzs
      if(abs(initp%sfcwater_mass(k)) > rk4min_sfcwater_mass)  then
           call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                       &
                    ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      elseif (k == 1) then
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
         initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
      else
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
         initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
      end if
   end do


   cpatch => csite%patch(ipa)

   !----- Looping over cohorts ------------------------------------------------------------!
   cohortloop: do ico=1,cpatch%ncohorts
      !----- Checking whether this is a prognostic cohort... ------------------------------!
      if (initp%solvable(ico)) then
         !----- Lastly we update leaf temperature and liquid fraction. --------------------!
         call qwtk8(initp%veg_energy(ico),initp%veg_water(ico),initp%hcapveg(ico)          &
                   ,initp%veg_temp(ico),initp%veg_fliq(ico))
      end if

   end do cohortloop

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
subroutine redistribute_snow(initp,csite,ipa)

   use rk4_coms      , only : rk4patchtype         & ! structure
                            , rk4min_sfcw_mass     & ! intent(in)
                            , rk4min_virt_water    & ! intent(in)
                            , rk4water_stab_thresh & ! intent(in)
                            , rk4min_sfcwater_mass & ! intent(in)
                            , rk4snowmin           & ! intent(in)
                            , newsnow              ! ! intent(in)
   use ed_state_vars , only : sitetype             & ! structure
                            , patchtype            ! ! structure
   use grid_coms     , only : nzs                  & ! intent(in)
                            , nzg                  ! ! intent(in)
   use soil_coms     , only : soil8                & ! intent(in)
                            , dslz8                & ! intent(in)
                            , dslzi8               & ! intent(in)
                            , thick                & ! intent(in)
                            , thicknet             ! ! intent(in)
   use consts_coms   , only : cice8                & ! intent(in)
                            , cliq8                & ! intent(in)
                            , t3ple8               & ! intent(in)
                            , wdns8                & ! intent(in)
                            , wdnsi8               & ! intent(in)
                            , tsupercool8          & ! intent(in)
                            , qliqt38              & ! intent(in)
                            , wdnsi8               ! ! intent(in)
   use therm_lib8    , only : qtk8                 & ! subroutine
                            , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp
   type(sitetype)         , target     :: csite
   integer                , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: kold
   integer                             :: newlayers
   integer                             :: nlayers
   integer                             :: ksn
   integer                             :: ksnnew
   integer                             :: k
   !----- Control variables ---------------------------------------------------------------!
   real(kind=8)                        :: wtold
   real(kind=8)                        :: wtnew
   real(kind=8), dimension(nzs)        :: newsfcw_mass
   real(kind=8), dimension(nzs)        :: newsfcw_energy
   real(kind=8), dimension(nzs)        :: newsfcw_depth
   real(kind=8)                        :: wdiff
   real(kind=8)                        :: totsnow
   real(kind=8)                        :: depthgain
   real(kind=8)                        :: wfree
   real(kind=8)                        :: qwfree
   real(kind=8)                        :: qw
   real(kind=8)                        :: w
   real(kind=8)                        :: dw
   real(kind=8)                        :: wtemp
   real(kind=8)                        :: wfliq
   real(kind=8)                        :: wfreeb
   real(kind=8)                        :: depthloss
   real(kind=8)                        :: snden
   real(kind=8)                        :: sndenmin
   real(kind=8)                        :: sndenmax
   real(kind=8)                        :: qwt
   real(kind=8)                        :: wt
   real(kind=8)                        :: soilhcap
   real(kind=8)                        :: free_surface_water_demand
   real(kind=8)                        :: Cr    ! snow waterholding capacity
   real(kind=8)                        :: gi    ! Partial density of ice
   integer                             :: nsoil
   !----- Constants -----------------------------------------------------------------------!
   logical                , parameter  :: debug   = .false.
   real(kind=8)           , parameter  :: Crmin   = 3.d-2
   real(kind=8)           , parameter  :: Crmax   = 1.d-1
   real(kind=8)           , parameter  :: ge      = 2.d2
   
   !---------------------------------------------------------------------------------------!


   !----- Initializing # of layers alias --------------------------------------------------!
   ksn       = initp%nlev_sfcwater

   if (ksn >= 1) then
      !------------------------------------------------------------------------------------!
      ! 1. There used to exist temporary water/snow layers here.  Check total mass to see  !
      !    whether there is still enough mass.                                             !
      !------------------------------------------------------------------------------------!
      totsnow = sum(initp%sfcwater_mass(1:ksn))
      if (totsnow < rk4min_sfcw_mass) then
         !----- Temporary layer is too negative, break it so the step can be rejected. ----!
         return
      elseif (totsnow <= rk4min_sfcwater_mass) then
         !---------------------------------------------------------------------------------!
         ! 1.a. Too little or negative mass.  Eliminate layers, ensuring that it will  not !
         !      leak mass or energy, by "stealing" them from the top soil layer.           !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_energy(1) = sum(initp%sfcwater_energy(1:ksn))
         initp%sfcwater_mass(1)   = sum(initp%sfcwater_mass(1:ksn))
         initp%soil_energy(nzg)   = initp%soil_energy(nzg)                                 &
                                  + initp%sfcwater_energy(1) * dslzi8(nzg)
         initp%soil_water(nzg)    = initp%soil_water(nzg)                                  &
                                  + initp%sfcwater_mass(1)   * dslzi8(nzg) * wdnsi8
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         initp%sfcwater_mass      = 0.d0
         initp%sfcwater_energy    = 0.d0
         initp%sfcwater_tempk     = initp%soil_tempk(nzg)
         initp%sfcwater_fracliq   = 0.d0
         initp%sfcwater_depth     = 0.d0       
         initp%nlev_sfcwater      = 0
         ksnnew                   = 0
      else
         !---------------------------------------------------------------------------------!
         ! 1.b.  Still something there, nothing changes at least not for the time being.   !
         !---------------------------------------------------------------------------------!
         ksnnew = ksn
         wfree               = 0.d0
         qwfree              = 0.d0
         depthgain           = 0.d0
      end if
   else
      !------------------------------------------------------------------------------------!
      ! 2.  No temporary layer, dealing with virtual layer.  Check whether the virtual     !
      !     layer would be thick enough to create a pond, otherwise skip the entire thing. !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water < rk4min_virt_water) then
         !----- Virtual layer is too negative, break it so the step can be rejected. ------!
         return
      elseif (initp%virtual_water <= rk4min_sfcwater_mass) then
         !---------------------------------------------------------------------------------!
         ! 2.a. Too little or negative mass in the virtual layer.  No layer will be creat- !
         !      ed, but before eliminating it, just make sure mass and energy will be      !
         !      conserved.                                                                 !
         !---------------------------------------------------------------------------------!
         ksnnew = 0
         initp%soil_energy(nzg)   = initp%soil_energy(nzg)                                 &
                                  + initp%virtual_heat * dslzi8(nzg)
         initp%soil_water(nzg)    = initp%soil_water(nzg)                                  &
                                  + initp%virtual_water * dslzi8(nzg) * wdnsi8
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         initp%virtual_water      = 0.d0
         initp%virtual_heat       = 0.d0
         initp%virtual_depth      = 0.d0
      else
         !---------------------------------------------------------------------------------!
         ! 2.b. No temporary layer, significant mass to add.  ksnnew will be at least one. !
         !      If there was no layer before, create one.                                  !
         !---------------------------------------------------------------------------------!
         wfree               = initp%virtual_water
         qwfree              = initp%virtual_heat
         depthgain           = initp%virtual_depth
         initp%virtual_water = 0.d0
         initp%virtual_heat  = 0.d0
         initp%virtual_depth = 0.d0
         ksnnew = 1
      end if
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 3. We now update the diagnostic variables, and ensure the layers are stable.  Loop    !
   !    over layers, from top to bottom this time.                                         !
   !---------------------------------------------------------------------------------------!
   totsnow = 0.d0
   do k = ksnnew,1,-1

      !----- Update current mass and energy of temporary layer ----------------------------!
      qw = initp%sfcwater_energy(k) + qwfree
      w  = initp%sfcwater_mass(k)   + wfree
      dw = initp%sfcwater_depth(k)  + depthgain

      !------------------------------------------------------------------------------------!
      !    Single layer, and this is a very thin one, which can cause numerical instabili- !
      ! ty. Force a fast heat exchange between this thin layer and the soil topmost level, !
      ! bringing both layers to a thermal equilibrium.                                     !
      !------------------------------------------------------------------------------------!
      if (ksnnew == 1 .and. initp%sfcwater_mass(k) < rk4water_stab_thresh) then

         !---------------------------------------------------------------------------------!
         !     Total internal energy and water of the combined system, in J/m² and kg/m²,  !
         ! respectively.                                                                   !
         !---------------------------------------------------------------------------------!
         qwt = qw + initp%soil_energy(nzg) * dslz8(nzg)
         wt  = w  + initp%soil_water(nzg)  * dslz8(nzg) * wdns8

         !----- Finding the equilibrium temperature and liquid/ice partition. -------------!
         soilhcap = soil8(csite%ntext_soil(nzg,ipa))%slcpd * dslz8(nzg)
         call qwtk8(qwt,wt,soilhcap,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))

         !---------------------------------------------------------------------------------!
         !    Computing internal energy of the temporary layer with the temperature and    !
         ! liquid/ice distribution we just found, for the mass the layer has.              !
         !---------------------------------------------------------------------------------!
         qw = w                                                                            &
            * (initp%sfcwater_fracliq(k) * cliq8 *(initp%sfcwater_tempk(k)-tsupercool8)    &
                  + (1.d0-initp%sfcwater_fracliq(k)) * cice8 * initp%sfcwater_tempk(k))

         !---------------------------------------------------------------------------------!
         !    Set the properties of top soil layer. Since internal energy is an extensive  !
         ! quantity, we can simply take the difference to be the soil internal energy,     !
         ! just remembering that we need to convert it back to J/m³. The other properties  !
         ! can be copied from the surface layer because we assumed phase and temperature   !
         ! equilibrium.                                                                    !
         !---------------------------------------------------------------------------------!
         initp%soil_energy(nzg)  = (qwt - qw) * dslzi8(nzg)
         initp%soil_tempk(nzg)   = initp%sfcwater_tempk(k)
         initp%soil_fracliq(nzg) = initp%sfcwater_fracliq(k)
      else
         !----- Layer is computationally stable, just update the temperature and phase ----!
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         call qtk8(qw/w,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      end if


      !------------------------------------------------------------------------------------!
      !    Shed liquid in excess of a 1:9 liquid-to-ice ratio through percolation.  Limit  !
      ! this shed amount (wfreeb) in lowest snow layer to amount top soil layer can hold.  !
      !------------------------------------------------------------------------------------!
      if (w > rk4min_sfcwater_mass) then

         !----- Find the fraction of liquid water of the shed amount. ---------------------!
         call qtk8(qw/w,wtemp,wfliq)

         if(newsnow) then

            !------------------------------------------------------------------------------!
            !    Alternative "free" water calculation.                                     !
            !    Anderson (1976), NOAA Tech Report NWS 19.                                 !
            !------------------------------------------------------------------------------!
            gi = w/dw * (1.d0 - wfliq)

            if (gi < ge) then
               Cr = Crmin + (Crmax - Crmin) * (ge - gi) / ge
            else
               Cr = Crmin
            end if

            if (wfliq > Cr / (1.d0 + Cr)) then
               wfreeb = w * (wfliq - Cr / (1.d0 + Cr))
            else
               wfreeb = 0.d0
            end if
         else
            wfreeb = max(0.d0, w * (wfliq-1.d-1)/9.d-1)
         end if
      else
         wfreeb = 0.0
      end if


      if (k == 1)then
           !----- Do "greedy" infiltration. -----------------------------------------------!
           nsoil = csite%ntext_soil(nzg,ipa)
           free_surface_water_demand = max(0.d0                                            &
                                          ,soil8(nsoil)%slmsts - initp%soil_water(nzg))    &
                                     * wdns8 * dslz8(nzg)
           wfreeb = min(wfreeb,free_surface_water_demand)
           qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(k)-tsupercool8)
           !----- Update topmost soil moisture and energy, updating temperature and phase -!
           initp%soil_water(nzg)  = initp%soil_water(nzg)  + wfreeb*wdnsi8*dslzi8(nzg)
           initp%soil_energy(nzg) = initp%soil_energy(nzg) + qwfree * dslzi8(nzg)
           soilhcap = soil8(nsoil)%slcpd
           call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                   &
                     ,soilhcap,initp%soil_tempk(nzg),initp%soil_fracliq(nzg))
      else
         !---- Not the first layer, just shed all free water, and compute its energy ------!
         qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(k)-tsupercool8)
      end if
      depthloss = wfreeb * wdnsi8

      !----- Remove water and internal energy losses due to percolation -------------------!
      initp%sfcwater_mass(k)  = w - wfreeb
      initp%sfcwater_depth(k) = initp%sfcwater_depth(k) + depthgain - depthloss
      if(initp%sfcwater_mass(k) > rk4min_sfcwater_mass) then
         initp%sfcwater_energy(k) = qw - qwfree
         call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k),initp%sfcwater_tempk(k) &
                 ,initp%sfcwater_fracliq(k))
      else
         initp%sfcwater_energy(k) = 0.d0
         initp%sfcwater_mass(k)   = 0.d0
         initp%sfcwater_depth(k)  = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end if

      !----- Integrate total "snow" -------------------------------------------------------!
      totsnow = totsnow + initp%sfcwater_mass(k)

      !----- Calculate density and depth of snow ------------------------------------------!
      snden    = initp%sfcwater_mass(k) / max(1.0d-6,initp%sfcwater_depth(k))
      sndenmax = wdns8
      sndenmin = max(3.d1, 2.d2 * (wfree + wfreeb)                                         &
               / max(rk4min_sfcwater_mass,initp%sfcwater_mass(k)))
      snden    = min(sndenmax, max(sndenmin,snden))
      initp%sfcwater_depth(k) = initp%sfcwater_mass(k) / snden

      !----- Set up input to next layer ---------------------------------------------------!
      wfree = wfreeb
      depthgain = depthloss
   end do

   !---------------------------------------------------------------------------------------!
   ! 4. Re-distribute snow layers to maintain prescribed distribution of mass.             !
   !---------------------------------------------------------------------------------------!
   if (totsnow <= rk4min_sfcwater_mass .or. ksnnew == 0) then
      initp%nlev_sfcwater = 0
      !----- Making sure that the unused layers have zero in everything -------------------!
      do k = 1, nzs
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do
   else
      !---- Check whether there is enough snow for a new layer. ---------------------------!
      nlayers   = ksnnew
      newlayers = 1
      do k = 1,nzs
         !----- Checking whether we need 
         if (      initp%sfcwater_mass(k)   >  rk4min_sfcwater_mass                        &
             .and. rk4snowmin * thicknet(k) <= totsnow                                     &
             .and. initp%sfcwater_energy(k) <  initp%sfcwater_mass(k)*qliqt38 ) then

            newlayers = newlayers + 1
         end if
      end do
      newlayers = min(newlayers, nzs, nlayers + 1)
      initp%nlev_sfcwater = newlayers
      kold  = 1
      wtnew = 1.d0
      wtold = 1.d0
      do k = 1,newlayers
         newsfcw_mass(k)   = totsnow * thick(k,newlayers)
         newsfcw_energy(k) = 0.d0
         newsfcw_depth(k)  = 0.d0
         !----- Finding new layer properties ----------------------------------------------!
         find_layer: do

            !----- Difference between old and new snow ------------------------------------!
            wdiff = wtnew * newsfcw_mass(k) - wtold * initp%sfcwater_mass(kold)  

            if (wdiff > 0.d0) then
               newsfcw_energy(k) = newsfcw_energy(k) + wtold * initp%sfcwater_energy(kold)
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtold * initp%sfcwater_depth(kold)
               wtnew  = wtnew - wtold * initp%sfcwater_mass(kold) / newsfcw_mass(k)
               kold   = kold + 1
               wtold  = 1.0
               if (kold > nlayers) exit find_layer
            else
               newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_energy(kold)                             &
                                 / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_depth(kold)                              &
                                 / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtold = wtold - wtnew * newsfcw_mass(k)                                     &
                             / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtnew = 1.
               exit find_layer
            end if
         end do find_layer
      end do

      !----- Updating the water/snow layer properties -------------------------------------!
      do k = 1,newlayers
         initp%sfcwater_mass(k)   = newsfcw_mass(k)
         initp%sfcwater_energy(k) = newsfcw_energy(k)
         initp%sfcwater_depth(k)  = newsfcw_depth(k)
         if (newsfcw_mass(k) > rk4min_sfcwater_mass) then
            call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                      &
                    ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
         elseif (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do

      !----- Making sure that the unused layers have zero in everything -------------------!
      do k = newlayers + 1, nzs
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do
   end if

   return
end subroutine redistribute_snow
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
   use rk4_coms             , only : rk4patchtype     & ! structure
                                   , rk4site          & ! intent(in)
                                   , rk4eps           & ! intent(in)
                                   , rk4min_sfcw_mass & ! intent(in)
                                   , rk4min_can_shv   & ! intent(in)
                                   , wcapcan          & ! intent(in)
                                   , wcapcani         ! ! intent(in)
   use ed_state_vars        , only : sitetype         & ! structure
                                   , patchtype        ! ! structure
   use consts_coms          , only : cice8            & ! intent(in)
                                   , cliq8            & ! intent(in)
                                   , alvl8            & ! intent(in)
                                   , alvi8            & ! intent(in)
                                   , alli8            & ! intent(in)
                                   , t3ple8           & ! intent(in)
                                   , wdns8            & ! intent(in)
                                   , idns8            & ! intent(in)
                                   , toodry8          & ! intent(in)
                                   , tsupercool8      & ! intent(in)
                                   , qliqt38          & ! intent(in)
                                   , wdnsi8           ! ! intent(in)
   use therm_lib8           , only : qwtk8            & ! subroutine
                                   , qtk8             ! ! subroutine
   use grid_coms            , only : nzg              ! ! intent(in)
   use soil_coms            , only : soil8            & ! intent(in)
                                   , dslzi8           & ! intent(in)
                                   , dslz8            ! ! intent(in)
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
   real(kind=8)                        :: virtual_temp
   real(kind=8)                        :: virtual_fliq
   logical                             :: slightlymoist
   logical                             :: slightlydry
   !---------------------------------------------------------------------------------------!

   !----- Inverse of time increment -------------------------------------------------------!
   hdidi = 1.d0 / hdid
   
   !----- Defining some aliases that will be often used during the integration. -----------!
   kt            = nzg
   nstop         = csite%ntext_soil(kt,ipa)

   !---------------------------------------------------------------------------------------!
   !     Now we check if the soil moisture is within safe bounds or seriously messed up.   !
   ! In boths case we won't bother adjusting.                                              !
   !---------------------------------------------------------------------------------------!
   if (initp%soil_water(kt) >= soil8(nstop)%soilcp .and.                                   &
       initp%soil_water(kt) <= soil8(nstop)%slmsts) then 
      !----- We are in the normal range, return... ----------------------------------------!
      return
   end if

   !----- Check whether we are just slightly off. -----------------------------------------!
   slightlymoist = initp%soil_water(kt) > soil8(nstop)%slmsts .and.                        &
                   initp%soil_water(kt) <= (1.d0 + rk4eps) * soil8(nstop)%slmsts
   slightlydry   = initp%soil_water(kt) < soil8(nstop)%soilcp .and.                        &
                   initp%soil_water(kt) >= (1.d0 - rk4eps) * soil8(nstop)%soilcp

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
         energy_needed = water_needed * initp%virtual_heat   / water_available
         depth_needed  = water_needed * initp%virtual_depth  / water_available

         !---------------------------------------------------------------------------------!
         !     Add the water and energy into the top layer, remove it from the surface,    !
         ! and quit.                                                                       !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)  = initp%soil_water(kt)  + water_needed  * dslzi8(kt)*wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)

         initp%virtual_depth   = initp%virtual_depth  - depth_needed
         initp%virtual_water   = initp%virtual_water  - water_needed
         initp%virtual_heat    = initp%virtual_heat   - energy_needed
         return
      elseif (water_available > 0.d0) then
         !---------------------------------------------------------------------------------!
         !     This water will not be enough, but we use it and seek for other sources to  !
         ! fill the remainder.                                                             !
         !---------------------------------------------------------------------------------!
         water_needed = water_needed - water_available
         initp%soil_water(kt)  = initp%soil_water(kt)                                      &
                               + initp%virtual_water * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt) = initp%soil_energy(kt) + initp%virtual_heat * dslzi8(kt)

         !----- Remove the virtual layer. -------------------------------------------------!
         initp%virtual_water  = 0.d0
         initp%virtual_heat   = 0.d0
         initp%virtual_depth  = 0.d0
      end if

      !------------------------------------------------------------------------------------!
      !    If we hit this point the temporary layers were not enough.  We now look for     !
      ! water in the layer immediately beneath the top soil layer.  We now find the amount !
      ! of water this layer can donate.  We will also need energy associated with the , so !
      ! we compute the temperature and liquid fraction.                                    !
      !------------------------------------------------------------------------------------!
      kb              = nzg -1 
      nsbeneath       = csite%ntext_soil(kb,ipa)
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
         initp%can_enthalpy    = initp%can_enthalpy   - energy_needed * wcapcani
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
         initp%can_enthalpy    = initp%can_enthalpy   - energy_available * wcapcani

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
                                     + (1.d0 - initp%soil_fracliq(kt)) * idns8)
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
      elseif (initp%virtual_water + water_excess > rk4min_sfcw_mass) then
         !---------------------------------------------------------------------------------!
         !     If the virtual layer will have some significant mass after adding the water !
         ! excess, we transfer the water to there.  It will likely become the new          !
         ! temporary surface water layer.                                                  !
         !---------------------------------------------------------------------------------!
         initp%soil_water(kt)     = initp%soil_water(kt)                                   &
                                  - water_excess * dslzi8(kt) * wdnsi8
         initp%soil_energy(kt)    = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
         
         initp%virtual_water      = initp%virtual_water     + water_excess
         initp%virtual_heat       = initp%virtual_heat      + energy_excess
         initp%virtual_depth      = initp%virtual_depth     + depth_excess

         return
      end if

      !------------------------------------------------------------------------------------!
      !    If we hit this point, then the amount of water is small but it can't go to the  !
      ! preferred destinations.  Adding on virtual pool wouldn't help because the water    !
      ! would be sent back to the top soil layer in redistribute_snow call.  Thus the      !
      ! excess now becomes the excess of water in the layer, plus any water left in the    !
      ! virtual layer...                                                                   !
      !    Anyway, we first eliminate any water left in the virtual layer by boiling it to !
      ! the atmosphere.  This is a tiny amount and even if supersaturation occurs, it      !
      ! shouldn't be enough to cause trouble.                                              !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water > 0.d0) then 
         call qtk8(initp%virtual_heat/initp%virtual_water,virtual_temp,virtual_fliq)
         initp%can_shv      = initp%can_shv      + initp%virtual_water  * wcapcani 

         initp%can_enthalpy = initp%can_enthalpy                                           &
                            + initp%virtual_water * (alvi8 - virtual_fliq * alli8)         &
                            * wcapcani

         initp%avg_vapor_gc = initp%avg_vapor_gc + initp%virtual_water  * hdidi
         
         !----- Say goodbye to the virtual layer... ---------------------------------------!
         initp%virtual_heat   = 0.d0
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
      nsbeneath  = csite%ntext_soil(kb,ipa)
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
         initp%soil_water(kb)  = initp%soil_water(kb)  + energy_excess * dslzi8(kb)
         !----- Update the fluxes too... --------------------------------------------------!
         initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) - water_excess * hdidi
         initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) + water_excess * hdidi
         !---------------------------------------------------------------------------------!

         return
      else
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
         initp%soil_water(kb)  = initp%soil_water(kb)  + energy_room * dslzi8(kb)
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
         initp%can_enthalpy    = initp%can_enthalpy    + energy_excess * wcapcani
         initp%can_shv         = initp%can_shv         + water_excess  * wcapcani
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
                                   , wcapcani           & ! intent(in)
                                   , rk4dry_veg_lwater  & ! intent(in)
                                   , rk4fullveg_lwater  ! ! intent(in)
   use ed_state_vars        , only : sitetype           & ! structure
                                   , patchtype          ! ! structure
   use consts_coms          , only : cice8              & ! intent(in)
                                   , cliq8              & ! intent(in)
                                   , alvl8              & ! intent(in)
                                   , alvi8              & ! intent(in)
                                   , alli8              & ! intent(in)
                                   , t3ple8             & ! intent(in)
                                   , wdns8              & ! intent(in)
                                   , idns8              & ! intent(in)
                                   , tsupercool8        & ! intent(in)
                                   , qliqt38            & ! intent(in)
                                   , wdnsi8             ! ! intent(in)
   use therm_lib8           , only : qwtk8              ! ! subroutine
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
   real(kind=8)                        :: min_leaf_water
   real(kind=8)                        :: max_leaf_water
   real(kind=8)                        :: veg_wshed
   real(kind=8)                        :: veg_qwshed
   real(kind=8)                        :: veg_dwshed
   real(kind=8)                        :: veg_dew
   real(kind=8)                        :: veg_qdew
   real(kind=8)                        :: hdidi
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)
   
   !----- Inverse of time increment -------------------------------------------------------!
   hdidi = 1.d0 / hdid

   !----- Looping over cohorts ------------------------------------------------------------!
   cohortloop: do ico=1,cpatch%ncohorts
      !----- Checking whether this is a prognostic cohort... ------------------------------!
      if (initp%solvable(ico)) then
         !---------------------------------------------------------------------------------!
         !   Now we find the maximum leaf water possible.                                  !
         !---------------------------------------------------------------------------------!
         rk4min_leaf_water = rk4min_veg_lwater * initp%tai(ico)
         min_leaf_water    = rk4dry_veg_lwater * initp%tai(ico)
         max_leaf_water    = rk4fullveg_lwater * initp%tai(ico)

         !------ Leaf water is too negative, break it so the step can be rejected. --------!
         if (initp%veg_water(ico) < rk4min_leaf_water) then
            return
         !----- Shedding excessive water to the ground ------------------------------------!
         elseif (initp%veg_water(ico) > max_leaf_water) then
            veg_wshed  = (initp%veg_water(ico)-max_leaf_water)
            veg_qwshed = veg_wshed                                                         &
                       * (initp%veg_fliq(ico) * cliq8 * (initp%veg_temp(ico)-tsupercool8)    &
                         + (1.d0-initp%veg_fliq(ico)) * cice8 * initp%veg_temp(ico))
            veg_dwshed = veg_wshed                                                         &
                       / (initp%veg_fliq(ico) * wdns8 + (1.d0-initp%veg_fliq(ico))*idns8)

            !----- Updating water mass and energy. ----------------------------------------!
            initp%veg_water(ico)  = initp%veg_water(ico)  - veg_wshed
            initp%veg_energy(ico) = initp%veg_energy(ico) - veg_qwshed
            
            !----- Updating virtual pool --------------------------------------------------!
            ksn = initp%nlev_sfcwater
            if (ksn > 0) then
               initp%sfcwater_mass(ksn)   = initp%sfcwater_mass(ksn)   + veg_wshed
               initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn) + veg_qwshed
               initp%sfcwater_depth(ksn)  = initp%sfcwater_depth(ksn)  + veg_dwshed
            else
               initp%virtual_water   = initp%virtual_water + veg_wshed
               initp%virtual_heat    = initp%virtual_heat  + veg_qwshed
               initp%virtual_depth   = initp%virtual_depth + veg_dwshed
            end if
            !----- Updating output fluxes -------------------------------------------------!
            initp%avg_wshed_vg     = initp%avg_wshed_vg     + veg_wshed  * hdidi
            initp%avg_qwshed_vg    = initp%avg_qwshed_vg    + veg_qwshed * hdidi
            initp%avg_intercepted  = initp%avg_intercepted  - veg_wshed  * hdidi
            initp%avg_qintercepted = initp%avg_qintercepted - veg_qwshed * hdidi

         !---------------------------------------------------------------------------------!
         !    If veg_water is tiny, exchange moisture with the air.  If the tiny number is !
         ! positive, then "donate" the remaining as "boiling" (fast evaporation or         !
         ! sublimation).  Likewise, if veg_water is tiny and negative, exchange moisture   !
         ! with the air, "stealing" moisture as fast "dew/frost" condensation.             !
         !---------------------------------------------------------------------------------!
         elseif (initp%veg_water(ico) < min_leaf_water) then
            veg_dew = - initp%veg_water(ico)
            veg_qdew = veg_dew * (alvi8 - initp%veg_fliq(ico) * alli8)

            !----- Updating state variables -----------------------------------------------!
            initp%veg_water(ico)  = 0.d0
            initp%veg_energy(ico) = initp%veg_energy(ico)  + veg_qdew
            initp%can_shv         = initp%can_shv          - veg_dew  * wcapcani
            initp%can_enthalpy    = initp%can_enthalpy     - veg_qdew * wcapcani

            !----- Updating output flux ---------------------------------------------------!
            initp%avg_vapor_vc    = initp%avg_vapor_vc   - veg_dew * hdidi
            initp%ebudget_latent  = initp%ebudget_latent - veg_qdew
         end if

         !----- Lastly we update leaf temperature and liquid fraction. --------------------!
         call qwtk8(initp%veg_energy(ico),initp%veg_water(ico),initp%hcapveg(ico)          &
                   ,initp%veg_temp(ico),initp%veg_fliq(ico))
      end if

   end do cohortloop

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

   errmax       = max(0.0,abs(yerr%can_enthalpy/yscal%can_enthalpy))
   troublemaker = large_error(yerr%can_enthalpy,yscal%can_enthalpy)
   write(unit=*,fmt=onefmt) 'CAN_ENTHALPY:',errmax,yerr%can_enthalpy,yscal%can_enthalpy    &
                                           ,troublemaker

   errmax       = max(errmax,abs(yerr%can_shv/yscal%can_shv))
   troublemaker = large_error(yerr%can_shv,yscal%can_shv)
   write(unit=*,fmt=onefmt) 'CAN_SHV:',errmax,yerr%can_shv,yscal%can_shv,troublemaker

   errmax = max(errmax,abs(yerr%can_co2/yscal%can_co2))
   troublemaker = large_error(yerr%can_co2,yscal%can_co2)
   write(unit=*,fmt=onefmt) 'CAN_CO2:',errmax,yerr%can_co2,yscal%can_co2,troublemaker

   errmax = max(errmax,abs(yerr%virtual_heat/yscal%virtual_heat))
   troublemaker = large_error(yerr%virtual_heat,yscal%virtual_heat)
   write(unit=*,fmt=onefmt) 'VIRTUAL_HEAT:',errmax,yerr%virtual_heat,yscal%virtual_heat    &
                                           ,troublemaker

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
   write(unit=*,fmt='(a)'      ) ' Cohort_level variables (only the solvable ones):'
   write(unit=*,fmt='(10(a,1x))')        'Name            ','   PFT','         LAI'        &
                                      ,'         WAI','         WPA','         TAI'        &
                                      ,'   Max.Error','   Abs.Error','       Scale'        &
                                      ,'Problem(T|F)'
   do ico = 1,cpatch%ncohorts
      if (y%solvable(ico)) then
         errmax       = max(errmax,abs(yerr%veg_water(ico)/yscal%veg_water(ico)))
         troublemaker = large_error(yerr%veg_water(ico),yscal%veg_water(ico))
         write(unit=*,fmt=cohfmt) 'VEG_WATER:',cpatch%pft(ico),y%lai(ico),y%wai(ico)       &
                                              ,y%wpa(ico),y%tai(ico),errmax                &
                                              ,yerr%veg_water(ico),yscal%veg_water(ico)    &
                                              ,troublemaker
              

         errmax       = max(errmax,abs(yerr%veg_energy(ico)/yscal%veg_energy(ico)))
         troublemaker = large_error(yerr%veg_energy(ico),yscal%veg_energy(ico))
         write(unit=*,fmt=cohfmt) 'VEG_ENERGY:',cpatch%pft(ico),cpatch%lai(ico),y%wai(ico) &
                                               ,y%wpa(ico),y%tai(ico),errmax               &
                                               ,yerr%veg_energy(ico),yscal%veg_energy(ico) &
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

      errmax = max(errmax                                                                  &
                  ,abs(yerr%ebudget_latent/yscal%ebudget_latent))
      troublemaker = large_error(yerr%ebudget_latent                                       &
                                ,yscal%ebudget_latent)
      write(unit=*,fmt=onefmt) 'EN_LATENT:',errmax,yerr%ebudget_latent                     &
                              ,yscal%ebudget_latent,troublemaker


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
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: ipa
   !----- Local variable ------------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   integer                      :: ico 
   integer                      :: k   
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (csite) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'UTC'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
 
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only the solvable ones shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),11(a12,1x))')                                              &
         '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH','       BDEAD'   &
                           &,'      BALIVE','  VEG_ENERGY','    VEG_TEMP','   VEG_WATER'   &
                           &,'     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),11(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico) &
              ,cpatch%nplant(ico),cpatch%lai(ico),cpatch%dbh(ico),cpatch%bdead(ico)        &
              ,cpatch%balive(ico),cpatch%veg_energy(ico),cpatch%veg_temp(ico)              &
              ,cpatch%veg_water(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '   DIST_TYPE','         AGE','        AREA'          &
                                   &,'          RH','AVGDAILY_TMP','     SUM_CHD'          &
                                   &,'     SUM_DGD'
   write (unit=*,fmt='(i12,1x,6(es12.4,1x))')  csite%dist_type(ipa),csite%age(ipa)         &
         ,csite%area(ipa),csite%rh(ipa),csite%avg_daily_temp(ipa),csite%sum_chd(ipa)       &
         ,csite%sum_dgd(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(11(a12,1x))') '  VEG_HEIGHT','   VEG_ROUGH','         LAI'          &
                                    ,'        HTRY','     CAN_CO2','    CAN_TEMP'          &
                                    ,'     CAN_SHV','    CAN_RHOS','    CAN_PRSS'          &
                                    ,'CAN_ENTHALPY','   CAN_DEPTH'
   write (unit=*,fmt='(11(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)         &
         ,csite%lai(ipa),csite%htry(ipa),csite%can_co2(ipa),csite%can_temp(ipa)            &
         ,csite%can_shv(ipa),csite%can_rhos(ipa),csite%can_prss(ipa)                       &
         ,csite%can_enthalpy(ipa),csite%can_depth(ipa) 

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                   &,'       TSTAR','     RLONG_G','    RSHORT_G'          &
                                   &,'     RLONG_S'
   write (unit=*,fmt='(7(es12.4,1x))') csite%ustar(ipa),csite%qstar(ipa),csite%cstar(ipa)  &
         ,csite%tstar(ipa),csite%rlong_g(ipa),csite%rshort_g(ipa),csite%rlong_s(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,a12)') '  PFT','       REPRO'
   do k=1,n_pft
      write (unit=*,fmt='(i5,1x,es12.4)') k,csite%repro(k,ipa)
   end do

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4site%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,csite%ntext_soil(k,ipa)              &
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
   use rk4_coms              , only : rk4patchtype         & ! structure
                                    , rk4site              & ! intent(in)
                                    , rk4min_sfcwater_mass ! ! intent(in)
   use ed_state_vars         , only : sitetype             & ! structure
                                    , patchtype            ! ! structure
   use grid_coms             , only : nzg                  & ! intent(in)
                                    , nzs                  ! ! intent(in)
   use ed_misc_coms             , only : current_time         ! ! intent(in)
   use therm_lib8            , only : qtk8                 & ! subroutine
                                    , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: y
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer    :: cpatch
   integer                         :: k
   integer                         :: ico
   real(kind=8)                    :: virtual_temp, virtual_fliq
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (rk4patch) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'s'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(80a)')         ('-',k=1,80)
   write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
   write (unit=*,fmt='(a,1x,es12.4)') ' Air temperature    : ',rk4site%atm_tmp
   write (unit=*,fmt='(a,1x,es12.4)') ' Air enthalpy       : ',rk4site%atm_enthalpy
   write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov mixing ratio  : ',rk4site%atm_shv
   write (unit=*,fmt='(a,1x,es12.4)') ' CO2  mixing ratio  : ',rk4site%atm_co2
   write (unit=*,fmt='(a,1x,es12.4)') ' Pressure           : ',rk4site%atm_prss
   write (unit=*,fmt='(a,1x,es12.4)') ' Exner function     : ',rk4site%atm_exner
   write (unit=*,fmt='(a,1x,es12.4)') ' Wind speed         : ',rk4site%vels
   write (unit=*,fmt='(a,1x,es12.4)') ' Height             : ',rk4site%geoht
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. mass  flux : ',rk4site%pcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. heat  flux : ',rk4site%qpcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. depth flux : ',rk4site%dpcpg

   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only those solvable are shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','        HITE','         DBH','       BDEAD'   &
                           &,'      BALIVE','     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdead(ico)       &
              ,cpatch%balive(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','         WPA','         TAI','  VEG_ENERGY'   &
             ,'   VEG_WATER','    VEG_HCAP','    VEG_TEMP','    VEG_FLIQ'
   do ico = 1,cpatch%ncohorts
      if (y%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
               ,y%lai(ico),y%wpa(ico),y%tai(ico),y%veg_energy(ico),y%veg_water(ico)        &
               ,y%hcapveg(ico),y%veg_temp(ico),y%veg_fliq(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(11(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','   PATCH_LAI'         &
                                     ,'     CAN_CO2','CAN_ENTHALPY','   CAN_THETA'         &
                                     ,'    CAN_TEMP','     CAN_SHV','    CAN_RHOS'         &
                                     ,'    CAN_PRSS','   CAN_DEPTH'
                                     
   write (unit=*,fmt='(11(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)         &
                                       ,csite%lai(ipa),y%can_co2,y%can_enthalpy            &
                                       ,y%can_theta,y%can_temp,y%can_shv,y%can_rhos        &
                                       ,y%can_prss,y%can_depth

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(5(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                    ,'       TSTAR','       ESTAR'
   write (unit=*,fmt='(5(es12.4,1x))') y%ustar,y%qstar,y%cstar,y%tstar,y%estar

   write (unit=*,fmt='(80a)') ('-',k=1,80)
   if (y%virtual_water /= 0.) then
      call qtk8(y%virtual_heat/y%virtual_water,virtual_temp,virtual_fliq)
   else
      virtual_temp = y%soil_tempk(nzg)
      virtual_fliq = y%soil_fracliq(nzg)
   end if


   write (unit=*,fmt='(5(a12,1x))')  'VIRTUAL_FLAG','VIRTUAL_HEAT','  VIRT_WATER'          &
                                   &,'VIRTUAL_TEMP','VIRTUAL_FLIQ'
   write (unit=*,fmt='(i12,1x,4(es12.4,1x))') y%virtual_flag,y%virtual_heat                &
                                             ,y%virtual_water,virtual_temp,virtual_fliq
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(4(a12,1x))')    '  GROUND_SHV',' SURFACE_SSH','SURFACE_TEMP'        &
                                      ,'SURFACE_FLIQ'
   write (unit=*,fmt='(4(es12.4,1x))') y%ground_shv, y%surface_ssh, y%surface_temp         &
                                      ,y%surface_fliq

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4site%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,csite%ntext_soil(k,ipa)              &
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
