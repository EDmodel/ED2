!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns an initial value of zero for most cohort-level variables.    !
! There will be a few variables that must be initialised with a value other than zero, in  !
! which case there will be a note explaining the reason.                                   !
!------------------------------------------------------------------------------------------!
subroutine init_ed_cohort_vars(cpatch,ico, lsl)
   use ed_state_vars , only : patchtype           ! ! structure
   use allometry     , only : dbh2krdepth         ! ! function
   use pft_coms      , only : phenology           & ! intent(in)
                            , leaf_turnover_rate  & ! intent(in)
                            , Vm0                 & ! intent(in)
                            , sla                 ! ! intent(in)
   use ed_misc_coms  , only : writing_long        & ! intent(in)
                            , writing_eorq        & ! intent(in)
                            , writing_dcyc        ! ! intent(in)
   use phenology_coms, only : vm0_tran            & ! intent(in)
                            , vm0_slope           & ! intent(in)
                            , vm0_amp             & ! intent(in)
                            , vm0_min             & ! intent(in)
                            , llspan_inf          ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(patchtype), target     :: cpatch     ! Current patch
   integer        , intent(in) :: ico        ! Index of the current cohort
   integer        , intent(in) :: lsl        ! Lowest soil level layer
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Set all cohorts to not have enough leaf area index or wood area index to be        !
   ! resolved.  The code will update this every time step, but we must assign an initial   !
   ! value so the debugger won't complain.                                                 !
   !---------------------------------------------------------------------------------------!
   cpatch%leaf_resolvable(ico) = .false.
   cpatch%wood_resolvable(ico) = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The leaf life span is initialised with the inverse of the turnover rate.  Notice  !
   ! that the turnover rate is in years, but the life span is in months.  Also, some PFTs  !
   ! do not define the turnover rate (temperate cold-deciduous for example), in which case !
   ! we assign a meaningless number just to make sure the variable is initialised.         !
   !---------------------------------------------------------------------------------------!
   if (leaf_turnover_rate(cpatch%pft(ico)) > 0.0) then
      cpatch%llspan(ico) = 12.0 / leaf_turnover_rate(cpatch%pft(ico))
   else
      cpatch%llspan(ico) = llspan_inf
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      The maximum capacity of Rubisco to perform the carboxylase function (Vm) and the !
   ! specific leaf area (SLA) must be assigned with the default values.  These numbers     !
   ! will change only if the PFT uses a light-controlled phenology.                        !
   !---------------------------------------------------------------------------------------!
   select case(phenology(cpatch%pft(ico)))
   case (3)
      cpatch%vm_bar(ico) = vm0_amp / (1.0 + (cpatch%llspan(ico)/vm0_tran)**vm0_slope)      &
                         + vm0_min
   case default
      cpatch%vm_bar(ico) = Vm0(cpatch%pft(ico))
   end select
   cpatch%sla(ico) = sla(cpatch%pft(ico))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Start the fraction of open stomata with 1., since this is the most likely value   !
   ! at night time.  FS_open is initialised with 0., though.                               !
   !---------------------------------------------------------------------------------------!
   cpatch%fsw    (ico) = 1.0
   cpatch%fsn    (ico) = 1.0
   cpatch%fs_open(ico) = 0.0
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Turnover amplitude must be set to 1 because it is used for scaling of the leaf    !
   ! life span in the light-controlled phenology.                                          !
   !---------------------------------------------------------------------------------------!
   cpatch%turnover_amp    (ico) = 1.0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The carbon balance must be initialised with a number other than zero (and better  !
   ! be positive), otherwise a massive mortality will happen at the first year.  The 13th  !
   ! element (current month integration) must be set to 0, though, otherwise the first     !
   ! month carbon balance will be incorrect.                                               !
   !---------------------------------------------------------------------------------------!
   cpatch%cb               (1:12,ico) = 1.0
   cpatch%cb_lightmax      (1:12,ico) = 1.0
   cpatch%cb_moistmax      (1:12,ico) = 1.0
   cpatch%cb_mlmax         (1:12,ico) = 1.0
   cpatch%cbr_bar               (ico) = 1.0
   cpatch%cb                 (13,ico) = 0.0
   cpatch%cb_lightmax        (13,ico) = 0.0
   cpatch%cb_moistmax        (13,ico) = 0.0
   cpatch%cb_mlmax           (13,ico) = 0.0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      The root depth should be the actual level for the roots.                         !
   !---------------------------------------------------------------------------------------!
   cpatch%krdepth(ico) = dbh2krdepth(cpatch%hite(ico),cpatch%dbh(ico),cpatch%pft(ico),lsl)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !       First census must be 1.  Not sure what this variable does, though.  Recruit_dbh !
   ! is a diagnostic variable that tells the cohort status regarding DBH recruitment, and  !
   ! census_status is similar to recruit_dbh, except that is updated only when there is a  !
   ! census.                                                                               !
   !---------------------------------------------------------------------------------------!
   cpatch%first_census  (ico) = 1
   cpatch%recruit_dbh   (ico) = 0
   cpatch%census_status (ico) = 0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Most variables start with zero.                                                  !
   !---------------------------------------------------------------------------------------!
   cpatch%today_leaf_resp       (ico) = 0.0
   cpatch%today_root_resp       (ico) = 0.0
   cpatch%today_gpp             (ico) = 0.0
   cpatch%today_nppleaf         (ico) = 0.0
   cpatch%today_nppfroot        (ico) = 0.0
   cpatch%today_nppsapwood      (ico) = 0.0
   cpatch%today_nppcroot        (ico) = 0.0
   cpatch%today_nppseeds        (ico) = 0.0
   cpatch%today_nppwood         (ico) = 0.0
   cpatch%today_nppdaily        (ico) = 0.0
   cpatch%today_gpp_pot         (ico) = 0.0
   cpatch%today_gpp_lightmax    (ico) = 0.0
   cpatch%today_gpp_moistmax    (ico) = 0.0
   cpatch%today_gpp_mlmax       (ico) = 0.0
   cpatch%light_level           (ico) = 0.0
   cpatch%light_level_beam      (ico) = 0.0
   cpatch%light_level_diff      (ico) = 0.0

   cpatch%par_level_beam        (ico) = 0.0
   cpatch%par_level_diffu       (ico) = 0.0
   cpatch%par_level_diffd       (ico) = 0.0

   cpatch%gpp                   (ico) = 0.0
   cpatch%leaf_respiration      (ico) = 0.0
   cpatch%root_respiration      (ico) = 0.0
   cpatch%leaf_growth_resp      (ico) = 0.0
   cpatch%root_growth_resp      (ico) = 0.0
   cpatch%sapa_growth_resp      (ico) = 0.0
   cpatch%sapb_growth_resp      (ico) = 0.0
   cpatch%leaf_storage_resp     (ico) = 0.0
   cpatch%root_storage_resp     (ico) = 0.0
   cpatch%sapa_storage_resp     (ico) = 0.0
   cpatch%sapb_storage_resp     (ico) = 0.0
   cpatch%monthly_dndt          (ico) = 0.0
   cpatch%monthly_dlnndt        (ico) = 0.0
   cpatch%mort_rate           (:,ico) = 0.0
   cpatch%dagb_dt               (ico) = 0.0
   cpatch%dlnagb_dt             (ico) = 0.0
   cpatch%dba_dt                (ico) = 0.0
   cpatch%dlnba_dt              (ico) = 0.0
   cpatch%ddbh_dt               (ico) = 0.0
   cpatch%dlndbh_dt             (ico) = 0.0
   cpatch%par_l                 (ico) = 0.0
   cpatch%par_l_beam            (ico) = 0.0
   cpatch%par_l_diffuse         (ico) = 0.0
   cpatch%rshort_l              (ico) = 0.0
   cpatch%rshort_l_beam         (ico) = 0.0
   cpatch%rshort_l_diffuse      (ico) = 0.0
   cpatch%rlong_l               (ico) = 0.0
   cpatch%rshort_w              (ico) = 0.0
   cpatch%rshort_w_beam         (ico) = 0.0
   cpatch%rshort_w_diffuse      (ico) = 0.0
   cpatch%rlong_w               (ico) = 0.0
   cpatch%rad_profile         (:,ico) = 0.0
   cpatch%leaf_gbh              (ico) = 0.0
   cpatch%leaf_gbw              (ico) = 0.0
   cpatch%wood_gbh              (ico) = 0.0
   cpatch%wood_gbw              (ico) = 0.0
   cpatch%A_open                (ico) = 0.0
   cpatch%A_closed              (ico) = 0.0
   cpatch%A_light               (ico) = 0.0
   cpatch%A_rubp                (ico) = 0.0
   cpatch%A_co2                 (ico) = 0.0
   cpatch%psi_open              (ico) = 0.0
   cpatch%psi_closed            (ico) = 0.0
   cpatch%water_supply          (ico) = 0.0
   cpatch%gsw_open              (ico) = 0.0
   cpatch%gsw_closed            (ico) = 0.0
   cpatch%leaf_gsw              (ico) = 0.0
   cpatch%leaf_maintenance      (ico) = 0.0
   cpatch%root_maintenance      (ico) = 0.0
   cpatch%leaf_drop             (ico) = 0.0
   cpatch%paw_avg               (ico) = 0.0
   cpatch%elongf                (ico) = 0.0
   cpatch%new_recruit_flag      (ico) = 0
   cpatch%bseeds                (ico) = 0.0
   cpatch%leaf_energy           (ico) = 0.
   cpatch%leaf_hcap             (ico) = 0.
   cpatch%leaf_temp             (ico) = 0.
   cpatch%leaf_temp_pv          (ico) = 0.
   cpatch%leaf_vpdef            (ico) = 0.
   cpatch%leaf_water            (ico) = 0.
   cpatch%leaf_fliq             (ico) = 0.
   cpatch%wood_energy           (ico) = 0.
   cpatch%wood_hcap             (ico) = 0.
   cpatch%wood_temp             (ico) = 0.
   cpatch%wood_temp_pv          (ico) = 0.
   cpatch%wood_water            (ico) = 0.
   cpatch%wood_fliq             (ico) = 0.
   cpatch%veg_wind              (ico) = 0.
   cpatch%lsfc_shv_open         (ico) = 0.
   cpatch%lsfc_shv_closed       (ico) = 0.
   cpatch%lsfc_co2_open         (ico) = 0.
   cpatch%lsfc_co2_closed       (ico) = 0.
   cpatch%lint_shv              (ico) = 0.
   cpatch%lint_co2_open         (ico) = 0.
   cpatch%lint_co2_closed       (ico) = 0.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Fast-averages.                                                                    !
   !---------------------------------------------------------------------------------------!
   cpatch%fmean_gpp               (ico) = 0.0
   cpatch%fmean_npp               (ico) = 0.0
   cpatch%fmean_leaf_resp         (ico) = 0.0
   cpatch%fmean_root_resp         (ico) = 0.0
   cpatch%fmean_leaf_growth_resp  (ico) = 0.0
   cpatch%fmean_root_growth_resp  (ico) = 0.0
   cpatch%fmean_sapa_growth_resp  (ico) = 0.0
   cpatch%fmean_sapb_growth_resp  (ico) = 0.0
   cpatch%fmean_leaf_storage_resp (ico) = 0.0
   cpatch%fmean_root_storage_resp (ico) = 0.0
   cpatch%fmean_sapa_storage_resp (ico) = 0.0
   cpatch%fmean_sapb_storage_resp (ico) = 0.0
   cpatch%fmean_plresp            (ico) = 0.0
   cpatch%fmean_leaf_energy       (ico) = 0.0
   cpatch%fmean_leaf_water        (ico) = 0.0
   cpatch%fmean_leaf_hcap         (ico) = 0.0
   cpatch%fmean_leaf_vpdef        (ico) = 0.0
   cpatch%fmean_leaf_temp         (ico) = 0.0
   cpatch%fmean_leaf_fliq         (ico) = 0.0
   cpatch%fmean_leaf_gsw          (ico) = 0.0
   cpatch%fmean_leaf_gbw          (ico) = 0.0
   cpatch%fmean_wood_energy       (ico) = 0.0
   cpatch%fmean_wood_water        (ico) = 0.0
   cpatch%fmean_wood_hcap         (ico) = 0.0
   cpatch%fmean_wood_temp         (ico) = 0.0
   cpatch%fmean_wood_fliq         (ico) = 0.0
   cpatch%fmean_wood_gbw          (ico) = 0.0
   cpatch%fmean_fs_open           (ico) = 0.0
   cpatch%fmean_fsw               (ico) = 0.0
   cpatch%fmean_fsn               (ico) = 0.0
   cpatch%fmean_a_open            (ico) = 0.0
   cpatch%fmean_a_closed          (ico) = 0.0
   cpatch%fmean_a_net             (ico) = 0.0
   cpatch%fmean_a_light           (ico) = 0.0
   cpatch%fmean_a_rubp            (ico) = 0.0
   cpatch%fmean_a_co2             (ico) = 0.0
   cpatch%fmean_psi_open          (ico) = 0.0
   cpatch%fmean_psi_closed        (ico) = 0.0
   cpatch%fmean_water_supply      (ico) = 0.0
   cpatch%fmean_light_level       (ico) = 0.0
   cpatch%fmean_light_level_beam  (ico) = 0.0
   cpatch%fmean_light_level_diff  (ico) = 0.0

   cpatch%fmean_par_level_beam    (ico) = 0.0
   cpatch%fmean_par_level_diffu   (ico) = 0.0
   cpatch%fmean_par_level_diffd   (ico) = 0.0

   cpatch%fmean_par_l             (ico) = 0.0
   cpatch%fmean_par_l_beam        (ico) = 0.0
   cpatch%fmean_par_l_diff        (ico) = 0.0
   cpatch%fmean_rshort_l          (ico) = 0.0
   cpatch%fmean_rlong_l           (ico) = 0.0
   cpatch%fmean_sensible_lc       (ico) = 0.0
   cpatch%fmean_vapor_lc          (ico) = 0.0
   cpatch%fmean_transp            (ico) = 0.0
   cpatch%fmean_intercepted_al    (ico) = 0.0
   cpatch%fmean_wshed_lg          (ico) = 0.0
   cpatch%fmean_rshort_w          (ico) = 0.0
   cpatch%fmean_rlong_w           (ico) = 0.0
   cpatch%fmean_rad_profile     (:,ico) = 0.0
   cpatch%fmean_sensible_wc       (ico) = 0.0
   cpatch%fmean_vapor_wc          (ico) = 0.0
   cpatch%fmean_intercepted_aw    (ico) = 0.0
   cpatch%fmean_wshed_wg          (ico) = 0.0

   cpatch%fmean_lai               (ico) = 0.0
   cpatch%fmean_bdead             (ico) = 0.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The daily means are allocated only when the user wants the daily or the monthly   !
   ! output, or the mean diurnal cycle.                                                    !
   !---------------------------------------------------------------------------------------!
   if (writing_long) then
      cpatch%dmean_nppleaf           (ico) = 0.0
      cpatch%dmean_nppfroot          (ico) = 0.0
      cpatch%dmean_nppsapwood        (ico) = 0.0
      cpatch%dmean_nppcroot          (ico) = 0.0
      cpatch%dmean_nppseeds          (ico) = 0.0
      cpatch%dmean_nppwood           (ico) = 0.0
      cpatch%dmean_nppdaily          (ico) = 0.0
      cpatch%dmean_gpp               (ico) = 0.0
      cpatch%dmean_npp               (ico) = 0.0
      cpatch%dmean_leaf_resp         (ico) = 0.0
      cpatch%dmean_root_resp         (ico) = 0.0
      cpatch%dmean_leaf_growth_resp  (ico) = 0.0
      cpatch%dmean_root_growth_resp  (ico) = 0.0
      cpatch%dmean_sapa_growth_resp  (ico) = 0.0
      cpatch%dmean_sapb_growth_resp  (ico) = 0.0
      cpatch%dmean_leaf_storage_resp (ico) = 0.0
      cpatch%dmean_root_storage_resp (ico) = 0.0
      cpatch%dmean_sapa_storage_resp (ico) = 0.0
      cpatch%dmean_sapb_storage_resp (ico) = 0.0
      cpatch%dmean_plresp            (ico) = 0.0
      cpatch%dmean_leaf_energy       (ico) = 0.0
      cpatch%dmean_leaf_water        (ico) = 0.0
      cpatch%dmean_leaf_hcap         (ico) = 0.0
      cpatch%dmean_leaf_vpdef        (ico) = 0.0
      cpatch%dmean_leaf_temp         (ico) = 0.0
      cpatch%dmean_leaf_fliq         (ico) = 0.0
      cpatch%dmean_leaf_gsw          (ico) = 0.0
      cpatch%dmean_leaf_gbw          (ico) = 0.0
      cpatch%dmean_wood_energy       (ico) = 0.0
      cpatch%dmean_wood_water        (ico) = 0.0
      cpatch%dmean_wood_hcap         (ico) = 0.0
      cpatch%dmean_wood_temp         (ico) = 0.0
      cpatch%dmean_wood_fliq         (ico) = 0.0
      cpatch%dmean_wood_gbw          (ico) = 0.0
      cpatch%dmean_fs_open           (ico) = 0.0
      cpatch%dmean_fsw               (ico) = 0.0
      cpatch%dmean_fsn               (ico) = 0.0
      cpatch%dmean_a_open            (ico) = 0.0
      cpatch%dmean_a_closed          (ico) = 0.0
      cpatch%dmean_a_net             (ico) = 0.0
      cpatch%dmean_a_light           (ico) = 0.0
      cpatch%dmean_a_rubp            (ico) = 0.0
      cpatch%dmean_a_co2             (ico) = 0.0
      cpatch%dmean_psi_open          (ico) = 0.0
      cpatch%dmean_psi_closed        (ico) = 0.0
      cpatch%dmean_water_supply      (ico) = 0.0
      cpatch%dmean_light_level       (ico) = 0.0
      cpatch%dmean_light_level_beam  (ico) = 0.0
      cpatch%dmean_light_level_diff  (ico) = 0.0

      cpatch%dmean_par_level_beam    (ico) = 0.0
      cpatch%dmean_par_level_diffu   (ico) = 0.0
      cpatch%dmean_par_level_diffd   (ico) = 0.0

      cpatch%dmean_par_l             (ico) = 0.0
      cpatch%dmean_par_l_beam        (ico) = 0.0
      cpatch%dmean_par_l_diff        (ico) = 0.0
      cpatch%dmean_rshort_l          (ico) = 0.0
      cpatch%dmean_rlong_l           (ico) = 0.0
      cpatch%dmean_sensible_lc       (ico) = 0.0
      cpatch%dmean_vapor_lc          (ico) = 0.0
      cpatch%dmean_transp            (ico) = 0.0
      cpatch%dmean_intercepted_al    (ico) = 0.0
      cpatch%dmean_wshed_lg          (ico) = 0.0
      cpatch%dmean_rshort_w          (ico) = 0.0
      cpatch%dmean_rlong_w           (ico) = 0.0
      cpatch%dmean_rad_profile     (:,ico) = 0.0
      cpatch%dmean_sensible_wc       (ico) = 0.0
      cpatch%dmean_vapor_wc          (ico) = 0.0
      cpatch%dmean_intercepted_aw    (ico) = 0.0
      cpatch%dmean_wshed_wg          (ico) = 0.0
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The monthly means are allocated only when the user wants the monthly output or    !
   ! the mean diurnal cycle.                                                               !
   !---------------------------------------------------------------------------------------!
   if (writing_eorq) then
      cpatch%mmean_gpp                 (ico) = 0.0
      cpatch%mmean_npp                 (ico) = 0.0
      cpatch%mmean_leaf_resp           (ico) = 0.0
      cpatch%mmean_root_resp           (ico) = 0.0
      cpatch%mmean_leaf_growth_resp    (ico) = 0.0
      cpatch%mmean_root_growth_resp    (ico) = 0.0
      cpatch%mmean_sapa_growth_resp    (ico) = 0.0
      cpatch%mmean_sapb_growth_resp    (ico) = 0.0
      cpatch%mmean_leaf_storage_resp   (ico) = 0.0
      cpatch%mmean_root_storage_resp   (ico) = 0.0
      cpatch%mmean_sapa_storage_resp   (ico) = 0.0
      cpatch%mmean_sapb_storage_resp   (ico) = 0.0
      cpatch%mmean_plresp              (ico) = 0.0
      cpatch%mmean_leaf_energy         (ico) = 0.0
      cpatch%mmean_leaf_water          (ico) = 0.0
      cpatch%mmean_leaf_hcap           (ico) = 0.0
      cpatch%mmean_leaf_vpdef          (ico) = 0.0
      cpatch%mmean_leaf_temp           (ico) = 0.0
      cpatch%mmean_leaf_fliq           (ico) = 0.0
      cpatch%mmean_leaf_gsw            (ico) = 0.0
      cpatch%mmean_leaf_gbw            (ico) = 0.0
      cpatch%mmean_wood_energy         (ico) = 0.0
      cpatch%mmean_wood_water          (ico) = 0.0
      cpatch%mmean_wood_hcap           (ico) = 0.0
      cpatch%mmean_wood_temp           (ico) = 0.0
      cpatch%mmean_wood_fliq           (ico) = 0.0
      cpatch%mmean_wood_gbw            (ico) = 0.0
      cpatch%mmean_fs_open             (ico) = 0.0
      cpatch%mmean_fsw                 (ico) = 0.0
      cpatch%mmean_fsn                 (ico) = 0.0
      cpatch%mmean_a_open              (ico) = 0.0
      cpatch%mmean_a_closed            (ico) = 0.0
      cpatch%mmean_a_net               (ico) = 0.0
      cpatch%mmean_a_light             (ico) = 0.0
      cpatch%mmean_a_rubp              (ico) = 0.0
      cpatch%mmean_a_co2               (ico) = 0.0
      cpatch%mmean_psi_open            (ico) = 0.0
      cpatch%mmean_psi_closed          (ico) = 0.0
      cpatch%mmean_water_supply        (ico) = 0.0
      cpatch%mmean_light_level         (ico) = 0.0
      cpatch%mmean_light_level_beam    (ico) = 0.0
      cpatch%mmean_light_level_diff    (ico) = 0.0

      cpatch%mmean_par_level_beam      (ico) = 0.0
      cpatch%mmean_par_level_diffu     (ico) = 0.0
      cpatch%mmean_par_level_diffd     (ico) = 0.0

      cpatch%mmean_par_l               (ico) = 0.0
      cpatch%mmean_par_l_beam          (ico) = 0.0
      cpatch%mmean_par_l_diff          (ico) = 0.0
      cpatch%mmean_rshort_l            (ico) = 0.0
      cpatch%mmean_rlong_l             (ico) = 0.0
      cpatch%mmean_sensible_lc         (ico) = 0.0
      cpatch%mmean_vapor_lc            (ico) = 0.0
      cpatch%mmean_transp              (ico) = 0.0
      cpatch%mmean_intercepted_al      (ico) = 0.0
      cpatch%mmean_wshed_lg            (ico) = 0.0
      cpatch%mmean_rshort_w            (ico) = 0.0
      cpatch%mmean_rlong_w             (ico) = 0.0
      cpatch%mmean_rad_profile       (:,ico) = 0.0
      cpatch%mmean_sensible_wc         (ico) = 0.0
      cpatch%mmean_vapor_wc            (ico) = 0.0
      cpatch%mmean_intercepted_aw      (ico) = 0.0
      cpatch%mmean_wshed_wg            (ico) = 0.0
      cpatch%mmean_lai                 (ico) = 0.0
      cpatch%mmean_bleaf               (ico) = 0.0
      cpatch%mmean_broot               (ico) = 0.0
      cpatch%mmean_bstorage            (ico) = 0.0
      cpatch%mmean_mort_rate         (:,ico) = 0.0
      cpatch%mmean_leaf_maintenance    (ico) = 0.0
      cpatch%mmean_root_maintenance    (ico) = 0.0
      cpatch%mmean_leaf_drop           (ico) = 0.0
      cpatch%mmean_cb                  (ico) = 0.0
      cpatch%mmean_nppleaf             (ico) = 0.0
      cpatch%mmean_nppfroot            (ico) = 0.0
      cpatch%mmean_nppsapwood          (ico) = 0.0
      cpatch%mmean_nppcroot            (ico) = 0.0
      cpatch%mmean_nppseeds            (ico) = 0.0
      cpatch%mmean_nppwood             (ico) = 0.0
      cpatch%mmean_nppdaily            (ico) = 0.0
      cpatch%mmsqu_gpp                 (ico) = 0.0
      cpatch%mmsqu_npp                 (ico) = 0.0
      cpatch%mmsqu_plresp              (ico) = 0.0
      cpatch%mmsqu_sensible_lc         (ico) = 0.0
      cpatch%mmsqu_vapor_lc            (ico) = 0.0
      cpatch%mmsqu_transp              (ico) = 0.0
      cpatch%mmsqu_sensible_wc         (ico) = 0.0
      cpatch%mmsqu_vapor_wc            (ico) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The daily means are allocated only when the user wants the daily or the monthly   !
   ! output, or the mean diurnal cycle.                                                    !
   !---------------------------------------------------------------------------------------!
   if (writing_dcyc) then
      cpatch%qmean_gpp               (:,ico) = 0.0
      cpatch%qmean_npp               (:,ico) = 0.0
      cpatch%qmean_leaf_resp         (:,ico) = 0.0
      cpatch%qmean_root_resp         (:,ico) = 0.0
      cpatch%qmean_leaf_growth_resp  (:,ico) = 0.0
      cpatch%qmean_root_growth_resp  (:,ico) = 0.0
      cpatch%qmean_sapa_growth_resp  (:,ico) = 0.0
      cpatch%qmean_sapb_growth_resp  (:,ico) = 0.0
      cpatch%qmean_leaf_storage_resp (:,ico) = 0.0
      cpatch%qmean_root_storage_resp (:,ico) = 0.0
      cpatch%qmean_sapa_storage_resp (:,ico) = 0.0
      cpatch%qmean_sapb_storage_resp (:,ico) = 0.0
      cpatch%qmean_plresp            (:,ico) = 0.0
      cpatch%qmean_leaf_energy       (:,ico) = 0.0
      cpatch%qmean_leaf_water        (:,ico) = 0.0
      cpatch%qmean_leaf_hcap         (:,ico) = 0.0
      cpatch%qmean_leaf_vpdef        (:,ico) = 0.0
      cpatch%qmean_leaf_temp         (:,ico) = 0.0
      cpatch%qmean_leaf_fliq         (:,ico) = 0.0
      cpatch%qmean_leaf_gsw          (:,ico) = 0.0
      cpatch%qmean_leaf_gbw          (:,ico) = 0.0
      cpatch%qmean_wood_energy       (:,ico) = 0.0
      cpatch%qmean_wood_water        (:,ico) = 0.0
      cpatch%qmean_wood_hcap         (:,ico) = 0.0
      cpatch%qmean_wood_temp         (:,ico) = 0.0
      cpatch%qmean_wood_fliq         (:,ico) = 0.0
      cpatch%qmean_wood_gbw          (:,ico) = 0.0
      cpatch%qmean_fs_open           (:,ico) = 0.0
      cpatch%qmean_fsw               (:,ico) = 0.0
      cpatch%qmean_fsn               (:,ico) = 0.0
      cpatch%qmean_a_open            (:,ico) = 0.0
      cpatch%qmean_a_closed          (:,ico) = 0.0
      cpatch%qmean_a_net             (:,ico) = 0.0
      cpatch%qmean_a_light           (:,ico) = 0.0
      cpatch%qmean_a_rubp            (:,ico) = 0.0
      cpatch%qmean_a_co2             (:,ico) = 0.0
      cpatch%qmean_psi_open          (:,ico) = 0.0
      cpatch%qmean_psi_closed        (:,ico) = 0.0
      cpatch%qmean_water_supply      (:,ico) = 0.0
      cpatch%qmean_light_level       (:,ico) = 0.0
      cpatch%qmean_light_level_beam  (:,ico) = 0.0
      cpatch%qmean_light_level_diff  (:,ico) = 0.0

      cpatch%qmean_par_level_beam    (:,ico) = 0.0
      cpatch%qmean_par_level_diffu   (:,ico) = 0.0
      cpatch%qmean_par_level_diffd   (:,ico) = 0.0

      cpatch%qmean_par_l             (:,ico) = 0.0
      cpatch%qmean_par_l_beam        (:,ico) = 0.0
      cpatch%qmean_par_l_diff        (:,ico) = 0.0
      cpatch%qmean_rshort_l          (:,ico) = 0.0
      cpatch%qmean_rlong_l           (:,ico) = 0.0
      cpatch%qmean_sensible_lc       (:,ico) = 0.0
      cpatch%qmean_vapor_lc          (:,ico) = 0.0
      cpatch%qmean_transp            (:,ico) = 0.0
      cpatch%qmean_intercepted_al    (:,ico) = 0.0
      cpatch%qmean_wshed_lg          (:,ico) = 0.0
      cpatch%qmean_rshort_w          (:,ico) = 0.0
      cpatch%qmean_rlong_w           (:,ico) = 0.0
      cpatch%qmean_rad_profile     (:,:,ico) = 0.0
      cpatch%qmean_sensible_wc       (:,ico) = 0.0
      cpatch%qmean_vapor_wc          (:,ico) = 0.0
      cpatch%qmean_intercepted_aw    (:,ico) = 0.0
      cpatch%qmean_wshed_wg          (:,ico) = 0.0
      cpatch%qmsqu_gpp               (:,ico) = 0.0
      cpatch%qmsqu_npp               (:,ico) = 0.0
      cpatch%qmsqu_plresp            (:,ico) = 0.0
      cpatch%qmsqu_sensible_lc       (:,ico) = 0.0
      cpatch%qmsqu_vapor_lc          (:,ico) = 0.0
      cpatch%qmsqu_transp            (:,ico) = 0.0
      cpatch%qmsqu_sensible_wc       (:,ico) = 0.0
      cpatch%qmsqu_vapor_wc          (:,ico) = 0.0
   end if
   !---------------------------------------------------------------------------------------!


   return
end subroutine init_ed_cohort_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine initialise a bunch of patch-level variables.  This should be called !
! whenever a group of new patches is created.                                              !
!------------------------------------------------------------------------------------------!
subroutine init_ed_patch_vars(csite,ipaa,ipaz,lsl)
   use ed_state_vars  , only : sitetype             ! ! structure
   use ed_max_dims    , only : n_pft                ! ! intent(in)
   use grid_coms      , only : nzs                  & ! intent(in)
                             , nzg                  ! ! intent(in)
   use soil_coms      , only : slz                  ! ! intent(in)
   use ed_misc_coms   , only : writing_long         & ! intent(in)
                             , writing_eorq         & ! intent(in)
                             , writing_dcyc         & ! intent(in)
                             , ied_init_mode        ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)   , target     :: csite
   integer          , intent(in) :: ipaa
   integer          , intent(in) :: ipaz
   integer          , intent(in) :: lsl
   !----- Local variables. ----------------------------------------------------------------!
   integer                       :: ncohorts
   integer                       :: ipa
   !---------------------------------------------------------------------------------------!



   !----- Current and previous time steps.  They cannot be set to zero... -----------------!
   csite%htry                      (ipaa:ipaz) = 1.0
   csite%hprev                     (ipaa:ipaz) = 0.1
   !---------------------------------------------------------------------------------------!


   !------ Set water table to the deepest soil layer. -------------------------------------!
   csite%watertable(ipaa:ipaz) = slz(lsl)
   !---------------------------------------------------------------------------------------!


   !------ Initialise hydrology variables. ------------------------------------------------!
   csite%moist_dz           (ipaa:ipaz) = 0.0
   csite%ksat               (ipaa:ipaz) = 0.0
   csite%soil_sat_energy    (ipaa:ipaz) = 0.0
   csite%soil_sat_water     (ipaa:ipaz) = 0.0
   csite%soil_sat_heat      (ipaa:ipaz) = 0.0
   csite%runoff_A         (:,ipaa:ipaz) = 0.0
   csite%runoff_rate        (ipaa:ipaz) = 0.0
   csite%runoff             (ipaa:ipaz) = 0.0
   csite%qrunoff            (ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!


   !------ Initialise soil state variables. -----------------------------------------------!
   if (ied_init_mode /= 7)then
      csite%soil_water              (1:nzg,ipaa:ipaz) = 0.0
   end if
   csite%soil_energy                (1:nzg,ipaa:ipaz) = 0.0
   csite%soil_mstpot                (1:nzg,ipaa:ipaz) = 0.0
   csite%soil_tempk                 (1:nzg,ipaa:ipaz) = 0.0
   csite%soil_fracliq               (1:nzg,ipaa:ipaz) = 0.0
   csite%rootdense                  (1:nzg,ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!



   !------ Initialize sfcwater state variables. -------------------------------------------!
   csite%sfcwater_mass              (1:nzs,ipaa:ipaz) = 0.0
   csite%sfcwater_energy            (1:nzs,ipaa:ipaz) = 0.0
   csite%sfcwater_depth             (1:nzs,ipaa:ipaz) = 0.0
   csite%sfcwater_tempk             (1:nzs,ipaa:ipaz) = 0.0
   csite%sfcwater_fracliq           (1:nzs,ipaa:ipaz) = 0.0
   csite%total_sfcw_depth                 (ipaa:ipaz) = 0.0
   csite%snowfac                          (ipaa:ipaz) = 0.0
   csite%runoff                           (ipaa:ipaz) = 0.0
   csite%qrunoff                          (ipaa:ipaz) = 0.0
   csite%rshort_s                       (:,ipaa:ipaz) = 0.0
   csite%rshort_s_beam                  (:,ipaa:ipaz) = 0.0
   csite%rshort_s_diffuse               (:,ipaa:ipaz) = 0.0
   csite%par_s                          (:,ipaa:ipaz) = 0.0
   csite%par_s_beam                     (:,ipaa:ipaz) = 0.0
   csite%par_s_diffuse                  (:,ipaa:ipaz) = 0.0
   csite%rlong_s                          (ipaa:ipaz) = 0.0
   !------ Number of pounding/snow layers. This is an integer... --------------------------!
   csite%nlev_sfcwater                    (ipaa:ipaz) = 0
   !---------------------------------------------------------------------------------------!


   !------ Initialize peat depth state variables. EJL------------------------------------------!
   csite%litter_depth             (1:nzg,ipaa:ipaz) = 0.0

   !------ Decomposition rates... ---------------------------------------------------------!
   csite%A_decomp                  (1:nzg,ipaa:ipaz) = 0.0
   csite%At_decomp                  (1:nzg,ipaa:ipaz) = 0.0
   csite%Aw_decomp                  (1:nzg,ipaa:ipaz) = 0.0
   csite%f_decomp                  (1:nzg,ipaa:ipaz) = 0.0
   csite%rh                        (1:nzg,ipaa:ipaz) = 0.0
   csite%cwd_rh                          (ipaa:ipaz) = 0.0
   csite%today_A_decomp            (1:nzg,ipaa:ipaz) = 0.0
   csite%today_Aw_decomp            (1:nzg,ipaa:ipaz) = 0.0
   csite%today_At_decomp            (1:nzg,ipaa:ipaz) = 0.0
   csite%today_Af_decomp           (1:nzg,ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!


   !------ Miscellaneous variables. -------------------------------------------------------!
   csite%repro                   (1:n_pft,ipaa:ipaz) = 0.0
   csite%avg_daily_temp                  (ipaa:ipaz) = 0.0
   csite%avg_monthly_gndwater            (ipaa:ipaz) = 0.0
   csite%avg_monthly_waterdef            (ipaa:ipaz) = 0.0
   csite%plant_ag_biomass                (ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!


   !----- Maximum light variables. --------------------------------------------------------!
   csite%A_o_max                 (1:n_pft,ipaa:ipaz) = 0.0
   csite%A_c_max                 (1:n_pft,ipaa:ipaz) = 0.0
   csite%par_l_max                       (ipaa:ipaz) = 0.0
   csite%par_l_beam_max                  (ipaa:ipaz) = 0.0
   csite%par_l_diffuse_max               (ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!

   csite%co2budget_gpp                   (ipaa:ipaz) = 0.0
   csite%co2budget_rh                    (ipaa:ipaz) = 0.0
   csite%co2budget_plresp                (ipaa:ipaz) = 0.0
   csite%co2budget_initialstorage        (ipaa:ipaz) = 0.0
   csite%co2budget_loss2atm              (ipaa:ipaz) = 0.0
   csite%co2budget_denseffect            (ipaa:ipaz) = 0.0
   csite%co2budget_residual              (ipaa:ipaz) = 0.0
   csite%wbudget_precipgain              (ipaa:ipaz) = 0.0
   csite%wbudget_loss2atm                (ipaa:ipaz) = 0.0
   csite%wbudget_loss2runoff             (ipaa:ipaz) = 0.0
   csite%wbudget_loss2drainage           (ipaa:ipaz) = 0.0
   csite%wbudget_denseffect              (ipaa:ipaz) = 0.0
   csite%wbudget_initialstorage          (ipaa:ipaz) = 0.0
   csite%wbudget_residual                (ipaa:ipaz) = 0.0
   csite%ebudget_precipgain              (ipaa:ipaz) = 0.0
   csite%ebudget_netrad                  (ipaa:ipaz) = 0.0
   csite%ebudget_loss2atm                (ipaa:ipaz) = 0.0
   csite%ebudget_loss2runoff             (ipaa:ipaz) = 0.0
   csite%ebudget_loss2drainage           (ipaa:ipaz) = 0.0
   csite%ebudget_denseffect              (ipaa:ipaz) = 0.0
   csite%ebudget_prsseffect              (ipaa:ipaz) = 0.0
   csite%ebudget_initialstorage          (ipaa:ipaz) = 0.0
   csite%ebudget_residual                (ipaa:ipaz) = 0.0
   csite%rshort_g                        (ipaa:ipaz) = 0.0
   csite%rshort_g_beam                   (ipaa:ipaz) = 0.0
   csite%rshort_g_diffuse                (ipaa:ipaz) = 0.0
   csite%par_g                           (ipaa:ipaz) = 0.0
   csite%par_g_beam                      (ipaa:ipaz) = 0.0
   csite%par_g_diffuse                   (ipaa:ipaz) = 0.0
   csite%par_b                           (ipaa:ipaz) = 0.0
   csite%par_b_beam                      (ipaa:ipaz) = 0.0
   csite%par_b_diffuse                   (ipaa:ipaz) = 0.0
   csite%nir_b                           (ipaa:ipaz) = 0.0
   csite%nir_b_beam                      (ipaa:ipaz) = 0.0
   csite%nir_b_diffuse                   (ipaa:ipaz) = 0.0
   csite%rlong_g                         (ipaa:ipaz) = 0.0
   csite%rlong_s                         (ipaa:ipaz) = 0.0
   csite%albedo                          (ipaa:ipaz) = 0.0
   csite%albedo_par                      (ipaa:ipaz) = 0.0
   csite%albedo_nir                      (ipaa:ipaz) = 0.0
   csite%rlong_albedo                    (ipaa:ipaz) = 0.0
   csite%rlongup                         (ipaa:ipaz) = 0.0
   csite%parup                           (ipaa:ipaz) = 0.0
   csite%nirup                           (ipaa:ipaz) = 0.0
   csite%rshortup                        (ipaa:ipaz) = 0.0
   csite%rnet                            (ipaa:ipaz) = 0.0
   csite%fsc_in                          (ipaa:ipaz) = 0.0
   csite%ssc_in                          (ipaa:ipaz) = 0.0
   csite%ssl_in                          (ipaa:ipaz) = 0.0
   csite%fsn_in                          (ipaa:ipaz) = 0.0
   csite%total_plant_nitrogen_uptake     (ipaa:ipaz) = 0.0
   csite%mineralized_N_loss        (1:nzg,ipaa:ipaz) = 0.0
   csite%mineralized_N_input       (1:nzg,ipaa:ipaz) = 0.0
   csite%ustar                           (ipaa:ipaz) = 0.0
   csite%tstar                           (ipaa:ipaz) = 0.0
   csite%qstar                           (ipaa:ipaz) = 0.0
   csite%cstar                           (ipaa:ipaz) = 0.0
   csite%zeta                            (ipaa:ipaz) = 0.0
   csite%ribulk                          (ipaa:ipaz) = 0.0
   csite%upwp                            (ipaa:ipaz) = 0.0
   csite%tpwp                            (ipaa:ipaz) = 0.0
   csite%qpwp                            (ipaa:ipaz) = 0.0
   csite%cpwp                            (ipaa:ipaz) = 0.0
   csite%wpwp                            (ipaa:ipaz) = 0.0
   csite%can_theiv                       (ipaa:ipaz) = 0.0
   csite%can_vpdef                       (ipaa:ipaz) = 0.0
   csite%can_temp                        (ipaa:ipaz) = 0.0
   csite%can_temp_pv                     (ipaa:ipaz) = 0.0
   csite%can_shv                         (ipaa:ipaz) = 0.0
   csite%can_co2                         (ipaa:ipaz) = 0.0
   csite%can_rhos                        (ipaa:ipaz) = 0.0
   csite%can_prss                        (ipaa:ipaz) = 0.0
   csite%can_theta                       (ipaa:ipaz) = 0.0
   csite%can_depth                       (ipaa:ipaz) = 0.0
   csite%opencan_frac                    (ipaa:ipaz) = 0.0
   csite%ground_shv                      (ipaa:ipaz) = 0.0
   csite%ground_ssh                      (ipaa:ipaz) = 0.0
   csite%ground_temp                     (ipaa:ipaz) = 0.0
   csite%ground_fliq                     (ipaa:ipaz) = 0.0
   csite%ggbare                          (ipaa:ipaz) = 0.0
   csite%ggveg                           (ipaa:ipaz) = 0.0
   csite%ggnet                           (ipaa:ipaz) = 0.0
   csite%ggsoil                          (ipaa:ipaz) = 0.0
   csite%rough                           (ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Fast average variables.                                                            !
   !---------------------------------------------------------------------------------------!
   csite%fmean_rh                  (1:nzg,ipaa:ipaz) = 0.0
   csite%fmean_cwd_rh                    (ipaa:ipaz) = 0.0
   csite%fmean_nep                       (ipaa:ipaz) = 0.0
   csite%fmean_rk4step                   (ipaa:ipaz) = 0.0
   csite%fmean_available_water           (ipaa:ipaz) = 0.0
   csite%fmean_can_theiv                 (ipaa:ipaz) = 0.0
   csite%fmean_can_theta                 (ipaa:ipaz) = 0.0
   csite%fmean_can_vpdef                 (ipaa:ipaz) = 0.0
   csite%fmean_can_temp                  (ipaa:ipaz) = 0.0
   csite%fmean_can_shv                   (ipaa:ipaz) = 0.0
   csite%fmean_can_co2                   (ipaa:ipaz) = 0.0
   csite%fmean_can_rhos                  (ipaa:ipaz) = 0.0
   csite%fmean_can_prss                  (ipaa:ipaz) = 0.0
   csite%fmean_gnd_temp                  (ipaa:ipaz) = 0.0
   csite%fmean_gnd_shv                   (ipaa:ipaz) = 0.0
   csite%fmean_can_ggnd                  (ipaa:ipaz) = 0.0
   csite%fmean_sfcw_depth                (ipaa:ipaz) = 0.0
   csite%fmean_sfcw_energy               (ipaa:ipaz) = 0.0
   csite%fmean_sfcw_mass                 (ipaa:ipaz) = 0.0
   csite%fmean_sfcw_temp                 (ipaa:ipaz) = 0.0
   csite%fmean_sfcw_fliq                 (ipaa:ipaz) = 0.0
   csite%fmean_lit_depth                (ipaa:ipaz) = 0.0
   csite%fmean_rshort_gnd                (ipaa:ipaz) = 0.0
   csite%fmean_par_gnd                   (ipaa:ipaz) = 0.0
   csite%fmean_rlong_gnd                 (ipaa:ipaz) = 0.0
   csite%fmean_rlongup                   (ipaa:ipaz) = 0.0
   csite%fmean_parup                     (ipaa:ipaz) = 0.0
   csite%fmean_nirup                     (ipaa:ipaz) = 0.0
   csite%fmean_rshortup                  (ipaa:ipaz) = 0.0
   csite%fmean_rnet                      (ipaa:ipaz) = 0.0
   csite%fmean_albedo                    (ipaa:ipaz) = 0.0
   csite%fmean_albedo_par                (ipaa:ipaz) = 0.0
   csite%fmean_albedo_nir                (ipaa:ipaz) = 0.0
   csite%fmean_rlong_albedo              (ipaa:ipaz) = 0.0
   csite%fmean_ustar                     (ipaa:ipaz) = 0.0
   csite%fmean_tstar                     (ipaa:ipaz) = 0.0
   csite%fmean_qstar                     (ipaa:ipaz) = 0.0
   csite%fmean_cstar                     (ipaa:ipaz) = 0.0
   csite%fmean_carbon_ac                 (ipaa:ipaz) = 0.0
   csite%fmean_carbon_st                 (ipaa:ipaz) = 0.0
   csite%fmean_vapor_gc                  (ipaa:ipaz) = 0.0
   csite%fmean_vapor_ac                  (ipaa:ipaz) = 0.0
   csite%fmean_throughfall               (ipaa:ipaz) = 0.0
   csite%fmean_runoff                    (ipaa:ipaz) = 0.0
   csite%fmean_drainage                  (ipaa:ipaz) = 0.0
   csite%fmean_sensible_gc               (ipaa:ipaz) = 0.0
   csite%fmean_sensible_ac               (ipaa:ipaz) = 0.0
   csite%fmean_qthroughfall              (ipaa:ipaz) = 0.0
   csite%fmean_qrunoff                   (ipaa:ipaz) = 0.0
   csite%fmean_qdrainage                 (ipaa:ipaz) = 0.0
   csite%fmean_soil_energy             (:,ipaa:ipaz) = 0.0
   csite%fmean_soil_mstpot             (:,ipaa:ipaz) = 0.0
   csite%fmean_soil_water              (:,ipaa:ipaz) = 0.0
   csite%fmean_soil_temp               (:,ipaa:ipaz) = 0.0
   csite%fmean_soil_fliq               (:,ipaa:ipaz) = 0.0
   csite%fmean_smoist_gg               (:,ipaa:ipaz) = 0.0
   csite%fmean_transloss               (:,ipaa:ipaz) = 0.0
   csite%fmean_sensible_gg             (:,ipaa:ipaz) = 0.0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Daily means.                                                                       !
   !---------------------------------------------------------------------------------------!
   if (writing_long) then
      csite%dmean_A_decomp             (:,ipaa:ipaz) = 0.0
      csite%dmean_Af_decomp            (:,ipaa:ipaz) = 0.0
      csite%dmean_co2_residual           (ipaa:ipaz) = 0.0
      csite%dmean_energy_residual        (ipaa:ipaz) = 0.0
      csite%dmean_water_residual         (ipaa:ipaz) = 0.0
      csite%dmean_rh                   (:,ipaa:ipaz) = 0.0
      csite%dmean_cwd_rh                 (ipaa:ipaz) = 0.0
      csite%dmean_nep                    (ipaa:ipaz) = 0.0
      csite%dmean_rk4step                (ipaa:ipaz) = 0.0
      csite%dmean_available_water        (ipaa:ipaz) = 0.0
      csite%dmean_can_theiv              (ipaa:ipaz) = 0.0
      csite%dmean_can_theta              (ipaa:ipaz) = 0.0
      csite%dmean_can_vpdef              (ipaa:ipaz) = 0.0
      csite%dmean_can_temp               (ipaa:ipaz) = 0.0
      csite%dmean_can_shv                (ipaa:ipaz) = 0.0
      csite%dmean_can_co2                (ipaa:ipaz) = 0.0
      csite%dmean_can_rhos               (ipaa:ipaz) = 0.0
      csite%dmean_can_prss               (ipaa:ipaz) = 0.0
      csite%dmean_gnd_temp               (ipaa:ipaz) = 0.0
      csite%dmean_gnd_shv                (ipaa:ipaz) = 0.0
      csite%dmean_can_ggnd               (ipaa:ipaz) = 0.0
      csite%dmean_sfcw_depth             (ipaa:ipaz) = 0.0
      csite%dmean_sfcw_energy            (ipaa:ipaz) = 0.0
      csite%dmean_sfcw_mass              (ipaa:ipaz) = 0.0
      csite%dmean_sfcw_temp              (ipaa:ipaz) = 0.0
      csite%dmean_sfcw_fliq              (ipaa:ipaz) = 0.0
      csite%dmean_rshort_gnd             (ipaa:ipaz) = 0.0
      csite%dmean_par_gnd                (ipaa:ipaz) = 0.0
      csite%dmean_rlong_gnd              (ipaa:ipaz) = 0.0
      csite%dmean_rlongup                (ipaa:ipaz) = 0.0
      csite%dmean_parup                  (ipaa:ipaz) = 0.0
      csite%dmean_nirup                  (ipaa:ipaz) = 0.0
      csite%dmean_rshortup               (ipaa:ipaz) = 0.0
      csite%dmean_rnet                   (ipaa:ipaz) = 0.0
      csite%dmean_albedo                 (ipaa:ipaz) = 0.0
      csite%dmean_albedo_par             (ipaa:ipaz) = 0.0
      csite%dmean_albedo_nir             (ipaa:ipaz) = 0.0
      csite%dmean_rlong_albedo           (ipaa:ipaz) = 0.0
      csite%dmean_ustar                  (ipaa:ipaz) = 0.0
      csite%dmean_tstar                  (ipaa:ipaz) = 0.0
      csite%dmean_qstar                  (ipaa:ipaz) = 0.0
      csite%dmean_cstar                  (ipaa:ipaz) = 0.0
      csite%dmean_carbon_ac              (ipaa:ipaz) = 0.0
      csite%dmean_carbon_st              (ipaa:ipaz) = 0.0
      csite%dmean_vapor_gc               (ipaa:ipaz) = 0.0
      csite%dmean_vapor_ac               (ipaa:ipaz) = 0.0
      csite%dmean_throughfall            (ipaa:ipaz) = 0.0
      csite%dmean_runoff                 (ipaa:ipaz) = 0.0
      csite%dmean_drainage               (ipaa:ipaz) = 0.0
      csite%dmean_sensible_gc            (ipaa:ipaz) = 0.0
      csite%dmean_sensible_ac            (ipaa:ipaz) = 0.0
      csite%dmean_qthroughfall           (ipaa:ipaz) = 0.0
      csite%dmean_qrunoff                (ipaa:ipaz) = 0.0
      csite%dmean_qdrainage              (ipaa:ipaz) = 0.0
      csite%dmean_soil_energy          (:,ipaa:ipaz) = 0.0
      csite%dmean_soil_mstpot          (:,ipaa:ipaz) = 0.0
      csite%dmean_soil_water           (:,ipaa:ipaz) = 0.0
      csite%dmean_soil_temp            (:,ipaa:ipaz) = 0.0
      csite%dmean_soil_fliq            (:,ipaa:ipaz) = 0.0
      csite%dmean_smoist_gg            (:,ipaa:ipaz) = 0.0
      csite%dmean_transloss            (:,ipaa:ipaz) = 0.0
      csite%dmean_sensible_gg          (:,ipaa:ipaz) = 0.0
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Monthly means.                                                                     !
   !---------------------------------------------------------------------------------------!
   if (writing_eorq) then
      csite%mmean_rh                   (:,ipaa:ipaz) = 0.0
      csite%mmean_cwd_rh                 (ipaa:ipaz) = 0.0
      csite%mmean_nep                    (ipaa:ipaz) = 0.0
      csite%mmean_rk4step                (ipaa:ipaz) = 0.0
      csite%mmean_available_water        (ipaa:ipaz) = 0.0
      csite%mmean_can_theiv              (ipaa:ipaz) = 0.0
      csite%mmean_can_theta              (ipaa:ipaz) = 0.0
      csite%mmean_can_vpdef              (ipaa:ipaz) = 0.0
      csite%mmean_can_temp               (ipaa:ipaz) = 0.0
      csite%mmean_can_shv                (ipaa:ipaz) = 0.0
      csite%mmean_can_co2                (ipaa:ipaz) = 0.0
      csite%mmean_can_rhos               (ipaa:ipaz) = 0.0
      csite%mmean_can_prss               (ipaa:ipaz) = 0.0
      csite%mmean_gnd_temp               (ipaa:ipaz) = 0.0
      csite%mmean_gnd_shv                (ipaa:ipaz) = 0.0
      csite%mmean_can_ggnd               (ipaa:ipaz) = 0.0
      csite%mmean_sfcw_depth             (ipaa:ipaz) = 0.0
      csite%mmean_sfcw_energy            (ipaa:ipaz) = 0.0
      csite%mmean_sfcw_mass              (ipaa:ipaz) = 0.0
      csite%mmean_sfcw_temp              (ipaa:ipaz) = 0.0
      csite%mmean_sfcw_fliq              (ipaa:ipaz) = 0.0
      csite%mmean_rshort_gnd             (ipaa:ipaz) = 0.0
      csite%mmean_par_gnd                (ipaa:ipaz) = 0.0
      csite%mmean_rlong_gnd              (ipaa:ipaz) = 0.0
      csite%mmean_rlongup                (ipaa:ipaz) = 0.0
      csite%mmean_parup                  (ipaa:ipaz) = 0.0
      csite%mmean_nirup                  (ipaa:ipaz) = 0.0
      csite%mmean_rshortup               (ipaa:ipaz) = 0.0
      csite%mmean_rnet                   (ipaa:ipaz) = 0.0
      csite%mmean_albedo                 (ipaa:ipaz) = 0.0
      csite%mmean_albedo_par             (ipaa:ipaz) = 0.0
      csite%mmean_albedo_nir             (ipaa:ipaz) = 0.0
      csite%mmean_rlong_albedo           (ipaa:ipaz) = 0.0
      csite%mmean_ustar                  (ipaa:ipaz) = 0.0
      csite%mmean_tstar                  (ipaa:ipaz) = 0.0
      csite%mmean_qstar                  (ipaa:ipaz) = 0.0
      csite%mmean_cstar                  (ipaa:ipaz) = 0.0
      csite%mmean_carbon_ac              (ipaa:ipaz) = 0.0
      csite%mmean_carbon_st              (ipaa:ipaz) = 0.0
      csite%mmean_vapor_gc               (ipaa:ipaz) = 0.0
      csite%mmean_vapor_ac               (ipaa:ipaz) = 0.0
      csite%mmean_throughfall            (ipaa:ipaz) = 0.0
      csite%mmean_runoff                 (ipaa:ipaz) = 0.0
      csite%mmean_drainage               (ipaa:ipaz) = 0.0
      csite%mmean_sensible_gc            (ipaa:ipaz) = 0.0
      csite%mmean_sensible_ac            (ipaa:ipaz) = 0.0
      csite%mmean_qthroughfall           (ipaa:ipaz) = 0.0
      csite%mmean_qrunoff                (ipaa:ipaz) = 0.0
      csite%mmean_qdrainage              (ipaa:ipaz) = 0.0
      csite%mmean_fast_soil_c            (ipaa:ipaz) = 0.0
      csite%mmean_slow_soil_c            (ipaa:ipaz) = 0.0
      csite%mmean_struct_soil_c          (ipaa:ipaz) = 0.0
      csite%mmean_struct_soil_l          (ipaa:ipaz) = 0.0
      csite%mmean_fast_soil_n            (ipaa:ipaz) = 0.0
      csite%mmean_mineral_soil_n         (ipaa:ipaz) = 0.0
      csite%mmean_A_decomp             (:,ipaa:ipaz) = 0.0
      csite%mmean_Af_decomp            (:,ipaa:ipaz) = 0.0
      csite%mmean_co2_residual           (ipaa:ipaz) = 0.0
      csite%mmean_energy_residual        (ipaa:ipaz) = 0.0
      csite%mmean_water_residual         (ipaa:ipaz) = 0.0
      csite%mmean_soil_energy          (:,ipaa:ipaz) = 0.0
      csite%mmean_soil_mstpot          (:,ipaa:ipaz) = 0.0
      csite%mmean_soil_water           (:,ipaa:ipaz) = 0.0
      csite%mmean_soil_temp            (:,ipaa:ipaz) = 0.0
      csite%mmean_soil_fliq            (:,ipaa:ipaz) = 0.0
      csite%mmean_smoist_gg            (:,ipaa:ipaz) = 0.0
      csite%mmean_transloss            (:,ipaa:ipaz) = 0.0
      csite%mmean_sensible_gg          (:,ipaa:ipaz) = 0.0
      csite%mmsqu_rh                     (ipaa:ipaz) = 0.0
      csite%mmsqu_cwd_rh                 (ipaa:ipaz) = 0.0
      csite%mmsqu_nep                    (ipaa:ipaz) = 0.0
      csite%mmsqu_rlongup                (ipaa:ipaz) = 0.0
      csite%mmsqu_parup                  (ipaa:ipaz) = 0.0
      csite%mmsqu_nirup                  (ipaa:ipaz) = 0.0
      csite%mmsqu_rshortup               (ipaa:ipaz) = 0.0
      csite%mmsqu_rnet                   (ipaa:ipaz) = 0.0
      csite%mmsqu_albedo                 (ipaa:ipaz) = 0.0
      csite%mmsqu_ustar                  (ipaa:ipaz) = 0.0
      csite%mmsqu_carbon_ac              (ipaa:ipaz) = 0.0
      csite%mmsqu_carbon_st              (ipaa:ipaz) = 0.0
      csite%mmsqu_vapor_gc               (ipaa:ipaz) = 0.0
      csite%mmsqu_vapor_ac               (ipaa:ipaz) = 0.0
      csite%mmsqu_sensible_gc            (ipaa:ipaz) = 0.0
      csite%mmsqu_sensible_ac            (ipaa:ipaz) = 0.0
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Mean diel.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (writing_dcyc) then
      csite%qmean_rh                   (:,ipaa:ipaz) = 0.0
      csite%qmean_cwd_rh               (:,ipaa:ipaz) = 0.0
      csite%qmean_nep                  (:,ipaa:ipaz) = 0.0
      csite%qmean_rk4step              (:,ipaa:ipaz) = 0.0
      csite%qmean_available_water      (:,ipaa:ipaz) = 0.0
      csite%qmean_can_theiv            (:,ipaa:ipaz) = 0.0
      csite%qmean_can_theta            (:,ipaa:ipaz) = 0.0
      csite%qmean_can_vpdef            (:,ipaa:ipaz) = 0.0
      csite%qmean_can_temp             (:,ipaa:ipaz) = 0.0
      csite%qmean_can_shv              (:,ipaa:ipaz) = 0.0
      csite%qmean_can_co2              (:,ipaa:ipaz) = 0.0
      csite%qmean_can_rhos             (:,ipaa:ipaz) = 0.0
      csite%qmean_can_prss             (:,ipaa:ipaz) = 0.0
      csite%qmean_gnd_temp             (:,ipaa:ipaz) = 0.0
      csite%qmean_gnd_shv              (:,ipaa:ipaz) = 0.0
      csite%qmean_can_ggnd             (:,ipaa:ipaz) = 0.0
      csite%qmean_sfcw_depth           (:,ipaa:ipaz) = 0.0
      csite%qmean_sfcw_energy          (:,ipaa:ipaz) = 0.0
      csite%qmean_sfcw_mass            (:,ipaa:ipaz) = 0.0
      csite%qmean_sfcw_temp            (:,ipaa:ipaz) = 0.0
      csite%qmean_sfcw_fliq            (:,ipaa:ipaz) = 0.0
      csite%qmean_rshort_gnd           (:,ipaa:ipaz) = 0.0
      csite%qmean_par_gnd              (:,ipaa:ipaz) = 0.0
      csite%qmean_rlong_gnd            (:,ipaa:ipaz) = 0.0
      csite%qmean_rlongup              (:,ipaa:ipaz) = 0.0
      csite%qmean_parup                (:,ipaa:ipaz) = 0.0
      csite%qmean_nirup                (:,ipaa:ipaz) = 0.0
      csite%qmean_rshortup             (:,ipaa:ipaz) = 0.0
      csite%qmean_rnet                 (:,ipaa:ipaz) = 0.0
      csite%qmean_albedo               (:,ipaa:ipaz) = 0.0
      csite%qmean_albedo_par           (:,ipaa:ipaz) = 0.0
      csite%qmean_albedo_nir           (:,ipaa:ipaz) = 0.0
      csite%qmean_rlong_albedo         (:,ipaa:ipaz) = 0.0
      csite%qmean_ustar                (:,ipaa:ipaz) = 0.0
      csite%qmean_tstar                (:,ipaa:ipaz) = 0.0
      csite%qmean_qstar                (:,ipaa:ipaz) = 0.0
      csite%qmean_cstar                (:,ipaa:ipaz) = 0.0
      csite%qmean_carbon_ac            (:,ipaa:ipaz) = 0.0
      csite%qmean_carbon_st            (:,ipaa:ipaz) = 0.0
      csite%qmean_vapor_gc             (:,ipaa:ipaz) = 0.0
      csite%qmean_vapor_ac             (:,ipaa:ipaz) = 0.0
      csite%qmean_throughfall          (:,ipaa:ipaz) = 0.0
      csite%qmean_runoff               (:,ipaa:ipaz) = 0.0
      csite%qmean_drainage             (:,ipaa:ipaz) = 0.0
      csite%qmean_sensible_gc          (:,ipaa:ipaz) = 0.0
      csite%qmean_sensible_ac          (:,ipaa:ipaz) = 0.0
      csite%qmean_qthroughfall         (:,ipaa:ipaz) = 0.0
      csite%qmean_qrunoff              (:,ipaa:ipaz) = 0.0
      csite%qmean_qdrainage            (:,ipaa:ipaz) = 0.0
      csite%qmsqu_rh                   (:,ipaa:ipaz) = 0.0
      csite%qmsqu_cwd_rh               (:,ipaa:ipaz) = 0.0
      csite%qmsqu_nep                  (:,ipaa:ipaz) = 0.0
      csite%qmsqu_rlongup              (:,ipaa:ipaz) = 0.0
      csite%qmsqu_parup                (:,ipaa:ipaz) = 0.0
      csite%qmsqu_nirup                (:,ipaa:ipaz) = 0.0
      csite%qmsqu_rshortup             (:,ipaa:ipaz) = 0.0
      csite%qmsqu_rnet                 (:,ipaa:ipaz) = 0.0
      csite%qmsqu_albedo               (:,ipaa:ipaz) = 0.0
      csite%qmsqu_ustar                (:,ipaa:ipaz) = 0.0
      csite%qmsqu_carbon_ac            (:,ipaa:ipaz) = 0.0
      csite%qmsqu_carbon_st            (:,ipaa:ipaz) = 0.0
      csite%qmsqu_vapor_gc             (:,ipaa:ipaz) = 0.0
      csite%qmsqu_vapor_ac             (:,ipaa:ipaz) = 0.0
      csite%qmsqu_sensible_gc          (:,ipaa:ipaz) = 0.0
      csite%qmsqu_sensible_ac          (:,ipaa:ipaz) = 0.0
      csite%qmean_soil_energy        (:,:,ipaa:ipaz) = 0.0
      csite%qmean_soil_mstpot        (:,:,ipaa:ipaz) = 0.0
      csite%qmean_soil_water         (:,:,ipaa:ipaz) = 0.0
      csite%qmean_soil_temp          (:,:,ipaa:ipaz) = 0.0
      csite%qmean_soil_fliq          (:,:,ipaa:ipaz) = 0.0
      csite%qmean_smoist_gg          (:,:,ipaa:ipaz) = 0.0
      csite%qmean_transloss          (:,:,ipaa:ipaz) = 0.0
      csite%qmean_sensible_gg        (:,:,ipaa:ipaz) = 0.0
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Update the number of cohorts in this site.                                        !
   !---------------------------------------------------------------------------------------!
   ncohorts = 0
   do ipa=1,csite%npatches
      ncohorts = ncohorts + csite%patch(ipa)%ncohorts
   end do
   csite%cohort_count = ncohorts
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ed_patch_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine initialises some site-level variables.                              !
!------------------------------------------------------------------------------------------!
subroutine init_ed_site_vars(cpoly)
   use ed_state_vars, only : polygontype      ! ! intent(in)
   use ed_max_dims  , only : n_pft            & ! intent(in)
                           , n_dbh            ! ! intent(in)
   use pft_coms     , only : agri_stock       & ! intent(in)
                           , plantation_stock ! ! intent(in)
   use ed_misc_coms , only : writing_long     & ! intent(in)
                           , writing_eorq     & ! intent(in)
                           , writing_dcyc     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(polygontype), target     :: cpoly
   !----- External functions. -------------------------------------------------------------!
   integer          , external   :: julday
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Size and PFT tables.                                                             !
   !---------------------------------------------------------------------------------------!
   cpoly%basal_area       (1:n_pft, 1:n_dbh, :) = 0.0
   cpoly%agb              (1:n_pft, 1:n_dbh, :) = 0.0

   cpoly%basal_area_growth(1:n_pft, 1:n_dbh, :) = 0.0
   cpoly%basal_area_mort  (1:n_pft, 1:n_dbh, :) = 0.0
   cpoly%basal_area_cut   (1:n_pft, 1:n_dbh, :) = 0.0

   cpoly%agb_growth       (1:n_pft, 1:n_dbh, :) = 0.0
   cpoly%agb_mort         (1:n_pft, 1:n_dbh, :) = 0.0
   cpoly%agb_cut          (1:n_pft, 1:n_dbh, :) = 0.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Radiation-related variables.                                                     !
   !---------------------------------------------------------------------------------------!
   cpoly%cosaoi    (:) = 0.0
   cpoly%daylight  (:) = 0.0
   cpoly%nighttime (:) = .true.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Hydrology-related variables.                                                     !
   !---------------------------------------------------------------------------------------!
   cpoly%runoff   (:) = 0.0
   cpoly%qrunoff  (:) = 0.0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Phenology variables (green_leaf_factor should be merged with elongf, as they     !
   ! represent the same thing).                                                            !
   !---------------------------------------------------------------------------------------!
   cpoly%green_leaf_factor(1:n_pft,:) = 1.0
   cpoly%leaf_aging_factor(1:n_pft,:) = 1.0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialise the minimum monthly temperature with a very large value, this is      !
   ! going to be reduced as the canopy temperature is updated.                             !
   !---------------------------------------------------------------------------------------!
   cpoly%min_monthly_temp(:) = huge(1.)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialise monthly rainfall with some arbitrary but high number.  This will      !
   ! probably prevent fires to happen at the first year, but all data will be replaced by  !
   ! actual rainfall after 12 months.  In the future we may initialise with climatological !
   ! rainfall.                                                                             !
   !---------------------------------------------------------------------------------------!
   cpoly%avg_monthly_pcpg(:,:) = 500.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialise several disturbance- and LU-related variables.                        !
   !---------------------------------------------------------------------------------------!
   cpoly%plantation                       (:) = 0
   cpoly%agri_stocking_pft                (:) = agri_stock
   cpoly%agri_stocking_density            (:) = 10.0
   cpoly%plantation_stocking_pft          (:) = plantation_stock
   cpoly%plantation_stocking_density      (:) = 4.0
   cpoly%primary_harvest_target           (:) = 0.0
   cpoly%secondary_harvest_target         (:) = 0.0
   cpoly%primary_harvest_memory           (:) = 0.0
   cpoly%secondary_harvest_memory         (:) = 0.0
   cpoly%ignition_rate                    (:) = 0.0
   cpoly%lambda_fire                    (:,:) = 0.0
   cpoly%disturbance_memory           (:,:,:) = 0.0
   cpoly%disturbance_rates            (:,:,:) = 0.0
   !---------------------------------------------------------------------------------------!



   !----- Initialise the mean radiation. --------------------------------------------------!
   cpoly%rad_avg                          (:) = 200.0
   !---------------------------------------------------------------------------------------!


   !----- Fast means. ---------------------------------------------------------------------!
   cpoly%fmean_atm_theiv                (:) = 0.0
   cpoly%fmean_atm_theta                (:) = 0.0
   cpoly%fmean_atm_temp                 (:) = 0.0
   cpoly%fmean_atm_vpdef                (:) = 0.0
   cpoly%fmean_atm_shv                  (:) = 0.0
   cpoly%fmean_atm_rshort               (:) = 0.0
   cpoly%fmean_atm_rshort_diff          (:) = 0.0
   cpoly%fmean_atm_par                  (:) = 0.0
   cpoly%fmean_atm_par_diff             (:) = 0.0
   cpoly%fmean_atm_rlong                (:) = 0.0
   cpoly%fmean_atm_vels                 (:) = 0.0
   cpoly%fmean_atm_rhos                 (:) = 0.0
   cpoly%fmean_atm_prss                 (:) = 0.0
   cpoly%fmean_atm_co2                  (:) = 0.0
   cpoly%fmean_pcpg                     (:) = 0.0
   cpoly%fmean_qpcpg                    (:) = 0.0
   cpoly%fmean_dpcpg                    (:) = 0.0
   !---------------------------------------------------------------------------------------!



   !----- Daily means. --------------------------------------------------------------------!
   if (writing_long) then
      cpoly%dmean_atm_theiv             (:) = 0.0
      cpoly%dmean_atm_theta             (:) = 0.0
      cpoly%dmean_atm_temp              (:) = 0.0
      cpoly%dmean_atm_vpdef             (:) = 0.0
      cpoly%dmean_atm_shv               (:) = 0.0
      cpoly%dmean_atm_rshort            (:) = 0.0
      cpoly%dmean_atm_rshort_diff       (:) = 0.0
      cpoly%dmean_atm_par               (:) = 0.0
      cpoly%dmean_atm_par_diff          (:) = 0.0
      cpoly%dmean_atm_rlong             (:) = 0.0
      cpoly%dmean_atm_vels              (:) = 0.0
      cpoly%dmean_atm_rhos              (:) = 0.0
      cpoly%dmean_atm_prss              (:) = 0.0
      cpoly%dmean_atm_co2               (:) = 0.0
      cpoly%dmean_pcpg                  (:) = 0.0
      cpoly%dmean_qpcpg                 (:) = 0.0
      cpoly%dmean_dpcpg                 (:) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !----- Monthly means. ------------------------------------------------------------------!
   if (writing_eorq) then
      cpoly%mmean_atm_theiv             (:) = 0.0
      cpoly%mmean_atm_theta             (:) = 0.0
      cpoly%mmean_atm_temp              (:) = 0.0
      cpoly%mmean_atm_vpdef             (:) = 0.0
      cpoly%mmean_atm_shv               (:) = 0.0
      cpoly%mmean_atm_rshort            (:) = 0.0
      cpoly%mmean_atm_rshort_diff       (:) = 0.0
      cpoly%mmean_atm_par               (:) = 0.0
      cpoly%mmean_atm_par_diff          (:) = 0.0
      cpoly%mmean_atm_rlong             (:) = 0.0
      cpoly%mmean_atm_vels              (:) = 0.0
      cpoly%mmean_atm_rhos              (:) = 0.0
      cpoly%mmean_atm_prss              (:) = 0.0
      cpoly%mmean_atm_co2               (:) = 0.0
      cpoly%mmean_pcpg                  (:) = 0.0
      cpoly%mmean_qpcpg                 (:) = 0.0
      cpoly%mmean_dpcpg                 (:) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !----- Monthly means. ------------------------------------------------------------------!
   if (writing_dcyc) then
      cpoly%qmean_atm_theiv           (:,:) = 0.0
      cpoly%qmean_atm_theta           (:,:) = 0.0
      cpoly%qmean_atm_temp            (:,:) = 0.0
      cpoly%qmean_atm_vpdef           (:,:) = 0.0
      cpoly%qmean_atm_shv             (:,:) = 0.0
      cpoly%qmean_atm_rshort          (:,:) = 0.0
      cpoly%qmean_atm_rshort_diff     (:,:) = 0.0
      cpoly%qmean_atm_par             (:,:) = 0.0
      cpoly%qmean_atm_par_diff        (:,:) = 0.0
      cpoly%qmean_atm_rlong           (:,:) = 0.0
      cpoly%qmean_atm_vels            (:,:) = 0.0
      cpoly%qmean_atm_rhos            (:,:) = 0.0
      cpoly%qmean_atm_prss            (:,:) = 0.0
      cpoly%qmean_atm_co2             (:,:) = 0.0
      cpoly%qmean_pcpg                (:,:) = 0.0
      cpoly%qmean_qpcpg               (:,:) = 0.0
      cpoly%qmean_dpcpg               (:,:) = 0.0
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ed_site_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_ed_poly_vars(cgrid)
   use structural_growth_module
   use ed_state_vars, only : edtype       & ! structure
                           , polygontype  & ! structure
                           , sitetype     ! ! structure
   use ed_misc_coms , only : writing_long & ! intent(in)
                           , writing_eorq & ! intent(in)
                           , writing_dcyc ! ! intent(in)
   use consts_coms  , only : day_sec      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   integer                    :: ipy
   integer                    :: isi
   real                       :: soil_C
   real                       :: soil_N
   real                       :: veg_C
   real                       :: veg_N
   real                       :: patchload
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   !
   !---------------------------------------------------------------------------------------!
   !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.  This  !
   ! works in off-line runs, but it causes memory leaks (and crashes) in the coupled runs  !
   ! over the ocean, where cgrid%npolygons can be 0 if one of the subdomains falls         !
   ! entirely over the ocean.  Thanks!                                                     !
   !---------------------------------------------------------------------------------------!
   ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
   !                 !      safe way of initialising the variable.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Define a nominal initial value of patch workload.  Normally we start with the    !
   ! RK4 time step to be 1 second, so each patch will contribute with 86400 time steps per !
   ! day.                                                                                  !
   !---------------------------------------------------------------------------------------!
   patchload = day_sec

   do ipy = 1,cgrid%npolygons

      !------------------------------------------------------------------------------------!
      !     This is the right and safe place to initialise polygon-level (cgrid) vari-     !
      ! ables, so in case npolygons is zero this will not cause memory leaks.  I know,     !
      ! this never happens in off-line runs, but it is quite common in coupled runs...     !
      ! Whenever one of the nodes receives a sub-domain where all the points are over the  !
      ! ocean, ED will not assign any polygon in that sub-domain, which means that that    !
      ! node will have 0 polygons, and the variables cannot be allocated.  If you try to   !
      ! access the polygon level variable outside the loop, then the model crashes due to  !
      ! segmentation violation (a bad thing), whereas by putting the variables here both   !
      ! the off-line model and the coupled runs will work, because this loop will be       !
      ! skipped when there is no polygon.                                                  !
      !------------------------------------------------------------------------------------!
      call compute_C_and_N_storage(cgrid,ipy,soil_C, soil_N, veg_C, veg_N)
      cgrid%cbudget_initialstorage(ipy) = soil_C + veg_C
      cgrid%nbudget_initialstorage(ipy) = soil_N + veg_N
      cgrid%cbudget_nep           (ipy) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      cgrid%avg_lai_ebalvars        (:,:,ipy) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Hydrology stuff.                                                             !
      !------------------------------------------------------------------------------------!
      !cgrid%wbar     (ipy) = 0.0
      !cgrid%Te       (ipy) = 0.0
      cgrid%zbar     (ipy) = 0.0
      cgrid%sheat    (ipy) = 0.0
      cgrid%baseflow (ipy) = 0.0
      cgrid%runoff   (ipy) = 0.0
      cgrid%qrunoff  (ipy) = 0.0
      cgrid%swliq    (ipy) = 0.0
      !------------------------------------------------------------------------------------!


      !----- Set all biomass and soil pools to zero. --------------------------------------!
      cgrid%total_agb                   (ipy) = 0.0
      cgrid%total_basal_area            (ipy) = 0.0
      cgrid%total_agb_growth            (ipy) = 0.0
      cgrid%total_agb_mort              (ipy) = 0.0
      cgrid%total_agb_recruit           (ipy) = 0.0
      cgrid%total_basal_area_growth     (ipy) = 0.0
      cgrid%total_basal_area_mort       (ipy) = 0.0
      cgrid%total_basal_area_recruit    (ipy) = 0.0
      cgrid%nplant                  (:,:,ipy) = 0.0
      cgrid%agb                     (:,:,ipy) = 0.0
      cgrid%lai                     (:,:,ipy) = 0.0
      cgrid%wai                     (:,:,ipy) = 0.0
      cgrid%basal_area              (:,:,ipy) = 0.0
      cgrid%bdead                   (:,:,ipy) = 0.0
      cgrid%balive                  (:,:,ipy) = 0.0
      cgrid%bleaf                   (:,:,ipy) = 0.0
      cgrid%broot                   (:,:,ipy) = 0.0
      cgrid%bsapwooda               (:,:,ipy) = 0.0
      cgrid%bsapwoodb               (:,:,ipy) = 0.0
      cgrid%bseeds                  (:,:,ipy) = 0.0
      cgrid%bstorage                (:,:,ipy) = 0.0
      cgrid%bdead_n                 (:,:,ipy) = 0.0
      cgrid%balive_n                (:,:,ipy) = 0.0
      cgrid%bleaf_n                 (:,:,ipy) = 0.0
      cgrid%broot_n                 (:,:,ipy) = 0.0
      cgrid%bsapwooda_n             (:,:,ipy) = 0.0
      cgrid%bsapwoodb_n             (:,:,ipy) = 0.0
      cgrid%bseeds_n                (:,:,ipy) = 0.0
      cgrid%bstorage_n              (:,:,ipy) = 0.0
      cgrid%leaf_maintenance        (:,:,ipy) = 0.0
      cgrid%root_maintenance        (:,:,ipy) = 0.0
      cgrid%leaf_drop               (:,:,ipy) = 0.0
      cgrid%fast_soil_c                 (ipy) = 0.0
      cgrid%slow_soil_c                 (ipy) = 0.0
      cgrid%struct_soil_c               (ipy) = 0.0
      cgrid%cwd_c                       (ipy) = 0.0
      cgrid%fast_soil_n                 (ipy) = 0.0
      cgrid%mineral_soil_n              (ipy) = 0.0
      cgrid%cwd_n                       (ipy) = 0.0
      !------------------------------------------------------------------------------------!




      !----- Count how many patches we have, and add to the workload. ---------------------!
      cgrid%workload(:,ipy)  = 0.0
      cpoly => cgrid%polygon(ipy)
      do isi = 1, cpoly%nsites
         csite => cpoly%site(isi)
         cgrid%workload(1:12,ipy) = cgrid%workload(1:12,ipy)                               &
                                  + real(csite%npatches) * patchload
      end do
      !------------------------------------------------------------------------------------!


      !----- Fast averages. ---------------------------------------------------------------!
      cgrid%fmean_gpp                  (ipy) = 0.0
      cgrid%fmean_npp                  (ipy) = 0.0
      cgrid%fmean_leaf_resp            (ipy) = 0.0
      cgrid%fmean_root_resp            (ipy) = 0.0
      cgrid%fmean_leaf_growth_resp     (ipy) = 0.0
      cgrid%fmean_root_growth_resp     (ipy) = 0.0
      cgrid%fmean_sapa_growth_resp     (ipy) = 0.0
      cgrid%fmean_sapb_growth_resp     (ipy) = 0.0
      cgrid%fmean_leaf_storage_resp    (ipy) = 0.0
      cgrid%fmean_root_storage_resp    (ipy) = 0.0
      cgrid%fmean_sapa_storage_resp    (ipy) = 0.0
      cgrid%fmean_sapb_storage_resp    (ipy) = 0.0
      cgrid%fmean_plresp               (ipy) = 0.0
      cgrid%fmean_leaf_energy          (ipy) = 0.0
      cgrid%fmean_leaf_water           (ipy) = 0.0
      cgrid%fmean_leaf_hcap            (ipy) = 0.0
      cgrid%fmean_leaf_vpdef           (ipy) = 0.0
      cgrid%fmean_leaf_temp            (ipy) = 0.0
      cgrid%fmean_leaf_fliq            (ipy) = 0.0
      cgrid%fmean_leaf_gsw             (ipy) = 0.0
      cgrid%fmean_leaf_gbw             (ipy) = 0.0
      cgrid%fmean_wood_energy          (ipy) = 0.0
      cgrid%fmean_wood_water           (ipy) = 0.0
      cgrid%fmean_wood_hcap            (ipy) = 0.0
      cgrid%fmean_wood_temp            (ipy) = 0.0
      cgrid%fmean_wood_fliq            (ipy) = 0.0
      cgrid%fmean_wood_gbw             (ipy) = 0.0
      cgrid%fmean_fs_open              (ipy) = 0.0
      cgrid%fmean_fsw                  (ipy) = 0.0
      cgrid%fmean_fsn                  (ipy) = 0.0
      cgrid%fmean_a_open               (ipy) = 0.0
      cgrid%fmean_a_closed             (ipy) = 0.0
      cgrid%fmean_a_net                (ipy) = 0.0
      cgrid%fmean_a_light              (ipy) = 0.0
      cgrid%fmean_a_rubp               (ipy) = 0.0
      cgrid%fmean_a_co2                (ipy) = 0.0
      cgrid%fmean_psi_open             (ipy) = 0.0
      cgrid%fmean_psi_closed           (ipy) = 0.0
      cgrid%fmean_water_supply         (ipy) = 0.0
      cgrid%fmean_par_l                (ipy) = 0.0
      cgrid%fmean_par_l_beam           (ipy) = 0.0
      cgrid%fmean_par_l_diff           (ipy) = 0.0
      cgrid%fmean_rshort_l             (ipy) = 0.0
      cgrid%fmean_rlong_l              (ipy) = 0.0
      cgrid%fmean_sensible_lc          (ipy) = 0.0
      cgrid%fmean_vapor_lc             (ipy) = 0.0
      cgrid%fmean_transp               (ipy) = 0.0
      cgrid%fmean_intercepted_al       (ipy) = 0.0
      cgrid%fmean_wshed_lg             (ipy) = 0.0
      cgrid%fmean_rshort_w             (ipy) = 0.0
      cgrid%fmean_rlong_w              (ipy) = 0.0
      cgrid%fmean_sensible_wc          (ipy) = 0.0
      cgrid%fmean_vapor_wc             (ipy) = 0.0
      cgrid%fmean_intercepted_aw       (ipy) = 0.0
      cgrid%fmean_wshed_wg             (ipy) = 0.0

      cgrid%fmean_lai                  (ipy) = 0.0
      cgrid%fmean_bdead                (ipy) = 0.0

      cgrid%fmean_rh                 (:,ipy) = 0.0
      cgrid%fmean_cwd_rh               (ipy) = 0.0
      cgrid%fmean_nep                  (ipy) = 0.0
      cgrid%fmean_rk4step              (ipy) = 0.0
      cgrid%fmean_available_water      (ipy) = 0.0
      cgrid%fmean_can_theiv            (ipy) = 0.0
      cgrid%fmean_can_theta            (ipy) = 0.0
      cgrid%fmean_can_vpdef            (ipy) = 0.0
      cgrid%fmean_can_temp             (ipy) = 0.0
      cgrid%fmean_can_shv              (ipy) = 0.0
      cgrid%fmean_can_co2              (ipy) = 0.0
      cgrid%fmean_can_rhos             (ipy) = 0.0
      cgrid%fmean_can_prss             (ipy) = 0.0
      cgrid%fmean_gnd_temp             (ipy) = 0.0
      cgrid%fmean_gnd_shv              (ipy) = 0.0
      cgrid%fmean_can_ggnd             (ipy) = 0.0
      cgrid%fmean_sfcw_depth           (ipy) = 0.0
      cgrid%fmean_sfcw_energy          (ipy) = 0.0
      cgrid%fmean_sfcw_mass            (ipy) = 0.0
      cgrid%fmean_sfcw_temp            (ipy) = 0.0
      cgrid%fmean_sfcw_fliq            (ipy) = 0.0
      cgrid%fmean_rshort_gnd           (ipy) = 0.0
      cgrid%fmean_par_gnd              (ipy) = 0.0
      cgrid%fmean_rlong_gnd            (ipy) = 0.0
      cgrid%fmean_rlongup              (ipy) = 0.0
      cgrid%fmean_parup                (ipy) = 0.0
      cgrid%fmean_nirup                (ipy) = 0.0
      cgrid%fmean_rshortup             (ipy) = 0.0
      cgrid%fmean_rnet                 (ipy) = 0.0
      cgrid%fmean_albedo               (ipy) = 0.0
      cgrid%fmean_albedo_par           (ipy) = 0.0
      cgrid%fmean_albedo_nir           (ipy) = 0.0
      cgrid%fmean_rlong_albedo         (ipy) = 0.0
      cgrid%fmean_ustar                (ipy) = 0.0
      cgrid%fmean_tstar                (ipy) = 0.0
      cgrid%fmean_qstar                (ipy) = 0.0
      cgrid%fmean_cstar                (ipy) = 0.0
      cgrid%fmean_carbon_ac            (ipy) = 0.0
      cgrid%fmean_carbon_st            (ipy) = 0.0
      cgrid%fmean_vapor_gc             (ipy) = 0.0
      cgrid%fmean_vapor_ac             (ipy) = 0.0
      cgrid%fmean_throughfall          (ipy) = 0.0
      cgrid%fmean_runoff               (ipy) = 0.0
      cgrid%fmean_drainage             (ipy) = 0.0
      cgrid%fmean_sensible_gc          (ipy) = 0.0
      cgrid%fmean_sensible_ac          (ipy) = 0.0
      cgrid%fmean_qthroughfall         (ipy) = 0.0
      cgrid%fmean_qrunoff              (ipy) = 0.0
      cgrid%fmean_qdrainage            (ipy) = 0.0
      cgrid%fmean_atm_theiv            (ipy) = 0.0
      cgrid%fmean_atm_theta            (ipy) = 0.0
      cgrid%fmean_atm_temp             (ipy) = 0.0
      cgrid%fmean_atm_vpdef            (ipy) = 0.0
      cgrid%fmean_atm_shv              (ipy) = 0.0
      cgrid%fmean_atm_rshort           (ipy) = 0.0
      cgrid%fmean_atm_rshort_diff      (ipy) = 0.0
      cgrid%fmean_atm_par              (ipy) = 0.0
      cgrid%fmean_atm_par_diff         (ipy) = 0.0
      cgrid%fmean_atm_rlong            (ipy) = 0.0
      cgrid%fmean_atm_vels             (ipy) = 0.0
      cgrid%fmean_atm_rhos             (ipy) = 0.0
      cgrid%fmean_atm_prss             (ipy) = 0.0
      cgrid%fmean_atm_co2              (ipy) = 0.0
      cgrid%fmean_pcpg                 (ipy) = 0.0
      cgrid%fmean_qpcpg                (ipy) = 0.0
      cgrid%fmean_dpcpg                (ipy) = 0.0
      cgrid%fmean_soil_wetness         (ipy) = 0.0
      cgrid%fmean_skin_temp            (ipy) = 0.0
      cgrid%fmean_soil_energy        (:,ipy) = 0.0
      cgrid%fmean_soil_mstpot        (:,ipy) = 0.0
      cgrid%fmean_soil_water         (:,ipy) = 0.0
      cgrid%fmean_soil_temp          (:,ipy) = 0.0
      cgrid%fmean_soil_fliq          (:,ipy) = 0.0
      cgrid%fmean_smoist_gg          (:,ipy) = 0.0
      cgrid%fmean_transloss          (:,ipy) = 0.0
      cgrid%fmean_sensible_gg        (:,ipy) = 0.0
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Daily means.                                                                   !
      !------------------------------------------------------------------------------------!
      if (writing_long) then
         cgrid%dmean_nppleaf              (ipy) = 0.0
         cgrid%dmean_nppfroot             (ipy) = 0.0
         cgrid%dmean_nppsapwood           (ipy) = 0.0
         cgrid%dmean_nppcroot             (ipy) = 0.0
         cgrid%dmean_nppseeds             (ipy) = 0.0
         cgrid%dmean_nppwood              (ipy) = 0.0
         cgrid%dmean_nppdaily             (ipy) = 0.0
         cgrid%dmean_A_decomp           (:,ipy) = 0.0
         cgrid%dmean_Af_decomp          (:,ipy) = 0.0
         cgrid%dmean_co2_residual         (ipy) = 0.0
         cgrid%dmean_energy_residual      (ipy) = 0.0
         cgrid%dmean_water_residual       (ipy) = 0.0
         cgrid%dmean_gpp                  (ipy) = 0.0
         cgrid%dmean_npp                  (ipy) = 0.0
         cgrid%dmean_leaf_resp            (ipy) = 0.0
         cgrid%dmean_root_resp            (ipy) = 0.0
         cgrid%dmean_leaf_growth_resp     (ipy) = 0.0
         cgrid%dmean_root_growth_resp     (ipy) = 0.0
         cgrid%dmean_sapa_growth_resp     (ipy) = 0.0
         cgrid%dmean_sapb_growth_resp     (ipy) = 0.0
         cgrid%dmean_leaf_storage_resp    (ipy) = 0.0
         cgrid%dmean_root_storage_resp    (ipy) = 0.0
         cgrid%dmean_sapa_storage_resp    (ipy) = 0.0
         cgrid%dmean_sapb_storage_resp    (ipy) = 0.0
         cgrid%dmean_plresp               (ipy) = 0.0
         cgrid%dmean_leaf_energy          (ipy) = 0.0
         cgrid%dmean_leaf_water           (ipy) = 0.0
         cgrid%dmean_leaf_hcap            (ipy) = 0.0
         cgrid%dmean_leaf_vpdef           (ipy) = 0.0
         cgrid%dmean_leaf_temp            (ipy) = 0.0
         cgrid%dmean_leaf_fliq            (ipy) = 0.0
         cgrid%dmean_leaf_gsw             (ipy) = 0.0
         cgrid%dmean_leaf_gbw             (ipy) = 0.0
         cgrid%dmean_wood_energy          (ipy) = 0.0
         cgrid%dmean_wood_water           (ipy) = 0.0
         cgrid%dmean_wood_hcap            (ipy) = 0.0
         cgrid%dmean_wood_temp            (ipy) = 0.0
         cgrid%dmean_wood_fliq            (ipy) = 0.0
         cgrid%dmean_wood_gbw             (ipy) = 0.0
         cgrid%dmean_fs_open              (ipy) = 0.0
         cgrid%dmean_fsw                  (ipy) = 0.0
         cgrid%dmean_fsn                  (ipy) = 0.0
         cgrid%dmean_a_open               (ipy) = 0.0
         cgrid%dmean_a_closed             (ipy) = 0.0
         cgrid%dmean_a_net                (ipy) = 0.0
         cgrid%dmean_a_light              (ipy) = 0.0
         cgrid%dmean_a_rubp               (ipy) = 0.0
         cgrid%dmean_a_co2                (ipy) = 0.0
         cgrid%dmean_psi_open             (ipy) = 0.0
         cgrid%dmean_psi_closed           (ipy) = 0.0
         cgrid%dmean_water_supply         (ipy) = 0.0
         cgrid%dmean_par_l                (ipy) = 0.0
         cgrid%dmean_par_l_beam           (ipy) = 0.0
         cgrid%dmean_par_l_diff           (ipy) = 0.0
         cgrid%dmean_rshort_l             (ipy) = 0.0
         cgrid%dmean_rlong_l              (ipy) = 0.0
         cgrid%dmean_sensible_lc          (ipy) = 0.0
         cgrid%dmean_vapor_lc             (ipy) = 0.0
         cgrid%dmean_transp               (ipy) = 0.0
         cgrid%dmean_intercepted_al       (ipy) = 0.0
         cgrid%dmean_wshed_lg             (ipy) = 0.0
         cgrid%dmean_rshort_w             (ipy) = 0.0
         cgrid%dmean_rlong_w              (ipy) = 0.0
         cgrid%dmean_sensible_wc          (ipy) = 0.0
         cgrid%dmean_vapor_wc             (ipy) = 0.0
         cgrid%dmean_intercepted_aw       (ipy) = 0.0
         cgrid%dmean_wshed_wg             (ipy) = 0.0
         cgrid%dmean_rh                 (:,ipy) = 0.0
         cgrid%dmean_cwd_rh               (ipy) = 0.0
         cgrid%dmean_nep                  (ipy) = 0.0
         cgrid%dmean_rk4step              (ipy) = 0.0
         cgrid%dmean_available_water      (ipy) = 0.0
         cgrid%dmean_can_theiv            (ipy) = 0.0
         cgrid%dmean_can_theta            (ipy) = 0.0
         cgrid%dmean_can_vpdef            (ipy) = 0.0
         cgrid%dmean_can_temp             (ipy) = 0.0
         cgrid%dmean_can_shv              (ipy) = 0.0
         cgrid%dmean_can_co2              (ipy) = 0.0
         cgrid%dmean_can_rhos             (ipy) = 0.0
         cgrid%dmean_can_prss             (ipy) = 0.0
         cgrid%dmean_gnd_temp             (ipy) = 0.0
         cgrid%dmean_gnd_shv              (ipy) = 0.0
         cgrid%dmean_can_ggnd             (ipy) = 0.0
         cgrid%dmean_sfcw_depth           (ipy) = 0.0
         cgrid%dmean_sfcw_energy          (ipy) = 0.0
         cgrid%dmean_sfcw_mass            (ipy) = 0.0
         cgrid%dmean_sfcw_temp            (ipy) = 0.0
         cgrid%dmean_sfcw_fliq            (ipy) = 0.0
         cgrid%dmean_rshort_gnd           (ipy) = 0.0
         cgrid%dmean_par_gnd              (ipy) = 0.0
         cgrid%dmean_rlong_gnd            (ipy) = 0.0
         cgrid%dmean_rlongup              (ipy) = 0.0
         cgrid%dmean_parup                (ipy) = 0.0
         cgrid%dmean_nirup                (ipy) = 0.0
         cgrid%dmean_rshortup             (ipy) = 0.0
         cgrid%dmean_rnet                 (ipy) = 0.0
         cgrid%dmean_albedo               (ipy) = 0.0
         cgrid%dmean_albedo_par           (ipy) = 0.0
         cgrid%dmean_albedo_nir           (ipy) = 0.0
         cgrid%dmean_rlong_albedo         (ipy) = 0.0
         cgrid%dmean_ustar                (ipy) = 0.0
         cgrid%dmean_tstar                (ipy) = 0.0
         cgrid%dmean_qstar                (ipy) = 0.0
         cgrid%dmean_cstar                (ipy) = 0.0
         cgrid%dmean_carbon_ac            (ipy) = 0.0
         cgrid%dmean_carbon_st            (ipy) = 0.0
         cgrid%dmean_vapor_gc             (ipy) = 0.0
         cgrid%dmean_vapor_ac             (ipy) = 0.0
         cgrid%dmean_throughfall          (ipy) = 0.0
         cgrid%dmean_runoff               (ipy) = 0.0
         cgrid%dmean_drainage             (ipy) = 0.0
         cgrid%dmean_sensible_gc          (ipy) = 0.0
         cgrid%dmean_sensible_ac          (ipy) = 0.0
         cgrid%dmean_qthroughfall         (ipy) = 0.0
         cgrid%dmean_qrunoff              (ipy) = 0.0
         cgrid%dmean_qdrainage            (ipy) = 0.0
         cgrid%dmean_atm_theiv            (ipy) = 0.0
         cgrid%dmean_atm_theta            (ipy) = 0.0
         cgrid%dmean_atm_temp             (ipy) = 0.0
         cgrid%dmean_atm_vpdef            (ipy) = 0.0
         cgrid%dmean_atm_shv              (ipy) = 0.0
         cgrid%dmean_atm_rshort           (ipy) = 0.0
         cgrid%dmean_atm_rshort_diff      (ipy) = 0.0
         cgrid%dmean_atm_par              (ipy) = 0.0
         cgrid%dmean_atm_par_diff         (ipy) = 0.0
         cgrid%dmean_atm_rlong            (ipy) = 0.0
         cgrid%dmean_atm_vels             (ipy) = 0.0
         cgrid%dmean_atm_rhos             (ipy) = 0.0
         cgrid%dmean_atm_prss             (ipy) = 0.0
         cgrid%dmean_atm_co2              (ipy) = 0.0
         cgrid%dmean_pcpg                 (ipy) = 0.0
         cgrid%dmean_qpcpg                (ipy) = 0.0
         cgrid%dmean_dpcpg                (ipy) = 0.0
         cgrid%dmean_soil_energy        (:,ipy) = 0.0
         cgrid%dmean_soil_mstpot        (:,ipy) = 0.0
         cgrid%dmean_soil_water         (:,ipy) = 0.0
         cgrid%dmean_soil_temp          (:,ipy) = 0.0
         cgrid%dmean_soil_fliq          (:,ipy) = 0.0
         cgrid%dmean_smoist_gg          (:,ipy) = 0.0
         cgrid%dmean_transloss          (:,ipy) = 0.0
         cgrid%dmean_sensible_gg        (:,ipy) = 0.0
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Monthly means.                                                                 !
      !------------------------------------------------------------------------------------!
      if (writing_eorq) then
         cgrid%mmean_gpp                  (ipy) = 0.0
         cgrid%mmean_npp                  (ipy) = 0.0
         cgrid%mmean_leaf_resp            (ipy) = 0.0
         cgrid%mmean_root_resp            (ipy) = 0.0
         cgrid%mmean_leaf_growth_resp     (ipy) = 0.0
         cgrid%mmean_root_growth_resp     (ipy) = 0.0
         cgrid%mmean_sapa_growth_resp     (ipy) = 0.0
         cgrid%mmean_sapb_growth_resp     (ipy) = 0.0
         cgrid%mmean_leaf_storage_resp    (ipy) = 0.0
         cgrid%mmean_root_storage_resp    (ipy) = 0.0
         cgrid%mmean_sapa_storage_resp    (ipy) = 0.0
         cgrid%mmean_sapb_storage_resp    (ipy) = 0.0
         cgrid%mmean_plresp               (ipy) = 0.0
         cgrid%mmean_leaf_energy          (ipy) = 0.0
         cgrid%mmean_leaf_water           (ipy) = 0.0
         cgrid%mmean_leaf_hcap            (ipy) = 0.0
         cgrid%mmean_leaf_vpdef           (ipy) = 0.0
         cgrid%mmean_leaf_temp            (ipy) = 0.0
         cgrid%mmean_leaf_fliq            (ipy) = 0.0
         cgrid%mmean_leaf_gsw             (ipy) = 0.0
         cgrid%mmean_leaf_gbw             (ipy) = 0.0
         cgrid%mmean_wood_energy          (ipy) = 0.0
         cgrid%mmean_wood_water           (ipy) = 0.0
         cgrid%mmean_wood_hcap            (ipy) = 0.0
         cgrid%mmean_wood_temp            (ipy) = 0.0
         cgrid%mmean_wood_fliq            (ipy) = 0.0
         cgrid%mmean_wood_gbw             (ipy) = 0.0
         cgrid%mmean_fs_open              (ipy) = 0.0
         cgrid%mmean_fsw                  (ipy) = 0.0
         cgrid%mmean_fsn                  (ipy) = 0.0
         cgrid%mmean_a_open               (ipy) = 0.0
         cgrid%mmean_a_closed             (ipy) = 0.0
         cgrid%mmean_a_net                (ipy) = 0.0
         cgrid%mmean_a_light              (ipy) = 0.0
         cgrid%mmean_a_rubp               (ipy) = 0.0
         cgrid%mmean_a_co2                (ipy) = 0.0
         cgrid%mmean_psi_open             (ipy) = 0.0
         cgrid%mmean_psi_closed           (ipy) = 0.0
         cgrid%mmean_water_supply         (ipy) = 0.0
         cgrid%mmean_par_l                (ipy) = 0.0
         cgrid%mmean_par_l_beam           (ipy) = 0.0
         cgrid%mmean_par_l_diff           (ipy) = 0.0
         cgrid%mmean_rshort_l             (ipy) = 0.0
         cgrid%mmean_rlong_l              (ipy) = 0.0
         cgrid%mmean_sensible_lc          (ipy) = 0.0
         cgrid%mmean_vapor_lc             (ipy) = 0.0
         cgrid%mmean_transp               (ipy) = 0.0
         cgrid%mmean_intercepted_al       (ipy) = 0.0
         cgrid%mmean_wshed_lg             (ipy) = 0.0
         cgrid%mmean_rshort_w             (ipy) = 0.0
         cgrid%mmean_rlong_w              (ipy) = 0.0
         cgrid%mmean_sensible_wc          (ipy) = 0.0
         cgrid%mmean_vapor_wc             (ipy) = 0.0
         cgrid%mmean_intercepted_aw       (ipy) = 0.0
         cgrid%mmean_wshed_wg             (ipy) = 0.0
         cgrid%mmean_rh                 (:,ipy) = 0.0
         cgrid%mmean_cwd_rh               (ipy) = 0.0
         cgrid%mmean_nep                  (ipy) = 0.0
         cgrid%mmean_rk4step              (ipy) = 0.0
         cgrid%mmean_available_water      (ipy) = 0.0
         cgrid%mmean_can_theiv            (ipy) = 0.0
         cgrid%mmean_can_theta            (ipy) = 0.0
         cgrid%mmean_can_vpdef            (ipy) = 0.0
         cgrid%mmean_can_temp             (ipy) = 0.0
         cgrid%mmean_can_shv              (ipy) = 0.0
         cgrid%mmean_can_co2              (ipy) = 0.0
         cgrid%mmean_can_rhos             (ipy) = 0.0
         cgrid%mmean_can_prss             (ipy) = 0.0
         cgrid%mmean_gnd_temp             (ipy) = 0.0
         cgrid%mmean_gnd_shv              (ipy) = 0.0
         cgrid%mmean_can_ggnd             (ipy) = 0.0
         cgrid%mmean_sfcw_depth           (ipy) = 0.0
         cgrid%mmean_sfcw_energy          (ipy) = 0.0
         cgrid%mmean_sfcw_mass            (ipy) = 0.0
         cgrid%mmean_sfcw_temp            (ipy) = 0.0
         cgrid%mmean_sfcw_fliq            (ipy) = 0.0
         cgrid%mmean_rshort_gnd           (ipy) = 0.0
         cgrid%mmean_par_gnd              (ipy) = 0.0
         cgrid%mmean_rlong_gnd            (ipy) = 0.0
         cgrid%mmean_rlongup              (ipy) = 0.0
         cgrid%mmean_parup                (ipy) = 0.0
         cgrid%mmean_nirup                (ipy) = 0.0
         cgrid%mmean_rshortup             (ipy) = 0.0
         cgrid%mmean_rnet                 (ipy) = 0.0
         cgrid%mmean_albedo               (ipy) = 0.0
         cgrid%mmean_albedo_par           (ipy) = 0.0
         cgrid%mmean_albedo_nir           (ipy) = 0.0
         cgrid%mmean_rlong_albedo         (ipy) = 0.0
         cgrid%mmean_ustar                (ipy) = 0.0
         cgrid%mmean_tstar                (ipy) = 0.0
         cgrid%mmean_qstar                (ipy) = 0.0
         cgrid%mmean_cstar                (ipy) = 0.0
         cgrid%mmean_carbon_ac            (ipy) = 0.0
         cgrid%mmean_carbon_st            (ipy) = 0.0
         cgrid%mmean_vapor_gc             (ipy) = 0.0
         cgrid%mmean_vapor_ac             (ipy) = 0.0
         cgrid%mmean_throughfall          (ipy) = 0.0
         cgrid%mmean_runoff               (ipy) = 0.0
         cgrid%mmean_drainage             (ipy) = 0.0
         cgrid%mmean_sensible_gc          (ipy) = 0.0
         cgrid%mmean_sensible_ac          (ipy) = 0.0
         cgrid%mmean_qthroughfall         (ipy) = 0.0
         cgrid%mmean_qrunoff              (ipy) = 0.0
         cgrid%mmean_qdrainage            (ipy) = 0.0
         cgrid%mmean_soil_energy        (:,ipy) = 0.0
         cgrid%mmean_soil_mstpot        (:,ipy) = 0.0
         cgrid%mmean_soil_water         (:,ipy) = 0.0
         cgrid%mmean_soil_temp          (:,ipy) = 0.0
         cgrid%mmean_soil_fliq          (:,ipy) = 0.0
         cgrid%mmean_smoist_gg          (:,ipy) = 0.0
         cgrid%mmean_transloss          (:,ipy) = 0.0
         cgrid%mmean_sensible_gg        (:,ipy) = 0.0
         cgrid%mmean_lai              (:,:,ipy) = 0.0
         cgrid%mmean_bleaf            (:,:,ipy) = 0.0
         cgrid%mmean_broot            (:,:,ipy) = 0.0
         cgrid%mmean_bstorage         (:,:,ipy) = 0.0
         cgrid%mmean_bleaf_n          (:,:,ipy) = 0.0
         cgrid%mmean_broot_n          (:,:,ipy) = 0.0
         cgrid%mmean_bstorage_n       (:,:,ipy) = 0.0
         cgrid%mmean_leaf_maintenance (:,:,ipy) = 0.0
         cgrid%mmean_root_maintenance (:,:,ipy) = 0.0
         cgrid%mmean_leaf_drop        (:,:,ipy) = 0.0
         cgrid%mmean_fast_soil_c          (ipy) = 0.0
         cgrid%mmean_slow_soil_c          (ipy) = 0.0
         cgrid%mmean_struct_soil_c        (ipy) = 0.0
         cgrid%mmean_struct_soil_l        (ipy) = 0.0
         cgrid%mmean_cwd_c                (ipy) = 0.0
         cgrid%mmean_fast_soil_n          (ipy) = 0.0
         cgrid%mmean_mineral_soil_n       (ipy) = 0.0
         cgrid%mmean_cwd_n                (ipy) = 0.0
         cgrid%mmean_nppleaf              (ipy) = 0.0
         cgrid%mmean_nppfroot             (ipy) = 0.0
         cgrid%mmean_nppsapwood           (ipy) = 0.0
         cgrid%mmean_nppcroot             (ipy) = 0.0
         cgrid%mmean_nppseeds             (ipy) = 0.0
         cgrid%mmean_nppwood              (ipy) = 0.0
         cgrid%mmean_nppdaily             (ipy) = 0.0
         cgrid%mmean_A_decomp           (:,ipy) = 0.0
         cgrid%mmean_Af_decomp          (:,ipy) = 0.0
         cgrid%mmean_co2_residual         (ipy) = 0.0
         cgrid%mmean_energy_residual      (ipy) = 0.0
         cgrid%mmean_water_residual       (ipy) = 0.0
         cgrid%mmean_atm_theiv            (ipy) = 0.0
         cgrid%mmean_atm_theta            (ipy) = 0.0
         cgrid%mmean_atm_temp             (ipy) = 0.0
         cgrid%mmean_atm_vpdef            (ipy) = 0.0
         cgrid%mmean_atm_shv              (ipy) = 0.0
         cgrid%mmean_atm_rshort           (ipy) = 0.0
         cgrid%mmean_atm_rshort_diff      (ipy) = 0.0
         cgrid%mmean_atm_par              (ipy) = 0.0
         cgrid%mmean_atm_par_diff         (ipy) = 0.0
         cgrid%mmean_atm_rlong            (ipy) = 0.0
         cgrid%mmean_atm_vels             (ipy) = 0.0
         cgrid%mmean_atm_rhos             (ipy) = 0.0
         cgrid%mmean_atm_prss             (ipy) = 0.0
         cgrid%mmean_atm_co2              (ipy) = 0.0
         cgrid%mmean_pcpg                 (ipy) = 0.0
         cgrid%mmean_qpcpg                (ipy) = 0.0
         cgrid%mmean_dpcpg                (ipy) = 0.0
         cgrid%mmsqu_gpp                  (ipy) = 0.0
         cgrid%mmsqu_npp                  (ipy) = 0.0
         cgrid%mmsqu_plresp               (ipy) = 0.0
         cgrid%mmsqu_sensible_lc          (ipy) = 0.0
         cgrid%mmsqu_vapor_lc             (ipy) = 0.0
         cgrid%mmsqu_transp               (ipy) = 0.0
         cgrid%mmsqu_sensible_wc          (ipy) = 0.0
         cgrid%mmsqu_vapor_wc             (ipy) = 0.0
         cgrid%mmsqu_rh                   (ipy) = 0.0
         cgrid%mmsqu_cwd_rh               (ipy) = 0.0
         cgrid%mmsqu_nep                  (ipy) = 0.0
         cgrid%mmsqu_rlongup              (ipy) = 0.0
         cgrid%mmsqu_parup                (ipy) = 0.0
         cgrid%mmsqu_nirup                (ipy) = 0.0
         cgrid%mmsqu_rshortup             (ipy) = 0.0
         cgrid%mmsqu_rnet                 (ipy) = 0.0
         cgrid%mmsqu_albedo               (ipy) = 0.0
         cgrid%mmsqu_ustar                (ipy) = 0.0
         cgrid%mmsqu_carbon_ac            (ipy) = 0.0
         cgrid%mmsqu_carbon_st            (ipy) = 0.0
         cgrid%mmsqu_vapor_gc             (ipy) = 0.0
         cgrid%mmsqu_vapor_ac             (ipy) = 0.0
         cgrid%mmsqu_sensible_gc          (ipy) = 0.0
         cgrid%mmsqu_sensible_ac          (ipy) = 0.0
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Mean diel.                                                                     !
      !------------------------------------------------------------------------------------!
      if (writing_dcyc) then
         cgrid%qmean_gpp                (:,ipy) = 0.0
         cgrid%qmean_npp                (:,ipy) = 0.0
         cgrid%qmean_leaf_resp          (:,ipy) = 0.0
         cgrid%qmean_root_resp          (:,ipy) = 0.0
         cgrid%qmean_leaf_growth_resp   (:,ipy) = 0.0
         cgrid%qmean_root_growth_resp   (:,ipy) = 0.0
         cgrid%qmean_sapa_growth_resp   (:,ipy) = 0.0
         cgrid%qmean_sapb_growth_resp   (:,ipy) = 0.0
         cgrid%qmean_leaf_storage_resp  (:,ipy) = 0.0
         cgrid%qmean_root_storage_resp  (:,ipy) = 0.0
         cgrid%qmean_sapa_storage_resp  (:,ipy) = 0.0
         cgrid%qmean_sapb_storage_resp  (:,ipy) = 0.0
         cgrid%qmean_plresp             (:,ipy) = 0.0
         cgrid%qmean_leaf_energy        (:,ipy) = 0.0
         cgrid%qmean_leaf_water         (:,ipy) = 0.0
         cgrid%qmean_leaf_hcap          (:,ipy) = 0.0
         cgrid%qmean_leaf_vpdef         (:,ipy) = 0.0
         cgrid%qmean_leaf_temp          (:,ipy) = 0.0
         cgrid%qmean_leaf_fliq          (:,ipy) = 0.0
         cgrid%qmean_leaf_gsw           (:,ipy) = 0.0
         cgrid%qmean_leaf_gbw           (:,ipy) = 0.0
         cgrid%qmean_wood_energy        (:,ipy) = 0.0
         cgrid%qmean_wood_water         (:,ipy) = 0.0
         cgrid%qmean_wood_hcap          (:,ipy) = 0.0
         cgrid%qmean_wood_temp          (:,ipy) = 0.0
         cgrid%qmean_wood_fliq          (:,ipy) = 0.0
         cgrid%qmean_wood_gbw           (:,ipy) = 0.0
         cgrid%qmean_fs_open            (:,ipy) = 0.0
         cgrid%qmean_fsw                (:,ipy) = 0.0
         cgrid%qmean_fsn                (:,ipy) = 0.0
         cgrid%qmean_a_open             (:,ipy) = 0.0
         cgrid%qmean_a_closed           (:,ipy) = 0.0
         cgrid%qmean_a_net              (:,ipy) = 0.0
         cgrid%qmean_a_light            (:,ipy) = 0.0
         cgrid%qmean_a_rubp             (:,ipy) = 0.0
         cgrid%qmean_a_co2              (:,ipy) = 0.0
         cgrid%qmean_psi_open           (:,ipy) = 0.0
         cgrid%qmean_psi_closed         (:,ipy) = 0.0
         cgrid%qmean_water_supply       (:,ipy) = 0.0
         cgrid%qmean_par_l              (:,ipy) = 0.0
         cgrid%qmean_par_l_beam         (:,ipy) = 0.0
         cgrid%qmean_par_l_diff         (:,ipy) = 0.0
         cgrid%qmean_rshort_l           (:,ipy) = 0.0
         cgrid%qmean_rlong_l            (:,ipy) = 0.0
         cgrid%qmean_sensible_lc        (:,ipy) = 0.0
         cgrid%qmean_vapor_lc           (:,ipy) = 0.0
         cgrid%qmean_transp             (:,ipy) = 0.0
         cgrid%qmean_intercepted_al     (:,ipy) = 0.0
         cgrid%qmean_wshed_lg           (:,ipy) = 0.0
         cgrid%qmean_rshort_w           (:,ipy) = 0.0
         cgrid%qmean_rlong_w            (:,ipy) = 0.0
         cgrid%qmean_sensible_wc        (:,ipy) = 0.0
         cgrid%qmean_vapor_wc           (:,ipy) = 0.0
         cgrid%qmean_intercepted_aw     (:,ipy) = 0.0
         cgrid%qmean_wshed_wg           (:,ipy) = 0.0
         cgrid%qmean_rh                 (:,ipy) = 0.0
         cgrid%qmean_cwd_rh             (:,ipy) = 0.0
         cgrid%qmean_nep                (:,ipy) = 0.0
         cgrid%qmean_rk4step            (:,ipy) = 0.0
         cgrid%qmean_available_water    (:,ipy) = 0.0
         cgrid%qmean_can_theiv          (:,ipy) = 0.0
         cgrid%qmean_can_theta          (:,ipy) = 0.0
         cgrid%qmean_can_vpdef          (:,ipy) = 0.0
         cgrid%qmean_can_temp           (:,ipy) = 0.0
         cgrid%qmean_can_shv            (:,ipy) = 0.0
         cgrid%qmean_can_co2            (:,ipy) = 0.0
         cgrid%qmean_can_rhos           (:,ipy) = 0.0
         cgrid%qmean_can_prss           (:,ipy) = 0.0
         cgrid%qmean_gnd_temp           (:,ipy) = 0.0
         cgrid%qmean_gnd_shv            (:,ipy) = 0.0
         cgrid%qmean_can_ggnd           (:,ipy) = 0.0
         cgrid%qmean_sfcw_depth         (:,ipy) = 0.0
         cgrid%qmean_sfcw_energy        (:,ipy) = 0.0
         cgrid%qmean_sfcw_mass          (:,ipy) = 0.0
         cgrid%qmean_sfcw_temp          (:,ipy) = 0.0
         cgrid%qmean_sfcw_fliq          (:,ipy) = 0.0
         cgrid%qmean_rshort_gnd         (:,ipy) = 0.0
         cgrid%qmean_par_gnd            (:,ipy) = 0.0
         cgrid%qmean_rlong_gnd          (:,ipy) = 0.0
         cgrid%qmean_rlongup            (:,ipy) = 0.0
         cgrid%qmean_parup              (:,ipy) = 0.0
         cgrid%qmean_nirup              (:,ipy) = 0.0
         cgrid%qmean_rshortup           (:,ipy) = 0.0
         cgrid%qmean_rnet               (:,ipy) = 0.0
         cgrid%qmean_albedo             (:,ipy) = 0.0
         cgrid%qmean_albedo_par         (:,ipy) = 0.0
         cgrid%qmean_albedo_nir         (:,ipy) = 0.0
         cgrid%qmean_rlong_albedo       (:,ipy) = 0.0
         cgrid%qmean_ustar              (:,ipy) = 0.0
         cgrid%qmean_tstar              (:,ipy) = 0.0
         cgrid%qmean_qstar              (:,ipy) = 0.0
         cgrid%qmean_cstar              (:,ipy) = 0.0
         cgrid%qmean_carbon_ac          (:,ipy) = 0.0
         cgrid%qmean_carbon_st          (:,ipy) = 0.0
         cgrid%qmean_vapor_gc           (:,ipy) = 0.0
         cgrid%qmean_vapor_ac           (:,ipy) = 0.0
         cgrid%qmean_throughfall        (:,ipy) = 0.0
         cgrid%qmean_runoff             (:,ipy) = 0.0
         cgrid%qmean_drainage           (:,ipy) = 0.0
         cgrid%qmean_sensible_gc        (:,ipy) = 0.0
         cgrid%qmean_sensible_ac        (:,ipy) = 0.0
         cgrid%qmean_qthroughfall       (:,ipy) = 0.0
         cgrid%qmean_qrunoff            (:,ipy) = 0.0
         cgrid%qmean_qdrainage          (:,ipy) = 0.0
         cgrid%qmean_soil_energy      (:,:,ipy) = 0.0
         cgrid%qmean_soil_mstpot      (:,:,ipy) = 0.0
         cgrid%qmean_soil_water       (:,:,ipy) = 0.0
         cgrid%qmean_soil_temp        (:,:,ipy) = 0.0
         cgrid%qmean_soil_fliq        (:,:,ipy) = 0.0
         cgrid%qmean_smoist_gg        (:,:,ipy) = 0.0
         cgrid%qmean_transloss        (:,:,ipy) = 0.0
         cgrid%qmean_sensible_gg      (:,:,ipy) = 0.0
         cgrid%qmean_atm_theiv          (:,ipy) = 0.0
         cgrid%qmean_atm_theta          (:,ipy) = 0.0
         cgrid%qmean_atm_temp           (:,ipy) = 0.0
         cgrid%qmean_atm_vpdef          (:,ipy) = 0.0
         cgrid%qmean_atm_shv            (:,ipy) = 0.0
         cgrid%qmean_atm_rshort         (:,ipy) = 0.0
         cgrid%qmean_atm_rshort_diff    (:,ipy) = 0.0
         cgrid%qmean_atm_par            (:,ipy) = 0.0
         cgrid%qmean_atm_par_diff       (:,ipy) = 0.0
         cgrid%qmean_atm_rlong          (:,ipy) = 0.0
         cgrid%qmean_atm_vels           (:,ipy) = 0.0
         cgrid%qmean_atm_rhos           (:,ipy) = 0.0
         cgrid%qmean_atm_prss           (:,ipy) = 0.0
         cgrid%qmean_atm_co2            (:,ipy) = 0.0
         cgrid%qmean_pcpg               (:,ipy) = 0.0
         cgrid%qmean_qpcpg              (:,ipy) = 0.0
         cgrid%qmean_dpcpg              (:,ipy) = 0.0
         cgrid%qmsqu_gpp                (:,ipy) = 0.0
         cgrid%qmsqu_npp                (:,ipy) = 0.0
         cgrid%qmsqu_plresp             (:,ipy) = 0.0
         cgrid%qmsqu_sensible_lc        (:,ipy) = 0.0
         cgrid%qmsqu_vapor_lc           (:,ipy) = 0.0
         cgrid%qmsqu_transp             (:,ipy) = 0.0
         cgrid%qmsqu_sensible_wc        (:,ipy) = 0.0
         cgrid%qmsqu_vapor_wc           (:,ipy) = 0.0
         cgrid%qmsqu_rh                 (:,ipy) = 0.0
         cgrid%qmsqu_cwd_rh             (:,ipy) = 0.0
         cgrid%qmsqu_nep                (:,ipy) = 0.0
         cgrid%qmsqu_rlongup            (:,ipy) = 0.0
         cgrid%qmsqu_parup              (:,ipy) = 0.0
         cgrid%qmsqu_nirup              (:,ipy) = 0.0
         cgrid%qmsqu_rshortup           (:,ipy) = 0.0
         cgrid%qmsqu_rnet               (:,ipy) = 0.0
         cgrid%qmsqu_albedo             (:,ipy) = 0.0
         cgrid%qmsqu_ustar              (:,ipy) = 0.0
         cgrid%qmsqu_carbon_ac          (:,ipy) = 0.0
         cgrid%qmsqu_carbon_st          (:,ipy) = 0.0
         cgrid%qmsqu_vapor_gc           (:,ipy) = 0.0
         cgrid%qmsqu_vapor_ac           (:,ipy) = 0.0
         cgrid%qmsqu_sensible_gc        (:,ipy) = 0.0
         cgrid%qmsqu_sensible_ac        (:,ipy) = 0.0
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ed_poly_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will assign the values of some diagnostic variables, such as soil    !
! and temporary layer temperature and liquid fraction, and the surface properties.         !
!------------------------------------------------------------------------------------------!
subroutine new_patch_sfc_props(csite,ipa,mzg,mzs,ntext_soil)
   use ed_state_vars , only : sitetype           & ! structure
                            , patchtype          ! ! structure
   use soil_coms     , only : soil               & ! intent(in), look-up table
                            , tiny_sfcwater_mass & ! intent(in)
                            , soil_rough         & ! intent(in)
                            , ny07_eq04_a        & ! intent(in)
                            , ny07_eq04_m        & ! intent(in)
                            , matric_potential   ! ! intent(in)
   use consts_coms   , only : wdns               & ! intent(in)
                            , fsdns              & ! intent(in)
                            , fsdnsi             ! ! intent(in)
   use therm_lib     , only : uextcm2tl          & ! subroutine
                            , uint2tl            ! ! subroutine
   use ed_therm_lib  , only : ed_grndvap         ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)                 , target     :: csite          ! Current site
   integer                        , intent(in) :: ipa            ! Current patch #
   integer                        , intent(in) :: mzg            ! # of soil layers
   integer                        , intent(in) :: mzs            ! # of sfc. water layers
   integer        , dimension(mzg), intent(in) :: ntext_soil     ! Soil texture
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: k              ! Layer counter
   integer                                     :: nsoil          ! Soil texture class
   real                                        :: tot_sfcw_mass  ! Total mass
   real                                        :: bulk_sfcw_dens ! Bulk density
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find soil temperature and liquid water fraction, and soil matric potential.      !
   !---------------------------------------------------------------------------------------!
   do k = 1, mzg
      nsoil = ntext_soil(k)
      call uextcm2tl(csite%soil_energy(k,ipa), csite%soil_water(k,ipa)*wdns                &
                    ,soil(nsoil)%slcpd, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa))
      csite%soil_mstpot(k,ipa) = matric_potential(nsoil,csite%soil_water(k,ipa))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Determine the number of temporary snow/surface water layers.  This is done by       !
   ! checking the mass.  In case there is a layer, we convert sfcwater_energy from J/m2 to !
   ! J/kg, and compute the temperature and liquid fraction.                                !
   !---------------------------------------------------------------------------------------!
   csite%nlev_sfcwater   (ipa) = 0
   tot_sfcw_mass               = 0.
   csite%total_sfcw_depth(ipa) = 0.
   snowloop: do k=1,mzs
      !----- Leave the loop if there is not enough mass in this layer... ------------------!
      if (csite%sfcwater_mass(k,ipa) <= tiny_sfcwater_mass)  exit snowloop
      csite%nlev_sfcwater(ipa) = k
      tot_sfcw_mass                = tot_sfcw_mass  + csite%sfcwater_mass (k,ipa)
      csite%total_sfcw_depth(ipa)  = csite%total_sfcw_depth(ipa)                           &
                                   + csite%sfcwater_depth(k,ipa)
      csite%sfcwater_energy(k,ipa) = csite%sfcwater_energy(k,ipa)                          &
                                   / csite%sfcwater_mass(k,ipa)
      call uint2tl(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa)                &
                  ,csite%sfcwater_fracliq(k,ipa))
   end do snowloop
   !---------------------------------------------------------------------------------------!
   !     Now, just to be safe, we will assign zeroes to layers above.                      !
   !---------------------------------------------------------------------------------------!
   do k=csite%nlev_sfcwater(ipa)+1,mzs
      csite%sfcwater_mass(k,ipa)   = 0.
      csite%sfcwater_energy(k,ipa) = 0.
      csite%sfcwater_depth(k,ipa)  = 0.
      if (k == 1) then
         csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(mzg,ipa)
         csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(mzg,ipa)
      else
         csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
         csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the fraction of the canopy covered in snow.  I could not find any            !
   ! reference for the original method (commented out), so I implemented the method used   !
   ! in CLM-4, which is based on:                                                          !
   !                                                                                       !
   ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover     !
   !    fraction and its evaluation over large North American river basins,                !
   !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                            !
   !---------------------------------------------------------------------------------------!
   ! csite%snowfac(ipa) = min(0.99, csite%total_sfcw_depth(ipa)/csite%veg_height(ipa))
   if (tot_sfcw_mass > tiny_sfcwater_mass) then
      bulk_sfcw_dens     = max( fsdns, min( wdns                                           &
                              , tot_sfcw_mass / csite%total_sfcw_depth(ipa)))
      csite%snowfac(ipa) = max( 0.0, min( 0.99                                             &
                              , tanh( csite%total_sfcw_depth(ipa)                          &
                                    / ( ny07_eq04_a * soil_rough                           &
                                      * (bulk_sfcw_dens * fsdnsi) ** ny07_eq04_m ) ) ) )
   else
      csite%snowfac(ipa) = 0.0
   end if
   !---------------------------------------------------------------------------------------!


   !----- Now we can compute the surface properties. --------------------------------------!
   k=max(1,csite%nlev_sfcwater(ipa))
   call ed_grndvap(csite%nlev_sfcwater(ipa),ntext_soil(mzg)                                &
                  ,csite%soil_water(mzg,ipa),csite%soil_tempk(mzg,ipa)                     &
                  ,csite%soil_fracliq(mzg,ipa),csite%sfcwater_tempk(k,ipa)                 &
                  ,csite%snowfac(ipa),csite%can_prss(ipa)                                  &
                  ,csite%can_shv(ipa),csite%ground_shv(ipa),csite%ground_ssh(ipa)          &
                  ,csite%ground_temp(ipa),csite%ground_fliq(ipa),csite%ggsoil(ipa))
   !---------------------------------------------------------------------------------------!

   return
end subroutine new_patch_sfc_props
!==========================================================================================!
!==========================================================================================!
