!==========================================================================================!
!==========================================================================================!
! MODULE: GROWTH_BALIVE
!
!> \brief   Various routines handling plant C and N use given size and allometry.
!> \details Essentially, this file contains dbalive_dt and their libraries of fns.
!> \author  Translated from ED1 by Ryan Knox and Marcos Longo
!> \author  31 Aug 2015 - Big refactoring in commit b8fb585, Daniel Scott
!> \author  6 Oct 2017 - Included bark as an "active" tissue, Marcos Longo
!------------------------------------------------------------------------------------------!
module growth_balive
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: DBALIVE_DT
   !
   !> \brief   Updates living biomass.
   !> \details Calls a variety of subroutines controlling plant carbon balances, C and N
   !>          xfers among plant pools, updating of storage and growth respiration,
   !>          leaf maintenance,and update mortality rates.
   !> \author  Translated from ED1 by Ryan Knox and Marcos Longo
   !> \author  31 Aug 2015 - Big refactoring in commit b8fb585, Daniel Scott
   !> \warning The order of the operations here affect the C/N budgets, so don't
   !>          change it unless you really know what you are doing.
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt(cgrid,gr_tfact0,year_o_day,veget_dyn_on)
      use ed_state_vars       , only : edtype                     & ! structure
                                     , polygontype                & ! structure
                                     , sitetype                   & ! structure
                                     , patchtype                  ! ! structure
      use met_driver_coms     , only : met_driv_state             ! ! structure
      use pft_coms            , only : plant_N_supply_scale       & ! intent(in)
                                     , is_grass                   ! ! intent(in)
      use physiology_coms     , only : N_plant_lim                ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap              & ! function
                                     , update_veg_energy_cweh     ! ! function
      use allometry           , only : area_indices               & ! subroutine
                                     , size2bt                    & ! subroutine
                                     , size2xb                    & ! subroutine
                                     , bl2h                       & ! function
                                     , h2dbh                      & ! function
                                     , ed_biomass                 ! ! function
      use mortality           , only : mortality_rates            ! ! subroutine
      use fuse_fiss_utils     , only : terminate_cohorts          & ! subroutine
                                     , sort_cohorts               ! ! subroutine
      use ed_misc_coms        , only : igrass                     ! ! intent(in)
      use consts_coms         , only : tiny_num                   & ! intent(in)
                                     , r_tol_trunc                ! ! intent(in)
      use stable_cohorts      , only : is_resolvable              ! ! function
      use budget_utils        , only : reset_cbudget_committed    ! ! sub-routine
      use update_derived_utils, only : update_patch_derived_props ! ! sub-routine
      use plant_hydro         , only : rwc2tw                     & ! sub-routine
                                     , twi2twe                    & ! sub-routine
                                     , update_plc                 ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)        , target     :: cgrid        !< the ed grid
      real                , intent(in) :: gr_tfact0    !< Corr. factor for steady growth
      real                , intent(in) :: year_o_day   !< "time factor" i.e. call frequency
      logical             , intent(in) :: veget_dyn_on !< Is vegetation dynamics on?
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)   , pointer    :: cpoly
      type(met_driv_state), pointer    :: cmet
      type(sitetype)      , pointer    :: csite
      type(patchtype)     , pointer    :: cpatch
      integer                          :: ipy
      integer                          :: isi
      integer                          :: ipa
      integer                          :: ico
      integer                          :: ipft
      integer                          :: phenstatus_in
      integer                          :: krdepth_in
      integer                          :: xfer_case
      real                             :: metnpp_actual
      real                             :: metnpp_pot
      real                             :: metnpp_lightmax
      real                             :: metnpp_moistmax
      real                             :: metnpp_mlmax
      real                             :: growresp_actual
      real                             :: growresp_pot
      real                             :: growresp_lightmax
      real                             :: growresp_moistmax
      real                             :: growresp_mlmax
      real                             :: npp_actual
      real                             :: npp_pot
      real                             :: tissue_maintenance
      real                             :: storage_maintenance
      real                             :: fgc_in_in
      real                             :: fsc_in_in
      real                             :: stgc_in_in
      real                             :: stsc_in_in
      real                             :: pat_balive_in
      real                             :: pat_bdead_in
      real                             :: pat_bstorage_in
      real                             :: pat_carbon_miss
      real                             :: pat_metnpp_actual
      real                             :: pat_npp_actual
      real                             :: pat_tissue_maintenance
      real                             :: bleaf_in
      real                             :: broot_in
      real                             :: bsapwooda_in
      real                             :: bsapwoodb_in
      real                             :: bbarka_in
      real                             :: bbarkb_in
      real                             :: balive_in
      real                             :: bdeada_in
      real                             :: bdeadb_in
      real                             :: hite_in
      real                             :: dbh_in
      real                             :: nplant_in
      real                             :: bstorage_in
      real                             :: agb_in
      real                             :: lai_in
      real                             :: wai_in
      real                             :: cai_in
      real                             :: ba_in
      real                             :: nitrogen_supply
      real                             :: dlnndt
      real                             :: old_leaf_hcap
      real                             :: old_wood_hcap
      real                             :: old_leaf_water
      real                             :: old_wood_water
      real                             :: old_leaf_water_im2
      real                             :: old_wood_water_im2
      real                             :: nitrogen_uptake
      real                             :: N_uptake_pot
      real                             :: tr_bleaf
      real                             :: tr_broot
      real                             :: tr_bbarka
      real                             :: tr_bbarkb
      real                             :: tr_bsapwooda
      real                             :: tr_bsapwoodb
      real                             :: tr_bstorage
      real                             :: carbon_debt
      real                             :: balive_aim
      real                             :: elim_nplant
      real                             :: elim_lai
      real                             :: carbon_miss
      logical                          :: flushing
      logical                          :: on_allometry
      !------------------------------------------------------------------------------------!




      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Save patch-level litter inputs before growth balive.  We use these   !
               ! variables to check carbon conservation at the patch level.                !
               !---------------------------------------------------------------------------!
               fgc_in_in              = csite%fgc_in (ipa)
               fsc_in_in              = csite%fsc_in (ipa)
               stgc_in_in             = csite%stgc_in(ipa)
               stsc_in_in             = csite%stsc_in(ipa)
               pat_balive_in          = 0.0
               pat_bdead_in           = 0.0
               pat_bstorage_in        = 0.0
               pat_carbon_miss        = 0.0
               pat_metnpp_actual      = 0.0
               pat_npp_actual         = 0.0
               pat_tissue_maintenance = 0.0
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Reset the patch-level storage and growth respiration, in kgC/m2/day.  !
               ! These variables are used in the fast time steps to release these          !
               ! committed respiration fluxes to the canopy air space.  We use these       !
               ! variables instead of the cohort-level variables because the cohorts may   !
               ! be terminated, in which case their committed carbon losses to CAS will    !
               ! never be accounted, violating carbon conservation.                        !
               !---------------------------------------------------------------------------!
               csite%commit_storage_resp(ipa) = 0.0
               csite%commit_growth_resp (ipa) = 0.0
               !---------------------------------------------------------------------------!


               !----- Loop over cohorts. --------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts


                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)
                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0

                  !------------------------------------------------------------------------!
                  !     This variable should be zero most of the time.  However there are  !
                  ! a few instances in which perfect carbon conservation is not attainable !
                  ! (i.e. plant with extremely negative NPP, no storage, and insufficient  !
                  ! leaf material).  These are extremely rare and may generate a warning,  !
                  ! but not a fatal error.                                                 !
                  !------------------------------------------------------------------------!
                  carbon_miss = 0.0
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Save variables before growth, so we can revert in case             !
                  ! ivegt_dynamics is zero.                                                !
                  !------------------------------------------------------------------------!
                  balive_in       = cpatch%balive          (ico)
                  bdeada_in       = cpatch%bdeada          (ico)
                  bdeadb_in       = cpatch%bdeadb          (ico)
                  bleaf_in        = cpatch%bleaf           (ico)
                  broot_in        = cpatch%broot           (ico)
                  bsapwooda_in    = cpatch%bsapwooda       (ico)
                  bsapwoodb_in    = cpatch%bsapwoodb       (ico)
                  bbarka_in       = cpatch%bbarka          (ico)
                  bbarkb_in       = cpatch%bbarkb          (ico)
                  hite_in         = cpatch%hite            (ico)
                  dbh_in          = cpatch%dbh             (ico)
                  nplant_in       = cpatch%nplant          (ico)
                  bstorage_in     = cpatch%bstorage        (ico)
                  agb_in          = cpatch%agb             (ico)
                  ba_in           = cpatch%basarea         (ico)
                  phenstatus_in   = cpatch%phenology_status(ico)
                  lai_in          = cpatch%lai             (ico)
                  wai_in          = cpatch%wai             (ico)
                  cai_in          = cpatch%crown_area      (ico)
                  krdepth_in      = cpatch%krdepth         (ico)
                  pat_balive_in   = pat_balive_in   + nplant_in * balive_in
                  pat_bdead_in    = pat_bdead_in    + nplant_in * (bdeada_in + bdeadb_in)
                  pat_bstorage_in = pat_bstorage_in + nplant_in * bstorage_in
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Save original heat capacities and water content (needed for the    !
                  ! energy and water budget checks).                                       !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                  old_wood_hcap      = cpatch%wood_hcap     (ico)
                  old_leaf_water     = cpatch%leaf_water    (ico)
                  old_wood_water     = cpatch%wood_water    (ico)
                  old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                  old_wood_water_im2 = cpatch%wood_water_im2(ico)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Compute and apply maintenance costs.                               !
                  !------------------------------------------------------------------------!
                  call update_maintenance(cpatch,ico,year_o_day,csite%avg_daily_temp(ipa)  &
                                         ,tissue_maintenance,storage_maintenance)
                  !------------------------------------------------------------------------!


                  !----- Update patch total storage maintenance costs. --------------------!
                  pat_tissue_maintenance         = pat_tissue_maintenance                  &
                                                 + cpatch%nplant(ico) * tissue_maintenance
                  csite%commit_storage_resp(ipa) = csite%commit_storage_resp(ipa)          &
                                                 + cpatch%nplant(ico) * storage_maintenance
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Find daily metabolic NPP (GPP minus leaf and fine root            !
                  ! respiration.  We quantify both the actual metabolic and the potential  !
                  ! rates.                                                                 !
                  !------------------------------------------------------------------------!
                  call get_metabolic_npp(cpatch,ico,metnpp_actual,metnpp_pot               &
                                        ,metnpp_lightmax,metnpp_moistmax,metnpp_mlmax)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  ! Find growth respiration.                                               !
                  ! MLO - Double check, but I changed the order so growth respiration is   !
                  !       updated before we update daily NPP.  Growth respiration          !
                  !       continues to be released in the day after, just like storage     !
                  !       respiration but the model accounts for the committed carbon      !
                  !       loss.  By updating this before, we can estimate the relative     !
                  !       carbon balance based on the same days for actual and potential   !
                  !       (before there was a one day shift, not a big deal but            !
                  !       inconsistent nonetheless).                                       !
                  !------------------------------------------------------------------------!
                  call get_growth_respiration(cpatch,ico,metnpp_actual,metnpp_pot          &
                                             ,metnpp_lightmax,metnpp_moistmax,metnpp_mlmax &
                                             ,growresp_actual,growresp_pot                 &
                                             ,growresp_lightmax,growresp_moistmax          &
                                             ,growresp_mlmax,csite%commit_growth_resp(ipa))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Find the total NPP and update the carbon balance variables (used  !
                  ! to estimate density-dependent mortality).                              !
                  !------------------------------------------------------------------------!
                  call update_carbon_balances(cpatch,ipa,ico,metnpp_actual,metnpp_pot      &
                                             ,metnpp_lightmax,metnpp_moistmax,metnpp_mlmax &
                                             ,growresp_actual,growresp_pot                 &
                                             ,growresp_lightmax,growresp_moistmax          &
                                             ,growresp_mlmax,tissue_maintenance            &
                                             ,storage_maintenance,npp_actual,npp_pot )
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      Find the percentage loss of xylem conductance if applicable (used !
                  ! to estimate hydraulic failure mortality).                              !
                  !------------------------------------------------------------------------!
                  call update_plc(cpatch,ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and bstorage.             !
                  !------------------------------------------------------------------------!
                  call get_c_xfers(csite,ipa,ico,npp_actual                                &
                                  ,cpoly%green_leaf_factor(ipft,isi),gr_tfact0             &
                                  ,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb,tr_bbarka   &
                                  ,tr_bbarkb,tr_bstorage,carbon_debt,flushing,balive_aim   &
                                  ,carbon_miss,xfer_case)

                  call apply_c_xfers(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda             &
                                    ,tr_bsapwoodb,tr_bbarka,tr_bbarkb,tr_bstorage)
                  
                  call update_today_npp_vars(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda     &
                                            ,tr_bsapwoodb,tr_bbarka,tr_bbarkb              &
                                            ,npp_actual)

                  call update_nitrogen(flushing,ipft,npp_actual,cpatch%nplant(ico)         &
                                       ,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb        &
                                       ,tr_bbarka,tr_bbarkb,tr_bstorage,nitrogen_uptake    &
                                       ,csite%fgn_in(ipa),csite%fsn_in(ipa))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Update transient carbon (metabolic NPP, total NPP and missed       !
                  ! carbon due to cohort extremely negative carbon balance.                !
                  !------------------------------------------------------------------------!
                  pat_metnpp_actual = pat_metnpp_actual + cpatch%nplant(ico) * metnpp_actual
                  pat_npp_actual    = pat_npp_actual    + cpatch%nplant(ico) * npp_actual
                  pat_carbon_miss   = pat_carbon_miss   + cpatch%nplant(ico) * carbon_miss
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !   In case this is new grass, we must update the grass size.            !
                  !------------------------------------------------------------------------!
                  if ( is_grass(ipft).and. igrass==1) then
                     !---------------------------------------------------------------------!
                     !    New grasses may update height and "DBH" every day.               !
                     !---------------------------------------------------------------------!
                     cpatch%hite(ico) = bl2h(cpatch%bleaf(ico), cpatch%sla(ico), ipft)
                     cpatch%dbh(ico)  = h2dbh(cpatch%hite(ico), ipft)
                     !---------------------------------------------------------------------!
                 else
                     !---------------------------------------------------------------------!
                     !                  Update the phenology status.                       !
                     !---------------------------------------------------------------------!
                     on_allometry = ( balive_aim - cpatch%balive(ico) ) <=                 &
                                    ( r_tol_trunc * balive_aim )
                     if (flushing .and. cpatch%elongf(ico) == 1.0 .and. on_allometry) then
                        cpatch%phenology_status(ico) = 0
                     elseif(cpatch%bleaf(ico) < tiny_num .and.                             &
                            (cpatch%phenology_status(ico) == 0 .or.                        &
                             cpatch%phenology_status(ico) == 1)) then
                        cpatch%phenology_status(ico) = 1
                     end if
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,npp_pot,N_uptake_pot               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_nitrogen_uptake(ipa) =                                 &
                                       csite%total_plant_nitrogen_uptake(ipa)              &
                                     + nitrogen_uptake * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!



                  !----- Calculate plant N limitation factor. -----------------------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     nitrogen_supply = plant_N_supply_scale * cpatch%broot(ico)            &
                                     * csite%mineralized_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply                                     &
                                     / (nitrogen_supply + N_uptake_pot)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Update mortality rates.  Notice that the only mortality rate that !
                  ! changes daily is the frost mortality, and the disturbance mortality is !
                  ! not updated here (it is updated in the main disturbance procedure).    !
                  !                                                                        !
                  !      We no longer include mortality rates due to disturbance in the    !
                  ! big-leaf simulations, this is now done at disturbance.f90.             !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ico,csite%avg_daily_temp(ipa),csite%age(ipa) &
                                      ,csite%dist_type(ipa))
                  dlnndt   = - sum(cpatch%mort_rate(1:5,ico))
                  !------------------------------------------------------------------------!

                  !----- Update monthly mortality rates [1/month]. ------------------------!
                  cpatch%monthly_dlnndt(ico) = cpatch%monthly_dlnndt(ico)                  &
                                             + dlnndt * year_o_day
                  !------------------------------------------------------------------------!


                  !----- Updating LAI, WAI, and CAI. --------------------------------------!
                  call area_indices(cpatch, ico)
                  !------------------------------------------------------------------------!



                  !----- Update above-ground biomass. -------------------------------------!
                  cpatch%agb(ico) = ed_biomass(cpatch, ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     In case vegetation dynamics is turned off, overwrite state         !
                  ! variables values with the original values.                             !
                  !------------------------------------------------------------------------!
                  if (.not. veget_dyn_on) then
                     cpatch%balive          (ico) = balive_in
                     cpatch%bdeada          (ico) = bdeada_in
                     cpatch%bdeadb          (ico) = bdeadb_in
                     cpatch%bleaf           (ico) = bleaf_in
                     cpatch%broot           (ico) = broot_in
                     cpatch%bsapwooda       (ico) = bsapwooda_in
                     cpatch%bsapwoodb       (ico) = bsapwoodb_in
                     cpatch%bbarka          (ico) = bbarka_in
                     cpatch%bbarkb          (ico) = bbarkb_in
                     cpatch%hite            (ico) = hite_in
                     cpatch%dbh             (ico) = dbh_in
                     cpatch%nplant          (ico) = nplant_in
                     cpatch%bstorage        (ico) = bstorage_in
                     cpatch%agb             (ico) = agb_in
                     cpatch%basarea         (ico) = ba_in
                     cpatch%phenology_status(ico) = phenstatus_in
                     cpatch%lai             (ico) = lai_in
                     cpatch%wai             (ico) = wai_in
                     cpatch%crown_area      (ico) = cai_in
                     cpatch%krdepth         (ico) = krdepth_in
                  end if
                  !------------------------------------------------------------------------!



                  !----- Update (commercial) timber biomass. ------------------------------!
                  cpatch%btimber(ico) = size2bt( cpatch%dbh(ico),cpatch%hite(ico)          &
                                               , cpatch%bdeada(ico),cpatch%bsapwooda(ico)  &
                                               , cpatch%bbarka(ico),cpatch%pft(ico) )
                  !------------------------------------------------------------------------!



                  !----- Update bark thickness. -------------------------------------------!
                  cpatch%thbark(ico)  = size2xb( cpatch%dbh(ico),cpatch%hite(ico)          &
                                               , cpatch%bbarka(ico),cpatch%bbarkb(ico)     &
                                               , cpatch%sla(ico),cpatch%pft(ico) )
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore, update           !
                  ! vegetation energy and heat capacity.                                   !
                  !------------------------------------------------------------------------!
                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdeada(ico)                 &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                  call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                    &
                             ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                      &
                             ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico)  &
                             ,cpatch%broot(ico),cpatch%dbh(ico),cpatch%pft(ico)            &
                             ,cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico))
                  call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)       &
                              ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)               &
                              ,cpatch%wood_water_im2(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap    &
                                             ,old_leaf_water,old_wood_water                &
                                             ,old_leaf_water_im2,old_wood_water_im2        &
                                             ,.true.,.false.)
                  !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico,.false.,.false.,'dbalive_dt')
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Before we let these plants carry on with their lives, we must check !
                  ! that we can account for all changes in carbon.  In case there is any-  !
                  ! thing missing or excessive, stop the simulation.                       !
                  !------------------------------------------------------------------------!
                  if (veget_dyn_on) then
                     call check_balive_cohort(csite,ipa,ico,bleaf_in,broot_in,bsapwooda_in &
                                             ,bsapwoodb_in,bbarka_in,bbarkb_in,balive_in   &
                                             ,bstorage_in,bdeada_in,bdeadb_in              &
                                             ,phenstatus_in,metnpp_actual,npp_actual       &
                                             ,growresp_actual,tissue_maintenance           &
                                             ,storage_maintenance,carbon_miss,xfer_case)
                  end if
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!


               !------ Update the committed carbon for the upcoming day. ------------------!
               call reset_cbudget_committed(csite,ipa,veget_dyn_on)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Terminate and sort cohorts in case we are using the new grass scheme,  !
               ! as height and biomass may change every day.                               !
               !---------------------------------------------------------------------------!
               if (veget_dyn_on .and. (igrass == 1)) then
                  call terminate_cohorts(csite,ipa,cmet,.false.,elim_nplant,elim_lai)
                  call sort_cohorts(cpatch)
               end if
               !---------------------------------------------------------------------------!


               !----- Update litter. ------------------------------------------------------!
               call update_litter_inputs(csite,ipa)
               !---------------------------------------------------------------------------!



               !----- Update patch LAI, WAI, height, roughness... -------------------------!
               call update_patch_derived_props(csite,ipa,.true.)
               !---------------------------------------------------------------------------!



               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Make sure that the patch did not try to smuggle or evade carbon.      !
               ! For this we compare the live stocks and the necromass inputs before and   !
               ! after updating living tissues.                                            !
               !---------------------------------------------------------------------------!
               if (veget_dyn_on) then
                  call check_balive_patch(csite,ipa,fgc_in_in,fsc_in_in,stgc_in_in         &
                                         ,stsc_in_in,pat_balive_in,pat_bdead_in            &
                                         ,pat_bstorage_in,pat_metnpp_actual,pat_npp_actual &
                                         ,pat_tissue_maintenance,pat_carbon_miss)
               end if
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine dbalive_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Obtain tissue and storage losses due to maintenance costs (aka turnover), and      !
   ! update the pools.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine update_maintenance(cpatch,ico,year_o_day,tempk                               &
                                ,tissue_maintenance,storage_maintenance)
      use ed_state_vars  , only : patchtype               ! ! structure
      use pft_coms       , only : phenology               & ! intent(in)
                                , leaf_turnover_rate      & ! intent(in)
                                , root_turnover_rate      & ! intent(in)
                                , bark_turnover_rate      & ! intent(in)
                                , storage_turnover_rate   ! ! intent(in)
      use physiology_coms, only : trait_plasticity_scheme ! ! intent(in)
      use consts_coms    , only : tiny_num                ! ! intent(in)
      use ed_misc_coms   , only : storage_resp_scheme     ! ! intent(in)
      use allometry      , only : ed_balive               ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(in)   :: year_o_day
      real           , intent(in)   :: tempk
      real           , intent(out)  :: tissue_maintenance
      real           , intent(out)  :: storage_maintenance
      !----- Local variables. -------------------------------------------------------------!
      integer                       :: ipft
      real                          :: maintenance_temp_dep
      real                          :: fp_turnover
      real                          :: stor_mco_o_balive
      logical                       :: dynamic_trait
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!




      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !                                 TISSUE MAINTENANCE                                 !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In case this simulation allows trait plasticity, compute the correction factor !
      ! for turnover.                                                                      !
      !------------------------------------------------------------------------------------!
      dynamic_trait =       ( trait_plasticity_scheme == 1  .or. phenology(ipft) == 3)     &
                      .and. ( leaf_turnover_rate(ipft) > 0. )
      if (dynamic_trait) then
         !----- Trait plasticity, or light-controlled phenology. --------------------------!
         fp_turnover = 12. / (cpatch%llspan(ico) * leaf_turnover_rate(ipft))
         !---------------------------------------------------------------------------------!
      else
         !----- No dynamic trait, set fp_turnover to 1 and use default values. ------------!
         fp_turnover = 1.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the maintenance costs.  This will depend on the type of phenology that   !
      ! the PFT has.   The year_o_day term applied converts the maintenance rates to       !
      ! [kgC/plant/day].                                                                   !
      !------------------------------------------------------------------------------------!
      cpatch%leaf_maintenance (ico) = leaf_turnover_rate (ipft) * cpatch%bleaf (ico)       &
                                    * year_o_day * fp_turnover
      cpatch%root_maintenance (ico) = root_turnover_rate (ipft) * cpatch%broot (ico)       &
                                    * year_o_day * fp_turnover
      cpatch%barka_maintenance(ico) = bark_turnover_rate(ipft)  * cpatch%bbarka(ico)       &
                                    * year_o_day * fp_turnover
      cpatch%barkb_maintenance(ico) = bark_turnover_rate(ipft)  * cpatch%bbarkb(ico)       &
                                    * year_o_day * fp_turnover

      select case (phenology(ipft))
      case (0)
         !---------------------------------------------------------------------------------!
         !     Evergreens, like pines.  The turnover rates will be adjusted by a function  !
         ! of temperature, which approaches 0 as the temperature goes down.                !
         !---------------------------------------------------------------------------------!
         !------ Find a temperature dependence adjustment. --------------------------------!
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
         !----- Scale maintenance by biomass and apply the temperature correction. --------!
         cpatch%leaf_maintenance (ico) = cpatch%leaf_maintenance (ico)*maintenance_temp_dep
         cpatch%root_maintenance (ico) = cpatch%root_maintenance (ico)*maintenance_temp_dep
         cpatch%barka_maintenance(ico) = cpatch%barka_maintenance(ico)*maintenance_temp_dep
         cpatch%barkb_maintenance(ico) = cpatch%barkb_maintenance(ico)*maintenance_temp_dep
         !---------------------------------------------------------------------------------!
      case default
         !----- Do nothing. ---------------------------------------------------------------!
         continue
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------ Accumulate total investment in tissue maintenance. --------------------------!
      tissue_maintenance  = cpatch%leaf_maintenance (ico) + cpatch%root_maintenance (ico)  &
                          + cpatch%barka_maintenance(ico) + cpatch%barkb_maintenance(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decrement tissue biomass to account for maintenance costs.                      !
      !------------------------------------------------------------------------------------!
      cpatch%bleaf (ico) = cpatch%bleaf (ico) - cpatch%leaf_maintenance (ico)
      cpatch%broot (ico) = cpatch%broot (ico) - cpatch%root_maintenance (ico)
      cpatch%bbarka(ico) = cpatch%bbarka(ico) - cpatch%barka_maintenance(ico)
      cpatch%bbarkb(ico) = cpatch%bbarkb(ico) - cpatch%barkb_maintenance(ico)
      cpatch%balive(ico) = ed_balive(cpatch,ico)
      !------------------------------------------------------------------------------------!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !                                STORAGE MAINTENANCE                                 !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     The commented line is an experimental and arbitrary test, borrowed from        !
      ! maintainence temperature dependency for needleleaves. [[MCD]]                      !
      !------------------------------------------------------------------------------------!
      ! maintenance_temp_dep = 1.0                                                         &
      !                      / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
      maintenance_temp_dep = 1.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the storage maintenance cost (turnover_rate), in kgC/plant/year, akin to   !
      ! leaf and fine root maintenance.  We will apply the storage respiration due to      !
      !  maintenance in subroutine apply_maintenance.                                      !
      !------------------------------------------------------------------------------------!
      storage_maintenance = cpatch%bstorage(ico) * storage_turnover_rate(ipft)             &
                          * year_o_day * maintenance_temp_dep
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Decrement storage biomass to account for maintenance costs.                     !
      !------------------------------------------------------------------------------------!
      cpatch%bstorage(ico) = cpatch%bstorage(ico) - storage_maintenance
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the contribution from each tissue to total storage respiration.           !
      !------------------------------------------------------------------------------------!
      select case (storage_resp_scheme)
      case (0)
         !---- Partition amongst tissues not required, send everything to AG sapwood. -----!
         cpatch%leaf_storage_resp (ico) = 0.0
         cpatch%root_storage_resp (ico) = 0.0
         cpatch%sapa_storage_resp (ico) = storage_maintenance
         cpatch%sapb_storage_resp (ico) = 0.0
         cpatch%barka_storage_resp(ico) = 0.0
         cpatch%barkb_storage_resp(ico) = 0.0
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     Partition storage respiration amongst all living tissues.   It is unlikely, !
         ! but in case living biomass is zero and storage turnover is not zero, we put all !
         ! the heterotrophic respiration in AG sapwood (like storage_resp_scheme = 0).     !
         !---------------------------------------------------------------------------------!
         if (cpatch%balive(ico) >= tiny_num) then
            stor_mco_o_balive              = storage_maintenance / cpatch%balive   (ico)
            cpatch%leaf_storage_resp (ico) = stor_mco_o_balive   * cpatch%bleaf    (ico)
            cpatch%root_storage_resp (ico) = stor_mco_o_balive   * cpatch%broot    (ico)
            cpatch%sapa_storage_resp (ico) = stor_mco_o_balive   * cpatch%bsapwooda(ico)
            cpatch%sapb_storage_resp (ico) = stor_mco_o_balive   * cpatch%bsapwoodb(ico)
            cpatch%barka_storage_resp(ico) = stor_mco_o_balive   * cpatch%bbarka   (ico)
            cpatch%barkb_storage_resp(ico) = stor_mco_o_balive   * cpatch%bbarkb   (ico)
         else
            !---- Partition amongst tissues not possible, send everything to AG sapwood. --!
            cpatch%leaf_storage_resp (ico) = 0.0
            cpatch%root_storage_resp (ico) = 0.0
            cpatch%sapa_storage_resp (ico) = storage_maintenance
            cpatch%sapb_storage_resp (ico) = 0.0
            cpatch%barka_storage_resp(ico) = 0.0
            cpatch%barkb_storage_resp(ico) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      return
   end subroutine update_maintenance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the daily metabolic NPP.  We also find the potential      !
   ! values: pot (no nitrogen limitation); lightmax (no light limitation); moistmax (no    !
   ! soil moisture limitation); mlmax (neither light nor moisture limitation).             !
   !---------------------------------------------------------------------------------------!
   subroutine get_metabolic_npp(cpatch,ico,metnpp_actual,metnpp_pot,metnpp_lightmax        &
                               ,metnpp_moistmax,metnpp_mlmax)
      use ed_state_vars  , only : patchtype          ! ! structure
      use consts_coms    , only : umol_2_kgC         & ! intent(in)
                                , day_sec            & ! intent(in)
                                , tiny_num           ! ! intent(in)
      use ed_max_dims    , only : n_pft              ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ico
      real                     , intent(out) :: metnpp_actual
      real                     , intent(out) :: metnpp_pot
      real                     , intent(out) :: metnpp_lightmax
      real                     , intent(out) :: metnpp_moistmax
      real                     , intent(out) :: metnpp_mlmax
      !----- Local variables. -------------------------------------------------------------!
      real                                   :: f_unitconv
      real                                   :: metab_resp
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the daily "metabolic" NPP, i.e. GPP minus leaf respiration minus fine    !
      ! root respiration.                                                                  !
      !------------------------------------------------------------------------------------!
      if (cpatch%nplant(ico) >= tiny_num) then
         !----- Find conversion from umol/m2/s to kgC/plant/day. --------------------------!
         f_unitconv = umol_2_kgC * day_sec / cpatch%nplant(ico)
         !---------------------------------------------------------------------------------!

         !----- Find the metabolic respiration, assumed to be the same for all NPP. -------!
         metab_resp    = cpatch%today_leaf_resp(ico)                                       & 
                       + cpatch%today_root_resp(ico)                                       &
                       + cpatch%today_stem_resp(ico)                                       
         !---------------------------------------------------------------------------------!


         !----- Find the metabolic NPP. ---------------------------------------------------!
         metnpp_actual   = ( cpatch%today_gpp         (ico) - metab_resp ) * f_unitconv
         metnpp_pot      = ( cpatch%today_gpp_pot     (ico) - metab_resp ) * f_unitconv
         metnpp_lightmax = ( cpatch%today_gpp_lightmax(ico) - metab_resp ) * f_unitconv
         metnpp_moistmax = ( cpatch%today_gpp_moistmax(ico) - metab_resp ) * f_unitconv
         metnpp_mlmax    = ( cpatch%today_gpp_mlmax   (ico) - metab_resp ) * f_unitconv
         !---------------------------------------------------------------------------------!
      else
         !----- nplant is somehow zero (this shouldn't happen). ---------------------------!
         metnpp_actual   = 0.0
         metnpp_pot      = 0.0
         metnpp_lightmax = 0.0
         metnpp_moistmax = 0.0
         metnpp_mlmax    = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine get_metabolic_npp
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     Find the growth respiration, which will be released to the canopy air space as    !
   ! the day progresses (although it refers to growth based on the previous day's NPP.     !
   !---------------------------------------------------------------------------------------!
   subroutine get_growth_respiration(cpatch,ico,metnpp_actual,metnpp_pot,metnpp_lightmax   &
                                    ,metnpp_moistmax,metnpp_mlmax,growresp_actual          &
                                    ,growresp_pot,growresp_lightmax,growresp_moistmax      &
                                    ,growresp_mlmax,pat_growth_resp)
      use ed_state_vars, only : patchtype             ! ! structure
      use consts_coms  , only : tiny_num              ! ! intent(in)
      use ed_misc_coms , only : growth_resp_scheme    ! ! intent(in)
      use pft_coms     , only : growth_resp_factor    ! ! intent(in)
       implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: metnpp_actual
      real           , intent(in)    :: metnpp_pot
      real           , intent(in)    :: metnpp_lightmax
      real           , intent(in)    :: metnpp_moistmax
      real           , intent(in)    :: metnpp_mlmax
      real           , intent(out)   :: growresp_actual
      real           , intent(out)   :: growresp_pot
      real           , intent(out)   :: growresp_lightmax
      real           , intent(out)   :: growresp_moistmax
      real           , intent(out)   :: growresp_mlmax
      real           , intent(inout) :: pat_growth_resp
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: grow_resp_o_balive
      !------------------------------------------------------------------------------------!



      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Total growth respiration (all of them, actual and the potential ones).  In    !
      ! case the cohort had a bad day and ended up with negative metabolic NPP, growth     !
      ! respiration becomes zero.                                                          !
      !------------------------------------------------------------------------------------!
      growresp_actual   = max(0.,metnpp_actual   * growth_resp_factor(ipft))
      growresp_pot      = max(0.,metnpp_pot      * growth_resp_factor(ipft))
      growresp_lightmax = max(0.,metnpp_lightmax * growth_resp_factor(ipft))
      growresp_moistmax = max(0.,metnpp_moistmax * growth_resp_factor(ipft))
      growresp_mlmax    = max(0.,metnpp_mlmax    * growth_resp_factor(ipft))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check whether actual growth respiration is positive.  In case it isn't, don't  !
      ! bother checking the partition amongst tissues.                                     !
      !------------------------------------------------------------------------------------!
      if (growresp_actual < tiny_num) then
         growresp_actual               = 0.0
         cpatch%leaf_growth_resp (ico) = 0.0
         cpatch%root_growth_resp (ico) = 0.0
         cpatch%sapa_growth_resp (ico) = 0.0
         cpatch%sapb_growth_resp (ico) = 0.0
         cpatch%barka_growth_resp(ico) = 0.0
         cpatch%barkb_growth_resp(ico) = 0.0
         return
      end if
      !------------------------------------------------------------------------------------!


      !----- Update the committed growth respiration at patch level. ----------------------!
      pat_growth_resp = pat_growth_resp + cpatch%nplant(ico) * growresp_actual
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the contribution from each tissue to total growth respiration.            !
      !------------------------------------------------------------------------------------!
      select case (growth_resp_scheme)
      case (0)
         !---- Partition amongst tissues not required, send everything to AG sapwood. -----!
         cpatch%leaf_growth_resp (ico) = 0.0
         cpatch%root_growth_resp (ico) = 0.0
         cpatch%sapa_growth_resp (ico) = growresp_actual
         cpatch%sapb_growth_resp (ico) = 0.0
         cpatch%barka_growth_resp(ico) = 0.0
         cpatch%barkb_growth_resp(ico) = 0.0
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     Partition growth respiration amongst all living tissues.   It is pretty     !
         ! much impossible, but in case the plant had no living tissues and still managed  !
         ! to assimilate carbon, we put all growth respiration to AG sapwood, like         !
         ! growth_resp_scheme 0.  Otherwise, allocation growth respiration proportionally  !
         ! to each tissue.                                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%balive(ico) >= tiny_num) then
            grow_resp_o_balive            = growresp_actual    / cpatch%balive   (ico)
            cpatch%leaf_growth_resp(ico)  = grow_resp_o_balive * cpatch%bleaf    (ico)
            cpatch%root_growth_resp(ico)  = grow_resp_o_balive * cpatch%broot    (ico)
            cpatch%sapa_growth_resp(ico)  = grow_resp_o_balive * cpatch%bsapwooda(ico)
            cpatch%sapb_growth_resp(ico)  = grow_resp_o_balive * cpatch%bsapwoodb(ico)
            cpatch%barka_growth_resp(ico) = grow_resp_o_balive * cpatch%bbarka   (ico)
            cpatch%barkb_growth_resp(ico) = grow_resp_o_balive * cpatch%bbarkb   (ico)
         else
            !---- Partition amongst tissues not possible, send everything to AG sapwood. --!
            cpatch%leaf_growth_resp (ico) = 0.0
            cpatch%root_growth_resp (ico) = 0.0
            cpatch%sapa_growth_resp (ico) = growresp_actual
            cpatch%sapb_growth_resp (ico) = 0.0
            cpatch%barka_growth_resp(ico) = 0.0
            cpatch%barkb_growth_resp(ico) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine get_growth_respiration
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the potential and total carbon balances for density-      !
   ! -dependent mortality rates.                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine update_carbon_balances(cpatch,ipa,ico,metnpp_actual,metnpp_pot               &
                                    ,metnpp_lightmax,metnpp_moistmax,metnpp_mlmax          &
                                    ,growresp_actual,growresp_pot,growresp_lightmax        &
                                    ,growresp_moistmax,growresp_mlmax,tissue_maintenance   &
                                    ,storage_maintenance,npp_actual,npp_pot )
      use ed_state_vars  , only : patchtype          ! ! structure
      use physiology_coms, only : iddmort_scheme     ! ! intent(in)
      use ed_misc_coms   , only : current_time       ! ! intent(in)
      use ed_max_dims    , only : n_pft              ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: metnpp_actual
      real                     , intent(in)  :: metnpp_pot
      real                     , intent(in)  :: metnpp_lightmax
      real                     , intent(in)  :: metnpp_moistmax
      real                     , intent(in)  :: metnpp_mlmax
      real                     , intent(in)  :: growresp_actual
      real                     , intent(in)  :: growresp_pot
      real                     , intent(in)  :: growresp_lightmax
      real                     , intent(in)  :: growresp_moistmax
      real                     , intent(in)  :: growresp_mlmax
      real                     , intent(in)  :: tissue_maintenance
      real                     , intent(in)  :: storage_maintenance
      real                     , intent(out) :: npp_actual
      real                     , intent(out) :: npp_pot
      !----- Local variables. -------------------------------------------------------------!
      integer                                :: ipft
      real                                   :: total_maintenance
      real                                   :: npp_lightmax
      real                                   :: npp_moistmax
      real                                   :: npp_mlmax
      !----- Local constants. -------------------------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !------------------------------------------------------------------------------------!


      !----- Alias for PFT type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the net primary productivity (metabolic NPP minus growth respiration),   !
      ! the actual and all the potential ones.                                             !
      !                                                                                    !
      ! MLO - This is a difference from the mainline, because the daily NPP uses the       !
      !       growth rate calculated from the previous day metabolic NPP, not the one from !
      !       two days ago.  The reason for this change is to make the carbon balance      !
      !       vectors consistent (so they are all based on the same metabolic              !
      !       respiration).                                                                !
      !------------------------------------------------------------------------------------!
      npp_actual   = metnpp_actual   - growresp_actual
      npp_pot      = metnpp_pot      - growresp_pot
      npp_lightmax = metnpp_lightmax - growresp_lightmax
      npp_moistmax = metnpp_moistmax - growresp_moistmax
      npp_mlmax    = metnpp_mlmax    - growresp_mlmax
      !------------------------------------------------------------------------------------!



      !----- Set total maintenance for DD mortality according to the selected method. -----!
      select case (iddmort_scheme)
      case (0) ! Storage is not accounted.
         total_maintenance = tissue_maintenance
      case (1) ! Storage is accounted.
         total_maintenance = tissue_maintenance + storage_maintenance
      end select
      !------------------------------------------------------------------------------------!



      !----- Update carbon balance variables (used for density-dependent mortality). ------!
      cpatch%cb         (13,ico) = cpatch%cb         (13,ico)                              &
                                 + npp_actual   - total_maintenance
      cpatch%cb_lightmax(13,ico) = cpatch%cb_lightmax(13,ico)                              &
                                 + npp_lightmax - total_maintenance
      cpatch%cb_moistmax(13,ico) = cpatch%cb_moistmax(13,ico)                              &
                                 + npp_moistmax - total_maintenance
      cpatch%cb_mlmax   (13,ico) = cpatch%cb_mlmax   (13,ico)                              &
                                 + npp_mlmax    - total_maintenance
      !------------------------------------------------------------------------------------!



      if (print_debug) then

         if (first_time(ipft)) then
            first_time(ipft) = .false.
            write (unit=30+ipft,fmt='(a10,29(1x,a18))')                                    &
               '      TIME','             PATCH','            COHORT','            NPLANT' &
                           ,'         TODAY_NPP','  LEAF_GROWTH_RESP','  ROOT_GROWTH_RESP' &
                           ,'  SAPA_GROWTH_RESP','  SAPB_GROWTH_RESP',' BARKA_GROWTH_RESP' &
                           ,' BARKB_GROWTH_RESP','         TODAY_GPP','TODAY_GPP_LIGHTMAX' &
                           ,'TODAY_GPP_MOISTMAX','   TODAY_GPP_MLMAX','   TODAY_LEAF_RESP' &
                           ,'   TODAY_ROOT_RESP','   TODAY_STEM_RESP','TODAY_NPP_LIGHTMAX' &
                           ,'TODAY_NPP_MOISTMAX','   TODAY_NPP_MLMAX','                CB' &
                           ,'       CB_LIGHTMAX','       CB_MOISTMAX','          CB_MLMAX' &
                           ,'  LEAF_MAINTENANCE','  ROOT_MAINTENANCE',' BARKA_MAINTENANCE' &
                           ,' BARKB_MAINTENANCE', ' STOR_MAINTENANCE'
         end if

         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i18),28(1x,es18.5))')               &
              current_time%month,'/',current_time%date,'/',current_time%year               &
             ,ipa,ico,cpatch%nplant(ico),npp_actual,cpatch%leaf_growth_resp(ico)           &
             ,cpatch%root_growth_resp(ico),cpatch%sapa_growth_resp(ico)                    &
             ,cpatch%sapb_growth_resp(ico),cpatch%barka_growth_resp(ico)                   &
             ,cpatch%barkb_growth_resp(ico),cpatch%today_gpp(ico)                          &
             ,cpatch%today_gpp_lightmax(ico),cpatch%today_gpp_moistmax(ico)                &
             ,cpatch%today_gpp_mlmax(ico),cpatch%today_leaf_resp(ico)                      &
             ,cpatch%today_root_resp(ico),cpatch%today_stem_resp(ico)                      &       
             ,npp_lightmax,npp_moistmax,npp_mlmax                                          &
             ,cpatch%cb(13,ico),cpatch%cb_lightmax(13,ico),cpatch%cb_moistmax(13,ico)      &
             ,cpatch%cb_mlmax(13,ico),cpatch%leaf_maintenance(ico)                         &
             ,cpatch%root_maintenance(ico),cpatch%barka_maintenance(ico)                   &
             ,cpatch%barkb_maintenance(ico),storage_maintenance
      end if

      return
   end subroutine update_carbon_balances
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: GET_C_XFERS
   !
   !> \brief   Calculates plant-internal C transfers for growth and maintainance.
   !> \details Uses phenology, allometry, carbon balance, and current C pool sizes to
   !>          determine (signed) transfers to leaf, root, sapwooda, sapwoodb, bark,
   !>          and storage.
   !> \warning The order of the operations here affect the C/N budgets, so don't
   !>          change it unless you really know what you are doing.
   !---------------------------------------------------------------------------------------!
   subroutine get_c_xfers(csite,ipa,ico,npp_actual,green_leaf_factor,gr_tfact0             &
                         ,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb,tr_bbarka,tr_bbarkb  &
                         ,tr_bstorage,carbon_debt,flushing,balive_aim,carbon_miss          &
                         ,xfer_case)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , qbark        & ! intent(in)
                               , agf_bs       & ! intent(in)
                               , is_grass     & ! intent(in)
                               , balive_crit  ! ! intent(in)
      use allometry     , only : size2bl      & ! function
                               , ba2h         ! ! function
      use consts_coms   , only : tiny_num     ! ! intent(in)
      use ed_misc_coms  , only : igrass       & ! intent(in)
                               , iallom       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite             !< Current Site
      integer        , intent(in)    :: ipa               !< Loop-Current Patch
      integer        , intent(in)    :: ico               !< Loop-Current Cohort
      real           , intent(in)    :: npp_actual        !< Plant net carbon uptake
      real           , intent(in)    :: green_leaf_factor !< Cohort leaf-age param.
      real           , intent(in)    :: gr_tfact0         !< Corr. factor for steady growth
      real           , intent(out)   :: tr_bleaf          !< Transfer to leaf C pool
      real           , intent(out)   :: tr_broot          !< Transfer to root C pool
      real           , intent(out)   :: tr_bsapwooda      !< Transfer to sapwooda C pool
      real           , intent(out)   :: tr_bsapwoodb      !< Transfer to sapwoodb C pool
      real           , intent(out)   :: tr_bbarka         !< Transfer to barka C pool
      real           , intent(out)   :: tr_bbarkb         !< Transfer to barkb C pool
      real           , intent(out)   :: tr_bstorage       !< Transfer to storage C pool
      real           , intent(out)   :: carbon_debt       !< Net cohort carbon uptake
      logical        , intent(out)   :: flushing          !< Flag for leaf flush
      real           , intent(out)   :: balive_aim        !< Desired cohort balive value
      real           , intent(inout) :: carbon_miss       !< Carbon from unaccounted source
      integer        , intent(out)   :: xfer_case         !< Transfer case (for debugging)
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bleaf_aim
      real                           :: broot_aim
      real                           :: bsapwooda_aim
      real                           :: bsapwoodb_aim
      real                           :: bbarka_aim
      real                           :: bbarkb_aim
      real                           :: balive_max
      real                           :: bleaf_max
      real                           :: bloss_max
      real                           :: height_aim
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwooda
      real                           :: delta_bsapwoodb
      real                           :: delta_bbarka
      real                           :: delta_bbarkb
      real                           :: delta_btotal
      real                           :: available_carbon
      real                           :: gtf_bleaf
      real                           :: gtf_broot
      real                           :: gtf_bbarka
      real                           :: gtf_bbarkb
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bbarka
      real                           :: f_bbarkb
      real                           :: f_bsapwooda
      real                           :: f_bsapwoodb
      logical                        :: time_to_flush
      integer                        :: phen_stat_in
      !logical          , parameter   :: printout = .false.
      !character(len=11), parameter   :: fracfile = 'cballoc.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      !logical          , save        :: first_time = .true.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! This could have been garbage collected out of the code or updated following the    !
      ! modularization of growth_balive.f90, but is being left in place as a template in   !
      ! case it should be maintained.                                                      !
      !----- First time, and the user wants to print the output.  Make a header. ----------!
      !if (first_time) then
      !   if (printout) then
      !      open (unit=66,file=fracfile,status='replace',action='write')
      !      write (unit=66,fmt='(24(a,1x))')                                              &
      !        ,'        YEAR','       MONTH','         DAY','         PFT','   PHENOLOGY' &
      !        ,'PHEN_STAT_IN','PHN_STAT_OUT','  FLUSH_TIME',' AVAILABLE_C','      ELONGF' &
      !        ,'  GREEN_LEAF','    ON_ALLOM',' DELTA_BLEAF',' DELTA_BROOT','   DELTA_BSA' &
      !        ,'   DELTA_BSB','DELTA_BBARKA','DELTA_BBARKB','    TR_BLEAF','    TR_BROOT' &
      !        ,'    TR_BSAPA','    TR_BSAPB','    TR_BARKA','    TR_BARKB'
      !      close (unit=66,status='keep')
      !   end if
      !   first_time = .false.
      !end if
      !------------------------------------------------------------------------------------!

      tr_bleaf     = 0.0
      tr_broot     = 0.0
      tr_bsapwooda = 0.0
      tr_bsapwoodb = 0.0
      tr_bbarka    = 0.0
      tr_bbarkb    = 0.0
      tr_bstorage  = 0.0

      cpatch => csite%patch(ipa)

      ipft = cpatch%pft(ico)
      phen_stat_in = cpatch%phenology_status(ico)
      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! npp_actual is negative, but the sum of npp_actual and bstorage is positive.  In    !
      ! this case, we allow plants to grow leaves.                                         !
      !------------------------------------------------------------------------------------!
      available_carbon = cpatch%bstorage(ico) + npp_actual
      time_to_flush    = npp_actual >= tiny_num .or.                                       &
                         ( available_carbon >= tiny_num .and.                              &
                           cpatch%phenology_status(ico) == 1 )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide maximum leaf based on life form.  Grasses may grow every day. 
      !------------------------------------------------------------------------------------!
      if (igrass == 1 .and. is_grass(ipft) .and. time_to_flush) then
         !---------------------------------------------------------------------------------!
         !    New grasses and they are in a vegetative growth phase, increase bleaf to                       !
         !---------------------------------------------------------------------------------!
         balive_max = min(balive_crit(ipft),cpatch%balive(ico) + available_carbon)
         height_aim = ba2h(balive_max,ipft)
         bleaf_max  = balive_max / (1.0 + q(ipft) + (qsw(ipft) + qbark(ipft)) * height_aim)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Maximum bleaf that the allometric relationship would allow.  If the plant   !
         ! is drought stressed (elongf<1), we down-regulate allocation to balive.          !
         !---------------------------------------------------------------------------------!
         bleaf_max     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico),ipft)
         height_aim    = cpatch%hite(ico)
         !---------------------------------------------------------------------------------!
      end if
      bleaf_aim     = bleaf_max * green_leaf_factor * cpatch%elongf(ico)
      broot_aim     = bleaf_aim * q    (ipft)
      bsapwooda_aim = bleaf_aim * qsw  (ipft) * height_aim * agf_bs(ipft)
      bsapwoodb_aim = bleaf_aim * qsw  (ipft) * height_aim * (1. - agf_bs(ipft))
      bbarka_aim    = bleaf_aim * qbark(ipft) * height_aim * agf_bs(ipft)
      bbarkb_aim    = bleaf_aim * qbark(ipft) * height_aim * (1. - agf_bs(ipft))
      balive_aim    = bleaf_aim     + broot_aim                                            &
                    + bsapwooda_aim + bsapwoodb_aim                                        &
                    + bbarka_aim    + bbarkb_aim
      !------------------------------------------------------------------------------------!
      


      !------------------------------------------------------------------------------------!
      !      Check whether to increase living tissue biomass or not.                       !
      !------------------------------------------------------------------------------------!
      flushing = .false.
      if (time_to_flush) then
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            !     There are leaves, we are not actively dropping leaves and we're off      !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!
            flushing = .true.
            !------------------------------------------------------------------------------!



            !---- Amount that bleaf, broot, and bsapwood are off allometry. ---------------!
            delta_bleaf     = max (0.0, bleaf_aim     - cpatch%bleaf    (ico))
            delta_broot     = max (0.0, broot_aim     - cpatch%broot    (ico))
            delta_bsapwooda = max (0.0, bsapwooda_aim - cpatch%bsapwooda(ico))
            delta_bsapwoodb = max (0.0, bsapwoodb_aim - cpatch%bsapwoodb(ico))
            delta_bbarka    = max (0.0, bbarka_aim    - cpatch%bbarka   (ico))
            delta_bbarkb    = max (0.0, bbarkb_aim    - cpatch%bbarkb   (ico))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            ! MLO: Find correction factors for growth.  Rationale: the original approach   !
            !      will try to fix allometry in the first day of the month, to catch up    !
            !      the increase in heartwood biomass (bdead).  This could lead to large    !
            !      CO2 fluxes due to growth respiration in one day of the month.  This is  !
            !      less of an issue in the old allometry because sapwood is so under-      !
            !      estimated, but with the correct ratios this flux can be large.          !
            !      To make the increments more steady throughout the month, we include     !
            !      a scaling factor that makes the increments the same throughout the      !
            !      month if available carbon is not limiting.  Currently this is only      !
            !      applied to iallom 3 so it doesn't break other people's code, although   !
            !      this may change in the future.                                          !
            !                                                                              !
            !      - gr_tfact0 splits the increment to bring back to allometry equal       !
            !        amounts every day if the tissue does not have maintenance costs.      !
            !      - gtf_xxxx corrects gr_tfact0 to account for the maintenance, so the    !
            !        net growth (i.e. increment from the value of btissue before           !
            !        maintenance was applied) is the same every day.                       !
            !------------------------------------------------------------------------------!
            select case (iallom)
            case (3,4,5)
               if (.not. (is_grass(ipft) .and. igrass == 1) ) then
                  if (delta_bleaf >= tiny_num) then
                     gtf_bleaf = ( cpatch%leaf_maintenance(ico)                            &
                                 + gr_tfact0 * (delta_bleaf-cpatch%leaf_maintenance(ico))) &
                               / delta_bleaf
                  else
                     gtf_bleaf = gr_tfact0
                  end if
                  if (delta_broot >= tiny_num) then
                     gtf_broot = ( cpatch%root_maintenance(ico)                            &
                                 + gr_tfact0 * (delta_broot-cpatch%root_maintenance(ico))) &
                               / delta_broot
                  else
                     gtf_broot = gr_tfact0
                  end if
                  if (delta_bbarka >= tiny_num) then
                     gtf_bbarka = ( cpatch%barka_maintenance(ico)                          &
                                  + gr_tfact0                                              &
                                  * (delta_bbarka - cpatch%barka_maintenance(ico)) )       &
                               / delta_bbarka
                  else
                     gtf_bbarka = gr_tfact0
                  end if
                  if (delta_bbarkb >= tiny_num) then
                     gtf_bbarkb = ( cpatch%barkb_maintenance(ico)                          &
                                  + gr_tfact0                                              &
                                  * (delta_bbarkb - cpatch%barkb_maintenance(ico)) )       &
                               / delta_bbarkb
                  else
                     gtf_bbarkb = gr_tfact0
                  end if
                  !----- Correct deltas based on the time of the month and turnover. ------!
                  delta_bleaf     = delta_bleaf     * gtf_bleaf
                  delta_broot     = delta_broot     * gtf_broot
                  delta_bsapwooda = delta_bsapwooda * gr_tfact0 ! sapwood turnover is zero.
                  delta_bsapwoodb = delta_bsapwoodb * gr_tfact0 ! sapwood turnover is zero.
                  delta_bbarka    = delta_bbarka    * gtf_bbarka
                  delta_bbarkb    = delta_bbarkb    * gtf_bbarkb
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!


            !---- Total sought transfer. --------------------------------------------------!
            delta_btotal    = delta_bleaf     + delta_broot     + delta_bsapwooda          &
                            + delta_bsapwoodb + delta_bbarka    + delta_bbarkb
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     In case the available carbon is less than what we need to get back to    !
            ! allometry, grow pools in proportion to demand.  Otherwise, put the excess    !
            ! carbon into bstorage.                                                        !
            !------------------------------------------------------------------------------!
            if (delta_btotal >= tiny_num) then
               f_total      = min(1.0, available_carbon / delta_btotal)
               tr_bleaf     = delta_bleaf      * f_total
               tr_broot     = delta_broot      * f_total
               tr_bsapwooda = delta_bsapwooda  * f_total
               tr_bsapwoodb = delta_bsapwoodb  * f_total
               tr_bbarka    = delta_bbarka     * f_total
               tr_bbarkb    = delta_bbarkb     * f_total
               xfer_case    = 2
            else
               tr_bleaf     = 0.
               tr_broot     = 0.
               tr_bsapwooda = 0.
               tr_bsapwoodb = 0.
               tr_bbarka    = 0.
               tr_bbarkb    = 0.
               xfer_case    = 1
            end if
            !------------------------------------------------------------------------------!


            !----- Change in storage is NPP minus transfer to other tissues. --------------!
            tr_bstorage = npp_actual   - tr_bleaf     - tr_broot     - tr_bsapwooda        &
                        - tr_bsapwoodb - tr_bbarka    - tr_bbarkb
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            tr_bstorage  = npp_actual
            xfer_case    = 0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         !   Carbon balance is negative, decide the source of carbon based on the          !
         ! phenology status.  If plants were already dropping leaves, then we don't take   !
         ! the carbon from storage unless there is no leaf or root biomass left.  If       !
         ! plants should be growing but they aren't, then we burn the storage first, and   !
         ! if the situation persists, then plants start destroying their living tissues.   !
         !---------------------------------------------------------------------------------!
         carbon_debt = -npp_actual
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            !    Plants should be growing or at their maximum, first we try to take all    !
            ! the carbon needed from storage.                                              !
            !------------------------------------------------------------------------------!
            if (cpatch%bstorage(ico) > carbon_debt) then
               !------ Storage loss will make up the carbon debt. -------------------------!
               tr_bstorage = -1.0 * carbon_debt
               carbon_debt =  0.0
               xfer_case   = -1
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Not enough carbon in storage.  Take everything then start destroying  !
               ! tissues.                                                                  !
               !---------------------------------------------------------------------------!
               carbon_debt = carbon_debt - cpatch%bstorage(ico)
               tr_bstorage = -1.0*cpatch%bstorage(ico)

               !---------------------------------------------------------------------------!
               !     Find total biomass that can be lost.  We take an amount proportional  !
               ! to the current biomass of each the active pools.                          !
               !---------------------------------------------------------------------------!
               bloss_max   = cpatch%bleaf (ico) + cpatch%broot (ico)                       &
                           + cpatch%bbarka(ico) + cpatch%bbarkb(ico)
               f_bleaf     = cpatch%bleaf    (ico) / bloss_max
               f_broot     = cpatch%broot    (ico) / bloss_max
               f_bbarka    = cpatch%bbarka   (ico) / bloss_max
               f_bbarkb    = cpatch%bbarkb   (ico) / bloss_max

               if (bloss_max > carbon_debt) then
                  !----- Remove biomass accordingly. --------------------------------------!
                  tr_bleaf     = -1.0 * carbon_debt * f_bleaf
                  tr_broot     = -1.0 * carbon_debt * f_broot
                  tr_bbarka    = -1.0 * carbon_debt * f_bbarka
                  tr_bbarkb    = -1.0 * carbon_debt * f_bbarkb
                  carbon_debt  = 0.0
                  xfer_case    = -2
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort had extremely large negative NPP, storage and active   !
                  ! tissues were not sufficient to close the negative carbon balance.      !
                  ! Eliminate all biomass from leaves, fine roots, and bark.               !
                  !------------------------------------------------------------------------!
                  carbon_debt  = carbon_debt - bloss_max
                  tr_bleaf     = -1.0 * cpatch%bleaf    (ico)
                  tr_broot     = -1.0 * cpatch%broot    (ico)
                  tr_bbarka    = -1.0 * cpatch%bbarka   (ico)
                  tr_bbarkb    = -1.0 * cpatch%bbarkb   (ico)
                  xfer_case    = -90
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     This situation also poses a dead end for this cohort, as it has    !
                  ! no way out (no storage to build active tissues, and no remaining       !
                  ! active tissues).  Flag the cohort as inviable and terminate it.        !
                  !------------------------------------------------------------------------!
                  cpatch%is_viable(ico) = .false.
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         case (-1,-2)
            !------------------------------------------------------------------------------!
            !      Plants were already shedding leaves.  We swap the order here and remove !
            ! living tissues first, and only if there is nothing left that we remove       !
            ! storage.  Bark is not 100% living tissue, but for simplicity we don't        !
            ! separate inner bark from outer bark.                                         !
            !------------------------------------------------------------------------------!
            bloss_max   = cpatch%bleaf (ico) + cpatch%broot (ico)                          &
                        + cpatch%bbarka(ico) + cpatch%bbarkb(ico)
            if (bloss_max >= tiny_num) then
               f_bleaf     = cpatch%bleaf    (ico) / bloss_max
               f_broot     = cpatch%broot    (ico) / bloss_max
               f_bbarka    = cpatch%bbarka   (ico) / bloss_max
               f_bbarkb    = cpatch%bbarkb   (ico) / bloss_max
            else
               f_bleaf     = 0.
               f_broot     = 0.
               f_bbarka    = 0.
               f_bbarkb    = 0.
            end if

            if (bloss_max > carbon_debt) then
               !----- Remove biomass accordingly. -----------------------------------------!
               tr_bleaf     = -1.0 * carbon_debt * f_bleaf
               tr_broot     = -1.0 * carbon_debt * f_broot
               tr_bbarka    = -1.0 * carbon_debt * f_bbarka
               tr_bbarkb    = -1.0 * carbon_debt * f_bbarkb
               carbon_debt  = 0.0
               xfer_case    = -11
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Not enough biomass, remove everything.                                !
               !---------------------------------------------------------------------------!
               carbon_debt  = carbon_debt - bloss_max
               tr_bleaf     = -1.0 * cpatch%bleaf    (ico)
               tr_broot     = -1.0 * cpatch%broot    (ico)
               tr_bbarka    = -1.0 * cpatch%bbarka   (ico)
               tr_bbarkb    = -1.0 * cpatch%bbarkb   (ico)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     The living tissues weren't enough to meet the demand, remove what is  !
               ! still needed from the storage.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%bstorage(ico) > carbon_debt) then
                  !----- Enough carbon in storage, take all carbon needed from there. -----!
                  tr_bstorage = -1.0 * carbon_debt
                  carbon_debt =  0.0
                  xfer_case    = -12
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort had extremely large negative NPP, storage and active   !
                  ! tissues were not sufficient to close the negative carbon balance.      !
                  ! Eliminate all storage.                                                 !
                  !------------------------------------------------------------------------!
                  tr_bstorage = -1.0*cpatch%bstorage(ico)
                  carbon_debt = carbon_debt - cpatch%bstorage(ico)
                  xfer_case   = -92
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     This situation also poses a dead end for this cohort, as it has    !
                  ! no way out (no storage to build active tissues, and no remaining       !
                  ! active tissues).  Flag the cohort as inviable and terminate it.        !
                  !------------------------------------------------------------------------!
                  cpatch%is_viable(ico) = .false.
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      In case of inviable cohorts, their carbon budget may not be closed.  To    !
         ! try to avoid unaccounted carbon, we also surrender sapwood.   In case the       !
         ! budget is still not closed, report the missing carbon.                          !
         !---------------------------------------------------------------------------------!
         if ( carbon_debt >= tiny_num) then
            !---- Check maximum amount that can be removed from sapwood. ------------------!
            bloss_max = cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)
            if (bloss_max >= tiny_num) then
               f_bsapwooda    = cpatch%bsapwooda(ico) / bloss_max
               f_bsapwoodb    = cpatch%bsapwoodb(ico) / bloss_max
            else
               f_bsapwooda    = 0.0
               f_bsapwoodb    = 0.0
            end if
            !------------------------------------------------------------------------------!

            if (bloss_max > carbon_debt) then
               !----- Remove biomass accordingly. -----------------------------------------!
               tr_bsapwooda = -1.0 * carbon_debt * f_bsapwooda
               tr_bsapwoodb = -1.0 * carbon_debt * f_bsapwoodb
               carbon_debt  = 0.0
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Not enough biomass, remove everything.                                !
               !---------------------------------------------------------------------------!
               carbon_debt  = carbon_debt - bloss_max
               tr_bsapwooda = -1.0 * cpatch%bsapwooda(ico)
               tr_bsapwoodb = -1.0 * cpatch%bsapwoodb(ico)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Report the missing carbon.                                            !
               !---------------------------------------------------------------------------!
               carbon_miss = carbon_miss + carbon_debt
               xfer_case   = -99
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
!      if (printout) then
!         open (unit=66,file=fracfile,status='old',position='append',action='write')
!         write(unit=66,fmt='(7(i12,1x),1(11x,l1,1x),3(f12.6,1x),1(11x,l1,1x),8(f12.8,1x))')&
!               current_time%year,current_time%month,current_time%date,ipft,phenology(ipft) &
!              ,phen_stat_in,cpatch%phenology_status(ico),time_to_flush,available_carbon    &
!              ,cpatch%elongf(ico),green_leaf_factor,on_allometry,delta_bleaf,delta_broot   &
!              ,delta_bsapwooda,delta_bsapwoodb,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb
!         close (unit=66,status='keep')
!      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine get_c_xfers
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   subroutine apply_c_xfers(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb         &
                           ,tr_bbarka,tr_bbarkb,tr_bstorage)
      use ed_state_vars , only : patchtype  ! ! structure
      use allometry     , only : ed_balive  ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: tr_bleaf
      real           , intent(in)    :: tr_broot
      real           , intent(in)    :: tr_bsapwooda
      real           , intent(in)    :: tr_bsapwoodb
      real           , intent(in)    :: tr_bbarka
      real           , intent(in)    :: tr_bbarkb
      real           , intent(in)    :: tr_bstorage
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update the carbon pools of the living tissues.                                 !
      !------------------------------------------------------------------------------------!
      cpatch%bleaf    (ico) = cpatch%bleaf    (ico) + tr_bleaf
      cpatch%broot    (ico) = cpatch%broot    (ico) + tr_broot
      cpatch%bsapwooda(ico) = cpatch%bsapwooda(ico) + tr_bsapwooda
      cpatch%bsapwoodb(ico) = cpatch%bsapwoodb(ico) + tr_bsapwoodb
      cpatch%bbarka   (ico) = cpatch%bbarka   (ico) + tr_bbarka
      cpatch%bbarkb   (ico) = cpatch%bbarkb   (ico) + tr_bbarkb
      cpatch%balive   (ico) = ed_balive(cpatch,ico)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the amount of carbon used to recover the tissues that were off-allometry,  !
      ! take that from the carbon balance first, then use some of the storage if needed    !
      ! be.                                                                                !
      !------------------------------------------------------------------------------------!
      cpatch%bstorage(ico) = cpatch%bstorage(ico) + tr_bstorage
      !------------------------------------------------------------------------------------!

      return
   end subroutine apply_c_xfers
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   subroutine update_today_npp_vars(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb &
                                   ,tr_bbarka,tr_bbarkb,npp_actual)
      use ed_state_vars , only : patchtype    ! ! structure
      use consts_coms   , only : tiny_num     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(in)   :: tr_bleaf
      real           , intent(in)   :: tr_broot
      real           , intent(in)   :: tr_bsapwooda
      real           , intent(in)   :: tr_bsapwoodb
      real           , intent(in)   :: tr_bbarka
      real           , intent(in)   :: tr_bbarkb
      real           , intent(in)   :: npp_actual
      !----- Local variables. -------------------------------------------------------------!
      real                          :: tr_bsapwood
      real                          :: tr_bbark
      !------------------------------------------------------------------------------------!

      tr_bsapwood = tr_bsapwooda + tr_bsapwoodb
      tr_bbark    = tr_bbarka    + tr_bbarkb

      !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------------!
      cpatch%today_nppleaf   (ico) = max(tr_bleaf    * cpatch%nplant(ico), 0.0)
      cpatch%today_nppfroot  (ico) = max(tr_broot    * cpatch%nplant(ico), 0.0)
      cpatch%today_nppsapwood(ico) = max(tr_bsapwood * cpatch%nplant(ico), 0.0)
      cpatch%today_nppbark   (ico) = max(tr_bbark    * cpatch%nplant(ico), 0.0)
      
      cpatch%today_nppdaily  (ico) = npp_actual      * cpatch%nplant(ico)
      !------------------------------------------------------------------------------------!


   end subroutine update_today_npp_vars
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   subroutine update_nitrogen(flushing,ipft,npp_actual,nplant,tr_bleaf,tr_broot            &
                             ,tr_bsapwooda,tr_bsapwoodb,tr_bbarka,tr_bbarkb,tr_bstorage    &
                             ,nitrogen_uptake,fgn_in,fsn_in)
      use consts_coms   , only : tiny_num                 ! ! intent(in)
      use pft_coms      , only : agf_bs                   & ! intent(in)
                               , c2n_storage              & ! intent(in)
                               , c2n_leaf                 & ! intent(in)
                               , c2n_stem                 & ! intent(in)
                               , f_labile_leaf            & ! intent(in)
                               , f_labile_stem            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      logical        , intent(in)    :: flushing
      integer        , intent(in)    :: ipft
      real           , intent(in)    :: npp_actual
      real           , intent(in)    :: nplant
      real           , intent(in)    :: tr_bleaf
      real           , intent(in)    :: tr_broot
      real           , intent(in)    :: tr_bsapwooda
      real           , intent(in)    :: tr_bsapwoodb
      real           , intent(in)    :: tr_bbarka
      real           , intent(in)    :: tr_bbarkb
      real           , intent(in)    :: tr_bstorage
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: fgn_in
      real           , intent(inout) :: fsn_in
      !----- Local variables. -------------------------------------------------------------!
      real                           :: n2c_labile_leaf
      real                           :: n2c_labile_stem
      real                           :: n2c_labile_storage
      real                           :: tr_a_bfast
      real                           :: tr_b_bfast
      real                           :: tr_a_bstruct
      real                           :: tr_b_bstruct
      real                           :: tr_balive
      real                           :: f_tr_bfast
      real                           :: f_tr_bstruct
      real                           :: tr_a_bstorage
      real                           :: tr_b_bstorage
      !------------------------------------------------------------------------------------!


      !------ Partition the components into labile and lignified (and AG and BG). ---------!
      n2c_labile_leaf    = (         f_labile_leaf(ipft)   / c2n_leaf(ipft)                &
                           + ( 1.0 - f_labile_leaf(ipft) ) / c2n_stem(ipft) )
      n2c_labile_stem    = (         f_labile_stem(ipft)   / c2n_leaf(ipft)                &
                           + ( 1.0 - f_labile_stem(ipft) ) / c2n_stem(ipft) )
      n2c_labile_storage = 1.0 / c2n_storage
      tr_a_bfast         = f_labile_leaf(ipft) * tr_bleaf                                  &
                         + f_labile_stem(ipft) * (tr_bsapwooda + tr_bbarka)
      tr_b_bfast         = f_labile_leaf(ipft) * tr_broot                                  &
                         + f_labile_stem(ipft) * (tr_bsapwoodb + tr_bbarkb)
      tr_a_bstruct       = (1.0 - f_labile_leaf(ipft)) * tr_bleaf                          &
                         + (1.0 - f_labile_stem(ipft)) * (tr_bsapwooda + tr_bbarka)
      tr_b_bstruct       = (1.0 - f_labile_leaf(ipft)) * tr_broot                          &
                         + (1.0 - f_labile_stem(ipft)) * (tr_bsapwoodb + tr_bbarkb)
      tr_balive          = tr_a_bfast + tr_b_bfast + tr_a_bstruct + tr_b_bstruct
      tr_a_bstorage      =       agf_bs(ipft)  * tr_bstorage
      tr_b_bstorage      = (1. - agf_bs(ipft)) * tr_bstorage
      !------------------------------------------------------------------------------------!

      if (abs(tr_balive) < tiny_num) then
         f_tr_bfast   = 0.5
         f_tr_bstruct = 0.5
      else
         f_tr_bfast   = max(0.,min(1.,( tr_a_bfast + tr_b_bfast)  / tr_balive))
         f_tr_bstruct = 1.0 - f_tr_bfast
      end if

      if (flushing) then
         !---------------------------------------------------------------------------------!
         !     Check whether there we're adding to or taking from storage.                 !
         !---------------------------------------------------------------------------------!
         if (tr_bstorage <= 0.0)  then
            !------------------------------------------------------------------------------!
            ! We are using all of daily C gain and some of bstorage.                       !
            ! Calculate N demand from using daily C gain.                                  !
            !------------------------------------------------------------------------------!
            if (npp_actual < 0.0) then
               nitrogen_uptake = nitrogen_uptake + npp_actual / c2n_storage
               nitrogen_uptake = nitrogen_uptake                                           &
                               + ( npp_actual      - tr_bstorage                           &
                                 - tr_a_bstruct    - tr_b_bstruct )                        &
                               * ( n2c_labile_leaf -  n2c_labile_storage )                 &
                               + ( npp_actual      - tr_bstorage                           &
                                 - tr_a_bfast      - tr_b_bfast  )                         &
                               * ( n2c_labile_stem -  n2c_labile_storage )
            else
               nitrogen_uptake = nitrogen_uptake                                           &
                               + ( npp_actual     - tr_a_bstruct - tr_b_bstruct )          &
                               * n2c_labile_leaf                                           &
                               + ( npp_actual     - tr_a_bfast   - tr_b_bfast   )          &
                               * n2c_labile_stem

               !---------------------------------------------------------------------------!
               ! Calculate additional N uptake from transfer of C from storage to balive.  !
               !---------------------------------------------------------------------------!
               nitrogen_uptake  = nitrogen_uptake                                          &
                                +  ( - f_tr_bfast   * tr_bstorage )                        &
                                * ( n2c_labile_leaf -  n2c_labile_storage )                &
                                +  ( - f_tr_bstruct * tr_bstorage )                        &
                                * ( n2c_labile_stem -  n2c_labile_storage )
            end if

         else
            !------------------------------------------------------------------------------!
            !     N uptake for fraction of daily C gain going to balive.                   !
            !------------------------------------------------------------------------------!
            nitrogen_uptake = nitrogen_uptake                                              &
                            + (npp_actual     - tr_a_bstruct - tr_b_bstruct - tr_bstorage) &
                            * n2c_labile_leaf                                              &
                            + (npp_actual     - tr_a_bfast   - tr_b_bfast   - tr_bstorage) &
                            * n2c_labile_stem
            !----------- N uptake for fraction of daily C gain going to bstorage. ---------!
            nitrogen_uptake = nitrogen_uptake + tr_bstorage * n2c_labile_storage
         end if
         !---------------------------------------------------------------------------------!
      else
         if (tr_bstorage > 0.0) then
            nitrogen_uptake = nitrogen_uptake + tr_bstorage * n2c_labile_storage
         else
            fgn_in = fgn_in - nplant  * ( tr_a_bstorage * n2c_labile_storage               &
                                        + tr_a_bfast    * n2c_labile_leaf                  &
                                        + tr_a_bstruct  * n2c_labile_stem )
            fsn_in = fsn_in - nplant  * ( tr_b_bstorage * n2c_labile_storage               &
                                        + tr_b_bfast    * n2c_labile_leaf                  &
                                        + tr_b_bstruct  * n2c_labile_stem )
         end if
      end if


   end subroutine update_nitrogen
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   subroutine potential_N_uptake(cpatch,ico,npp_pot,N_uptake_pot,green_leaf_factor)
      use ed_state_vars , only : patchtype     ! ! structure
      use pft_coms      , only : c2n_storage   & ! intent(in)
                               , c2n_leaf      & ! intent(in)
                               , c2n_stem      & ! intent(in)
                               , f_labile_leaf ! ! intent(in)
      use allometry     , only : size2bl       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: npp_pot
      real           , intent(inout) :: N_uptake_pot
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      !------------------------------------------------------------------------------------!

      ipft = cpatch%pft(ico)

      if ( cpatch%phenology_status(ico) == 0 .and. npp_pot > 0.0 ) then

         !----- Positive carbon balance with plants fully flushed. ------------------------!
         N_uptake_pot = N_uptake_pot + npp_pot / c2n_storage
         !---------------------------------------------------------------------------------!

      elseif (cpatch%phenology_status(ico) == 1) then
         !---------------------------------------------------------------------------------!
         !   This calculation of bl_max is wrong for grass, but they should not have       !
         ! phenology_status=1 yet.                                                         !
         ! MLO - I don't see problems as long as phenology(grass) is evergreen.            !
         !---------------------------------------------------------------------------------!
         bl_max = size2bl(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico),ipft)           &
                * green_leaf_factor * cpatch%elongf(ico)
         bl_pot = cpatch%bleaf(ico) + npp_pot

         if (bl_pot > bl_max) then
            !------------------------------------------------------------------------------!
            !     This increment would take us over the limit, so we assign all that can   !
            ! go for leaves to them, and put the remainder in storage.                     !
            !------------------------------------------------------------------------------!
            increment    = npp_pot - (bl_max-cpatch%bleaf(ico))
            N_uptake_pot = N_uptake_pot + increment / c2n_storage
            increment    = bl_max-cpatch%bleaf(ico)
            N_uptake_pot = N_uptake_pot + increment                                        &
                         * (        f_labile_leaf(ipft)  / c2n_leaf(ipft)                  &
                           + (1.0 - f_labile_leaf(ipft)) / c2n_stem(ipft) )
            !------------------------------------------------------------------------------!
         elseif (npp_pot > 0.0) then

            !------------------------------------------------------------------------------!
            !      This increment did not exceed the limit, put everything in leaves.  We  !
            ! don't compute the uptake if carbon balance is negative, just because there   !
            ! will be no uptake...                                                         !
            !------------------------------------------------------------------------------!
            N_uptake_pot = N_uptake_pot + npp_pot                                          &
                         * (        f_labile_leaf(ipft)  / c2n_leaf(ipft)                  &
                           + (1.0 - f_labile_leaf(ipft)) / c2n_stem(ipft) )
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine potential_N_uptake
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine updates litter carbon and nitrogen inputs.                       !
   !---------------------------------------------------------------------------------------!
   subroutine update_litter_inputs(csite,ipa)

      use ed_state_vars, only : patchtype      & ! structure
                              , sitetype       ! ! structure
      use pft_coms     , only : c2n_leaf       & ! intent(in)
                              , c2n_stem       & ! intent(in)
                              , l2n_stem       & ! intent(in)
                              , c2n_storage    & ! intent(in)
                              , f_labile_leaf  & ! intent(in)
                              , f_labile_stem  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: fg_litter
      real                         :: fs_litter
      real                         :: sg_litter
      real                         :: ss_litter
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !      Add fine root, leaf, and bark turnover to the litter.                         !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !----- Split litter into fast (labile) and structural. ---------------------------!
         fg_litter = cpatch%nplant(ico)                                                    &
                  * ( f_labile_leaf(ipft) * cpatch%leaf_maintenance(ico)                   &
                    + f_labile_stem(ipft) * cpatch%barka_maintenance(ico) )
         fs_litter = cpatch%nplant(ico)                                                    &
                  * ( f_labile_leaf(ipft) * cpatch%root_maintenance(ico)                   &
                    + f_labile_stem(ipft) * cpatch%barkb_maintenance(ico) )
         sg_litter = cpatch%nplant(ico)                                                    &
                  * ( (1.0 - f_labile_leaf(ipft)) * cpatch%leaf_maintenance(ico)           &
                    + (1.0 - f_labile_stem(ipft)) * cpatch%barka_maintenance(ico) )
         ss_litter = cpatch%nplant(ico)                                                    &
                  * ( (1.0 - f_labile_leaf(ipft)) * cpatch%root_maintenance(ico)           &
                    + (1.0 - f_labile_stem(ipft)) * cpatch%barkb_maintenance(ico) )
         !---------------------------------------------------------------------------------!



         !------ Update soil carbon and soil nitrogen. ------------------------------------!
         csite%fgc_in (ipa) = csite%fgc_in (ipa) + fg_litter
         csite%fsc_in (ipa) = csite%fsc_in (ipa) + fs_litter
         csite%fgn_in (ipa) = csite%fgn_in (ipa) + fg_litter / c2n_leaf(ipft)
         csite%fsn_in (ipa) = csite%fsn_in (ipa) + fs_litter / c2n_leaf(ipft)
         csite%stgc_in(ipa) = csite%stgc_in(ipa) + sg_litter
         csite%stsc_in(ipa) = csite%stsc_in(ipa) + ss_litter
         csite%stgl_in(ipa) = csite%stgl_in(ipa) + sg_litter * l2n_stem / c2n_stem(ipft)
         csite%stsl_in(ipa) = csite%stsl_in(ipa) + ss_litter * l2n_stem / c2n_stem(ipft)
         csite%stgn_in(ipa) = csite%stgn_in(ipa) + sg_litter / c2n_stem(ipft)
         csite%stsn_in(ipa) = csite%stsn_in(ipa) + ss_litter / c2n_stem(ipft)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     When storage carbon is lost, allow the associated nitrogen to go to litter  !
         ! in order to maintain prescribed C2N ratio.  Carbon does not go to litter        !
         ! because it becomes CO2.                                                         !
         !---------------------------------------------------------------------------------!
         csite%fsn_in(ipa) = csite%fsn_in(ipa)                                             &
                           + csite%commit_storage_resp(ipa) / c2n_storage
         !---------------------------------------------------------------------------------!
      end do
      return
   end subroutine update_litter_inputs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the cohort level.  There are  !
   ! a few, extreme instances in which carbon conservation is not attainable (e.g.         !
   ! severely negative daily NPP when plant has almost no carbon left or minor truncation  !
   ! errors), and we account for them.  Otherwise, we stop any cohort that is attempting   !
   ! to smuggle or to evade carbon.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine check_balive_cohort(csite,ipa,ico,bleaf_in,broot_in,bsapwooda_in             &
                                 ,bsapwoodb_in,bbarka_in,bbarkb_in,balive_in,bstorage_in   &
                                 ,bdeada_in,bdeadb_in,phenstatus_in,metnpp_actual          &
                                 ,npp_actual,growresp_actual,tissue_maintenance            &
                                 ,storage_maintenance,carbon_miss,xfer_case)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use allometry    , only : size2bl            ! ! function
      use budget_utils , only : tol_carbon_budget  ! ! intent(in)
      use pft_coms     , only : q                  & ! intent(in)
                              , qsw                & ! intent(in)
                              , qbark              & ! intent(in)
                              , agf_bs             & ! intent(in)
                              , min_dbh            & ! intent(in)
                              , hgt_min            ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target        :: csite
      integer         , intent(in)    :: ipa
      integer         , intent(in)    :: ico
      integer         , intent(in)    :: phenstatus_in
      integer         , intent(in)    :: xfer_case
      real            , intent(in)    :: bleaf_in
      real            , intent(in)    :: broot_in
      real            , intent(in)    :: bsapwooda_in
      real            , intent(in)    :: bsapwoodb_in
      real            , intent(in)    :: bbarka_in
      real            , intent(in)    :: bbarkb_in
      real            , intent(in)    :: bdeada_in
      real            , intent(in)    :: bdeadb_in
      real            , intent(in)    :: balive_in
      real            , intent(in)    :: bstorage_in
      real            , intent(in)    :: metnpp_actual
      real            , intent(in)    :: npp_actual
      real            , intent(in)    :: growresp_actual
      real            , intent(in)    :: tissue_maintenance
      real            , intent(in)    :: storage_maintenance
      real            , intent(inout) :: carbon_miss
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ipft
      real                         :: bleaf_ok_min
      real                         :: broot_ok_min
      real                         :: bsapwooda_ok_min
      real                         :: bsapwoodb_ok_min
      real                         :: bbarka_ok_min
      real                         :: bbarkb_ok_min
      real                         :: bstorage_ok_min
      real                         :: btotal_in
      real                         :: btotal_fn
      real                         :: delta_btotal
      real                         :: resid_btotal
      real                         :: bstorage_ref
      logical                      :: neg_biomass
      logical                      :: btotal_violation
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=09), parameter :: fmtl='(a,1x,l1)'
      character(len=13), parameter :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Handy aliases. ---------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!


      !----- First, find the minimum possible scale for each pool. ------------------------!
      bleaf_ok_min     = size2bl(min_dbh(ipft),hgt_min(ipft),cpatch%sla(ico),ipft)
      broot_ok_min     = q(ipft) * bleaf_ok_min
      bsapwooda_ok_min =     agf_bs(ipft)  * qsw  (ipft) * cpatch%hite(ico) * bleaf_ok_min
      bsapwoodb_ok_min = (1.-agf_bs(ipft)) * qsw  (ipft) * cpatch%hite(ico) * bleaf_ok_min
      bbarka_ok_min    =     agf_bs(ipft)  * qbark(ipft) * cpatch%hite(ico) * bleaf_ok_min
      bbarkb_ok_min    = (1.-agf_bs(ipft)) * qbark(ipft) * cpatch%hite(ico) * bleaf_ok_min
      bstorage_ok_min  = bleaf_ok_min
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Because the changes may be more closely related with thedaily NPP than with    !
      ! the initial storage pool when storage is almost zero but NPP is not, we combine    !
      ! storage and absolute NPP as the reference scale for storage, to avoid false alarms !
      ! of unacceptable negative biomass when storage is nearly zero.                      !
      !------------------------------------------------------------------------------------!
      bstorage_ref = bstorage_in + abs(npp_actual)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Then, if possible, set the minimum acceptable deviation based on the input, to !
      ! avoid false alarms due to truncation errors when the pool is much greater than the !
      ! minimum size but loses all its stocks.   Storage does not follow a stable          !
      ! allometric relationship.  Because the changes may be more closely related with the !
      ! daily NPP than with the initial storage pool when storage is almost zero but NPP   !
      ! is not, we make sure that to scale the tolerance in such way that does not create  !
      ! false alarms.                                                                      !
      !------------------------------------------------------------------------------------!
      bleaf_ok_min     = - tol_carbon_budget * max(bleaf_in    , bleaf_ok_min    )
      broot_ok_min     = - tol_carbon_budget * max(broot_in    , broot_ok_min    )
      bsapwooda_ok_min = - tol_carbon_budget * max(bsapwooda_in, bsapwooda_ok_min)
      bsapwoodb_ok_min = - tol_carbon_budget * max(bsapwoodb_in, bsapwoodb_ok_min)
      bbarka_ok_min    = - tol_carbon_budget * max(bbarka_in   , bbarka_ok_min   )
      bbarkb_ok_min    = - tol_carbon_budget * max(bbarkb_in   , bbarkb_ok_min   )
      bstorage_ok_min  = - tol_carbon_budget * max(bstorage_ref, bstorage_ok_min )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check every tissue and storage, to make sure they have sensible numbers.  Tiny  !
      ! negative stocks will be tolerated, but don't fix anything in case any of the pools !
      ! are too negative.                                                                  !
      !------------------------------------------------------------------------------------!
      neg_biomass    = cpatch%bleaf    (ico) < bleaf_ok_min     .or.                       &
                       cpatch%broot    (ico) < broot_ok_min     .or.                       &
                       cpatch%bsapwooda(ico) < bsapwooda_ok_min .or.                       &
                       cpatch%bsapwoodb(ico) < bsapwoodb_ok_min .or.                       &
                       cpatch%bbarka   (ico) < bbarka_ok_min    .or.                       &
                       cpatch%bbarkb   (ico) < bbarkb_ok_min    .or.                       &
                       cpatch%bstorage (ico) < bstorage_ok_min
      if (.not. neg_biomass) then
         !----- Account for any potential violation of carbon stocks. ---------------------!
         carbon_miss = carbon_miss                                                         &
                     - min(cpatch%bleaf    (ico),0.0) - min(cpatch%broot    (ico),0.0)     &
                     - min(cpatch%bsapwooda(ico),0.0) - min(cpatch%bsapwoodb(ico),0.0)     &
                     - min(cpatch%bbarka   (ico),0.0) - min(cpatch%bbarkb   (ico),0.0)     &
                     - min(cpatch%bstorage (ico),0.0)
         !---------------------------------------------------------------------------------!


         !----- Make sure that all pools are non-negative. --------------------------------!
         cpatch%bleaf    (ico) = max(cpatch%bleaf    (ico),0.0)
         cpatch%broot    (ico) = max(cpatch%broot    (ico),0.0)
         cpatch%bsapwooda(ico) = max(cpatch%bsapwooda(ico),0.0)
         cpatch%bsapwoodb(ico) = max(cpatch%bsapwoodb(ico),0.0)
         cpatch%bbarka   (ico) = max(cpatch%bbarka   (ico),0.0)
         cpatch%bbarkb   (ico) = max(cpatch%bbarkb   (ico),0.0)
         cpatch%bstorage (ico) = max(cpatch%bstorage (ico),0.0)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after active tissue growth.                     !
      !------------------------------------------------------------------------------------!
      btotal_in         = balive_in + bdeada_in + bdeadb_in + bstorage_in
      btotal_fn         = cpatch%balive(ico) + cpatch%bdeada  (ico)                        &
                        + cpatch%bdeadb(ico) + cpatch%bstorage(ico)
      delta_btotal      = npp_actual      - tissue_maintenance - storage_maintenance
      resid_btotal      = btotal_fn - btotal_in - delta_btotal - carbon_miss
      btotal_violation  = abs(resid_btotal) > (tol_carbon_budget * btotal_in)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( neg_biomass .or. btotal_violation ) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|       !!!   Cohort Balive budget failed   !!!      |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
         write(unit=*,fmt=fmti )  ' COHORT              : ',ico
         write(unit=*,fmt=fmti )  ' IPFT                : ',ipft
         write(unit=*,fmt=fmtf )  ' DBH                 : ',cpatch%dbh(ico)
         write(unit=*,fmt=fmtf )  ' HITE                : ',cpatch%hite(ico)
         write(unit=*,fmt=fmtl )  ' NEG_BIOMASS         : ',neg_biomass
         write(unit=*,fmt=fmtl )  ' BTOTAL_VIOLATION    : ',btotal_violation
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmti )  ' PHENSTATUS_IN       : ',phenstatus_in
         write(unit=*,fmt=fmtf )  ' BLEAF_IN            : ',bleaf_in
         write(unit=*,fmt=fmtf )  ' BROOT_IN            : ',broot_in
         write(unit=*,fmt=fmtf )  ' BSAPWOODA_IN        : ',bsapwooda_in
         write(unit=*,fmt=fmtf )  ' BSAPWOODB_IN        : ',bsapwoodb_in
         write(unit=*,fmt=fmtf )  ' BBARKA_IN           : ',bbarka_in
         write(unit=*,fmt=fmtf )  ' BBARKB_IN           : ',bbarkb_in
         write(unit=*,fmt=fmtf )  ' BDEADA_IN           : ',bdeada_in
         write(unit=*,fmt=fmtf )  ' BDEADB_IN           : ',bdeadb_in
         write(unit=*,fmt=fmtf )  ' BALIVE_IN           : ',balive_in
         write(unit=*,fmt=fmtf )  ' BSTORAGE_IN         : ',bstorage_in
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmti )  ' PHENSTATUS_FN       : ',cpatch%phenology_status(ico)
         write(unit=*,fmt=fmtf )  ' BLEAF_FN            : ',cpatch%bleaf           (ico)
         write(unit=*,fmt=fmtf )  ' BROOT_FN            : ',cpatch%broot           (ico)
         write(unit=*,fmt=fmtf )  ' BSAPWOODA_FN        : ',cpatch%bsapwooda       (ico)
         write(unit=*,fmt=fmtf )  ' BSAPWOODB_FN        : ',cpatch%bsapwoodb       (ico)
         write(unit=*,fmt=fmtf )  ' BBARKA_FN           : ',cpatch%bbarka          (ico)
         write(unit=*,fmt=fmtf )  ' BBARKB_FN           : ',cpatch%bbarkb          (ico)
         write(unit=*,fmt=fmtf )  ' BDEADA_FN           : ',cpatch%bdeada          (ico)
         write(unit=*,fmt=fmtf )  ' BDEADB_FN           : ',cpatch%bdeadb          (ico)
         write(unit=*,fmt=fmtf )  ' BSTORAGE_FN         : ',cpatch%bstorage        (ico)
         write(unit=*,fmt=fmtf )  ' BALIVE_FN           : ',cpatch%balive          (ico)
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' TOL_CARBON_BUDGET   : ',tol_carbon_budget
         write(unit=*,fmt=fmtf )  ' BLEAF_OK_MIN        : ',bleaf_ok_min
         write(unit=*,fmt=fmtf )  ' BROOT_OK_MIN        : ',broot_ok_min
         write(unit=*,fmt=fmtf )  ' BSAPWOODA_OK_MIN    : ',bsapwooda_ok_min
         write(unit=*,fmt=fmtf )  ' BSAPWOODB_OK_MIN    : ',bsapwoodb_ok_min
         write(unit=*,fmt=fmtf )  ' BBARKA_OK_MIN       : ',bbarka_ok_min
         write(unit=*,fmt=fmtf )  ' BBARKB_OK_MIN       : ',bbarkb_ok_min
         write(unit=*,fmt=fmtf )  ' BSTORAGE_OK_MIN     : ',bstorage_ok_min
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmti )  ' XFER_CASE           : ',xfer_case
         write(unit=*,fmt=fmtf )  ' METABOLIC_NPP       : ',metnpp_actual
         write(unit=*,fmt=fmtf )  ' NPP                 : ',npp_actual
         write(unit=*,fmt=fmtf )  ' GROWTH_RESPIRATION  : ',growresp_actual
         write(unit=*,fmt=fmtf )  ' LEAF_MAINTENANCE    : ',cpatch%leaf_maintenance(ico)
         write(unit=*,fmt=fmtf )  ' ROOT_MAINTENANCE    : ',cpatch%root_maintenance(ico)
         write(unit=*,fmt=fmtf )  ' BALIVE_MAINTENANCE  : ',tissue_maintenance
         write(unit=*,fmt=fmtf )  ' STORAGE_MAINTENANCE : ',storage_maintenance
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BTOTAL_IN           : ',btotal_in
         write(unit=*,fmt=fmtf )  ' BTOTAL_FN           : ',btotal_fn
         write(unit=*,fmt=fmtf )  ' DELTA_BTOTAL        : ',delta_btotal
         write(unit=*,fmt=fmtf )  ' CARBON_MISS         : ',carbon_miss
         write(unit=*,fmt=fmtf )  ' RESIDUAL_BTOTAL     : ',resid_btotal
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  ' '


         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_balive_cohort','growth_balive.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_balive_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the patch level.  There are   !
   ! a few, extreme instances in which carbon conservation is not attainable (e.g.         !
   ! severely negative daily NPP when plant has almost no carbon left or minor truncation  !
   ! errors), and we account for them.  Otherwise, we stop any patch that is attempting    !
   ! to smuggle or to evade carbon.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine check_balive_patch(csite,ipa,fgc_in_in,fsc_in_in,stgc_in_in,stsc_in_in       &
                                ,pat_balive_in,pat_bdead_in,pat_bstorage_in                &
                                ,pat_metnpp_actual,pat_npp_actual                          &
                                ,pat_tissue_maintenance,pat_carbon_miss)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use budget_utils , only : tol_carbon_budget  ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target        :: csite
      integer         , intent(in)    :: ipa
      real            , intent(in)    :: fgc_in_in
      real            , intent(in)    :: fsc_in_in
      real            , intent(in)    :: stgc_in_in
      real            , intent(in)    :: stsc_in_in
      real            , intent(in)    :: pat_balive_in
      real            , intent(in)    :: pat_bdead_in
      real            , intent(in)    :: pat_bstorage_in
      real            , intent(in)    :: pat_metnpp_actual
      real            , intent(in)    :: pat_npp_actual
      real            , intent(in)    :: pat_tissue_maintenance
      real            , intent(in)    :: pat_carbon_miss
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer       :: cpatch
      integer                         :: ico
      real                            :: fgc_in_fn
      real                            :: fsc_in_fn
      real                            :: stgc_in_fn
      real                            :: stsc_in_fn
      real                            :: pat_balive_fn
      real                            :: pat_bdead_fn
      real                            :: pat_bstorage_fn
      real                            :: pat_btotal_in
      real                            :: pat_btotal_fn
      real                            :: pat_biomass_in
      real                            :: pat_biomass_fn
      real                            :: pat_storage_maintenance
      real                            :: soilc_in_in
      real                            :: soilc_in_fn
      real                            :: resid_pat_btotal
      logical                         :: pat_btotal_violation
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=13), parameter :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Handy aliases. ---------------------------------------------------------------!
      cpatch                  => csite%patch  (ipa)
      fgc_in_fn               =  csite%fgc_in (ipa)
      fsc_in_fn               =  csite%fsc_in (ipa)
      stgc_in_fn              =  csite%stgc_in(ipa)
      stsc_in_fn              =  csite%stsc_in(ipa)
      pat_storage_maintenance = csite%commit_storage_resp(ipa)
      !------------------------------------------------------------------------------------!



      !----- Count current stocks. --------------------------------------------------------!
      pat_balive_fn = 0.0
      pat_bdead_fn  = 0.0
      pat_bstorage_fn = 0.0
      do ico=1,cpatch%ncohorts
         pat_balive_fn   = pat_balive_fn   + cpatch%nplant(ico) * cpatch%balive  (ico)
         pat_bdead_fn    = pat_bdead_fn                                                    &
                         + cpatch%nplant(ico) * ( cpatch%bdeada(ico) + cpatch%bdeadb(ico) )
         pat_bstorage_fn = pat_bstorage_fn + cpatch%nplant(ico) * cpatch%bstorage(ico)
      end do
      !------------------------------------------------------------------------------------!


      !------ Summary of the carbon stocks. -----------------------------------------------!
      pat_biomass_in = pat_balive_in + pat_bdead_in + pat_bstorage_in
      pat_biomass_fn = pat_balive_fn + pat_bdead_fn + pat_bstorage_fn
      soilc_in_in    = fgc_in_in + fsc_in_in + stgc_in_in + stsc_in_in
      soilc_in_fn    = fgc_in_fn + fsc_in_fn + stgc_in_fn + stsc_in_fn
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after active tissue growth.  Storage            !
      ! maintenance is positive in the residual calculation because it is a committed loss !
      ! that doesn't go to the soil carbon, but to the canopy air space.                   !
      !------------------------------------------------------------------------------------!
      pat_btotal_in        = pat_biomass_in + soilc_in_in
      pat_btotal_fn        = pat_biomass_fn + soilc_in_fn
      resid_pat_btotal     = pat_btotal_fn - pat_btotal_in                                 &
                           + pat_storage_maintenance - pat_npp_actual - pat_carbon_miss
      pat_btotal_violation = abs(resid_pat_btotal) > (tol_carbon_budget * pat_btotal_in)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( pat_btotal_violation ) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|       !!!   Patch Balive budget failed   !!!       |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
         write(unit=*,fmt=fmti )  ' DIST_TYPE           : ',csite%dist_type(ipa)
         write(unit=*,fmt=fmti )  ' NCOHORTS            : ',cpatch%ncohorts
         write(unit=*,fmt=fmtf )  ' AGE                 : ',csite%age(ipa)
         write(unit=*,fmt=fmtf )  ' VEG_HEIGHT          : ',csite%veg_height(ipa)
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BALIVE_IN           : ',pat_balive_in
         write(unit=*,fmt=fmtf )  ' BDEAD_IN            : ',pat_bdead_in
         write(unit=*,fmt=fmtf )  ' BSTORAGE_IN         : ',pat_bstorage_in
         write(unit=*,fmt=fmtf )  ' FGC_IN_IN           : ',fgc_in_in
         write(unit=*,fmt=fmtf )  ' FSC_IN_IN           : ',fsc_in_in
         write(unit=*,fmt=fmtf )  ' STGC_IN_IN          : ',stgc_in_in
         write(unit=*,fmt=fmtf )  ' STSC_IN_IN          : ',stsc_in_in
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BALIVE_FN           : ',pat_balive_fn
         write(unit=*,fmt=fmtf )  ' BDEAD_FN            : ',pat_bdead_fn
         write(unit=*,fmt=fmtf )  ' BSTORAGE_FN         : ',pat_bstorage_fn
         write(unit=*,fmt=fmtf )  ' FGC_IN_FN           : ',fgc_in_fn
         write(unit=*,fmt=fmtf )  ' FSC_IN_FN           : ',fsc_in_fn
         write(unit=*,fmt=fmtf )  ' STGC_IN_FN          : ',stgc_in_fn
         write(unit=*,fmt=fmtf )  ' STSC_IN_FN          : ',stsc_in_fn
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BIOMASS_IN          : ',pat_biomass_in
         write(unit=*,fmt=fmtf )  ' SOILC_IN_IN         : ',soilc_in_in
         write(unit=*,fmt=fmtf )  ' BIOMASS_FN          : ',pat_biomass_fn
         write(unit=*,fmt=fmtf )  ' SOILC_IN_FN         : ',soilc_in_fn
         write(unit=*,fmt=fmtf )  ' METABOLIC_NPP       : ',pat_metnpp_actual
         write(unit=*,fmt=fmtf )  ' NPP                 : ',pat_npp_actual
         write(unit=*,fmt=fmtf )  ' TISSUE_MAINTENANCE  : ',pat_tissue_maintenance
         write(unit=*,fmt=fmtf )  ' STORAGE_MAINTENANCE : ',pat_storage_maintenance
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BTOTAL_IN           : ',pat_btotal_in
         write(unit=*,fmt=fmtf )  ' BTOTAL_FN           : ',pat_btotal_fn
         write(unit=*,fmt=fmtf )  ' CARBON_MISS         : ',pat_carbon_miss
         write(unit=*,fmt=fmtf )  ' RESIDUAL_BTOTAL     : ',resid_pat_btotal
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  ' '


         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_balive_patch','growth_balive.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_balive_patch
   !=======================================================================================!
   !=======================================================================================!
end module growth_balive
!==========================================================================================!
!==========================================================================================!

