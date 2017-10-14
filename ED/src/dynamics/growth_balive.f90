!==========================================================================================!
!==========================================================================================!
! MODULE: GROWTH_BALIVE
!
!> \brief   Various routines handling plant C and N use given size and allometry.
!> \details Essentially, this file contains dbalive_dt[_eq0] and their libraries of fns.
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
   subroutine dbalive_dt(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , qbark                  & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , is_grass               & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , size2bt                & ! subroutine
                                 , size2xb                & ! subroutine
                                 , ed_biomass             ! ! function
      use mortality       , only : mortality_rates        ! ! subroutine
      use fuse_fiss_utils , only : sort_cohorts           ! ! subroutine
      use ed_misc_coms    , only : igrass                 & ! intent(in)
                                 , growth_resp_scheme     & ! intent(in)
                                 , storage_resp_scheme    & ! intent(in)
                                 , ibigleaf               ! ! intent(in)
      use budget_utils    , only : update_budget          ! ! sub-routine
      use consts_coms   , only : tiny_num     ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid !< the ed grid
      real             , intent(in) :: tfact !< "time factor" i.e. call frequency
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_lightmax
      real                          :: carbon_balance_moistmax
      real                          :: carbon_balance_mlmax
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: dlnndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      real                          :: temp_dep
      real                          :: growth_resp_int          ! Growth resp / balive
      real                          :: storage_resp_int         ! Growth resp / balive
      real                          :: tr_bleaf
      real                          :: tr_broot
      real                          :: tr_bbark
      real                          :: tr_bsapwooda
      real                          :: tr_bsapwoodb
      real                          :: tr_bstorage
      real                          :: cb_decrement
      real                          :: carbon_debt
      real                          :: balive_aim
      logical                       :: flushing
      logical                       :: on_allometry
      !------------------------------------------------------------------------------------!


      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0

               !----- Loop over cohorts. --------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0
                  
                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + q(ipft) + (qsw(ipft) + qbark(ipft)) * cpatch%hite(ico)
                  salloci = 1.0 / salloc

                  !------------------------------------------------------------------------!
                  !     Compute and apply maintenance costs and get daily C gain.          !
                  !------------------------------------------------------------------------!
                  call get_maintenance(cpatch,ico,tfact,csite%avg_daily_temp(ipa))
                  call get_daily_C_gain(cpatch,ico,daily_C_gain)
                  call apply_maintenance(cpatch,ico,tfact,cb_decrement)

                  call update_cb(cpatch,ico,cb_decrement)

                  !------------------------------------------------------------------------!
                  !    Storage respiration/turnover_rate.                                  !
                  !    Calculate in same way as leaf and root turnover in kgC/plant/year.  !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The commented line is an experimental and arbitrary test, borrowed !
                  ! from maintainence temperature dependency. [[MCD]]                      !
                  !------------------------------------------------------------------------!
                  ! temp_dep = 1.0                                                         &
                  !          / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
                  temp_dep = 1.0
                  !------------------------------------------------------------------------!
                  select case(storage_resp_scheme)
                  case(0)
                     cpatch%leaf_storage_resp(ico) = 0.0
                     cpatch%root_storage_resp(ico) = 0.0
                     cpatch%sapa_storage_resp(ico) = cpatch%bstorage(ico)                  &
                                                   * storage_turnover_rate(ipft)           &
                                                   * tfact * temp_dep
                     cpatch%sapb_storage_resp(ico) = 0.0
                     cpatch%bark_storage_resp(ico) = 0.0

                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                            - cpatch%sapa_storage_resp(ico)
                  case(1)
                     storage_resp_int = cpatch%bstorage(ico) / cpatch%balive(ico)          &
                                      * storage_turnover_rate(ipft) * tfact * temp_dep

                     cpatch%leaf_storage_resp(ico) = storage_resp_int *cpatch%bleaf(ico)
                     cpatch%root_storage_resp(ico) = storage_resp_int *cpatch%broot(ico)
                     cpatch%sapa_storage_resp(ico) = storage_resp_int *cpatch%bsapwooda(ico)
                     cpatch%sapb_storage_resp(ico) = storage_resp_int *cpatch%bsapwoodb(ico)
                     cpatch%bark_storage_resp(ico) = storage_resp_int *cpatch%bbark(ico)

                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                          - cpatch%leaf_storage_resp(ico)                  &
                                          - cpatch%root_storage_resp(ico)                  &
                                          - cpatch%sapa_storage_resp(ico)                  &
                                          - cpatch%sapb_storage_resp(ico)                  &
                                          - cpatch%bark_storage_resp(ico)
                  end select

                  !------------------------------------------------------------------------!
                  !     When storage carbon is lost, allow the associated nitrogen to go   !
                  ! to litter in order to maintain prescribed C2N ratio.                   !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa) = csite%fsn_in(ipa) + (cpatch%leaf_storage_resp(ico)   &
                                                        +  cpatch%root_storage_resp(ico)   &
                                                        +  cpatch%sapa_storage_resp(ico)   &
                                                        +  cpatch%sapb_storage_resp(ico)   &
                                                        +  cpatch%bark_storage_resp(ico))  &
                                                        / c2n_storage * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_lightmax    &
                                            ,carbon_balance_moistmax,carbon_balance_mlmax)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  select case(growth_resp_scheme)
                  case(0)
                     cpatch%sapa_growth_resp(ico) = max(0.0, daily_C_gain                  &
                                                             * growth_resp_factor(ipft))
                     cpatch%leaf_growth_resp(ico) = 0.0
                     cpatch%root_growth_resp(ico) = 0.0
                     cpatch%sapb_growth_resp(ico) = 0.0
                     cpatch%bark_growth_resp(ico) = 0.0
                  case(1)
                     growth_resp_int = max(0.0, daily_C_gain * growth_resp_factor(ipft)    &
                                                             / cpatch%balive(ico))
                     
                     cpatch%leaf_growth_resp(ico) = growth_resp_int * cpatch%bleaf    (ico)
                     cpatch%root_growth_resp(ico) = growth_resp_int * cpatch%broot    (ico)
                     cpatch%sapa_growth_resp(ico) = growth_resp_int * cpatch%bsapwooda(ico)
                     cpatch%sapb_growth_resp(ico) = growth_resp_int * cpatch%bsapwoodb(ico)
                     cpatch%bark_growth_resp(ico) = growth_resp_int * cpatch%bbark    (ico)
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and bstorage.             !
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)

                  if(is_grass(ipft).and. igrass==1) then
                      call alloc_plant_c_balance_grass(csite,ipa,ico,salloc,salloci        &
                                                ,carbon_balance,nitrogen_uptake            &
                                                ,cpoly%green_leaf_factor(ipft,isi))
                      call sort_cohorts(cpatch)
                  else
                     call get_c_xfers(csite,ipa,ico,carbon_balance                         &
                                     ,cpoly%green_leaf_factor(ipft,isi),tr_bleaf,tr_broot  &
                                     ,tr_bsapwooda,tr_bsapwoodb,tr_bbark,tr_bstorage       &
                                     ,carbon_debt,flushing,balive_aim)

                     call apply_c_xfers(cpatch,ico,carbon_balance,tr_bleaf,tr_broot        &
                                       ,tr_bsapwooda,tr_bsapwoodb,tr_bbark,tr_bstorage)
                     
                     call update_today_npp_vars(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda  &
                                                ,tr_bsapwoodb,tr_bbark,carbon_balance)

                     call update_nitrogen(flushing,ipft,carbon_balance,cpatch%nplant(ico)  &
                                          ,tr_bleaf,tr_broot,tr_bstorage,tr_bbark          &
                                          ,nitrogen_uptake,csite%fsn_in(ipa))

                     !---------------------------------------------------------------------!
                     !     Update the phenology status.                                    !
                     !---------------------------------------------------------------------!
                     on_allometry = balive_aim - cpatch%balive(ico) <= 0.000001*balive_aim
                     if (flushing .and. cpatch%elongf(ico) == 1.0 .and. on_allometry) then
                        cpatch%phenology_status(ico) = 0
                     elseif(cpatch%bleaf(ico) < tiny_num .and.                             &
                            (cpatch%phenology_status(ico) == 0 .or.                        &
                             cpatch%phenology_status(ico) == 1)) then
                        cpatch%phenology_status(ico) = 1
                     end if
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,balive_in,carbon_balance_pot       &
                                            ,N_uptake_pot                                  &
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
                  !      Update mortality rates.  Notice that the only mortality rate that !
                  ! changes daily is the frost mortality, and the disturbance mortality is !
                  ! not updated here (it is updated in the main disturbance procedure).    !
                  !                                                                        !
                  !      We no longer include mortality rates due to disturbance in the    !
                  ! big-leaf simulations, this is now done at disturbance.f90.             !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ico,csite%avg_daily_temp(ipa),csite%age(ipa) &
                                      ,csite%dist_type(ipa))
                  dlnndt   = - sum(cpatch%mort_rate(1:4,ico))
                  dndt     = dlnndt * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!

                  !----- Update monthly mortality rates [plants/m2/month and 1/month]. ----!
                  cpatch%monthly_dndt  (ico) = cpatch%monthly_dndt  (ico) + dndt   * tfact
                  cpatch%monthly_dlnndt(ico) = cpatch%monthly_dlnndt(ico) + dlnndt * tfact
                  !------------------------------------------------------------------------!


                  !----- Updating LAI, WAI, and CAI. --------------------------------------!
                  call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%dbh(ico)   &
                                   ,cpatch%hite(ico) ,cpatch%pft(ico),cpatch%sla(ico)      &
                                   ,cpatch%lai(ico),cpatch%wai(ico),cpatch%crown_area(ico))
                  !------------------------------------------------------------------------!



                  !----- Update above-ground biomass. -------------------------------------!
                  cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)         &
                                              ,cpatch%bsapwooda(ico),cpatch%bbark(ico)     &
                                              ,cpatch%pft(ico))
                  !------------------------------------------------------------------------!



                  !----- Update (commercial) timber biomass. ------------------------------!
                  cpatch%btimber(ico) = size2bt( cpatch%dbh(ico),cpatch%hite(ico)          &
                                               , cpatch%bdead(ico),cpatch%bsapwooda(ico)   &
                                               , cpatch%bbark(ico),cpatch%pft(ico) )
                  !------------------------------------------------------------------------!



                  !----- Update bark thickness. -------------------------------------------!
                  cpatch%thbark(ico)  = size2xb( cpatch%dbh(ico),cpatch%hite(ico)          &
                                               , cpatch%bbark(ico),cpatch%pft(ico) )
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore, update           !
                  ! vegetation energy and heat capacity.                                   !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap         = cpatch%leaf_hcap(ico)
                  old_wood_hcap         = cpatch%wood_hcap(ico)
                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdead(ico)                  &
                                    ,cpatch%bsapwooda(ico),cpatch%bbark(ico)               &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
                  !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico)
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!

               !----- Update litter. ------------------------------------------------------!
               call litter(csite,ipa)
               !---------------------------------------------------------------------------!

               !----- Update patch LAI, WAI, height, roughness... -------------------------!
               call update_patch_derived_props(csite,ipa)
               !---------------------------------------------------------------------------!

               !----- Recalculate storage terms (for budget assessment). ------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)
               !---------------------------------------------------------------------------!

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
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
   !     This subroutine will compute the terms relative to growth of the living tissues,  !
   ! without actually updating it.                                                         !
   ! IMPORTANT: The order of the operations here affect the C/N budgets, so don't change   !
   !            the order of the operations unless you really know what you are doing.     !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt_eq_0(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , qbark                  & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use ed_misc_coms    , only : ibigleaf               & ! intent(in)
                                 , growth_resp_scheme     & ! intent(in)
                                 , storage_resp_scheme    ! ! intent(in)
      use mortality       , only : mortality_rates        ! ! subroutine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_lightmax
      real                          :: carbon_balance_moistmax
      real                          :: carbon_balance_mlmax
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: dlnndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      real                          :: temp_dep
      real                          :: growth_resp_int         ! Growth resp / balive
      real                          :: storage_resp_int        ! Growth resp / balive
      !------------------------------------------------------------------------------------!


      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0

               !----- Loop over cohorts. --------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)


                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0
                  
                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + q(ipft) + (qsw(ipft) + qbark(ipft)) * cpatch%hite(ico)
                  salloci = 1.0 / salloc

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs using actual pools.                      !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,cpatch%broot(ico),cpatch%bleaf(ico)    &
                                        ,cpatch%bbark(ico),tfact,daily_C_gain              &
                                        ,csite%avg_daily_temp(ipa))

                  !------------------------------------------------------------------------!
                  !    For the no vegetation dynamics case, we update the carbon balance   !
                  ! but we do NOT update the living tissues.                               !
                  !------------------------------------------------------------------------!
                  cpatch%cb         (13,ico) = cpatch%cb                  (13,ico)         &
                                             - cpatch%leaf_maintenance       (ico)         &
                                             - cpatch%root_maintenance       (ico)         &
                                             - cpatch%bark_maintenance       (ico)
                  cpatch%cb_lightmax(13,ico) = cpatch%cb_lightmax         (13,ico)         &
                                             - cpatch%leaf_maintenance       (ico)         &
                                             - cpatch%root_maintenance       (ico)         &
                                             - cpatch%bark_maintenance       (ico)
                  cpatch%cb_moistmax(13,ico) = cpatch%cb_moistmax         (13,ico)         &
                                             - cpatch%leaf_maintenance       (ico)         &
                                             - cpatch%root_maintenance       (ico)         &
                                             - cpatch%bark_maintenance       (ico)
                  cpatch%cb_mlmax(13,ico)    = cpatch%cb_mlmax            (13,ico)         &
                                             - cpatch%leaf_maintenance       (ico)         &
                                             - cpatch%root_maintenance       (ico)         &
                                             - cpatch%bark_maintenance       (ico)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Storage respriation/turnover_rate.                                  !
                  !    Calculate in same way as leaf and root turnover in kgC/plant/year.  !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The commented line is an experimental and arbitrary test, borrowed !
                  ! from maintainence temperature dependency. [[MCD]]                      !
                  !------------------------------------------------------------------------!
                  ! temp_dep = 1.0                                                         &
                  !          / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
                  temp_dep = 1.0
                  !------------------------------------------------------------------------!
                  select case(storage_resp_scheme)
                  case(0)
                     cpatch%sapa_storage_resp(ico) = cpatch%bstorage(ico)                  &
                                                   * storage_turnover_rate(ipft)           &
                                                   * tfact * temp_dep
                     cpatch%leaf_storage_resp(ico) = 0.0
                     cpatch%root_storage_resp(ico) = 0.0
                     cpatch%sapb_storage_resp(ico) = 0.0
                     cpatch%bark_storage_resp(ico) = 0.0

                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                          - cpatch%sapa_storage_resp(ico)
                  case(1)
                     storage_resp_int = cpatch%bstorage(ico) * storage_turnover_rate(ipft) &
                                      * tfact * temp_dep

                     cpatch%leaf_storage_resp(ico) = storage_resp_int *cpatch%bleaf(ico)
                     cpatch%root_storage_resp(ico) = storage_resp_int *cpatch%broot(ico)
                     cpatch%sapa_storage_resp(ico) = storage_resp_int *cpatch%bsapwooda(ico)
                     cpatch%sapb_storage_resp(ico) = storage_resp_int *cpatch%bsapwoodb(ico)
                     cpatch%bark_storage_resp(ico) = storage_resp_int *cpatch%bbark(ico)

                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                            - cpatch%leaf_storage_resp(ico)                &
                                            - cpatch%root_storage_resp(ico)                &
                                            - cpatch%sapa_storage_resp(ico)                &
                                            - cpatch%sapb_storage_resp(ico)                &
                                            - cpatch%bark_storage_resp(ico)
                  end select

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_lightmax    &
                                            ,carbon_balance_moistmax,carbon_balance_mlmax)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  select case(growth_resp_scheme)
                  case(0)
                     cpatch%sapa_growth_resp(ico) = max(0.0, daily_C_gain                  &
                                                             * growth_resp_factor(ipft))
                     cpatch%leaf_growth_resp(ico) = 0.0
                     cpatch%root_growth_resp(ico) = 0.0
                     cpatch%sapb_growth_resp(ico) = 0.0
                     cpatch%bark_growth_resp(ico) = 0.0
                  case(1)
                     growth_resp_int = max(0.0, daily_C_gain * growth_resp_factor(ipft)    &
                                                             / cpatch%balive(ico))
                     
                     cpatch%leaf_growth_resp(ico) = growth_resp_int * cpatch%bleaf(ico)
                     cpatch%root_growth_resp(ico) = growth_resp_int * cpatch%broot(ico)
                     cpatch%sapa_growth_resp(ico) = growth_resp_int * cpatch%bsapwooda(ico)
                     cpatch%sapb_growth_resp(ico) = growth_resp_int * cpatch%bsapwoodb(ico)
                     cpatch%bark_growth_resp(ico) = growth_resp_int * cpatch%bbark(ico)
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Find the potential allocation to the living tissues, but don't    !
                  ! update them.                                                           !
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)
                  call alloc_plant_c_balance_eq_0(csite,ipa,ico,carbon_balance             &
                                                 ,nitrogen_uptake                          &
                                                 ,cpoly%green_leaf_factor(ipft,isi))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,balive_in,carbon_balance_pot       &
                                            ,N_uptake_pot                                  &
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
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ico,csite%avg_daily_temp(ipa),csite%age(ipa) &
                                      ,csite%dist_type(ipa))
                  select case (ibigleaf)
                  case (0)
                     dlnndt   = - sum(cpatch%mort_rate(1:4,ico))
                     dndt     = dlnndt * cpatch%nplant(ico)
                  case (1)
                     dlnndt   = - sum(cpatch%mort_rate(1:5,ico))
                     dndt     = dlnndt * cpatch%nplant(ico)
                  end select
                  !------------------------------------------------------------------------!


                  !----- Update monthly mortality rates [plants/m2/month and 1/month]. ----!
                  cpatch%monthly_dndt  (ico) = cpatch%monthly_dndt  (ico) + dndt   * tfact
                  cpatch%monthly_dlnndt(ico) = cpatch%monthly_dlnndt(ico) + dlnndt * tfact
                  !------------------------------------------------------------------------!

               end do cohortloop
               !---------------------------------------------------------------------------!

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine dbalive_dt_eq_0
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   subroutine get_maintenance(cpatch,ico,tfact,tempk)
      use ed_state_vars, only : patchtype             ! ! structure
      use pft_coms     , only : phenology             & ! intent(in)
                              , leaf_turnover_rate    & ! intent(in)
                              , root_turnover_rate    & ! intent(in)
                              , bark_turnover_rate    & ! intent(in)
                              , storage_turnover_rate ! ! intent(in)
      use ed_misc_coms , only : storage_resp_scheme   ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(in)   :: tfact
      real           , intent(in)   :: tempk
      !----- Local variables. -------------------------------------------------------------!
      integer                       :: ipft
      real                          :: maintenance_temp_dep
      real                          :: temp_dep
      real                          :: storage_resp_int
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      !      Find the maintenance costs.  This will depend on the type of phenology that   !
      ! the PFT has.   The tfact term applied converts the maintenance rates to            !
      ! [kgC/plant/day].                                                                   !
      !------------------------------------------------------------------------------------!
      cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) *cpatch%bleaf(ico) *tfact
      cpatch%root_maintenance(ico) = root_turnover_rate(ipft) *cpatch%broot(ico) *tfact
      cpatch%bark_maintenance(ico) = bark_turnover_rate(ipft) *cpatch%bbark(ico) *tfact
      
      select case (phenology(ipft))
      case (0)
         !---------------------------------------------------------------------------------!
         !     Evergreens, like pines.  The turnover rates will be adjusted by a function  !
         ! of temperature, which approaches 0 as the temperature goes down.                !
         !---------------------------------------------------------------------------------!
         !------ Find a temperature dependence adjustment. --------------------------------!
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
         !----- Scale maintenance by biomass and apply the temperature correction. --------!
         cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico)* maintenance_temp_dep
         cpatch%root_maintenance(ico) = cpatch%root_maintenance(ico)* maintenance_temp_dep
         cpatch%bark_maintenance(ico) = cpatch%bark_maintenance(ico)* maintenance_temp_dep
         !---------------------------------------------------------------------------------!

      case (3)
         !---------------------------------------------------------------------------------!
         !      Light phenology.  Leaf turnover rate will be adjusted according to the     !
         ! amplitude that comes from the dependence on the radiation (running average).    !
         ! Roots are the same as the other plants that don't depend on temperature.        !
         !---------------------------------------------------------------------------------!
         cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico)                       &
                                      * cpatch%turnover_amp(ico)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      
      return
   end subroutine get_maintenance
   !=======================================================================================!
   !=======================================================================================!

   
   
   !=======================================================================================!
   !=======================================================================================!
   subroutine apply_maintenance(cpatch,ico,tfact,cb_decrement)
      use ed_state_vars, only : patchtype             ! ! structure
      use pft_coms     , only : storage_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC            & ! intent(in)
                              , day_sec               ! ! intent(in)
      use ed_misc_coms , only : storage_resp_scheme   ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(in)   :: tfact
      real           , intent(out)  :: cb_decrement
      !----- Local Vars -------------------------------------------------------------------!
      real                          :: rresp
      !------------------------------------------------------------------------------------!
      cb_decrement = 0.0
      
      !------------------------------------------------------------------------------------!
      ! Apply the standard update.                                                         !
      !------------------------------------------------------------------------------------!
                                    
      cpatch%bleaf (ico) = cpatch%bleaf(ico) - cpatch%leaf_maintenance(ico)
      cpatch%broot (ico) = cpatch%broot(ico) - cpatch%root_maintenance(ico)
      cpatch%bbark (ico) = cpatch%bbark(ico) - cpatch%bark_maintenance(ico)
      cpatch%balive(ico) = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bsapwooda(ico)   &
                         + cpatch%bsapwoodb(ico) + cpatch%bbark(ico)
      cb_decrement       = cpatch%leaf_maintenance(ico) + cpatch%root_maintenance(ico)     &
                         + cpatch%bark_maintenance(ico)
      !------------------------------------------------------------------------------------!
      
   end subroutine apply_maintenance
   !=======================================================================================!
   !=======================================================================================!
   


   !=======================================================================================!
   !=======================================================================================!
   subroutine get_daily_C_gain(cpatch,ico,daily_C_gain)
      use ed_state_vars, only : patchtype          ! ! structure
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(out)  :: daily_C_gain
      !------------------------------------------------------------------------------------!

      if(cpatch%nplant(ico) > tiny(1.0)) then
         daily_C_gain = umol_2_kgC * day_sec * ( cpatch%today_gpp(ico)                     &
                                               - cpatch%today_leaf_resp(ico)               &
                                               - cpatch%today_root_resp(ico))              &
                                             / cpatch%nplant(ico)
      else
         daily_C_gain = 0.0
      end if

      end subroutine get_daily_C_gain
   !=======================================================================================!
   !=======================================================================================!
   

   
   
   !=======================================================================================!
   !=======================================================================================!
   subroutine update_cb(cpatch,ico,cb_decrement)
      use ed_state_vars, only : patchtype             ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      real           , intent(in)   :: cb_decrement
      !------------------------------------------------------------------------------------!
   
      cpatch%cb         (13,ico) = cpatch%cb(13,ico)          - cb_decrement
      cpatch%cb_lightmax(13,ico) = cpatch%cb_lightmax(13,ico) - cb_decrement
      cpatch%cb_moistmax(13,ico) = cpatch%cb_moistmax(13,ico) - cb_decrement
      cpatch%cb_mlmax   (13,ico) = cpatch%cb_mlmax   (13,ico) - cb_decrement


   end subroutine update_cb
   !=======================================================================================!
   !=======================================================================================!
   

   
   !=======================================================================================!
   !=======================================================================================!
   subroutine update_today_npp_vars(cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb &
                                   ,tr_bbark,carbon_balance)
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
      real           , intent(in)   :: tr_bbark
      real           , intent(in)   :: carbon_balance
      !----- Local variables. -------------------------------------------------------------!
      real                          :: tr_bsapwood
      !------------------------------------------------------------------------------------!
      
      tr_bsapwood = tr_bsapwooda + tr_bsapwoodb

      !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
      cpatch%today_nppleaf   (ico) = max(tr_bleaf    * cpatch%nplant(ico), 0.0)
      cpatch%today_nppfroot  (ico) = max(tr_broot    * cpatch%nplant(ico), 0.0)
      cpatch%today_nppsapwood(ico) = max(tr_bsapwood * cpatch%nplant(ico), 0.0)
      cpatch%today_nppbark   (ico) = max(tr_bbark    * cpatch%nplant(ico), 0.0)
      
      cpatch%today_nppdaily(ico)  = carbon_balance * cpatch%nplant(ico)
      !------------------------------------------------------------------------------!


   end subroutine update_today_npp_vars
   !=======================================================================================!
   !=======================================================================================!
   
   
   


   !=======================================================================================!
   !=======================================================================================!
   subroutine update_nitrogen(flushing,ipft,carbon_balance,nplant,tr_bleaf,tr_broot        &
                              ,tr_bstorage,tr_bbark,nitrogen_uptake,fsn_in)
      use consts_coms   , only : tiny_num                 ! ! intent(in)
      use decomp_coms   , only : f_labile                 ! ! intent(in)
      use pft_coms      , only : c2n_storage              & ! intent(in)
                               , c2n_leaf                 & ! intent(in)
                               , c2n_stem                 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      logical        , intent(in)    :: flushing
      integer        , intent(in)    :: ipft
      real           , intent(in)    :: carbon_balance
      real           , intent(in)    :: nplant
      real           , intent(in)    :: tr_bleaf
      real           , intent(in)    :: tr_broot
      real           , intent(in)    :: tr_bstorage
      real           , intent(in)    :: tr_bbark
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: fsn_in
      !----- Local variables. -------------------------------------------------------------!
      real                           :: n2c_labile
      !------------------------------------------------------------------------------------!
      
      n2c_labile = (f_labile(ipft)/c2n_leaf(ipft) + (1.0-f_labile(ipft))/c2n_stem(ipft))

      if (flushing) then
         !---------------------------------------------------------------------------------!
         !     Check whether there we're adding to or taking from storage.                 !
         !---------------------------------------------------------------------------------!
         if (tr_bstorage <= 0.0)  then
            !------------------------------------------------------------------------------!
            ! We are using all of daily C gain and some of bstorage.                       !
            ! Calculate N demand from using daily C gain.                                  !
            !------------------------------------------------------------------------------!
            if (carbon_balance < 0.0) then
               nitrogen_uptake = nitrogen_uptake + carbon_balance / c2n_storage
               nitrogen_uptake = nitrogen_uptake                                           &
                               + (carbon_balance - tr_bstorage)                            &
                               * ( f_labile(ipft) / c2n_leaf(ipft)                         &
                                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft)                 &
                                 -  1.0 / c2n_storage)
               
            else
               nitrogen_uptake = nitrogen_uptake + carbon_balance                          &
                              * ( f_labile(ipft) / c2n_leaf(ipft)                          &
                                + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )

               !---------------------------------------------------------------------------!
               ! Calculate additional N uptake from transfer of C from storage to balive.  !
               !---------------------------------------------------------------------------!
               nitrogen_uptake  = nitrogen_uptake +  ( - 1.0 * tr_bstorage )               &
                                * ( f_labile(ipft)  / c2n_leaf(ipft)                       &
                                  + (1.0 - f_labile(ipft)) / c2n_stem(ipft)                &
                                  -  1.0 / c2n_storage)
            end if

         else
            !------------------------------------------------------------------------------!
            !     N uptake for fraction of daily C gain going to balive.                   !
            !------------------------------------------------------------------------------!
            nitrogen_uptake = nitrogen_uptake + (carbon_balance - tr_bstorage)             &
                            * ( f_labile(ipft) / c2n_leaf(ipft)                            &
                              + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
            !----------- N uptake for fraction of daily C gain going to bstorage. ---------!
            nitrogen_uptake = nitrogen_uptake + tr_bstorage / c2n_storage
         end if
         !---------------------------------------------------------------------------------!
      else
         if (tr_bstorage > 0.0) then
            nitrogen_uptake = nitrogen_uptake + tr_bstorage / c2n_storage
         else
            !------------------------------------------------------------------------------!
            !    MLO: Is this correct? It seems that most places treat leaves, fine roots  !
            ! and sapwood with the same f_labile (I was planning to do the same for bark   !
            ! although I can see bark and sapwood having different fractions).             !
            !------------------------------------------------------------------------------!
            fsn_in = fsn_in                                                                &
                   - tr_bstorage           / c2n_storage    * nplant                       &
                   - (tr_bleaf + tr_broot) * n2c_labile     * nplant                       &
                   - tr_bbark              / c2n_stem(ipft) * nplant
            !------------------------------------------------------------------------------!
         end if
      end if      

      
   end subroutine update_nitrogen
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
   subroutine get_c_xfers(csite,ipa,ico,carbon_balance,green_leaf_factor,tr_bleaf,tr_broot &
                         ,tr_bsapwooda,tr_bsapwoodb,tr_bbark,tr_bstorage,carbon_debt       &
                         ,flushing,balive_aim)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : phenology    ! ! intent(in)
      use pft_coms      , only : q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , qbark        & ! intent(in)
                               , agf_bs       ! ! intent(in)
      use allometry     , only : size2bl      ! ! function
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use consts_coms   , only : tiny_num     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite               !< Current Site
      integer        , intent(in)    :: ipa                 !< Loop-Current Patch
      integer        , intent(in)    :: ico                 !< Loop-Current Cohort
      real           , intent(in)    :: carbon_balance      !< Plant net carbon uptake
      real           , intent(in)    :: green_leaf_factor   !< Cohort leaf-age param.
      real           , intent(out)   :: tr_bleaf            !< Transfer to leaf C pool
      real           , intent(out)   :: tr_broot            !< Transfer to root C pool
      real           , intent(out)   :: tr_bsapwooda        !< Transfer to sapwooda C pool
      real           , intent(out)   :: tr_bsapwoodb        !< Transfer to sapwoodb C pool
      real           , intent(out)   :: tr_bbark            !< Transfer to bark C pool
      real           , intent(out)   :: tr_bstorage         !< Transfer to storage C pool
      real           , intent(out)   :: carbon_debt         !< Net cohort carbon uptake
      logical        , intent(out)   :: flushing            !< Flag for leaf flush
      real           , intent(out)   :: balive_aim          !< Desired cohort balive value
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bleaf_aim
      real                           :: broot_aim
      real                           :: bsapwooda_aim
      real                           :: bsapwoodb_aim
      real                           :: bbark_aim
      real                           :: bleaf_max
      real                           :: balive_max
      real                           :: bloss_max
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwooda
      real                           :: delta_bsapwoodb
      real                           :: delta_bbark
      real                           :: available_carbon
      real                           :: increment
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwooda
      real                           :: f_bsapwoodb
      real                           :: f_bbark
      real                           :: f_bstorage
      logical                        :: on_allometry
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
      !      write (unit=66,fmt='(21(a,1x))')                                               &
      !        ,'        YEAR','       MONTH','         DAY','         PFT','   PHENOLOGY'  &
      !        ,'PHEN_STAT_IN','PHN_STAT_OUT','  FLUSH_TIME',' AVAILABLE_C','      ELONGF'  &
      !        ,'  GREEN_LEAF','    ON_ALLOM',' DELTA_BLEAF',' DELTA_BROOT','   DELTA_BSA'  &
      !        ,'   DELTA_BSB',' DELTA_BBARK','    TR_BLEAF','    TR_BROOT','      TR_BSA'  &
      !        ,'      TR_BSB'
      !      close (unit=66,status='keep')
      !   end if
      !   first_time = .false.
      !end if
      !------------------------------------------------------------------------------------!

      tr_bleaf     = 0.0
      tr_broot     = 0.0
      tr_bsapwooda = 0.0
      tr_bsapwoodb = 0.0
      tr_bbark     = 0.0
      tr_bstorage  = 0.0

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 
      phen_stat_in = cpatch%phenology_status(ico)
      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!
      available_carbon = cpatch%bstorage(ico) + carbon_balance
      time_to_flush    = carbon_balance > 0.0 .or.                                         &
                         ( available_carbon > 0.0 .and. cpatch%phenology_status(ico) == 1 )
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Maximum bleaf that the allometric relationship would allow.  If the plant is   !
      ! drought stressed (elongf < 1), we don't allow it to get back to full allometry.    !
      !------------------------------------------------------------------------------------!
      bleaf_max     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
      bleaf_aim     = bleaf_max * green_leaf_factor * cpatch%elongf(ico)
      broot_aim     = bleaf_aim * q    (ipft)
      bsapwooda_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * agf_bs(ipft)
      bsapwoodb_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * (1. - agf_bs(ipft))
      bbark_aim     = bleaf_aim * qbark(ipft) * cpatch%hite(ico)
      balive_aim    = bleaf_aim + broot_aim + bsapwooda_aim + bsapwoodb_aim + bbark_aim
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Check whether to increase living tissue biomass or not.                       !
      !------------------------------------------------------------------------------------!
      flushing = .false.
      if (time_to_flush) then 
         select case (cpatch%phenology_status(ico))
         case (0,1)
            flushing = .true.
            !------------------------------------------------------------------------------!
            !     There are leaves, we are not actively dropping leaves and we're off      !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!

            !---- Amount that bleaf, broot, and bsapwood are off allometry. ---------------!
            delta_bleaf     = max (0.0, bleaf_aim     - cpatch%bleaf    (ico))
            delta_broot     = max (0.0, broot_aim     - cpatch%broot    (ico))
            delta_bsapwooda = max (0.0, bsapwooda_aim - cpatch%bsapwooda(ico))
            delta_bsapwoodb = max (0.0, bsapwoodb_aim - cpatch%bsapwoodb(ico))
            delta_bbark     = max (0.0, bbark_aim     - cpatch%bbark    (ico))
            !------------------------------------------------------------------------------!

            if(cpatch%elongf(ico) < tiny_num) then
               
               write(*,'(a)')' ============================================'
               write(*,'(a)')' LINE 990 growth_balive.f90'
               write(*,'(a)')' subroutine alloc_plant_c_balance'
               write(*,'(a)')' '
               write(*,'(a)')' An elongation factor of effectively zero'
               write(*,'(a)')' has been detected during the transfer'
               write(*,'(a)')' of storage carbon back to active leaf pool.'
               write(*,'(a)')' This routine is expecting a non-zero '
               write(*,'(a)')' elongation as status leaves exist.'
               write(*,'(a)')' This is a minor bug that appears to trigger'
               write(*,'(a)')' in rare cases when veg dynamics are off and'
               write(*,'(a)')' drought stress is high.'
               write(*,'(a)')' '
               write(*,'(a)')' Continuing with 0 storage transfer.'
               write(*,'(a)')' ============================================'

               f_total=0.0
            else

               !---------------------------------------------------------------------------!
               !     If the available carbon is less than what we need to get back to      !
               ! allometry.  Grow pools in proportion to demand.  If we have enough        !
               ! carbon, we'll put the extra into bstorage.                                !
               !---------------------------------------------------------------------------!
               f_bleaf     = delta_bleaf     / bleaf_aim
               f_broot     = delta_broot     / broot_aim
               f_bsapwooda = delta_bsapwooda / bsapwooda_aim
               f_bsapwoodb = delta_bsapwoodb / bsapwoodb_aim
               f_bbark     = delta_bbark     / bbark_aim
               f_total     = f_bleaf + f_broot + f_bsapwooda + f_bsapwoodb + f_bbark
               !---------------------------------------------------------------------------!
            end if

            !------------------------------------------------------------------------------!
            !     We only allow transfer from storage to living tissues if there is need   !
            ! to transfer.                                                                 !
            !------------------------------------------------------------------------------!
            if (f_total > 0.0) then
               tr_bleaf     = min(delta_bleaf    , f_bleaf     / f_total * available_carbon)
               tr_broot     = min(delta_broot    , f_broot     / f_total * available_carbon)
               tr_bsapwooda = min(delta_bsapwooda, f_bsapwooda / f_total * available_carbon)
               tr_bsapwoodb = min(delta_bsapwoodb, f_bsapwoodb / f_total * available_carbon)
               tr_bbark     = min(delta_bbark    , f_bbark     / f_total * available_carbon)
            end if
            !------------------------------------------------------------------------------!

            tr_bstorage = carbon_balance - tr_bleaf - tr_broot - tr_bsapwooda              &
                        - tr_bsapwoodb - tr_bbark
         case default
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            tr_bstorage  = carbon_balance
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
         carbon_debt = -carbon_balance
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            !    Plants should be growing or at their maximum, first we try to take all    !
            ! the carbon needed from storage.                                              !
            !------------------------------------------------------------------------------!
            if (cpatch%bstorage(ico) > carbon_debt) then
               tr_bstorage = carbon_balance
            else
               !---------------------------------------------------------------------------!
               !     Not enough carbon in storage.  Take everything then start destroying  !
               ! tissues.                                                                  !
               !---------------------------------------------------------------------------!
               carbon_debt = carbon_debt - cpatch%bstorage(ico)
               tr_bstorage = -1.0*cpatch%bstorage(ico)

               !---------------------------------------------------------------------------!
               !     Find total biomass that can be lost.  We take an amount proportional  !
               ! to the current biomass of each the pools.                                 !
               !---------------------------------------------------------------------------!
               bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bbark(ico)
               f_bleaf     = cpatch%bleaf    (ico) / bloss_max
               f_broot     = cpatch%broot    (ico) / bloss_max
               f_bbark     = cpatch%bbark    (ico) / bloss_max

               if (bloss_max > carbon_debt) then
                  !----- Remove biomass accordingly. --------------------------------------!
                  tr_bleaf = f_bleaf * -1.0 * carbon_debt
                  tr_broot = f_broot * -1.0 * carbon_debt
                  tr_bbark = f_bbark * -1.0 * carbon_debt
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort didn't know how to save carbon during its life, and    !
                  ! has spent everything it had and now it is sunk in huge debt that it    !
                  ! can't afford.  It is with profound sadness that we announce that this  !
                  ! cohort is going to fertilizer business.                                !
                  !------------------------------------------------------------------------!
                  carbon_debt = bloss_max
                  tr_bleaf = -1.0*cpatch%bleaf(ico);
                  tr_broot = -1.0*cpatch%broot(ico);
                  tr_bbark = -1.0*cpatch%bbark(ico);
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
            bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bbark(ico)
            f_bleaf     = cpatch%bleaf    (ico) / bloss_max
            f_broot     = cpatch%broot    (ico) / bloss_max
            f_bbark     = cpatch%bbark    (ico) / bloss_max

            if (bloss_max > carbon_debt) then
               !----- Remove biomass accordingly. -----------------------------------------!
               tr_bleaf = f_bleaf * -1.0 * carbon_debt
               tr_broot = f_broot * -1.0 * carbon_debt
               tr_bbark = f_bbark * -1.0 * carbon_debt
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Not enough biomass, remove everything.                                !
               !---------------------------------------------------------------------------!
               carbon_debt = carbon_debt - bloss_max
               tr_bleaf = -1.0*cpatch%bleaf(ico);
               tr_broot = -1.0*cpatch%broot(ico);
               tr_bbark = -1.0*cpatch%bbark(ico);
               !---------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------!
               !     The living tissues weren't enough to meet the demand, remove what is  !
               ! still needed from the storage.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%bstorage(ico) > carbon_debt) then
                  !----- Enough carbon in storage, take all carbon needed from there. -----!
                  tr_bstorage = -1.0 * carbon_debt
                  !cpatch%bstorage(ico) = cpatch%bstorage(ico) - carbon_debt
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort didn't know how to save carbon during its life, and    !
                  ! has spent everything it had and now it is sunk in huge debt that it    !
                  ! can't afford.  It is with profound sadness that we announce that this  !
                  ! cohort is going to fertilizer business.                                !
                  !------------------------------------------------------------------------!
                  tr_bstorage = -1.0*cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!
               end if
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end if

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
   subroutine apply_c_xfers(cpatch,ico,carbon_balance,tr_bleaf,tr_broot,tr_bsapwooda       &
                           ,tr_bsapwoodb,tr_bbark,tr_bstorage)
      use ed_state_vars , only : patchtype  ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: carbon_balance
      real           , intent(in)    :: tr_bleaf
      real           , intent(in)    :: tr_broot
      real           , intent(in)    :: tr_bsapwooda
      real           , intent(in)    :: tr_bsapwoodb
      real           , intent(in)    :: tr_bbark
      real           , intent(in)    :: tr_bstorage
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update the carbon pools of the living tissues.                                 !
      !------------------------------------------------------------------------------------!
      cpatch%bleaf    (ico) = cpatch%bleaf    (ico) + tr_bleaf
      cpatch%broot    (ico) = cpatch%broot    (ico) + tr_broot
      cpatch%bsapwooda(ico) = cpatch%bsapwooda(ico) + tr_bsapwooda
      cpatch%bsapwoodb(ico) = cpatch%bsapwoodb(ico) + tr_bsapwoodb
      cpatch%bbark    (ico) = cpatch%bbark    (ico) + tr_bbark

      cpatch%balive   (ico) = cpatch%bleaf(ico)     + cpatch%broot(ico)                    &
                            + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)                &
                            + cpatch%bbark(ico)
      !------------------------------------------------------------------------------------!


      !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------------!
      cpatch%today_nppleaf(ico)    = tr_bleaf       * cpatch%nplant(ico)
      cpatch%today_nppfroot(ico)   = tr_broot       * cpatch%nplant(ico)
      cpatch%today_nppsapwood(ico) = (tr_bsapwooda + tr_bsapwoodb)* cpatch%nplant(ico)
      cpatch%today_nppbark(ico)    = tr_bbark       * cpatch%nplant(ico)
      cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
      !------------------------------------------------------------------------------------!
      

      !------------------------------------------------------------------------------------!
      !    Find the amount of carbon used to recover the tissues that were off-allometry,  !
      ! take that from the carbon balance first, then use some of the storage if needed    !
      ! be.                                                                                !
      !------------------------------------------------------------------------------------!
      cpatch%bstorage(ico) = max(0.0, cpatch%bstorage(ico) + tr_bstorage)
      !------------------------------------------------------------------------------------!

      return
   end subroutine apply_c_xfers
   !=======================================================================================!
   !=======================================================================================!
   
   
   


   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_maintenance(cpatch,ico,br,bl,bb,tfact,daily_C_gain,tempk)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : phenology          & ! intent(in)
                              , leaf_turnover_rate & ! intent(in)
                              , root_turnover_rate & ! intent(in)
                              , bark_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            & ! intent(in)
                              , tiny_num           ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: br
      real           , intent(in)    :: bl
      real           , intent(in)    :: bb
      real           , intent(in)    :: tfact
      real           , intent(in)    :: tempk
      real           , intent(out)   :: daily_C_gain
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: maintenance_temp_dep
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      !      Find the maintenance costs.  This will depend on the type of phenology that   !
      ! the PFT has.   The tfact term applied converts the maintenance rates to            !
      ! [kgC/plant/day].                                                                   !
      !------------------------------------------------------------------------------------!
      select case (phenology(ipft))
      case (0)
         !---------------------------------------------------------------------------------!
         !     Evergreens, like pines.  The turnover rates will be adjusted by a function  !
         ! of temperature, which approaches 0 as the temperature goes down.                !
         !---------------------------------------------------------------------------------!
         !------ Find a temperature dependence adjustment. --------------------------------!
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
         !----- Scale maintenance by biomass and apply the temperature correction. --------!
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl                      &
                                      * maintenance_temp_dep * tfact
         cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br                      &
                                      * maintenance_temp_dep * tfact
         cpatch%bark_maintenance(ico) = bark_turnover_rate(ipft) * bb                      &
                                      * maintenance_temp_dep * tfact
         !---------------------------------------------------------------------------------!

      case (3)
         !---------------------------------------------------------------------------------!
         !      Light phenology.  Leaf turnover rate will be adjusted according to the     !
         ! amplitude that comes from the dependence on the radiation (running average).    !
         ! Roots are the same as the other plants that don't depend on temperature.        !
         !---------------------------------------------------------------------------------!
         cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br * tfact
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl                      &
                                      * cpatch%turnover_amp(ico) * tfact
         cpatch%bark_maintenance(ico) = bark_turnover_rate(ipft) * bb * tfact
         !---------------------------------------------------------------------------------!


      case default
         !---------------------------------------------------------------------------------!
         !      Ohter phenologies, use the standard turnover rates, scaled by biomass      !
         ! only.                                                                           !
         !---------------------------------------------------------------------------------!
         cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br * tfact
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl * tfact
         cpatch%bark_maintenance(ico) = bark_turnover_rate(ipft) * bb * tfact
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      In case maintenance is too small, flush values to zero.                       !
      !------------------------------------------------------------------------------------!
      if (cpatch%leaf_maintenance(ico) < tiny_num) cpatch%leaf_maintenance(ico) = 0.0
      if (cpatch%root_maintenance(ico) < tiny_num) cpatch%root_maintenance(ico) = 0.0
      if (cpatch%bark_maintenance(ico) < tiny_num) cpatch%bark_maintenance(ico) = 0.0
      !------------------------------------------------------------------------------------!



      !----- Compute daily C uptake [kgC/plant/day]. --------------------------------------!
      if(cpatch%nplant(ico) > tiny(1.0)) then
         daily_C_gain = umol_2_kgC * day_sec * ( cpatch%today_gpp(ico)                     &
                                               - cpatch%today_leaf_resp(ico)               &
                                               - cpatch%today_root_resp(ico))              &
                                             / cpatch%nplant(ico)
      else
         daily_C_gain = 0.0
      end if

      return
   end subroutine plant_maintenance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance             &
                                   ,carbon_balance_pot,carbon_balance_lightmax             &
                                   ,carbon_balance_moistmax,carbon_balance_mlmax)
      use ed_state_vars  , only : patchtype          ! ! structure
      use pft_coms       , only : growth_resp_factor ! ! intent(in)
      use consts_coms    , only : umol_2_kgC         & ! intent(in)
                                , day_sec            ! ! intent(in)
      use ed_misc_coms   , only : current_time       ! ! intent(in)
      use ed_max_dims    , only : n_pft              ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: daily_C_gain
      real                     , intent(out) :: carbon_balance
      real                     , intent(out) :: carbon_balance_pot
      real                     , intent(out) :: carbon_balance_lightmax
      real                     , intent(out) :: carbon_balance_moistmax
      real                     , intent(out) :: carbon_balance_mlmax
      !----- Local variables. -------------------------------------------------------------!
      real                                   :: daily_C_gain_pot
      real                                   :: daily_C_gain_lightmax
      real                                   :: daily_C_gain_moistmax
      real                                   :: daily_C_gain_mlmax
      real                                   :: growth_respiration_pot
      real                                   :: growth_respiration_lightmax
      real                                   :: growth_respiration_moistmax
      real                                   :: growth_respiration_mlmax
      integer                                :: ipft
      !----- Local constants. -------------------------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !------------------------------------------------------------------------------------!

      !----- Alias for PFT type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Calculate actual daily carbon balance: kgC/plant/day.                        !
      !------------------------------------------------------------------------------------!
      carbon_balance = daily_C_gain - cpatch%leaf_growth_resp (ico)                        &
                                    - cpatch%root_growth_resp (ico)                        &
                                    - cpatch%sapa_growth_resp (ico)                        &
                                    - cpatch%sapb_growth_resp (ico)                        &
                                    - cpatch%bark_growth_resp (ico)
      !------------------------------------------------------------------------------------!

      if (cpatch%nplant(ico) > tiny(1.0)) then

         !---------------------------------------------------------------------------------!
         !      Calculate potential carbon balance (used for nitrogen demand function).    !
         ! [kgC/plant/day].                                                                !
         !---------------------------------------------------------------------------------!
         daily_C_gain_pot       = umol_2_kgC * day_sec * ( cpatch%today_gpp_pot(ico)       &
                                                         - cpatch%today_leaf_resp(ico)     &
                                                         - cpatch%today_root_resp(ico))    &
                                                       / cpatch%nplant(ico)
         growth_respiration_pot = max(0.0, daily_C_gain_pot * growth_resp_factor(ipft))
         carbon_balance_pot     = daily_C_gain_pot - growth_respiration_pot
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find carbon balance under full light and full soil moisture.  They will be !
         ! used for density-dependent mortality.  Units: [kgC/plant/day].                  !
         !---------------------------------------------------------------------------------!
         !------ Full light. --------------------------------------------------------------!
         daily_C_gain_lightmax       = umol_2_kgC * day_sec                                &
                                     * ( cpatch%today_gpp_lightmax(ico)                    &
                                       - cpatch%today_leaf_resp   (ico)                    &
                                       - cpatch%today_root_resp   (ico) )                  &
                                     / cpatch%nplant(ico)
         growth_respiration_lightmax = max(0.0, daily_C_gain_lightmax                      &
                                              * growth_resp_factor(ipft) )
         carbon_balance_lightmax     = daily_C_gain_lightmax - growth_respiration_lightmax
         !------ Full soil moisture. ------------------------------------------------------!
         daily_C_gain_moistmax       = umol_2_kgC * day_sec                                &
                                     * ( cpatch%today_gpp_moistmax(ico)                    &
                                       - cpatch%today_leaf_resp   (ico)                    &
                                       - cpatch%today_root_resp   (ico) )                  &
                                     / cpatch%nplant(ico)
         growth_respiration_moistmax = max(0.0, daily_C_gain_moistmax                      &
                                              * growth_resp_factor(ipft) )
         carbon_balance_moistmax     = daily_C_gain_moistmax - growth_respiration_moistmax
         !------ Full soil moisture and light. --------------------------------------------!
         daily_C_gain_mlmax          = umol_2_kgC * day_sec                                &
                                     * ( cpatch%today_gpp_mlmax(ico)                       &
                                       - cpatch%today_leaf_resp   (ico)                    &
                                       - cpatch%today_root_resp   (ico) )                  &
                                     / cpatch%nplant(ico)
         growth_respiration_mlmax    = max(0.0, daily_C_gain_mlmax                         &
                                              * growth_resp_factor(ipft) )
         carbon_balance_mlmax        = daily_C_gain_mlmax - growth_respiration_mlmax
         !---------------------------------------------------------------------------------!
      else
         carbon_balance_pot      = 0.0
         carbon_balance_lightmax = 0.0
         carbon_balance_moistmax = 0.0
         carbon_balance_mlmax    = 0.0
      end if

      !----- Carbon balances for mortality. -----------------------------------------------!
      cpatch%cb         (13,ico) = cpatch%cb         (13,ico) + carbon_balance
      cpatch%cb_lightmax(13,ico) = cpatch%cb_lightmax(13,ico) + carbon_balance_lightmax
      cpatch%cb_moistmax(13,ico) = cpatch%cb_moistmax(13,ico) + carbon_balance_moistmax
      cpatch%cb_mlmax   (13,ico) = cpatch%cb_mlmax   (13,ico) + carbon_balance_mlmax

      if (print_debug) then

         if (first_time(ipft)) then
            first_time(ipft) = .false.
            write (unit=30+ipft,fmt='(a10,25(1x,a18))')                                    &
               '      TIME','             PATCH','            COHORT','            NPLANT' &
                           ,'          CB_TODAY','  LEAF_GROWTH_RESP','  ROOT_GROWTH_RESP' &
                           ,'  SAPA_GROWTH_RESP','  SAPB_GROWTH_RESP','  BARK_GROWTH_RESP' &
                           ,'         TODAY_GPP','TODAY_GPP_LIGHTMAX','TODAY_GPP_MOISTMAX' &
                           ,'   TODAY_GPP_MLMAX','   TODAY_LEAF_RESP','   TODAY_ROOT_RESP' &
                           ,' CB_LIGHTMAX_TODAY',' CB_MOISTMAX_TODAY','    CB_MLMAX_TODAY' &
                           ,'                CB','       CB_LIGHTMAX','       CB_MOISTMAX' &
                           ,'          CB_MLMAX','  LEAF_MAINTENANCE','  ROOT_MAINTENANCE' &
                           ,'  BARK_MAINTENANCE'
         end if

         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i18),24(1x,es18.5))')               &
              current_time%month,'/',current_time%date,'/',current_time%year               &
             ,ipa,ico,cpatch%nplant(ico),carbon_balance,cpatch%leaf_growth_resp(ico)       &
             ,cpatch%root_growth_resp(ico),cpatch%sapa_growth_resp(ico)                    &
             ,cpatch%sapb_growth_resp(ico),cpatch%bark_growth_resp(ico)                    &
             ,cpatch%today_gpp(ico),cpatch%today_gpp_lightmax(ico)                         &
             ,cpatch%today_gpp_moistmax(ico),cpatch%today_gpp_mlmax(ico)                   &
             ,cpatch%today_leaf_resp(ico),cpatch%today_root_resp(ico)                      &
             ,carbon_balance_lightmax,carbon_balance_moistmax,carbon_balance_mlmax         &
             ,cpatch%cb(13,ico),cpatch%cb_lightmax(13,ico),cpatch%cb_moistmax(13,ico)      &
             ,cpatch%cb_mlmax(13,ico),cpatch%leaf_maintenance(ico)                         &
             ,cpatch%root_maintenance(ico),cpatch%bark_maintenance(ico)
      end if

      return
   end subroutine plant_carbon_balances
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_plant_c_balance(csite,ipa,ico,carbon_balance,nitrogen_uptake           &
                                   ,green_leaf_factor)
      use ed_state_vars , only : sitetype                 & ! structure
                               , patchtype                ! ! structure
      use pft_coms      , only : phenology                & ! intent(in)
                               , c2n_storage              & ! intent(in)
                               , c2n_leaf                 & ! intent(in)
                               , sla                      & ! intent(in)
                               , q                        & ! intent(in)
                               , qsw                      & ! intent(in)
                               , qbark                    & ! intent(in)
                               , agf_bs                   & ! intent(in)
                               , c2n_stem                 ! ! intent(in)
      use decomp_coms   , only : f_labile                 ! ! intent(in)
      use allometry     , only : size2bl                  ! ! function
      use phenology_coms, only : elongf_min               ! ! intent(in)
      use consts_coms   , only : tiny_num                 ! ! intent(in)
      use ed_misc_coms  , only : current_time             ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bleaf_aim
      real                           :: broot_aim
      real                           :: bsapwooda_aim
      real                           :: bsapwoodb_aim
      real                           :: bbark_aim
      real                           :: balive_aim
      real                           :: bleaf_max
      real                           :: balive_max
      real                           :: bloss_max
      real                           :: carbon_debt
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwooda
      real                           :: delta_bsapwoodb
      real                           :: delta_bbark
      real                           :: available_carbon
      real                           :: increment
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwooda
      real                           :: f_bsapwoodb
      real                           :: f_bbark
      real                           :: f_bstorage
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwooda
      real                           :: tr_bsapwoodb
      real                           :: tr_bbark
      logical                        :: on_allometry
      logical                        :: time_to_flush
      integer                        :: phen_stat_in 
      logical          , parameter   :: printout = .false.
      character(len=11), parameter   :: fracfile = 'cballoc.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save        :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then
         if (printout) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(22(a,1x))')                                               &
              ,'        YEAR','       MONTH','         DAY','         PFT','   PHENOLOGY'  &
              ,'PHEN_STAT_IN','PHN_STAT_OUT','  FLUSH_TIME',' AVAILABLE_C','      ELONGF'  &
              ,'  GREEN_LEAF','    ON_ALLOM',' DELTA_BLEAF',' DELTA_BROOT','   DELTA_BSA'  &
              ,'   DELTA_BSB',' DELTA_BBARK','    TR_BLEAF','    TR_BROOT','      TR_BSA'  &
              ,'      TR_BSB','    TR_BBARK'
            close (unit=66,status='keep')
         end if
         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      phen_stat_in = cpatch%phenology_status(ico)
      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!
      available_carbon = cpatch%bstorage(ico) + carbon_balance
      time_to_flush    = carbon_balance > 0.0 .or.                                         &
                         ( available_carbon > 0.0 .and. cpatch%phenology_status(ico) == 1 )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Check whether to increase living tissue biomass or not.                       !
      !------------------------------------------------------------------------------------!
      if (time_to_flush) then 
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            !     There are leaves, we are not actively dropping leaves and we're off      !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Maximum bleaf that the allometric relationship would allow.  If the      !
            ! plant is drought stress (elongf < 1), we do not allow the plant to get back  !
            ! to full allometry.                                                           !
            !------------------------------------------------------------------------------!
            bleaf_max     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
            bleaf_aim     = bleaf_max * green_leaf_factor * cpatch%elongf(ico)
            broot_aim     = bleaf_aim * q(ipft)
            bsapwooda_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * agf_bs(ipft)
            bsapwoodb_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * (1. - agf_bs(ipft))
            bbark_aim     = bleaf_aim * qbark(ipft) * cpatch%hite(ico)
            balive_aim    = bleaf_aim + broot_aim + bsapwooda_aim + bsapwoodb_aim          &
                          + bbark_aim
            !---- Amount that bleaf, broot, and bsapwood are off allometry. ---------------!
            delta_bleaf     = max (0.0, bleaf_aim     - cpatch%bleaf    (ico))
            delta_broot     = max (0.0, broot_aim     - cpatch%broot    (ico))
            delta_bsapwooda = max (0.0, bsapwooda_aim - cpatch%bsapwooda(ico))
            delta_bsapwoodb = max (0.0, bsapwoodb_aim - cpatch%bsapwoodb(ico))
            delta_bbark     = max (0.0, bbark_aim     - cpatch%bbark    (ico))
            !------------------------------------------------------------------------------!

            if(cpatch%elongf(ico) < tiny_num) then
               
               write(*,'(a)')' ============================================'
               write(*,'(a)')' LINE 990 growth_balive.f90'
               write(*,'(a)')' subroutine alloc_plant_c_balance'
               write(*,'(a)')' '
               write(*,'(a)')' An elongation factor of effectively zero'
               write(*,'(a)')' has been detected during the transfer'
               write(*,'(a)')' of storage carbon back to active leaf pool.'
               write(*,'(a)')' This routine is expecting a non-zero '
               write(*,'(a)')' elongation as status leaves exist.'
               write(*,'(a)')' This is a minor bug that appears to trigger'
               write(*,'(a)')' in rare cases when veg dynamics are off and'
               write(*,'(a)')' drought stress is high.'
               write(*,'(a)')' '
               write(*,'(a)')' Continuing with 0 storage transfer.'
               write(*,'(a)')' ============================================'

               f_total=0.0
            else

               !---------------------------------------------------------------------------!
               !     If the available carbon is less than what we need to get back to      !
               ! allometry.  Grow pools in proportion to demand.  If we have enough        !
               ! carbon, we'll put the extra into bstorage.                                !
               !---------------------------------------------------------------------------!
               f_bleaf     = delta_bleaf     / bleaf_aim
               f_broot     = delta_broot     / broot_aim
               f_bsapwooda = delta_bsapwooda / bsapwooda_aim
               f_bsapwoodb = delta_bsapwoodb / bsapwoodb_aim
               f_bbark     = delta_bbark     / bbark_aim
               f_total     = f_bleaf + f_broot + f_bsapwooda + f_bsapwoodb + f_bbark
               !---------------------------------------------------------------------------!
            end if

            !------------------------------------------------------------------------------!
            !     We only allow transfer from storage to living tissues if there is need   !
            ! to transfer.                                                                 !
            !------------------------------------------------------------------------------!
            if (f_total > 0.0) then
               tr_bleaf     = min(delta_bleaf    , f_bleaf     / f_total * available_carbon)
               tr_broot     = min(delta_broot    , f_broot     / f_total * available_carbon)
               tr_bsapwooda = min(delta_bsapwooda, f_bsapwooda / f_total * available_carbon)
               tr_bsapwoodb = min(delta_bsapwoodb, f_bsapwoodb / f_total * available_carbon)
               tr_bbark     = min(delta_bbark    , f_bbark     / f_total * available_carbon)
            else
               tr_bleaf     = 0.
               tr_broot     = 0.
               tr_bsapwooda = 0.
               tr_bsapwoodb = 0.
               tr_bbark     = 0.
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Update the carbon pools of the living tissues.                           !
            !------------------------------------------------------------------------------!
            cpatch%bleaf    (ico) = cpatch%bleaf    (ico) + tr_bleaf
            cpatch%broot    (ico) = cpatch%broot    (ico) + tr_broot
            cpatch%bsapwooda(ico) = cpatch%bsapwooda(ico) + tr_bsapwooda
            cpatch%bsapwoodb(ico) = cpatch%bsapwoodb(ico) + tr_bsapwoodb
            cpatch%bbark    (ico) = cpatch%bbark    (ico) + tr_bbark

            cpatch%balive   (ico) = cpatch%bleaf(ico) + cpatch%broot(ico)                  &
                                  + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)          &
                                  + cpatch%bbark(ico)
            !------------------------------------------------------------------------------!


            !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
            cpatch%today_nppleaf(ico)    = tr_bleaf       * cpatch%nplant(ico)
            cpatch%today_nppfroot(ico)   = tr_broot       * cpatch%nplant(ico)
            cpatch%today_nppsapwood(ico) = (tr_bsapwooda + tr_bsapwoodb)* cpatch%nplant(ico)
            cpatch%today_nppbark(ico)    = tr_bbark       * cpatch%nplant(ico)
            cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
            !------------------------------------------------------------------------------!
            

            !------------------------------------------------------------------------------!
            !    Find the amount of carbon used to recover the tissues that were off-      !
            ! -allometry, take that from the carbon balance first, then use some of the    !
            ! storage if needed be.                                                        !
            !------------------------------------------------------------------------------!
            increment = carbon_balance -  tr_bleaf - tr_broot - tr_bsapwooda               &
                      - tr_bsapwoodb - tr_bbark
            cpatch%bstorage(ico) = max(0.0, cpatch%bstorage(ico) + increment)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check whether there is still some carbon to go to storage or if we are   !
            ! burning some of the storage.                                                 !
            !------------------------------------------------------------------------------!
            if (increment <= 0.0)  then
               !---------------------------------------------------------------------------!
               !    We are using up all of daily C gain and some of bstorage.  First       !
               ! calculate N demand from using daily C gain.                               !
               !---------------------------------------------------------------------------!
               if (carbon_balance < 0.0) then
                  nitrogen_uptake = nitrogen_uptake + carbon_balance / c2n_storage
                  nitrogen_uptake = nitrogen_uptake                                        &
                                  + (carbon_balance - increment)                           &
                                  * ( f_labile(ipft) / c2n_leaf(ipft)                      &
                                    + (1.0 - f_labile(ipft)) / c2n_stem(ipft)              &
                                    -  1.0 / c2n_storage)
                  
               else
                  nitrogen_uptake = nitrogen_uptake + carbon_balance                       &
                                 * ( f_labile(ipft) / c2n_leaf(ipft)                       &
                                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )

                  !------------------------------------------------------------------------!
                  !     Now calculate additional N uptake required from transfer of C from !
                  ! storage to balive.                                                     !
                  !------------------------------------------------------------------------!
                  nitrogen_uptake  = nitrogen_uptake +  ( - 1.0 * increment )              &
                                   * ( f_labile(ipft)  / c2n_leaf(ipft)                    &
                                     + (1.0 - f_labile(ipft)) / c2n_stem(ipft)             &
                                     -  1.0 / c2n_storage)
               end if

            else
               !---------------------------------------------------------------------------!
               !     N uptake for fraction of daily C gain going to balive.                !
               !---------------------------------------------------------------------------!
               nitrogen_uptake = nitrogen_uptake + (carbon_balance - increment)            &
                               * ( f_labile(ipft) / c2n_leaf(ipft)                         &
                                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
               !----- N uptake for fraction of daily C gain going to bstorage. ------------!
               nitrogen_uptake = nitrogen_uptake + increment / c2n_storage
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check whether we are on allometry or not.                                !
            !------------------------------------------------------------------------------!
            on_allometry = (balive_aim - cpatch%balive(ico)) <= 0.000001*balive_aim
            if (cpatch%elongf(ico) == 1.0 .and. on_allometry) then
               !---------------------------------------------------------------------------!
               !     We're back to allometry, change phenology_status.                     !
               !---------------------------------------------------------------------------!
               cpatch%phenology_status(ico) = 0
            end if
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
            nitrogen_uptake      = nitrogen_uptake      + carbon_balance / c2n_storage
                                 
            !----- NPP allocation in diff pools in Kg C/m2/day. ---------------------------!
            cpatch%today_nppleaf(ico)    = 0.0
            cpatch%today_nppfroot(ico)   = 0.0
            cpatch%today_nppsapwood(ico) = 0.0
            cpatch%today_nppbark(ico)    = 0.0
            cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
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
         carbon_debt = -carbon_balance
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            !    Plants should be growing or at their maximum, first we try to take all    !
            ! the carbon needed from storage.                                              !
            !------------------------------------------------------------------------------!
            if (cpatch%bstorage(ico) > carbon_debt) then
               !----- Enough carbon in storage, take all carbon needed from there. --------!
               cpatch%bstorage(ico) = cpatch%bstorage(ico) - carbon_debt
               csite%fsn_in(ipa)    = csite%fsn_in(ipa)                                    &
                                    + carbon_debt / c2n_storage * cpatch%nplant(ico)
            else
               !---------------------------------------------------------------------------!
               !     Not enough carbon in storage.  Take everything then start destroying  !
               ! tissues.                                                                  !
               !---------------------------------------------------------------------------!
               carbon_debt          = carbon_debt - cpatch%bstorage(ico)
               csite%fsn_in(ipa)    = csite%fsn_in(ipa) + cpatch%bstorage(ico)             &
                                                        / c2n_storage * cpatch%nplant(ico)
               cpatch%bstorage(ico) = 0.0

               !---------------------------------------------------------------------------!
               !     Find total biomass that can be lost.  We take an amount proportional  !
               ! to the current biomass of each the pools.  Bark is not a 100% living      !
               ! tissue, but for simplicity we don't separate inner and outer bark.        !
               !---------------------------------------------------------------------------!
               bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bbark(ico)
               f_bleaf     = cpatch%bleaf    (ico) / bloss_max
               f_broot     = cpatch%broot    (ico) / bloss_max
               f_bbark     = cpatch%bbark    (ico) / bloss_max

               if (bloss_max > carbon_debt) then
                  !----- Remove biomass accordingly. --------------------------------------!
                  cpatch%bleaf    (ico) = cpatch%bleaf    (ico) - f_bleaf     * carbon_debt
                  cpatch%broot    (ico) = cpatch%broot    (ico) - f_broot     * carbon_debt
                  cpatch%bbark    (ico) = cpatch%bbark    (ico) - f_bbark     * carbon_debt
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort didn't know how to save carbon during its life, and    !
                  ! has spent everything it had and now it is sunk in huge debt that it    !
                  ! can't afford.  It is with profound sadness that we announce that this  !
                  ! cohort is going to fertilizer business.                                !
                  !------------------------------------------------------------------------!
                  carbon_debt           = bloss_max
                  cpatch%bleaf    (ico) = 0.0
                  cpatch%broot    (ico) = 0.0
                  cpatch%bbark    (ico) = 0.0
                  !------------------------------------------------------------------------!
               end if
               cpatch%phenology_status(ico) = 1
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Update living biomass.                                                !
               !---------------------------------------------------------------------------!
               cpatch%balive(ico) = cpatch%bleaf(ico) + cpatch%broot(ico)                  &
                                  + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)          &
                                  + cpatch%bbark(ico)
               !---------------------------------------------------------------------------!




               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !                  NOT SURE IF THIS IS CORRECT N ACCOUNTING                 !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               csite%fsn_in(ipa)    = csite%fsn_in(ipa) + carbon_debt                      &
                                    * ( f_labile(ipft) / c2n_leaf(ipft)                    &
                                      + (1.0 - f_labile(ipft)) / c2n_stem(ipft))           &
                                    * cpatch%nplant(ico)   
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         case (-1,-2)
            !------------------------------------------------------------------------------!
            !      Plants were already shedding leaves.  We swap the order here and remove !
            ! living tissues first, and only if there is nothing left that we remove       !
            ! storage.                                                                     !
            !------------------------------------------------------------------------------!
            bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico)
            f_bleaf     = cpatch%bleaf    (ico) / bloss_max
            f_broot     = cpatch%broot    (ico) / bloss_max
            f_bbark     = cpatch%bbark    (ico) / bloss_max

            if (bloss_max > carbon_debt) then
               !----- Remove biomass accordingly. -----------------------------------------!
               cpatch%bleaf    (ico) = cpatch%bleaf    (ico) - f_bleaf     * carbon_debt
               cpatch%broot    (ico) = cpatch%broot    (ico) - f_broot     * carbon_debt
               cpatch%bbark    (ico) = cpatch%bbark    (ico) - f_bbark     * carbon_debt
               !---------------------------------------------------------------------------!



               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !                  NOT SURE IF THIS IS CORRECT N ACCOUNTING                 !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               csite%fsn_in(ipa)    = csite%fsn_in(ipa) + carbon_debt                      &
                                    * ( f_labile(ipft) / c2n_leaf(ipft)                    &
                                      + (1.0 - f_labile(ipft)) / c2n_stem(ipft))           &
                                    * cpatch%nplant(ico)   
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Not enough biomass, remove everything.                                !
               !---------------------------------------------------------------------------!
               carbon_debt           = carbon_debt - bloss_max
               cpatch%bleaf    (ico) = 0.0
               cpatch%broot    (ico) = 0.0
               cpatch%bbark    (ico) = 0.0
               !---------------------------------------------------------------------------!



               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !                  NOT SURE IF THIS IS CORRECT N ACCOUNTING                 !
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               csite%fsn_in(ipa)    = csite%fsn_in(ipa) + bloss_max                        &
                                    * ( f_labile(ipft) / c2n_leaf(ipft)                    &
                                      + (1.0 - f_labile(ipft)) / c2n_stem(ipft))           &
                                    * cpatch%nplant(ico)   
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     The living tissues weren't enough to meet the demand, remove what is  !
               ! still needed from the storage.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%bstorage(ico) > carbon_debt) then
                  !----- Enough carbon in storage, take all carbon needed from there. -----!
                  cpatch%bstorage(ico) = cpatch%bstorage(ico) - carbon_debt
                  csite%fsn_in(ipa)    = csite%fsn_in(ipa)                                 &
                                       + carbon_debt / c2n_storage * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     This cohort didn't know how to save carbon during its life, and    !
                  ! has spent everything it had and now it is sunk in huge debt that it    !
                  ! can't afford.  It is with profound sadness that we announce that this  !
                  ! cohort is going to fertilizer business.                                !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa)    = csite%fsn_in(ipa) + cpatch%bstorage(ico)          &
                                       / c2n_storage * cpatch%nplant(ico)
                  cpatch%bstorage(ico) = 0.0
                  !------------------------------------------------------------------------!
               end if
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Update living biomass.                                                   !
            !------------------------------------------------------------------------------!
            cpatch%balive(ico) = cpatch%bleaf(ico) + cpatch%broot(ico)                     &
                               + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)             &
                               + cpatch%bbark(ico)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !---- NPP allocation in diff pools in KgC/m2/day. --------------------------------!
         cpatch%today_nppleaf(ico)    = 0.0
         cpatch%today_nppfroot(ico)   = 0.0
         cpatch%today_nppsapwood(ico) = 0.0
         cpatch%today_nppbark (ico)   = 0.0
         cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
         !---------------------------------------------------------------------------------!
      end if

      !------------------------------------------------------------------------------------!
      if (printout) then
         open (unit=66,file=fracfile,status='old',position='append',action='write')
         write(unit=66                                                                     &
              ,fmt='(7(i12,1x),1(11x,l1,1x),3(f12.6,1x),1(11x,l1,1x),10(f12.8,1x))')       &
               current_time%year,current_time%month,current_time%date,ipft,phenology(ipft) &
              ,phen_stat_in,cpatch%phenology_status(ico),time_to_flush,available_carbon    &
              ,cpatch%elongf(ico),green_leaf_factor,on_allometry,delta_bleaf,delta_broot   &
              ,delta_bsapwooda,delta_bsapwoodb,delta_bbark,tr_bleaf,tr_broot,tr_bsapwooda  &
              ,tr_bsapwoodb,tr_bbark
         close (unit=66,status='keep')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine alloc_plant_c_balance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Alternative subroutine for grasses.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_plant_c_balance_grass(csite,ipa,ico,salloc,salloci,carbon_balance      &
                                   ,nitrogen_uptake,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , qbark        & ! intent(in)
                               , agf_bs       & ! intent(in)
                               , r_fract      & ! intent(in)
                               , hgt_max      & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : bl2h         & ! function
                               , h2dbh        & ! function
                               , dbh2h        ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max
      real                           :: balive_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: old_status
      real                           :: delta_balive
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwooda
      real                           :: delta_bsapwoodb
      real                           :: delta_bbark
      real                           :: available_carbon
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwooda
      real                           :: f_bsapwoodb
      real                           :: f_bbark
      real                           :: f_resp
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwooda
      real                           :: tr_bsapwoodb
      real                           :: tr_bbark
      real                           :: bl
      real                           :: dbh_to_height
      logical                        :: on_allometry
      logical                        :: time_to_flush
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!
      increment     = cpatch%bstorage(ico) + carbon_balance
      time_to_flush = (carbon_balance <= 0.0) .and. (increment > 0.0) .and.                &
                      (cpatch%phenology_status(ico) == 1) 

      if (carbon_balance > 0.0 .or. time_to_flush) then 
         if ((cpatch%hite(ico)*(1 + 1e-4)) < hgt_max(ipft)) then ! - could use repro_min_h here instead
            !------------------------------------------------------------------------------!
            ! The grass is in a vegetative growth phase, put carbon into growth.           !
            !------------------------------------------------------------------------------!
            !--allow grass to use carbon from that day and from storage to grow
            available_carbon = carbon_balance + cpatch%bstorage(ico)

            !--scale maximum growth by elongf (currently grass is "evergreen" so elongf=1)
            delta_balive = available_carbon * cpatch%elongf(ico)
            increment    = available_carbon * (1. - cpatch%elongf(ico))

            delta_bleaf     = delta_balive * salloci * green_leaf_factor
            delta_broot     = delta_balive * salloci * q  (ipft)
            delta_bsapwooda = delta_balive * salloci * qsw(ipft)   * cpatch%hite(ico)      &
                            * agf_bs(ipft)
            delta_bsapwoodb = delta_balive * salloci * qsw(ipft)   * cpatch%hite(ico)      &
                            * (1.-agf_bs(ipft))
            delta_bbark     = delta_balive * salloci * qbark(ipft) * cpatch%hite(ico)



            cpatch%bleaf(ico)     = cpatch%bleaf(ico)     + delta_bleaf
            cpatch%broot(ico)     = cpatch%broot(ico)     + delta_broot
            cpatch%bsapwooda(ico) = cpatch%bsapwooda(ico) + delta_bsapwooda
            cpatch%bsapwoodb(ico) = cpatch%bsapwoodb(ico) + delta_bsapwoodb
            cpatch%bbark(ico)     = cpatch%bbark(ico)     + delta_bbark

            cpatch%balive(ico)   = cpatch%bleaf(ico) + cpatch%broot(ico)                   &
                                 + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)           &
                                 + cpatch%bbark(ico)

            !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
            cpatch%today_nppleaf    (ico) = delta_bleaf       * cpatch%nplant(ico)
            cpatch%today_nppfroot   (ico) = delta_broot       * cpatch%nplant(ico)
            cpatch%today_nppsapwood (ico) = (delta_bsapwooda + delta_bsapwoodb)            &
                                                              * cpatch%nplant(ico)
            cpatch%today_nppbark    (ico) = delta_bbark       * cpatch%nplant(ico)
            cpatch%today_nppdaily   (ico) = carbon_balance    * cpatch%nplant(ico)

            !----- update height for grasses to match new leaf mass -----------------------!
            cpatch%hite(ico) = min(hgt_max(ipft), bl2h(cpatch%bleaf(ico), ipft))  !limit by maximum height
            cpatch%dbh(ico)  = h2dbh(cpatch%hite(ico), ipft) !--effective_dbh value for grasses
                
            
            !----- put remaining carbon in the storage pool -------------------------------!
            cpatch%bstorage(ico) = max(0.0, cpatch%bstorage(ico) + increment)
            !------------------------------------------------------------------------------!
            

            if (increment <= 0.0)  then
               !---------------------------------------------------------------------------!
               !    We are using up all of daily C gain and some of bstorage.  First       !
               ! calculate N demand from using daily C gain.                               !
               !---------------------------------------------------------------------------!

               nitrogen_uptake = nitrogen_uptake + carbon_balance                       &
                              * ( f_labile(ipft) / c2n_leaf(ipft)                       &
                                + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )

               !------------------------------------------------------------------------!
               !     Now calculate additional N uptake required from transfer of C from !
               ! storage to balive.                                                     !
               !------------------------------------------------------------------------!
               nitrogen_uptake  = nitrogen_uptake +  ( - 1.0 * increment )              &
                                * ( f_labile(ipft)  / c2n_leaf(ipft)                    &
                                  + (1.0 - f_labile(ipft)) / c2n_stem(ipft)             &
                                  -  1.0 / c2n_storage)

            else
               !---------------------------------------------------------------------------!
               !     N uptake for fraction of daily C gain going to balive.                !
               !---------------------------------------------------------------------------!
               nitrogen_uptake = nitrogen_uptake + (carbon_balance - increment)            &
                               * ( f_labile(ipft) / c2n_leaf(ipft)                         &
                                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
               !----- N uptake for fraction of daily C gain going to bstorage. ------------!
               nitrogen_uptake = nitrogen_uptake + increment / c2n_storage
            end if

         else  !-- plant is at the maximum height.
            !------------------------------------------------------------------------------!
            !     Grass is at its maximum height.  Put carbon gain into storage.           !
            !   For agriculture, put carbon into grain (and storage?)                      !
            !------------------------------------------------------------------------------!
            
            !--test here if pft is agriculture, if so put most carbon into grain and ------!
            !- maybe a little into storage -- STILL TO BE WRITTEN! ------------------------!
            !----- Consider adding allocation to harvest pool here ------------------------!
            !---ALS=== Agriculture
            increment = carbon_balance ! subtract the part that goes into grain for ag here

            cpatch%bstorage(ico) = cpatch%bstorage(ico) + increment
            nitrogen_uptake      = nitrogen_uptake      + increment / c2n_storage
                                 
            
            
            !----- NPP allocation in diff pools in Kg C/m2/day. ---------------------------!
            cpatch%today_nppleaf(ico)    = 0.0
            cpatch%today_nppfroot(ico)   = 0.0
            cpatch%today_nppsapwood(ico) = 0.0
            cpatch%today_nppbark(ico)    = 0.0
            cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
         end if   !-- end height loop for cb>0


      else !-- carbon_balance <0
         !---------------------------------------------------------------------------------!
         !   Carbon balance is negative, take it out of storage.                           !
         !---------------------------------------------------------------------------------!
         increment =  cpatch%bstorage(ico) + carbon_balance
         


         if (increment <= 0.0)  then  !-- carbon_balance > bstorage
            !----- Take carbon from the Storage pool first --------------------------------!
            increment            =  - increment
            cpatch%bstorage(ico) = 0.0
            csite%fsn_in(ipa)    = csite%fsn_in(ipa) + cpatch%bstorage(ico) / c2n_storage  &
                                                     * cpatch%nplant(ico)

            !-- Take the remaining carbon from balive -------------------------------------!
            cpatch%balive(ico)    = cpatch%balive(ico) - increment
            cpatch%bleaf(ico)     = cpatch%balive(ico) * salloci * green_leaf_factor
            cpatch%broot(ico)     = cpatch%balive(ico) * salloci * q(ipft)
            cpatch%bsapwooda(ico) = cpatch%balive(ico) * salloci * qsw(ipft)               &
                                  * cpatch%hite(ico)   * agf_bs(ipft)
            cpatch%bsapwoodb(ico) = cpatch%balive(ico) * salloci * qsw(ipft)               &
                                  * cpatch%hite(ico)   * (1. - agf_bs(ipft))
            cpatch%bbark(ico)     = cpatch%balive(ico) * salloci * qbark(ipft)             &
                                  * cpatch%hite(ico)


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !                  NOT SURE IF THIS IS CORRECT N ACCOUNTING                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            csite%fsn_in(ipa) = csite%fsn_in(ipa) + increment                              &
                              * ( f_labile(ipft) / c2n_leaf(ipft)                          &
                                + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )                &
                              * cpatch%nplant(ico)

         else  !-- carbon_balance < bstorage
            !------ Remove the necessary carbon from the storage pool.  -------------------!
            !------ Dont' forget the nitrogen. --------------------------------------------!
            cpatch%bstorage(ico) =  cpatch%bstorage(ico) + carbon_balance
            csite%fsn_in(ipa)    = csite%fsn_in(ipa) - carbon_balance                      &
                                 * ( f_labile(ipft) / c2n_leaf(ipft)                       &
                                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft))              &
                                 * cpatch%nplant(ico)   
         end if

         !---- NPP allocation in diff pools in KgC/m2/day. --------------------------------!
         cpatch%today_nppleaf(ico)    = 0.0
         cpatch%today_nppfroot(ico)   = 0.0
         cpatch%today_nppsapwood(ico) = 0.0
         cpatch%today_nppbark(ico)    = 0.0
         cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
      end if

            
      return
   end subroutine alloc_plant_c_balance_grass
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     Find the allocation terms, but don't really allocate them.                        !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_plant_c_balance_eq_0(csite,ipa,ico,carbon_balance,nitrogen_uptake      &
                                        ,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , qbark        & ! intent(in)
                               , agf_bs       & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : size2bl      ! ! function
      use consts_coms   , only : tiny_num     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bleaf_aim
      real                           :: broot_aim
      real                           :: bsapwooda_aim
      real                           :: bsapwoodb_aim
      real                           :: bbark_aim
      real                           :: balive_aim
      real                           :: bleaf_max
      real                           :: increment
      real                           :: carbon_loss
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwooda
      real                           :: delta_bsapwoodb
      real                           :: delta_bbark
      real                           :: available_carbon
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwooda
      real                           :: f_bsapwoodb
      real                           :: f_bbark
      real                           :: f_resp
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwooda
      real                           :: tr_bsapwoodb
      real                           :: tr_bbark
      logical                        :: on_allometry
      logical                        :: time_to_flush
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!
      increment     = cpatch%bstorage(ico) + carbon_balance
      time_to_flush = (carbon_balance <= 0.0) .and. (increment > 0.0) .and.                &
                      (cpatch%phenology_status(ico) == 1) 

      if (carbon_balance > 0.0 .or. time_to_flush) then 
         select case (cpatch%phenology_status(ico))
         case (0,1)
            !------------------------------------------------------------------------------!
            ! There are leaves, we are not actively dropping leaves and we're off          !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!
            available_carbon = cpatch%bstorage(ico) + carbon_balance

            !------------------------------------------------------------------------------!
            !     Maximum bleaf that the allometric relationship would allow.  If the      !
            ! plant is drought stress (elongf < 1), we do not allow the plant to get back  !
            ! to full allometry.                                                           !
            !------------------------------------------------------------------------------!
            bleaf_max     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
            bleaf_aim     = bleaf_max * green_leaf_factor * cpatch%elongf(ico)
            broot_aim     = bleaf_aim * q    (ipft) * cpatch%elongf(ico)
            bsapwooda_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * agf_bs(ipft)
            bsapwoodb_aim = bleaf_aim * qsw  (ipft) * cpatch%hite(ico) * (1. - agf_bs(ipft))
            bbark_aim     = bleaf_aim * qbark(ipft) * cpatch%hite(ico)
            balive_aim    = bleaf_aim + broot_aim + bsapwooda_aim + bsapwoodb_aim          &
                          + bbark_aim
            !---- Amount that bleaf, broot, and bsapwood are off allometry. ---------------!
            delta_bleaf     = max (0.0, bleaf_aim     - cpatch%bleaf    (ico))
            delta_broot     = max (0.0, broot_aim     - cpatch%broot    (ico))
            delta_bsapwooda = max (0.0, bsapwooda_aim - cpatch%bsapwooda(ico))
            delta_bsapwoodb = max (0.0, bsapwoodb_aim - cpatch%bsapwoodb(ico))
            delta_bbark     = max (0.0, bbark_aim     - cpatch%bbark    (ico))
            !------------------------------------------------------------------------------!

            if(cpatch%elongf(ico) < tiny_num) then

               write(*,'(a)')' ============================================'
               write(*,'(a)')' LINE 1660 growth_balive.f90'
               write(*,'(a)')' subroutine alloc_plant_c_balance_eq_0'
               write(*,'(a)')' '
               write(*,'(a)')' An elongation factor of effectively zero'
               write(*,'(a)')' has been detected during the transfer'
               write(*,'(a)')' of storage carbon back to active leaf pool.'
               write(*,'(a)')' This routine is expecting a non-zero '
               write(*,'(a)')' elongation as status is flushing.'
               write(*,'(a)')' This is a minor bug that appears to trigger'
               write(*,'(a)')' in rare cases when veg dynamics are off and'
               write(*,'(a)')' drought stress is high.'
               write(*,'(a)')' '
               write(*,'(a)')' Continuing with 0 storage transfer.'
               write(*,'(a)')' ============================================'

               f_total=0.0
            else
               
               !---------------------------------------------------------------------------!
               !     If the available carbon is less than what we need to get back to      !
               ! allometry.  Grow pools in proportion to demand.  If we have enough        !
               ! carbon, we'll put the extra into bstorage.                                !
               !---------------------------------------------------------------------------!
               f_bleaf     = delta_bleaf     / bleaf_aim
               f_broot     = delta_broot     / broot_aim
               f_bsapwooda = delta_bsapwooda / bsapwooda_aim
               f_bsapwoodb = delta_bsapwoodb / bsapwoodb_aim
               f_bbark     = delta_bbark     / bbark_aim
               f_total     = f_bleaf + f_broot + f_bsapwooda + f_bsapwoodb + f_bbark
               !---------------------------------------------------------------------------!
               
            end if

            !------------------------------------------------------------------------------!
            !     We only allow transfer from storage to living tissues if there is need   !
            ! to transfer.                                                                 !
            !------------------------------------------------------------------------------!
            if (f_total > 0.0) then
               tr_bleaf     = min(delta_bleaf    , (f_bleaf/f_total)     * available_carbon)
               tr_broot     = min(delta_broot    , (f_broot/f_total)     * available_carbon)
               tr_bsapwooda = min(delta_bsapwooda, (f_bsapwooda/f_total) * available_carbon)
               tr_bsapwoodb = min(delta_bsapwoodb, (f_bsapwoodb/f_total) * available_carbon)
               tr_bbark     = min(delta_bbark    , (f_bbark/f_total)     * available_carbon)
            else
               tr_bleaf     = 0.
               tr_broot     = 0.
               tr_bsapwooda = 0.
               tr_bsapwoodb = 0.
               tr_bbark     = 0.
            end if
            !------------------------------------------------------------------------------!

            cpatch%bleaf    (ico) = cpatch%bleaf    (ico) + tr_bleaf
            cpatch%broot    (ico) = cpatch%broot    (ico) + tr_broot
            cpatch%bsapwooda(ico) = cpatch%bsapwooda(ico) + tr_bsapwooda
            cpatch%bsapwoodb(ico) = cpatch%bsapwoodb(ico) + tr_bsapwoodb
            cpatch%bbark    (ico) = cpatch%bbark    (ico) + tr_bbark

            cpatch%balive(ico)    = cpatch%bleaf    (ico) + cpatch%broot    (ico)          &
                                  + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico)          &
                                  + cpatch%bbark    (ico)
            
            !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
            cpatch%today_nppleaf   (ico)= tr_bleaf                      * cpatch%nplant(ico)
            cpatch%today_nppfroot  (ico)= tr_broot                      * cpatch%nplant(ico)
            cpatch%today_nppsapwood(ico)= (tr_bsapwooda + tr_bsapwoodb) * cpatch%nplant(ico)
            cpatch%today_nppbark   (ico)= tr_bbark                      * cpatch%nplant(ico)
            cpatch%today_nppdaily  (ico)= carbon_balance * cpatch%nplant(ico)
            
            !------------------------------------------------------------------------------!
            !    Find the amount of carbon used to recover the tissues that were off-      !
            ! -allometry, take that from the carbon balance first, then use some of the    !
            ! storage if needed be.                                                        !
            !------------------------------------------------------------------------------!
            increment = carbon_balance -  tr_bleaf - tr_broot - tr_bsapwooda               &
                      - tr_bsapwoodb - tr_bbark
            cpatch%bstorage(ico) = max(0.0, cpatch%bstorage(ico) + increment)
            !------------------------------------------------------------------------------!

            if (increment <= 0.0)  then
               !---------------------------------------------------------------------------!
               !    We are using up all of daily C gain and some of bstorage.  First       !
               ! calculate N demand from using daily C gain.                               !
               !---------------------------------------------------------------------------!
               if (carbon_balance < 0.0) then
                  nitrogen_uptake = nitrogen_uptake + carbon_balance / c2n_storage
                  nitrogen_uptake = nitrogen_uptake                                        &
                                  + (carbon_balance - increment)                           &
                                  * ( f_labile(ipft) / c2n_leaf(ipft)                      &
                                    + (1.0 - f_labile(ipft)) / c2n_stem(ipft)              &
                                    -  1.0 / c2n_storage)
                  
               else
                  nitrogen_uptake = nitrogen_uptake + carbon_balance                       &
                                 * ( f_labile(ipft) / c2n_leaf(ipft)                       &
                                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )

                  !------------------------------------------------------------------------!
                  !     Now calculate additional N uptake required from transfer of C from !
                  ! storage to balive.                                                     !
                  !------------------------------------------------------------------------!
                  nitrogen_uptake  = nitrogen_uptake +  ( - 1.0 * increment )              &
                                   * ( f_labile(ipft)  / c2n_leaf(ipft)                    &
                                     + (1.0 - f_labile(ipft)) / c2n_stem(ipft)             &
                                     -  1.0 / c2n_storage)
               end if

            else
               !---------------------------------------------------------------------------!
               !     N uptake for fraction of daily C gain going to balive.                !
               !---------------------------------------------------------------------------!
               nitrogen_uptake = nitrogen_uptake + (carbon_balance - increment)            &
                               * ( f_labile(ipft) / c2n_leaf(ipft)                         &
                                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
               !----- N uptake for fraction of daily C gain going to bstorage. ------------!
               nitrogen_uptake = nitrogen_uptake + increment / c2n_storage
            end if

            on_allometry = 2.0 * abs(balive_aim - cpatch%balive(ico))                      &
                         / (balive_aim + cpatch%balive(ico))          < 1.e-6
            if (cpatch%elongf(ico) == 1.0 .and. on_allometry) then
               !---------------------------------------------------------------------------!
               !     We're back to allometry, change phenology_status.                     !
               !---------------------------------------------------------------------------!
               cpatch%phenology_status(ico) = 0
            end if
         case default
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
            nitrogen_uptake      = nitrogen_uptake      + carbon_balance / c2n_storage
            !------------------------------------------------------------------------------!


            !----- NPP allocation in diff pools in Kg C/m2/day. ---------------------------!
            cpatch%today_nppleaf(ico)    = 0.0
            cpatch%today_nppfroot(ico)   = 0.0
            cpatch%today_nppsapwood(ico) = 0.0
            cpatch%today_nppbark(ico)    = 0.0
            cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


      else
         !---------------------------------------------------------------------------------!
         !   Carbon balance is negative, take it out of storage.                           !
         !---------------------------------------------------------------------------------!
         carbon_loss = - (cpatch%bstorage(ico) + carbon_balance)

         if (carbon_loss >= 0.0)  then
            !----- Use Storage pool first then take out of balive. ------------------------!
            cpatch%bstorage(ico) = 0.0
            csite%fsn_in(ipa)    = csite%fsn_in(ipa)
         else
            !------ Burn the storage pool.  Dont' forget the nitrogen. --------------------!
            cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
            csite%fsn_in(ipa)    = csite%fsn_in(ipa)
         end if

         !---- NPP allocation in diff pools in KgC/m2/day. --------------------------------!
         cpatch%today_nppleaf(ico)    = 0.0
         cpatch%today_nppfroot(ico)   = 0.0
         cpatch%today_nppsapwood(ico) = 0.0
         cpatch%today_nppbark(ico)    = 0.0
         cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine alloc_plant_c_balance_eq_0
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   subroutine potential_N_uptake(cpatch,ico,balive_in,carbon_balance_pot,N_uptake_pot      &
                                ,green_leaf_factor)
      use ed_state_vars , only : patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : size2bl      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: balive_in
      real           , intent(in)    :: carbon_balance_pot
      real           , intent(inout) :: N_uptake_pot
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      !------------------------------------------------------------------------------------!

      ipft = cpatch%pft(ico) 

      if ( cpatch%phenology_status(ico) == 0 .and. carbon_balance_pot > 0.0 ) then

         !----- Positive carbon balance with plants fully flushed. ------------------------!
         N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage

      elseif (cpatch%phenology_status(ico) == 1) then
         ! this calculation of bl_max is wrong for grass, but they should not have phenology_status=1 yet
         bl_max = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)                           &
                * green_leaf_factor * cpatch%elongf(ico)
         bl_pot = cpatch%bleaf(ico) + carbon_balance_pot

         if (bl_pot > bl_max) then
            !------------------------------------------------------------------------------!
            !     This increment would take us over the limit, so we assign all that can   !
            ! go for leaves to them, and put the remainder in storage.                     !
            !------------------------------------------------------------------------------!
            increment    = carbon_balance_pot - (bl_max-cpatch%bleaf(ico))
            N_uptake_pot = N_uptake_pot + increment / c2n_storage
            increment    = bl_max-cpatch%bleaf(ico)
            N_uptake_pot = N_uptake_pot + increment                                        &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
         elseif (carbon_balance_pot > 0.0) then

            !------------------------------------------------------------------------------!
            !      This increment did not exceed the limit, put everything in leaves.  We  !
            ! don't compute the uptake if carbon balance is negative, just because there   !
            ! will be no uptake...                                                         !
            !------------------------------------------------------------------------------!
            N_uptake_pot = N_uptake_pot + carbon_balance_pot                               &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
         end if
      end if

      return
   end subroutine potential_N_uptake
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine litter(csite,ipa)

      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! ! structure
      use pft_coms     , only : c2n_leaf  & ! intent(in)
                              , c2n_stem  & ! intent(in)
                              , l2n_stem  ! ! intent(in)
      use decomp_coms  , only : f_labile  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !      Add fine root, leaf, and bark turnover to the litter.                         !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         plant_litter   = cpatch%nplant(ico)                                               &
                        * ( cpatch%leaf_maintenance(ico) + cpatch%root_maintenance(ico) )
         plant_litter_f = plant_litter * f_labile(ipft)
         plant_litter_s = plant_litter * (1.0 - f_labile(ipft))                            &
                        + cpatch%nplant(ico) * cpatch%bark_maintenance(ico)

         csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f
         csite%fsn_in(ipa) = csite%fsn_in(ipa) + plant_litter_f / c2n_leaf(ipft)

         csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
         csite%ssl_in(ipa) = csite%ssl_in(ipa) + plant_litter_s * l2n_stem / c2n_stem(ipft)
      end do
      return
   end subroutine litter
   !=======================================================================================!
   !=======================================================================================!
end module growth_balive
!==========================================================================================!
!==========================================================================================!
