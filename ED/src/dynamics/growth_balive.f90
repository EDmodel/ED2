!============================================================================!
!============================================================================!
module growth_balive
   !=========================================================================!
   !=========================================================================!


   contains



   !=========================================================================!
   !=========================================================================!
   !     This subroutine will update the alive biomass, and compute the      !
   ! respiration terms other than leaf respiration.                          !
   ! IMPORTANT: The order of the operations here affect the C/N budgests, so !
   !            don't change the order of the operations unless you really   !
   !            know what you are doing.                                     !
   !-------------------------------------------------------------------------!
   subroutine dbalive_dt(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_hcapveg           & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           ! ! subroutine
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -----------------------------------------------!
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
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_hcapveg
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      real                          :: temp_dep
      !----- External functions. --------------------------------------------!
      real             , external   :: mortality_rates


      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -----------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0

               !----- Loop over cohorts. ------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -----------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Initialize cohort nitrogen uptake. -----------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0
                  
                  ! Set allocation factors
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc
                  
                  !----- Subtract preceding day's storage respiration. ------!
                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                &
                                       - cpatch%storage_respiration(ico)

                  !----------------------------------------------------------!
                  !     When storage carbon is lost, allow the associated    !
                  ! nitrogen to go to litter in order to maintain prescribed !
                  ! C2N ratio.                                               !
                  !----------------------------------------------------------!
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                      &
                                    + cpatch%storage_respiration(ico)        &
                                    / c2n_storage * cpatch%nplant(ico)

                  !----------------------------------------------------------!
                  !     If plants are off allometry, move carbon from        !
                  ! bstorage to balive.                                      !
                  !----------------------------------------------------------!
                  call transfer_C_from_storage(cpatch,ico,salloc             &
                                              ,nitrogen_uptake,N_uptake_pot)
                  
                  !----- Calculate leaf, fine root biomass. -----------------!
                  if (cpatch%phenology_status(ico) /= 2) then
                     bl = cpoly%green_leaf_factor(ipft,isi)                  &
                        * cpatch%balive(ico) * salloci
                  else
                     bl = 0.0
                  end if

                  br = q(ipft) * cpatch%balive(ico) * salloci 

                  !----------------------------------------------------------!
                  !     Compute maintenance costs.                           !
                  !----------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain &
                                        ,csite%avg_daily_temp(ipa))
                  
                  !----- Subtract maintenance costs from balive. ------------!
                  cpatch%balive(ico)    = cpatch%balive(ico)                 &
                                        - cpatch%leaf_maintenance_costs(ico) &
                                        - cpatch%root_maintenance_costs(ico)
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                  &
                                        - cpatch%leaf_maintenance_costs(ico) &
                                        - cpatch%root_maintenance_costs(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)              &
                                        - cpatch%leaf_maintenance_costs(ico) &
                                        - cpatch%root_maintenance_costs(ico)

                  !----------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon      !
                  ! balances.                                                !
                  !----------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain     &
                                            ,carbon_balance                  &
                                            ,carbon_balance_pot              &
                                            ,carbon_balance_max)

                  !----------------------------------------------------------!
                  !      Compute respiration rates for coming day            !
                  ! [kgC/plant/day].                                         !
                  !----------------------------------------------------------!
                  cpatch%growth_respiration(ico) =                           &
                          max(0.0, daily_C_gain * growth_resp_factor(ipft))
                  temp_dep = 1.0!! / (1.0 + exp(0.4 * (278.15 - csite%avg_daily_temp(ipa))))  !! experimental and arbitrary (borrowed from maintainence temp dep) [[MCD]]
                  cpatch%storage_respiration(ico) =                          &
                          cpatch%bstorage(ico) * storage_turnover_rate(ipft) &
                        * tfact * temp_dep
                  cpatch%vleaf_respiration(ico) =                            &
                          (1.0 - cpoly%green_leaf_factor(ipft,isi))          &
                        / (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))     &
                        * cpatch%balive(ico) * storage_turnover_rate(ipft)   &
                        * tfact * temp_dep


                  !----------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and         !
                  ! bstorage.                                                !
                  !----------------------------------------------------------!
                  balive_in = cpatch%balive(ico)
                  call alloc_plant_c_balance(csite,ipa,ico,salloc,salloci    &
                                  ,carbon_balance,nitrogen_uptake            &
                                  ,cpoly%green_leaf_factor(ipft,isi))

                  !----------------------------------------------------------!
                  !     Do a shadow calculation to see what would have       !
                  ! happened if stomata were open.  This is used to          !
                  ! calculate potential nitrogen uptake, N_uptake_pot.       !
                  !----------------------------------------------------------!

                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,salloc,salloci       &
                                  ,balive_in,carbon_balance_pot,N_uptake_pot &
                                  ,cpoly%green_leaf_factor(ipft,isi))
                  end if

                  !----------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.    !
                  !----------------------------------------------------------!
                  csite%total_plant_nitrogen_uptake(ipa) =                   &
                         csite%total_plant_nitrogen_uptake(ipa)              &
                       + nitrogen_uptake * cpatch%nplant(ico)

                  !----- Calculate plant N limitation factor. ---------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     br              = q(ipft) * cpatch%balive(ico) * salloci
                     nitrogen_supply = plant_N_supply_scale * br             &
                                     * csite%mineralized_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply                       &
                                     / (nitrogen_supply + N_uptake_pot)
                  end if
                  
                  !----------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality     !
                  ! changes daily.                                           !
                  !----------------------------------------------------------!
                  dndt = - mortality_rates(cpatch,ipa,ico                    &
                                          ,csite%avg_daily_temp(ipa))        &
                       * cpatch%nplant(ico) * tfact
                  
                  !----- Update monthly mortality rate [plants/m2/month]. ---!
                  cpatch%monthly_dndt(ico) = cpatch%monthly_dndt(ico) + dndt

     
                  !----- Updating LAI, WPA, and WAI. ------------------------!
                  call area_indices( cpatch%nplant(ico), cpatch%bleaf(ico)   &
                                   , cpatch%bdead(ico) , cpatch%balive(ico)  &
                                   , cpatch%dbh(ico)   , cpatch%hite(ico)    &
                                   , cpatch%pft(ico)   , cpatch%sla(ico)     &
                                   , cpatch%lai(ico)   , cpatch%wpa(ico)     &
                                   , cpatch%wai(ico)   )

                  !----------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore,    !
                  ! update vegetation energy and heat capacity.              !
                  !----------------------------------------------------------!
                  old_hcapveg = cpatch%hcapveg(ico)
                  cpatch%hcapveg(ico) =                                      &
                         calc_hcapveg(cpatch%bleaf(ico) ,cpatch%bdead(ico)   &
                                     ,cpatch%balive(ico),cpatch%nplant(ico)  &
                                     ,cpatch%hite(ico)  ,cpatch%pft(ico)     &
                                     ,cpatch%phenology_status(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_hcapveg)
                  !----- Likewise, the total heat capacity must be updated. -!
                  csite%hcapveg(ipa) = csite%hcapveg(ipa)                    &
                                     + cpatch%hcapveg(ico) - old_hcapveg
                  !----------------------------------------------------------!
               end do
               
               !----- Update litter. ----------------------------------------!
               call litter(csite,ipa)
               
               !----- Update patch LAI, WAI, height, roughness... -----------!
               call update_patch_derived_props(csite,cpoly%lsl(isi)          &
                                              ,cpoly%met(isi)%prss,ipa)

               !----- Recalculate storage terms (for budget assessment). ----!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

               !----- It's a new day, reset average daily temperature. ------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do

      return
   end subroutine dbalive_dt
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   !    This subroutine will transfer some of the stored carbon to balive in !
   ! order to put the plant back on allometry.                               !
   !-------------------------------------------------------------------------!
   subroutine transfer_C_from_storage(cpatch,ico,salloc,nitrogen_uptake      &
                                     ,N_uptake_pot)
      use ed_state_vars , only : patchtype
      use pft_coms      , only : c2n_leaf    & ! intent(in)
                               , c2n_storage & ! intent(in)
                               , c2n_stem    ! ! intent(in)
      use decomp_coms   , only : f_labile    ! ! intent(in)
      use allometry     , only : dbh2bl      ! ! function
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: N_uptake_pot
      !----- Local variables. -----------------------------------------------!
      integer                        :: ipft
      real                           :: off_allometry_cb
      real                           :: increment
      !----------------------------------------------------------------------!


      !----------------------------------------------------------------------!
      !     Only do the transfer if leaves exist.                            !
      !----------------------------------------------------------------------!
      if (cpatch%phenology_status(ico) == 2) return
     
      !----- Alias for pft type. --------------------------------------------!
      ipft = cpatch%pft(ico)
     
      !----- Determine how much biomass we need to go back to allometry. ----!
      off_allometry_cb = dbh2bl(cpatch%dbh(ico),ipft) * salloc               &
                       - cpatch%balive(ico)

      !----- If plants have storage, transfer it to balive. -----------------!
      increment            = max(0.0,min(max(0.0, off_allometry_cb)          &
                                ,cpatch%bstorage(ico)))
      cpatch%balive(ico)   = cpatch%balive(ico) + increment
      cpatch%bstorage(ico) = cpatch%bstorage(ico) - increment

      !----------------------------------------------------------------------!
      !      N uptake is required since c2n_leaf < c2n_storage.  Units are   !
      ! kgN/plant/day.                                                       !
      !----------------------------------------------------------------------!
      nitrogen_uptake = increment                                            &
                      * (        f_labile(ipft)  / c2n_leaf(ipft)            &
                        + (1.0 - f_labile(ipft)) / c2n_stem(ipft)            &
                        -  1.0 / c2n_storage)
      N_uptake_pot    = nitrogen_uptake

      return
   end subroutine transfer_C_from_storage
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   subroutine plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain,tempk)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : phenology          & ! intent(in)
                              , root_turnover_rate & ! intent(in)
                              , leaf_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: br
      real           , intent(in)    :: bl
      real           , intent(in)    :: tfact
      real           , intent(in)    :: tempk
      real           , intent(out)   :: daily_C_gain
      !----- Local variables. -----------------------------------------------!
      integer                        :: ipft
      real                           :: maintenance_temp_dep
      !----------------------------------------------------------------------!

      !------ Alias for plant functional type. ------------------------------!
      ipft = cpatch%pft(ico)

      !------ Get the temperature dependence. -------------------------------!
      if (phenology(ipft) == 0) then
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
      else
         maintenance_temp_dep = 1.0
      end if

      !----- Calculate maintenance demand (kgC/plant/year). -----------------!
      cpatch%root_maintenance_costs(ico) = root_turnover_rate(ipft) * br * maintenance_temp_dep
      if(phenology(ipft) /= 3)then
         cpatch%leaf_maintenance_costs(ico) = leaf_turnover_rate(ipft) * bl * maintenance_temp_dep
      else
         cpatch%leaf_maintenance_costs(ico) = leaf_turnover_rate(ipft) * bl * cpatch%turnover_amp(ico) * maintenance_temp_dep
      endif


      !----- Convert units of maintenance to [kgC/plant/day]. ---------------!
      cpatch%root_maintenance_costs(ico) = cpatch%root_maintenance_costs(ico) * tfact
      cpatch%leaf_maintenance_costs(ico) = cpatch%leaf_maintenance_costs(ico) * tfact

      !----- Compute daily C uptake [kgC/plant/day]. ------------------------!
      if(cpatch%nplant(ico) > tiny(1.0)) then
         daily_C_gain = umol_2_kgC * day_sec                                 &
                      * ( cpatch%dmean_gpp(ico)                              &
                        - cpatch%dmean_leaf_resp(ico)                        &
                        - cpatch%dmean_root_resp(ico))                       &
                      / cpatch%nplant(ico)
      else
         daily_C_gain = 0.0
      end if

      return
   end subroutine plant_maintenance
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   subroutine plant_carbon_balances(cpatch,ipa,ico,daily_C_gain              &
                                   ,carbon_balance,carbon_balance_pot        &
                                   ,carbon_balance_max)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : growth_resp_factor ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      use ed_max_dims  , only : n_pft              ! ! intent(in)
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: daily_C_gain
      real                     , intent(out) :: carbon_balance
      real                     , intent(out) :: carbon_balance_pot
      real                     , intent(out) :: carbon_balance_max
      !----- Local variables. -----------------------------------------------!
      real                                   :: daily_C_gain_pot
      real                                   :: daily_C_gain_max
      real                                   :: growth_respiration_pot
      real                                   :: growth_respiration_max
      integer                                :: ipft
      !----- Local constants. -----------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. ---------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !----------------------------------------------------------------------!

      !----- Alias for PFT type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Calculate actual daily carbon balance: kgC/plant/day. ---------!
      carbon_balance = daily_C_gain - cpatch%growth_respiration(ico)         &
                                    - cpatch%vleaf_respiration(ico)

      if (cpatch%nplant(ico) > tiny(1.0)) then

         !-------------------------------------------------------------------!
         !      Calculate potential carbon balance (used for nitrogen demand !
         ! function).  [kgC/plant/day].                                      !
         !-------------------------------------------------------------------!
         daily_C_gain_pot       = umol_2_kgC * day_sec                       &
                                * ( cpatch%dmean_gpp_pot(ico)                &
                                  - cpatch%dmean_leaf_resp(ico)              &
                                  - cpatch%dmean_root_resp(ico))             &
                                / cpatch%nplant(ico)
         growth_respiration_pot = max(0.0, daily_C_gain_pot                  &
                                         * growth_resp_factor(ipft) )
         carbon_balance_pot     = daily_C_gain_pot - growth_respiration_pot  &
                                - cpatch%vleaf_respiration(ico)

         !----- Calculate maximum carbon balance (used for mortality). ------!
         daily_C_gain_max       = umol_2_kgC * day_sec                       &
                                * ( cpatch%dmean_gpp_max(ico)                &
                                  - cpatch%dmean_leaf_resp(ico)              &
                                  - cpatch%dmean_root_resp(ico) )            &
                                / cpatch%nplant(ico)
         growth_respiration_max = max(0.0, daily_C_gain_max                  &
                                         * growth_resp_factor(ipft))
         carbon_balance_max     = daily_C_gain_max                           &
                                - growth_respiration_max                     &
                                - cpatch%vleaf_respiration(ico)
      else
         carbon_balance_max = 0.0
         carbon_balance_pot = 0.0
      end if

      !----- Carbon balances for mortality. ---------------------------------!
      cpatch%cb(13,ico)     = cpatch%cb(13,ico) + carbon_balance
      cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) + carbon_balance_max

      if (print_debug) then

         if (first_time(ipft)) then
            first_time(ipft) = .false.
            write (unit=30+ipft,fmt='(a10,15(1x,a12))')                      &
                '      TIME','       PATCH','      COHORT','      NPLANT'    &
                            ,'    CB_TODAY',' GROWTH_RESP','  VLEAF_RESP'    &
                            ,'   DMEAN_GPP','DMEAN_GPPMAX','  DMEAN_LEAF'    &
                            ,'  DMEAN_ROOT',' CBMAX_TODAY','          CB'    &
                            ,'       CBMAX','  LEAF_MAINT','  ROOT_MAINT'
         end if

         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),13(1x,es12.5))') &
              current_time%month,'/',current_time%date,'/',current_time%year &
             ,ipa,ico,cpatch%nplant(ico),carbon_balance                      &
             ,cpatch%growth_respiration(ico),cpatch%vleaf_respiration(ico)   &
             ,cpatch%dmean_gpp(ico),cpatch%dmean_gpp_max(ico)                &
             ,cpatch%dmean_leaf_resp(ico),cpatch%dmean_root_resp(ico)        &
             ,carbon_balance_max,cpatch%cb(13,ico),cpatch%cb_max(13,ico)     &
             ,cpatch%leaf_maintenance_costs(ico) ,cpatch%root_maintenance_costs(ico)
      end if

      return
   end subroutine plant_carbon_balances
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   subroutine alloc_plant_c_balance(csite,ipa,ico,salloc,salloci             &
                                   ,carbon_balance,nitrogen_uptake           &
                                   ,green_leaf_factor)
      use ed_state_vars, only : sitetype     & ! structure
                              , patchtype    ! ! structure
      use pft_coms     , only : c2n_storage  & ! intent(in)
                              , c2n_leaf     & ! intent(in)
                              , sla          & ! intent(in)
                              , c2n_stem     ! ! intent(in)
      use decomp_coms  , only : f_labile     ! ! intent(in)
      use allometry    , only : dbh2bl       ! ! function
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -----------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: old_status
      !----------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      if (cpatch%phenology_status(ico) == 0 .and. carbon_balance > 0.0) then
         !-------------------------------------------------------------------!
         !      Simply update monthly carbon gain.  This will be used for    !
         ! structural growth at the end of the month.                        !
         !-------------------------------------------------------------------!
         cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
         nitrogen_uptake      = nitrogen_uptake + carbon_balance             &
                                                / c2n_storage
         cpatch%bleaf(ico)    = cpatch%balive(ico) * salloci                 &
                              * green_leaf_factor

      elseif (cpatch%phenology_status(ico) < 2) then
         !-------------------------------------------------------------------!
         !      There are leaves.  Here we will compute the maximum amount   !
         ! that can go to leaves, and put any excess in storage.             !
         !-------------------------------------------------------------------!

         !----- Maximum bleaf that the allometric relationship would allow. -!
         bl_max = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor

         !----- Maximum bleaf if all gain went to leaf biomass. -------------!
         bl_pot = green_leaf_factor                                          &
                * (cpatch%balive(ico) + carbon_balance) * salloci

         if (bl_pot >= bl_max) then
            !----------------------------------------------------------------!
            !     The carbon gain is more than what can go to leaves.  Add   !
            ! all that can to go leaves to leaf biomass, and put the remain- !
            ! der in the storage.                                            !
            !----------------------------------------------------------------!
            increment = carbon_balance                                       &
                      - ( dbh2bl(cpatch%dbh(ico),ipft) * salloc              &
                        - cpatch%balive(ico))
            cpatch%bstorage(ico) = cpatch%bstorage(ico) + increment

            !----------------------------------------------------------------!
            !    Compute nitrogen uptake in order to conserve C/N ratio.     !
            !----------------------------------------------------------------!
            nitrogen_uptake    = nitrogen_uptake + increment / c2n_storage
            increment          = dbh2bl(cpatch%dbh(ico),ipft) * salloc       &
                               - cpatch%balive(ico)
            cpatch%balive(ico) = cpatch%balive(ico) + increment
            nitrogen_uptake    = nitrogen_uptake                             &
                               + increment * ( f_labile(ipft)                &
                                             / c2n_leaf(ipft)                &
                                             + (1.0 - f_labile(ipft))        &
                                             / c2n_stem(ipft))
            cpatch%bleaf(ico)  = bl_max
            !----------------------------------------------------------------!
            !    That is the maximum leaf biomass that this cohort can have. !
            ! Therefore, we set phenology status to 0 (leaves not growing).  !
            !----------------------------------------------------------------!
            cpatch%phenology_status(ico) = 0

         else

            !----------------------------------------------------------------!
            !    It will not exceed limit, so just add to balive.  There is  !
            ! still some room for leaf growing, so phenology status is set   !
            ! to 1.                                                          !
            !----------------------------------------------------------------!
            cpatch%balive(ico) = max(0.0,cpatch%balive(ico) + carbon_balance)
            cpatch%phenology_status(ico) = 1
            cpatch%bleaf(ico) = cpatch%balive(ico) * salloci                 &
                              * green_leaf_factor

            !----------------------------------------------------------------!
            !     Update nitrogen uptake/soil nitrogen to preserve the C/N   !
            ! ratio.                                                         !
            !----------------------------------------------------------------!
            if (carbon_balance < 0.0) then
               csite%fsn_in(ipa) = csite%fsn_in(ipa)                         &
                                 - carbon_balance                            &
                                 * ( f_labile(ipft) / c2n_leaf(ipft)         &
                                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft))     &
                                 * cpatch%nplant(ico)
            else
               nitrogen_uptake = nitrogen_uptake                             &
                               + carbon_balance                              &
                               * ( f_labile(ipft) / c2n_leaf(ipft)           &
                                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
            end if
         end if

      else
         !-------------------------------------------------------------------!
         !      In this case, either carbon balance is negative.  Simply add !
         ! the carbon balance.                                               !
         !-------------------------------------------------------------------!
         cpatch%balive(ico) = max(0.0,cpatch%balive(ico) + carbon_balance)
         cpatch%bleaf(ico)  = 0.0
         csite%fsn_in(ipa)  = csite%fsn_in(ipa)                              &
                            - carbon_balance                                 &
                            * ( f_labile(ipft) / c2n_leaf(ipft)              &
                              + (1.0 - f_labile(ipft)) / c2n_stem(ipft))&
                            * cpatch%nplant(ico)

      end if

      return
   end subroutine alloc_plant_c_balance
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   subroutine potential_N_uptake(cpatch,ico,salloc,salloci,balive_in         &
                                ,carbon_balance_pot,N_uptake_pot             &
                                ,green_leaf_factor)
      use ed_state_vars , only : patchtype   ! ! structure
      use pft_coms      , only : c2n_storage & ! intent(in)
                               , c2n_leaf    & ! intent(in)
                               , c2n_stem    ! ! intent(in)
      use decomp_coms   , only : f_labile    ! ! intent(in)
      use allometry     , only : dbh2bl      ! ! intent(in)
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: balive_in
      real           , intent(in)    :: carbon_balance_pot
      real           , intent(inout) :: N_uptake_pot
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -----------------------------------------------!
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      !----------------------------------------------------------------------!

      ipft = cpatch%pft(ico) 

      if (cpatch%phenology_status(ico) == 0 .and.                            &
          carbon_balance_pot > 0.0 ) then

         !----- Positive carbon balance with plants fully flushed. ----------!
         N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage

      elseif (cpatch%phenology_status(ico) < 2) then

         !-------------------------------------------------------------------! 
         !    There are at least some leaves.  First we compute the maximum  !
         ! possible bleaf and the how much leaf biomass we could attain if   !
         ! all carbon balance went to leaf biomass.  We then decide what to  !
         ! do based on whether we would exceed the limit or not.             !
         !-------------------------------------------------------------------! 
         bl_max = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor
         bl_pot = green_leaf_factor * (balive_in + carbon_balance_pot)       &
                * salloci

         if (bl_pot > bl_max) then
            !----------------------------------------------------------------!
            !     This increment would take us over the limit, so we assign  !
            ! all that can go for leaves to them, and put the remainder in   !
            ! in storage.                                                    !
            !----------------------------------------------------------------!
            increment    = carbon_balance_pot                                &
                         - (dbh2bl(cpatch%dbh(ico),ipft)*salloc - balive_in)
            N_uptake_pot = N_uptake_pot + increment / c2n_storage
            increment    = dbh2bl(cpatch%dbh(ico),ipft)*salloc - balive_in
            N_uptake_pot = N_uptake_pot                                      &
                         + increment * ( f_labile(ipft) / c2n_leaf(ipft)     &
                                       + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
         elseif (carbon_balance_pot > 0.0) then

            !----------------------------------------------------------------!
            !      This increment did not exceed the limit, put everything   !
            ! in leaves.  We don't compute the uptake if carbon balance is   !
            ! negative, just because there will be no uptake...              !
            !----------------------------------------------------------------!
            N_uptake_pot = N_uptake_pot                                      &
                         + carbon_balance_pot                                &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                 &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))

         end if
      end if

      return
   end subroutine potential_N_uptake
   !=========================================================================!
   !=========================================================================!






   !=========================================================================!
   !=========================================================================!
   subroutine litter(csite,ipa)

      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! ! structure
      use pft_coms     , only : c2n_leaf  & ! intent(in)
                              , c2n_stem  & ! intent(in)
                              , l2n_stem  ! ! intent(in)
      use decomp_coms  , only : f_labile  ! ! intent(in)
      implicit none
      !----- Arguments. -----------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -----------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
      !----------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !----------------------------------------------------------------------!
      !      Add fine root and leaf turnover to the litter.                  !
      !----------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         plant_litter   = (cpatch%leaf_maintenance_costs(ico) &
          + cpatch%root_maintenance_costs(ico)) * cpatch%nplant(ico)
         plant_litter_f = plant_litter * f_labile(ipft)
         plant_litter_s = plant_litter - plant_litter_f

         csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f
         csite%fsn_in(ipa) = csite%fsn_in(ipa)                               &
                           + plant_litter_f / c2n_leaf(ipft)

         csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
         csite%ssl_in(ipa) = csite%ssl_in(ipa)                               &
                           + plant_litter_s * l2n_stem / c2n_stem(ipft)
      end do
      return
   end subroutine litter
   !=========================================================================!
   !=========================================================================!
end module growth_balive
!============================================================================!
!============================================================================!
