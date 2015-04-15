!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the photosynthesis scheme (Farquar and Leuning).  This  !
! is called every step, but not every sub-step.                                            !
!------------------------------------------------------------------------------------------!
subroutine canopy_photosynthesis(csite,cmet,mzg,ipa,lsl,ntext_soil                         &
                                ,leaf_aging_factor,green_leaf_factor)
   use ed_state_vars  , only : sitetype           & ! structure
                             , patchtype          ! ! structure
   use ed_max_dims    , only : n_pft              ! ! intent(in)
   use pft_coms       , only : leaf_width         & ! intent(in)
                             , water_conductance  & ! intent(in)
                             , include_pft        & ! intent(in)
                             , vm0                & ! intent(in)
                             , leaf_turnover_rate ! ! intent(in)
   use soil_coms      , only : soil               & ! intent(in)
                             , slz                & ! intent(in)
                             , slzt               & ! intent(in)
                             , dslz               ! ! intent(in)
   use consts_coms    , only : t00                & ! intent(in)
                             , epi                & ! intent(in)
                             , wdnsi              & ! intent(in)
                             , wdns               & ! intent(in)
                             , umols_2_kgCyr      & ! intent(in)
                             , yr_day             & ! intent(in)
                             , lnexp_min          & ! intent(in)
                             , tiny_num           ! ! intent(in)
   use ed_misc_coms   , only : current_time       & ! intent(in)
                             , dtlsm              & ! intent(in)
                             , frqsum             ! ! intent(in)
   use met_driver_coms, only : met_driv_state     ! ! structure
   use physiology_coms, only : print_photo_debug  & ! intent(in)
                             , h2o_plant_lim      ! ! intent(in)
   use phenology_coms , only : llspan_inf         ! ! intent(in)
   use farq_leuning   , only : lphysiol_full      ! ! sub-routine
   use allometry      , only : h2crownbh          ! ! function
   use therm_lib      , only : qslif              ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite             ! Current site
   type(met_driv_state)      , target      :: cmet              ! Current met. conditions.
   integer                   , intent(in)  :: ipa               ! Current patch #
   integer                   , intent(in)  :: lsl               ! Lowest soil level
   integer                   , intent(in)  :: mzg               ! Number of soil layers
   integer, dimension(mzg)   , intent(in)  :: ntext_soil        ! Soil class
   real   , dimension(n_pft) , intent(in)  :: leaf_aging_factor ! 
   real   , dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch             ! Current site
   integer                                 :: ico                ! Current cohort #
   integer                                 :: tuco               ! Tallest used cohort
   integer                                 :: ipft
   integer                                 :: tpft
   integer                                 :: k
   integer                                 :: kroot
   integer                                 :: nsoil
   integer                                 :: limit_flag
   logical, dimension(mzg+1)               :: root_depth_indices ! 
   logical                                 :: las
   real   , dimension(:)    , allocatable  :: avail_h2o_coh
   real                                    :: leaf_par
   real                                    :: leaf_resp
   real                                    :: d_A_light_max
   real                                    :: d_A_rubp_max
   real                                    :: d_A_co2_max
   real                                    :: d_gsw_open
   real                                    :: d_gsw_closed
   real                                    :: d_lsfc_shv_open
   real                                    :: d_lsfc_shv_closed
   real                                    :: d_lsfc_co2_open
   real                                    :: d_lsfc_co2_closed
   real                                    :: d_lint_co2_open
   real                                    :: d_lint_co2_closed
   real                                    :: swp
   real                                    :: vm
   real                                    :: mcheight
   real                                    :: compp
   real                                    :: broot_tot
   real                                    :: broot_loc
   real                                    :: water_demand
   real                                    :: psiplusz
   real                                    :: avail_h2o_lyr
   real                                    :: wilting_factor
   real                                    :: pss_available_water
   real                                    :: vm0_tuco
   real                                    :: llspan_tuco
   real                                    :: can_ssh
   integer, dimension(n_pft)               :: tuco_pft
   !----- Locally saved variables. --------------------------------------------------------!
   real                          , save    :: dtlsm_o_frqsum
   logical                       , save    :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Assign the constant scaling factor. ---------------------------------------------!
   if (first_time) then
      first_time     = .false.
      dtlsm_o_frqsum = dtlsm / frqsum
   end if
   !---------------------------------------------------------------------------------------!


   !----- Point to the cohort structures --------------------------------------------------!
   cpatch => csite%patch(ipa)
   !---------------------------------------------------------------------------------------!



   !----- Allocate the available water function for plants. -------------------------------!
   if (cpatch%ncohorts > 0) then
      allocate (avail_h2o_coh(cpatch%ncohorts))
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Calculate liquid water available for transpiration.   The way this is done        !
   ! depends on how the water limitation is to be solved.                                  !
   !---------------------------------------------------------------------------------------!
   select case (h2o_plant_lim)
   case (0,1)
      !------------------------------------------------------------------------------------!
      !     Available water is defined as the soil moisture (mass) above wilting point,    !
      ! scaled by liquid water fraction.                                                   !
      !------------------------------------------------------------------------------------!
      do ico = 1, cpatch%ncohorts
         !---- Aliases for rooting depth and PFT. -----------------------------------------!
         kroot = cpatch%krdepth(ico)
         ipft  = cpatch%pft(ico)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the available water for each layer.                                    !
         !---------------------------------------------------------------------------------!
         avail_h2o_coh(ico) = 0.
         do k=mzg,kroot,-1
            !----- Alias for soil type. ---------------------------------------------------!
            nsoil = ntext_soil(k)
            !------------------------------------------------------------------------------!



            !----- Find the available water factor for this layer. ------------------------!
            avail_h2o_lyr = max(0.0, (csite%soil_water(k,ipa) - soil(nsoil)%soilwp))       &
                          * csite%soil_fracliq(k,ipa) * wdns * dslz(k)
            !------------------------------------------------------------------------------!



            !----- Add the factor from this layer to the integral. ------------------------!
            avail_h2o_coh(ico) = avail_h2o_coh(ico) + avail_h2o_lyr
            !------------------------------------------------------------------------------!
         end do
      end do
      !------------------------------------------------------------------------------------!

   case (2)
      !------------------------------------------------------------------------------------!
      !     The available water factor is the soil moisture at field capacity minus wilt-  !
      ! ing point, scaled by the wilting factor, defined as a function of soil potential   !
      ! and height between roots and mid-crown.                                            !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         !---- Aliases for rooting depth and PFT. -----------------------------------------!
         kroot = cpatch%krdepth(ico)
         ipft  = cpatch%pft(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the mean height of the crown (to represent the distance between        !
         ! the ground and the leaves.                                                      !
         !---------------------------------------------------------------------------------!
         mcheight = 0.5 * (cpatch%hite(ico) + h2crownbh(cpatch%hite(ico),ipft))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the available water for each layer.                                    !
         !---------------------------------------------------------------------------------!
         avail_h2o_coh(ico) = 0.
         do k = mzg, kroot, -1
            !----- Alias for soil type. ---------------------------------------------------!
            nsoil = ntext_soil(k)
            !------------------------------------------------------------------------------!


            !----- Find the potential for this layer. -------------------------------------!
            psiplusz = slzt(k) - mcheight + csite%soil_mstpot(k,ipa)
            !------------------------------------------------------------------------------!


            !----- Find the available water factor for this layer. ------------------------!
            wilting_factor   = (psiplusz - soil(nsoil)%slpotwp)                            &
                             / (soil(nsoil)%slpotfc - soil(nsoil)%slpotwp)
            avail_h2o_lyr    = min( 1.0, max( 0.0, wilting_factor ) )                      &
                             * csite%soil_fracliq(k,ipa)                                   &
                             * ( soil(nsoil)%sfldcap - soil(nsoil)%soilwp )                &
                             * wdns * dslz(k)
            !------------------------------------------------------------------------------!



            !----- Add the factor from this layer to the integral. ------------------------!
            avail_h2o_coh(ico) = avail_h2o_coh(ico) + avail_h2o_lyr
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialize the array of maximum photosynthesis rates used in the mortality        !
   ! function.                                                                             !
   !---------------------------------------------------------------------------------------!
   csite%A_o_max(1:n_pft,ipa) = 0.0
   csite%A_c_max(1:n_pft,ipa) = 0.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the tallest cohort amongst all cohorts in this patch that is resolved        !
   ! (tuco).  In addition, we must find the tallest cohort for each PFT, so in case the    !
   ! we are using the light phenology, we use that value for Vm0 and leaf life span.       !
   !---------------------------------------------------------------------------------------!
   las         = .false.
   tuco_pft(:) = 0
   do ico = 1,cpatch%ncohorts
      ipft = cpatch%pft(ico)

      !----- If this is the tallest cohort to be used, we save its index. -----------------!
      if (.not. las .and. cpatch%leaf_resolvable(ico)) then
         las  = .true.
         tuco = ico
      end if
      !------------------------------------------------------------------------------------!



      !----- If this is the tallest cohort for this specific PFT, we save the index too. --!
      if (tuco_pft(ipft) == 0 .and. cpatch%leaf_resolvable(ico)) then
         tuco_pft(ipft) = ico
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    There is at least one cohort that meet requirements.  And this is tallest one, so  !
   ! we can use it to compute the maximum photosynthetic rates, i.e., the rate the cohort  !
   ! would have if it were at the top of the canopy.  This is used for the mortality       !
   ! function.                                                                             !
   !---------------------------------------------------------------------------------------!
   if (las) then
      !----- We now loop over PFTs, not cohorts, skipping those we are not using. ---------!
      do ipft = 1, n_pft

         if (include_pft(ipft)) then

            !------------------------------------------------------------------------------!
            !      Find the tallest cohort for this PFT.  In case the patch no longer has  !
            ! the PFT, then we just the default Vm0 and leaf life span.  This only matters !
            ! for light-controlled phenology, not the standard cases.                      !
            !------------------------------------------------------------------------------!
            tpft = tuco_pft(ipft)
            if (tpft == 0) then
               !---------------------------------------------------------------------------!
               !    This patch doesn't have any cohort of this PFT left, use default       !
               ! values.                                                                   !
               !---------------------------------------------------------------------------!
               vm0_tuco    = Vm0(ipft)
               if (leaf_turnover_rate(ipft) == 0.) then
                  llspan_tuco = llspan_inf
               else
                  llspan_tuco = 12. / leaf_turnover_rate(ipft)
               end if
            else
               !---------------------------------------------------------------------------!
               !    Use Vm0 and leaf life span of the tallest cohort of this PFT, so we    !
               ! avoid punishing or helping the plants too much in case the PFTs don't     !
               ! match.                                                                    !
               !---------------------------------------------------------------------------!
               vm0_tuco    = cpatch%vm_bar(tpft)
               llspan_tuco = cpatch%llspan(tpft)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Find the 100% relative humidity.  This is a temporary test to make the    !
            ! maximum carbon balance less negative.                                        !
            !------------------------------------------------------------------------------!
            can_ssh = csite%can_shv(ipa)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Call the photosynthesis for maximum photosynthetic rates.  The units      !
            ! of the input and output are the standard in most of ED modules, but many of  !
            ! them are converted inside the photosynthesis model.                          !
            !    Notice that the units that are per unit area are per m² of leaf, not the  !
            ! patch area.                                                                  !
            !------------------------------------------------------------------------------!
            call lphysiol_full(            & !
               csite%can_prss(ipa)         & ! Canopy air pressure              [       Pa]
             , csite%can_rhos(ipa)         & ! Canopy air density               [    kg/m³]
             , can_ssh                     & ! Canopy air sp. humidity          [    kg/kg]
             , csite%can_co2(ipa)          & ! Canopy air CO2 mixing ratio      [ µmol/mol]
             , ipft                        & ! Plant functional type            [      ---]
             , csite%par_l_max(ipa)        & ! Absorbed photos. active rad.     [ W/m²leaf]
             , cpatch%leaf_temp(tuco)      & ! Leaf temperature                 [        K]
             , cpatch%lint_shv(tuco)       & ! Leaf intercellular spec. hum.    [    kg/kg]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , llspan_tuco                 & ! Leaf life span                   [       yr]
             , vm0_tuco                    & ! Average Vm function              [µmol/m²/s]
             , cpatch%leaf_gbw(tuco)       & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , csite%A_o_max(ipft,ipa)     & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , csite%A_c_max(ipft,ipa)     & ! Photosynthesis rate     (closed) [µmol/m²/s]
             , d_A_light_max               & ! Photosynthesis rate     (light)  [µmol/m²/s]
             , d_A_rubp_max                & ! Photosynthesis rate     (RuBP)   [µmol/m²/s]
             , d_A_co2_max                 & ! Photosynthesis rate     (CO2)    [µmol/m²/s]
             , d_gsw_open                  & ! Stom. condct. of water  (open)   [  kg/m²/s]
             , d_gsw_closed                & ! Stom. condct. of water  (closed) [  kg/m²/s]
             , d_lsfc_shv_open             & ! Leaf sfc. sp. humidity  (open)   [    kg/kg]
             , d_lsfc_shv_closed           & ! Leaf sfc. sp. humidity  (closed) [    kg/kg]
             , d_lsfc_co2_open             & ! Leaf sfc. CO2 mix. rat. (open)   [ µmol/mol]
             , d_lsfc_co2_closed           & ! Leaf sfc. CO2 mix. rat. (closed) [ µmol/mol]
             , d_lint_co2_open             & ! Intercellular CO2       (open)   [ µmol/mol]
             , d_lint_co2_closed           & ! Intercellular CO2       (closed) [ µmol/mol]
             , leaf_resp                   & ! Leaf respiration rate            [µmol/m²/s]
             , vm                          & ! Max. capacity of Rubisco         [µmol/m²/s]
             , compp                       & ! Gross photo. compensation point  [ µmol/mol]
             , limit_flag                  & ! Photosynthesis limitation flag   [      ---]
             )
         end if
      end do
         
   else
      !---- There is no "active" cohort. --------------------------------------------------!
      csite%A_o_max(1:n_pft,ipa) = 0.0
      csite%A_c_max(1:n_pft,ipa) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Initialize some variables.                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Total root biomass (in kgC/m2) and patch sum available water. -------------------!
   pss_available_water = 0.0
   broot_tot           = 0.0
   !----- Initialize variables for transpiration calculation. -----------------------------!
   root_depth_indices(:) = .false.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts, from tallest to shortest.                                   !
   !---------------------------------------------------------------------------------------!
   cohortloop: do ico = 1,cpatch%ncohorts
         
      !------------------------------------------------------------------------------------!
      !     Only need to worry about photosyn if radiative transfer has been  done for     !
      ! this cohort.                                                                       !
      !------------------------------------------------------------------------------------!
      if (cpatch%leaf_resolvable(ico)) then

            !----- Alias for PFT and root layer. ------------------------------------------!
            ipft  = cpatch%pft(ico)
            kroot = cpatch%krdepth(ico)

            !------------------------------------------------------------------------------!
            !    Scale photosynthetically active radiation per unit of leaf.               !
            !------------------------------------------------------------------------------!
            leaf_par = cpatch%par_l(ico) / cpatch%lai(ico) 
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Call the photosynthesis for actual photosynthetic rates.  The units       !
            ! of the input and output are the standard in most of ED modules, but many of  !
            ! them are converted inside the photosynthesis model.                          !
            !    Notice that the units that are per unit area are per m² of leaf, not the  !
            ! patch area.                                                                  !
            !------------------------------------------------------------------------------!
            call lphysiol_full(            & !
               csite%can_prss(ipa)         & ! Canopy air pressure              [       Pa]
             , csite%can_rhos(ipa)         & ! Canopy air density               [    kg/m³]
             , csite%can_shv(ipa)          & ! Canopy air sp. humidity          [    kg/kg]
             , csite%can_co2(ipa)          & ! Canopy air CO2 mixing ratio      [ µmol/mol]
             , ipft                        & ! Plant functional type            [      ---]
             , leaf_par                    & ! Absorbed photos. active rad.     [ W/m²leaf]
             , cpatch%leaf_temp(ico)       & ! Leaf temperature                 [        K]
             , cpatch%lint_shv(ico)        & ! Leaf intercellular spec. hum.    [    kg/kg]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , cpatch%llspan(ico)          & ! Leaf life span                   [       yr]
             , cpatch%vm_bar(ico)          & ! Average Vm function              [µmol/m²/s]
             , cpatch%leaf_gbw(ico)        & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , cpatch%A_open(ico)          & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , cpatch%A_closed(ico)        & ! Photosynthesis rate     (closed) [µmol/m²/s]
             , cpatch%A_light(ico)         & ! Photosynthesis rate     (light)  [µmol/m²/s]
             , cpatch%A_rubp(ico)          & ! Photosynthesis rate     (RuBP)   [µmol/m²/s]
             , cpatch%A_co2(ico)           & ! Photosynthesis rate     (CO2)    [µmol/m²/s]
             , cpatch%gsw_open(ico)        & ! Stom. condct. of water  (open)   [  kg/m²/s]
             , cpatch%gsw_closed(ico)      & ! Stom. condct. of water  (closed) [  kg/m²/s]
             , cpatch%lsfc_shv_open(ico)   & ! Leaf sfc. sp. humidity  (open)   [    kg/kg] 
             , cpatch%lsfc_shv_closed(ico) & ! Leaf sfc. sp. humidity  (closed) [    kg/kg]
             , cpatch%lsfc_co2_open(ico)   & ! Leaf sfc. CO2 mix. rat. (open)   [ µmol/mol]
             , cpatch%lsfc_co2_closed(ico) & ! Leaf sfc. CO2 mix. rat. (closed) [ µmol/mol]
             , cpatch%lint_co2_open(ico)   & ! Intercellular CO2       (open)   [ µmol/mol]
             , cpatch%lint_co2_closed(ico) & ! Intercellular CO2       (closed) [ µmol/mol]
             , leaf_resp                   & ! Leaf respiration rate            [µmol/m²/s]
             , vm                          & ! Max. capacity of Rubisco         [µmol/m²/s]
             , compp                       & ! Gross photo. compensation point  [ µmol/mol]
             , limit_flag                  & ! Photosynthesis limitation flag   [      ---]
             )

            !----- Convert leaf respiration to [µmol/m²ground/s] --------------------------!
            cpatch%leaf_respiration(ico) = leaf_resp * cpatch%lai (ico)
            cpatch%today_leaf_resp(ico)  = cpatch%today_leaf_resp (ico)                    &
                                         + cpatch%leaf_respiration(ico)
            !----- The output variable must be in [kgC/plant/yr]. -------------------------!
            cpatch%fmean_leaf_resp(ico)  = cpatch%fmean_leaf_resp (ico)                    &
                                         + cpatch%leaf_respiration(ico)                    &
                                         * dtlsm_o_frqsum * umols_2_kgCyr                  &
                                         / cpatch%nplant          (ico)

            !----- Root biomass [kg/m2]. --------------------------------------------------!
            broot_loc = cpatch%broot(ico)  * cpatch%nplant(ico)

            !----- Supply of water. -------------------------------------------------------!
            cpatch%water_supply      (ico) = water_conductance       (ipft) * broot_loc    &
                                           * avail_h2o_coh            (ico)
            cpatch%fmean_water_supply(ico) = cpatch%fmean_water_supply(ico)                &
                                           + cpatch%water_supply      (ico)                &
                                           * dtlsm_o_frqsum

            root_depth_indices(kroot) = .true.
            broot_tot                 = broot_tot + broot_loc
            pss_available_water       = pss_available_water                                &
                                      + avail_h2o_coh(ico) * broot_loc

            !------------------------------------------------------------------------------!
            !     Determine the fraction of open stomata due to water limitation.          !
            ! This is a function of the ratio between the potential water demand           !
            ! (cpatch%psi_open, which is the average over the last time step), and the     !
            ! supply (cpatch%water_supply).                                                !
            !------------------------------------------------------------------------------!
            select case (h2o_plant_lim)
            case (0)
               !---- No water limitation, fsw is always 1.0. ------------------------------!
               cpatch%fsw(ico) = 1.0

            case (1,2)
               water_demand    = cpatch%psi_open(ico) * cpatch%lai(ico)
               if (cpatch%water_supply (ico) < tiny_num) then
                  cpatch%fsw(ico) = 0.0
               else
                  cpatch%fsw(ico) = 1.0 / (1.0 + water_demand / cpatch%water_supply(ico))
               end if
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Photorespiration can become important at high temperatures.  If so,     !
            ! close down the stomata.                                                      !
            !------------------------------------------------------------------------------!
            if (cpatch%A_open(ico) < cpatch%A_closed(ico)) then
               cpatch%fs_open(ico) = 0.0
            else
               cpatch%fs_open(ico) = cpatch%fsw(ico) * cpatch%fsn(ico)
            end if

            !----- Net stomatal conductance. ----------------------------------------------!
            cpatch%leaf_gsw(ico) =        cpatch%fs_open(ico)  * cpatch%gsw_open(ico)      &
                                 + (1.0 - cpatch%fs_open(ico)) * cpatch%gsw_closed(ico)
            !------------------------------------------------------------------------------!


            !----- GPP, averaged over frqstate. -------------------------------------------!
            cpatch%gpp(ico) = max(0., cpatch%lai(ico)                                      &
                                    * ( cpatch%fs_open(ico) * cpatch%A_open(ico)           &
                                      + (1. - cpatch%fs_open(ico)) * cpatch%A_closed(ico)) &
                                    + cpatch%leaf_respiration(ico) )
            !----- The average must be in [kgC/plant/yr]. ---------------------------------!
            cpatch%fmean_gpp(ico) = cpatch%fmean_gpp(ico)                                  &
                                  + cpatch%gpp      (ico) * umols_2_kgCyr * dtlsm_o_frqsum &
                                  / cpatch%nplant(ico)
            !------------------------------------------------------------------------------!


            !----- GPP, summed over 1 day. [µmol/m²ground] --------------------------------!
            cpatch%today_gpp(ico) = cpatch%today_gpp(ico) + cpatch%gpp(ico)
            !------------------------------------------------------------------------------!


            !----- Potential GPP if no N limitation. [µmol/m²ground] ----------------------!
            cpatch%today_gpp_pot(ico) = cpatch%today_gpp_pot(ico)                          &
                                      + cpatch%lai(ico)                                    &
                                      * ( cpatch%fsw(ico) * cpatch%A_open(ico)             &
                                        + (1.0 - cpatch%fsw(ico)) * cpatch%A_closed(ico))  &
                                      + cpatch%leaf_respiration(ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find the maximum productivities:                                         !
            !                                                                              !
            !     - today_gpp_lightmax: productivity of this cohort if it were at the top  !
            !                           of the canopy (full light), with the actual fsw.   !
            !     - today_gpp_moistmax: productivity of this cohort if the soil moisture   !
            !                           was such that fsw would be 1 (full moisture), with !
            !                           the actual light.                                  !
            !     - today_gpp_mlmax:    productivity of this cohort if it was at the top   !
            !                           of the canopy (full light) AND the soil moisture   !
            !                           was such that fsw would be 1 (full moisture).      !
            !                                                                              !
            !     These productivites are used to scale the relative carbon balance, which !
            ! will control density-dependent mortality.                                    !
            !------------------------------------------------------------------------------!
            cpatch%today_gpp_lightmax(ico) = cpatch%today_gpp_lightmax(ico)                &
                                           + cpatch%lai(ico)                               &
                                           * ( cpatch%fs_open(ico)                         &
                                             * csite%A_o_max(ipft,ipa)                     &
                                             + (1.0 - cpatch%fs_open(ico))                 &
                                             * csite%A_c_max(ipft,ipa) )                   &
                                           + cpatch%leaf_respiration(ico)
            cpatch%today_gpp_moistmax(ico) = cpatch%today_gpp_moistmax(ico)                &
                                           + cpatch%lai(ico) * cpatch%A_open(ico)          &
                                           + cpatch%leaf_respiration(ico)
            cpatch%today_gpp_mlmax(ico)    = cpatch%today_gpp_mlmax(ico)                   &
                                           + cpatch%lai(ico) * csite%A_o_max(ipft,ipa)     &
                                           + cpatch%leaf_respiration(ico)
            !------------------------------------------------------------------------------!

      else
         !----- If the cohort wasn't solved, we must assign some zeroes. ------------------!
         cpatch%A_open(ico)               = 0.0
         cpatch%A_closed(ico)             = 0.0
         cpatch%psi_open(ico)             = 0.0
         cpatch%psi_closed(ico)           = 0.0
         cpatch%water_supply(ico)         = 0.0
         cpatch%gsw_open(ico)             = 0.0
         cpatch%gsw_closed(ico)           = 0.0
         cpatch%leaf_gbh(ico)             = 0.0
         cpatch%leaf_gbw(ico)             = 0.0
         cpatch%leaf_gsw(ico)             = 0.0
         cpatch%gpp(ico)                  = 0.0
         cpatch%leaf_respiration(ico)     = 0.0
         vm                               = 0.0
         limit_flag                       = 0
      end if
      
      !------------------------------------------------------------------------------------!
      !    Not really a part of the photosynthesis scheme, but this will do it.  We must   !
      ! integrate the "mean" of the remaining respiration terms, except for the root one.  !
      ! This is done regardless on whether the cohort is doing photosynthesis.             !
      !                                                                                    !
      !    The "_respiration(ico) terms are in kgC/plant/day, so we must also multiply     !
      ! them by the number of years per day so the output is in kgC/plant/yr.  High time   !
      ! we switched everything to SI...                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%fmean_growth_resp (ico) = cpatch%fmean_growth_resp  (ico)                     &
                                     + cpatch%growth_respiration (ico) * dtlsm_o_frqsum    &
                                     * yr_day
      cpatch%fmean_storage_resp(ico) = cpatch%fmean_storage_resp (ico)                     &
                                     + cpatch%storage_respiration(ico) * dtlsm_o_frqsum    &
                                     * yr_day
      cpatch%fmean_vleaf_resp  (ico) = cpatch%fmean_vleaf_resp   (ico)                     &
                                     + cpatch%vleaf_respiration  (ico) * dtlsm_o_frqsum    &
                                     * yr_day
      !------------------------------------------------------------------------------------!

      if (print_photo_debug) then
         call print_photo_details(cmet,csite,ipa,ico,limit_flag,vm,compp)
      end if
   end do cohortloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Add the contribution of this time step to the average available water.  This is   !
   ! done only if there is some cohort transpiring.                                        !
   !---------------------------------------------------------------------------------------!
   if (broot_tot > 1.e-20) then
      csite%fmean_available_water(ipa) = csite%fmean_available_water(ipa)                  &
                                       + pss_available_water * dtlsm_o_frqsum / broot_tot
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    De-allocate the temporary vector.                                                  !
   !---------------------------------------------------------------------------------------!
   if (cpatch%ncohorts > 0) then
      deallocate(avail_h2o_coh)
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine canopy_photosynthesis
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine prints some extra information on the photosynthesis driver in a     !
! convenient ascii file for debugging purposes.                                            !
!------------------------------------------------------------------------------------------!
subroutine print_photo_details(cmet,csite,ipa,ico,limit_flag,vm,compp)
   use ed_max_dims    , only : str_len            ! ! intent(in)
   use ed_state_vars  , only : sitetype           & ! structure
                             , patchtype          ! ! structure
   use met_driver_coms, only : met_driv_state     ! ! structure
   use physiology_coms, only : photo_prefix       ! ! intent(in)
   use ed_misc_coms   , only : current_time       ! ! intent(in)
   use consts_coms    , only : Watts_2_Ein        & ! intent(in)
                             , mol_2_umol         & ! intent(in)
                             , t008               ! ! intent(in)
   use pft_coms       , only : quantum_efficiency & ! intent(in)
                             , photosyn_pathway   ! ! intent(in)
   use physiology_coms, only : quantum_efficiency_T ! ! intent(in)
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite           ! Current site
   type(met_driv_state)      , target      :: cmet            ! Current met. conditions.
   integer                   , intent(in)  :: ipa             ! Current patch number
   integer                   , intent(in)  :: ico             ! Current cohort number
   integer                   , intent(in)  :: limit_flag      ! Limitation flag
   real                      , intent(in)  :: vm              ! Maximum Rubisco capacity
   real                      , intent(in)  :: compp           ! GPP compensation point
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype)           , pointer     :: jpatch          ! Current site
   type(patchtype)           , pointer     :: cpatch          ! Current site
   character(len=str_len)                  :: photo_fout      ! File with the cohort info
   integer                                 :: ipft
   integer                                 :: jpa
   integer                                 :: jco
   logical                                 :: isthere
   real                                    :: leaf_resp
   real                                    :: par_area
   real                                    :: nir_area
   real                                    :: parv
   real                                    :: nirv
   real                                    :: util_parv
   real                                    :: alpha
   !----- Local constants. ----------------------------------------------------------------!
   character(len=10), parameter :: hfmt='(63(a,1x))'
   character(len=48), parameter :: bfmt='(3(i13,1x),1(es13.6,1x),2(i13,1x),57(es13.6,1x))'
   !----- Locally saved variables. --------------------------------------------------------!
   logical                   , save        :: first_time=.true.
   !---------------------------------------------------------------------------------------!


   !----- Make some aliases. --------------------------------------------------------------!
   cpatch      => csite%patch(ipa)

   ipft        =  cpatch%pft             (ico)
   leaf_resp   =  cpatch%leaf_respiration(ico)
   !---------------------------------------------------------------------------------------!

   if (cpatch%leaf_resolvable(ico)) then
      par_area   = cpatch%par_l(ico) * Watts_2_Ein * mol_2_umol
      parv       = par_area / cpatch%lai(ico)
      nir_area   = (cpatch%rshort_l(ico) - cpatch%par_l(ico)) * Watts_2_Ein * mol_2_umol
      nirv       = nir_area / cpatch%lai(ico)
      
      !------------------------------------------------------------------------------------!
      !    Is alpha (quantum efficiency) temperature dependent?  If so, calculate after    !
      !    Ehlringer and Ollebjorkman 1977, if not use default value from ed_params                                                   !
      !------------------------------------------------------------------------------------!
      select case(quantum_efficiency_T)
      case(1)
           select case (photosyn_pathway(ipft))
           case (4)
               alpha         = dble(quantum_efficiency(ipft))       
           case (3)       
               alpha         = dble(-0.0016*(dble(cpatch%leaf_temp(ico))-t008)+0.1040)
           end select
      case default
            alpha         = dble(quantum_efficiency(ipft))      
      end select
      
      util_parv  = alpha * parv
   else
      par_area  = 0.0
      parv      = 0.0
      nir_area  = 0.0
      nirv      = 0.0
      util_parv = 0.0
   end if

   !---------------------------------------------------------------------------------------!
   !     First time here.  Delete all files.                                               !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      do jpa = 1, csite%npatches
         jpatch => csite%patch(jpa)
         do jco = 1, jpatch%ncohorts
            write (photo_fout,fmt='(a,2(a,i4.4),a)')                                       &
                  trim(photo_prefix),'patch_',jpa,'_cohort_',jco,'.txt'
            inquire(file=trim(photo_fout),exist=isthere)
            if (isthere) then
               !---- Open the file to delete when closing. --------------------------------!
               open (unit=57,file=trim(photo_fout),status='old',action='write')
               close(unit=57,status='delete')
            end if
         end do
      end do
      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!




   !----- Create the file name. -----------------------------------------------------------!
   write (photo_fout,fmt='(a,2(a,i4.4),a)') trim(photo_prefix),'patch_',ipa                &
                                                              ,'_cohort_',ico,'.txt'
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check whether the file exists or not.  In case it doesn't, create it and add the   !
   ! header.                                                                               !
   !---------------------------------------------------------------------------------------!
   inquire(file=trim(photo_fout),exist=isthere)
   if (.not. isthere) then
      open  (unit=57,file=trim(photo_fout),status='replace',action='write')
      write (unit=57,fmt=hfmt)   '         YEAR', '        MONTH', '          DAY'         &
                               , '         TIME', '          PFT', '   LIMIT_FLAG'         &
                               , '       HEIGHT', '       NPLANT', '        BLEAF'         &
                               , '          LAI', '    LEAF_HCAP', '   LEAF_WATER'         &
                               , '    LEAF_TEMP', '    WOOD_TEMP', '     CAN_TEMP'         &
                               , '     ATM_TEMP', '  GROUND_TEMP', '      CAN_SHV'         &
                               , '      ATM_SHV', '   GROUND_SHV', 'LSFC_SHV_OPEN'         &
                               , 'LSFC_SHV_CLOS', '     LINT_SHV', '     ATM_PRSS'         &
                               , '     CAN_PRSS', '         PCPG', '     CAN_RHOS'         &
                               , '      ATM_CO2', '      CAN_CO2', 'LSFC_CO2_OPEN'         &
                               , 'LSFC_CO2_CLOS', 'LINT_CO2_OPEN', 'LINT_CO2_CLOS'         &
                               , '        COMPP', '     PAR_AREA', '         PARV'         &
                               , '    UTIL_PARV', '     NIR_AREA', '         NIRV'         &
                               , '          GPP', '    LEAF_RESP', '     LEAF_GBH'         &
                               , '     LEAF_GBW', '     WOOD_GBH', '     WOOD_GBW'         &
                               , '     LEAF_GSW', '       A_OPEN', '       A_CLOS'         &
                               , '      A_LIGHT', '       A_RUBP', '        A_CO2'         &
                               , '     GSW_OPEN', '     GSW_CLOS', '     PSI_OPEN'         &
                               , '     PSI_CLOS', '   H2O_SUPPLY', '          FSW'         &
                               , '          FSN', '      FS_OPEN', '     ATM_WIND'         &
                               , '     VEG_WIND', '        USTAR', '           VM'
                               
                              
      close (unit=57,status='keep')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Re-open the file at the last line, and include the current status.                !
   !---------------------------------------------------------------------------------------!
   open (unit=57,file=trim(photo_fout),status='old',action='write',position='append')
   write(unit=57,fmt=bfmt)                                                                 &
     current_time%year          , current_time%month         , current_time%date           &
   , current_time%time          , cpatch%pft(ico)            , limit_flag                  &
   , cpatch%hite(ico)           , cpatch%nplant(ico)         , cpatch%bleaf(ico)           &
   , cpatch%lai(ico)            , cpatch%leaf_hcap(ico)      , cpatch%leaf_water(ico)      &
   , cpatch%leaf_temp(ico)      , cpatch%wood_temp(ico)      , csite%can_temp(ipa)         &
   , cmet%atm_tmp               , csite%ground_temp(ipa)     , csite%can_shv(ipa)          &
   , cmet%atm_shv               , csite%ground_shv(ipa)      , cpatch%lsfc_shv_open(ico)   &
   , cpatch%lsfc_shv_closed(ico), cpatch%lint_shv(ico)       , cmet%prss                   &
   , csite%can_prss(ipa)        , cmet%pcpg                  , csite%can_rhos(ipa)         &
   , cmet%atm_co2               , csite%can_co2(ipa)         , cpatch%lsfc_co2_open(ico)   &
   , cpatch%lsfc_co2_closed(ico), cpatch%lint_co2_open(ico)  , cpatch%lint_co2_closed(ico) &
   , compp                      , par_area                   , parv                        &
   , util_parv                  , nir_area                   , nirv                        &
   , cpatch%gpp(ico)            , leaf_resp                  , cpatch%leaf_gbh(ico)        &
   , cpatch%leaf_gbw(ico)       , cpatch%wood_gbh(ico)       , cpatch%wood_gbw(ico)        &
   , cpatch%leaf_gsw(ico)       , cpatch%A_open(ico)         , cpatch%A_closed(ico)        &
   , cpatch%A_light(ico)        , cpatch%A_rubp(ico)         , cpatch%A_co2   (ico)        &
   , cpatch%gsw_open(ico)       , cpatch%gsw_closed(ico)     , cpatch%psi_open(ico)        &
   , cpatch%psi_closed(ico)     , cpatch%water_supply(ico)   , cpatch%fsw(ico)             &
   , cpatch%fsn(ico)            , cpatch%fs_open(ico)        , cmet%vels                   &
   , cpatch%veg_wind(ico)       , csite%ustar(ipa)           , vm

   close(unit=57,status='keep')
   !---------------------------------------------------------------------------------------!

   return
end subroutine print_photo_details
!==========================================================================================!
!==========================================================================================!
