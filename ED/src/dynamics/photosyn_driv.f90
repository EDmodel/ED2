!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the photosynthesis scheme (Farquar and Leuning).  This  !
! is called every step, but not every sub-step.                                            !
!------------------------------------------------------------------------------------------!
subroutine canopy_photosynthesis(csite,cmet,mzg,ipa,ed_ktrans,lsl                          &
                                ,leaf_aging_factor,green_leaf_factor)
   use ed_state_vars  , only : sitetype          & ! structure
                             , patchtype         ! ! structure
   use ed_max_dims    , only : n_pft             ! ! intent(in)
   use pft_coms       , only : leaf_width        & ! intent(in)
                             , water_conductance & ! intent(in)
                             , q                 & ! intent(in)
                             , qsw               & ! intent(in)
                             , include_pft       ! ! intent(in)
   use soil_coms      , only : soil              & ! intent(in)
                             , dslz              ! ! intent(in)
   use consts_coms    , only : t00               & ! intent(in)
                             , epi               & ! intent(in)
                             , wdnsi             & ! intent(in)
                             , wdns              & ! intent(in)
                             , kgCday_2_umols    ! ! intent(in)
   use met_driver_coms, only : met_driv_state    ! ! structure
   use physiology_coms, only : print_photo_debug & ! intent(in)
                             , new_fsw_method    ! ! intent(in)
   use farq_leuning   , only : lphysiol_full     ! ! sub-routine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite             ! Current site
   type(met_driv_state)      , target      :: cmet              ! Current met. conditions.
   integer                   , intent(in)  :: ipa               ! Current patch #
   integer                   , intent(in)  :: lsl               ! Lowest soil level
   integer                   , intent(in)  :: mzg               ! Number of soil layers
   real   , dimension(n_pft) , intent(in)  :: leaf_aging_factor ! 
   real   , dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
   integer, dimension(mzg)   , intent(out) :: ed_ktrans         ! 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch             ! Current site
   integer                                 :: ico                ! Current cohort #
   integer                                 :: tuco               ! Tallest used cohort
   integer                                 :: ipft
   integer                                 :: k1
   integer                                 :: k2
   integer                                 :: nsoil
   integer                                 :: limit_flag
   logical, dimension(mzg)                 :: root_depth_indices ! 
   logical                                 :: las
   real   , dimension(mzg)                 :: available_liquid_water
   real                                    :: leaf_par
   real                                    :: leaf_resp
   real                                    :: d_gsw_open
   real                                    :: d_gsw_closed
   real                                    :: d_lsfc_shv_open
   real                                    :: d_lsfc_shv_closed
   real                                    :: d_lsfc_co2_open
   real                                    :: d_lsfc_co2_closed
   real                                    :: d_lint_co2_open
   real                                    :: d_lint_co2_closed
   real                                    :: slpotv
   real                                    :: swp
   real                                    :: vm
   real                                    :: compp
   real                                    :: broot_tot
   real                                    :: broot_loc
   real                                    :: pss_available_water
   !---------------------------------------------------------------------------------------!


   !----- Pointing to the cohort structures -----------------------------------------------!
   cpatch => csite%patch(ipa)

   !----- Finding the patch-level Total Leaf and Wood Area Index. -------------------------!
   csite%lai(ipa) = 0.0
   csite%wpa(ipa) = 0.0
   csite%wai(ipa) = 0.0
   do ico=1,cpatch%ncohorts
      csite%lai(ipa)  = csite%lai(ipa)  + cpatch%lai(ico)
      csite%wpa(ipa)  = csite%wpa(ipa)  + cpatch%wpa(ico)
      csite%wai(ipa)  = csite%wai(ipa)  + cpatch%wai(ico)
   end do


   !----- Calculate liquid water available for transpiration. -----------------------------!
   nsoil = csite%ntext_soil(mzg,ipa)
   available_liquid_water(mzg) = wdns * dslz(mzg) * csite%soil_fracliq(mzg,ipa)            &
                               * max(0.0, csite%soil_water(mzg,ipa) - soil(nsoil)%soilwp )
   do k1 = mzg-1, lsl, -1
      nsoil = csite%ntext_soil(k1,ipa)
      available_liquid_water(k1) = available_liquid_water(k1+1)                            &
                                 + wdns * dslz(k1) * csite%soil_fracliq(k1,ipa)            &
                                 * max(0.0, csite%soil_water(k1,ipa) - soil(nsoil)%soilwp )
   end do
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Initialize the array of maximum photosynthesis rates used in the mortality        !
   ! function.                                                                             !
   !---------------------------------------------------------------------------------------!
   csite%A_o_max(1:n_pft,ipa) = 0.0
   csite%A_c_max(1:n_pft,ipa) = 0.0

   !---------------------------------------------------------------------------------------!
   !     Find the tallest cohort with TAI above minimum, sufficient heat capacity, and not !
   ! buried in snow.  The first two conditions are redundant, but we will keep them for    !
   ! the time being, so it is going to be safer.                                           !
   !---------------------------------------------------------------------------------------!
   las = .false.
   do ico = 1,cpatch%ncohorts
      !----- If this is the tallest cohort to be used, we save its index. -----------------!
      if (.not. las .and. cpatch%solvable(ico)) then
         las  = .true.
         tuco = ico
      end if
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
         if (include_pft(ipft) == 1)then

            !------------------------------------------------------------------------------!
            !    Scale photosynthetically active radiation per unit of leaf.               !
            !------------------------------------------------------------------------------!
            leaf_par = cpatch%par_v(tuco) / cpatch%lai(tuco)

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
             , csite%can_shv(ipa)          & ! Canopy air sp. humidity          [    kg/kg]
             , csite%can_co2(ipa)          & ! Canopy air CO2 mixing ratio      [ µmol/mol]
             , ipft                        & ! Plant functional type            [      ---]
             , leaf_par                    & ! Absorbed photos. active rad.     [     W/m²]
             , cpatch%veg_temp(tuco)       & ! Vegetation temperature           [        K]
             , cpatch%lint_shv(tuco)       & ! Leaf intercellular spec. hum.    [    kg/kg]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , cpatch%llspan(tuco)         & ! Leaf life span                   [       yr]
             , cpatch%vm_bar(tuco)         & ! Average Vm function              [µmol/m²/s]
             , cpatch%gbw(tuco)            & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , csite%A_o_max(ipft,ipa)     & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , csite%A_c_max(ipft,ipa)     & ! Photosynthesis rate     (closed) [µmol/m²/s]
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
             , csite%old_stoma_data_max(ipft,ipa) & ! Previous state            [      ---]
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
      if (cpatch%solvable(ico)) then

            !----- Alias for PFT ----------------------------------------------------------!
            ipft = cpatch%pft(ico)

            !------------------------------------------------------------------------------!
            !    Scale photosynthetically active radiation per unit of leaf.               !
            !------------------------------------------------------------------------------!
            leaf_par = cpatch%par_v(ico) / cpatch%lai(ico) 


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
             , leaf_par                    & ! Absorbed photos. active rad.     [     W/m²]
             , cpatch%veg_temp(ico)        & ! Vegetation temperature           [        K]
             , cpatch%lint_shv(ico)        & ! Leaf intercellular spec. hum.    [    kg/kg]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , cpatch%llspan(ico)          & ! Leaf life span                   [       yr]
             , cpatch%vm_bar(ico)          & ! Average Vm function              [µmol/m²/s]
             , cpatch%gbw(ico)             & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , cpatch%A_open(ico)          & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , cpatch%A_closed(ico)        & ! Photosynthesis rate     (closed) [µmol/m²/s]
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
             , csite%old_stoma_data_max(ipft,ipa) & ! Previous state            [      ---]
             )

            !----- Convert leaf respiration to [µmol/m²ground/s] --------------------------!
            cpatch%leaf_respiration(ico) = leaf_resp * cpatch%lai(ico)
            cpatch%mean_leaf_resp(ico)   = cpatch%mean_leaf_resp(ico)                      &
                                         + cpatch%leaf_respiration(ico)
            cpatch%today_leaf_resp(ico)  = cpatch%today_leaf_resp(ico)                     &
                                         + cpatch%leaf_respiration(ico)

            !----- Root biomass [kg/m2]. --------------------------------------------------!
            broot_loc = q(ipft) * cpatch%balive(ico)                                       &
                      / (1.0 + q(ipft) + cpatch%hite(ico) * qsw(ipft) )                    &
                      * cpatch%nplant(ico)

            !----- Supply of water. -------------------------------------------------------!
            cpatch%water_supply(ico) = water_conductance(ipft)                             &
                                     * available_liquid_water(cpatch%krdepth(ico))         &
                                     * broot_loc

            root_depth_indices(cpatch%krdepth(ico)) = .true.
            broot_tot = broot_tot + broot_loc
            pss_available_water = pss_available_water                                      &
                                + available_liquid_water(cpatch%krdepth(ico)) * broot_loc

            !------------------------------------------------------------------------------!
            !     Determine the fraction of open stomata due to water limitation.          !
            ! This is a function of the ratio between the potential water demand           !
            ! (cpatch%psi_open, which is the average over the last time step), and the     !
            ! supply (cpatch%water_supply).                                                !
            !------------------------------------------------------------------------------!
            if (new_fsw_method) then
               if (cpatch%water_supply(ico) > 1.e-20) then
                  cpatch%fsw(ico) = 0.5 * ( 1. - tanh( 2. * (cpatch%psi_open(ico)          &
                                                            /cpatch%water_supply(ico)-2.)))
               else
                  cpatch%fsw(ico) = 0.0
               end if
            else
               cpatch%fsw(ico) = cpatch%water_supply(ico)                                  &
                               / max( 1.0e-20                                              &
                                    , cpatch%water_supply(ico) + cpatch%psi_open(ico))
            end if
            ! cpatch%fsw(ico) = 1.0
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
            cpatch%stomatal_conductance(ico) =  cpatch%fs_open(ico) *cpatch%gsw_open(ico)  &
                                             + (1.0 - cpatch%fs_open(ico))                 &
                                             * cpatch%gsw_closed(ico)

            !----- GPP, averaged over frqstate. -------------------------------------------!
            cpatch%gpp(ico)       = cpatch%lai(ico)                                        &
                                  * ( cpatch%fs_open(ico) * cpatch%A_open(ico)             &
                                    + (1.0 - cpatch%fs_open(ico)) * cpatch%A_closed(ico) ) &
                                  + cpatch%leaf_respiration(ico)
            cpatch%mean_gpp(ico)  = cpatch%mean_gpp(ico) + cpatch%gpp(ico)

            !----- GPP, summed over 1 day. [µmol/m²ground] --------------------------------!
            cpatch%today_gpp(ico) = cpatch%today_gpp(ico) + cpatch%gpp(ico)

            !----- Potential GPP if no N limitation. [µmol/m²ground] ----------------------!
            cpatch%today_gpp_pot(ico) = cpatch%today_gpp_pot(ico)                          &
                                      + cpatch%lai(ico)                                    &
                                      * ( cpatch%fsw(ico) * cpatch%A_open(ico)             &
                                        + (1.0 - cpatch%fsw(ico)) * cpatch%A_closed(ico))  &
                                      + cpatch%leaf_respiration(ico)

            !----- Maximum GPP if at the top of the canopy [µmol/m²ground] ----------------!
            cpatch%today_gpp_max(ico) = cpatch%today_gpp_max(ico)                          &
                                      + cpatch%lai(ico)                                    &
                                      * ( cpatch%fs_open(ico) * csite%A_o_max(ipft,ipa)    &
                                        + (1.0 - cpatch%fs_open(ico))                      &
                                          * csite%A_c_max(ipft,ipa))                       &
                                      + cpatch%leaf_respiration(ico)

      else
         !----- If the cohort wasn't solved, we must assign some zeroes. ------------------!
         cpatch%A_open(ico)               = 0.0
         cpatch%A_closed(ico)             = 0.0
         cpatch%psi_open(ico)             = 0.0
         cpatch%psi_closed(ico)           = 0.0
         cpatch%water_supply(ico)         = 0.0
         cpatch%gsw_open(ico)             = 0.0
         cpatch%gsw_closed(ico)           = 0.0
         cpatch%gbh(ico)                  = 0.0
         cpatch%gbw(ico)                  = 0.0
         cpatch%stomatal_conductance(ico) = 0.0
         cpatch%gpp(ico)                  = 0.0
         cpatch%leaf_respiration(ico)     = 0.0
         vm                               = 0.0
         limit_flag                       = 0
      end if
      
      !------------------------------------------------------------------------------------!
      !    Not really a part of the photosynthesis scheme, but this will do it.  We must   !
      ! integrate the "mean" of the remaining respiration terms, except for the root one.  !
      ! This is done regardless on whether the cohort is doing photosynthesis.  Also, we   !
      ! convert units so all fast respiration terms are in [µmol/m²ground/s].              !
      !------------------------------------------------------------------------------------!
      cpatch%mean_growth_resp (ico) = cpatch%mean_growth_resp (ico)                        &
                                    + cpatch%growth_respiration (ico) * kgCday_2_umols     &
                                    * cpatch%nplant(ico)
      cpatch%mean_storage_resp(ico) = cpatch%mean_storage_resp(ico)                        &
                                    + cpatch%storage_respiration(ico) * kgCday_2_umols     &
                                    * cpatch%nplant(ico)
      cpatch%mean_vleaf_resp  (ico) = cpatch%mean_vleaf_resp  (ico)                        &
                                    + cpatch%vleaf_respiration  (ico) * kgCday_2_umols     &
                                    * cpatch%nplant(ico)                                    
      !------------------------------------------------------------------------------------!

      if (print_photo_debug) then
         call print_photo_details(cmet,csite,ipa,ico,limit_flag,vm,compp)
      end if
   end do cohortloop

   !---------------------------------------------------------------------------------------!
   !     Add the contribution of this time step to the average available water.            !
   !---------------------------------------------------------------------------------------!
   if (broot_tot > 1.e-20) then
      csite%avg_available_water(ipa) = csite%avg_available_water(ipa)                      &
                                     + pss_available_water / broot_tot
   !else
   !  Add nothing, the contribution of this time is zero since no cohort can transpire... 
   end if

   !---------------------------------------------------------------------------------------!
   !     For plants of a given rooting depth, determine soil level from which transpired   !
   ! water is to be extracted.                                                             !
   !---------------------------------------------------------------------------------------!
   ed_ktrans(:) = 0
   do k1 = lsl, mzg
      !---- Assign a very large negative, so it will update it at least once. -------------!
      swp = -huge(1.)
      if (root_depth_indices(k1)) then
         do k2 = k1, mzg
            nsoil = csite%ntext_soil(k2,ipa)
            !------------------------------------------------------------------------------!
            !      Find slpotv using the available liquid water, since ice is unavailable  !
            ! for transpiration.                                                           !
            !------------------------------------------------------------------------------!
            slpotv = soil(nsoil)%slpots * csite%soil_fracliq(k2,ipa)                       &
                   * (soil(nsoil)%slmsts / csite%soil_water(k2,ipa)) ** soil(nsoil)%slbs

            !------------------------------------------------------------------------------!
            !      Find layer in root zone with highest slpotv AND soil_water above        !
            ! minimum soilwp.  Set ktrans to this layer.                                   !
            !------------------------------------------------------------------------------!
            if (slpotv > swp .and. csite%soil_water(k2,ipa) > soil(nsoil)%soilwp) then
               swp = slpotv
               ed_ktrans(k1) = k2
            end if
         end do
      end if
   end do

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
                             , mol_2_umol         ! ! intent(in)
   use pft_coms       , only : quantum_efficiency ! ! intent(in)
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
   real                                    :: stom_condct
   real                                    :: par_area
   real                                    :: parv
   real                                    :: util_parv
   !----- Local constants. ----------------------------------------------------------------!
   character(len=10), parameter :: hfmt='(55(a,1x))'
   character(len=48), parameter :: bfmt='(3(i13,1x),1(es13.6,1x),2(i13,1x),49(es13.6,1x))'
   !----- Locally saved variables. --------------------------------------------------------!
   logical                   , save        :: first_time=.true.
   !---------------------------------------------------------------------------------------!


   !----- Make some aliases. --------------------------------------------------------------!
   cpatch      => csite%patch(ipa)
   ipft        =  cpatch%pft(ico)
   leaf_resp   =  cpatch%leaf_respiration(ico)
   stom_condct =  cpatch%stomatal_conductance(ico)
   !---------------------------------------------------------------------------------------!

   if (cpatch%solvable(ico)) then
      par_area   = cpatch%par_v(ico) * Watts_2_Ein * mol_2_umol
      parv       = par_area / cpatch%lai(ico)
      util_parv  = quantum_efficiency(ipft) * parv
   else
      par_area  = 0.0
      parv      = 0.0
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
                               , '          LAI', '      HCAPVEG', '    VEG_WATER'         &
                               , '     VEG_TEMP', '     CAN_TEMP', '     ATM_TEMP'         &
                               , '  GROUND_TEMP', '      CAN_SHV', '      ATM_SHV'         &
                               , '   GROUND_SHV', 'LSFC_SHV_OPEN', 'LSFC_SHV_CLOS'         &
                               , '     LINT_SHV', '     ATM_PRSS', '     CAN_PRSS'         &
                               , '         PCPG', '     CAN_RHOS', '      ATM_CO2'         &
                               , '      CAN_CO2', 'LSFC_CO2_OPEN', 'LSFC_CO2_CLOS'         &
                               , 'LINT_CO2_OPEN', 'LINT_CO2_CLOS', '        COMPP'         &
                               , '     PAR_AREA', '         PARV', '    UTIL_PARV'         &
                               , '          GPP', '    LEAF_RESP', '          GBH'         &
                               , '          GBW', '  STOM_CONDCT', '       A_OPEN'         &
                               , '       A_CLOS', '     GSW_OPEN', '     GSW_CLOS'         &
                               , '     PSI_OPEN', '     PSI_CLOS', '   H2O_SUPPLY'         &
                               , '          FSW', '          FSN', '      FS_OPEN'         &
                               , '     ATM_WIND', '     VEG_WIND', '        USTAR'         &
                               , '           VM'
                              
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
   , cpatch%lai(ico)            , cpatch%hcapveg(ico)        , cpatch%veg_water(ico)       &
   , cpatch%veg_temp(ico)       , csite%can_temp(ipa)        , cmet%atm_tmp                &
   , csite%ground_temp(ipa)     , csite%can_shv(ipa)         , cmet%atm_shv                &
   , csite%ground_shv(ipa)      , cpatch%lsfc_shv_open(ico)  , cpatch%lsfc_shv_closed(ico) &
   , cpatch%lint_shv(ico)       , cmet%prss                  , csite%can_prss(ipa)         &
   , cmet%pcpg                  , csite%can_rhos(ipa)        , cmet%atm_co2                &
   , csite%can_co2(ipa)         , cpatch%lsfc_co2_open(ico)  , cpatch%lsfc_co2_closed(ico) &
   , cpatch%lint_co2_open(ico)  , cpatch%lint_co2_closed(ico), compp                       &
   , par_area                   , parv                       , util_parv                   &
   , cpatch%gpp(ico)            , leaf_resp                  , cpatch%gbh(ico)             &
   , cpatch%gbw(ico)            , stom_condct                , cpatch%A_open(ico)          &
   , cpatch%A_closed(ico)       , cpatch%gsw_open(ico)       , cpatch%gsw_closed(ico)      &
   , cpatch%psi_open(ico)       , cpatch%psi_closed(ico)     , cpatch%water_supply(ico)    &
   , cpatch%fsw(ico)            , cpatch%fsn(ico)            , cpatch%fs_open(ico)         &
   , cmet%vels                  , cpatch%veg_wind(ico)       , csite%ustar(ipa)            &
   , vm
   

   close(unit=57,status='keep')
   !---------------------------------------------------------------------------------------!

   return
end subroutine print_photo_details
!==========================================================================================!
!==========================================================================================!
