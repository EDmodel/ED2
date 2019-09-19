!==========================================================================================!
!==========================================================================================!
!     Module that controls changes in leaf biomass due to phenology.                       !
!------------------------------------------------------------------------------------------!
module phenology_driv
   contains

   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine controls the changes in leaf biomass due to phenology.             !
   !---------------------------------------------------------------------------------------!
   subroutine phenology_driver(cgrid, doy, month, tfact,veget_dyn_on)
      use ed_state_vars  , only : edtype                & ! structure
                                , polygontype           & ! structure
                                , sitetype              ! ! structure
      use phenology_coms , only : iphen_scheme          ! ! intent(in)
      use ed_misc_coms   , only : current_time          ! ! intent(in)
      use phenology_aux  , only : prescribed_leaf_state & ! subroutine
                                , update_thermal_sums   & ! subroutine
                                , update_turnover       ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target      :: cgrid
      integer           , intent(in)  :: doy
      integer           , intent(in)  :: month
      real              , intent(in)  :: tfact
      logical           , intent(in)  :: veget_dyn_on
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer     :: cpoly
      type(sitetype)    , pointer     :: csite
      integer                         :: ipy
      integer                         :: isi
      integer                         :: ipa
      !------------------------------------------------------------------------------------!

      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !     Get the patch-level average daily temperature, which is needed for       !
            ! mortality, recruitment and some phenology schemes.                           !
            !------------------------------------------------------------------------------!
            do ipa = 1,csite%npatches
               csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) * tfact
            end do
            
            select case (iphen_scheme)
            case (-1,0,2)
               !---------------------------------------------------------------------------!
               !     Default predictive scheme (Botta et al.) or the modified drought      !
               ! deciduous phenology for broadleaf PFTs.                                   !
               !---------------------------------------------------------------------------!
               call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)
               
            case (1)
               !----- Use prescribed phenology. -------------------------------------------!
               call prescribed_leaf_state(cgrid%lat(ipy), current_time%month               &
                                         ,current_time%year, doy                           &
                                         ,cpoly%green_leaf_factor(:,isi)                   &
                                         ,cpoly%leaf_aging_factor(:,isi)                   &
                                         ,cpoly%phen_pars(isi) ) 
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)


            case (3,4)
               !----- Light-controlled predictive phenology scheme. -----------------------!
               call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
               call update_turnover(cpoly,isi)
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)
            end select
         end do
      end do

      return
   end subroutine phenology_driver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine update_phenology(doy, cpoly, isi, lat,veget_dyn_on)
      use ed_state_vars  , only : polygontype              & ! structure
                                , sitetype                 & ! structure
                                , patchtype                ! ! structure
      use grid_coms      , only : nzg                      ! ! intent(in)
      use pft_coms       , only : phenology                & ! intent(in)
                                , c2n_leaf                 & ! intent(in)
                                , q                        & ! intent(in)
                                , leaf_psi_tlp             & ! intent(in)
                                , high_psi_threshold       & ! intent(in)
                                , low_psi_threshold        & ! intent(in)
                                , leaf_shed_rate           & ! intent(in)
                                , leaf_grow_rate           & ! intent(in)
                                , l2n_stem                 & ! intent(in)
                                , c2n_stem                 & ! intent(in)
                                , c2n_storage              & ! intent(in)
                                , f_labile_leaf            ! ! intent(in)
      use phenology_coms , only : retained_carbon_fraction & ! intent(in)
                                , root_phen_factor         & ! intent(in)
                                , iphen_scheme             & ! intent(in)
                                , elongf_min               & ! intent(in)
                                , elongf_flush             ! ! intent(in)
      use consts_coms    , only : t3ple                    & ! intent(in)
                                , cice                     & ! intent(in)
                                , cliq                     & ! intent(in)
                                , alli                     ! ! intent(in)
      use ed_therm_lib   , only : calc_veg_hcap            & ! function
                                , update_veg_energy_cweh   ! ! subroutine
      use therm_lib      , only : tq2enthalpy              ! ! function
      use ed_max_dims    , only : n_pft                    ! ! intent(in)
      use ed_misc_coms   , only : current_time             ! ! intent(in)
      use allometry      , only : area_indices             & ! subroutine
                                , ed_balive                & ! function
                                , ed_biomass               & ! function
                                , size2bl                  ! ! function
      use phenology_aux  , only : daylength                ! ! function
      use plant_hydro    , only : rwc2tw                   & ! sub-routine
                                , twi2twe                  ! ! sub-routine
      use stable_cohorts , only : is_resolvable            ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)        , target     :: cpoly
      real                     , intent(in) :: lat
      integer                  , intent(in) :: doy  ! Day of year (1=Jan 1, 365/366=Dec 31)
      integer                  , intent(in) :: isi
      logical                  , intent(in) :: veget_dyn_on
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)           , pointer    :: csite
      type(patchtype)          , pointer    :: cpatch
      integer                               :: ipa
      integer                               :: ico
      integer                               :: isoil_lev
      integer                               :: kroot
      integer                               :: ipft
      logical                               :: leaf_out_cold
      logical                               :: drop_cold
      real                                  :: daylight
      real                                  :: delta_bleaf
      real                                  :: bl_max
      real                                  :: bl_full
      real                                  :: delta_broot
      real                                  :: br_max
      real                                  :: old_leaf_hcap
      real                                  :: old_wood_hcap
      real                                  :: old_leaf_water_im2
      real                                  :: old_wood_water_im2
      real                                  :: veg_energy_im2_in
      real                                  :: veg_water_im2_in
      real                                  :: elongf_try
      real                                  :: elongf_grow
      real                                  :: bleaf_in
      real                                  :: broot_in
      real                                  :: bstorage_in
      real                                  :: carbon_miss
      real                                  :: fast_n_drop
      real                                  :: struct_n_drop
      real                                  :: fgc_in_in
      real                                  :: stgc_in_in
      real                                  :: fsc_in_in
      real                                  :: stsc_in_in
      real                                  :: pat_bleaf_in
      real                                  :: pat_broot_in
      real                                  :: pat_bstorage_in
      real                                  :: pat_carbon_miss
      !----- Variables used only for debugging purposes. ----------------------------------!
      logical                  , parameter  :: printphen=.false.
      logical, dimension(n_pft), save       :: first_time=.true.
      !------------------------------------------------------------------------------------!


      !----- Level to evaluate the soil temperature. --------------------------------------!
      isoil_lev = nzg 
      !------------------------------------------------------------------------------------!

      !----- Calculate daylength for this gridcell. ---------------------------------------!
      daylight = daylength(lat, doy) 
      !------------------------------------------------------------------------------------!

      !----- Loop over patches. -----------------------------------------------------------!
      csite => cpoly%site(isi)
      !------------------------------------------------------------------------------------!

      patchloop: do ipa = 1,csite%npatches
         cpatch => csite%patch(ipa)


         !---------------------------------------------------------------------------------!
         !      Save patch-level litter inputs before growth balive.  We use these vari-   !
         ! ables to check carbon conservation at the patch level.                          !
         !---------------------------------------------------------------------------------!
         fgc_in_in       = csite%fgc_in (ipa)
         stgc_in_in      = csite%stgc_in(ipa)
         fsc_in_in       = csite%fsc_in (ipa)
         stsc_in_in      = csite%stsc_in(ipa)
         pat_bleaf_in    = 0.0
         pat_broot_in    = 0.0
         pat_bstorage_in = 0.0
         pat_carbon_miss = 0.0
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Litter inputs (fsc_in, stsc_in and alikes) used to be reset here.  This is  !
         ! not an intuitive place to reset them, plus carbon transfers could be lost       !
         ! because some carbon may go to the inputs after update_C_and_N_pools and this    !
         ! point (mostly because of cohort termination).  The input litter pools are now   !
         ! reset in update_C_and_N_pools, immediately after they are transferred to the    !
         ! litter pools.                                                                   !
         !---------------------------------------------------------------------------------!

         !----- Determine what phenology thresholds have been crossed. --------------------!
         call phenology_thresholds(daylight,csite%soil_tempk(isoil_lev,ipa)                &
                                  ,csite%sum_chd(ipa),csite%sum_dgd(ipa)                   &
                                  ,drop_cold,leaf_out_cold)

         cohortloop: do ico = 1,cpatch%ncohorts
            ipft    = cpatch%pft(ico)
            kroot   = cpatch%krdepth(ico)
            
            !----- Save input leaf and storage biomass for when dynamics is disabled. -----!
            bleaf_in          = cpatch%bleaf   (ico)
            broot_in          = cpatch%broot   (ico)
            bstorage_in       = cpatch%bstorage(ico)
            pat_bleaf_in      = pat_bleaf_in    + cpatch%nplant(ico) * cpatch%bleaf   (ico)
            pat_broot_in      = pat_broot_in    + cpatch%nplant(ico) * cpatch%broot   (ico)
            pat_bstorage_in   = pat_bstorage_in + cpatch%nplant(ico) * cpatch%bstorage(ico)
            veg_water_im2_in  = cpatch%leaf_water_im2(ico) + cpatch%wood_water_im2(ico)
            veg_energy_im2_in = cpatch%leaf_water_im2(ico)                                 &
                              * tq2enthalpy(cpatch%leaf_temp(ico),1.0,.true.)              &
                              + cpatch%wood_water_im2(ico)                                 &
                              * tq2enthalpy(cpatch%wood_temp(ico),1.0,.true.)
            !------------------------------------------------------------------------------!


            !----- Initially, we assume all leaves and roots stay. ------------------------!
            cpatch%leaf_drop(ico) = 0.0
            cpatch%root_drop(ico) = 0.0
            !------------------------------------------------------------------------------!

            !----- Find cohort-specific thresholds. ---------------------------------------!
            select case (iphen_scheme)
            case (1)
               !----- Get cohort-specific thresholds for prescribed phenology. ------------!
               call assign_prescribed_phen(cpoly%green_leaf_factor(ipft,isi)               &
                                          ,cpoly%leaf_aging_factor(ipft,isi)               &
                                          ,cpatch%dbh(ico),cpatch%hite(ico),ipft           &
                                          ,drop_cold,leaf_out_cold,bl_max)
            case default
               !----- Drop_cold is computed in phenology_thresholds for Botta scheme. -----!
               if (drop_cold) bl_max = 0.0
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Here we decide what to do depending on the phenology habit. There are    !
            ! five different types:                                                        !
            ! 0. Evergreen         - neither cold nor drought makes these plants to drop   !
            !                        their leaves;                                         !
            ! 1. Drought deciduous - these plants will drop all leaves when drought        !
            !                        conditions happen. By drought conditions we mean a    !
            !                        time when the available water drops below a threshold !
            ! 2. Cold deciduous    - these plants will drop their leaves when cold         !
            !                        conditions happen.                                    !
            ! 3. Light phenology   - these plants will control their leaf structure with   !
            !                        the light history (last 10 days).  They are also      !
            !                        drought-deciduous, similar to phenology 4;            !
            ! 4. Drought deciduous - similar to one, but the threshold is compared against !
            !                        a 10-day running average rather than the instant-     !
            !                        aneous value.                                         !
            !------------------------------------------------------------------------------!
            select case (phenology(ipft))
            case (0)
               !---------------------------------------------------------------------------!
               !    Evergreen, there is nothing to be done here except to assign maximum   !
               ! elongation factor.                                                        !
               !---------------------------------------------------------------------------!
               cpatch%elongf(ico)    = 1.0
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !     Drought deciduous.  Now we must check whether the plants still have   !
               ! enough water or it is too dry, or if there is some drought relief so      !
               ! leaves can start to grow again.                                           !
               !---------------------------------------------------------------------------!
               !----- This is the first guess for the new elongation factor. --------------!
               elongf_try  = max(0.0, min (1.0, cpatch%paw_avg(ico)))
               !---------------------------------------------------------------------------!



               if (elongf_try < 1.0 .and. cpatch%phenology_status(ico) /= -2) then


                  !------------------------------------------------------------------------!
                  !     Environment is drier, shed all leaves. Find leaf drop in terms of  !
                  ! carbon and nitrogen.  Nitrogen drop must account for different C:N     !
                  ! ratios between leaves and storage.                                     !
                  !------------------------------------------------------------------------!
                  delta_bleaf           = cpatch%bleaf(ico)
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  fast_n_drop           = f_labile_leaf(ipft) * delta_bleaf                &
                                        * ( 1./c2n_leaf(ipft)                              &
                                          - retained_carbon_fraction   / c2n_storage )
                  struct_n_drop         = (1.0 - f_labile_leaf(ipft)) * delta_bleaf        &
                                        * ( 1. / c2n_stem(ipft)                            &
                                          - retained_carbon_fraction   / c2n_storage )
                  !------------------------------------------------------------------------!



                  !----- Update storage only if vegetation dynamics is on. ----------------!
                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                              &
                                       + delta_bleaf * retained_carbon_fraction
                  !------------------------------------------------------------------------!



                  !----- Update plant carbon pools. ---------------------------------------!
                  cpatch%bleaf           (ico) = 0.0
                  cpatch%elongf          (ico) = 0.0
                  cpatch%phenology_status(ico) = -2
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Send the dropped leaves to soil carbon and nitrogen pools.         !
                  ! To conserve nitrogen, we assume all the lost nitrogen goes to the fast !
                  ! pool.                                                                  !
                  !------------------------------------------------------------------------!
                  csite%fgc_in(ipa)  = csite%fgc_in(ipa)                                   &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * f_labile_leaf(ipft)
                  csite%fgn_in(ipa)  = csite%fgn_in(ipa)                                   &
                                     + cpatch%nplant(ico) * fast_n_drop
                  csite%stgc_in(ipa) = csite%stgc_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))
                  csite%stgl_in(ipa) = csite%stgl_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))                         &
                                     * l2n_stem / c2n_stem(ipft)
                  csite%stgn_in(ipa) = csite%stgn_in(ipa)                                  &
                                     + cpatch%nplant(ico) * struct_n_drop
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico)  = cpatch%cb          (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_lightmax (13,ico)  = cpatch%cb_lightmax (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_moistmax (13,ico)  = cpatch%cb_moistmax (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_mlmax (13,ico)     = cpatch%cb_mlmax (13,ico)                  &
                                               - cpatch%leaf_drop (ico)
                  !------------------------------------------------------------------------!

               elseif(elongf_try > 1.0 .and. cpatch%phenology_status(ico) == -2) then
                  !------------------------------------------------------------------------!
                  !      It is time to flush.  Change phenology_status will update carbon  !
                  ! pools in growth_balive.                                                !
                  !------------------------------------------------------------------------!
                  cpatch%phenology_status(ico) = 1
                  cpatch%elongf          (ico) = 1.0
                  !------------------------------------------------------------------------!
               end if  ! critical moisture

            case (2)
               !---------------------------------------------------------------------------!
               !    Cold deciduous.  Here we must check two possibilities:                 !
               !                                                                           !
               ! 1. It is cold, and the plants have leaves, so we flag them with           !
               !    phenology_status=0 (leaves not growing) and the plants will start      !
               !    losing their leaves;                                                   !
               ! 2. The plant has no leaves, but the temperature and light conditions are  !
               !    okay again, and leaves can start growing.                              !
               !---------------------------------------------------------------------------!
               if (cpatch%phenology_status(ico) /= -2 .and. drop_cold) then
                  if (cpoly%green_leaf_factor(ipft,isi) < elongf_min) then
                      bl_max = 0.0
                  end if
                  delta_bleaf = cpatch%bleaf(ico) - bl_max

                  if (delta_bleaf > 0.0) then
                     !----- Set flag to -1 (dropping leaves). -----------------------------!
                     cpatch%phenology_status(ico) = -1
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Find leaf drop in terms of carbon and nitrogen.  Nitrogen drop  !
                     ! must account for different C:N ratios between leaves and storage.   !
                     !---------------------------------------------------------------------!
                     cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                     fast_n_drop           = f_labile_leaf(ipft) * delta_bleaf             &
                                           * ( 1./c2n_leaf(ipft)                           &
                                             - retained_carbon_fraction   / c2n_storage )
                     struct_n_drop         = (1.0 - f_labile_leaf(ipft)) * delta_bleaf     &
                                           * ( 1. / c2n_stem(ipft)                         &
                                             - retained_carbon_fraction   / c2n_storage )
                     !---------------------------------------------------------------------!



                     !----- Adjust plant carbon pools. ------------------------------------!
                     cpatch%bleaf(ico)    = cpatch%bleaf(ico) - delta_bleaf
                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                          + retained_carbon_fraction * delta_bleaf
                     !------------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Send the dropped leaves to soil carbon and nitrogen pools.      !
                     ! To conserve nitrogen, we assume all the lost nitrogen goes to the   !
                     ! fast pool.                                                          !
                     !---------------------------------------------------------------------!
                     csite%fgc_in(ipa)  = csite%fgc_in(ipa)                                &
                                        + cpatch%nplant(ico) * cpatch%leaf_drop(ico)       &
                                        * f_labile_leaf(ipft)
                     csite%fgn_in(ipa)  = csite%fgn_in(ipa)                                &
                                        + cpatch%nplant(ico) * fast_n_drop
                     csite%stgc_in(ipa) = csite%stgc_in(ipa)                               &
                                        + cpatch%nplant(ico) * cpatch%leaf_drop(ico)       &
                                        * (1.0 - f_labile_leaf(ipft))
                     csite%stgl_in(ipa) = csite%stgl_in(ipa)                               &
                                        + cpatch%nplant(ico) * cpatch%leaf_drop(ico)       &
                                        * (1.0 - f_labile_leaf(ipft))                      &
                                        * l2n_stem / c2n_stem(ipft)
                     csite%stgn_in(ipa) = csite%stgn_in(ipa)                               &
                                        + cpatch%nplant(ico) * struct_n_drop
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Deduct the leaf drop from the carbon balance.                  !
                     !---------------------------------------------------------------------!
                     cpatch%cb          (13,ico) = cpatch%cb         (13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax(13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax(13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax (13,ico) &
                                                 - cpatch%leaf_drop (ico)


                     !---------------------------------------------------------------------!
                  end if

                  !----- Set status flag. -------------------------------------------------!
                  if (bl_max == 0.0) then
                     cpatch%phenology_status(ico) = -2
                     cpatch%elongf(ico) = 0.
                  else
                     cpatch%elongf(ico) = 1.0 ! It should become green_leaf_factor...
                  end if
                  !------------------------------------------------------------------------!

               elseif (.not. drop_cold .and. cpatch%phenology_status(ico) == -2            &
                       .and. leaf_out_cold) then
                  !------------------------------------------------------------------------!
                  !      Update the phenology status (1 means that leaves are growing),    !
                  !------------------------------------------------------------------------!
                  cpatch%phenology_status(ico) = 1
                  cpatch%elongf          (ico) = 1.0 
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


            case (3,4) 
               !---------------------------------------------------------------------------!
               !    Drought deciduous.  Here we must check two possibilities:              !
               !                                                                           !
               ! 1. The soil has been dry recently, and the plants have leaves, so we flag !
               !    them with phenology_status=0 (leaves not growing) and the plants will  !
               !    start losing their leaves;                                             !
               ! 2. The plant has no leaves, but the soil has started to come back to more !
               !    moist conditions.  Given this situation, leaves can start growing      !
               !    again.                                                                 !
               !---------------------------------------------------------------------------!
               !----- This is the first guess for the new elongation factor. --------------!
               elongf_try  = max(0.0, min (1.0, cpatch%paw_avg(ico)))
               !----- If extremely dry, force the cohort to shed all leaves... ------------!
               if (elongf_try < elongf_min) elongf_try = 0.0
               !---------------------------------------------------------------------------!



               !----- Find the maximum allowed leaf biomass. ------------------------------!
               bl_max = elongf_try * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Delta_bleaf is the difference between the current leaf biomass and    !
               ! the maximum permitted given the soil moisture conditions.  If delta_bleaf !
               ! is positive, it means that the plant has more leaves than it should.      !
               !---------------------------------------------------------------------------!
               delta_bleaf = cpatch%bleaf(ico) - bl_max
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Check whether drought is becoming more or less severe.                !
               !---------------------------------------------------------------------------!
               if (delta_bleaf > 0.0) then
                  !------------------------------------------------------------------------!
                  !    Drought conditions are becoming more severe, drop leaves.           !
                  !------------------------------------------------------------------------!
                  if (elongf_try >= elongf_min) then
                     cpatch%phenology_status(ico) = -1
                  else
                     cpatch%phenology_status(ico) = -2
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find leaf drop in terms of carbon and nitrogen.  Nitrogen drop     !
                  ! must account for different C:N ratios between leaves and storage.      !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  fast_n_drop           = f_labile_leaf(ipft) * delta_bleaf                &
                                        * ( 1./c2n_leaf(ipft)                              &
                                          - retained_carbon_fraction   / c2n_storage )
                  struct_n_drop         = (1.0 - f_labile_leaf(ipft)) * delta_bleaf        &
                                        * ( 1. / c2n_stem(ipft)                            &
                                          - retained_carbon_fraction   / c2n_storage )
                  !------------------------------------------------------------------------!




                  !----- Adjust plant carbon pools and elongation factor. -----------------!
                  cpatch%bleaf     (ico) = cpatch%bleaf(ico) - delta_bleaf
                  cpatch%bstorage  (ico) = cpatch%bstorage(ico)                            &
                                         + retained_carbon_fraction * delta_bleaf
                  cpatch%elongf    (ico) = elongf_try
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Send the dropped leaves to soil carbon and nitrogen pools.         !
                  ! To conserve nitrogen, we assume all the lost nitrogen goes to the      !
                  ! fast pool.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fgc_in(ipa)  = csite%fgc_in(ipa)                                   &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * f_labile_leaf(ipft)
                  csite%fgn_in(ipa)  = csite%fgn_in(ipa)                                   &
                                     + cpatch%nplant(ico) * fast_n_drop
                  csite%stgc_in(ipa) = csite%stgc_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))
                  csite%stgl_in(ipa) = csite%stgl_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))                         &
                                     * l2n_stem / c2n_stem(ipft)
                  csite%stgn_in(ipa) = csite%stgn_in(ipa)                                  &
                                     + cpatch%nplant(ico) * struct_n_drop
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico) = cpatch%cb          (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax    (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  !------------------------------------------------------------------------!

               elseif (cpatch%phenology_status(ico) /= 0) then
                  !------------------------------------------------------------------------!
                  !       Elongation factor could increase, but we first check whether it  !
                  ! is safe to do so based on the phenology status.                        !
                  !------------------------------------------------------------------------!
                  select case(cpatch%phenology_status(ico))
                  case (1)
                     !----- Leaves were already growing, keep growing. --------------------!
                     cpatch%elongf          (ico) = elongf_try
                     !---------------------------------------------------------------------!
                  case (-1,-2)
                     !---------------------------------------------------------------------!
                     !     Leaves were dropping or gone, we first check that conditions    !
                     ! are really improving before we turn on leaf production.             !
                     !---------------------------------------------------------------------!
                     elongf_grow = min(1.0,max(elongf_flush,cpatch%elongf(ico)+0.02))
                     if (elongf_try >= elongf_grow) then
                        cpatch%elongf          (ico) = elongf_try
                        cpatch%phenology_status(ico) = 1
                     end if
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            case (5) 
               !---------------------------------------------------------------------------!
               !    Drought deciduous driven by plant hydrodynamics.  We track the number  !
               ! of consecutive wet days and dry days.  We then modify the phenology       !
               ! status and elongf whenever these numbers cross the wet- or dry-day        !
               ! thresholds.                                                               !
               !---------------------------------------------------------------------------!


               !----- Update consecutive dry days. ----------------------------------------!
               if (cpatch%dmax_leaf_psi(ico) < leaf_psi_tlp(ipft)) then
                  !---- Another dry day. --------------------------------------------------!
                  cpatch%low_leaf_psi_days(ico) = cpatch%low_leaf_psi_days(ico) + 1
                  !------------------------------------------------------------------------!
               else
                  !---- A relatively wet day.  Reset the number of dry days. --------------!
                  cpatch%low_leaf_psi_days(ico) = 0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !----- Update consecutive wet days. ----------------------------------------!
               if (cpatch%dmax_leaf_psi(ico) >= 0.5 * leaf_psi_tlp(ipft)) then
                  !---- Another wet day. --------------------------------------------------!
                  cpatch%high_leaf_psi_days(ico) = cpatch%high_leaf_psi_days(ico) + 1
                  !------------------------------------------------------------------------!
               else
                  !---- A relatively dry day.  Reset the number of wet days. --------------!
                  cpatch%high_leaf_psi_days(ico) = 0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !----- Modify elongf and phenology_status whenever necessary. --------------!
               if (cpatch%low_leaf_psi_days(ico) >= low_psi_threshold(ipft)) then
                  !----- Too many dry days, decrease elongation factor. -------------------!
                  elongf_try = max(0., cpatch%elongf(ico) - leaf_shed_rate(ipft))
                  !------------------------------------------------------------------------!
               else if (cpatch%high_leaf_psi_days(ico) >= high_psi_threshold(ipft)) then
                  !----- Good sequence of wet days, increase elongation factor. -----------!
                  select case (cpatch%phenology_status(ico))
                  case (-2)
                     !---------------------------------------------------------------------!
                     !    Currently without any leaves.  Let the cohort grow a small       !
                     ! fraction of leaves to check whether the condition has indeed become !
                     ! better.                                                             !
                     !---------------------------------------------------------------------!
                     elongf_try = elongf_min + 0.01
                     !---------------------------------------------------------------------!
                  case default
                     !----- Apply typical increase factor. --------------------------------!
                     elongf_try = min(1.0, cpatch%elongf(ico) + leaf_grow_rate(ipft))
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!
               else
                  !---- Conditions did not change much, keep the same elongation factor. --!
                  elongf_try = cpatch%elongf(ico)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !----- If extremely dry, force the cohort to shed all leaves... ------------!
               if (elongf_try < elongf_min) elongf_try = 0.0
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Find the maximum allowed leaf/root biomass.  The senescing rate for    !
               ! fine roots, relative to leaf shedding, is controlled by root_phen_factor. !
               ! It is normally a good idea to set this parameter to be greater than 1, to !
               ! allow for some fine roots to persist even if all leaves have been shed,   !
               ! otherwise plants will be unable to extract water once the rains return.   !
               !---------------------------------------------------------------------------!
               bl_full = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
               bl_max = elongf_try * bl_full
               if (root_phen_factor > 0.) then
                  br_max = bl_full * q(ipft)                                               &
                         * (elongf_try + root_phen_factor - 1.) / root_phen_factor
               else
                  br_max = bl_full * q(ipft)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Delta_bleaf is the difference between the current leaf biomass and    !
               ! the maximum permitted given the soil moisture conditions.  If delta_bleaf !
               ! is positive, it means that the plant has more leaves than it should.      !
               !---------------------------------------------------------------------------!
               delta_bleaf = cpatch%bleaf(ico) - bl_max
               delta_broot = cpatch%broot(ico) - br_max
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Check whether drought is becoming more or less severe.                !
               !---------------------------------------------------------------------------!
               if (delta_bleaf > 0.0) then
                  !------------------------------------------------------------------------!
                  !    Drought conditions are becoming more severe, drop leaves.           !
                  !------------------------------------------------------------------------!
                  if (elongf_try >= elongf_min) then
                     cpatch%phenology_status(ico) = -1
                  else
                     cpatch%phenology_status(ico) = -2
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find leaf drop in terms of carbon and nitrogen.  Nitrogen drop     !
                  ! must account for different C:N ratios between leaves and storage.      !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  fast_n_drop           = f_labile_leaf(ipft) * delta_bleaf                &
                                        * ( 1./c2n_leaf(ipft)                              &
                                          - retained_carbon_fraction   / c2n_storage )
                  struct_n_drop         = (1.0 - f_labile_leaf(ipft)) * delta_bleaf        &
                                        * ( 1. / c2n_stem(ipft)                            &
                                          - retained_carbon_fraction   / c2n_storage )
                  !------------------------------------------------------------------------!




                  !----- Adjust plant carbon pools and elongation factor. -----------------!
                  cpatch%bleaf     (ico) = cpatch%bleaf(ico) - delta_bleaf
                  cpatch%bstorage  (ico) = cpatch%bstorage(ico)                            &
                                         + retained_carbon_fraction * delta_bleaf
                  cpatch%elongf    (ico) = elongf_try
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Send the dropped leaves to soil carbon and nitrogen pools.         !
                  ! To conserve nitrogen, we assume all the lost nitrogen goes to the      !
                  ! fast pool.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fgc_in(ipa)  = csite%fgc_in(ipa)                                   &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * f_labile_leaf(ipft)
                  csite%fgn_in(ipa)  = csite%fgn_in(ipa)                                   &
                                     + cpatch%nplant(ico) * fast_n_drop
                  csite%stgc_in(ipa) = csite%stgc_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))
                  csite%stgl_in(ipa) = csite%stgl_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%leaf_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))                         &
                                     * l2n_stem / c2n_stem(ipft)
                  csite%stgn_in(ipa) = csite%stgn_in(ipa)                                  &
                                     + cpatch%nplant(ico) * struct_n_drop
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico) = cpatch%cb          (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax    (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  !------------------------------------------------------------------------!
               elseif (cpatch%phenology_status(ico) /= 0) then
                  !------------------------------------------------------------------------!
                  !       Elongation factor could increase. No need to check the safety    !
                  ! for flushing leaves, because this is done during the calculation of    !
                  ! high_leaf_psi_days.                                                    !
                  !------------------------------------------------------------------------!
                  cpatch%elongf          (ico) = elongf_try
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Only modify phenology_status when delta_bleaf < 0. This avoids      !
                  ! changing phenology_status = 1 while elongf = 0.                        !
                  !------------------------------------------------------------------------!
                  if (delta_bleaf < 0.) cpatch%phenology_status(ico) = 1
                  !------------------------------------------------------------------------!
               end if
               !------------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Check whether drought is becoming more or less severe.                !
               !---------------------------------------------------------------------------!
               if (delta_broot > 0.0) then


                  !------------------------------------------------------------------------!
                  !     Find leaf drop in terms of carbon and nitrogen.  Nitrogen drop     !
                  ! must account for different C:N ratios between leaves and storage.      !
                  !------------------------------------------------------------------------!
                  cpatch%root_drop(ico) = (1.0 - retained_carbon_fraction) * delta_broot
                  fast_n_drop           = f_labile_leaf(ipft) * delta_broot                &
                                        * ( 1./c2n_leaf(ipft)                              &
                                          - retained_carbon_fraction   / c2n_storage )
                  struct_n_drop         = (1.0 - f_labile_leaf(ipft)) * delta_broot        &
                                        * ( 1. / c2n_stem(ipft)                            &
                                          - retained_carbon_fraction   / c2n_storage )
                  !------------------------------------------------------------------------!




                  !----- Adjust plant carbon pools and elongation factor. -----------------!
                  cpatch%bleaf     (ico) = cpatch%bleaf(ico) - delta_broot
                  cpatch%bstorage  (ico) = cpatch%bstorage(ico)                            &
                                         + retained_carbon_fraction * delta_broot
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Send the dropped leaves to soil carbon and nitrogen pools.         !
                  ! To conserve nitrogen, we assume all the lost nitrogen goes to the      !
                  ! fast pool.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fsc_in(ipa)  = csite%fsc_in(ipa)                                   &
                                     + cpatch%nplant(ico) * cpatch%root_drop(ico)          &
                                     * f_labile_leaf(ipft)
                  csite%fsn_in(ipa)  = csite%fsn_in(ipa)                                   &
                                     + cpatch%nplant(ico) * fast_n_drop
                  csite%stsc_in(ipa) = csite%stsc_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%root_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))
                  csite%stsl_in(ipa) = csite%stsl_in(ipa)                                  &
                                     + cpatch%nplant(ico) * cpatch%root_drop(ico)          &
                                     * (1.0 - f_labile_leaf(ipft))                         &
                                     * l2n_stem / c2n_stem(ipft)
                  csite%stsn_in(ipa) = csite%stsn_in(ipa)                                  &
                                     + cpatch%nplant(ico) * struct_n_drop
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico) = cpatch%cb          (13,ico)                &
                                              - cpatch%root_drop      (ico)
                  cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax (13,ico)                &
                                              - cpatch%root_drop      (ico)
                  cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax (13,ico)                &
                                              - cpatch%root_drop      (ico)
                  cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax    (13,ico)                &
                                              - cpatch%root_drop      (ico)
                  !------------------------------------------------------------------------!
               end if
               !------------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     In case vegetation dynamics is turned off, we replace leaf biomass with  !
            ! the equilibrium leaf biomass, whilst we keep the same storage pool as the    !
            ! input.                                                                       !
            !------------------------------------------------------------------------------!
            if (.not. veget_dyn_on) then
               elongf_try           = cpoly%green_leaf_factor(ipft,isi) * cpatch%elongf(ico)
               cpatch%bleaf(ico)    = elongf_try                                           &
                                    * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
               cpatch%bstorage(ico) = bstorage_in
               !----- Fix phenology status according to elongation factor. ----------------!
               if (elongf_try == 1.0) then
                  cpatch%phenology_status(ico) = 0
               elseif (elongf_try < elongf_min) then
                  cpatch%phenology_status(ico) = -2
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Adjust root biomass in case phenology is 5 (drought-deciduous driven !
               ! by hydrodynamics).                                                        !
               !---------------------------------------------------------------------------!
               if (phenology(ipft) == 5 .and. root_phen_factor > 0.) then
                  cpatch%broot(ico) = q(ipft) * (elongf_try + root_phen_factor - 1.)       &
                                    * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)       &
                                    / root_phen_factor
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- Update balive. ---------------------------------------------------------!
            cpatch%balive(ico) = ed_balive(cpatch,ico)
            !------------------------------------------------------------------------------!


            !----- Update LAI, WAI, and CAI accordingly. ----------------------------------!
            call area_indices(cpatch, ico)
            !------------------------------------------------------------------------------!




            !----- Update above-ground biomass. -------------------------------------------!
            cpatch%agb(ico) = ed_biomass(cpatch, ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    The leaf biomass of the cohort has changed, update the vegetation energy  !
            ! using a constant temperature assumption.                                     !
            !------------------------------------------------------------------------------!
            old_leaf_hcap       = cpatch%leaf_hcap(ico)
            old_wood_hcap       = cpatch%wood_hcap(ico)
            old_leaf_water_im2  = cpatch%leaf_water_im2(ico)
            old_wood_water_im2  = cpatch%wood_water_im2(ico)
            call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico),cpatch%bsapwooda(ico)  &
                              ,cpatch%bbarka(ico),cpatch%nplant(ico),cpatch%pft(ico)       &
                              ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
            call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                          &
                       ,cpatch%bleaf(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)      &
                       ,cpatch%bdeada(ico),cpatch%bdeadb(ico),cpatch%broot(ico)            &
                       ,cpatch%dbh(ico),cpatch%pft(ico),cpatch%leaf_water_int(ico)         &
                       ,cpatch%wood_water_int(ico))
            !----- Convert total water content (kgW/plant) to extensive (kgW/m2). ---------!
            call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)             &
                        ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)                     &
                        ,cpatch%wood_water_im2(ico))
            call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap          &
                                       ,old_leaf_water_im2,old_wood_water_im2)
            call is_resolvable(csite,ipa,ico)
            !------------------------------------------------------------------------------!


            !----- Print some debugging stuff if the code is set for it. ------------------!
            if (printphen) then
               ipft=cpatch%pft(ico)
               if (first_time(ipft)) then
                  first_time(ipft) = .false.
                  write (unit=40+ipft,fmt='(a10,7(1x,a12))')                               &
                      '      TIME','       PATCH','      COHORT','      NPLANT'            &
                                  ,'   LEAF_DROP','   ROOT_DROP','     PAW_AVG'            &
                                  ,'      ELONGF'
               end if
            
               write (unit=40+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),6(1x,es12.5))')         &
                    current_time%month,'/',current_time%date,'/',current_time%year,ipa,ico &
                   ,cpatch%nplant(ico),cpatch%leaf_drop(ico),cpatch%root_drop(ico)         &
                   ,cpatch%paw_avg(ico),cpatch%elongf(ico)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Before we move to the next cohort, we make sure that carbon is being     !
            ! conserved by this cohort.                                                    !
            !------------------------------------------------------------------------------!
            if (veget_dyn_on) then
               call check_bphen_cohort(csite,ipa,ico,bleaf_in,broot_in,bstorage_in         &
                                      ,cpoly%green_leaf_factor(ipft,isi),carbon_miss)
               pat_carbon_miss = pat_carbon_miss + cpatch%nplant(ico) * carbon_miss
            end if
            !------------------------------------------------------------------------------!

         end do cohortloop
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!


      end do patchloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine update_phenology
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine establishes whether it's time to drop leaves or start flushing    !
   ! them for cold deciduous or temperate drought deciduous.                               !
   ! MLO. Shouldn't we have a similar criterion for both tropical and temperate, based on  !
   ! a long term dry condition?                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine phenology_thresholds(daylight,soil_temp,sum_chd,sum_dgd,drop_cold            &
                                  ,leaf_out_cold)
      use phenology_coms, only : dl_tr        & ! intent(in)
                               , st_tr1       & ! intent(in)
                               , st_tr2       & ! intent(in)
                               , phen_a       & ! intent(in)
                               , phen_b       & ! intent(in)
                               , phen_c       & ! intent(in)
                               , iphen_scheme ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real                   , intent(in)    :: daylight      ! Daytime Length
      real                   , intent(in)    :: soil_temp     ! 
      real                   , intent(inout) :: sum_dgd       !
      real                   , intent(inout) :: sum_chd       !
      logical                , intent(out)   :: drop_cold     !
      logical                , intent(out)   :: leaf_out_cold !
      !----- Local variables --------------------------------------------------------------!
      real                                   :: gdd_threshold
      !------------------------------------------------------------------------------------!

      !----- Initialize variables. --------------------------------------------------------!
      drop_cold     = .false.
      leaf_out_cold = .false.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check which phenology scheme to use.                                           !
      !------------------------------------------------------------------------------------!
      select case (iphen_scheme)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Prescribed phenology, skip this and fill with the prescribed phenology.     !
         !---------------------------------------------------------------------------------!
         drop_cold     = .false.
         leaf_out_cold = .false.
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !    Predicted phenology, check if this is the time to drop leaves or to start    !
         ! flushing.                                                                       !
         !---------------------------------------------------------------------------------!

         !----- Too cold or too dark, time to shed leaves... ------------------------------!
         drop_cold = (daylight <= dl_tr .and. soil_temp < st_tr1) .or.  soil_temp < st_tr2
         !---------------------------------------------------------------------------------!


         !----- Warm again, time to flush leaves. -----------------------------------------!
         gdd_threshold = phen_a + phen_b * exp(phen_c * sum_chd)
         leaf_out_cold = sum_dgd >= gdd_threshold
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine phenology_thresholds
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine assigns the prescribed phenology.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine assign_prescribed_phen(green_leaf_factor,leaf_aging_factor,dbh,height,pft    &
                                    ,drop_cold,leaf_out_cold,bl_max)
      use allometry     , only : size2bl
      use phenology_coms, only : elongf_min
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      logical, intent(out) :: drop_cold
      logical, intent(out) :: leaf_out_cold
      real   , intent(out) :: bl_max
      real   , intent(in)  :: green_leaf_factor
      real   , intent(in)  :: leaf_aging_factor
      real   , intent(in)  :: dbh
      real   , intent(in)  :: height
      integer, intent(in)  :: pft
      !------------------------------------------------------------------------------------!

      drop_cold     = green_leaf_factor /= leaf_aging_factor
      leaf_out_cold = green_leaf_factor > elongf_min .and. (.not. drop_cold)
      bl_max        = green_leaf_factor * size2bl(dbh, height, pft)

      return
   end subroutine assign_prescribed_phen
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the cohort level.  Minor      !
   ! truncation errors ma cause slightly negative pools.  In this case, we fix them and    !
   ! account for the missed carbon source.  Otherwise, we stop any cohort that is          !
   ! attempting to smuggle or to evade carbon.                                             !
   !---------------------------------------------------------------------------------------!
   subroutine check_bphen_cohort(csite,ipa,ico,bleaf_in,broot_in,bstorage_in               &
                                ,green_leaf_factor,carbon_miss)
      use ed_state_vars  , only : sitetype                 & ! structure
                                , patchtype                ! ! structure
      use allometry      , only : size2bl                  ! ! function
      use budget_utils   , only : tol_carbon_budget        ! ! intent(in)
      use pft_coms       , only : min_dbh                  & ! intent(in)
                                , hgt_min                  & ! intent(in)
                                , q                        & ! intent(in)
                                , phenology                ! ! intent(in)
      use phenology_coms , only : retained_carbon_fraction & ! intent(in)
                                , iphen_scheme             ! ! intent(in)
      use ed_misc_coms   , only : current_time             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target      :: csite
      integer         , intent(in)  :: ipa
      integer         , intent(in)  :: ico
      real            , intent(in)  :: bleaf_in
      real            , intent(in)  :: broot_in
      real            , intent(in)  :: bstorage_in
      real            , intent(in)  :: green_leaf_factor
      real            , intent(out) :: carbon_miss
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer     :: cpatch
      integer                       :: ipft
      real                          :: bleaf_ok_min
      real                          :: broot_ok_min
      real                          :: bstorage_ok_min
      real                          :: btotal_in
      real                          :: btotal_fn
      real                          :: delta_btotal
      real                          :: resid_btotal
      logical                       :: neg_biomass
      logical                       :: btotal_violation
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter  :: fmti='(a,1x,i14)'
      character(len=13), parameter  :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter  :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Handy aliases. ---------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      ipft   =  cpatch%pft (ico)
      !------------------------------------------------------------------------------------!


      !----- Find the minimum acceptable biomass. -----------------------------------------!
      bleaf_ok_min     = - tol_carbon_budget * size2bl(min_dbh(ipft),hgt_min(ipft),ipft)
      broot_ok_min     = q(ipft) * bleaf_ok_min
      bstorage_ok_min  = bleaf_ok_min + broot_ok_min
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check leaves and storage, to make sure they have sensible numbers.  Tiny        !
      ! negative stocks will be tolerated because of floating point truncation, but don't  !
      ! fix anything in case any of the pools is too negative.                             !
      !------------------------------------------------------------------------------------!
      neg_biomass    = cpatch%bleaf    (ico) < bleaf_ok_min     .or.                       &
                       cpatch%broot    (ico) < broot_ok_min     .or.                       &
                       cpatch%bstorage (ico) < bstorage_ok_min
      if (.not. neg_biomass) then
         !----- Account for any potential violation of carbon stocks. ---------------------!
         carbon_miss = - min(cpatch%bleaf   (ico),0.0)                                     &
                       - min(cpatch%broot   (ico),0.0)                                     &
                       - min(cpatch%bstorage(ico),0.0)
         !---------------------------------------------------------------------------------!


         !----- Make sure that all pools are non-negative. --------------------------------!
         cpatch%bleaf    (ico) = max(cpatch%bleaf    (ico),0.0)
         cpatch%broot    (ico) = max(cpatch%broot    (ico),0.0)
         cpatch%bstorage (ico) = max(cpatch%bstorage (ico),0.0)
         !---------------------------------------------------------------------------------!
      else
         !----- Set missing carbon to zero so the code works with debugging. --------------!
         carbon_miss = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after leaf phenology.                           !
      !------------------------------------------------------------------------------------!
      btotal_in         = bleaf_in          + broot_in          + bstorage_in
      btotal_fn         = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bstorage(ico)
      delta_btotal      = - cpatch%leaf_drop(ico) - cpatch%root_drop(ico)
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
         write(unit=*,fmt='(a)')  '|       !!!   Cohort Bphen budget failed   !!!       |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
         write(unit=*,fmt=fmti )  ' COHORT              : ',ico
         write(unit=*,fmt=fmti )  ' IPFT                : ',ipft
         write(unit=*,fmt=fmtf )  ' DBH                 : ',cpatch%dbh   (ico)
         write(unit=*,fmt=fmtf )  ' HITE                : ',cpatch%hite  (ico)
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmti )  ' IPHEN_SCHEME        : ',iphen_scheme
         write(unit=*,fmt=fmti )  ' PHENOLOGY_STRATEGY  : ',phenology(ipft)
         write(unit=*,fmt=fmti )  ' PHENOLOGY_STATUS    : ',cpatch%phenology_status(ico)
         write(unit=*,fmt=fmtf )  ' ELONGF              : ',cpatch%elongf(ico)
         write(unit=*,fmt=fmtf )  ' GREEN_LEAF_FACTOR   : ',green_leaf_factor
         write(unit=*,fmt=fmtf )  ' RETAINED_C_FRACTION : ',retained_carbon_fraction
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BLEAF_IN            : ',bleaf_in
         write(unit=*,fmt=fmtf )  ' BROOT_IN            : ',broot_in
         write(unit=*,fmt=fmtf )  ' BSTORAGE_IN         : ',bstorage_in
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BLEAF_FN            : ',cpatch%bleaf    (ico)
         write(unit=*,fmt=fmtf )  ' BROOT_FN            : ',cpatch%broot    (ico)
         write(unit=*,fmt=fmtf )  ' BSTORAGE_FN         : ',cpatch%bstorage (ico)
         write(unit=*,fmt=fmtf )  ' LEAF_DROP_FN        : ',cpatch%leaf_drop(ico)
         write(unit=*,fmt=fmtf )  ' ROOT_DROP_FN        : ',cpatch%root_drop(ico)
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
                         ,'check_bphen_cohort','phenology_driv.f90')
      end if
      !------------------------------------------------------------------------------------!
      return
   end subroutine check_bphen_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the patch level.  Minor       !
   ! truncation errors may cause slightly negative pools, which are accounted.  Otherwise, !
   ! we stop any patch that is attempting to smuggle or to evade carbon.                   !
   !---------------------------------------------------------------------------------------!
   subroutine check_bphen_patch(csite,ipa,fgc_in_in,stgc_in_in,fsc_in_in,stsc_in_in        &
                               ,pat_bleaf_in,pat_broot_in,pat_bstorage_in,pat_carbon_miss)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use budget_utils , only : tol_carbon_budget  ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target        :: csite
      integer         , intent(in)    :: ipa
      real            , intent(in)    :: fgc_in_in
      real            , intent(in)    :: stgc_in_in
      real            , intent(in)    :: fsc_in_in
      real            , intent(in)    :: stsc_in_in
      real            , intent(in)    :: pat_bleaf_in
      real            , intent(in)    :: pat_broot_in
      real            , intent(in)    :: pat_bstorage_in
      real            , intent(in)    :: pat_carbon_miss
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer       :: cpatch
      integer                         :: ico
      integer                         :: ipft
      real                            :: fgc_in_fn
      real                            :: stgc_in_fn
      real                            :: fsc_in_fn
      real                            :: stsc_in_fn
      real                            :: pat_bleaf_fn
      real                            :: pat_broot_fn
      real                            :: pat_bstorage_fn
      real                            :: pat_leaf_drop
      real                            :: pat_root_drop
      real                            :: pat_biomass_in
      real                            :: pat_biomass_fn
      real                            :: pat_btotal_in
      real                            :: pat_btotal_fn
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
      cpatch     => csite%patch  (ipa)
      fgc_in_fn  =  csite%fgc_in (ipa)
      stgc_in_fn =  csite%stgc_in(ipa)
      fsc_in_fn  =  csite%fsc_in (ipa)
      stsc_in_fn =  csite%stsc_in(ipa)
      !------------------------------------------------------------------------------------!


      !----- Count current stocks. --------------------------------------------------------!
      pat_bleaf_fn    = 0.0
      pat_broot_fn    = 0.0
      pat_bstorage_fn = 0.0
      pat_leaf_drop   = 0.0
      pat_root_drop   = 0.0
      do ico=1,cpatch%ncohorts
         ipft            = cpatch%pft(ico)
         pat_bleaf_fn    = pat_bleaf_fn    + cpatch%nplant(ico) * cpatch%bleaf    (ico)
         pat_broot_fn    = pat_broot_fn    + cpatch%nplant(ico) * cpatch%broot    (ico)
         pat_bstorage_fn = pat_bstorage_fn + cpatch%nplant(ico) * cpatch%bstorage (ico)
         pat_leaf_drop   = pat_leaf_drop   + cpatch%nplant(ico) * cpatch%leaf_drop(ico)
         pat_root_drop   = pat_root_drop   + cpatch%nplant(ico) * cpatch%root_drop(ico)
      end do
      !------------------------------------------------------------------------------------!


      !------ Summary of the carbon stocks. -----------------------------------------------!
      pat_biomass_in = pat_bleaf_in    + pat_broot_in    + pat_bstorage_in
      pat_biomass_fn = pat_bleaf_fn    + pat_broot_fn    + pat_bstorage_fn
      soilc_in_in    = fgc_in_in       + stgc_in_in      + fsc_in_in       + stsc_in_in
      soilc_in_fn    = fgc_in_fn       + stgc_in_fn      + fsc_in_fn       + stsc_in_fn
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after active tissue growth.  Yield is positive  !
      ! in the residual calculation because it is a committed permanent loss of carbon     !
      ! (shipped off-site).                                                                !
      !------------------------------------------------------------------------------------!
      pat_btotal_in        = pat_biomass_in + soilc_in_in
      pat_btotal_fn        = pat_biomass_fn + soilc_in_fn
      resid_pat_btotal     = pat_btotal_fn  - pat_btotal_in - pat_carbon_miss
      pat_btotal_violation = abs(resid_pat_btotal) > (tol_carbon_budget * pat_btotal_in)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( pat_btotal_violation ) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|        !!!   Patch Bphen budget failed   !!!       |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BLEAF_IN            : ',pat_bleaf_in
         write(unit=*,fmt=fmtf )  ' BROOT_IN            : ',pat_broot_in
         write(unit=*,fmt=fmtf )  ' BSTORAGE_IN         : ',pat_bstorage_in
         write(unit=*,fmt=fmtf )  ' FGC_IN_IN           : ',fgc_in_in
         write(unit=*,fmt=fmtf )  ' FSC_IN_IN           : ',fsc_in_in
         write(unit=*,fmt=fmtf )  ' STGC_IN_IN          : ',stgc_in_in
         write(unit=*,fmt=fmtf )  ' STSC_IN_IN          : ',stsc_in_in
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BLEAF_FN            : ',pat_bleaf_fn
         write(unit=*,fmt=fmtf )  ' BROOT_FN            : ',pat_broot_fn
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
         write(unit=*,fmt=fmtf )  ' LEAF_DROP           : ',pat_leaf_drop
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BTOTAL_IN           : ',pat_btotal_in
         write(unit=*,fmt=fmtf )  ' BTOTAL_FN           : ',pat_btotal_fn
         write(unit=*,fmt=fmtf )  ' CARBON_MISS         : ',pat_carbon_miss
         write(unit=*,fmt=fmtf )  ' RESIDUAL_BTOTAL     : ',resid_pat_btotal
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  ' '


         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_bphen_patch','phenology_driv.f90')
      end if
      !------------------------------------------------------------------------------------!





      return
   end subroutine check_bphen_patch
   !=======================================================================================!
   !=======================================================================================!
end module phenology_driv
!==========================================================================================!
!==========================================================================================!
