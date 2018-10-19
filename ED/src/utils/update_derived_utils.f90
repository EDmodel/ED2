module update_derived_utils
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the update of derived properties.                      !
   !---------------------------------------------------------------------------------------!
   subroutine update_derived_props(cgrid)
      use ed_state_vars , only : edtype      & ! structure
                               , polygontype & ! structure
                               , sitetype    ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer :: cpoly
      type(sitetype)    , pointer :: csite
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      !------------------------------------------------------------------------------------!
      
      do ipy = 1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           csite => cpoly%site(isi)

           do ipa = 1,csite%npatches
              call update_patch_derived_props(csite,ipa,.false.)
           end do

           call update_site_derived_props(cpoly, 0, isi)
        end do

        call update_polygon_derived_props(cgrid)
      end do

      return
   end subroutine update_derived_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will assign values derived from the basic properties of a given   !
   ! cohort.                                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_derived_props(cpatch,ico,lsl,new_year)

      use ed_state_vars  , only : patchtype               ! ! structure
      use pft_coms       , only : is_grass                ! ! function
      use allometry      , only : bd2dbh                  & ! function
                                , dbh2h                   & ! function
                                , size2krdepth            & ! function
                                , bl2dbh                  & ! function
                                , bl2h                    & ! function
                                , size2bl                 & ! function
                                , size2bt                 & ! function
                                , size2xb                 & ! function
                                , ed_balive               & ! function
                                , ed_biomass              & ! function
                                , area_indices            ! ! subroutine
      use physiology_coms, only : trait_plasticity_scheme ! ! intent(in)
      use consts_coms    , only : pio4                    ! ! intent(in)
      use ed_misc_coms   , only : igrass                  & ! intent(in)
                                , current_time            ! ! intent(in)
      use detailed_coms  , only : dt_census               & ! intent(in)
                                , yr1st_census            & ! intent(in)
                                , mon1st_census           & ! intent(in)
                                , min_recruit_dbh         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target     :: cpatch
      integer        , intent(in) :: ico
      integer        , intent(in) :: lsl
      logical        , intent(in) :: new_year
      !----- Local variables --------------------------------------------------------------!
      real                        :: bleaf_max
      integer                     :: ipft
      integer                     :: elapsed_months
      logical                     :: census_time
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the number of elapsed months since the first census, and decide whether    !
      ! there was a census last month or not.  It is absolutely fine to be negative, al-   !
      ! though none of the cohorts will be flagged as measured until the first census.     !
      !                                                                                    !
      ! IMPORTANT: the flag will be updated only AFTER the census month, because the time  !
      !            step is at the beginning of the month, likely to be just before the     !
      !            census...                                                               !
      !------------------------------------------------------------------------------------!
      elapsed_months = (current_time%year-yr1st_census-1)*12                               &
                     + current_time%month + (12 - mon1st_census - 1)
      census_time    = elapsed_months >= 0 .and. mod(elapsed_months,dt_census) == 0
      !------------------------------------------------------------------------------------!

      ipft    = cpatch%pft(ico)

      !----- Get DBH and height -----------------------------------------------------------!
      if (is_grass(ipft) .and. igrass == 1) then
          !---- New grasses get dbh_effective and height from bleaf. ----------------------!
          cpatch%dbh(ico)  = bl2dbh(cpatch%bleaf(ico), ipft)
          cpatch%hite(ico) = bl2h  (cpatch%bleaf(ico), ipft)
      else
          !---- Trees and old grasses get dbh from bdead. ---------------------------------!
          cpatch%dbh(ico)  = bd2dbh(ipft, cpatch%bdeada(ico), cpatch%bdeadb(ico))
          cpatch%hite(ico) = dbh2h (ipft, cpatch%dbh   (ico))
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update the recruitment flag regarding DBH if needed.                           !
      !------------------------------------------------------------------------------------!
      if (cpatch%dbh(ico) >= min_recruit_dbh) then
         cpatch%recruit_dbh(ico) = min(2,cpatch%recruit_dbh(ico) + 1)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update the census status if this is the time to do so.                         !
      !------------------------------------------------------------------------------------!
      if ( cpatch%dbh(ico) >= min_recruit_dbh .and. census_time ) then
         cpatch%census_status(ico) = min(2,cpatch%census_status(ico) + 1)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Because DBH may have increased, the maximum leaf biomass may be different,     !
      ! which will put plants off allometry even if they were on-allometry before.  Here   !
      ! we check whether this is the case.                                                 !
      !------------------------------------------------------------------------------------!
      if ((.not. is_grass(ipft)) .or. igrass /= 1) then
         select case (cpatch%phenology_status(ico))
         case (0)
            bleaf_max = size2bl(cpatch%dbh(ico),cpatch%hite(ico),cpatch%pft(ico))
            if (cpatch%bleaf(ico) < bleaf_max) cpatch%phenology_status(ico) = 1
         end select
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update plastic traits (SLA, Vm0).  This must be done before calculating LAI.   !
      !------------------------------------------------------------------------------------!
      select case (trait_plasticity_scheme)
      case (-1,1) ! Update trait every year
         if (new_year) call update_cohort_plastic_trait(cpatch,ico)
      case (-2,2) ! Update trait every month
         call update_cohort_plastic_trait(cpatch,ico)
      end select
      !------------------------------------------------------------------------------------!



      !----- Update LAI, WAI, and CAI. ----------------------------------------------------!
      call area_indices(cpatch, ico)
      !------------------------------------------------------------------------------------!


      !----- Update derived properties (AGB, BA, timber stocks and bark thickness). -------!
      cpatch%balive (ico) = ed_balive (cpatch, ico)
      cpatch%basarea(ico) = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
      cpatch%agb    (ico) = ed_biomass(cpatch, ico)
      cpatch%btimber(ico) = size2bt(cpatch%dbh(ico),cpatch%hite(ico),cpatch%bdeada(ico)    &
                                   ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)               &
                                   ,cpatch%pft(ico))
      cpatch%thbark(ico)  = size2xb(cpatch%dbh(ico),cpatch%hite(ico),cpatch%bbarka(ico)    &
                                   ,cpatch%bbarkb(ico),cpatch%pft(ico))
      !------------------------------------------------------------------------------------!


      !----- Update rooting depth ---------------------------------------------------------!
      cpatch%krdepth(ico) = size2krdepth(cpatch%hite(ico),cpatch%dbh(ico),ipft,lsl)
      !if new root depth is smaller keep the old one
      return
   end subroutine update_cohort_derived_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE UPDATE_COHORT_PLASTIC_TRAIT 
   !< \brief This subroutine will assign values for plastic functional traits driven by
   !< local light environment and thus depending on vertical structure of the canopy.
   !< \warning This function should be called after update_derived_cohort_props
   !< \details Refs: \n
   !<      Lloyd J, et al. 2010. Optimisation of photosynthetic carbon gain and
   !< within-canopy gradients of associated foliar traits for Amazon forest
   !< trees. Biogesciences, 7(6):1833-1859. doi:10.5194/bg-7-1833-2010.\n
   !=======================================================================================!
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_plastic_trait(cpatch,ico)
      use ed_state_vars  , only : patchtype               ! ! structure
      use pft_coms       , only : SLA                     & ! intent(in)
                                , kplastic_vm0            & ! intent(in)
                                , kplastic_sla            & ! intent(in)
                                , kplastic_ll             & ! intent(in)
                                , eplastic_vm0            & ! intent(in)
                                , eplastic_sla            & ! intent(in)
                                , laimax_plastic          & ! intent(in)
                                , lma_slope               & ! intent(in)
                                , Vm0                     & ! intent(in)
                                , leaf_turnover_rate      & ! intent(in)
                                , is_tropical             & ! intent(in)
                                , is_grass                ! ! intent(in)
      use consts_coms    , only : lnexp_min               & ! intent(in)
                                , lnexp_max               ! ! intent(in)
      use allometry      , only : size2bl                 ! ! function
      use physiology_coms, only : trait_plasticity_scheme ! ! intent(in)
      use phenology_coms , only : llspan_inf              ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target     :: cpatch       ! Current patch
      integer        , intent(in) :: ico          ! Cohort index
      !----- Local variables --------------------------------------------------------------!
      integer                     :: ipft         ! Alias for current PFT
      integer                     :: jco          ! Cohort count
      real                        :: max_cum_lai  ! Potential cumulative LAI
      real                        :: bl_max       ! Maximum attainable leaf biomass
      real                        :: lnexp        ! FPE-safe exponential test
      real                        :: new_sla      ! Updated SLA
      real                        :: sla_scaler   ! Scaling factor for SLA
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !  for now, only update the trait for tropical trees. However, It can be applied to  !
      !  temperate forests as well because the parameters come from a meta-analysis        !
      !  by Lloyd et al. 2010                                                              !
      !------------------------------------------------------------------------------------!
      ipft    = cpatch%pft(ico)
      if (is_grass(ipft) .or. (.not. is_tropical(ipft))) return
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 1. Find the maximum cumulative lai above the current cohort using the current SLA. !
      !------------------------------------------------------------------------------------!
      max_cum_lai = 0.  ! Set cumulative LAI as zero, and update only when # cohorts > 1.
      if (ico > 1) then
         !----- Accumulate LAI from the top cohort to current cohort. ---------------------!
         do jco = 1,ico-1
            bl_max      = size2bl(cpatch%dbh(jco),cpatch%hite(jco),cpatch%pft(jco))
            max_cum_lai = max_cum_lai + bl_max * cpatch%sla(jco) * cpatch%nplant(jco)
         end do
         !---------------------------------------------------------------------------------!
         max_cum_lai = min(laimax_plastic(ipft),max_cum_lai)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 2.  Update Vm0.  This should be defined at the top of canopy [sun-lit leaves].     !
      !     Note that the sign of kplastic_vm0 is typically negative, so this should       !
      !     reduce Vm0.                                                                    !
      !------------------------------------------------------------------------------------!
      lnexp              = max(lnexp_min,kplastic_vm0(ipft) * max_cum_lai)
      cpatch%vm_bar(ico) = Vm0(ipft) * exp(lnexp)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 3.  Update SLA.  Decide whether to use the bottom or top of the canopy as the      !
      !     reference.                                                                     !
      !------------------------------------------------------------------------------------!
      select case (trait_plasticity_scheme)
      case ( 1, 2)
         !------ SLA is defined at the top of canopy, use LAI to change SLA. --------------!
         lnexp   = max(lnexp_min,min(lnexp_max,kplastic_sla(ipft) * max_cum_lai))
         new_sla = SLA(ipft) * exp(lnexp)
         !---------------------------------------------------------------------------------!
      case (-1,-2)
         !------ SLA is defined at the bottom of canopy, use height to change SLA. --------!
         new_sla = SLA(ipft) / (1. + lma_slope(ipft) * cpatch%hite(ico))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 4.  Here we also need to retrospectively change leaf level state variables because !
      !     the leaf area has changed while we want to keep the flux the same. This is     !
      !     necessary for plant hydraulic calculations, which uses the water fluxes from   !
      !     'Last Timestep'.  For now we only update psi_open and psi_closed, which will   !
      !     be used in plant_hydro_driver. We will leave A_open and A_closed unchanged     !
      !     because growth of the day has already happen at this time point in the model.  !
      !------------------------------------------------------------------------------------!
      sla_scaler             = cpatch%sla(ico) / new_sla
      cpatch%sla       (ico) = new_sla
      cpatch%psi_open  (ico) = cpatch%psi_open  (ico) * sla_scaler
      cpatch%psi_closed(ico) = cpatch%psi_closed(ico) * sla_scaler
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 5.  Update leaf life span.                                                         !
      !------------------------------------------------------------------------------------!
      if (leaf_turnover_rate(ipft) > 0.0) then
         !---------------------------------------------------------------------------------!
         !     Decide which method to employ based on TRAIT_PLASTICITY_SCHEME.             !
         !---------------------------------------------------------------------------------!
         select case (trait_plasticity_scheme)
         case ( 1, 2)
            !------------------------------------------------------------------------------!
            !    Use a leaf longevity extinction/expansion factor derived from digitised   !
            ! data (Fig 5a of RK16).                                                       !
            !------------------------------------------------------------------------------!
            lnexp              = max(lnexp_min,min(lnexp_max,kplastic_ll(ipft)*max_cum_lai))
            cpatch%llspan(ico) = 12. / leaf_turnover_rate(ipft) * exp(lnexp)
            !------------------------------------------------------------------------------!
         case (-1,-2)
            !------------------------------------------------------------------------------!
            !    Follow Eqn. 1 of X17, assuming that the net assimilation is linearly      !
            ! correlated with Vcmax25, and find the ratio between canopy values and the    !
            ! plastic values.  For the term b we use X17 slope from Table S1               !
            ! (log10(b) ~ log10(Vcmax_m)).                                                 !
            !------------------------------------------------------------------------------!
            cpatch%llspan(ico) = 12. / leaf_turnover_rate(ipft)                            &
                               * ( cpatch%vm_bar(ico) / Vm0(ipft) ) ** eplastic_vm0(ipft)  &
                               * ( cpatch%sla   (ico) / SLA(ipft) ) ** eplastic_sla(ipft)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !---- Nothing lasts forever, so impose a maximum life span. ----------------------!
         cpatch%llspan(ico) = llspan_inf
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_cohort_plastic_trait
   !==========================================================================================!
   !==========================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will compute the growth and mortality rates.                       !
   !---------------------------------------------------------------------------------------!
   subroutine update_vital_rates(cpatch,ico,dbh_in,nplant_in,agb_in,ba_in,area,basal_area  &
                                 ,agb,basal_area_growth,agb_growth,basal_area_mort,agb_mort)

      use ed_state_vars , only : patchtype    ! ! structure
      use ed_max_dims   , only : n_pft        & ! intent(in)
                               , n_dbh        ! ! intent(in)
      use ed_misc_coms  , only : ddbhi        ! ! intent(in)
      use consts_coms   , only : pio4         ! ! intent(in)
      use pft_coms      , only : is_grass     ! ! function
      use allometry     , only : ed_biomass   ! ! function
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype)              , target        :: cpatch
      real                         , intent(in)    :: dbh_in
      real                         , intent(in)    :: nplant_in
      real                         , intent(in)    :: agb_in
      real                         , intent(in)    :: ba_in
      real                         , intent(in)    :: area
      integer                      , intent(in)    :: ico
      real, dimension(n_pft, n_dbh), intent(inout) :: basal_area
      real, dimension(n_pft, n_dbh), intent(inout) :: agb
      real, dimension(n_pft, n_dbh), intent(inout) :: basal_area_growth
      real, dimension(n_pft, n_dbh), intent(inout) :: agb_growth
      real, dimension(n_pft, n_dbh), intent(inout) :: basal_area_mort
      real, dimension(n_pft, n_dbh), intent(inout) :: agb_mort
      !----- Local variables --------------------------------------------------------------!
      integer                                      :: ipft
      integer                                      :: idbh

      !------------------------------------------------------------------------------------!


      !----- Make the alias for PFT type. -------------------------------------------------!
      ipft = cpatch%pft(ico)

      !----- Find the DBH bin. ------------------------------------------------------------!
      idbh = max(1,min(n_dbh,ceiling(dbh_in*ddbhi)))

      !----- Find the new basal area and above-ground biomass. ----------------------------!
      cpatch%basarea(ico)    = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
      cpatch%agb(ico)        = ed_biomass(cpatch, ico)

      !------------------------------------------------------------------------------------!
      !     Change the agb growth to kgC/plant/year, basal area to cm2/plant/year, and DBH !
      ! growth to cm/year.                                                                 !
      !------------------------------------------------------------------------------------!
      cpatch%dagb_dt   (ico) =    (cpatch%agb(ico)     - agb_in ) * 12.0
      cpatch%dlnagb_dt (ico) = log(cpatch%agb(ico)     / agb_in ) * 12.0
      cpatch%dba_dt    (ico) =    (cpatch%basarea(ico) - ba_in  ) * 12.0
      cpatch%dlnba_dt  (ico) = log(cpatch%basarea(ico) / ba_in  ) * 12.0
      cpatch%ddbh_dt   (ico) =    (cpatch%dbh(ico)     - dbh_in ) * 12.0
      cpatch%dlndbh_dt (ico) = log(cpatch%dbh(ico)     / dbh_in ) * 12.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     These are polygon-level variable, so they are done in kgC/m2.  Update the      !
      ! current basal area and above-ground biomass.                                       !
      !------------------------------------------------------------------------------------!
      if (is_grass(ipft)) return
      basal_area(ipft, idbh) = basal_area(ipft, idbh)                                      &
                             + area * cpatch%nplant(ico) * cpatch%basarea(ico)
      agb(ipft, idbh)        = agb(ipft, idbh)                                             &
                             + area * cpatch%nplant(ico) * cpatch%agb(ico)

      !------------------------------------------------------------------------------------!
      !    The growth and mortality census are applied only on those cohorts present on    !
      ! the first census.                                                                  !
      !------------------------------------------------------------------------------------!
      if (cpatch%first_census(ico) /= 1) return

      !------------------------------------------------------------------------------------!
      !   Computed for plants alive both at past census and current census.  These will be !
      ! given in cm2/m2/yr and kgC/m2/yr, respectively.                                    !
      !------------------------------------------------------------------------------------!
      basal_area_growth(ipft,idbh) = basal_area_growth(ipft,idbh)                          &
                                   + area * cpatch%nplant(ico) * pio4                      &
                                   * (cpatch%dbh(ico) * cpatch%dbh(ico) - dbh_in * dbh_in)
      agb_growth(ipft,idbh)        = agb_growth(ipft,idbh)                                 &
                                   + area * cpatch%nplant(ico)                             &
                                   * (cpatch%agb(ico) - agb_in)

      !------------------------------------------------------------------------------------!
      !    Computed for plants alive at past census but dead at current census.  These     !
      ! variables are also given in cm2/m2/yr and kgC/m2/yr, respectively.                 !
      !------------------------------------------------------------------------------------!
      basal_area_mort(ipft,idbh) = basal_area_mort(ipft,idbh)                              &
                                 + area * (nplant_in - cpatch%nplant(ico)) * ba_in

      !----- Calculation based on mort_litter includes TOTAL biomass, not AGB [[mcd]]. ----!
      agb_mort(ipft,idbh)        = agb_mort(ipft,idbh)                                     &
                                 + area * (nplant_in - cpatch%nplant(ico)) * agb_in

      return
   end subroutine update_vital_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine will take care of derived patch-level structural quantities.     !
   ! These depend on the results from reproduction, which in turn depends on structural    !
   ! growth results from all patches.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine update_patch_derived_props(csite,ipa,update_zcaneff)
     
      use ed_state_vars       , only : sitetype                   & ! structure
                                     , patchtype                  ! ! structure
      use allometry           , only : ed_biomass                 ! ! function
      use canopy_air_coms     , only : veg_height_min             & ! intent(in)
                                     , minimum_canopy_depth       & ! intent(in)
                                     , vh2vr                      & ! intent(in)
                                     , vh2dh                      ! ! intent(in)
      use soil_coms           , only : soil_rough                 & ! intent(in)
                                     , ny07_eq04_a                & ! intent(in)
                                     , ny07_eq04_m                & ! intent(in)
                                     , tiny_sfcwater_mass         ! ! intent(in)
      use consts_coms         , only : wdns                       & ! intent(in)
                                     , fsdns                      & ! intent(in)
                                     , fsdnsi                     & ! intent(in)
                                     , mmdryi                     & ! intent(in)
                                     , umol_2_kgC                 ! ! intent(in)
      use therm_lib           , only : tq2enthalpy                ! ! function
      use ed_misc_coms        , only : frqsumi                    ! ! intent(in)

      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      logical         , intent(in) :: update_zcaneff
      !----- Local variables --------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      real                         :: weight
      real                         :: weight_sum
      real                         :: total_sfcw_mass
      real                         :: bulk_sfcw_dens
      real                         :: fdelta_storage
      real                         :: old_can_depth
      real                         :: can_enthalpy
      integer                      :: ico
      integer                      :: k
      integer                      :: ksn
      integer                      :: ipft
      !------------------------------------------------------------------------------------!


      !----- Find the total snow depth. ---------------------------------------------------!
      ksn = csite%nlev_sfcwater(ipa)
      csite%total_sfcw_depth(ipa) = 0.
      total_sfcw_mass             = 0.
      do k=1,ksn
         csite%total_sfcw_depth(ipa) = csite%total_sfcw_depth(ipa)                         &
                                     + csite%sfcwater_depth(k,ipa)
         total_sfcw_mass             = total_sfcw_mass                                     &
                                     + csite%sfcwater_mass(k,ipa)
      end do
      !------------------------------------------------------------------------------------!

      !----- Reset properties. ------------------------------------------------------------!
      csite%veg_height      (ipa) = 0.0
      weight_sum                  = 0.0
      csite%opencan_frac    (ipa) = 1.0
      !------------------------------------------------------------------------------------!


      cpatch => csite%patch(ipa)

      !----- Loop over cohorts and integrate the patch-level properties. ------------------!
      do ico = 1,cpatch%ncohorts

         ipft = cpatch%pft(ico)


         !----- Compute the patch-level above-ground biomass
         csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                         &
                                     + ed_biomass(cpatch, ico) * cpatch%nplant(ico)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Compute average vegetation height, weighting using basal area.  We add the  !
         ! cohorts only until when the canopy is closed, this way we will not bias the     !
         ! vegetation height or the canopy depth towards the cohorts that live in the      !
         ! understorey.  Also, we must take into account the depth of the temporary        !
         ! surface water or snow, because this will make the plants "shorter".             !
         !---------------------------------------------------------------------------------!
         if (csite%opencan_frac(ipa) > 0.0) then
            weight                  = cpatch%nplant(ico) * cpatch%basarea(ico)
            weight_sum              = weight_sum + weight
            csite%veg_height(ipa)   = csite%veg_height(ipa) + cpatch%hite(ico) * weight
            csite%opencan_frac(ipa) = csite%opencan_frac(ipa)                              &
                                    * (1.0 - cpatch%crown_area(ico))
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !----- Normalise the vegetation height, making sure that it is above the minimum. ---!
      if (weight_sum > tiny(1.0)) then
         csite%veg_height(ipa)  = max(veg_height_min,csite%veg_height(ipa) / weight_sum)
      else
         csite%veg_height(ipa)  = veg_height_min
      end if
      !------------------------------------------------------------------------------------!



      !----- Find the patch roughness due to vegetation. ----------------------------------!
      csite%veg_rough(ipa) = vh2vr * csite%veg_height(ipa)
      !------------------------------------------------------------------------------------!



      !----- Find the 0-plane displacement due to vegetation. -----------------------------!
      csite%veg_displace(ipa) = vh2dh * csite%veg_height(ipa)
      !------------------------------------------------------------------------------------!



      !----- Update the canopy depth, and impose the minimum if needed be. ----------------!
      old_can_depth        = csite%can_depth(ipa)
      csite%can_depth(ipa) = max(csite%veg_height(ipa), minimum_canopy_depth)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the changes in canopy air space storage due to the change in canopy       !
      ! depth.  We update the budget variables during the patch and cohort dynamics step,  !
      ! but we skip this step during the initialisation, when we are creating a new patch  !
      ! or when we are updating (e.g. disturbance).                                        !
      !------------------------------------------------------------------------------------!
      if (update_zcaneff) then
         !----- Compute specific enthalpy, which will be used to find changes in storage. -!
         can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
         !---------------------------------------------------------------------------------!


         !------ Find the volume change of the canopy air space. --------------------------!
         fdelta_storage                  = frqsumi * csite%can_rhos(ipa)                   &
                                         * (csite%can_depth(ipa) - old_can_depth)
         !---------------------------------------------------------------------------------!


         !------ Update the change in storage due to change in CAS capacity. --------------!
         csite%co2budget_zcaneffect(ipa) = csite%co2budget_zcaneffect(ipa)                 &
                                         + fdelta_storage * csite%can_co2(ipa)             &
                                         * mmdryi
         csite%cbudget_zcaneffect  (ipa) = csite%cbudget_zcaneffect(ipa)                   &
                                         + fdelta_storage * csite%can_co2(ipa)             &
                                         * mmdryi * umol_2_kgC
         csite%wbudget_zcaneffect  (ipa) = csite%wbudget_zcaneffect(ipa)                   &
                                         + fdelta_storage * csite%can_shv(ipa)
         csite%ebudget_zcaneffect  (ipa) = csite%ebudget_zcaneffect(ipa)                   &
                                         + fdelta_storage * can_enthalpy
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the fraction of the canopy covered in snow.  I could not find any         !
      ! reference for the original method (commented out), so I implemented the method     !
      ! used in CLM-4, which is based on:                                                  !
      !                                                                                    !
      ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover  !
      !    fraction and its evaluation over large North American river basins,             !
      !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                         !
      !------------------------------------------------------------------------------------!
      ! csite%snowfac(ipa) = min(0.99, csite%total_sfcw_depth(ipa)/csite%veg_height(ipa))
      if (total_sfcw_mass > tiny_sfcwater_mass) then
         bulk_sfcw_dens     = max( fsdns, min( wdns                                        &
                                 , total_sfcw_mass / csite%total_sfcw_depth(ipa)))
         csite%snowfac(ipa) = max( 0.0, min( 0.99                                          &
                                 , tanh( csite%total_sfcw_depth(ipa)                       &
                                       / ( ny07_eq04_a * soil_rough                        &
                                         * (bulk_sfcw_dens * fsdnsi) ** ny07_eq04_m ) ) ) )
      else
         csite%snowfac(ipa) = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !----- Find the PFT-dependent size distribution of this patch. ----------------------!
      call patch_pft_size_profile(csite,ipa)
      !------------------------------------------------------------------------------------!


      !----- Update the cohort count (may be redundant as well...) ------------------------!
      csite%cohort_count(ipa) = cpatch%ncohorts
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_patch_derived_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine will take care of some diagnostic thermodynamic properties.      !
   !---------------------------------------------------------------------------------------!
   subroutine update_patch_thermo_props(csite,ipaa,ipaz,mzg,mzs,ntext_soil)
     
      use ed_state_vars, only : sitetype         ! ! structure
      use therm_lib    , only : idealdenssh      & ! function
                              , press2exner      & ! function
                              , extheta2temp     & ! function
                              , uextcm2tl        & ! function
                              , uint2tl          ! ! function
      use consts_coms  , only : t00              & ! intent(in)
                              , wdns             ! ! intent(in)
      use soil_coms    , only : soil             & ! intent(in)
                              , matric_potential ! ! function
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)                , target     :: csite
      integer                       , intent(in) :: ipaa
      integer                       , intent(in) :: ipaz
      integer                       , intent(in) :: mzg
      integer                       , intent(in) :: mzs
      integer       , dimension(mzg), intent(in) :: ntext_soil
      !----- Local variables. -------------------------------------------------------------!
      integer                                    :: ipa
      integer                                    :: nsoil
      integer                                    :: ksn
      integer                                    :: k
      real                                       :: soilhcap
      real                                       :: can_exner
      !------------------------------------------------------------------------------------!


      do ipa=ipaa,ipaz

         !----- Canopy air temperature and density. ---------------------------------------!
         can_exner           = press2exner (csite%can_prss(ipa))
         csite%can_temp(ipa) = extheta2temp(can_exner,csite%can_theta(ipa))
         csite%can_rhos(ipa) = idealdenssh ( csite%can_prss  (ipa)                         &
                                           , csite%can_temp  (ipa)                         &
                                           , csite%can_shv   (ipa)                         )
         !---------------------------------------------------------------------------------!


         !----- Update soil temperature and liquid water fraction. ------------------------!
         do k = 1, mzg
            nsoil    = ntext_soil(k)
            soilhcap = soil(nsoil)%slcpd
            call uextcm2tl(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*wdns,soilhcap  &
                          ,csite%soil_tempk(k,ipa),csite%soil_fracliq(k,ipa))
            csite%soil_mstpot(k,ipa) = matric_potential(nsoil,csite%soil_water(k,ipa))
         end do
         !---------------------------------------------------------------------------------!



         !----- Update temporary surface water temperature and liquid water fraction. -----!
         ksn = csite%nlev_sfcwater(ipa)
         csite%total_sfcw_depth(ipa) = 0.
         do k = 1, ksn
            call uint2tl(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa)          &
                        ,csite%sfcwater_fracliq(k,ipa))
            csite%total_sfcw_depth(ipa) =  csite%total_sfcw_depth(ipa)                     &
                                        +  csite%sfcwater_depth(k,ipa)
         end do
         do k = ksn+1,mzs
            if (k == 1) then
               csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk  (mzg,ipa)
               csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(mzg,ipa)
            else
               csite%sfcwater_tempk  (k,ipa) = csite%sfcwater_tempk  (k-1,ipa)
               csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            end if
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_patch_thermo_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine will update the fast mean properties, similarly to the routine   !
   ! above.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine update_patch_thermo_fmean(csite,ipaa,ipaz,mzg,ntext_soil)
     
      use ed_state_vars, only : sitetype           ! ! structure
      use therm_lib    , only : idealdenssh        & ! function
                              , press2exner        & ! function
                              , extheta2temp       & ! function
                              , uextcm2tl          & ! function
                              , uint2tl            ! ! function
      use consts_coms  , only : t00                & ! intent(in)
                              , wdns               ! ! intent(in)
      use soil_coms    , only : soil               & ! intent(in)
                              , tiny_sfcwater_mass & ! intent(in)
                              , matric_potential   ! ! function
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)                , target     :: csite
      integer                       , intent(in) :: ipaa
      integer                       , intent(in) :: ipaz
      integer                       , intent(in) :: mzg
      integer       , dimension(mzg), intent(in) :: ntext_soil
      !----- Local variables. -------------------------------------------------------------!
      integer                                    :: ipa
      integer                                    :: nsoil
      integer                                    :: k
      real                                       :: soilhcap
      real                                       :: can_exner
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      do ipa=ipaa,ipaz

         !----- Canopy air temperature and density. ---------------------------------------!
         can_exner                 = press2exner (csite%fmean_can_prss(ipa))
         csite%fmean_can_temp(ipa) = extheta2temp(can_exner,csite%fmean_can_theta(ipa))
         csite%fmean_can_rhos(ipa) = idealdenssh ( csite%fmean_can_prss  (ipa)             &
                                                 , csite%fmean_can_temp  (ipa)             &
                                                 , csite%fmean_can_shv   (ipa)             )
         !---------------------------------------------------------------------------------!


         !----- Update soil temperature and liquid water fraction. ------------------------!
         do k = 1, mzg
            nsoil    = ntext_soil(k)
            soilhcap = soil(nsoil)%slcpd
            call uextcm2tl( csite%fmean_soil_energy(k,ipa)                                 &
                          , csite%fmean_soil_water (k,ipa) * wdns                          &
                          , soilhcap                                                       &
                          , csite%fmean_soil_temp  (k,ipa)                                 &
                          , csite%fmean_soil_fliq  (k,ipa) )
            csite%fmean_soil_mstpot(k,ipa) =                                               &
                                    matric_potential(nsoil,csite%fmean_soil_water(k,ipa))
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !   If the patch had some temporary snow/pounding layer, convert the mean energy  !
         ! to J/kg, then find the mean temperature and liquid fraction.  Otherwise, set    !
         ! them to either zero or default values.                                          !
         !---------------------------------------------------------------------------------!
         if (csite%fmean_sfcw_mass(ipa) > tiny_sfcwater_mass) then
            csite%fmean_sfcw_energy(ipa) = csite%fmean_sfcw_energy(ipa)                    &
                                         / csite%fmean_sfcw_mass(ipa)
            call uint2tl(csite%fmean_sfcw_energy(ipa),csite%fmean_sfcw_temp(ipa)           &
                        ,csite%fmean_sfcw_fliq(ipa))
         else
            csite%fmean_sfcw_mass  (ipa)  = 0.
            csite%fmean_sfcw_depth (ipa)  = 0.
            csite%fmean_sfcw_energy(ipa)  = 0.
            csite%fmean_sfcw_temp  (ipa)  = csite%fmean_soil_temp(mzg,ipa)
            csite%fmean_sfcw_fliq  (ipa)  = csite%fmean_soil_fliq(mzg,ipa)
         end if
         !---------------------------------------------------------------------------------!
      end do
      return
   end subroutine update_patch_thermo_fmean
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will update the derived properties at the site level.             !
   !---------------------------------------------------------------------------------------!
   subroutine update_site_derived_props(cpoly,census_flag,isi)
     
      use ed_state_vars , only : polygontype  & ! structure
                               , sitetype     & ! structure
                               , patchtype    ! ! structure
      use consts_coms   , only : pio4         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype) , target     :: cpoly
      integer           , intent(in) :: census_flag
      integer           , intent(in) :: isi
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: bdbh
      integer                        :: ipa
      integer                        :: ico
      integer                        :: ipft
      !------------------------------------------------------------------------------------!
      
      !----- Initialise the variables before looping. -------------------------------------!
      cpoly%basal_area(:,:,isi) = 0.0
      cpoly%agb       (:,:,isi) = 0.0

      csite => cpoly%site(isi)

      !----- Loop over patches. -----------------------------------------------------------!
      do ipa = 1,csite%npatches

         cpatch => csite%patch(ipa)

         !----- Loop over cohorts. --------------------------------------------------------!
         do ico = 1,cpatch%ncohorts
            ipft = cpatch%pft(ico)

            !----- Update basal area and above-ground biomass. ----------------------------!
            if (census_flag == 0 .or. cpatch%first_census(ico) == 1) then
               bdbh = max(0,min( int(cpatch%dbh(ico) * 0.1), 10)) + 1

               cpoly%basal_area(ipft,bdbh,isi) = cpoly%basal_area(ipft, bdbh,isi)          &
                                               + cpatch%basarea(ico) * cpatch%nplant(ico)  &
                                               * csite%area(ipa)   
               cpoly%agb(ipft,bdbh,isi)        = cpoly%agb(ipft, bdbh,isi)                 &
                                               + cpatch%agb(ico)     * cpatch%nplant(ico)  &
                                               * csite%area(ipa)
            end if
         end do
      end do
      
      return
   end subroutine update_site_derived_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following subroutine finds the polygon averages from site-, patch-, and       !
   ! cohort--level properties whose time step is longer than DTLSM (days, months, years).  !
   ! Fluxes, meteorological input, thermodynamic properties, and radiation are aggregated  !
   ! in  sub-routine aggregate_polygon_fmean.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine update_polygon_derived_props(cgrid)
      use ed_state_vars         , only : edtype             & ! structure
                                       , polygontype        & ! structure
                                       , sitetype           & ! structure
                                       , patchtype          ! ! structure
      use soil_coms             , only : dslz               ! ! intent(in)
      use grid_coms             , only : nzg                ! ! intent(in)
      use ed_max_dims           , only : n_dbh              ! ! intent(in)
      use ed_misc_coms          , only : ddbhi              ! ! intent(in)
      use pft_coms              , only : c2n_leaf           & ! intent(in)
                                       , c2n_stem           & ! intent(in)
                                       , c2n_storage        & ! intent(in)
                                       , c2n_recruit        ! ! intent(in)
      use consts_coms           , only : tiny_num           ! ! intent(in)

      implicit none
      !----- Arguments.      --------------------------------------------------------------!
      type(edtype)         , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)    , pointer :: cpoly
      type(sitetype)       , pointer :: csite
      type(patchtype)      , pointer :: cpatch
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: p
      integer                        :: d
      integer                        :: k
      real                           :: poly_area_i
      real                           :: site_area_i
      real                           :: site_wgt
      real                           :: patch_wgt
      real                           :: rdepth
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      ! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   !
      !------------------------------------------------------------------------------------!
      !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.     !
      ! This works in off-line runs, but it causes memory leaks (and crashes) in the       !
      ! coupled runs over the ocean, where cgrid%npolygons can be 0 if one of the sub-     !
      ! -domains falls entirely over the ocean.  Thanks!                                   !
      !------------------------------------------------------------------------------------!
      ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
      !                 !      safe way of initialising the variable.
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !                                                                                 !
         !     Initialise properties.                                                      !
         !                                                                                 !
         !     This is the right and safe place to initialise polygon-level (cgrid) vari-  !
         ! ables, so in case npolygons is zero this will not cause memory leaks.  I know,  !
         ! this never happens in off-line runs, but it is quite common in coupled runs...  !
         ! Whenever one of the nodes receives a sub-domain where all the points are over   !
         ! the ocean, ED will not assign any polygon in that sub-domain, which means that  !
         ! that node will have 0 polygons, and the variables cannot be allocated.  If you  !
         ! try to access the polygon level variable outside the loop, then the model       !
         ! crashes due to segmentation violation (a bad thing), whereas by putting the     !
         ! variables here both the off-line model and the coupled runs will work, because  !
         ! this loop will be skipped when there is no polygon.                             !
         !---------------------------------------------------------------------------------!
         ! cgrid%blah(ipy) = 0.0 ! <<- This way works for all cases. 
         !---------------------------------------------------------------------------------!
         cgrid%nplant              (:,:,ipy) = 0.0
         cgrid%agb                 (:,:,ipy) = 0.0
         cgrid%lai                 (:,:,ipy) = 0.0
         cgrid%wai                 (:,:,ipy) = 0.0
         cgrid%basal_area          (:,:,ipy) = 0.0
         cgrid%thbark              (:,:,ipy) = 0.0
         cgrid%bdeada              (:,:,ipy) = 0.0
         cgrid%bdeadb              (:,:,ipy) = 0.0
         cgrid%btimber             (:,:,ipy) = 0.0
         cgrid%balive              (:,:,ipy) = 0.0
         cgrid%bleaf               (:,:,ipy) = 0.0
         cgrid%broot               (:,:,ipy) = 0.0
         cgrid%bsapwooda           (:,:,ipy) = 0.0
         cgrid%bsapwoodb           (:,:,ipy) = 0.0
         cgrid%bbarka              (:,:,ipy) = 0.0
         cgrid%bbarkb              (:,:,ipy) = 0.0
         cgrid%bseeds              (:,:,ipy) = 0.0
         cgrid%byield              (:,:,ipy) = 0.0
         cgrid%bstorage            (:,:,ipy) = 0.0
         cgrid%bdeada_n            (:,:,ipy) = 0.0
         cgrid%bdeadb_n            (:,:,ipy) = 0.0
         cgrid%balive_n            (:,:,ipy) = 0.0
         cgrid%bleaf_n             (:,:,ipy) = 0.0
         cgrid%broot_n             (:,:,ipy) = 0.0
         cgrid%bsapwooda_n         (:,:,ipy) = 0.0
         cgrid%bsapwoodb_n         (:,:,ipy) = 0.0
         cgrid%bbarka_n            (:,:,ipy) = 0.0
         cgrid%bbarkb_n            (:,:,ipy) = 0.0
         cgrid%bseeds_n            (:,:,ipy) = 0.0
         cgrid%bstorage_n          (:,:,ipy) = 0.0
         cgrid%leaf_maintenance    (:,:,ipy) = 0.0
         cgrid%root_maintenance    (:,:,ipy) = 0.0
         cgrid%barka_maintenance   (:,:,ipy) = 0.0
         cgrid%barkb_maintenance   (:,:,ipy) = 0.0
         cgrid%leaf_drop           (:,:,ipy) = 0.0
         cgrid%fast_grnd_c             (ipy) = 0.0
         cgrid%fast_soil_c             (ipy) = 0.0
         cgrid%microbe_soil_c          (ipy) = 0.0
         cgrid%slow_soil_c             (ipy) = 0.0
         cgrid%passive_soil_c          (ipy) = 0.0
         cgrid%struct_grnd_c           (ipy) = 0.0
         cgrid%struct_grnd_l           (ipy) = 0.0
         cgrid%struct_soil_c           (ipy) = 0.0
         cgrid%struct_soil_l           (ipy) = 0.0
         cgrid%fast_grnd_n             (ipy) = 0.0
         cgrid%fast_soil_n             (ipy) = 0.0
         cgrid%struct_grnd_n           (ipy) = 0.0
         cgrid%struct_soil_n           (ipy) = 0.0
         cgrid%mineral_soil_n          (ipy) = 0.0
         !---------------------------------------------------------------------------------!
         !     Some of these variables are redundant with variables above.  Perhaps we     !
         ! should find a single way to report them.                                        !
         !---------------------------------------------------------------------------------!
         cgrid%Nbiomass_uptake         (ipy) = 0.0
         cgrid%Cleaf_litter_flux       (ipy) = 0.0
         cgrid%Croot_litter_flux       (ipy) = 0.0
         cgrid%Nleaf_litter_flux       (ipy) = 0.0
         cgrid%Nroot_litter_flux       (ipy) = 0.0
         cgrid%Ngross_min              (ipy) = 0.0
         cgrid%Nnet_min                (ipy) = 0.0
         !---------------------------------------------------------------------------------!




         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop over sites.                                                            !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i=1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop over patches.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!


               !----- Integrate soil properties. ------------------------------------------!
               cgrid%fast_grnd_c   (ipy) = cgrid%fast_grnd_c        (ipy)                  &
                                         + csite%fast_grnd_c        (ipa)                  &
                                         * patch_wgt
               cgrid%fast_soil_c   (ipy) = cgrid%fast_soil_c        (ipy)                  &
                                         + csite%fast_soil_c        (ipa)                  &
                                         * patch_wgt
               cgrid%struct_grnd_c (ipy) = cgrid%struct_grnd_c      (ipy)                  &
                                         + csite%structural_grnd_c  (ipa)                  &
                                         * patch_wgt
               cgrid%struct_soil_c (ipy) = cgrid%struct_soil_c      (ipy)                  &
                                         + csite%structural_soil_c  (ipa)                  &
                                         * patch_wgt
               cgrid%struct_grnd_l (ipy) = cgrid%struct_grnd_l      (ipy)                  &
                                         + csite%structural_grnd_l  (ipa)                  &
                                         * patch_wgt
               cgrid%struct_soil_l (ipy) = cgrid%struct_soil_l      (ipy)                  &
                                         + csite%structural_soil_l  (ipa)                  &
                                         * patch_wgt
               cgrid%microbe_soil_c(ipy) = cgrid%microbe_soil_c     (ipy)                  &
                                         + csite%microbial_soil_c   (ipa)                  &
                                         * patch_wgt
               cgrid%slow_soil_c   (ipy) = cgrid%slow_soil_c        (ipy)                  &
                                         + csite%slow_soil_c        (ipa)                  &
                                         * patch_wgt
               cgrid%passive_soil_c(ipy) = cgrid%passive_soil_c     (ipy)                  &
                                         + csite%passive_soil_c     (ipa)                  &
                                         * patch_wgt
               cgrid%fast_grnd_n   (ipy) = cgrid%fast_grnd_n        (ipy)                  &
                                         + csite%fast_grnd_n        (ipa)                  &
                                         * patch_wgt
               cgrid%fast_soil_n   (ipy) = cgrid%fast_soil_n        (ipy)                  &
                                         + csite%fast_soil_n        (ipa)                  &
                                         * patch_wgt
               cgrid%struct_grnd_n (ipy) = cgrid%struct_grnd_n      (ipy)                  &
                                         + csite%structural_grnd_n  (ipa)                  &
                                         * patch_wgt
               cgrid%struct_soil_n (ipy) = cgrid%struct_soil_n      (ipy)                  &
                                         + csite%structural_soil_n  (ipa)                  &
                                         * patch_wgt
               cgrid%mineral_soil_n(ipy) = cgrid%mineral_soil_n     (ipy)                  &
                                         + csite%mineralized_soil_n (ipa)                  &
                                         * patch_wgt
               !---------------------------------------------------------------------------!



               !----- Zero the root fraction (patch-level diagnostic). --------------------!
               csite%rootdense(:,ipa) = 0.
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Loop over cohorts.                                                    !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !----- Find the PFT and DBH class to which this cohort belongs. ---------!
                  p = cpatch%pft(ico)
                  d = max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbhi)))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Rooting fraction: step 1, find root biomass per cubic meter         !
                  !    broot*nplant/rooting_depth   [kg/plant]*[plant/m2]/[m]              !
                  !------------------------------------------------------------------------!
                  rdepth = sum(dslz(cpatch%krdepth(ico):nzg))
                  do k=cpatch%krdepth(ico),nzg
                     csite%rootdense(k,ipa) = csite%rootdense(k,ipa)                       &
                                            + cpatch%broot(ico)*cpatch%nplant(ico) / rdepth
                  end do
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Integrate cohort-based properties.  Make sure that the polygon-     !
                  ! -level gets the right units (i.e., no /plant, but /m2).                !
                  !------------------------------------------------------------------------!
                  cgrid%nplant           (p,d,ipy) = cgrid%nplant            (p,d,ipy)     &
                                                   + cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%lai              (p,d,ipy) = cgrid%lai               (p,d,ipy)     &
                                                   + cpatch%lai                  (ico)     &
                                                   * patch_wgt
                  cgrid%wai              (p,d,ipy) = cgrid%wai               (p,d,ipy)     &
                                                   + cpatch%wai                  (ico)     &
                                                   * patch_wgt
                  cgrid%agb              (p,d,ipy) = cgrid%agb               (p,d,ipy)     &
                                                   + cpatch%agb                  (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%basal_area       (p,d,ipy) = cgrid%basal_area        (p,d,ipy)     &
                                                   + cpatch%basarea              (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bdeada           (p,d,ipy) = cgrid%bdeada            (p,d,ipy)     &
                                                   + cpatch%bdeada               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bdeadb           (p,d,ipy) = cgrid%bdeadb            (p,d,ipy)     &
                                                   + cpatch%bdeadb               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%btimber          (p,d,ipy) = cgrid%btimber           (p,d,ipy)     &
                                                   + cpatch%btimber              (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%balive           (p,d,ipy) = cgrid%balive            (p,d,ipy)     &
                                                   + cpatch%balive               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bleaf            (p,d,ipy) = cgrid%bleaf             (p,d,ipy)     &
                                                   + cpatch%bleaf                (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%broot            (p,d,ipy) = cgrid%broot             (p,d,ipy)     &
                                                   + cpatch%broot                (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bsapwooda        (p,d,ipy) = cgrid%bsapwooda         (p,d,ipy)     &
                                                   + cpatch%bsapwooda            (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bsapwoodb        (p,d,ipy) = cgrid%bsapwoodb         (p,d,ipy)     &
                                                   + cpatch%bsapwoodb            (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bbarka           (p,d,ipy) = cgrid%bbarka            (p,d,ipy)     &
                                                   + cpatch%bbarka               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bbarkb           (p,d,ipy) = cgrid%bbarkb            (p,d,ipy)     &
                                                   + cpatch%bbarkb               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bseeds           (p,d,ipy) = cgrid%bseeds            (p,d,ipy)     &
                                                   + cpatch%bseeds               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%byield           (p,d,ipy) = cgrid%byield            (p,d,ipy)     &
                                                   + cpatch%byield               (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bstorage         (p,d,ipy) = cgrid%bstorage          (p,d,ipy)     &
                                                   + cpatch%bstorage             (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bdeada_n         (p,d,ipy) = cgrid%bdeada_n          (p,d,ipy)     &
                                                   + cpatch%bdeada               (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bdeadb_n         (p,d,ipy) = cgrid%bdeadb_n          (p,d,ipy)     &
                                                   + cpatch%bdeadb               (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%balive_n         (p,d,ipy) = cgrid%balive_n          (p,d,ipy)     &
                                                   + ( ( cpatch%bleaf            (ico)     &
                                                       + cpatch%broot            (ico) )   &
                                                       / c2n_leaf                  (p)     &
                                                     + ( cpatch%bsapwooda        (ico)     &
                                                       + cpatch%bsapwoodb        (ico)     &
                                                       + cpatch%bbarka           (ico)     &
                                                       + cpatch%bbarkb           (ico) )   &
                                                       / c2n_stem                  (p)   ) &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bleaf_n          (p,d,ipy) = cgrid%bleaf_n           (p,d,ipy)     &
                                                   + cpatch%bleaf                (ico)     &
                                                   / c2n_leaf                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%broot_n          (p,d,ipy) = cgrid%broot_n           (p,d,ipy)     &
                                                   + cpatch%broot                (ico)     &
                                                   / c2n_leaf                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bsapwooda_n      (p,d,ipy) = cgrid%bsapwooda_n       (p,d,ipy)     &
                                                   + cpatch%bsapwooda            (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bsapwoodb_n      (p,d,ipy) = cgrid%bsapwoodb_n       (p,d,ipy)     &
                                                   + cpatch%bsapwoodb            (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bbarka_n         (p,d,ipy) = cgrid%bbarka_n          (p,d,ipy)     &
                                                   + cpatch%bbarka               (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bbarkb_n         (p,d,ipy) = cgrid%bbarkb_n          (p,d,ipy)     &
                                                   + cpatch%bbarkb               (ico)     &
                                                   / c2n_stem                      (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bseeds_n         (p,d,ipy) = cgrid%bseeds_n          (p,d,ipy)     &
                                                   + cpatch%bseeds               (ico)     &
                                                   / c2n_recruit                   (p)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%bstorage_n       (p,d,ipy) = cgrid%bstorage_n        (p,d,ipy)     &
                                                   + cpatch%bstorage             (ico)     &
                                                   / c2n_storage                           &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%leaf_maintenance (p,d,ipy) = cgrid%leaf_maintenance  (p,d,ipy)     &
                                                   + cpatch%leaf_maintenance     (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%root_maintenance (p,d,ipy) = cgrid%root_maintenance  (p,d,ipy)     &
                                                   + cpatch%root_maintenance     (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%barka_maintenance(p,d,ipy) = cgrid%barka_maintenance (p,d,ipy)     &
                                                   + cpatch%barka_maintenance    (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%barkb_maintenance(p,d,ipy) = cgrid%barkb_maintenance (p,d,ipy)     &
                                                   + cpatch%barkb_maintenance    (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  cgrid%leaf_drop        (p,d,ipy) = cgrid%leaf_drop         (p,d,ipy)     &
                                                   + cpatch%leaf_drop            (ico)     &
                                                   * cpatch%nplant               (ico)     &
                                                   * patch_wgt
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Bark thickness is a weighted average.  Here we integrate the       !
                  ! thickness, multiplied by basal area, then outside the siteloop we      !
                  ! normalise by total basal area.                                         !
                  !------------------------------------------------------------------------!
                  cgrid%thbark          (p,d,ipy) = cgrid%thbark          (p,d,ipy)        &
                                                  + cpatch%thbark             (ico)        &
                                                  * cpatch%nplant             (ico)        &
                                                  * cpatch%basarea            (ico)        &
                                                  * patch_wgt
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !  CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! !
                  !------------------------------------------------------------------------!
                  !   These variables were originally in integrate_ed_daily_output_flux, I !
                  ! moved them to here just to be consistent (they are not "dmean"         !
                  ! variables).  However, they are slightly different than the original,   !
                  ! which were missing nplant (so the results were in not in units of      !
                  ! kg/m2/yr, but kg/plant/yr).                                            !
                  !------------------------------------------------------------------------!
                  cgrid%Cleaf_litter_flux(ipy) = cgrid%Cleaf_litter_flux(ipy)              &
                                               + cpatch%leaf_maintenance(ico)              &
                                               * cpatch%nplant          (ico)              &
                                               * patch_wgt
                  cgrid%Croot_litter_flux(ipy) = cgrid%Croot_litter_flux(ipy)              &
                                               + cpatch%root_maintenance(ico)              &
                                               * cpatch%nplant          (ico)              &
                                               * patch_wgt
                  cgrid%Nleaf_litter_flux(ipy) = cgrid%Cleaf_litter_flux(ipy)              &
                                               + cpatch%leaf_maintenance(ico) / c2n_leaf(p)&
                                               * cpatch%nplant          (ico)              &
                                               * patch_wgt
                  cgrid%Nroot_litter_flux(ipy) = cgrid%Croot_litter_flux(ipy)              &
                                               + cpatch%root_maintenance(ico) / c2n_leaf(p)&
                                               * cpatch%nplant          (ico)              &
                                               * patch_wgt
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Update the patch-related N budget variables.                           !
               !---------------------------------------------------------------------------!
               cgrid%Ngross_min     (ipy) = cgrid%Ngross_min(ipy)                          &
                                          + csite%mineralized_N_input   (ipa)   * patch_wgt
               cgrid%Ngross_min     (ipy) = cgrid%Ngross_min(ipy)                          &
                                          + ( csite%mineralized_N_input (ipa)              &
                                            - csite%mineralized_N_loss  (ipa) ) * patch_wgt
               cgrid%Nbiomass_uptake(ipy) = cgrid%Ngross_min(ipy)                          &
                                          + csite%total_plant_nitrogen_uptake(ipa)         &
                                          * patch_wgt
               !---------------------------------------------------------------------------!

            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Normalise bark thickness.                                                   !
         !---------------------------------------------------------------------------------!
         where (cgrid%basal_area(:,:,ipy) > tiny_num)
            cgrid%thbark(:,:,ipy) = cgrid%thbark(:,:,ipy) / cgrid%basal_area(:,:,ipy)
         elsewhere
            cgrid%thbark(:,:,ipy) = 0.
         end where
         !---------------------------------------------------------------------------------!


      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine update_polygon_derived_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will read the regular soil moisture and temperature dataset.       !
   !---------------------------------------------------------------------------------------!
   subroutine read_soil_moist_temp(cgrid)

      use ed_state_vars , only : edtype       & ! structure
                               , polygontype  & ! structure
                               , sitetype     & ! structure
                               , patchtype    ! ! structure
      use soil_coms     , only : soilstate_db & ! intent(in)
                               , soil         & ! intent(in)
                               , slz          ! ! intent(in)
      use consts_coms   , only : wdns         & ! intent(in)
                               , t3ple        ! ! intent(in)
      use therm_lib     , only : cmtl2uext    ! ! function
      use grid_coms     , only : nzg          & ! intent(in)
                               , nzs          ! ! intent(in)
      use ed_therm_lib  , only : ed_grndvap   ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid         ! Alias for current ED grid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly          ! Alias for current polygon
      type(sitetype)    , pointer    :: csite          ! Alias for current site
      type(patchtype)   , pointer    :: cpatch         ! Alias for current patch
      integer                        :: ntext          !
      integer                        :: ilat           !
      integer                        :: ilon           !
      integer                        :: ilatf          !
      integer                        :: ilonf          !
      integer                        :: nls            !
      integer                        :: nlsw1          !
      integer                        :: k              !
      integer                        :: ipy            !
      integer                        :: isi            !
      integer                        :: ipa            !
      logical                        :: l1             !
      real                           :: glat           !
      real                           :: glon           !
      real                           :: tmp1           !
      real                           :: tmp2           !
      real                           :: soilw1         !
      real                           :: soilw2         !
      !----- Local constants.  ------------------------------------------------------------!
      logical           , parameter  :: harvard_override = .false.
      integer           , parameter  :: nlon = 144
      integer           , parameter  :: nlat = 73
      real              , parameter  :: dlon = 2.5
      real              , parameter  :: dlat = 2.5
      !------------------------------------------------------------------------------------!

      !----- Check whether the dataset exists and crash the run if it doesn't. ------------!
      inquire(file=trim(soilstate_db),exist=l1)
      if (.not.l1) then
         write(unit=*,fmt='(a)') ' ISOILSTATEINIT is set to read the initial'
         write(unit=*,fmt='(a)') ' soil moisture and temperature from a file.  The file'
         write(unit=*,fmt='(a)') ' specified by SOILSTATE_DB, however, doesn''t exist!'
         call fatal_error('Soil database '//trim(soilstate_db)//' not found!'              &
                        &,'read_soil_moist_temp','update_derived_props.f90')
      end if

      open (unit=12,file=trim(soilstate_db),form='formatted',status='old',position='rewind')

      !---- Loop over latitude levels, from north pole then southwards. -------------------!
      latloop: do ilatf = 1,nlat

         !----- Loop over longitude levels, from Greenwich Meridian then eastwards. -------!
         lonloop: do ilonf = 1,nlon 
    
            !------------------------------------------------------------------------------!
            !     Read in reanalysis: two temperatures and moistures, corresponding to     !
            ! different depths.                                                            !
            ! + soilw1, soilw2 are relative porosities and thus range from [0-1].          !
            ! + tmp1, tmp2 are temperature in kelvin.                                      !
            !------------------------------------------------------------------------------!
            read (unit=12,fmt=*) tmp1,tmp2,soilw1,soilw2

            !----- Make sure the numbers make sense... ------------------------------------!
            if (tmp1 > 0.0 .and. tmp2 > 0.0 .and. soilw1 > 0.0 .and. soilw2 > 0.0) then

               !----- Loop over land points. ----------------------------------------------!
               polyloop: do ipy=1,cgrid%npolygons
                  cpoly => cgrid%polygon(ipy)

                  !----- Land point lat, lon. ---------------------------------------------!
                  glat = cgrid%lat(ipy)
                  glon = cgrid%lon(ipy)
                  
                  if(glon < 0.0) glon = glon + 360.0
                  
                  !----- Find reanalysis point corresponding to this land point. ----------!
                  if(glat >= 0.0)then
                     ilat = nint((90.0 - glat)/dlat) + 1
                  else
                     ilat = nlat - nint((90.0 - abs(glat))/dlat)
                  end if
                  ilon = int(glon/dlon) + 1
                  
                  !----- If we are at the right point, fill the array. --------------------!
                  if(ilat == ilatf .and. ilon == ilonf)then

                     !------ Loop over sites and patches. ---------------------------------!
                     siteloop: do isi=1,cpoly%nsites
                        csite => cpoly%site(isi)

                        patchloop: do ipa=1,csite%npatches
                           cpatch => csite%patch(ipa)

                           do k=1,nzg
                              ntext = cpoly%ntext_soil(k,isi)

                              if(abs(slz(k)) < 0.1)then
                                 csite%soil_tempk(k,ipa) = tmp1
                                 csite%soil_water(k,ipa) = max(soil(ntext)%soilcp          &
                                                              ,soilw1 * soil(ntext)%slmsts)
                              else
                                 csite%soil_tempk(k,ipa) = tmp2
                                 csite%soil_water(k,ipa) = max(soil(ntext)%soilcp          &
                                                              ,soilw2 * soil(ntext)%slmsts)
                              endif
                              if (csite%soil_tempk(k,ipa) > t3ple) then
                                 csite%soil_fracliq(k,ipa) = 1.0
                              elseif (csite%soil_tempk(k,ipa) < t3ple) then
                                 csite%soil_fracliq(k,ipa) = 0.0
                              else
                                 csite%soil_fracliq(k,ipa) = 0.5
                              end if
                              csite%soil_energy(k,ipa) =                                   &
                                               cmtl2uext( soil(ntext)%slcpd                &
                                                        , csite%soil_water(k,ipa)* wdns    &
                                                        , csite%soil_tempk(k,ipa)          &
                                                        , csite%soil_fracliq(k,ipa))
                           end do


                          !----- Initial condition is with no snow/pond. ------------------!
                          csite%nlev_sfcwater(ipa)    = 0
                          csite%total_sfcw_depth(ipa) = 0.
                           do k=1,nzs
                              csite%sfcwater_energy (k,ipa) = 0.
                              csite%sfcwater_depth  (k,ipa) = 0.
                              csite%sfcwater_mass   (k,ipa) = 0.
                              csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
                              csite%sfcwater_fracliq(k,ipa) =                              &
                                                            csite%sfcwater_fracliq(nzg,ipa)
                           end do

                           if(harvard_override)then
                              csite%soil_tempk(1,ipa)     = 277.6
                              csite%soil_tempk(2:4,ipa)   = 276.0
                              csite%soil_energy(1,ipa)    =   1.5293664e8
                              csite%soil_energy(2,ipa)    =   1.4789957e8
                              csite%soil_energy(3:4,ipa)  =   1.4772002e8
                              csite%soil_water(1:4,ipa)   =   0.41595e+0
                              csite%soil_fracliq(1:4,ipa) =   1.0
                           endif
                           
                           !----- Compute the ground specific humidity. -------------------!
                           ntext = cpoly%ntext_soil(k,isi)
                           nls   = csite%nlev_sfcwater(ipa)
                           nlsw1 = max(1,nls)
                           call ed_grndvap(nls,ntext,csite%soil_water(nzg,ipa)             &
                                          ,csite%soil_tempk(nzg,ipa)                       &
                                          ,csite%soil_fracliq(nzg,ipa)                     &
                                          ,csite%sfcwater_tempk(nlsw1,ipa)                 &
                                          ,csite%snowfac(ipa),csite%can_prss(ipa)          &
                                          ,csite%can_shv(ipa),csite%ground_shv(ipa)        &
                                          ,csite%ground_ssh(ipa),csite%ground_temp(ipa)    &
                                          ,csite%ground_fliq(ipa),csite%ggsoil(ipa))

                        end do patchloop
                     end do siteloop
                  end if
               end do polyloop
            end if
         end do lonloop
      end do latloop

      close(unit=12,status='keep')
      return

   end subroutine read_soil_moist_temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will convert the integrated number of time steps in steps/day,    !
   ! then it will update the monthly mean workload.                                        !
   !---------------------------------------------------------------------------------------!
   subroutine update_workload(cgrid)
      use ed_state_vars, only : edtype        ! ! structure
      use ed_misc_coms , only : current_time  & ! intent(in)
                              , simtime       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target     :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(simtime)            :: lastmonth
      integer                  :: lmon
      integer                  :: ipy
      real                     :: ndaysi
      !------------------------------------------------------------------------------------!

      !----- Find last month information. -------------------------------------------------!
      call lastmonthdate(current_time,lastmonth,ndaysi)
      lmon = lastmonth%month

      !------------------------------------------------------------------------------------!
      !     Loop over all polygons, normalise the workload, then copy it to the cor-       !
      ! responding month.  Then copy the scratch column (13) to the appropriate month, and !
      ! reset it.                                                                          !
      !------------------------------------------------------------------------------------!
      do ipy=1,cgrid%npolygons
         cgrid%workload(lmon,ipy) = cgrid%workload(13,ipy) * ndaysi
         cgrid%workload(13,ipy)   = 0.
      end do

      return
   end subroutine update_workload
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine updates the model time.  It was moved to here to reduce the      !
   ! number of routines that must be written twice (off-line and coupled model).           !
   !---------------------------------------------------------------------------------------!
   subroutine update_model_time_dm(ctime,dtlong)

      use ed_misc_coms, only : simtime ! ! variable type
      use consts_coms , only : day_sec & ! intent(in)
                             , hr_sec  & ! intent(in)
                             , min_sec ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(simtime)         , intent(inout) :: ctime  ! Current time
      real                  , intent(in)    :: dtlong ! Model time step
      !----- Local variables. -------------------------------------------------------------!
      integer               , dimension(12) :: daymax
      !----- External functions. ----------------------------------------------------------!
      logical               , external      :: isleap ! This function check for leap years
      !------------------------------------------------------------------------------------!



      !----- Assume that the year is not leap. --------------------------------------------!
      daymax=(/31,28,31,30,31,30,31,31,30,31,30,31/)
      !------------------------------------------------------------------------------------!


      !----- Update the time. -------------------------------------------------------------!
      ctime%time = ctime%time + dtlong
      !------------------------------------------------------------------------------------!

      if (ctime%time >= day_sec)then

         !----- If time is greater than one day, update the day. --------------------------!
         ctime%time = ctime%time - day_sec
         ctime%date = ctime%date + 1

         !----- If the year is leap, correct February's number of days. -------------------!
         if (isleap(ctime%year)) daymax(2) = 29

         !----- If we have reached the end of the month, update the month. ----------------!
         if (ctime%date > daymax(ctime%month)) then
            ctime%date  = 1
            ctime%month = ctime%month + 1

            !------ If we have reached the end of the year, update the year. --------------!
            if (ctime%month == 13) then
               ctime%month = 1
               ctime%year = ctime%year + 1
            end if
         end if

      elseif (ctime%time < 0.0) then
         !----- Time is negative, go back one day. ----------------------------------------!
         ctime%time = ctime%time + day_sec
         ctime%date = ctime%date - 1

         !----- Day is zero, which means we must go to the previous month. ----------------!
         if (ctime%date == 0) then
            ctime%month = ctime%month - 1

            !----- Month is zero, which means that we must go to the previous year. -------!
            if (ctime%month == 0) then
               ctime%month = 12
               ctime%year = ctime%year - 1
            end if

            !------------------------------------------------------------------------------!
            !     Fix the month in case it it a leap year, and make the day the last of    !
            ! the previous month.                                                          !
            !------------------------------------------------------------------------------!
            if (isleap(ctime%year)) daymax(2) = 29
            ctime%date = daymax(ctime%month)
         end if
      end if
      !------------------------------------------------------------------------------------!



      !----- Update the hours, minutes, and seconds. --------------------------------------!
      ctime%hour = floor(ctime%time / hr_sec)
      ctime%min  = floor((ctime%time - real(ctime%hour)*hr_sec)/min_sec)
      ctime%sec  = floor(ctime%time - real(ctime%hour)*hr_sec - real(ctime%min)*min_sec)
      !------------------------------------------------------------------------------------!


      return
   end subroutine update_model_time_dm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       This subroutine scales all the "extensive" cohort properties when cohort area   !
   ! or population changes.                                                                !
   ! IMPORTANT: Only cohort-level variables that have units per area (m2 ground) should be !
   !            rescaled.  Variables whose units are per plant, m2 leaf, or m2 wood        !
   !            SHOULD NOT be included here.                                               !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_extensive_props(cpatch,aco,zco,mult)
      use ed_state_vars, only : patchtype    ! ! structure
      use ed_misc_coms , only : writing_long & ! intent(in)
                              , writing_eorq & ! intent(in)
                              , writing_dcyc ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target     :: cpatch    ! Current patch
      integer        , intent(in) :: aco       ! First cohort to be rescaled
      integer        , intent(in) :: zco       ! Last  cohort to be rescaled
      real           , intent(in) :: mult      ! Scale factor
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: ico       ! Cohort counter
      real                        :: mult_2    ! Square of the scale factor
      !------------------------------------------------------------------------------------!


      !----- Set up the scale factor for monthly mean sum of squares. ---------------------!
      mult_2 = mult * mult
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over all cohorts.                                                         !
      !------------------------------------------------------------------------------------!
      cohortloop: do ico=aco,zco
         cpatch%lai                (ico) = cpatch%lai                (ico) * mult
         cpatch%wai                (ico) = cpatch%wai                (ico) * mult
         cpatch%nplant             (ico) = cpatch%nplant             (ico) * mult
         cpatch%today_gpp          (ico) = cpatch%today_gpp          (ico) * mult
         cpatch%today_nppleaf      (ico) = cpatch%today_nppleaf      (ico) * mult
         cpatch%today_nppfroot     (ico) = cpatch%today_nppfroot     (ico) * mult
         cpatch%today_nppsapwood   (ico) = cpatch%today_nppsapwood   (ico) * mult
         cpatch%today_nppbark      (ico) = cpatch%today_nppbark      (ico) * mult
         cpatch%today_nppcroot     (ico) = cpatch%today_nppcroot     (ico) * mult
         cpatch%today_nppseeds     (ico) = cpatch%today_nppseeds     (ico) * mult
         cpatch%today_nppwood      (ico) = cpatch%today_nppwood      (ico) * mult
         cpatch%today_nppdaily     (ico) = cpatch%today_nppdaily     (ico) * mult
         cpatch%today_gpp_pot      (ico) = cpatch%today_gpp_pot      (ico) * mult
         cpatch%today_gpp_lightmax (ico) = cpatch%today_gpp_lightmax (ico) * mult
         cpatch%today_gpp_moistmax (ico) = cpatch%today_gpp_moistmax (ico) * mult
         cpatch%today_gpp_mlmax    (ico) = cpatch%today_gpp_mlmax    (ico) * mult
         cpatch%today_leaf_resp    (ico) = cpatch%today_leaf_resp    (ico) * mult
         cpatch%today_root_resp    (ico) = cpatch%today_root_resp    (ico) * mult
         cpatch%gpp                (ico) = cpatch%gpp                (ico) * mult
         cpatch%leaf_respiration   (ico) = cpatch%leaf_respiration   (ico) * mult
         cpatch%root_respiration   (ico) = cpatch%root_respiration   (ico) * mult
         cpatch%monthly_dndt       (ico) = cpatch%monthly_dndt       (ico) * mult
         cpatch%leaf_water         (ico) = cpatch%leaf_water         (ico) * mult
         cpatch%leaf_hcap          (ico) = cpatch%leaf_hcap          (ico) * mult
         cpatch%leaf_energy        (ico) = cpatch%leaf_energy        (ico) * mult
         cpatch%wood_water         (ico) = cpatch%wood_water         (ico) * mult
         cpatch%wood_hcap          (ico) = cpatch%wood_hcap          (ico) * mult
         cpatch%wood_energy        (ico) = cpatch%wood_energy        (ico) * mult
         !----- Crown area shall not exceed 1. --------------------------------------------!
         cpatch%crown_area         (ico) = min(1.,cpatch%crown_area  (ico) * mult)
         !----- Fast-scale means. ---------------------------------------------------------!
         cpatch%fmean_leaf_energy   (ico) = cpatch%fmean_leaf_energy   (ico) * mult
         cpatch%fmean_leaf_water    (ico) = cpatch%fmean_leaf_water    (ico) * mult
         cpatch%fmean_leaf_hcap     (ico) = cpatch%fmean_leaf_hcap     (ico) * mult
         cpatch%fmean_wood_energy   (ico) = cpatch%fmean_wood_energy   (ico) * mult
         cpatch%fmean_wood_water    (ico) = cpatch%fmean_wood_water    (ico) * mult
         cpatch%fmean_wood_hcap     (ico) = cpatch%fmean_wood_hcap     (ico) * mult
         cpatch%fmean_water_supply  (ico) = cpatch%fmean_water_supply  (ico) * mult
         cpatch%fmean_par_l         (ico) = cpatch%fmean_par_l         (ico) * mult
         cpatch%fmean_par_l_beam    (ico) = cpatch%fmean_par_l_beam    (ico) * mult
         cpatch%fmean_par_l_diff    (ico) = cpatch%fmean_par_l_diff    (ico) * mult
         cpatch%fmean_rshort_l      (ico) = cpatch%fmean_rshort_l      (ico) * mult
         cpatch%fmean_rlong_l       (ico) = cpatch%fmean_rlong_l       (ico) * mult
         cpatch%fmean_sensible_lc   (ico) = cpatch%fmean_sensible_lc   (ico) * mult
         cpatch%fmean_vapor_lc      (ico) = cpatch%fmean_vapor_lc      (ico) * mult
         cpatch%fmean_transp        (ico) = cpatch%fmean_transp        (ico) * mult
         cpatch%fmean_intercepted_al(ico) = cpatch%fmean_intercepted_al(ico) * mult
         cpatch%fmean_wshed_lg      (ico) = cpatch%fmean_wshed_lg      (ico) * mult
         cpatch%fmean_rshort_w      (ico) = cpatch%fmean_rshort_w      (ico) * mult
         cpatch%fmean_rlong_w       (ico) = cpatch%fmean_rlong_w       (ico) * mult
         cpatch%fmean_rad_profile (:,ico) = cpatch%fmean_rad_profile (:,ico) * mult
         cpatch%fmean_sensible_wc   (ico) = cpatch%fmean_sensible_wc   (ico) * mult
         cpatch%fmean_vapor_wc      (ico) = cpatch%fmean_vapor_wc      (ico) * mult
         cpatch%fmean_intercepted_aw(ico) = cpatch%fmean_intercepted_aw(ico) * mult
         cpatch%fmean_wshed_wg      (ico) = cpatch%fmean_wshed_wg      (ico) * mult
         !----- Daily means. --------------------------------------------------------------!
         if (writing_long) then
            cpatch%dmean_leaf_energy   (ico) = cpatch%dmean_leaf_energy   (ico) * mult
            cpatch%dmean_leaf_water    (ico) = cpatch%dmean_leaf_water    (ico) * mult
            cpatch%dmean_leaf_hcap     (ico) = cpatch%dmean_leaf_hcap     (ico) * mult
            cpatch%dmean_wood_energy   (ico) = cpatch%dmean_wood_energy   (ico) * mult
            cpatch%dmean_wood_water    (ico) = cpatch%dmean_wood_water    (ico) * mult
            cpatch%dmean_wood_hcap     (ico) = cpatch%dmean_wood_hcap     (ico) * mult
            cpatch%dmean_water_supply  (ico) = cpatch%dmean_water_supply  (ico) * mult
            cpatch%dmean_par_l         (ico) = cpatch%dmean_par_l         (ico) * mult
            cpatch%dmean_par_l_beam    (ico) = cpatch%dmean_par_l_beam    (ico) * mult
            cpatch%dmean_par_l_diff    (ico) = cpatch%dmean_par_l_diff    (ico) * mult
            cpatch%dmean_rshort_l      (ico) = cpatch%dmean_rshort_l      (ico) * mult
            cpatch%dmean_rlong_l       (ico) = cpatch%dmean_rlong_l       (ico) * mult
            cpatch%dmean_sensible_lc   (ico) = cpatch%dmean_sensible_lc   (ico) * mult
            cpatch%dmean_vapor_lc      (ico) = cpatch%dmean_vapor_lc      (ico) * mult
            cpatch%dmean_transp        (ico) = cpatch%dmean_transp        (ico) * mult
            cpatch%dmean_intercepted_al(ico) = cpatch%dmean_intercepted_al(ico) * mult
            cpatch%dmean_wshed_lg      (ico) = cpatch%dmean_wshed_lg      (ico) * mult
            cpatch%dmean_rshort_w      (ico) = cpatch%dmean_rshort_w      (ico) * mult
            cpatch%dmean_rlong_w       (ico) = cpatch%dmean_rlong_w       (ico) * mult
            cpatch%dmean_rad_profile (:,ico) = cpatch%dmean_rad_profile (:,ico) * mult
            cpatch%dmean_sensible_wc   (ico) = cpatch%dmean_sensible_wc   (ico) * mult
            cpatch%dmean_vapor_wc      (ico) = cpatch%dmean_vapor_wc      (ico) * mult
            cpatch%dmean_intercepted_aw(ico) = cpatch%dmean_intercepted_aw(ico) * mult
            cpatch%dmean_wshed_wg      (ico) = cpatch%dmean_wshed_wg      (ico) * mult
         end if
         !----- Monthly means. ------------------------------------------------------------!
         if (writing_eorq) then
            cpatch%mmean_lai           (ico) = cpatch%mmean_lai           (ico) * mult
            cpatch%mmean_leaf_energy   (ico) = cpatch%mmean_leaf_energy   (ico) * mult
            cpatch%mmean_leaf_water    (ico) = cpatch%mmean_leaf_water    (ico) * mult
            cpatch%mmean_leaf_hcap     (ico) = cpatch%mmean_leaf_hcap     (ico) * mult
            cpatch%mmean_wood_energy   (ico) = cpatch%mmean_wood_energy   (ico) * mult
            cpatch%mmean_wood_water    (ico) = cpatch%mmean_wood_water    (ico) * mult
            cpatch%mmean_wood_hcap     (ico) = cpatch%mmean_wood_hcap     (ico) * mult
            cpatch%mmean_water_supply  (ico) = cpatch%mmean_water_supply  (ico) * mult
            cpatch%mmean_par_l         (ico) = cpatch%mmean_par_l         (ico) * mult
            cpatch%mmean_par_l_beam    (ico) = cpatch%mmean_par_l_beam    (ico) * mult
            cpatch%mmean_par_l_diff    (ico) = cpatch%mmean_par_l_diff    (ico) * mult
            cpatch%mmean_rshort_l      (ico) = cpatch%mmean_rshort_l      (ico) * mult
            cpatch%mmean_rlong_l       (ico) = cpatch%mmean_rlong_l       (ico) * mult
            cpatch%mmean_sensible_lc   (ico) = cpatch%mmean_sensible_lc   (ico) * mult
            cpatch%mmean_vapor_lc      (ico) = cpatch%mmean_vapor_lc      (ico) * mult
            cpatch%mmean_transp        (ico) = cpatch%mmean_transp        (ico) * mult
            cpatch%mmean_intercepted_al(ico) = cpatch%mmean_intercepted_al(ico) * mult
            cpatch%mmean_wshed_lg      (ico) = cpatch%mmean_wshed_lg      (ico) * mult
            cpatch%mmean_rshort_w      (ico) = cpatch%mmean_rshort_w      (ico) * mult
            cpatch%mmean_rlong_w       (ico) = cpatch%mmean_rlong_w       (ico) * mult
            cpatch%mmean_rad_profile (:,ico) = cpatch%mmean_rad_profile (:,ico) * mult
            cpatch%mmean_sensible_wc   (ico) = cpatch%mmean_sensible_wc   (ico) * mult
            cpatch%mmean_vapor_wc      (ico) = cpatch%mmean_vapor_wc      (ico) * mult
            cpatch%mmean_intercepted_aw(ico) = cpatch%mmean_intercepted_aw(ico) * mult
            cpatch%mmean_wshed_wg      (ico) = cpatch%mmean_wshed_wg      (ico) * mult
            cpatch%mmsqu_gpp           (ico) = cpatch%mmsqu_gpp           (ico) * mult_2
            cpatch%mmsqu_npp           (ico) = cpatch%mmsqu_npp           (ico) * mult_2
            cpatch%mmsqu_plresp        (ico) = cpatch%mmsqu_plresp        (ico) * mult_2
            cpatch%mmsqu_sensible_lc   (ico) = cpatch%mmsqu_sensible_lc   (ico) * mult_2
            cpatch%mmsqu_vapor_lc      (ico) = cpatch%mmsqu_vapor_lc      (ico) * mult_2
            cpatch%mmsqu_transp        (ico) = cpatch%mmsqu_transp        (ico) * mult_2
            cpatch%mmsqu_sensible_wc   (ico) = cpatch%mmsqu_sensible_wc   (ico) * mult_2
            cpatch%mmsqu_vapor_wc      (ico) = cpatch%mmsqu_vapor_wc      (ico) * mult_2
         end if
         !----- Mean diel. ----------------------------------------------------------------!
         if (writing_dcyc) then
            cpatch%qmean_leaf_energy   (:,ico) = cpatch%qmean_leaf_energy   (:,ico)*mult
            cpatch%qmean_leaf_water    (:,ico) = cpatch%qmean_leaf_water    (:,ico)*mult
            cpatch%qmean_leaf_hcap     (:,ico) = cpatch%qmean_leaf_hcap     (:,ico)*mult
            cpatch%qmean_wood_energy   (:,ico) = cpatch%qmean_wood_energy   (:,ico)*mult
            cpatch%qmean_wood_water    (:,ico) = cpatch%qmean_wood_water    (:,ico)*mult
            cpatch%qmean_wood_hcap     (:,ico) = cpatch%qmean_wood_hcap     (:,ico)*mult
            cpatch%qmean_water_supply  (:,ico) = cpatch%qmean_water_supply  (:,ico)*mult
            cpatch%qmean_par_l         (:,ico) = cpatch%qmean_par_l         (:,ico)*mult
            cpatch%qmean_par_l_beam    (:,ico) = cpatch%qmean_par_l_beam    (:,ico)*mult
            cpatch%qmean_par_l_diff    (:,ico) = cpatch%qmean_par_l_diff    (:,ico)*mult
            cpatch%qmean_rshort_l      (:,ico) = cpatch%qmean_rshort_l      (:,ico)*mult
            cpatch%qmean_rlong_l       (:,ico) = cpatch%qmean_rlong_l       (:,ico)*mult
            cpatch%qmean_sensible_lc   (:,ico) = cpatch%qmean_sensible_lc   (:,ico)*mult
            cpatch%qmean_vapor_lc      (:,ico) = cpatch%qmean_vapor_lc      (:,ico)*mult
            cpatch%qmean_transp        (:,ico) = cpatch%qmean_transp        (:,ico)*mult
            cpatch%qmean_intercepted_al(:,ico) = cpatch%qmean_intercepted_al(:,ico)*mult
            cpatch%qmean_wshed_lg      (:,ico) = cpatch%qmean_wshed_lg      (:,ico)*mult
            cpatch%qmean_rshort_w      (:,ico) = cpatch%qmean_rshort_w      (:,ico)*mult
            cpatch%qmean_rlong_w       (:,ico) = cpatch%qmean_rlong_w       (:,ico)*mult
            cpatch%qmean_rad_profile (:,:,ico) = cpatch%qmean_rad_profile (:,:,ico)*mult
            cpatch%qmean_sensible_wc   (:,ico) = cpatch%qmean_sensible_wc   (:,ico)*mult
            cpatch%qmean_vapor_wc      (:,ico) = cpatch%qmean_vapor_wc      (:,ico)*mult
            cpatch%qmean_intercepted_aw(:,ico) = cpatch%qmean_intercepted_aw(:,ico)*mult
            cpatch%qmean_wshed_wg      (:,ico) = cpatch%qmean_wshed_wg      (:,ico)*mult
            cpatch%qmsqu_gpp           (:,ico) = cpatch%qmsqu_gpp           (:,ico)*mult_2
            cpatch%qmsqu_npp           (:,ico) = cpatch%qmsqu_npp           (:,ico)*mult_2
            cpatch%qmsqu_plresp        (:,ico) = cpatch%qmsqu_plresp        (:,ico)*mult_2
            cpatch%qmsqu_sensible_lc   (:,ico) = cpatch%qmsqu_sensible_lc   (:,ico)*mult_2
            cpatch%qmsqu_vapor_lc      (:,ico) = cpatch%qmsqu_vapor_lc      (:,ico)*mult_2
            cpatch%qmsqu_transp        (:,ico) = cpatch%qmsqu_transp        (:,ico)*mult_2
            cpatch%qmsqu_sensible_wc   (:,ico) = cpatch%qmsqu_sensible_wc   (:,ico)*mult_2
            cpatch%qmsqu_vapor_wc      (:,ico) = cpatch%qmsqu_vapor_wc      (:,ico)*mult_2
         end if
         !---------------------------------------------------------------------------------!

      end do cohortloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_cohort_extensive_props
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine patch_pft_size_profile(csite,ipa)
      use ed_state_vars        , only : sitetype   & ! structure
                                      , patchtype  ! ! structure
      use fusion_fission_coms  , only : ff_nhgt    & ! intent(in)
                                      , hgt_class  ! ! intent(in)
      use allometry            , only : size2bl    ! ! intent(in)
      use ed_max_dims          , only : n_pft      ! ! intent(in)
      use pft_coms             , only : hgt_min    & ! intent(in)
                                      , dbh_crit   & ! intent(in)
                                      , is_grass   ! ! intent(in)
      use ed_misc_coms         , only : igrass     ! ! intent(in)
      use canopy_radiation_coms, only : ihrzrad    & ! intent(in)
                                      , cci_hmax   ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target     :: csite     ! Current site
      integer                , intent(in) :: ipa       ! Current patch index
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch    ! Current patch
      integer                             :: ipft      ! PFT index
      integer                             :: ihgt      ! Height class index
      integer                             :: ico       ! Counters
      real(kind=4)                        :: lai_pot   ! Potential LAI
      real(kind=4)                        :: hgt_eff   ! Effective height
      real(kind=4)                        :: sz_fact   ! Size correction factor 
      !------------------------------------------------------------------------------------!


      !----- Reset all bins to zero. ------------------------------------------------------!
      do ipft=1,n_pft
         do ihgt=1,ff_nhgt
            csite%cumlai_profile(ipft,ihgt,ipa) = 0.0
         end do
      end do
      !------------------------------------------------------------------------------------!



      !----- Update bins ------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cohortloop: do ico = 1,cpatch%ncohorts

         !----- Find the PFT class. -------------------------------------------------------!
         ipft    = cpatch%pft(ico)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Check whether this cohort is almost at the minimum height given its PFT.    !
         ! If it is, then we will skip it.                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < hgt_min(ipft) + 0.2) cycle cohortloop
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Decide whether to use actual height or effective height (to account for     !
         ! emergent trees).                                                                !
         !---------------------------------------------------------------------------------!
         select case (ihrzrad)
         case (2,4)
            sz_fact = max(dbh_crit(ipft),cpatch%dbh(ico))/dbh_crit(ipft)
            hgt_eff = min(cci_hmax, cpatch%hite(ico) * sz_fact * sz_fact)
            ihgt    = min(ff_nhgt,max(1,count(hgt_class < hgt_eff)))
         case default
            ihgt    = min(ff_nhgt,max(1,count(hgt_class < cpatch%hite(ico))))
         end select
         !---------------------------------------------------------------------------------!


         !----- Find the potential (on-allometry) leaf area index. ------------------------!
         if (is_grass(ipft) .and. igrass==1) then
             !--use actual bleaf for grass
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico) * cpatch%bleaf(ico)
         else
             !--use dbh for trees
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico)                                &
                     * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         end if
         !---------------------------------------------------------------------------------!


         !----- Add the potential LAI to the bin. -----------------------------------------!
         csite%cumlai_profile(ipft,ihgt,ipa) = csite%cumlai_profile(ipft,ihgt,ipa)         &
                                             + lai_pot
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!



      !----- Integrate the leaf area index from top to bottom. ----------------------------!
      do ihgt=ff_nhgt-1,1,-1
         do ipft=1,n_pft
            csite%cumlai_profile(ipft,ihgt,ipa) = csite%cumlai_profile(ipft,ihgt  ,ipa)    &
                                                + csite%cumlai_profile(ipft,ihgt+1,ipa)
         end do
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine patch_pft_size_profile
   !=======================================================================================!
   !=======================================================================================!
end module update_derived_utils
